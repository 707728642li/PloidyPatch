from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
from pathlib import Path
from pathlib import PurePosixPath
import shutil

import pytest

from ploidypatch.holdout_contract import (
    ArtifactSource,
    HoldoutContract,
    ReferenceContract,
    TargetResolvedParameters,
    staged_relative_path,
)


ROOT = Path(__file__).parents[1]
SCRIPTS = {
    name: (ROOT / "scripts" / name).read_text(encoding="utf-8")
    for name in (
        "run_actinidia_miniprot_upstream_v0.5.sh",
        "build_actinidia_method_trio_candidate_pools_v0.5.sh",
        "run_actinidia_blind_union_self_wgd_v0.5.sh",
        "score_actinidia_candidates_blind_v0.5.sh",
        "run_actinidia_blind_pipeline_v0.5.sh",
    )
}
VERIFIER_PATH = ROOT / "scripts/verify_external_holdout_blind_context_v0.5.py"
VERIFIER_SPEC = importlib.util.spec_from_file_location("blind_context_v05", VERIFIER_PATH)
assert VERIFIER_SPEC and VERIFIER_SPEC.loader
VERIFIER = importlib.util.module_from_spec(VERIFIER_SPEC)
VERIFIER_SPEC.loader.exec_module(VERIFIER)


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


@pytest.mark.parametrize("name", sorted(SCRIPTS))
def test_candidate_scripts_bind_generic_freezes_and_three_blind_roles(name: str) -> None:
    script = SCRIPTS[name]
    for marker in (
        "PLOIDYPATCH_BLIND_RUNNER",
        "PLOIDYPATCH_NETWORK_ACCESS",
        "PLOIDYPATCH_STAGED_INPUT_ROOT",
        "PLOIDYPATCH_BLIND_BENCHMARK_ROOT",
        "PLOIDYPATCH_HOLDOUT_CONTRACT",
        "PLOIDYPATCH_PROTOCOL_FREEZE",
        "PLOIDYPATCH_EXECUTION_FREEZE",
        "role_manifest.tsv",
        "role_contract.json",
        "verify_external_holdout_blind_context_v0.5.py",
        "--expected-holdout-id actinidia_red5_v0.5",
        "--expected-primary-chromosomes 29",
    ):
        assert marker in script, (name, marker)
    assert "/nas_data" in script
    assert "evaluator_only" in script
    assert "target_complete" in script or "/complete" in script


@pytest.mark.parametrize("name", sorted(SCRIPTS))
def test_candidate_scripts_have_no_mechanical_species_or_model_residue(name: str) -> None:
    folded = SCRIPTS[name].casefold()
    for forbidden in (
        "populus",
        "salix",
        "ptr_v4",
        "20260930",
        "salicoid",
        "ploidypatch_ranker_v0.5",
        "ploidypatch.composite_ranker.v0.5",
        "untouched_confirmatory_external",
        "external_inputs/actinidia_v0.5",
        "preflight_input_manifest.tsv",
    ):
        assert forbidden not in folded, (name, forbidden)


def test_miniprot_upstream_uses_red5_29_and_provider_specific_hierarchies() -> None:
    script = SCRIPTS["run_actinidia_miniprot_upstream_v0.5.sh"]
    assert "actinidia_chinensis_red5.tsv" in script
    assert "len(allowed) != 29" in script
    assert "exactly 29 seqids" in script
    assert "shared_gene_transcript_id" in script
    assert "gene_with_direct_children" in script
    assert "normalize_maker_transcript_hierarchy_v0.5.py" in script
    assert "synthesize_missing_transcript_exons.py" in script
    assert "--repair-parent-bounds" in script
    assert "ploidypatch.missing_transcript_exon_compat.v2" in script
    assert 'report.get("child_coordinate_or_cds_changes") is not False' in script
    assert 'counts.get("parent_bounds_repaired") != len(repairs)' in script
    assert 'cds.get("input") != cds.get("output")' in script
    assert 'counts.get("input_cds_records") != counts.get("output_cds_records")' in script
    target_binding = script.index("normalized Red5 primary genome differs")
    hierarchy = script.index("normalize_maker_transcript_hierarchy_v0.5.py", target_binding)
    primary = script.index("normalize primary-annotation", hierarchy)
    synthesize = script.index("synthesize_missing_transcript_exons.py", primary)
    projection = script.index('"$miniprot_bin" -I', synthesize)
    assert target_binding < hierarchy < primary < synthesize < projection


def test_two_candidate_references_contribute_one_vote_per_method_family() -> None:
    script = SCRIPTS["build_actinidia_method_trio_candidate_pools_v0.5.sh"]
    assert 'for method in gemoma lifton; do' in script
    assert script.count('--candidate "actinidia_eriantha=${raw[') == 1
    assert script.count('--candidate "actinidia_rufa=${raw[') == 1
    assert script.count('--candidate "miniprot=$working_root/methods/miniprot') == 1
    assert script.count('--candidate "gemoma=$working_root/methods/gemoma') == 1
    assert script.count('--candidate "lifton=$working_root/methods/lifton') == 1
    assert "within_method_reference_vote_count\\t1" in script
    assert "one_method_family_vote" not in script or "within_method_reference_vote_count" in script


def test_self_wgd_is_blind_candidate_union_and_actinidia_event_only() -> None:
    script = SCRIPTS["run_actinidia_blind_union_self_wgd_v0.5.sh"]
    assert "consensus/primary_union/blind/candidate.gff3" in script
    assert "actinidia_specific_ad_alpha_blind_recomputed" in script
    assert "candidate_candidate_pair_policy\\treject_as_circular" in script
    assert "--min-block-pairs 20" in script
    assert "--max-target-seqs 20" in script
    assert "--evalue 1e-5" in script


def test_scorer_is_exact_frozen_v04_with_production_guard() -> None:
    script = SCRIPTS["score_actinidia_candidates_blind_v0.5.sh"]
    assert "PLOIDYPATCH_COMPOSITE_MODEL_FREEZE" in script
    assert "protocol_artifacts/config/actinidia_external_validation_policy_v0.5.tsv" in script
    assert 'policy.get("model_version") != "PloidyPatch_ranker_v0.4"' in script
    assert 'manifest.get("schema_version") != "ploidypatch.composite_ranker.v0.4"' in script
    assert "score-support-conditioned-candidates" in script
    assert "apply-conflict-winner-guard" in script
    assert "v04_automatic_approval" in script
    assert "automatic_approved\") != 0" in script
    assert "winner_mismatch_count\") != 0" in script


def test_blind_pipeline_has_one_positional_root_and_no_truth_mount_fallback() -> None:
    script = SCRIPTS["run_actinidia_blind_pipeline_v0.5.sh"]
    assert "if [[ $# -ne 1 ]]" in script
    assert "usage: $0 ISOLATED_PROJECT_ROOT" in script
    assert "PLOIDYPATCH_BLIND_OUTPUT_ROOT" in script
    assert "pipeline_blind_context.json" in script
    assert "for reference in actinidia_eriantha actinidia_rufa" in script
    assert "run_gemoma_homology.sh" in script
    assert "run_lifton_transfer.sh" in script
    assert "build_actinidia_method_trio_candidate_pools_v0.5.sh" in script
    assert "run_actinidia_blind_union_self_wgd_v0.5.sh" in script
    assert "score_actinidia_candidates_blind_v0.5.sh" in script
    target_gate = script.index("normalized Red5 target differs from sealed blind benchmark")
    first_homology = script.index("run_gemoma_homology.sh", target_gate)
    assert script.index("blind_manifest.json", target_gate - 500) < target_gate < first_homology


def test_pipeline_writes_atomic_truth_free_frozen_command_log() -> None:
    script = SCRIPTS["run_actinidia_blind_pipeline_v0.5.sh"]
    assert "command_log=$project_root/pipeline_commands.tsv" in script
    assert "command_log_working=${command_log}.working" in script
    assert "implementation_manifest.tsv" in script
    assert "frozen_relative_script" in script
    assert "script_sha256" in script
    assert 'mv "$command_log_working" "$command_log"' in script
    command_rows = script.split("<<'COMMANDS'\n", 1)[1].split("\nCOMMANDS", 1)[0]
    assert len(command_rows.splitlines()) == 10
    folded = command_rows.casefold()
    for forbidden in ("nas_data", "evaluator", "truth", "labels", "target_complete"):
        assert forbidden not in folded
    assert script.index('mv "$command_log_working" "$command_log"') < script.index(
        'bash "$code_root/scripts/run_actinidia_miniprot_upstream_v0.5.sh"'
    )


def test_formal_blind_output_relative_paths_are_materialized() -> None:
    pool = SCRIPTS["build_actinidia_method_trio_candidate_pools_v0.5.sh"]
    score = SCRIPTS["score_actinidia_candidates_blind_v0.5.sh"]
    assert "results/copy_collapse/external/actinidia_v0.5_method_trio" in pool
    assert "consensus/primary_union/blind/decisions.tsv" in score
    assert "consensus/primary_union/blind/candidate.gff3.manifest.json" in score
    assert "results/copy_collapse/external/actinidia/v0.5" not in score
    assert "results/copy_collapse/external/actinidia_v0.5_blind_rankings" in score
    assert 'scores/v04.tsv"' in score
    assert 'scores/v04.tsv.manifest.json"' in score


def test_shared_context_verifier_is_header_aware_and_fail_closed() -> None:
    verifier = VERIFIER_PATH.read_text(encoding="utf-8")
    assert "csv.DictReader" in verifier
    assert '"artifact"' in verifier
    assert '"staged_relative_path"' in verifier
    assert "Role manifest header differs" in verifier
    assert "len(rows) != 15" in verifier
    assert "visible_candidate_safe_artifacts" in verifier
    assert "len(visible_artifacts) != 7" in verifier
    assert '"ploidypatch.blind_benchmark_input.v0.5"' in verifier
    assert 'benchmark.get("truth_access") is not False' in verifier
    assert 'execution.get("network_access_in_blind_runner") is not False' in verifier
    assert '"ploidypatch.composite_ranker.v0.4"' in verifier


def test_shared_context_verifier_accepts_only_matching_synthetic_blind_roles(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    input_root = tmp_path / "holdout"
    protocol = tmp_path / "protocol"
    execution = tmp_path / "execution"
    benchmark = tmp_path / "blind_benchmark"
    for root in (input_root, protocol, execution, benchmark):
        root.mkdir()

    specifications = (
        ("target", "Target", "target_red5", "red5"),
        ("candidate_reference", "Eriantha", "candidate_eriantha", "aer"),
        ("candidate_reference", "Rufa", "candidate_rufa", "aru"),
        ("evaluator_reference", "Rhododendron", "evaluator_rhododendron", "rhs"),
        ("evaluator_reference", "Diospyros", "evaluator_diospyros", "dol"),
    )
    references = []
    rows = []
    for role, species, bundle, prefix in specifications:
        artifacts = {}
        for artifact_name in ("genome", "gff3", "protein"):
            payload = f"{species}:{artifact_name}\n".encode()
            artifact = ArtifactSource(
                source_relative_path=PurePosixPath(f"source/{species}/{artifact_name}.dat"),
                bytes=len(payload),
                sha256=hashlib.sha256(payload).hexdigest(),
            )
            artifacts[artifact_name] = artifact
        reference = ReferenceContract(
            role=role,
            species_id=species,
            release=f"{species}_release",
            bundle_id=bundle,
            wgdi_prefix=prefix,
            primary_seqid_table=PurePosixPath(f"config/{species}.tsv"),
            **artifacts,
        )
        references.append(reference)
        for artifact_name, artifact in reference.artifact_items():
            relative = staged_relative_path(reference, artifact_name).as_posix()
            visible = role == "candidate_reference" or (
                role == "target" and artifact_name == "genome"
            )
            if visible:
                destination = input_root / relative
                destination.parent.mkdir(parents=True, exist_ok=True)
                destination.write_bytes(f"{species}:{artifact_name}\n".encode())
            rows.append(
                {
                    "role": role,
                    "species_id": species,
                    "release": reference.release,
                    "bundle_id": bundle,
                    "wgdi_prefix": prefix,
                    "artifact": artifact_name,
                    "bytes": str(artifact.bytes),
                    "sha256": artifact.sha256,
                    "source_relative_path": artifact.source_relative_path.as_posix(),
                    "staged_relative_path": relative,
                    "staged_sha256": artifact.sha256,
                }
            )
    contract = HoldoutContract(
        schema_version="ploidypatch.external_holdout_contract.v0.5",
        holdout_id="actinidia_red5_v0.5",
        policy_id="ploidypatch_actinidia_external_validation_v0.5",
        test_role="target_level_predeclared_untouched_secondary_replication",
        model_version="PloidyPatch_ranker_v0.4",
        references=tuple(references),
        seeds=None,  # type: ignore[arg-type]
        target_resolved_parameters=TargetResolvedParameters(29, 0.75, 22),
        scientific_parameters=None,  # type: ignore[arg-type]
        truth_blind={},
    )
    contract_path = protocol / "contract.json"
    contract_path.write_text("{}\n", encoding="utf-8")
    role_manifest = input_root / "role_manifest.tsv"
    with role_manifest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=VERIFIER.ROLE_FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    role_contract = {
        "schema_version": VERIFIER.STAGE_SCHEMA,
        "holdout_id": contract.holdout_id,
        "policy_id": contract.policy_id,
        "test_role": contract.test_role,
        "model_version": contract.model_version,
        "code_commit": "a" * 40,
        "contract": {"sha256": _sha(contract_path)},
        "truth_blind": {},
        "target_resolved_parameters": {
            "primary_chromosome_count": 29,
            "minimum_target_chromosomes_fraction": 0.75,
            "minimum_target_chromosomes": 22,
        },
        "role_boundaries": {
            "candidate_only": "candidate_generation_only",
            "evaluator_only": "complete_target_annotation_and_evaluator_truth_references",
        },
    }
    (input_root / "role_contract.json").write_text(
        json.dumps(role_contract) + "\n", encoding="utf-8"
    )
    shutil.copyfile(role_manifest, protocol / "role_manifest.tsv")
    shutil.copyfile(input_root / "role_contract.json", protocol / "role_contract.json")
    (protocol / "SHA256SUMS").write_text("protocol freeze sentinel\n", encoding="utf-8")
    protocol_manifest = {
        "schema_version": VERIFIER.PROTOCOL_SCHEMA,
        "holdout_id": contract.holdout_id,
        "policy_id": contract.policy_id,
        "test_role": contract.test_role,
        "model_version": contract.model_version,
        "code_commit": "a" * 40,
        "contract_sha256": _sha(contract_path),
        "truth_access": False,
        "freeze_stage": "post_metadata_pre_pair_pre_candidate_pre_label",
        "wgd_pairs_enumerated": False,
        "candidate_counts_computed": False,
        "truth_labels_accessed": False,
        "staged_role_contract_sha256": _sha(input_root / "role_contract.json"),
        "staged_input_SHA256SUMS_sha256": "e" * 64,
    }
    (protocol / "protocol_manifest.json").write_text(
        json.dumps(protocol_manifest) + "\n", encoding="utf-8"
    )
    (execution / "SHA256SUMS").write_text("execution freeze sentinel\n", encoding="utf-8")
    execution_manifest = {
        "schema_version": VERIFIER.EXECUTION_SCHEMA,
        "holdout_id": contract.holdout_id,
        "policy_id": contract.policy_id,
        "test_role": contract.test_role,
        "model_version": contract.model_version,
        "code_commit": "a" * 40,
        "contract_sha256": _sha(contract_path),
        "freeze_stage": "post_metadata_pre_pair_pre_candidate_pre_label",
        "created_before": {
            "wgd_pair_enumeration": True,
            "candidate_generation": True,
            "candidate_labels": True,
            "candidate_scores": True,
        },
        "protocol_SHA256SUMS_sha256": _sha(protocol / "SHA256SUMS"),
        "network_access_in_blind_runner": False,
        "nas_data_mount_in_blind_runner": False,
        "complete_target_annotation_mount_in_blind_runner": False,
        "evaluator_only_mount_in_blind_runner": False,
        "truth_or_label_mount_in_blind_runner": False,
    }
    (execution / "execution_manifest.json").write_text(
        json.dumps(execution_manifest) + "\n", encoding="utf-8"
    )
    perturbed = benchmark / "perturbed.gff3"
    perturbed.write_text("##gff-version 3\n", encoding="utf-8")
    blind_manifest = {
        "schema_version": "ploidypatch.blind_benchmark_input.v0.5",
        "truth_access": False,
        "complete_target_annotation_access": False,
        "perturbed_annotation": {"file_name": perturbed.name, "sha256": _sha(perturbed)},
        "target_genome": {"mount_role": "shared_target_genome", "sha256": "a" * 64},
    }
    (benchmark / "blind_manifest.json").write_text(
        json.dumps(blind_manifest) + "\n", encoding="utf-8"
    )
    (benchmark / "SHA256SUMS").write_text(
        f"{_sha(benchmark / 'blind_manifest.json')}  blind_manifest.json\n"
        f"{_sha(perturbed)}  perturbed.gff3\n",
        encoding="utf-8",
    )
    monkeypatch.setenv("PLOIDYPATCH_BLIND_RUNNER", "1")
    monkeypatch.setenv("PLOIDYPATCH_NETWORK_ACCESS", "none")
    monkeypatch.setenv("PLOIDYPATCH_CODE_COMMIT", "a" * 40)
    monkeypatch.setattr(VERIFIER, "load_holdout_contract", lambda _: contract)
    # This unit test invokes the verifier directly on the host.  Production
    # invokes it inside bubblewrap, where the host /nas_data mount is absent.
    # Preserve the production fail-closed check while modelling that namespace
    # boundary explicitly in the synthetic unit test.
    original_exists = Path.exists

    def namespace_exists(path: Path) -> bool:
        if path == Path("/nas_data"):
            return False
        return original_exists(path)

    monkeypatch.setattr(Path, "exists", namespace_exists)
    report = VERIFIER.verify_context(
        input_root=input_root,
        contract_path=contract_path,
        protocol_root=protocol,
        execution_root=execution,
        blind_benchmark_root=benchmark,
        expected_holdout_id=contract.holdout_id,
        expected_primary_chromosomes=29,
    )
    assert len(report["visible_candidate_safe_artifacts"]) == 7
    assert report["network_access"] is False

    patch_commit = "b" * 40
    patch = {
        "schema_version": "ploidypatch.external_holdout_execution_patch.v0.5",
        "freeze_stage": (
            "post_evaluator_truth_failed_blind_pre_candidate_pre_score_pre_label_execution_patch"
        ),
        "base_code_commit": "a" * 40,
        "patch_code_commit": patch_commit,
        "superseded_execution_SHA256SUMS_sha256": "f" * 64,
        "base_protocol_SHA256SUMS_sha256": _sha(protocol / "SHA256SUMS"),
        "contract_sha256": _sha(contract_path),
        "staged_input_SHA256SUMS_sha256": "e" * 64,
        "scientific_protocol_changed": False,
        "contract_or_policy_changed": False,
        "model_or_threshold_changed": False,
        "staged_inputs_changed": False,
        "truth_or_benchmark_regenerated": False,
        "evaluator_truth_construction_completed_before_patch": True,
        "blind_candidate_wgd_completed_before_patch": False,
        "candidate_generation_completed_before_patch": False,
        "formal_scores_generated_before_patch": False,
        "truth_labels_accessed_before_patch": False,
        "automatic_approval": False,
        "changed_files": [
            {
                "status": "M",
                "relative_path": "src/ploidypatch/gff_compat.py",
                "patch_sha256": "c" * 64,
            }
        ],
        "failed_attempt": {"exit_status": 17, "tree_sha256": "d" * 64},
    }
    execution_manifest.update(
        {
            "code_commit": patch_commit,
            "freeze_stage": patch["freeze_stage"],
            "created_before": {
                "wgd_pair_enumeration": False,
                "candidate_generation": True,
                "candidate_labels": True,
                "candidate_scores": True,
            },
            "execution_patch": patch,
        }
    )
    (execution / "execution_manifest.json").write_text(
        json.dumps(execution_manifest) + "\n", encoding="utf-8"
    )
    (execution / "superseded_failed_attempt_manifest.tsv").write_text(
        "relative_path\tbytes\tsha256\nexit_status.txt\t3\t" + "d" * 64 + "\n",
        encoding="utf-8",
    )
    (execution / "patch_reason.md").write_text("reason\n", encoding="utf-8")
    monkeypatch.setenv("PLOIDYPATCH_CODE_COMMIT", patch_commit)
    patched_report = VERIFIER.verify_context(
        input_root=input_root,
        contract_path=contract_path,
        protocol_root=protocol,
        execution_root=execution,
        blind_benchmark_root=benchmark,
        expected_holdout_id=contract.holdout_id,
        expected_primary_chromosomes=29,
    )
    assert patched_report["execution_patch"]["active"] is True
    assert patched_report["execution_patch"]["base_code_commit"] == "a" * 40
    assert patched_report["execution_patch"]["patch_code_commit"] == patch_commit

    execution_manifest["network_access_in_blind_runner"] = True
    (execution / "execution_manifest.json").write_text(
        json.dumps(execution_manifest) + "\n", encoding="utf-8"
    )
    with pytest.raises(ValueError, match="blind context contract"):
        VERIFIER.verify_context(
            input_root=input_root,
            contract_path=contract_path,
            protocol_root=protocol,
            execution_root=execution,
            blind_benchmark_root=benchmark,
            expected_holdout_id=contract.holdout_id,
            expected_primary_chromosomes=29,
        )
