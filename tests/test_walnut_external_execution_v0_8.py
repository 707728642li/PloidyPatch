from __future__ import annotations

import importlib.util
import io
import json
from pathlib import Path
import subprocess
import sys
import tarfile
from types import SimpleNamespace

import pytest

from ploidypatch.artifact_manifest import sha256_file, write_sha256sums
from ploidypatch.holdout_contract import CORE_H1_MODEL_VERSION
from ploidypatch.walnut_h1_framework import (
    BLIND_OUTPUTS,
    EVALUATION_SCHEMA,
    EVALUATION_STATUS_SCHEMA,
    PIPELINE_ENTRIES,
    POLICY_ID,
    RAW_PREDICTION_TREE_KEYS,
    REQUIRED_ENVIRONMENTS,
    REVEAL_STATUS_SCHEMA,
    load_json,
    validate_status,
    verify_execution,
    verify_protocol,
)


ROOT = Path(__file__).parents[1]


def load(relative: str, name: str):
    path = ROOT / relative
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


FREEZE = load(
    "scripts/freeze_walnut_external_execution_v0.8.py", "walnut_execution_freeze_v0_8"
)
ROLE = load(
    "scripts/build_walnut_blind_role_root_v0.8.py", "walnut_blind_roles_v0_8"
)
CUSTODY = load(
    "scripts/finalize_walnut_blind_custody_v0.8.py", "walnut_blind_custody_v0_8"
)
REVEAL = load(
    "scripts/run_walnut_external_reveal_v0.8.py", "walnut_external_reveal_v0_8"
)


def write_json(path: Path, value: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def make_protocol_and_stage(tmp_path: Path) -> tuple[Path, Path]:
    stage = tmp_path / "stage"
    (stage / "shared_target").mkdir(parents=True)
    (stage / "candidate_only").mkdir()
    (stage / "shared_target/target.fa").write_text(">chr1\nACGT\n", encoding="utf-8")
    (stage / "candidate_only/reference.fa").write_text(">ref\nACGT\n", encoding="utf-8")
    fields = (
        "role", "species_id", "release", "bundle_id", "wgdi_prefix", "artifact",
        "bytes", "sha256", "source_relative_path", "staged_relative_path",
        "staged_sha256",
    )
    lines = ["\t".join(fields)]
    for bundle_number, (bundle, species, prefix) in enumerate(
        (
            ("candidate_mandshurica", "Juglans_mandshurica", "jma"),
            ("candidate_carya", "Carya_illinoinensis", "cil"),
        ),
        start=1,
    ):
        for artifact_number, artifact in enumerate(("genome", "gff3", "protein"), start=1):
            digest = f"{bundle_number * 10 + artifact_number:064x}"
            lines.append(
                "\t".join(
                    (
                        "candidate_reference", species, "synthetic", bundle, prefix,
                        artifact, "1", digest, f"source/{bundle}/{artifact}",
                        f"candidate_only/{bundle}/{artifact}", digest,
                    )
                )
            )
    (stage / "role_manifest.tsv").write_text("\n".join(lines) + "\n", encoding="utf-8")
    write_json(stage / "role_contract.json", {"holdout_id": "walnut_walnut2_v0.8"})
    write_sha256sums(stage)

    protocol = tmp_path / "protocol"
    protocol.mkdir()
    contract = ROOT / "config/holdouts/walnut_walnut2_v0.8/contract.json"
    (protocol / "contract.json").write_bytes(contract.read_bytes())
    (protocol / "role_manifest.tsv").write_bytes((stage / "role_manifest.tsv").read_bytes())
    (protocol / "role_contract.json").write_bytes((stage / "role_contract.json").read_bytes())
    write_json(
        protocol / "protocol_manifest.json",
        {
            "schema_version": "ploidypatch.walnut_core_h1_protocol_freeze.v0.8",
            "holdout_id": "walnut_walnut2_v0.8",
            "policy_id": POLICY_ID,
            "model_version": CORE_H1_MODEL_VERSION,
            "contract_sha256": sha256_file(protocol / "contract.json"),
            "staged_input_SHA256SUMS_sha256": sha256_file(stage / "SHA256SUMS"),
            "code_commit": "0" * 40,
            "formal_runner_frozen": False,
            "ranker_enabled": False,
            "h2_or_topology_ranking_enabled": False,
            "truth_access": False,
            "wgd_pairs_enumerated": False,
            "candidate_counts_computed": False,
            "truth_labels_accessed": False,
            "all_arm_collateral_loss_maximum": 0,
        },
    )
    write_sha256sums(protocol)
    return protocol, stage


def make_execution(tmp_path: Path, protocol: Path) -> Path:
    execution = tmp_path / "execution"
    execution.mkdir()
    environments = []
    for name in sorted(REQUIRED_ENVIRONMENTS):
        prefix = tmp_path / "envs" / name
        (prefix / "conda-meta").mkdir(parents=True)
        (prefix / "conda-meta/history").write_text("synthetic\n", encoding="utf-8")
        environments.append({"name": name, "host_prefix": str(prefix.resolve())})
    write_json(
        execution / "execution_manifest.json",
        {
            "schema_version": "ploidypatch.walnut_core_h1_execution_freeze.v0.8",
            "holdout_id": "walnut_walnut2_v0.8",
            "policy_id": POLICY_ID,
            "code_commit": "1" * 40,
            "protocol_SHA256SUMS_sha256": sha256_file(protocol / "SHA256SUMS"),
            "contract_sha256": sha256_file(protocol / "contract.json"),
            "pipeline_entries": PIPELINE_ENTRIES,
            "blind_outputs": BLIND_OUTPUTS,
            "ranker_or_model_execution": False,
            "h2_or_topology_ranking_enabled": False,
            "network_access_in_blind_runner": False,
            "nas_data_mount_in_blind_runner": False,
            "complete_target_annotation_mount_in_blind_runner": False,
            "evaluator_only_mount_in_blind_runner": False,
            "truth_or_label_mount_in_blind_runner": False,
            "environments": environments,
        },
    )
    write_sha256sums(execution)
    return execution


def test_execution_patch_requires_exact_diff_and_preserves_failed_attempt(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    protocol, _stage = make_protocol_and_stage(tmp_path)
    superseded = make_execution(tmp_path, protocol)
    code_root = tmp_path / "code"
    reason = code_root / "docs/patch.md"
    reason.parent.mkdir(parents=True)
    reason.write_text("source-only eligibility repair\n", encoding="utf-8")
    failed_log = tmp_path / "failed.log"
    failed_log.write_text("Hidden CDS chains remain\n", encoding="utf-8")
    changed = [
        ("M", "src/ploidypatch/copy_pair_sampling.py"),
        ("A", "docs/patch.md"),
    ]
    monkeypatch.setattr(FREEZE, "git_changed_files", lambda *_args: changed)

    audit = FREEZE.validate_execution_patch(
        code_root=code_root,
        code_commit="2" * 40,
        protocol=protocol,
        superseded_execution=superseded,
        failed_attempt_log=failed_log,
        failed_stage="structure_holdout_sentinel",
        patch_reason=reason,
        allowed_changed_files=[relative for _status, relative in changed],
    )
    assert audit["base_commit"] == "1" * 40
    assert audit["reason_relative"] == "docs/patch.md"
    assert audit["changed"] == changed
    assert audit["previous_patch_depth"] == 0

    prior = load_json(superseded / "execution_manifest.json")
    prior["execution_patch"] = {
        "schema_version": "ploidypatch.walnut_execution_patch.v0.8",
        "chain_depth": 1,
        "patch_code_commit": "1" * 40,
        "blind_candidate_generation_completed_before_patch": False,
        "formal_scores_generated_before_patch": False,
        "truth_labels_accessed_before_patch": False,
    }
    (superseded / "SHA256SUMS").unlink()
    write_json(superseded / "execution_manifest.json", prior)
    write_sha256sums(superseded)
    chained = FREEZE.validate_execution_patch(
        code_root=code_root,
        code_commit="2" * 40,
        protocol=protocol,
        superseded_execution=superseded,
        failed_attempt_log=failed_log,
        failed_stage="blind_normalization",
        patch_reason=reason,
        allowed_changed_files=[relative for _status, relative in changed],
    )
    assert chained["previous_patch_depth"] == 1

    pool_failure = FREEZE.validate_execution_patch(
        code_root=code_root,
        code_commit="2" * 40,
        protocol=protocol,
        superseded_execution=superseded,
        failed_attempt_log=failed_log,
        failed_stage="candidate_pools",
        patch_reason=reason,
        allowed_changed_files=[relative for _status, relative in changed],
    )
    assert pool_failure["failed_stage"] == "candidate_pools"

    with pytest.raises(ValueError, match="exact changed-file whitelist"):
        FREEZE.validate_execution_patch(
            code_root=code_root,
            code_commit="2" * 40,
            protocol=protocol,
            superseded_execution=superseded,
            failed_attempt_log=failed_log,
            failed_stage="candidate_methods",
            patch_reason=reason,
            allowed_changed_files=["docs/patch.md"],
        )


def make_failed_candidate_publication(tmp_path: Path) -> Path:
    attempt = tmp_path / "run6.working"
    attempt.mkdir()
    (attempt / "exit_status.txt").write_text("72\n", encoding="utf-8")
    (attempt / "stdout.log").write_text("blind run\n", encoding="utf-8")
    (attempt / "stderr.log").write_text(
        "missing exact Walnut blind output: raw_predictions.manifest.json\n",
        encoding="utf-8",
    )
    (attempt / "bwrap_command.txt").write_text(
        "bwrap --unshare-all --unshare-net --clearenv\n", encoding="utf-8"
    )
    write_json(
        attempt / "mount_manifest.json",
        {
            "mounts": [
                {"role": role, "host_path": f"/safe/{role}", "namespace_path": f"/holdout/{role}"}
                for role in ("shared_target", "candidate_only", "blind_benchmark")
            ]
        },
    )
    write_json(
        attempt / "namespace_role_validation.json",
        {
            "shared_target_visible": True,
            "candidate_only_visible": True,
            "blind_benchmark_visible": True,
            "evaluator_only_visible": False,
            "truth_visible": False,
            "complete_target_annotation_visible": False,
            "nas_data_visible": False,
        },
    )
    project = attempt / "project"
    project.mkdir()
    (project / "pipeline_commands.tsv").write_text("step\tcommand\n", encoding="utf-8")
    wrong = project / "results/copy_collapse/external/walnut/v0.8_h1"
    wrong.mkdir(parents=True)
    canonical_prefix = "results/copy_collapse/external/walnut_v0.8_h1/"
    for name, relative in BLIND_OUTPUTS.items():
        if name == "command_log":
            continue
        assert relative.startswith(canonical_prefix)
        path = wrong / relative.removeprefix(canonical_prefix)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(f"{name}\n", encoding="utf-8")
    write_sha256sums(wrong)
    return attempt


def test_post_pool_publication_patch_binds_full_failure_and_forbids_formal_output(
    tmp_path: Path,
) -> None:
    attempt = make_failed_candidate_publication(tmp_path)
    audit = FREEZE.audit_failed_candidate_publication(attempt, BLIND_OUTPUTS)
    assert audit["exit_status"] == 72
    assert audit["canonical_outputs_present"] is False
    assert audit["noncanonical_pool_checksum_verified"] is True
    assert audit["truth_labels_accessed"] is False
    assert audit["reuse_permitted"] is False
    assert audit["files"] == len(audit["rows"])

    formal = attempt / "project" / BLIND_OUTPUTS["retain_pool"]
    formal.parent.mkdir(parents=True, exist_ok=True)
    formal.write_text("forbidden\n", encoding="utf-8")
    with pytest.raises(ValueError, match="already contains canonical output"):
        FREEZE.audit_failed_candidate_publication(attempt, BLIND_OUTPUTS)


def make_benchmark(tmp_path: Path) -> Path:
    benchmark = tmp_path / "blind_benchmark"
    benchmark.mkdir()
    (benchmark / "perturbed.gff3").write_text("##gff-version 3\n", encoding="utf-8")
    write_json(
        benchmark / "blind_manifest.json",
        {
            "schema_version": "ploidypatch.walnut_h1_blind_benchmark_input.v0.8",
            "truth_access": False,
            "complete_target_annotation_access": False,
            "ranker_access": False,
            "h2_or_topology_ranking_access": False,
            "perturbed_annotation": {
                "file_name": "perturbed.gff3",
                "sha256": sha256_file(benchmark / "perturbed.gff3"),
            },
            "target_genome": {"mount_role": "shared_target_genome", "sha256": "a" * 64},
        },
    )
    write_sha256sums(benchmark)
    return benchmark


def test_protocol_execution_contract_is_exact_and_no_ranker(tmp_path: Path) -> None:
    protocol, _ = make_protocol_and_stage(tmp_path)
    execution = make_execution(tmp_path, protocol)
    protocol_manifest, contract = verify_protocol(protocol)
    execution_manifest, _, _ = verify_execution(execution, protocol)
    assert contract.model_version == CORE_H1_MODEL_VERSION
    assert protocol_manifest["ranker_enabled"] is False
    assert execution_manifest["ranker_or_model_execution"] is False
    assert set(BLIND_OUTPUTS) == {
        "raw_predictions_manifest", "retain_pool", "retain_decisions",
        "retain_manifest", "suppress_pool", "suppress_decisions",
        "suppress_manifest", "command_log",
    }
    assert all("score" not in value.casefold() for value in BLIND_OUTPUTS.values())


def test_pipeline_file_validation_rejects_ranker_h2_and_accepts_exact_closure(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    files = (
        set(FREEZE.ENGINEERING_FILES)
        | set(FREEZE.SCIENTIFIC_FILES)
        | set(PIPELINE_ENTRIES.values())
    )
    for relative in files:
        path = tmp_path / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("safe\n", encoding="utf-8")
    monkeypatch.setattr(FREEZE, "ENGINEERING_FILES", frozenset())
    monkeypatch.setattr(FREEZE, "SCIENTIFIC_FILES", frozenset())
    assert set(FREEZE.validate_pipeline_files(tmp_path, sorted(files))) == files
    bad = tmp_path / "scripts/run_h2_ranker.py"
    bad.write_text("unsafe\n", encoding="utf-8")
    with pytest.raises(ValueError, match="Ranker/H2"):
        FREEZE.validate_pipeline_files(tmp_path, [*files, "scripts/run_h2_ranker.py"])


def test_source_archive_is_regular_only_and_traversal_safe(tmp_path: Path) -> None:
    good = tmp_path / "good.tar"
    with tarfile.open(good, "w") as archive:
        payload = b"safe\n"
        member = tarfile.TarInfo("scripts/entry.py")
        member.size = len(payload)
        archive.addfile(member, io.BytesIO(payload))
    destination = tmp_path / "good"
    destination.mkdir()
    FREEZE.safe_extract_archive(good, destination)
    assert (destination / "scripts/entry.py").read_bytes() == b"safe\n"

    bad = tmp_path / "bad.tar"
    with tarfile.open(bad, "w") as archive:
        member = tarfile.TarInfo("../escape")
        member.size = 1
        archive.addfile(member, io.BytesIO(b"x"))
    with pytest.raises(ValueError, match="Unsafe source archive member"):
        FREEZE.safe_extract_archive(bad, tmp_path / "bad")


def test_tool_only_conda_environment_has_explicit_lock_without_fake_pip_lock(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    prefix = tmp_path / "tool_only"
    (prefix / "conda-meta").mkdir(parents=True)
    (prefix / "conda-meta/history").write_text("created\n", encoding="utf-8")

    def fake_run(command: list[str], **_kwargs: object) -> SimpleNamespace:
        assert command[:2] == ["conda", "list"]
        return SimpleNamespace(stdout=b"@EXPLICIT\nhttps://example.invalid/tool.tar.bz2\n")

    monkeypatch.setattr(FREEZE.subprocess, "run", fake_run)
    explicit, pip = FREEZE.environment_locks(prefix)
    assert explicit.startswith(b"@EXPLICIT\n")
    assert pip is None


def test_role_root_is_exact_three_roles_atomic_and_nonoverwriting(tmp_path: Path) -> None:
    protocol, stage = make_protocol_and_stage(tmp_path)
    benchmark = make_benchmark(tmp_path)
    output = tmp_path / "sealed_roles"
    ROLE.build_role_root(
        staged_inputs=stage, blind_benchmark=benchmark, protocol=protocol, output=output
    )
    manifest = load_json(output / "role_manifest.json")
    assert manifest["roles"] == ["shared_target", "candidate_only", "blind_benchmark"]
    assert manifest["ranker_or_model_present"] is False
    assert {path.name for path in output.iterdir()} == {
        "shared_target", "candidate_only", "blind_benchmark", "role_manifest.json",
        "SHA256SUMS",
    }
    execution = make_execution(tmp_path, protocol)
    raw_hashes = CUSTODY.expected_raw_input_hashes(
        role_root=output, protocol=protocol, execution=execution
    )
    assert set(raw_hashes) == {
        "staged_input_SHA256SUMS_sha256", "blind_benchmark_SHA256SUMS_sha256",
        "protocol_SHA256SUMS_sha256", "execution_SHA256SUMS_sha256",
        "target_genome_sha256", "perturbed_gff3_sha256",
        "candidate_mandshurica_genome_sha256", "candidate_mandshurica_gff3_sha256",
        "candidate_mandshurica_protein_sha256", "candidate_carya_genome_sha256",
        "candidate_carya_gff3_sha256", "candidate_carya_protein_sha256",
    }
    with pytest.raises(FileExistsError, match="overwrite"):
        ROLE.build_role_root(
            staged_inputs=stage, blind_benchmark=benchmark, protocol=protocol, output=output
        )


def test_role_root_rejects_evaluator_named_input(tmp_path: Path) -> None:
    protocol, stage = make_protocol_and_stage(tmp_path)
    benchmark = make_benchmark(tmp_path)
    (stage / "SHA256SUMS").unlink()
    forbidden = stage / "candidate_only/evaluator_only"
    forbidden.mkdir()
    (forbidden / "reference.gff3").write_text("unsafe\n", encoding="utf-8")
    write_sha256sums(stage)
    # Keep protocol lineage valid so the semantic firewall is the rejecting gate.
    protocol_manifest = load_json(protocol / "protocol_manifest.json")
    protocol_manifest["staged_input_SHA256SUMS_sha256"] = sha256_file(stage / "SHA256SUMS")
    (protocol / "SHA256SUMS").unlink()
    write_json(protocol / "protocol_manifest.json", protocol_manifest)
    write_sha256sums(protocol)
    with pytest.raises(ValueError, match="Forbidden evaluator/truth path"):
        ROLE.build_role_root(
            staged_inputs=stage, blind_benchmark=benchmark, protocol=protocol,
            output=tmp_path / "roles",
        )


def test_raw_prediction_tree_hash_is_permutation_stable_and_mutation_bound(
    tmp_path: Path,
) -> None:
    project = tmp_path / "project"
    trees: dict[str, dict[str, object]] = {}
    for number, name in enumerate(sorted(RAW_PREDICTION_TREE_KEYS)):
        root = project / "raw" / name
        root.mkdir(parents=True)
        (root / "z.gff3").write_text(f"z{number}\n", encoding="utf-8")
        (root / "a.gff3").write_text(f"a{number}\n", encoding="utf-8")
        digest, count, byte_count = CUSTODY.tree_digest(root)
        trees[name] = {
            "relative_path": root.relative_to(project).as_posix(),
            "file_count": count,
            "bytes": byte_count,
            "tree_sha256": digest,
        }
    raw = project / BLIND_OUTPUTS["raw_predictions_manifest"]
    write_json(
        raw,
        {
            "schema_version": CUSTODY.RAW_SCHEMA,
            "holdout_id": "walnut_walnut2_v0.8",
            "policy_id": POLICY_ID,
            "truth_access": False,
            "ranker_access": False,
            "method_families": ["miniprot", "gemoma", "lifton"],
            "candidate_references": ["candidate_mandshurica", "candidate_carya"],
            "within_method_reference_vote_count": 1,
            "input_hashes": {"synthetic": "a" * 64},
            "tree_hash_algorithm":
                "sha256_of_sorted_sha256_two_space_relative_path_newline",
            "raw_prediction_trees": dict(reversed(list(trees.items()))),
        },
    )
    CUSTODY.validate_raw_manifest(
        raw, project, expected_input_hashes={"synthetic": "a" * 64}
    )
    first = sorted(RAW_PREDICTION_TREE_KEYS)[0]
    (project / "raw" / first / "a.gff3").write_text("mutated\n", encoding="utf-8")
    with pytest.raises(ValueError, match="tree differs"):
        CUSTODY.validate_raw_manifest(
            raw, project, expected_input_hashes={"synthetic": "a" * 64}
        )


def make_evaluation(root: Path, *, extra: dict[str, object] | None = None) -> None:
    root.mkdir()
    write_json(
        root / "status.json",
        {"schema_version": EVALUATION_STATUS_SCHEMA, "status": "ready", "reason_codes": []},
    )
    value: dict[str, object] = {
        "schema_version": EVALUATION_SCHEMA,
        "status": "ready",
        "reason_codes": [],
        "ranker_or_model_executed": False,
        "h2_or_topology_ranking_executed": False,
        "bootstrap_replicates": 20_000,
        "bootstrap_unit": "paired_event",
        "metric": "event_exact_phased_CDS_recall_retain_distinct_minus_suppress_overlap",
        "all_arm_collateral_loss": 0,
    }
    value.update(extra or {})
    write_json(root / "evaluation.json", value)
    write_sha256sums(root)


def test_h1_evaluation_accepts_fixed_estimand_and_rejects_ranker_semantics(
    tmp_path: Path,
) -> None:
    good = tmp_path / "good"
    make_evaluation(good)
    assert REVEAL.validate_evaluation(good)["status"] == "ready"
    bad = tmp_path / "bad"
    make_evaluation(bad, extra={"average_precision": 0.99})
    with pytest.raises(ValueError, match="forbidden ranker/H2"):
        REVEAL.validate_evaluation(bad)


def test_tri_state_and_not_evaluable_reason_class_are_fail_closed() -> None:
    assert validate_status(
        {"schema_version": REVEAL_STATUS_SCHEMA, "status": "ready", "reason_codes": []},
        expected_schema=REVEAL_STATUS_SCHEMA,
    ) == "ready"
    REVEAL.validate_reason_class(
        "not_evaluable", ["formal_event_count_below_500"]
    )
    with pytest.raises(ValueError, match="non-data-gate"):
        REVEAL.validate_reason_class("not_evaluable", ["adapter_failed"])
    with pytest.raises(ValueError, match="Malformed tri-state"):
        validate_status(
            {"schema_version": REVEAL_STATUS_SCHEMA, "status": "invalid", "reason_codes": []},
            expected_schema=REVEAL_STATUS_SCHEMA,
        )


def test_reveal_not_evaluable_is_sealed_without_invoking_evaluator(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    protocol = tmp_path / "protocol"
    execution = tmp_path / "execution"
    blind = tmp_path / "blind"
    evaluator_only = tmp_path / "evaluator_only"
    for root in (protocol, execution, blind, evaluator_only):
        root.mkdir()
    (protocol / "contract.json").write_text("{}\n", encoding="utf-8")
    (protocol / "protocol.txt").write_text("frozen\n", encoding="utf-8")
    write_sha256sums(protocol)
    (execution / "execution.txt").write_text("frozen\n", encoding="utf-8")
    write_sha256sums(execution)
    write_json(blind / "custody_manifest.json", {"sealed": True})
    write_sha256sums(blind)
    (evaluator_only / "complete.gff3").write_text("##gff-version 3\n", encoding="utf-8")
    write_sha256sums(evaluator_only)

    fake_execution = {
        "environments": [{"name": "ploidypatch-dev", "host_prefix": str(tmp_path / "dev")}]
    }
    monkeypatch.setattr(
        REVEAL,
        "verify_execution",
        lambda _execution, _protocol: (
            fake_execution,
            {},
            SimpleNamespace(holdout_id="walnut_walnut2_v0.8"),
        ),
    )
    monkeypatch.setattr(REVEAL, "validate_custody", lambda *_args: {"sealed": True})
    calls: list[str] = []

    def fake_run(**kwargs) -> int:
        calls.append(Path(kwargs["entry"]).name)
        reveal_root = Path(kwargs["environment"]["PLOIDYPATCH_REVEAL_INPUTS_OUTPUT"])
        reveal_root.mkdir()
        write_json(
            reveal_root / "status.json",
            {
                "schema_version": REVEAL_STATUS_SCHEMA,
                "status": "not_evaluable",
                "reason_codes": ["formal_event_count_below_500"],
            },
        )
        write_sha256sums(reveal_root)
        Path(kwargs["stdout"]).write_text("builder\n", encoding="utf-8")
        Path(kwargs["stderr"]).write_text("", encoding="utf-8")
        return 0

    monkeypatch.setattr(REVEAL, "run_frozen_entry", fake_run)
    output = tmp_path / "reveal"
    REVEAL.orchestrate_reveal(
        project_root=tmp_path / "project",
        protocol=protocol,
        execution=execution,
        blind_run=blind,
        evaluator_only=evaluator_only,
        output=output,
    )
    assert calls == [PIPELINE_ENTRIES["reveal_input_builder"].split("/")[-1]]
    status = load_json(output / "status.json")
    assert status["status"] == "not_evaluable"
    assert status["evaluator_invoked"] is False
    assert not (output / "evaluation").exists()
    authorization = load_json(output / "reveal_authorization.json")
    assert authorization["truth_reveal_authorized"] is True
    assert authorization["ranker_or_model_authorized"] is False


def test_blind_launcher_has_network_namespace_firewall_and_no_model_mount() -> None:
    launcher = ROOT / "scripts/run_walnut_blind_isolated_v0.8.sh"
    text = launcher.read_text(encoding="utf-8")
    for token in (
        "--unshare-all", "--unshare-net", "--clearenv",
        "PLOIDYPATCH_NETWORK_ACCESS none", "/holdout/shared_target",
        "/holdout/candidate_only", "/holdout/blind_benchmark",
    ):
        assert token in text
    assert "--ro-bind $model" not in text
    assert "--setenv PLOIDYPATCH_COMPOSITE_MODEL_FREEZE" not in text
    if sys.platform != "win32":
        subprocess.run(["bash", "-n", str(launcher)], check=True)
