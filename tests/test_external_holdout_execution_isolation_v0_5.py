from __future__ import annotations

import csv
import importlib.util
import io
import json
import os
from pathlib import Path
from types import SimpleNamespace
import subprocess
import sys
import tarfile

import pytest

from ploidypatch.artifact_manifest import (
    sha256_file,
    verify_sha256sums,
    write_sha256sums,
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


PROTOCOL = load(
    "scripts/freeze_external_holdout_protocol_v0.5.py", "holdout_protocol_v0_5"
)
EXECUTION = load(
    "scripts/freeze_external_holdout_execution_v0.5.py", "holdout_execution_v0_5"
)
ROLE_ROOT = load(
    "scripts/build_external_holdout_blind_role_root_v0.5.py", "blind_role_root_v0_5"
)
CUSTODY = load(
    "scripts/finalize_external_holdout_blind_custody_v0.5.py", "blind_custody_v0_5"
)


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def test_actinidia_too_late_paths_pipeline_closure_and_blind_outputs_are_exact() -> None:
    contract = SimpleNamespace(holdout_id="actinidia_red5_v0.5")
    forbidden = set(PROTOCOL.default_forbidden_paths(contract))
    assert "data/derived/external_inputs/actinidia/v0.5" in forbidden
    assert "data/derived/holdout_inputs/actinidia_v0.5" not in forbidden
    assert {
        "data/derived/external_evaluator/actinidia_v0.5_wgdi_inputs",
        "results/evaluator/actinidia/v0.5",
        "benchmark/structure/copy_collapse_v0.5/red5_ps1_1.69.0/annotation_copy_collapse_seed20261010",
        "results/baselines/actinidia_v0.5",
        "results/copy_collapse/external/actinidia_v0.5_method_trio",
        "results/copy_collapse/external/actinidia_v0.5_blind_self_wgd",
        "results/copy_collapse/external/actinidia_v0.5_blind_rankings",
    } <= forbidden

    required = PROTOCOL.HOLDOUT_REQUIRED_PIPELINE_FILES[contract.holdout_id]
    assert {
        "scripts/prepare_actinidia_external_normalized_inputs_v0.5.sh",
        "scripts/prepare_actinidia_evaluator_wgdi_inputs_v0.5.sh",
        "scripts/run_actinidia_evaluator_wgdi_v0.5.sh",
        "scripts/infer_actinidia_external_pairs_v0.5.sh",
        "scripts/run_actinidia_copy_collapse_benchmark_v0.5.sh",
        "scripts/build_actinidia_complete_control_reveal_inputs_v0.5.sh",
        "scripts/evaluate_actinidia_external_v0.5.py",
        "scripts/evaluate_external_v0.5.py",
    } <= required
    assert "scripts/publish_gemoma_working.sh" in required
    assert PROTOCOL.HOLDOUT_REQUIRED_PIPELINE_ENTRIES[contract.holdout_id] == {
        "blind_pipeline": "scripts/run_actinidia_blind_pipeline_v0.5.sh",
        "reveal_input_builder": "scripts/build_actinidia_complete_control_reveal_inputs_v0.5.sh",
        "evaluator": "scripts/evaluate_actinidia_external_v0.5.py",
    }
    assert PROTOCOL.HOLDOUT_REQUIRED_BLIND_OUTPUTS[contract.holdout_id] == {
        "scores": "results/copy_collapse/external/actinidia_v0.5_blind_rankings/scores/v04.tsv",
        "score_manifest": "results/copy_collapse/external/actinidia_v0.5_blind_rankings/scores/v04.tsv.manifest.json",
        "pool_decisions": "results/copy_collapse/external/actinidia_v0.5_method_trio/consensus/primary_union/blind/decisions.tsv",
        "pool_manifest": "results/copy_collapse/external/actinidia_v0.5_method_trio/consensus/primary_union/blind/candidate.gff3.manifest.json",
        "command_log": "pipeline_commands.tsv",
    }


def test_source_archive_rejects_traversal_and_links(tmp_path: Path) -> None:
    for name, make_member in (
        ("traversal.tar", lambda: tarfile.TarInfo("../escape")),
        ("link.tar", lambda: tarfile.TarInfo("unsafe-link")),
    ):
        archive = tmp_path / name
        with tarfile.open(archive, "w") as handle:
            member = make_member()
            if name == "link.tar":
                member.type = tarfile.SYMTYPE
                member.linkname = "/etc/passwd"
                handle.addfile(member)
            else:
                payload = b"bad"
                member.size = len(payload)
                handle.addfile(member, io.BytesIO(payload))
        with pytest.raises(ValueError, match="Unsafe member"):
            EXECUTION.extract_archive_safely(archive, tmp_path / f"extract-{name}")


def test_strict_manifest_rejects_dot_paths_and_accepts_canonical_paths(
    tmp_path: Path,
) -> None:
    artifact = tmp_path / "nested/output.tsv"
    artifact.parent.mkdir()
    artifact.write_text("field\tvalue\n", encoding="utf-8")
    checksum = tmp_path / "SHA256SUMS"
    checksum.write_text(
        f"{sha256_file(artifact)}  ./nested/output.tsv\n", encoding="utf-8"
    )
    with pytest.raises(ValueError, match="Unsafe manifest path"):
        verify_sha256sums(tmp_path, ignore_checksum_file=True)
    checksum.unlink()
    write_sha256sums(tmp_path)
    assert verify_sha256sums(tmp_path, ignore_checksum_file=True) == {
        "nested/output.tsv": sha256_file(artifact)
    }


def test_execution_patch_requires_exact_changed_file_whitelist_and_reason() -> None:
    changed = [
        ("M", "src/ploidypatch/gff_compat.py"),
        ("A", "docs/ACTINIDIA_EXECUTION_PATCH_attempt1.md"),
        ("M", "tests/test_external_holdout_execution_isolation_v0_5.py"),
    ]
    allowed = [path for _, path in changed]
    EXECUTION.verify_changed_whitelist(
        changed, allowed, "docs/ACTINIDIA_EXECUTION_PATCH_attempt1.md"
    )
    with pytest.raises(ValueError, match="exactly equal"):
        EXECUTION.verify_changed_whitelist(
            changed,
            allowed[:-1],
            "docs/ACTINIDIA_EXECUTION_PATCH_attempt1.md",
        )
    with pytest.raises(ValueError, match="Non-implementation"):
        EXECUTION.verify_changed_whitelist(
            [("M", "config/scientific_policy.tsv"), ("A", allowed[1])],
            ["config/scientific_policy.tsv", allowed[1]],
            allowed[1],
        )


def make_failed_attempt(tmp_path: Path) -> Path:
    attempt = tmp_path / "attempt1.working"
    attempt.mkdir()
    (attempt / "exit_status.txt").write_text("17\n", encoding="utf-8")
    (attempt / "stdout.log").write_text("candidate-safe startup\n", encoding="utf-8")
    (attempt / "stderr.log").write_text("adapter failed early\n", encoding="utf-8")
    (attempt / "bwrap_command.txt").write_text(
        "bwrap --unshare-all --unshare-net --clearenv\n", encoding="utf-8"
    )
    (attempt / "mount_manifest.json").write_text(
        json.dumps(
            {
                "schema_version": "ploidypatch.blind_mount_manifest.v0.5",
                "mounts": [
                    {
                        "role": role,
                        "host_path": f"/sealed/{role}",
                        "namespace_path": f"/holdout/{role}",
                        "read_only": True,
                    }
                    for role in ("shared_target", "candidate_only", "blind_benchmark")
                ],
            }
        ),
        encoding="utf-8",
    )
    (attempt / "namespace_role_validation.json").write_text(
        json.dumps(
            {
                "schema_version": "ploidypatch.blind_namespace_validation.v0.5",
                "shared_target_visible": True,
                "candidate_only_visible": True,
                "blind_benchmark_visible": True,
                "evaluator_only_visible": False,
                "truth_visible": False,
                "complete_target_annotation_visible": False,
                "nas_data_visible": False,
            }
        ),
        encoding="utf-8",
    )
    (attempt / "project/results/baselines/raw").mkdir(parents=True)
    (attempt / "project/results/baselines/raw/prediction.gff3").write_text(
        "##gff-version 3\n", encoding="utf-8"
    )
    return attempt


def test_failed_attempt_patch_manifest_is_pre_candidate_pre_score_pre_label(
    tmp_path: Path,
) -> None:
    attempt = make_failed_attempt(tmp_path)
    output = tmp_path / "failed.tsv"
    summary = EXECUTION.freeze_failed_attempt_manifest(
        attempt,
        output,
        {
            "scores": "results/ranking/scores/v04.tsv",
            "score_manifest": "results/ranking/scores/v04.tsv.manifest.json",
            "pool_decisions": "results/pool/decisions.tsv",
            "pool_manifest": "results/pool/candidate.gff3.manifest.json",
            "command_log": "pipeline_commands.tsv",
        },
    )
    assert summary["exit_status"] == 17
    assert len(summary["tree_sha256"]) == 64
    (attempt / "project/results/pool").mkdir(parents=True)
    (attempt / "project/results/pool/candidate.gff3").write_text(
        "##gff-version 3\n", encoding="utf-8"
    )
    with pytest.raises(ValueError, match="pre-candidate"):
        EXECUTION.freeze_failed_attempt_manifest(
            attempt,
            tmp_path / "rejected.tsv",
            {"scores": "results/ranking/scores/v04.tsv"},
        )


def test_custody_accepts_only_strict_execution_patch_provenance(tmp_path: Path) -> None:
    execution = tmp_path / "execution"
    protocol = tmp_path / "protocol"
    model = tmp_path / "model"
    for root in (execution, protocol, model):
        root.mkdir()
    (protocol / "contract.json").write_text("{}\n", encoding="utf-8")
    (protocol / "SHA256SUMS").write_text("base protocol\n", encoding="utf-8")
    (model / "SHA256SUMS").write_text("model\n", encoding="utf-8")
    failed_row = {
        "relative_path": "exit_status.txt",
        "bytes": 3,
        "sha256": "d" * 64,
    }
    write_tsv(execution / "superseded_failed_attempt_manifest.tsv", [failed_row])
    (execution / "patch_reason.md").write_text("reason\n", encoding="utf-8")
    protocol_manifest = {"code_commit": "a" * 40}
    patch = {
        "schema_version": "ploidypatch.external_holdout_execution_patch.v0.5",
        "freeze_stage": EXECUTION.PATCH_STAGE,
        "base_code_commit": "a" * 40,
        "patch_code_commit": "b" * 40,
        "base_protocol_SHA256SUMS_sha256": sha256_file(protocol / "SHA256SUMS"),
        "superseded_execution_SHA256SUMS_sha256": "e" * 64,
        "contract_sha256": sha256_file(protocol / "contract.json"),
        "composite_model_SHA256SUMS_sha256": sha256_file(model / "SHA256SUMS"),
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
        "failed_attempt": {
            "exit_status": 17,
            "tree_sha256": EXECUTION._tree_digest([failed_row]),
            "files": 1,
            "bytes": 3,
        },
    }
    manifest = {
        "code_commit": "b" * 40,
        "freeze_stage": EXECUTION.PATCH_STAGE,
        "created_before": {
            "wgd_pair_enumeration": False,
            "candidate_generation": True,
            "candidate_labels": True,
            "candidate_scores": True,
        },
        "execution_patch": patch,
    }
    CUSTODY.validate_execution_patch(
        manifest, protocol_manifest, execution, protocol, model
    )
    patch["model_or_threshold_changed"] = True
    with pytest.raises(ValueError, match="firewall"):
        CUSTODY.validate_execution_patch(
            manifest, protocol_manifest, execution, protocol, model
        )


def make_role_fixture(tmp_path: Path) -> tuple[Path, Path, Path, SimpleNamespace]:
    stage = tmp_path / "stage"
    (stage / "shared_target/target_red5").mkdir(parents=True)
    (stage / "candidate_only/candidate_eriantha").mkdir(parents=True)
    target = stage / "shared_target/target_red5/genome.fa.gz"
    candidate = stage / "candidate_only/candidate_eriantha/genome.fa.gz"
    target.write_bytes(b"sealed target genome\n")
    candidate.write_bytes(b"sealed candidate genome\n")
    rows = [
        {
            "role": "target",
            "species_id": "Actinidia_chinensis_Red5",
            "release": "Red5",
            "bundle_id": "target_red5",
            "wgdi_prefix": "red5",
            "artifact": "genome",
            "bytes": target.stat().st_size,
            "sha256": sha256_file(target),
            "source_relative_path": "raw/target/genome.fa.gz",
            "staged_relative_path": "shared_target/target_red5/genome.fa.gz",
            "staged_sha256": sha256_file(target),
        },
        {
            "role": "candidate_reference",
            "species_id": "Actinidia_eriantha",
            "release": "White",
            "bundle_id": "candidate_eriantha",
            "wgdi_prefix": "aer",
            "artifact": "genome",
            "bytes": candidate.stat().st_size,
            "sha256": sha256_file(candidate),
            "source_relative_path": "raw/candidate/genome.fa.gz",
            "staged_relative_path": "candidate_only/candidate_eriantha/genome.fa.gz",
            "staged_sha256": sha256_file(candidate),
        },
    ]
    write_tsv(stage / "role_manifest.tsv", rows)
    (stage / "role_contract.json").write_text("{}\n", encoding="utf-8")
    write_sha256sums(stage)

    benchmark = tmp_path / "blind"
    benchmark.mkdir()
    perturbed = benchmark / "perturbed.gff3"
    perturbed.write_text("##gff-version 3\n", encoding="utf-8")
    (benchmark / "blind_manifest.json").write_text(
        json.dumps(
            {
                "schema_version": "ploidypatch.blind_benchmark_input.v0.5",
                "truth_access": False,
                "complete_target_annotation_access": False,
                "perturbed_annotation": {
                    "file_name": "perturbed.gff3",
                    "sha256": sha256_file(perturbed),
                },
                "target_genome": {
                    "mount_role": "shared_target_genome",
                    "sha256": "a" * 64,
                },
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    write_sha256sums(benchmark)

    contract = SimpleNamespace(
        holdout_id="actinidia_red5_v0.5",
        policy_id="ploidypatch_actinidia_external_validation_v0.5",
        model_version="PloidyPatch_ranker_v0.4",
        target=SimpleNamespace(species_id="Actinidia_chinensis_Red5"),
    )
    protocol = tmp_path / "protocol"
    protocol.mkdir()
    (protocol / "contract.json").write_text("{}\n", encoding="utf-8")
    (protocol / "role_manifest.tsv").write_bytes((stage / "role_manifest.tsv").read_bytes())
    (protocol / "role_contract.json").write_bytes((stage / "role_contract.json").read_bytes())
    (protocol / "protocol_manifest.json").write_text(
        json.dumps(
            {
                "schema_version": ROLE_ROOT.PROTOCOL_SCHEMA,
                "holdout_id": contract.holdout_id,
                "code_commit": "a" * 40,
                "contract_sha256": sha256_file(protocol / "contract.json"),
                "staged_input_SHA256SUMS_sha256": sha256_file(stage / "SHA256SUMS"),
            }
        )
        + "\n",
        encoding="utf-8",
    )
    write_sha256sums(protocol)
    return stage, benchmark, protocol, contract


def test_role_root_builder_is_atomic_truth_free_and_exact(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    stage, benchmark, protocol, contract = make_role_fixture(tmp_path)
    monkeypatch.setattr(ROLE_ROOT, "load_holdout_contract", lambda _: contract)
    monkeypatch.setattr(ROLE_ROOT, "chmod_read_only", lambda _: None)
    output = tmp_path / "blind-role"
    ROLE_ROOT.build_role_root(
        staged_inputs=stage,
        blind_benchmark=benchmark,
        protocol=protocol,
        output=output,
    )
    verify_sha256sums(output, ignore_checksum_file=True)
    manifest = json.loads((output / "role_manifest.json").read_text(encoding="utf-8"))
    assert manifest["roles"] == [
        "shared_target",
        "candidate_only",
        "blind_benchmark",
    ]
    assert manifest["truth_access"] is False
    assert manifest["complete_target_annotation_present"] is False
    assert manifest["evaluator_references_present"] is False
    assert not (output / "evaluator_only").exists()
    assert not Path(str(output) + ".partial").exists()
    with pytest.raises(FileExistsError, match="overwrite"):
        ROLE_ROOT.build_role_root(
            staged_inputs=stage,
            blind_benchmark=benchmark,
            protocol=protocol,
            output=output,
        )


@pytest.mark.skipif(os.name == "nt", reason="bubblewrap mount namespaces are Linux-only")
def test_mount_contract_accepts_exact_three_roles_and_rejects_evaluator(
    tmp_path: Path,
) -> None:
    role_root = tmp_path / "roles"
    execution = tmp_path / "execution"
    protocol = tmp_path / "protocol"
    model = tmp_path / "model"
    output = tmp_path / "blind-output"
    for path in (role_root, execution, protocol, model, output, execution / "source"):
        path.mkdir(parents=True, exist_ok=True)
    for role in ("shared_target", "candidate_only", "blind_benchmark"):
        (role_root / role).mkdir()
    environments = {}
    for name in sorted(EXECUTION.REQUIRED_ENVIRONMENTS):
        prefix = tmp_path / "envs" / name
        prefix.mkdir(parents=True)
        environments[name] = prefix
    (execution / "execution_manifest.json").write_text(
        json.dumps(
            {
                "environments": [
                    {"name": name, "host_prefix": str(prefix)}
                    for name, prefix in sorted(environments.items())
                ]
            }
        ),
        encoding="utf-8",
    )
    (role_root / "role_manifest.json").write_text("{}\n", encoding="utf-8")
    (protocol / "role_manifest.tsv").write_text("role\n", encoding="utf-8")
    (protocol / "role_contract.json").write_text("{}\n", encoding="utf-8")
    mounts = [
        {
            "role": role,
            "host_path": str(role_root / role),
            "namespace_path": f"/holdout/{role}",
            "read_only": True,
        }
        for role in ("shared_target", "candidate_only", "blind_benchmark")
    ] + [
        {
            "role": "frozen_execution",
            "host_path": str(execution),
            "namespace_path": "/frozen/execution",
            "read_only": True,
        },
        {
            "role": "frozen_source",
            "host_path": str(execution / "source"),
            "namespace_path": "/frozen/source",
            "read_only": True,
        },
        {
            "role": "frozen_source",
            "host_path": str(execution / "source"),
            "namespace_path": "/run/blind-run/project/code",
            "read_only": True,
        },
        {
            "role": "frozen_protocol",
            "host_path": str(protocol),
            "namespace_path": "/frozen/protocol",
            "read_only": True,
        },
        {
            "role": "frozen_model",
            "host_path": str(model),
            "namespace_path": "/frozen/model",
            "read_only": True,
        },
        {
            "role": "blind_output",
            "host_path": str(output),
            "namespace_path": "/run/blind-run",
            "read_only": False,
        },
    ] + [
        {
            "role": f"frozen_environment:{name}",
            "host_path": str(prefix),
            "namespace_path": namespace,
            "read_only": True,
        }
        for name, prefix in sorted(environments.items())
        for namespace in (
            str(prefix.resolve()),
            f"/frozen/envs/{name}",
            f"/run/blind-run/project/envs/{name}",
        )
    ] + [
        {
            "role": "system_role_metadata",
            "host_path": str(host),
            "namespace_path": namespace,
            "read_only": True,
        }
        for host, namespace in (
            (role_root / "role_manifest.json", "/holdout/blind_role_manifest.json"),
            (protocol / "role_manifest.tsv", "/holdout/role_manifest.tsv"),
            (protocol / "role_contract.json", "/holdout/role_contract.json"),
        )
    ]
    manifest_path = tmp_path / "mounts.json"
    payload = {
        "schema_version": CUSTODY.MOUNT_SCHEMA,
        "holdout_id": "actinidia_red5_v0.5",
        "mounts": mounts,
    }
    manifest_path.write_text(json.dumps(payload), encoding="utf-8")
    CUSTODY.validate_mount_manifest(
        manifest_path,
        payload["holdout_id"],
        role_root,
        execution,
        protocol,
        model,
        output,
    )
    payload["mounts"].append(
        {
            "role": "evaluator_only",
            "host_path": str(tmp_path / "evaluator_only"),
            "namespace_path": "/holdout/evaluator_only",
            "read_only": True,
        }
    )
    manifest_path.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(ValueError, match="Unsafe blind mount|Unexpected blind mount"):
        CUSTODY.validate_mount_manifest(
            manifest_path,
            payload["holdout_id"],
            role_root,
            execution,
            protocol,
            model,
            output,
        )


def test_custody_candidate_tables_reject_truth_or_label_columns(tmp_path: Path) -> None:
    clean = tmp_path / "clean.tsv"
    write_tsv(clean, [{"candidate_digest": "a", "score": "0.2"}])
    assert CUSTODY.score_digests(clean) == {"a"}
    leaked = tmp_path / "leaked.tsv"
    write_tsv(leaked, [{"candidate_digest": "a", "truth_label": "1"}])
    with pytest.raises(ValueError, match="truth/label"):
        CUSTODY.score_digests(leaked)


def test_bwrap_and_reveal_scripts_encode_fail_closed_isolation_and_barrier() -> None:
    blind = (ROOT / "scripts/run_external_holdout_blind_isolated_v0.5.sh").read_text(
        encoding="utf-8"
    )
    assert "--unshare-all --unshare-net" in blind
    assert "--clearenv" in blind
    assert "--setenv PLOIDYPATCH_NETWORK_ACCESS none" in blind
    assert '--ro-bind "$role_root/shared_target" /holdout/shared_target' in blind
    assert '--ro-bind "$role_root/candidate_only" /holdout/candidate_only' in blind
    assert '--ro-bind "$role_root/blind_benchmark" /holdout/blind_benchmark' in blind
    assert '--ro-bind "$role_root"' not in blind
    assert "--ro-bind /nas_data" not in blind
    assert '--ro-bind "$role_root/evaluator_only"' not in blind
    assert '"$lock_python" - "$mount_tsv"' in blind
    assert "python3 -" not in blind
    assert "find . -type f ! -path ./SHA256SUMS -printf '%P\\0'" in blind

    reveal = (ROOT / "scripts/run_external_holdout_reveal_v0.5.sh").read_text(
        encoding="utf-8"
    )
    custody_gate = reveal.index("This is the hard truth-access barrier")
    authorization = reveal.index("external_holdout_reveal_authorization.v0.5")
    evaluator_resolution = reveal.index('evaluator_only=$(realpath "$evaluator_only_argument")')
    builder = reveal.index('PLOIDYPATCH_REVEAL_INPUTS_OUTPUT="$reveal_inputs"')
    evaluator = reveal.index('"$model_python" "$evaluator"')
    assert custody_gate < authorization < evaluator_resolution < builder < evaluator
    assert "ready_for_evaluation" in reveal
    assert "not_evaluable_without_rule_relaxation" in reveal
    assert "invalid_run" in reveal
    assert '"evaluator_invoked": False' in reveal
    assert '--blind-run "$blind_run/project"' in reveal
    assert "find . -type f ! -path ./SHA256SUMS -printf '%P\\0'" in reveal
    patch_stage = (
        "post_evaluator_truth_failed_blind_pre_candidate_pre_score_pre_label_execution_patch"
    )
    assert patch_stage == EXECUTION.PATCH_STAGE
    assert patch_stage in blind
    assert patch_stage in reveal
    assert "external_holdout_execution_patch.v0.5" in blind
    assert "external_holdout_execution_patch.v0.5" in reveal


def test_shell_skeletons_are_syntactically_valid() -> None:
    for relative in (
        "scripts/run_external_holdout_blind_isolated_v0.5.sh",
        "scripts/run_external_holdout_reveal_v0.5.sh",
    ):
        completed = subprocess.run(
            ["bash", "-n", relative],
            check=False,
            capture_output=True,
            text=True,
            cwd=ROOT,
        )
        assert completed.returncode == 0, completed.stderr
