from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
import subprocess
from types import SimpleNamespace

import pytest

from ploidypatch.artifact_manifest import sha256_file, write_sha256sums
from ploidypatch.coffea_h1_framework import (
    BLIND_OUTPUTS,
    PIPELINE_ENTRIES,
    RAW_PREDICTION_TREE_KEYS,
)
from ploidypatch.core_h1_evaluation import evaluate_core_h1_scores


ROOT = Path(__file__).resolve().parents[1]


def score(path: Path, values: list[bool]) -> None:
    recovered = sum(values)
    payload = {
        "schema_version": "ploidypatch.annotation_repair_score.v5",
        "quality_gate": {"grade": "pass"},
        "collateral_changes": {
            "baseline_transcript_structures_missing_from_candidate": 0
        },
        "event_recovery": {
            "events": len(values),
            "complete_cds_chain_recovery": recovered,
            "complete_cds_chain_recall": recovered / len(values),
        },
        "event_details": [
            {
                "event_id": f"event-{index}",
                "event_type": "annotation_copy_collapse",
                "complete_cds_chain_recovery": value,
            }
            for index, value in enumerate(values, start=1)
        ],
    }
    path.write_text(json.dumps(payload) + "\n", encoding="utf-8")


def test_coffea_blind_output_universe_is_exact() -> None:
    assert len(BLIND_OUTPUTS) == 20
    assert set(PIPELINE_ENTRIES) == {
        "blind_pipeline", "reveal_input_builder", "evaluator"
    }
    assert len(RAW_PREDICTION_TREE_KEYS) == 6
    for scope in ("combined", "bua_only", "mauritiana_only"):
        for arm in ("retain_distinct", "suppress_overlap"):
            prefix = f"{scope}_{arm}"
            assert {f"{prefix}_pool", f"{prefix}_decisions", f"{prefix}_manifest"} <= set(
                BLIND_OUTPUTS
            )


def test_core_h1_evaluates_combined_only_and_reports_secondary(tmp_path: Path) -> None:
    paths: dict[str, Path] = {}
    for scope in ("combined", "bua_only", "mauritiana_only"):
        for arm in ("retain_distinct", "suppress_overlap"):
            name = f"{scope}_{arm}"
            path = tmp_path / f"{name}.json"
            values = [True] * 8 if arm == "retain_distinct" else [False] * 8
            score(path, values)
            paths[name] = path
    result = evaluate_core_h1_scores(
        score_paths=paths,
        holdout_id="coffea_et39_v1.0",
        policy_id="ploidypatch_coffea_external_core_h1_v1.0",
        schema_version="test.coffea",
        primary_scope="combined",
        bootstrap_seed=20260912,
        bootstrap_replicates=20_000,
        output=tmp_path / "evaluation.json",
    )
    assert result["success"] is True
    assert result["formal_outcome"] == "formal_positive_external_result"
    assert result["paired_event_bootstrap"]["observed_delta"] == 1.0
    assert result["paired_event_bootstrap"]["ci_lower"] == 1.0
    assert set(result["arms"]) == {"combined", "bua_only", "mauritiana_only"}
    assert result["ranker_or_model_executed"] is False


def test_core_h1_rejects_collateral_loss(tmp_path: Path) -> None:
    paths: dict[str, Path] = {}
    for scope in ("combined", "bua_only", "mauritiana_only"):
        for arm in ("retain_distinct", "suppress_overlap"):
            name = f"{scope}_{arm}"
            path = tmp_path / f"{name}.json"
            score(path, [True, False])
            paths[name] = path
    value = json.loads(paths["bua_only_retain_distinct"].read_text())
    value["collateral_changes"][
        "baseline_transcript_structures_missing_from_candidate"
    ] = 1
    paths["bua_only_retain_distinct"].write_text(json.dumps(value))
    with pytest.raises(ValueError, match="collateral"):
        evaluate_core_h1_scores(
            score_paths=paths,
            holdout_id="x",
            policy_id="y",
            schema_version="z",
            primary_scope="combined",
            bootstrap_seed=20260912,
            bootstrap_replicates=20_000,
            output=tmp_path / "evaluation.json",
        )


def test_coffea_entries_and_blind_scripts_exist_and_parse() -> None:
    python_files = [
        "scripts/build_coffea_blind_role_root_v1.0.py",
        "scripts/build_coffea_evaluator_role_root_v1.0.py",
        "scripts/finalize_coffea_blind_custody_v1.0.py",
        "scripts/freeze_coffea_external_execution_v1.0.py",
        "scripts/build_coffea_complete_control_reveal_inputs_v1.0.py",
        "scripts/evaluate_coffea_external_h1_v1.0.py",
        "scripts/run_coffea_external_reveal_v1.0.py",
    ]
    subprocess.run(
        ["python", "-m", "py_compile", *python_files], cwd=ROOT, check=True
    )
    shell_files = [
        "scripts/run_coffea_blind_pipeline_v1.0.sh",
        "scripts/run_coffea_blind_isolated_v1.0.sh",
        "scripts/run_coffea_candidate_methods_v1.0.sh",
        "scripts/run_coffea_evaluator_wgdi_v1.0.sh",
        "scripts/run_coffea_evaluator_pipeline_v1.0.sh",
    ]
    subprocess.run(["bash", "-n", *shell_files], cwd=ROOT, check=True)
    blind_text = "\n".join((ROOT / path).read_text() for path in shell_files[:3]).casefold()
    assert "ploidypatch-model" not in blind_text
    assert "score_candidates" not in blind_text
    assert "topology_features" not in blind_text


def test_coffea_wgdi_runner_allows_conda_executable_symlinks_only() -> None:
    text = (ROOT / "scripts/run_coffea_evaluator_wgdi_v1.0.sh").read_text()
    assert 'for executable in "$python" "$diamond" "$wgdi" "$parallel"' in text
    assert '[[ -x $executable ]]' in text
    assert 'resolved=$(realpath "$executable")' in text
    assert 'for path in "$ks_helper" "$input/SHA256SUMS"' in text
    assert '[[ -s $path && ! -L $path ]]' in text


def test_coffea_blind_launcher_mounts_one_canonical_role_manifest() -> None:
    text = (ROOT / "scripts/run_coffea_blind_isolated_v1.0.sh").read_text()
    assert "/holdout/blind_role_manifest.json" not in text
    assert text.count(
        '--ro-bind "$role_root/role_manifest.json" /holdout/role_manifest.json'
    ) == 1
    assert "'host_role_manifest_sha256':sha('/holdout/role_manifest.json')" in text


def test_coffea_execution_patch_chain_is_fail_closed_before_labels() -> None:
    freezer = (ROOT / "scripts/freeze_coffea_external_execution_v1.0.py").read_text()
    framework = (ROOT / "src/ploidypatch/coffea_h1_framework.py").read_text()
    assert "post_evaluator_truth_failed_blind_pre_candidate_pre_label_execution_patch" in freezer
    assert "post_evaluator_truth_failed_blind_partial_candidate_pre_label_execution_patch_2" in freezer
    assert "post_evaluator_truth_two_complete_blind_runs_pre_label_" in freezer
    assert "post_blind_custody_reveal_authorized_pre_evaluator_environment_patch_4" in freezer
    assert "post_blind_custody_pre_truth_authorization_custody_lineage_patch_5" in freezer
    assert '"candidate_generation_completed_before_patch": patch_sequence >= 3' in freezer
    assert '"partial_candidate_generation_before_patch": patch_sequence == 2' in freezer
    assert '"truth_labels_accessed_before_patch": False' in freezer
    assert '"scientific_rules_or_thresholds_changed": patch_sequence == 3' in freezer
    assert '"biological_rules_or_performance_thresholds_changed": False' in freezer
    assert all(
        value in freezer
        for value in (
            "PATCH1_ALLOWED_PATHS", "PATCH2_ALLOWED_PATHS", "PATCH3_ALLOWED_PATHS",
            "PATCH4_ALLOWED_PATHS",
            "PATCH5_ALLOWED_PATHS",
        )
    )
    assert "Unexpected Coffea failed-attempt universe" in freezer
    assert "primary_combined_replay_outputs" in freezer
    assert "candidate_generation_completed_before_patch" in framework
    assert "scientific_rules_or_thresholds_changed" in framework
    assert (ROOT / "docs/COFFEA_EXECUTION_PATCH_SINGLE_REFERENCE_SCOPE_v1.0.md").is_file()
    assert (ROOT / "docs/COFFEA_BLIND_REPRODUCIBILITY_AMENDMENT_v1.1.md").is_file()
    assert (ROOT / "docs/COFFEA_REVEAL_EXECUTION_PATCH_v1.2.md").is_file()
    assert (ROOT / "docs/COFFEA_CUSTODY_LINEAGE_PATCH_v1.3.md").is_file()


def test_coffea_patch4_accepts_only_reveal_authorization(tmp_path: Path) -> None:
    import importlib.util

    script = ROOT / "scripts/freeze_coffea_external_execution_v1.0.py"
    spec = importlib.util.spec_from_file_location("coffea_freezer_patch4", script)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    root = tmp_path / "failed_reveal"
    root.mkdir()
    authorization = {
        "schema_version": "ploidypatch.coffea_h1_reveal_authorization.v1.0",
        "holdout_id": module.HOLDOUT_ID,
        "truth_reveal_authorized": True,
        "ranker_or_model_authorized": False,
        "h2_or_topology_ranking_authorized": False,
    }
    (root / "reveal_authorization.json").write_text(
        json.dumps(authorization), encoding="utf-8"
    )
    status, rows, evidence, replay = module.failed_attempt_files(
        root, patch_sequence=4
    )
    assert status == 1
    assert rows == evidence
    assert rows[0]["relative_path"] == "reveal_authorization.json"
    assert replay is None
    (root / "unexpected.txt").write_text("no\n", encoding="utf-8")
    with pytest.raises(ValueError, match="reveal failure universe differs"):
        module.failed_attempt_files(root, patch_sequence=4)


def test_coffea_patch5_accepts_only_pre_truth_custody_failure(tmp_path: Path) -> None:
    import importlib.util

    script = ROOT / "scripts/freeze_coffea_external_execution_v1.0.py"
    spec = importlib.util.spec_from_file_location("coffea_freezer_patch5", script)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    root = tmp_path / "failed_reveal"
    root.mkdir()
    error_log = root / "error.log"
    error_log.write_text(
        "ValueError: Coffea custody failed pre-truth authorization\n",
        encoding="utf-8",
    )
    status, rows, evidence, replay = module.failed_attempt_files(
        root, patch_sequence=5
    )
    assert status == 1
    assert rows == evidence
    assert rows[0]["sha256"] == sha256_file(error_log)
    assert replay is None
    error_log.write_text("different failure\n", encoding="utf-8")
    with pytest.raises(ValueError, match="reveal failure differs"):
        module.failed_attempt_files(root, patch_sequence=5)


def test_coffea_reveal_uses_explicit_patch5_custody_lineage() -> None:
    reveal = (ROOT / "scripts/run_coffea_external_reveal_v1.0.py").read_text()
    framework = (ROOT / "src/ploidypatch/coffea_h1_framework.py").read_text()
    assert 'patch.get("patch_sequence") == 5' in reveal
    assert "blind_custody_execution_SHA256SUMS_sha256" in reveal
    assert 'patch.get("blind_custody_manifest_sha256")' in reveal
    assert "custody_execution_chain_validation_fixed" in framework


def test_coffea_reveal_accepts_only_contained_conda_python_symlink(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    import importlib.util

    script = ROOT / "scripts/run_coffea_external_reveal_v1.0.py"
    spec = importlib.util.spec_from_file_location("coffea_reveal", script)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    prefix = tmp_path / "env"
    real_python = prefix / "libexec/python-real"
    real_python.parent.mkdir(parents=True)
    real_python.write_text("python\n", encoding="utf-8")
    real_python.chmod(0o755)
    python = prefix / "bin/python"
    python.parent.mkdir()
    try:
        python.symlink_to(real_python)
    except OSError:
        pytest.skip("symlink creation is unavailable")
    explicit = b"@EXPLICIT\nmock\n"
    execution = tmp_path / "execution"
    (execution / "environment_locks").mkdir(parents=True)
    lock = execution / "environment_locks/dev.explicit.txt"
    lock.write_bytes(explicit)
    manifest = {
        "environments": [
            {
                "name": "ploidypatch-dev",
                "host_prefix": str(prefix),
                "explicit_relative_path": "environment_locks/dev.explicit.txt",
                "explicit_sha256": hashlib.sha256(explicit).hexdigest(),
            }
        ]
    }
    (execution / "execution_manifest.json").write_text(
        json.dumps(manifest), encoding="utf-8"
    )
    monkeypatch.setattr(
        module.subprocess,
        "run",
        lambda *args, **kwargs: SimpleNamespace(stdout=explicit),
    )
    assert module.dev_environment(execution) == real_python.resolve()
    python.unlink()
    escaped = tmp_path / "outside-python"
    escaped.write_text("python\n", encoding="utf-8")
    escaped.chmod(0o755)
    python.symlink_to(escaped)
    with pytest.raises(ValueError, match="escapes its environment"):
        module.dev_environment(execution)


def test_coffea_root_manifests_use_canonical_nested_writer() -> None:
    for relative in (
        "scripts/run_coffea_blind_isolated_v1.0.sh",
        "scripts/run_coffea_evaluator_pipeline_v1.0.sh",
        "scripts/run_coffea_reconciliation_isolated_v1.1.sh",
    ):
        text = (ROOT / relative).read_text()
        assert "find . -type f ! -name SHA256SUMS" not in text
        assert "write_sha256sums" in text
        assert "verify_sha256sums" in text


def test_coffea_patch2_failed_attempt_binds_primary_combined_outputs(
    tmp_path: Path,
) -> None:
    import importlib.util

    script = ROOT / "scripts/freeze_coffea_external_execution_v1.0.py"
    spec = importlib.util.spec_from_file_location("coffea_freezer", script)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    root = tmp_path / "failed"
    root.mkdir()
    (root / "bwrap_command.txt").write_text("bwrap\n", encoding="utf-8")
    (root / "exit_status.txt").write_text("1\n", encoding="utf-8")
    (root / "stderr.log").write_text("single-reference failure\n", encoding="utf-8")
    (root / "stdout.log").write_text("blind projections complete\n", encoding="utf-8")
    (root / "mount_manifest.json").write_text(
        json.dumps(
            {
                "schema_version": "ploidypatch.coffea_core_h1_mount_manifest.v1.0",
                "network_access": False,
                "ranker_or_model_mounted": False,
            }
        ),
        encoding="utf-8",
    )
    (root / "namespace_role_validation.json").write_text(
        json.dumps(
            {
                "schema_version": "ploidypatch.coffea_core_h1_namespace_validation.v1.0",
                "complete_target_annotation_visible": False,
                "evaluator_only_visible": False,
                "truth_visible": False,
                "labels_visible": False,
                "nas_data_visible": False,
                "model_visible": False,
                "ranker_visible": False,
            }
        ),
        encoding="utf-8",
    )
    (root / "project").mkdir()
    (root / "project/pipeline_commands.tsv").write_text(
        "stage\tentry\n", encoding="utf-8"
    )
    failure = root / "project/results/copy_collapse/external/coffea_v1.0_h1.invalid_run"
    relative_files = {
        "raw_predictions.manifest.json": "raw\n",
        "combined/retain_distinct/blind/candidate.gff3": "retain pool\n",
        "combined/retain_distinct/blind/decisions.tsv": "retain decisions\n",
        "combined/retain_distinct/blind/candidate.gff3.manifest.json": "retain manifest\n",
        "combined/suppress_overlap/blind/candidate.gff3": "legacy pool\n",
        "combined/suppress_overlap/blind/decisions.tsv": "legacy decisions\n",
        "combined/suppress_overlap/blind/candidate.gff3.manifest.json": "legacy manifest\n",
        "scopes/bua_only/methods/miniprot/stderr.log": "requires two inputs\n",
        "scopes/bua_only/methods/miniprot/stdout.json": "\n",
    }
    for relative, text in relative_files.items():
        path = failure / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text, encoding="utf-8")

    status, rows, evidence, replay = module.failed_attempt_files(
        root, patch_sequence=2
    )
    assert status == 1
    assert len(rows) == 16
    assert len(evidence) == 12
    assert replay == {
        "combined_retain_distinct_pool": sha256_file(
            failure / "combined/retain_distinct/blind/candidate.gff3"
        ),
        "combined_retain_distinct_decisions": sha256_file(
            failure / "combined/retain_distinct/blind/decisions.tsv"
        ),
        "combined_suppress_overlap_pool": sha256_file(
            failure / "combined/suppress_overlap/blind/candidate.gff3"
        ),
        "combined_suppress_overlap_decisions": sha256_file(
            failure / "combined/suppress_overlap/blind/decisions.tsv"
        ),
    }


def test_coffea_patch2_custody_requires_exact_primary_replay(tmp_path: Path) -> None:
    import importlib.util

    script = ROOT / "scripts/finalize_coffea_blind_custody_v1.0.py"
    spec = importlib.util.spec_from_file_location("coffea_custody", script)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    keys = {
        "combined_retain_distinct_pool",
        "combined_retain_distinct_decisions",
        "combined_suppress_overlap_pool",
        "combined_suppress_overlap_decisions",
    }
    paths: dict[str, Path] = {}
    for key in keys:
        path = tmp_path / key
        path.write_text(key + "\n", encoding="utf-8")
        paths[key] = path
    manifest = {
        "freeze_stage": module.PATCH2_STAGE,
        "execution_patch": {
            "primary_combined_replay_outputs": {
                key: sha256_file(path) for key, path in paths.items()
            }
        },
    }
    module.validate_patch2_primary_replay(manifest, paths)
    paths["combined_retain_distinct_pool"].write_text("changed\n", encoding="utf-8")
    with pytest.raises(ValueError, match="changed a primary combined blind output"):
        module.validate_patch2_primary_replay(manifest, paths)


def test_coffea_patch3_accepts_only_complete_pre_custody_blind_run(
    tmp_path: Path,
) -> None:
    import importlib.util

    script = ROOT / "scripts/freeze_coffea_external_execution_v1.0.py"
    spec = importlib.util.spec_from_file_location("coffea_freezer_patch3", script)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    root = tmp_path / "run_b"
    root.mkdir()
    for relative, value in {
        "bwrap_command.txt": "bwrap\n",
        "exit_status.txt": "0\n",
        "stderr.log": "\n",
        "stdout.log": "blind complete\n",
        "project/pipeline_commands.tsv": "stage\tentry\n",
    }.items():
        path = root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(value, encoding="utf-8")
    (root / "mount_manifest.json").write_text(
        json.dumps(
            {
                "schema_version": "ploidypatch.coffea_core_h1_mount_manifest.v1.0",
                "network_access": False,
                "ranker_or_model_mounted": False,
            }
        ),
        encoding="utf-8",
    )
    (root / "namespace_role_validation.json").write_text(
        json.dumps(
            {
                "schema_version": "ploidypatch.coffea_core_h1_namespace_validation.v1.0",
                **{
                    key: False
                    for key in (
                        "complete_target_annotation_visible",
                        "evaluator_only_visible",
                        "truth_visible",
                        "labels_visible",
                        "nas_data_visible",
                        "model_visible",
                        "ranker_visible",
                    )
                },
            }
        ),
        encoding="utf-8",
    )
    completed = root / "project/results/copy_collapse/external/coffea_v1.0_h1"
    (completed / "raw_predictions.manifest.json").parent.mkdir(parents=True)
    (completed / "raw_predictions.manifest.json").write_text("{}\n", encoding="utf-8")
    for scope in ("combined", "bua_only", "mauritiana_only"):
        for arm in ("retain_distinct", "suppress_overlap"):
            blind = completed / scope / arm / "blind"
            blind.mkdir(parents=True)
            for name in (
                "candidate.gff3",
                "decisions.tsv",
                "candidate.gff3.manifest.json",
            ):
                (blind / name).write_text(name + "\n", encoding="utf-8")
    status, rows, evidence, replay = module.failed_attempt_files(
        root, patch_sequence=3
    )
    assert status == 0
    assert replay is None
    assert len(rows) == 26
    assert len(evidence) == 14
    (root / "custody_manifest.json").write_text("{}\n", encoding="utf-8")
    with pytest.raises(ValueError, match="unexpectedly completed custody"):
        module.failed_attempt_files(root, patch_sequence=3)


def test_filtered_role_manifest_keeps_exact_candidate_safe_rows(tmp_path: Path) -> None:
    import importlib.util

    script = ROOT / "scripts/build_coffea_blind_role_root_v1.0.py"
    spec = importlib.util.spec_from_file_location("coffea_role", script)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    source = tmp_path / "roles.tsv"; output = tmp_path / "filtered.tsv"
    fields = ["role", "species_id", "bundle_id", "artifact", "staged_relative_path", "bytes", "sha256"]
    rows = [
        {"role": "target", "species_id": "Coffea_arabica_ET39", "bundle_id": "target_et39", "artifact": artifact,
         "staged_relative_path": f"x/{artifact}", "bytes": "1", "sha256": "a" * 64}
        for artifact in ("genome", "gff3", "protein")
    ]
    for species, bundle in (("Coffea_eugenioides_BuA", "candidate_bua"), ("Coffea_mauritiana", "candidate_mauritiana")):
        rows.extend(
            {"role": "candidate_reference", "species_id": species, "bundle_id": bundle, "artifact": artifact,
             "staged_relative_path": f"y/{artifact}", "bytes": "1", "sha256": "b" * 64}
            for artifact in ("genome", "gff3", "protein")
        )
    rows.append({"role": "evaluator_reference", "species_id": "Gardenia", "bundle_id": "evaluator", "artifact": "genome",
                 "staged_relative_path": "z", "bytes": "1", "sha256": "c" * 64})
    with source.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    module.filtered_role_manifest(source, output)
    with output.open(encoding="utf-8", newline="") as handle:
        kept = list(csv.DictReader(handle, delimiter="\t"))
    assert len(kept) == 7
    assert {(row["role"], row["artifact"]) for row in kept} == {
        ("target", "genome"),
        ("candidate_reference", "genome"),
        ("candidate_reference", "gff3"),
        ("candidate_reference", "protein"),
    }
