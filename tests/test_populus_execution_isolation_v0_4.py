from __future__ import annotations

import csv
import importlib.util
import io
import json
import tarfile
from pathlib import Path

import pytest


ROOT = Path(__file__).parents[1]
FREEZE_PATH = ROOT / "scripts/freeze_populus_execution_implementation_v0.4.py"
CUSTODY_PATH = ROOT / "scripts/finalize_populus_blind_custody_v0.4.py"


def load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


FREEZE = load(FREEZE_PATH, "populus_execution_freeze")
CUSTODY = load(CUSTODY_PATH, "populus_blind_custody")


def write_checksums(root: Path, names: list[str]) -> None:
    with (root / "SHA256SUMS").open("x", encoding="utf-8", newline="") as handle:
        for name in names:
            handle.write(f"{CUSTODY.sha256(root / name)}  {name}\n")


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def test_execution_contract_lists_every_pipeline_family_and_primary_seqids() -> None:
    files = set(FREEZE.implementation_files(ROOT))
    expected = {
        "scripts/prepare_populus_external_normalized_inputs_v0.4.sh",
        "scripts/run_populus_evaluator_wgdi_v0.4.sh",
        "scripts/run_populus_copy_collapse_benchmark_v0.4.sh",
        "scripts/run_populus_blind_pipeline_v0.4.sh",
        "scripts/synthesize_missing_transcript_exons.py",
        "scripts/build_populus_complete_control_reveal_inputs_v0.4.sh",
        "scripts/evaluate_external_v0.4.py",
        "src/ploidypatch/conflict_guard.py",
        "src/ploidypatch/support_ranker.py",
    }
    assert expected <= files
    assert {
        "config/primary_seqids/populus_trichocarpa_v4.0.tsv",
        "config/primary_seqids/salix_purpurea_v5.0.tsv",
        "config/primary_seqids/salix_suchowensis_GCA_017552425.1.tsv",
        "config/primary_seqids/manihot_esculenta_v6.tsv",
        "config/primary_seqids/ricinus_communis_wild_castor.tsv",
    } <= files
    assert "data/derived/external_inputs/populus_v0.4" in FREEZE.FORBIDDEN_PRE_FREEZE_PATHS
    assert all(
        relative in FREEZE.FORBIDDEN_PRE_FREEZE_PATHS
        for relative in (
            "results/evaluator/populus/v0.4/truth_pairs",
            "results/copy_collapse/external/populus_v0.4_blind_rankings",
            "benchmark/structure/copy_collapse_v0.4/ptr_v4/annotation_copy_collapse_seed20260930",
        )
    )


def test_source_archive_rejects_path_traversal(tmp_path: Path) -> None:
    archive = tmp_path / "bad.tar"
    with tarfile.open(archive, "w") as handle:
        member = tarfile.TarInfo("../escape")
        payload = b"bad"
        member.size = len(payload)
        handle.addfile(member, io.BytesIO(payload))
    with pytest.raises(ValueError, match="Unsafe member"):
        FREEZE.extract_archive_safely(archive, tmp_path / "extract")


def test_execution_patch_binds_nonzero_failed_attempt_without_reveal(
    tmp_path: Path,
) -> None:
    attempt = tmp_path / "attempt1"
    attempt.mkdir()
    (attempt / "exit_status.txt").write_text("1\n", encoding="utf-8")
    (attempt / "stderr.log").write_text("candidate-safe failure\n", encoding="utf-8")
    output = tmp_path / "failed.tsv"

    status = FREEZE.freeze_failed_attempt_manifest(attempt, output)

    assert status == 1
    with output.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert {row["relative_path"] for row in rows} == {
        "exit_status.txt",
        "stderr.log",
    }
    assert {row["kind"] for row in rows} == {"file"}
    (attempt / "v04.tsv").write_text("candidate_digest\n", encoding="utf-8")
    with pytest.raises(ValueError, match="pre-score/pre-reveal"):
        FREEZE.freeze_failed_attempt_manifest(attempt, tmp_path / "rejected.tsv")


def test_execution_freeze_is_atomic_complete_and_non_overwriting(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    project = tmp_path / "project"
    code = project / "code"
    implementation = code / "scripts/only.py"
    implementation.parent.mkdir(parents=True)
    implementation.write_text("print('frozen')\n", encoding="utf-8")
    protocol = project / "protocol"
    model = project / "model"
    protocol.mkdir()
    model.mkdir()
    (protocol / "protocol_manifest.json").write_text(
        json.dumps({"schema_version": FREEZE.PROTOCOL_SCHEMA}), encoding="utf-8"
    )
    write_checksums(protocol, ["protocol_manifest.json"])
    (model / "model.json").write_text("{}\n", encoding="utf-8")
    write_checksums(model, ["model.json"])

    monkeypatch.setattr(FREEZE, "implementation_files", lambda _: ["scripts/only.py"])
    monkeypatch.setattr(FREEZE, "verify_git_state", lambda *_: None)
    monkeypatch.setattr(FREEZE, "environment_lock", lambda _: b"@EXPLICIT\n")
    monkeypatch.setattr(FREEZE, "pip_lock", lambda _: b"pkg==1\n")
    monkeypatch.setattr(FREEZE, "chmod_read_only", lambda _: None)

    def fake_archive(_: Path, __: str, output: Path) -> None:
        with tarfile.open(output, "w") as handle:
            handle.add(implementation, arcname="scripts/only.py")

    monkeypatch.setattr(FREEZE, "create_git_archive", fake_archive)
    environments = [
        "ploidypatch-dev",
        "ploidypatch-model",
        "ploidypatch-baseline",
        "ploidypatch-synteny",
        "ploidypatch-syngap",
        "ploidypatch-gemoma",
        "ploidypatch-lifton",
    ]
    env_args: list[str] = []
    for name in environments:
        prefix = project / "envs" / name
        prefix.mkdir(parents=True)
        env_args.extend(("--environment", f"{name}={prefix}"))
    output = project / "freeze"
    arguments = [
        "--project-root",
        str(project),
        "--code-root",
        str(code),
        "--protocol-freeze",
        str(protocol),
        "--composite-model-freeze",
        str(model),
        *env_args,
        "--code-commit",
        "a" * 40,
        "--output-dir",
        str(output),
    ]
    assert FREEZE.main(arguments) == 0
    assert not Path(str(output) + ".partial").exists()
    FREEZE.verify_sha256sums(output)
    manifest = json.loads((output / "execution_manifest.json").read_text(encoding="utf-8"))
    assert manifest["created_before"] == {
        "wgd_pair_enumeration": True,
        "candidate_generation": True,
        "candidate_labels": True,
        "candidate_scores": True,
    }
    assert len(manifest["environments"]) == 7
    with pytest.raises(FileExistsError, match="overwrite"):
        FREEZE.main(arguments)


def make_custody_fixture(tmp_path: Path) -> dict[str, Path]:
    execution = tmp_path / "execution"
    protocol = tmp_path / "protocol"
    model = tmp_path / "model"
    for root in (execution, protocol, model):
        root.mkdir()
    environment = tmp_path / "env"
    environment.mkdir()
    (execution / "source").mkdir()
    (execution / "execution_manifest.json").write_text(
        json.dumps(
            {
                "schema_version": CUSTODY.EXECUTION_SCHEMA,
                "environments": [
                    {"name": "ploidypatch-dev", "host_prefix": str(environment)}
                ],
            }
        ),
        encoding="utf-8",
    )
    write_checksums(execution, ["execution_manifest.json"])
    (protocol / "policy.tsv").write_text("field\tvalue\n", encoding="utf-8")
    write_checksums(protocol, ["policy.tsv"])
    (model / "model.json").write_text("{}\n", encoding="utf-8")
    write_checksums(model, ["model.json"])

    project = tmp_path / "blind_project"
    scores = project / CUSTODY.SCORES_RELATIVE
    decisions = project / CUSTODY.POOL_DECISIONS_RELATIVE
    scores.parent.mkdir(parents=True)
    decisions.parent.mkdir(parents=True)
    write_tsv(
        decisions,
        [
            {
                "consensus_digest": "a",
                "status": "accepted",
                "conflict_set_digest": "x",
            }
        ],
    )
    pool_manifest = project / CUSTODY.POOL_MANIFEST_RELATIVE
    pool_manifest.write_text(
        json.dumps(
            {
                "schema_version": CUSTODY.POOL_SCHEMA,
                "outputs": {"decisions": {"sha256": CUSTODY.sha256(decisions)}},
            }
        ),
        encoding="utf-8",
    )
    write_tsv(scores, [{"candidate_digest": "a", "v04_primary_rank_score": 1.0}])
    score_manifest = project / CUSTODY.SCORE_MANIFEST_RELATIVE
    score_manifest.write_text(
        json.dumps(
            {
                "schema_version": CUSTODY.SCORE_SCHEMA,
                "truth_access": False,
                "inputs": {
                    "pool_decisions": CUSTODY.sha256(decisions),
                    "pool_manifest": CUSTODY.sha256(pool_manifest),
                },
                "outputs": {"scores": {"sha256": CUSTODY.sha256(scores)}},
            }
        ),
        encoding="utf-8",
    )
    mounts = {}
    for key, name in (
        ("shared", "shared_target"),
        ("candidate", "candidate_only"),
        ("blind", "blind"),
    ):
        path = tmp_path / name
        path.mkdir()
        mounts[key] = path
    command = tmp_path / "command.txt"
    command.write_text("bwrap --unshare-all --unshare-net\n", encoding="utf-8")
    blind_sums = tmp_path / "BLIND_SHA256SUMS"
    blind_sums.write_text("0" * 64 + "  shared_target/genome.fa\n", encoding="utf-8")
    return {
        "execution": execution,
        "protocol": protocol,
        "model": model,
        "project": project,
        "scores": scores,
        "shared": mounts["shared"],
        "candidate": mounts["candidate"],
        "blind": mounts["blind"],
        "command": command,
        "blind_sums": blind_sums,
    }


def custody_args(paths: dict[str, Path], output: Path) -> list[str]:
    return [
        "--blind-project-root",
        str(paths["project"]),
        "--execution-freeze",
        str(paths["execution"]),
        "--protocol-freeze",
        str(paths["protocol"]),
        "--composite-model-freeze",
        str(paths["model"]),
        "--shared-target",
        str(paths["shared"]),
        "--candidate-only",
        str(paths["candidate"]),
        "--blind-annotation-root",
        str(paths["blind"]),
        "--blind-role-checksums",
        str(paths["blind_sums"]),
        "--runner-command",
        str(paths["command"]),
        "--bwrap-version",
        "bubblewrap 0.11",
        "--output",
        str(output),
    ]


def test_custody_seals_exact_scores_pool_and_negative_mount_claims(tmp_path: Path) -> None:
    paths = make_custody_fixture(tmp_path)
    output = tmp_path / "custody.json"
    assert CUSTODY.main(custody_args(paths, output)) == 0
    manifest = json.loads(output.read_text(encoding="utf-8"))
    assert manifest["schema_version"] == "ploidypatch.blind_run_custody.v1"
    assert manifest["network_access"] is False
    assert manifest["nas_data_mounted"] is False
    assert manifest["truth_mounted"] is False
    assert manifest["blind_outputs"]["scores_sha256"] == CUSTODY.sha256(paths["scores"])
    assert manifest["bubblewrap"]["required_flags"] == [
        "--unshare-all",
        "--unshare-net",
    ]


def test_custody_rejects_label_leakage_in_blind_scores(tmp_path: Path) -> None:
    paths = make_custody_fixture(tmp_path)
    scores = paths["scores"]
    scores.unlink()
    write_tsv(scores, [{"candidate_digest": "a", "truth_label": 1}])
    score_manifest = paths["project"] / CUSTODY.SCORE_MANIFEST_RELATIVE
    manifest = json.loads(score_manifest.read_text(encoding="utf-8"))
    manifest["outputs"]["scores"]["sha256"] = CUSTODY.sha256(scores)
    score_manifest.write_text(json.dumps(manifest), encoding="utf-8")
    with pytest.raises(ValueError, match="truth/label"):
        CUSTODY.main(custody_args(paths, tmp_path / "custody.json"))


def test_bwrap_launcher_mounts_only_blind_roles_and_disables_network() -> None:
    script = (ROOT / "scripts/run_populus_blind_isolated_v0.4.sh").read_text(
        encoding="utf-8"
    )
    assert "--unshare-all --unshare-net" in script
    assert '--ro-bind "$shared_target"' in script
    assert '--ro-bind "$candidate_only"' in script
    assert '--ro-bind "$blind_annotation"' in script
    assert '--ro-bind "$stage"' not in script
    assert '--ro-bind /nas_data' not in script
    assert 'PLOIDYPATCH_BLIND_RUNNER 1' in script
    assert '--setenv PYTHONPATH "$internal/code/src"' in script
    assert "--symlink usr/bin /bin" in script
    assert "[EXECUTION_FREEZE]" in script
    assert "seven frozen environments are required" in script
    assert "external_inputs/populus/v0.4" not in script


def test_reveal_barrier_precedes_truth_builder_and_uses_model_environment() -> None:
    script = (ROOT / "scripts/run_populus_external_reveal_v0.4.sh").read_text(
        encoding="utf-8"
    )
    barrier = script.index("populus_reveal_authorization.v0.4")
    builder = script.index('PLOIDYPATCH_BLIND_RUN_ROOT="$blind"')
    evaluator = script.index('"$model_python" "$execution/source/scripts/evaluate_external_v0.4.py"')
    assert barrier < builder < evaluator
    assert '"not_evaluable_without_rule_relaxation"' in script
    assert "if [[ $reveal_status != ready_for_evaluation ]]" in script
    not_evaluable_branch = script.index("if [[ $reveal_status != ready_for_evaluation ]]")
    assert not_evaluable_branch < evaluator
    assert '"evaluator_invoked": False' in script
    assert "if [[ $reveal_status == invalid_run ]]" in script
    invalid_branch = script.index("if [[ $reveal_status == invalid_run ]]")
    invalid_exit = script.index('echo "Populus evaluator sentinel/custody run is invalid', invalid_branch)
    next_publish = script.index('mv "$working" "$result"', invalid_exit)
    assert invalid_exit < next_publish
    assert "exit 1" in script[invalid_exit:next_publish]
    assert 'PLOIDYPATCH_BLIND_RUN_ROOT="$blind"' in script
    assert 'PLOIDYPATCH_EXECUTION_FREEZE_OVERRIDE="$execution"' in script
    assert "external_inputs/populus/v0.4" not in script
    assert 'model_prefix=$(awk' in script
    assert "export PYTHONPATH=$execution/source/src" in script
