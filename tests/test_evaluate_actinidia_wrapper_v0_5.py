from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import subprocess

import pytest


SCRIPT = Path(__file__).parents[1] / "scripts/evaluate_actinidia_external_v0.5.py"
SPEC = importlib.util.spec_from_file_location("evaluate_actinidia_v05", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def _json(path: Path, value: object) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value), encoding="utf-8")
    return path


def _text(path: Path, value: str = "x\n") -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(value, encoding="utf-8")
    return path


def _freeze(root: Path) -> None:
    with (root / "SHA256SUMS").open("x", encoding="utf-8", newline="") as handle:
        for path in sorted(root.rglob("*")):
            if path.is_file() and path.name != "SHA256SUMS":
                handle.write(f"{MODULE.sha256_file(path)}  {path.relative_to(root).as_posix()}\n")


def _fixture(tmp_path: Path) -> dict[str, Path]:
    model = tmp_path / "model"
    model.mkdir()
    _text(model / "model.bin")
    _freeze(model)

    protocol = tmp_path / "protocol"
    policy = _text(
        protocol
        / "protocol_artifacts/config/actinidia_external_validation_policy_v0.5.tsv",
        "field\tvalue\npolicy_id\tploidypatch_actinidia_external_validation_v0.5\n",
    )
    _json(
        protocol / "protocol_manifest.json",
        {
            "schema_version": MODULE.PROTOCOL_SCHEMA,
            "holdout_id": MODULE.HOLDOUT_ID,
            "policy_id": MODULE.POLICY_ID,
            "model_version": "PloidyPatch_ranker_v0.4",
            "composite_model_SHA256SUMS_sha256": MODULE.sha256_file(
                model / "SHA256SUMS"
            ),
        },
    )
    _freeze(protocol)

    execution = tmp_path / "execution"
    execution.mkdir()
    _json(
        execution / "execution_manifest.json",
        {
            "schema_version": MODULE.EXECUTION_SCHEMA,
            "holdout_id": MODULE.HOLDOUT_ID,
            "protocol_SHA256SUMS_sha256": MODULE.sha256_file(
                protocol / "SHA256SUMS"
            ),
        },
    )
    _freeze(execution)

    blind = tmp_path / "blind"
    scores = _text(blind / "project/scores.tsv")
    score_manifest = _json(blind / "project/scores.tsv.manifest.json", {"ok": True})
    decisions = _text(blind / "project/decisions.tsv")
    pool_manifest = _json(blind / "project/pool.manifest.json", {"ok": True})
    custody = _json(
        tmp_path / "custody.json",
        {
            "schema_version": MODULE.CUSTODY_SCHEMA,
            "holdout_id": MODULE.HOLDOUT_ID,
            "policy_id": MODULE.POLICY_ID,
            "frozen_inputs": {
                "execution_SHA256SUMS_sha256": MODULE.sha256_file(
                    execution / "SHA256SUMS"
                ),
                "protocol_SHA256SUMS_sha256": MODULE.sha256_file(
                    protocol / "SHA256SUMS"
                ),
                "composite_model_SHA256SUMS_sha256": MODULE.sha256_file(
                    model / "SHA256SUMS"
                ),
            },
            "blind_outputs": {
                name: {
                    "relative_path": path.relative_to(blind).as_posix(),
                    "sha256": MODULE.sha256_file(path),
                }
                for name, path in {
                    "scores": scores,
                    "score_manifest": score_manifest,
                    "pool_decisions": decisions,
                    "pool_manifest": pool_manifest,
                }.items()
            },
        },
    )

    reveal = tmp_path / "reveal"
    inputs = {
        "labels": _text(reveal / "labels.tsv"),
        "primary_pool_score": _json(reveal / "primary.json", {}),
        "legacy_pool_score": _json(reveal / "legacy.json", {}),
        "secondary:miniprot": _json(reveal / "miniprot.json", {}),
    }
    evaluability = _json(reveal / "evaluability.json", {})
    _json(
        reveal / "reveal_input_manifest.json",
        {
            "schema_version": MODULE.REVEAL_SCHEMA,
            "holdout_id": MODULE.HOLDOUT_ID,
            "policy_id": MODULE.POLICY_ID,
            "formal_status": "ready_for_evaluation",
            "custody_manifest_sha256": MODULE.sha256_file(custody),
            "evaluability": {
                "relative_path": evaluability.relative_to(reveal).as_posix(),
                "sha256": MODULE.sha256_file(evaluability),
            },
            "evaluation_inputs": {
                name: {
                    "relative_path": path.relative_to(reveal).as_posix(),
                    "sha256": MODULE.sha256_file(path),
                }
                for name, path in inputs.items()
            },
        },
    )
    _freeze(reveal)
    return {
        "execution": execution,
        "protocol": protocol,
        "model": model,
        "blind": blind,
        "custody": custody,
        "reveal": reveal,
        "policy": policy,
        "output": tmp_path / "output",
    }


def test_wrapper_resolves_only_hash_bound_inputs(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    paths = _fixture(tmp_path)
    called: list[list[str]] = []

    def fake_run(command: list[str], *, check: bool) -> subprocess.CompletedProcess[str]:
        assert check is True
        called.append(command)
        return subprocess.CompletedProcess(command, 0)

    monkeypatch.setattr(MODULE.subprocess, "run", fake_run)
    MODULE.evaluate(
        execution=paths["execution"],
        protocol=paths["protocol"],
        model=paths["model"],
        blind_run=paths["blind"],
        custody_path=paths["custody"],
        reveal_inputs=paths["reveal"],
        output=paths["output"],
    )
    assert len(called) == 1
    command = called[0]
    assert "--secondary-score" in command
    assert any(value.startswith("miniprot=") for value in command)
    assert str(paths["policy"]) in command
    roots = [
        command[index + 1]
        for index, value in enumerate(command)
        if value == "--input-root"
    ]
    assert roots == [
        f"blind_run={paths['blind']}",
        f"reveal_inputs={paths['reveal']}",
        f"protocol_freeze={paths['protocol']}",
    ]


def test_wrapper_rejects_reveal_mutation_before_evaluator(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    paths = _fixture(tmp_path)
    (paths["reveal"] / "labels.tsv").write_text("mutated\n", encoding="utf-8")
    monkeypatch.setattr(
        MODULE.subprocess,
        "run",
        lambda *args, **kwargs: pytest.fail("evaluator must not run"),
    )
    with pytest.raises(ValueError, match="SHA256"):
        MODULE.evaluate(
            execution=paths["execution"],
            protocol=paths["protocol"],
            model=paths["model"],
            blind_run=paths["blind"],
            custody_path=paths["custody"],
            reveal_inputs=paths["reveal"],
            output=paths["output"],
        )
