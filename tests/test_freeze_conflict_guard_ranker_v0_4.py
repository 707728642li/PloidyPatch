from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest


MODULE_PATH = (
    Path(__file__).parents[1] / "scripts" / "freeze_conflict_guard_ranker_v0.4.py"
)
SPEC = importlib.util.spec_from_file_location("guard_freeze", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_guard_policy_has_no_probability_threshold_or_auto_approval() -> None:
    commit = "a" * 40
    policy = MODULE.guard_policy(commit)
    assert policy["code_commit"] == commit
    assert policy["policy"]["automatic_approval"] is False
    assert policy["policy"]["calibrated_probability"] is False
    assert policy["policy"]["decision_threshold"] is None
    assert policy["production_replay"]["ranking_exact"] is True


def test_verify_sha256sums_rejects_modified_input(tmp_path: Path) -> None:
    artifact = tmp_path / "artifact.txt"
    artifact.write_text("frozen\n", encoding="utf-8")
    (tmp_path / "SHA256SUMS").write_text(
        f"{MODULE.sha256(artifact)}  artifact.txt\n", encoding="utf-8"
    )
    assert MODULE.verify_sha256sums(tmp_path) == {
        "artifact.txt": MODULE.sha256(artifact)
    }
    artifact.write_text("changed\n", encoding="utf-8")
    with pytest.raises(ValueError, match="verification failed"):
        MODULE.verify_sha256sums(tmp_path)


def test_development_fixture_expresses_required_exact_winner_gate(
    tmp_path: Path,
) -> None:
    fixture = {
        "schema_version": "ploidypatch.conflict_winner_guard_evaluation.v1",
        "code_commit": "b" * 40,
        "all_development_gates_pass": True,
        "datasets": {
            "toy": {
                "conflict_winner_identity": {
                    "mismatch_count": 0,
                    "identical_to_baseline": True,
                },
                "review_budgets": {
                    "baseline": {"top_250": {}},
                    "v04_guard": {"top_250": {}},
                },
            }
        },
    }
    path = tmp_path / "evaluation.json"
    path.write_text(json.dumps(fixture), encoding="utf-8")
    loaded = MODULE.load_json(path)
    identity = loaded["datasets"]["toy"]["conflict_winner_identity"]
    assert identity["mismatch_count"] == 0
    assert identity["identical_to_baseline"] is True
