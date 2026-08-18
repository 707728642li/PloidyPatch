from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

np = pytest.importorskip("numpy")
pytest.importorskip("scipy")
pytest.importorskip("sklearn")


SCRIPT = Path(__file__).parents[1] / "scripts" / "evaluate_apple_external_v0.3.py"
SPEC = importlib.util.spec_from_file_location("evaluate_apple_external_v0_3", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def test_event_bootstrap_is_paired_and_reproducible() -> None:
    primary = {
        "event_details": [
            {"event_id": f"e{index}", "complete_cds_chain_recovery": value}
            for index, value in enumerate((1, 1, 1, 0, 0))
        ]
    }
    legacy = {
        "event_details": [
            {"event_id": f"e{index}", "complete_cds_chain_recovery": value}
            for index, value in enumerate((1, 0, 0, 0, 0))
        ]
    }
    first, values = MODULE.event_bootstrap_delta(
        primary, legacy, replicates=1000, seed=17
    )
    second, repeated = MODULE.event_bootstrap_delta(
        primary, legacy, replicates=1000, seed=17
    )
    assert first["observed_delta"] == pytest.approx(0.4)
    assert first == second
    assert np.array_equal(values, repeated)


def test_chromosome_bootstrap_and_safety_metrics() -> None:
    labels = np.asarray([1, 0, 1, 0, 1, 0], dtype=np.uint8)
    primary = np.asarray([6, 1, 5, 2, 4, 3], dtype=float)
    baseline = np.asarray([1, 6, 2, 5, 3, 4], dtype=float)
    groups = np.asarray(["c1", "c1", "c2", "c2", "c3", "c3"])
    report, deltas = MODULE.chromosome_bootstrap_delta(
        labels, primary, baseline, groups, replicates=1000, seed=19
    )
    assert report["observed_delta"] > 0
    assert report["replicates_valid"] == 1000
    assert len(deltas) == 1000

    digests = [f"d{index}" for index in range(6)]
    review = MODULE.review_metrics(labels, primary, digests)
    assert review["top_1pct"]["true_positive"] == 1
    conflict = MODULE.conflict_metrics(
        labels, primary, digests, ["x", "x", "y", "y", "z", "z"]
    )
    assert conflict["evaluable_exactly_one_positive"] == 3
    assert conflict["top1_accuracy"] == 1.0


def test_named_secondary_paths_and_pool_summary(tmp_path: Path) -> None:
    first = tmp_path / "first.json"
    second = tmp_path / "second.json"
    parsed = MODULE.parse_named_paths(
        [f"miniprot={first}", f"support2={second}"]
    )
    assert parsed == {"miniprot": first, "support2": second}
    with pytest.raises(ValueError, match="unique NAME=PATH"):
        MODULE.parse_named_paths([f"same={first}", f"same={second}"])

    score = {
        "event_recovery": {
            "events": 800,
            "complete_cds_chain_recovery": 640,
            "complete_cds_chain_recall": 0.8,
        },
        "strict_cds_chain": {"precision": 0.7, "recall": 0.8, "f1": 0.7467},
        "background_subtraction": {"differential_candidate_cds_chains": 900},
        "collateral_changes": {
            "baseline_transcript_structures_missing_from_candidate": 0
        },
        "quality_gate": {"grade": "pass"},
    }
    assert MODULE.score_collateral_gate(score)
    summary = MODULE.summarize_pool_score(score)
    assert summary["complete_cds_chain_recovery"] == 640
    assert summary["differential_candidate_cds_chains"] == 900
    assert summary["quality_grade"] == "pass"
