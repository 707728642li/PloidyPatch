from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np


MODULE_PATH = (
    Path(__file__).parents[1]
    / "scripts"
    / "evaluate_conflict_winner_guard_v0.4.py"
)
SPEC = importlib.util.spec_from_file_location("conflict_guard_evaluator", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_evaluator_uses_production_guard_and_exact_winner_identity() -> None:
    rows = [
        {
            "candidate_digest": "a",
            "seqid": "1",
            "label": "1",
            "conflict_set": "x",
            "baseline": "1",
            "primary": "0",
        },
        {
            "candidate_digest": "b",
            "seqid": "1",
            "label": "0",
            "conflict_set": "x",
            "baseline": "0",
            "primary": "2",
        },
    ]
    baseline, primary, guard, conflicts, guarded_sets, audit = (
        MODULE.apply_production_guard(rows)
    )
    assert baseline.tolist() == [1.0, 0.0]
    assert primary.tolist() == [0.0, 2.0]
    assert guard.tolist() == baseline.tolist()
    assert conflicts == {"x": [0, 1]}
    assert guarded_sets == {"x"}
    assert audit["winner_mismatch_count"] == 0
    assert (
        audit["baseline_winner_mapping_sha256"]
        == audit["guard_winner_mapping_sha256"]
    )


def test_review_metrics_freeze_all_six_budgets_and_selection_digest() -> None:
    labels = np.asarray([1] + [0] * 299, dtype=np.uint8)
    scores = np.arange(300, dtype=float)[::-1]
    digests = [f"d{index:03d}" for index in range(300)]
    review = MODULE.review_metrics(labels, scores, digests)
    assert list(review) == [
        "top_0.5pct",
        "top_1pct",
        "top_2pct",
        "top_100",
        "top_250",
        "top_500",
    ]
    assert review["top_250"]["reviewed"] == 250
    assert review["top_500"]["reviewed"] == 300
    assert len(str(review["top_250"]["selected_digest_sha256"])) == 64
    assert review["top_250"]["positive_candidate_recall"] == 1.0
