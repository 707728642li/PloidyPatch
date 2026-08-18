from __future__ import annotations

import numpy as np

from ploidypatch.external_evaluation import review_metrics


def test_review_metrics_reports_all_frozen_budgets_with_digest_tie_break() -> None:
    labels = np.zeros(300, dtype=np.uint8)
    labels[-1] = 1
    scores = np.zeros(300, dtype=float)
    digests = [f"d{index:03d}" for index in reversed(range(300))]

    review = review_metrics(labels, scores, digests)

    assert list(review) == [
        "top_0.5pct",
        "top_1pct",
        "top_2pct",
        "top_100",
        "top_250",
        "top_500",
    ]
    assert [entry["reviewed"] for entry in review.values()] == [2, 3, 6, 100, 250, 300]
    # Every score is tied, so the digest ordering puts d000 (the positive) first.
    assert all(entry["true_positive"] == 1 for entry in review.values())
