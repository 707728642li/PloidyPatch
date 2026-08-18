from __future__ import annotations

import json
from pathlib import Path

import pytest

from ploidypatch.bootstrap import (
    SCORE_SCHEMA_VERSION,
    independent_event_bootstrap,
)
from ploidypatch.confusion_bootstrap import independent_confusion_bootstrap


def write_score(
    path: Path,
    prefix: str,
    values: tuple[bool, ...],
    *,
    false_positive: int,
) -> Path:
    true_positive = sum(values)
    false_negative = len(values) - true_positive
    precision = true_positive / (true_positive + false_positive)
    recall = true_positive / (true_positive + false_negative)
    f1 = 2 * true_positive / (
        2 * true_positive + false_positive + false_negative
    )
    report = {
        "schema_version": SCORE_SCHEMA_VERSION,
        "quality_gate": {"grade": "pass"},
        "strict_cds_chain": {
            "true_positive": true_positive,
            "false_positive": false_positive,
            "false_negative": false_negative,
            "precision": precision,
            "recall": recall,
            "f1": f1,
        },
        "event_details": [
            {
                "event_id": f"{prefix}-{index}",
                "event_type": "annotation_copy_collapse",
                "complete_cds_chain_recovery": value,
            }
            for index, value in enumerate(values)
        ],
    }
    path.write_text(json.dumps(report), encoding="utf-8")
    return path


def test_independent_event_bootstrap_accepts_distinct_event_ids(
    tmp_path: Path,
) -> None:
    left = write_score(tmp_path / "left.json", "left", (True, True, False), false_positive=1)
    right = write_score(tmp_path / "right.json", "right", (True, False), false_positive=2)
    report = independent_event_bootstrap(
        score_inputs=(("left", left), ("right", right)),
        output_json_path=tmp_path / "events.json",
        replicates=500,
        seed=13,
    )

    methods = {row["label"]: row for row in report["methods"]}
    assert methods["left"]["events"] == 3
    assert methods["right"]["events"] == 2
    assert report["parameters"]["paired_method_resampling"] is False
    assert report["independent_differences"][0]["observed_delta"] == pytest.approx(1 / 6)


def test_confusion_bootstrap_reports_all_metrics(tmp_path: Path) -> None:
    left = write_score(tmp_path / "left.json", "left", (True, True, False), false_positive=1)
    right = write_score(tmp_path / "right.json", "right", (True, False), false_positive=2)
    report = independent_confusion_bootstrap(
        score_inputs=(("left", left), ("right", right)),
        output_json_path=tmp_path / "confusion.json",
        replicates=500,
        seed=19,
    )

    methods = {row["label"]: row for row in report["methods"]}
    assert methods["left"]["counts"] == {
        "true_positive": 2,
        "false_positive": 1,
        "false_negative": 1,
    }
    assert set(methods["left"]["metrics"]) == {"precision", "recall", "f1"}
    assert report["parameters"]["resampling_unit"] == (
        "matched_or_unmatched_confusion_outcome"
    )


def test_confusion_bootstrap_rejects_inconsistent_report(tmp_path: Path) -> None:
    score = write_score(tmp_path / "score.json", "x", (True, False), false_positive=1)
    report = json.loads(score.read_text(encoding="utf-8"))
    report["strict_cds_chain"]["precision"] = 0.99
    score.write_text(json.dumps(report), encoding="utf-8")

    with pytest.raises(ValueError, match="Inconsistent"):
        independent_confusion_bootstrap(
            score_inputs=(("x", score),),
            output_json_path=tmp_path / "bad.json",
            replicates=100,
        )
