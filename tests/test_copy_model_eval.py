from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import pytest

from ploidypatch.copy_features import COPY_LABEL_SCHEMA_VERSION
from ploidypatch.copy_model import COPY_SCORE_SCHEMA_VERSION
from ploidypatch.copy_model_eval import evaluate_copy_candidate_scores


def _write_table(path: Path, fields: tuple[str, ...], rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_copy_ranking_evaluation_uses_tie_correct_metrics_and_frozen_decisions(
    tmp_path: Path,
) -> None:
    scores = tmp_path / "scores.tsv"
    labels = tmp_path / "labels.tsv"
    score_rows = [
        {"candidate_digest": "a", "seqid": "1", "support_methods": "m", "wgd_existing_partner": "1", "model_calibrated_probability": "0.9", "model_review_decision": "1", "model_high_confidence_decision": "0"},
        {"candidate_digest": "b", "seqid": "1", "support_methods": "m", "wgd_existing_partner": "1", "model_calibrated_probability": "0.8", "model_review_decision": "1", "model_high_confidence_decision": "0"},
        {"candidate_digest": "c", "seqid": "2", "support_methods": "n", "wgd_existing_partner": "0", "model_calibrated_probability": "0.8", "model_review_decision": "0", "model_high_confidence_decision": "0"},
        {"candidate_digest": "d", "seqid": "2", "support_methods": "n", "wgd_existing_partner": "0", "model_calibrated_probability": "0.1", "model_review_decision": "0", "model_high_confidence_decision": "0"},
    ]
    label_rows = [
        {"candidate_digest": "a", "label_exact_cds": "1"},
        {"candidate_digest": "b", "label_exact_cds": "0"},
        {"candidate_digest": "c", "label_exact_cds": "1"},
        {"candidate_digest": "d", "label_exact_cds": "0"},
    ]
    _write_table(scores, tuple(score_rows[0]), score_rows)
    _write_table(labels, tuple(label_rows[0]), label_rows)
    Path(str(scores) + ".manifest.json").write_text(
        json.dumps({"schema_version": COPY_SCORE_SCHEMA_VERSION, "truth_access": False, "outputs": {"scores": {"sha256": _sha(scores)}}, "thresholds": {"review": {"value": 0.5}}}),
        encoding="utf-8",
    )
    Path(str(labels) + ".manifest.json").write_text(
        json.dumps({"schema_version": COPY_LABEL_SCHEMA_VERSION, "evaluator_only": True, "outputs": {"labels": {"sha256": _sha(labels)}}}),
        encoding="utf-8",
    )
    report = evaluate_copy_candidate_scores(
        scored_tsv_path=scores,
        labeled_feature_tsv_path=labels,
        output_json_path=tmp_path / "report.json",
    )
    assert report["probability_metrics"]["roc_auc"] == pytest.approx(0.875)
    assert report["probability_metrics"]["average_precision"] == pytest.approx(5 / 6)
    assert report["frozen_policies"]["review"]["true_positive"] == 1
    assert report["frozen_policies"]["review"]["false_positive"] == 1
    assert report["counts"] == {"candidates": 4, "positive_exact_cds": 2, "negative_candidates": 2}
