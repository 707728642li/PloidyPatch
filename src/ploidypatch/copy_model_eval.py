from __future__ import annotations

import csv
import json
import math
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

from . import __version__
from .baseline import _file_sha256
from .copy_features import COPY_LABEL_SCHEMA_VERSION
from .copy_model import COPY_SCORE_SCHEMA_VERSION


COPY_RANKING_EVALUATION_SCHEMA_VERSION = "ploidypatch.copy_ranking_evaluation.v1"


def _read_rows(path: Path, required: set[str]) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or not required <= set(reader.fieldnames):
            raise ValueError(f"{path} lacks required copy-ranking columns")
        rows = list(reader)
    digests = [row["candidate_digest"] for row in rows]
    if any(not digest for digest in digests) or len(digests) != len(set(digests)):
        raise ValueError(f"{path} has empty or duplicate candidate digests")
    return rows


def _load_manifest(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"Manifest is not a JSON object: {path}")
    return value


def _parse_probability(value: str) -> float:
    try:
        probability = float(value)
    except ValueError as exc:
        raise ValueError(f"Invalid copy-model probability: {value!r}") from exc
    if not math.isfinite(probability) or not 0 <= probability <= 1:
        raise ValueError(f"Copy-model probability lies outside [0, 1]: {value!r}")
    return probability


def _confusion(labels: Sequence[int], decisions: Sequence[int]) -> dict[str, Any]:
    tp = sum(label == 1 and decision == 1 for label, decision in zip(labels, decisions, strict=True))
    fp = sum(label == 0 and decision == 1 for label, decision in zip(labels, decisions, strict=True))
    fn = sum(label == 1 and decision == 0 for label, decision in zip(labels, decisions, strict=True))
    tn = sum(label == 0 and decision == 0 for label, decision in zip(labels, decisions, strict=True))
    precision = tp / (tp + fp) if tp + fp else 0.0
    recall = tp / (tp + fn) if tp + fn else 0.0
    f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
    return {
        "true_positive": tp,
        "false_positive": fp,
        "false_negative": fn,
        "true_negative": tn,
        "precision": precision,
        "recall": recall,
        "f1": f1,
        "selected": tp + fp,
    }


def _average_precision(labels: Sequence[int], probabilities: Sequence[float]) -> float:
    total_positive = sum(labels)
    if total_positive == 0:
        raise ValueError("Average precision requires at least one positive")
    grouped: dict[float, list[int]] = defaultdict(list)
    for label, probability in zip(labels, probabilities, strict=True):
        grouped[probability].append(label)
    cumulative_positive = 0
    cumulative_rows = 0
    average_precision = 0.0
    for probability in sorted(grouped, reverse=True):
        group_labels = grouped[probability]
        group_positive = sum(group_labels)
        cumulative_positive += group_positive
        cumulative_rows += len(group_labels)
        average_precision += (
            group_positive / total_positive * cumulative_positive / cumulative_rows
        )
    return average_precision


def _roc_auc(labels: Sequence[int], probabilities: Sequence[float]) -> float:
    positives = sum(labels)
    negatives = len(labels) - positives
    if positives == 0 or negatives == 0:
        raise ValueError("ROC AUC requires both label classes")
    grouped: dict[float, list[int]] = defaultdict(list)
    for label, probability in zip(labels, probabilities, strict=True):
        grouped[probability].append(label)
    negatives_below = 0
    concordant = 0.0
    for probability in sorted(grouped):
        group = grouped[probability]
        group_positive = sum(group)
        group_negative = len(group) - group_positive
        concordant += group_positive * (negatives_below + 0.5 * group_negative)
        negatives_below += group_negative
    return concordant / (positives * negatives)


def _calibration_bins(
    labels: Sequence[int],
    probabilities: Sequence[float],
    *,
    bins: int,
    strategy: str,
) -> dict[str, Any]:
    records: list[dict[str, Any]] = []
    if strategy == "equal_frequency":
        ordered = sorted(range(len(labels)), key=lambda index: probabilities[index])
        partitions = [
            ordered[start * len(ordered) // bins : (start + 1) * len(ordered) // bins]
            for start in range(bins)
        ]
    elif strategy == "equal_width":
        partitions = [[] for _ in range(bins)]
        for index, probability in enumerate(probabilities):
            bin_index = min(int(probability * bins), bins - 1)
            partitions[bin_index].append(index)
    else:
        raise ValueError(f"Unknown calibration-bin strategy: {strategy}")
    ece = 0.0
    for bin_index, indices in enumerate(partitions):
        if not indices:
            continue
        mean_probability = sum(probabilities[index] for index in indices) / len(indices)
        observed = sum(labels[index] for index in indices) / len(indices)
        gap = abs(mean_probability - observed)
        ece += len(indices) / len(labels) * gap
        records.append(
            {
                "bin": bin_index,
                "rows": len(indices),
                "minimum_probability": min(probabilities[index] for index in indices),
                "maximum_probability": max(probabilities[index] for index in indices),
                "mean_probability": mean_probability,
                "observed_fraction": observed,
                "absolute_gap": gap,
            }
        )
    return {"strategy": strategy, "bins_requested": bins, "ece": ece, "bins": records}


def _probability_metrics(labels: Sequence[int], probabilities: Sequence[float]) -> dict[str, Any]:
    epsilon = 1e-15
    brier = sum((probability - label) ** 2 for label, probability in zip(labels, probabilities, strict=True)) / len(labels)
    cross_entropy = -sum(
        label * math.log(min(max(probability, epsilon), 1 - epsilon))
        + (1 - label) * math.log(min(max(1 - probability, epsilon), 1 - epsilon))
        for label, probability in zip(labels, probabilities, strict=True)
    ) / len(labels)
    return {
        "average_precision": _average_precision(labels, probabilities),
        "roc_auc": _roc_auc(labels, probabilities),
        "brier_score": brier,
        "log_loss": cross_entropy,
        "calibration_equal_frequency_20": _calibration_bins(
            labels, probabilities, bins=20, strategy="equal_frequency"
        ),
        "calibration_equal_width_10": _calibration_bins(
            labels, probabilities, bins=10, strategy="equal_width"
        ),
    }


def _stratified_metrics(rows: Sequence[Mapping[str, Any]], field: str) -> dict[str, Any]:
    grouped: dict[str, list[Mapping[str, Any]]] = defaultdict(list)
    for row in rows:
        grouped[str(row[field])].append(row)
    result: dict[str, Any] = {}
    for value, group in sorted(grouped.items()):
        labels = [int(row["label"]) for row in group]
        probabilities = [float(row["probability"]) for row in group]
        decisions = [int(row["review_decision"]) for row in group]
        record = {
            "rows": len(group),
            "positives": sum(labels),
            "review": _confusion(labels, decisions),
        }
        if 0 < sum(labels) < len(labels):
            record["average_precision"] = _average_precision(labels, probabilities)
            record["roc_auc"] = _roc_auc(labels, probabilities)
        result[value] = record
    return result


def evaluate_copy_candidate_scores(
    *,
    scored_tsv_path: str | Path,
    labeled_feature_tsv_path: str | Path,
    output_json_path: str | Path,
) -> dict[str, Any]:
    """Evaluate a frozen truth-blind score table against evaluator-only labels."""

    scored_path = Path(scored_tsv_path)
    score_manifest_path = Path(str(scored_path) + ".manifest.json")
    labeled_path = Path(labeled_feature_tsv_path)
    label_manifest_path = Path(str(labeled_path) + ".manifest.json")
    output_path = Path(output_json_path)
    for required in (scored_path, score_manifest_path, labeled_path, label_manifest_path):
        if not required.is_file() or required.stat().st_size == 0:
            raise ValueError(f"Missing or empty copy-ranking evaluation input: {required}")
    if output_path.exists():
        raise FileExistsError(f"Refusing to overwrite copy-ranking evaluation: {output_path}")
    score_manifest = _load_manifest(score_manifest_path)
    label_manifest = _load_manifest(label_manifest_path)
    if (
        score_manifest.get("schema_version") != COPY_SCORE_SCHEMA_VERSION
        or score_manifest.get("truth_access") is not False
        or score_manifest.get("outputs", {}).get("scores", {}).get("sha256")
        != _file_sha256(scored_path)
    ):
        raise ValueError("Copy-score manifest fails schema, truth, or checksum gate")
    if (
        label_manifest.get("schema_version") != COPY_LABEL_SCHEMA_VERSION
        or label_manifest.get("evaluator_only") is not True
        or label_manifest.get("outputs", {}).get("labels", {}).get("sha256")
        != _file_sha256(labeled_path)
    ):
        raise ValueError("Copy-label manifest fails schema, evaluator, or checksum gate")
    scores = _read_rows(
        scored_path,
        {
            "candidate_digest",
            "seqid",
            "support_methods",
            "wgd_existing_partner",
            "model_calibrated_probability",
            "model_review_decision",
            "model_high_confidence_decision",
        },
    )
    labels = _read_rows(labeled_path, {"candidate_digest", "label_exact_cds"})
    label_by_digest = {row["candidate_digest"]: row for row in labels}
    if {row["candidate_digest"] for row in scores} != set(label_by_digest):
        raise ValueError("Copy-score and label candidate universes differ")
    joined: list[dict[str, Any]] = []
    for row in scores:
        label = label_by_digest[row["candidate_digest"]]["label_exact_cds"]
        if label not in {"0", "1"}:
            raise ValueError("Copy label must be zero or one")
        review = row["model_review_decision"]
        high = row["model_high_confidence_decision"]
        if review not in {"0", "1"} or high not in {"0", "1"}:
            raise ValueError("Copy model decisions must be zero or one")
        joined.append(
            {
                "candidate_digest": row["candidate_digest"],
                "seqid": row["seqid"],
                "support_methods": row["support_methods"],
                "wgd_existing_partner": row["wgd_existing_partner"],
                "probability": _parse_probability(row["model_calibrated_probability"]),
                "review_decision": int(review),
                "high_confidence_decision": int(high),
                "label": int(label),
            }
        )
    binary_labels = [int(row["label"]) for row in joined]
    probabilities = [float(row["probability"]) for row in joined]
    review_decisions = [int(row["review_decision"]) for row in joined]
    high_decisions = [int(row["high_confidence_decision"]) for row in joined]
    report: dict[str, Any] = {
        "schema_version": COPY_RANKING_EVALUATION_SCHEMA_VERSION,
        "generator": {"name": "PloidyPatch", "version": __version__},
        "evaluation_scope": "evaluator_only_external_or_development_copy_candidate_ranking",
        "inputs": {
            "scores": {"file_name": scored_path.name, "sha256": _file_sha256(scored_path)},
            "score_manifest": {"file_name": score_manifest_path.name, "sha256": _file_sha256(score_manifest_path)},
            "labels": {"file_name": labeled_path.name, "sha256": _file_sha256(labeled_path)},
            "label_manifest": {"file_name": label_manifest_path.name, "sha256": _file_sha256(label_manifest_path)},
        },
        "counts": {
            "candidates": len(joined),
            "positive_exact_cds": sum(binary_labels),
            "negative_candidates": len(joined) - sum(binary_labels),
        },
        "probability_metrics": _probability_metrics(binary_labels, probabilities),
        "frozen_policies": {
            "review": _confusion(binary_labels, review_decisions),
            "high_confidence": _confusion(binary_labels, high_decisions),
        },
        "stratified": {
            "seqid": _stratified_metrics(joined, "seqid"),
            "support_methods": _stratified_metrics(joined, "support_methods"),
            "wgd_existing_partner": _stratified_metrics(joined, "wgd_existing_partner"),
        },
        "thresholds": score_manifest.get("thresholds"),
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("x", encoding="utf-8", newline="") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return report
