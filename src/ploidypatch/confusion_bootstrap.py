from __future__ import annotations

import json
import math
import random
from pathlib import Path
from typing import Any, Iterable

from . import __version__
from .bootstrap import SCORE_SCHEMA_VERSION, _file_sha256, _quantile


CONFUSION_BOOTSTRAP_SCHEMA_VERSION = "ploidypatch.confusion_bootstrap.v1"
SUPPORTED_CONFUSION_SECTIONS = frozenset(
    {"strict_cds_chain", "strict_transcript_structure"}
)
CONFUSION_METRICS = ("precision", "recall", "f1")


def _metrics(true_positive: int, false_positive: int, false_negative: int) -> dict[str, float]:
    precision_denominator = true_positive + false_positive
    recall_denominator = true_positive + false_negative
    precision = true_positive / precision_denominator if precision_denominator else 0.0
    recall = true_positive / recall_denominator if recall_denominator else 0.0
    f1_denominator = 2 * true_positive + false_positive + false_negative
    f1 = 2 * true_positive / f1_denominator if f1_denominator else 0.0
    return {"precision": precision, "recall": recall, "f1": f1}


def _read_confusion_score(
    label: str, path: str | Path, section: str
) -> tuple[dict[str, int], dict[str, float], dict[str, Any]]:
    score_path = Path(path)
    report = json.loads(score_path.read_text(encoding="utf-8"))
    if report.get("schema_version") != SCORE_SCHEMA_VERSION:
        raise ValueError(f"Unsupported score schema for {label}: {score_path}")
    if report.get("quality_gate", {}).get("grade") != "pass":
        raise ValueError(f"Score quality gate did not pass for {label}")
    section_report = report.get(section)
    if not isinstance(section_report, dict):
        raise ValueError(f"Score lacks {section} for {label}")
    counts: dict[str, int] = {}
    for field in ("true_positive", "false_positive", "false_negative"):
        value = section_report.get(field)
        if not isinstance(value, int) or isinstance(value, bool) or value < 0:
            raise ValueError(f"Invalid {section}.{field} for {label}")
        counts[field] = value
    if sum(counts.values()) == 0:
        raise ValueError(f"Empty confusion table for {label}")
    observed = _metrics(**counts)
    for metric, expected in observed.items():
        reported = section_report.get(metric)
        if not isinstance(reported, (int, float)) or isinstance(reported, bool):
            raise ValueError(f"Missing {section}.{metric} for {label}")
        if not math.isclose(float(reported), expected, rel_tol=1e-12, abs_tol=1e-12):
            raise ValueError(f"Inconsistent {section}.{metric} for {label}")
    return counts, observed, {
        "label": label,
        "file_name": score_path.name,
        "sha256": _file_sha256(score_path),
        "section": section,
        **counts,
    }


def independent_confusion_bootstrap(
    *,
    score_inputs: Iterable[tuple[str, str | Path]],
    output_json_path: str | Path,
    section: str = "strict_cds_chain",
    replicates: int = 10_000,
    seed: int = 20260807,
    alpha: float = 0.05,
) -> dict[str, Any]:
    """Bootstrap independent TP/FP/FN match outcomes and metric contrasts."""

    if section not in SUPPORTED_CONFUSION_SECTIONS:
        raise ValueError(f"Unsupported confusion section: {section}")
    if replicates < 100:
        raise ValueError("At least 100 bootstrap replicates are required")
    if not 0 < alpha < 1:
        raise ValueError("alpha must be between zero and one")
    output = Path(output_json_path)
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite bootstrap output: {output}")
    inputs = list(score_inputs)
    if not inputs:
        raise ValueError("At least one score input is required")
    labels = [label for label, _ in inputs]
    if any(not label for label in labels) or len(set(labels)) != len(labels):
        raise ValueError("Bootstrap labels must be non-empty and unique")

    counts_by_label: dict[str, dict[str, int]] = {}
    observed_by_label: dict[str, dict[str, float]] = {}
    input_manifests: list[dict[str, Any]] = []
    for label, path in inputs:
        counts, observed, manifest = _read_confusion_score(label, path, section)
        counts_by_label[label] = counts
        observed_by_label[label] = observed
        input_manifests.append(manifest)

    rng = random.Random(seed)
    method_samples = {
        label: {metric: [] for metric in CONFUSION_METRICS} for label in labels
    }
    pairs = [
        (left, right)
        for left_index, left in enumerate(labels)
        for right in labels[left_index + 1 :]
    ]
    delta_samples = {
        pair: {metric: [] for metric in CONFUSION_METRICS} for pair in pairs
    }
    for _ in range(replicates):
        replicate_metrics: dict[str, dict[str, float]] = {}
        for label in labels:
            counts = counts_by_label[label]
            true_positive = counts["true_positive"]
            false_positive = counts["false_positive"]
            total = sum(counts.values())
            outcomes = rng.choices(
                (0, 1, 2),
                cum_weights=(true_positive, true_positive + false_positive, total),
                k=total,
            )
            sampled_metrics = _metrics(
                true_positive=outcomes.count(0),
                false_positive=outcomes.count(1),
                false_negative=outcomes.count(2),
            )
            replicate_metrics[label] = sampled_metrics
            for metric in CONFUSION_METRICS:
                method_samples[label][metric].append(sampled_metrics[metric])
        for pair in pairs:
            for metric in CONFUSION_METRICS:
                delta_samples[pair][metric].append(
                    replicate_metrics[pair[0]][metric]
                    - replicate_metrics[pair[1]][metric]
                )

    lower_probability = alpha / 2
    upper_probability = 1 - alpha / 2
    method_reports: list[dict[str, Any]] = []
    for label in labels:
        intervals = {}
        for metric in CONFUSION_METRICS:
            samples = method_samples[label][metric]
            intervals[metric] = {
                "observed": observed_by_label[label][metric],
                "bootstrap_mean": sum(samples) / replicates,
                "ci_lower": _quantile(samples, lower_probability),
                "ci_upper": _quantile(samples, upper_probability),
            }
        method_reports.append(
            {"label": label, "counts": counts_by_label[label], "metrics": intervals}
        )
    difference_reports: list[dict[str, Any]] = []
    for left, right in pairs:
        metrics = {}
        for metric in CONFUSION_METRICS:
            samples = delta_samples[(left, right)][metric]
            metrics[metric] = {
                "observed_delta": observed_by_label[left][metric]
                - observed_by_label[right][metric],
                "bootstrap_mean_delta": sum(samples) / replicates,
                "ci_lower": _quantile(samples, lower_probability),
                "ci_upper": _quantile(samples, upper_probability),
                "probability_delta_gt_zero": sum(value > 0 for value in samples)
                / replicates,
                "probability_delta_eq_zero": sum(value == 0 for value in samples)
                / replicates,
            }
        difference_reports.append(
            {
                "left": left,
                "right": right,
                "delta_definition": "left_minus_right",
                "metrics": metrics,
            }
        )

    report = {
        "schema_version": CONFUSION_BOOTSTRAP_SCHEMA_VERSION,
        "generator": {"name": "PloidyPatch", "version": __version__},
        "section": section,
        "parameters": {
            "replicates": replicates,
            "seed": seed,
            "alpha": alpha,
            "interval": "percentile",
            "resampling_unit": "matched_or_unmatched_confusion_outcome",
            "paired_method_resampling": False,
            "independent_score_resampling": True,
        },
        "inputs": input_manifests,
        "methods": method_reports,
        "independent_differences": difference_reports,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x", encoding="utf-8", newline="") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return report
