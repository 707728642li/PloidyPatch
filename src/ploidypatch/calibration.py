from __future__ import annotations

import csv
import hashlib
import itertools
import json
import math
from pathlib import Path
from typing import Any

from . import __version__
from .perturb import _file_sha256


CALIBRATION_SCHEMA_VERSION = "ploidypatch.synteny_tier_calibration.v1"


RULE_GRID: tuple[tuple[str, str, tuple[float, ...]], ...] = (
    ("model_rank", "max", (1, 2, 20)),
    ("model_identity", "min", (0.5, 0.8, 0.9, 0.95, 0.98)),
    ("model_query_coverage", "min", (0.5, 0.8, 0.9, 0.95, 1.0)),
    (
        "baseline_existing_cds_overlap_fraction",
        "max",
        (0.0, 0.05, 0.1, 0.2),
    ),
    ("supporting_block_count", "min", (1, 2, 3)),
    ("target_gap_genes", "max", (1, 2, 5)),
    ("locus_span_bp", "max", (25_000, 100_000, 250_000, 500_000)),
    ("best_block_pairs", "min", (2, 25, 100, 500)),
    ("best_block_pvalue", "max", (0.001, 0.01, 0.05, 0.2)),
)
PROJECTION_RULE_GRID: tuple[tuple[str, str, tuple[float, ...]], ...] = (
    ("model_identity", "min", (0.5, 0.8, 0.9, 0.95, 0.98)),
    ("model_query_coverage", "min", (0.5, 0.8, 0.9, 0.95, 1.0)),
    (
        "baseline_existing_cds_overlap_fraction",
        "max",
        (0.0, 0.05, 0.1, 0.2),
    ),
    ("support_model_count", "min", (1, 2, 3, 5)),
    ("support_query_count", "min", (1, 2, 3)),
    ("support_source_count", "min", (1, 2)),
    ("support_min_identity", "min", (0.5, 0.8, 0.9, 0.95)),
    ("support_min_query_coverage", "min", (0.5, 0.8, 0.9, 0.95)),
)


def _rate(numerator: int, denominator: int) -> float | None:
    return numerator / denominator if denominator else None


def _f1(precision: float | None, recall: float | None) -> float | None:
    if precision is None or recall is None or precision + recall == 0:
        return None
    return 2 * precision * recall / (precision + recall)


def _metrics(prediction: int, labels: int, split: int) -> dict[str, Any]:
    selected = (prediction & split).bit_count()
    positives = (labels & split).bit_count()
    true_positive = (prediction & labels & split).bit_count()
    false_positive = selected - true_positive
    false_negative = positives - true_positive
    precision = _rate(true_positive, selected)
    recall = _rate(true_positive, positives)
    return {
        "selected": selected,
        "positives": positives,
        "true_positive": true_positive,
        "false_positive": false_positive,
        "false_negative": false_negative,
        "precision": precision,
        "recall": recall,
        "f1": _f1(precision, recall),
    }


def _condition_mask(
    values: list[float], operation: str, threshold: float
) -> int:
    mask = 0
    for index, value in enumerate(values):
        keep = value >= threshold if operation == "min" else value <= threshold
        if keep:
            mask |= 1 << index
    return mask


def _rule_complexity(
    rule: dict[str, float],
    rule_grid: tuple[tuple[str, str, tuple[float, ...]], ...],
) -> int:
    broadest = {
        field: (min(values) if operation == "min" else max(values))
        for field, operation, values in rule_grid
    }
    return sum(rule[field] != broadest[field] for field in rule)


def calibrate_synteny_tiers(
    *,
    labeled_tsv_path: str | Path,
    output_json_path: str | Path,
    label_column: str,
    high_precision_floor: float = 0.1,
    high_precision_min_selected: int = 20,
    eligible_column: str | None = None,
    feature_set: str = "synteny",
) -> dict[str, Any]:
    """Fit interpretable evidence tiers on chromosome-disjoint development data."""

    if not 0 < high_precision_floor <= 1:
        raise ValueError("high_precision_floor must be in (0, 1]")
    if high_precision_min_selected < 1:
        raise ValueError("high_precision_min_selected must be positive")
    if feature_set == "synteny":
        rule_grid = RULE_GRID
    elif feature_set == "projection_support":
        rule_grid = PROJECTION_RULE_GRID
    else:
        raise ValueError(f"Unsupported calibration feature set: {feature_set}")
    output = Path(output_json_path)
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite calibration report: {output}")

    with Path(labeled_tsv_path).open(
        "r", encoding="utf-8", newline=""
    ) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = set(reader.fieldnames or [])
        if "query_seqid" in fieldnames:
            split_column = "query_seqid"
        elif "model_seqid" in fieldnames:
            split_column = "model_seqid"
        else:
            split_column = "query_seqid"
        required = {"model_id", split_column, label_column} | {
            field for field, _, _ in rule_grid
        }
        if eligible_column is not None:
            required.add(eligible_column)
        missing = required - fieldnames
        if reader.fieldnames is None or missing:
            raise ValueError(
                "Labeled TSV is missing column(s): " + ", ".join(sorted(missing))
            )
        input_rows = list(reader)
    rows = []
    for row in input_rows:
        if eligible_column is None:
            rows.append(row)
            continue
        try:
            eligible = int(row[eligible_column])
        except ValueError as exc:
            raise ValueError(
                f"Eligibility column must be binary: {row.get('model_id', '')}"
            ) from exc
        if eligible not in {0, 1}:
            raise ValueError(
                f"Eligibility column must be binary: {row.get('model_id', '')}"
            )
        if eligible:
            rows.append(row)
    if not rows:
        raise ValueError("Labeled TSV contains no rows")

    model_ids: set[str] = set()
    numeric: dict[str, list[float]] = {field: [] for field, _, _ in rule_grid}
    labels = 0
    train_mask = 0
    validation_mask = 0
    split_by_seqid: dict[str, str] = {}
    for index, row in enumerate(rows):
        model_id = row["model_id"]
        if not model_id or model_id in model_ids:
            raise ValueError(f"Empty or duplicate model ID: {model_id}")
        model_ids.add(model_id)
        try:
            label = int(row[label_column])
            values = {
                field: float(row[field]) for field, _, _ in rule_grid
            }
        except ValueError as exc:
            raise ValueError(f"Non-numeric calibration value for {model_id}") from exc
        if label not in {0, 1}:
            raise ValueError(f"Calibration label must be binary: {model_id}")
        for field, value in values.items():
            numeric[field].append(value)
        if label:
            labels |= 1 << index
        seqid = row[split_column]
        split = split_by_seqid.setdefault(
            seqid,
            "validation"
            if int(hashlib.sha256(seqid.encode("utf-8")).hexdigest()[:8], 16) % 10
            < 3
            else "train",
        )
        if split == "validation":
            validation_mask |= 1 << index
        else:
            train_mask |= 1 << index
    if not train_mask or not validation_mask:
        raise ValueError("Chromosome hash split produced an empty partition")
    if not (labels & train_mask) or not (labels & validation_mask):
        raise ValueError("Chromosome hash split produced a partition without positives")

    condition_masks: dict[tuple[str, float], int] = {}
    univariate: dict[str, list[dict[str, Any]]] = {}
    all_mask = (1 << len(rows)) - 1
    for field, operation, thresholds in rule_grid:
        summaries = []
        for threshold in thresholds:
            mask = _condition_mask(numeric[field], operation, threshold)
            condition_masks[(field, threshold)] = mask
            summaries.append(
                {
                    "operation": operation,
                    "threshold": threshold,
                    "train": _metrics(mask, labels, train_mask),
                    "validation": _metrics(mask, labels, validation_mask),
                }
            )
        univariate[field] = summaries

    baseline = {
        "train": _metrics(all_mask, labels, train_mask),
        "validation": _metrics(all_mask, labels, validation_mask),
    }
    best_balanced: tuple[tuple[float, ...], dict[str, Any]] | None = None
    best_high_precision: tuple[tuple[float, ...], dict[str, Any]] | None = None
    evaluated_rules = 0
    threshold_sets = [values for _, _, values in rule_grid]
    for thresholds in itertools.product(*threshold_sets):
        prediction = all_mask
        rule = {}
        for (field, _, _), threshold in zip(rule_grid, thresholds, strict=True):
            prediction &= condition_masks[(field, threshold)]
            rule[field] = threshold
        train = _metrics(prediction, labels, train_mask)
        if train["selected"] == 0:
            continue
        evaluated_rules += 1
        validation = _metrics(prediction, labels, validation_mask)
        result = {
            "rule": rule,
            "complexity": _rule_complexity(rule, rule_grid),
            "train": train,
            "validation": validation,
        }
        balanced_key = (
            train["f1"] or 0.0,
            train["precision"] or 0.0,
            train["recall"] or 0.0,
            -result["complexity"],
        )
        if best_balanced is None or balanced_key > best_balanced[0]:
            best_balanced = (balanced_key, result)
        if (
            train["selected"] >= high_precision_min_selected
            and (train["precision"] or 0.0) >= high_precision_floor
        ):
            high_precision_key = (
                train["recall"] or 0.0,
                train["precision"] or 0.0,
                -result["complexity"],
                train["selected"],
            )
            if (
                best_high_precision is None
                or high_precision_key > best_high_precision[0]
            ):
                best_high_precision = (high_precision_key, result)

    report: dict[str, Any] = {
        "schema_version": CALIBRATION_SCHEMA_VERSION,
        "generator": {"name": "PloidyPatch", "version": __version__},
        "access": "evaluator_only_development",
        "input": {
            "file_name": Path(labeled_tsv_path).name,
            "sha256": _file_sha256(labeled_tsv_path),
            "rows": len(rows),
            "input_rows_before_eligibility_filter": len(input_rows),
            "label_column": label_column,
            "eligible_column": eligible_column,
            "feature_set": feature_set,
        },
        "split": {
            "unit": "chromosome_seqid",
            "column": split_column,
            "algorithm": "sha256_seqid_mod10_lt3_validation",
            "assignments": dict(sorted(split_by_seqid.items())),
            "train_rows": train_mask.bit_count(),
            "validation_rows": validation_mask.bit_count(),
        },
        "grid": {
            "declared_conditions": [
                {"field": field, "operation": operation, "thresholds": list(values)}
                for field, operation, values in rule_grid
            ],
            "cartesian_rules": math.prod(len(values) for values in threshold_sets),
            "nonempty_train_rules": evaluated_rules,
        },
        "baseline": baseline,
        "selected_policies": {
            "balanced_train_f1": best_balanced[1] if best_balanced else None,
            "high_precision": (
                best_high_precision[1] if best_high_precision else None
            ),
        },
        "high_precision_constraint": {
            "minimum_train_precision": high_precision_floor,
            "minimum_train_selected": high_precision_min_selected,
        },
        "univariate": univariate,
        "warning": (
            "Rules are selected on the development chromosomes only. Validation "
            "chromosomes must not influence policy selection, and an independent "
            "species/release remains required."
        ),
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x", encoding="utf-8", newline="") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return report
