#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import platform
import shutil
import sys
from collections import Counter
from pathlib import Path
from statistics import NormalDist
from typing import Any, Iterable, Mapping, Sequence

import numpy as np
import scipy
import sklearn
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    average_precision_score,
    brier_score_loss,
    log_loss,
    roc_auc_score,
)
from sklearn.model_selection import StratifiedGroupKFold

from ploidypatch import __version__ as ploidypatch_version
from ploidypatch.copy_features import FEATURE_COLUMNS
from ploidypatch.copy_model import (
    COPY_MODEL_SCHEMA_VERSION,
    FEATURE_SETS,
    fit_copy_feature_contract,
    transform_copy_feature_row,
)


TRAINING_SCHEMA_VERSION = "ploidypatch.copy_model_training.v1"
DEFAULT_C_GRID = (0.03, 0.1, 0.3, 1.0, 3.0)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def _read_labeled_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = set(FEATURE_COLUMNS) | {"label_exact_cds"}
        if reader.fieldnames is None or not required <= set(reader.fieldnames):
            raise ValueError("Labeled copy-feature table lacks required columns")
        rows = list(reader)
    if not rows:
        raise ValueError("Labeled copy-feature table has no rows")
    digests = [row["candidate_digest"] for row in rows]
    if any(not value for value in digests) or len(digests) != len(set(digests)):
        raise ValueError("Candidate digests must be nonempty and unique")
    if any(row["label_exact_cds"] not in {"0", "1"} for row in rows):
        raise ValueError("label_exact_cds must contain only zero and one")
    groups = {row["seqid"] for row in rows}
    if "" in groups or len(groups) < 5:
        raise ValueError("At least five nonempty chromosome groups are required")
    return rows


def _matrix(
    rows: Sequence[Mapping[str, str]], contract: Mapping[str, Any]
) -> np.ndarray:
    return np.asarray(
        [transform_copy_feature_row(row, contract) for row in rows], dtype=float
    )


def _fit_estimator(
    rows: Sequence[Mapping[str, str]],
    labels: np.ndarray,
    *,
    feature_set: str,
    c_value: float,
) -> tuple[dict[str, Any], LogisticRegression]:
    contract = fit_copy_feature_contract(rows, feature_set=feature_set)
    design = _matrix(rows, contract)
    estimator = LogisticRegression(
        C=c_value,
        penalty="l2",
        solver="lbfgs",
        class_weight=None,
        max_iter=5000,
        random_state=0,
    )
    estimator.fit(design, labels)
    if estimator.n_iter_[0] >= estimator.max_iter:
        raise RuntimeError("Logistic regression did not converge")
    return contract, estimator


def _validated_splits(
    labels: np.ndarray,
    groups: np.ndarray,
    *,
    n_splits: int,
    seed: int,
) -> list[tuple[np.ndarray, np.ndarray]]:
    unique_groups = np.unique(groups)
    if n_splits > len(unique_groups):
        raise ValueError("Grouped CV fold count exceeds chromosome count")
    splitter = StratifiedGroupKFold(
        n_splits=n_splits, shuffle=True, random_state=seed
    )
    placeholder = np.zeros(len(labels), dtype=np.uint8)
    splits = list(splitter.split(placeholder, labels, groups))
    for train, validation in splits:
        if len(np.unique(groups[train])) + len(np.unique(groups[validation])) != len(
            unique_groups
        ):
            raise AssertionError("Chromosome group leakage detected")
        if set(groups[train]) & set(groups[validation]):
            raise AssertionError("Chromosome group occurs in both CV partitions")
        if len(np.unique(labels[train])) != 2 or len(np.unique(labels[validation])) != 2:
            raise ValueError("Each grouped CV partition must contain both label classes")
    return splits


def _choose_c(
    rows: Sequence[Mapping[str, str]],
    labels: np.ndarray,
    groups: np.ndarray,
    *,
    feature_set: str,
    c_grid: Sequence[float],
    n_splits: int,
    seed: int,
) -> tuple[float, list[dict[str, Any]]]:
    splits = _validated_splits(
        labels, groups, n_splits=n_splits, seed=seed
    )
    records: list[dict[str, Any]] = []
    for c_value in c_grid:
        fold_scores: list[float] = []
        for train, validation in splits:
            train_rows = [rows[index] for index in train]
            validation_rows = [rows[index] for index in validation]
            contract, estimator = _fit_estimator(
                train_rows,
                labels[train],
                feature_set=feature_set,
                c_value=c_value,
            )
            probability = estimator.predict_proba(
                _matrix(validation_rows, contract)
            )[:, 1]
            fold_scores.append(
                float(average_precision_score(labels[validation], probability))
            )
        records.append(
            {
                "C": c_value,
                "fold_average_precision": fold_scores,
                "mean_average_precision": float(np.mean(fold_scores)),
            }
        )
    selected = max(
        records,
        key=lambda record: (record["mean_average_precision"], -record["C"]),
    )
    return float(selected["C"]), records


def _fit_sigmoid(logits: np.ndarray, labels: np.ndarray) -> tuple[float, float]:
    calibrator = LogisticRegression(
        penalty=None,
        solver="lbfgs",
        max_iter=5000,
        random_state=0,
    )
    calibrator.fit(logits.reshape(-1, 1), labels)
    if calibrator.n_iter_[0] >= calibrator.max_iter:
        raise RuntimeError("Sigmoid calibration did not converge")
    return float(calibrator.coef_[0, 0]), float(calibrator.intercept_[0])


def _sigmoid(values: np.ndarray) -> np.ndarray:
    clipped = np.clip(values, -700.0, 700.0)
    return 1.0 / (1.0 + np.exp(-clipped))


def _cross_fit_sigmoid(
    logits: np.ndarray,
    labels: np.ndarray,
    fold_ids: np.ndarray,
) -> tuple[np.ndarray, list[dict[str, Any]]]:
    calibrated = np.full(len(labels), np.nan, dtype=float)
    records: list[dict[str, Any]] = []
    for fold in sorted(np.unique(fold_ids)):
        validation = fold_ids == fold
        train = ~validation
        slope, intercept = _fit_sigmoid(logits[train], labels[train])
        calibrated[validation] = _sigmoid(
            slope * logits[validation] + intercept
        )
        records.append(
            {
                "outer_fold": int(fold),
                "slope": slope,
                "intercept": intercept,
                "calibration_rows": int(np.sum(train)),
                "heldout_rows": int(np.sum(validation)),
            }
        )
    if not np.all(np.isfinite(calibrated)):
        raise AssertionError("Cross-fitted calibration left non-finite predictions")
    return calibrated, records


def _confusion(labels: np.ndarray, decisions: np.ndarray) -> dict[str, Any]:
    labels = labels.astype(bool)
    decisions = decisions.astype(bool)
    tp = int(np.sum(labels & decisions))
    fp = int(np.sum(~labels & decisions))
    fn = int(np.sum(labels & ~decisions))
    tn = int(np.sum(~labels & ~decisions))
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


def _ece(
    labels: np.ndarray, probabilities: np.ndarray, bins: int = 10
) -> dict[str, Any]:
    ordered = np.argsort(probabilities, kind="stable")
    bin_records: list[dict[str, Any]] = []
    total = len(labels)
    ece = 0.0
    for indices in np.array_split(ordered, bins):
        if not len(indices):
            continue
        observed = float(np.mean(labels[indices]))
        predicted = float(np.mean(probabilities[indices]))
        absolute_gap = abs(observed - predicted)
        ece += len(indices) / total * absolute_gap
        bin_records.append(
            {
                "rows": len(indices),
                "minimum_probability": float(np.min(probabilities[indices])),
                "maximum_probability": float(np.max(probabilities[indices])),
                "mean_probability": predicted,
                "observed_fraction": observed,
                "absolute_gap": absolute_gap,
            }
        )
    return {"equal_frequency_bins": bin_records, "ece": ece}


def _probability_metrics(
    labels: np.ndarray, probabilities: np.ndarray
) -> dict[str, Any]:
    return {
        "average_precision": float(average_precision_score(labels, probabilities)),
        "roc_auc": float(roc_auc_score(labels, probabilities)),
        "brier_score": float(brier_score_loss(labels, probabilities)),
        "log_loss": float(log_loss(labels, probabilities, labels=[0, 1])),
        "calibration": _ece(labels, probabilities),
    }


def _wilson_lower(successes: int, trials: int, confidence: float) -> float:
    if trials == 0:
        return 0.0
    z_value = NormalDist().inv_cdf(confidence)
    proportion = successes / trials
    denominator = 1.0 + z_value * z_value / trials
    centre = proportion + z_value * z_value / (2.0 * trials)
    margin = z_value * math.sqrt(
        proportion * (1.0 - proportion) / trials
        + z_value * z_value / (4.0 * trials * trials)
    )
    return max(0.0, (centre - margin) / denominator)


def _threshold_candidates(
    labels: np.ndarray, probabilities: np.ndarray
) -> Iterable[dict[str, Any]]:
    order = np.argsort(-probabilities, kind="stable")
    sorted_probability = probabilities[order]
    sorted_labels = labels[order]
    cumulative_tp = np.cumsum(sorted_labels)
    total_positive = int(np.sum(labels))
    for index in range(len(order)):
        if index + 1 < len(order) and sorted_probability[index + 1] == sorted_probability[index]:
            continue
        selected = index + 1
        tp = int(cumulative_tp[index])
        fp = selected - tp
        fn = total_positive - tp
        precision = tp / selected
        recall = tp / total_positive if total_positive else 0.0
        f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
        yield {
            "value": float(sorted_probability[index]),
            "true_positive": tp,
            "false_positive": fp,
            "false_negative": fn,
            "selected": selected,
            "precision": precision,
            "recall": recall,
            "f1": f1,
        }


def _select_review_threshold(
    labels: np.ndarray, probabilities: np.ndarray
) -> dict[str, Any]:
    candidates = list(_threshold_candidates(labels, probabilities))
    if not candidates:
        raise ValueError("Cannot select a review threshold without candidates")
    selected = max(
        candidates,
        key=lambda row: (row["f1"], row["precision"], row["value"]),
    )
    return {**selected, "selection_rule": "maximize_F1_then_precision_then_threshold"}


def _select_high_confidence_threshold(
    labels: np.ndarray,
    probabilities: np.ndarray,
    *,
    minimum_precision_lower_bound: float,
    confidence: float,
    minimum_selected: int,
) -> dict[str, Any]:
    eligible: list[dict[str, Any]] = []
    for candidate in _threshold_candidates(labels, probabilities):
        lower = _wilson_lower(
            candidate["true_positive"], candidate["selected"], confidence
        )
        candidate = {**candidate, "precision_wilson_lower": lower}
        if (
            candidate["selected"] >= minimum_selected
            and lower >= minimum_precision_lower_bound
        ):
            eligible.append(candidate)
    policy = {
        "selection_rule": "maximize_recall_subject_to_one_sided_Wilson_precision_bound",
        "minimum_precision_lower_bound": minimum_precision_lower_bound,
        "one_sided_confidence": confidence,
        "minimum_selected": minimum_selected,
    }
    if not eligible:
        return {**policy, "value": None, "available": False}
    selected = max(
        eligible,
        key=lambda row: (
            row["recall"],
            row["precision_wilson_lower"],
            row["precision"],
            row["value"],
        ),
    )
    return {**policy, **selected, "available": True}


def _cross_fitted_threshold_policy(
    labels: np.ndarray,
    probabilities: np.ndarray,
    fold_ids: np.ndarray,
    *,
    high_precision: float,
    high_confidence: float,
    high_minimum_selected: int,
) -> dict[str, Any]:
    review_decisions = np.zeros(len(labels), dtype=bool)
    high_decisions = np.zeros(len(labels), dtype=bool)
    folds: list[dict[str, Any]] = []
    high_available = 0
    for fold in sorted(np.unique(fold_ids)):
        validation = fold_ids == fold
        train = ~validation
        review = _select_review_threshold(labels[train], probabilities[train])
        high = _select_high_confidence_threshold(
            labels[train],
            probabilities[train],
            minimum_precision_lower_bound=high_precision,
            confidence=high_confidence,
            minimum_selected=high_minimum_selected,
        )
        review_decisions[validation] = probabilities[validation] >= review["value"]
        if high["value"] is not None:
            high_available += 1
            high_decisions[validation] = probabilities[validation] >= high["value"]
        folds.append(
            {
                "heldout_fold": int(fold),
                "review_threshold": review["value"],
                "high_confidence_threshold": high["value"],
                "heldout_rows": int(np.sum(validation)),
            }
        )
    return {
        "leakage_control": "threshold_for_each_fold_selected_only_from_other_outer_folds",
        "folds": folds,
        "review": _confusion(labels, review_decisions),
        "high_confidence": {
            **_confusion(labels, high_decisions),
            "available_training_folds": high_available,
            "total_training_folds": len(folds),
        },
    }


def _run_nested_oof(
    rows: Sequence[Mapping[str, str]],
    labels: np.ndarray,
    groups: np.ndarray,
    outer_splits: Sequence[tuple[np.ndarray, np.ndarray]],
    *,
    feature_set: str,
    c_grid: Sequence[float],
    inner_folds: int,
    seed: int,
) -> dict[str, Any]:
    logits = np.full(len(rows), np.nan, dtype=float)
    raw_probability = np.full(len(rows), np.nan, dtype=float)
    fold_ids = np.full(len(rows), -1, dtype=int)
    fold_records: list[dict[str, Any]] = []
    for fold, (train, validation) in enumerate(outer_splits):
        train_rows = [rows[index] for index in train]
        validation_rows = [rows[index] for index in validation]
        c_value, tuning = _choose_c(
            train_rows,
            labels[train],
            groups[train],
            feature_set=feature_set,
            c_grid=c_grid,
            n_splits=inner_folds,
            seed=seed + 1000 + fold,
        )
        contract, estimator = _fit_estimator(
            train_rows,
            labels[train],
            feature_set=feature_set,
            c_value=c_value,
        )
        design = _matrix(validation_rows, contract)
        logits[validation] = estimator.decision_function(design)
        raw_probability[validation] = estimator.predict_proba(design)[:, 1]
        fold_ids[validation] = fold
        fold_records.append(
            {
                "outer_fold": fold,
                "training_chromosomes": sorted(set(groups[train])),
                "heldout_chromosomes": sorted(set(groups[validation])),
                "training_rows": len(train),
                "heldout_rows": len(validation),
                "training_positives": int(np.sum(labels[train])),
                "heldout_positives": int(np.sum(labels[validation])),
                "selected_C": c_value,
                "inner_tuning": tuning,
            }
        )
    if not np.all(np.isfinite(logits)) or np.any(fold_ids < 0):
        raise AssertionError("Nested grouped CV did not produce one OOF prediction per row")
    calibrated, calibration_records = _cross_fit_sigmoid(logits, labels, fold_ids)
    return {
        "logits": logits,
        "raw_probability": raw_probability,
        "calibrated_probability": calibrated,
        "fold_ids": fold_ids,
        "fold_records": fold_records,
        "calibration_records": calibration_records,
    }


def _group_bootstrap(
    labels: np.ndarray,
    probabilities: np.ndarray,
    groups: np.ndarray,
    *,
    threshold: float,
    replicates: int,
    seed: int,
) -> dict[str, Any]:
    unique_groups = sorted(np.unique(groups))
    group_indices = {group: np.flatnonzero(groups == group) for group in unique_groups}
    rng = np.random.default_rng(seed)
    metrics = {name: [] for name in ("average_precision", "brier_score", "precision", "recall", "f1")}
    for _ in range(replicates):
        sampled = rng.choice(unique_groups, size=len(unique_groups), replace=True)
        indices = np.concatenate([group_indices[group] for group in sampled])
        y_value = labels[indices]
        p_value = probabilities[indices]
        if len(np.unique(y_value)) != 2:
            continue
        confusion = _confusion(y_value, p_value >= threshold)
        metrics["average_precision"].append(
            float(average_precision_score(y_value, p_value))
        )
        metrics["brier_score"].append(float(brier_score_loss(y_value, p_value)))
        for name in ("precision", "recall", "f1"):
            metrics[name].append(confusion[name])
    if min(len(values) for values in metrics.values()) < replicates * 0.99:
        raise RuntimeError("Too many invalid chromosome-bootstrap replicates")
    summary: dict[str, Any] = {}
    for name, values in metrics.items():
        array = np.asarray(values)
        summary[name] = {
            "mean": float(np.mean(array)),
            "median": float(np.median(array)),
            "ci_95_percentile": [
                float(np.quantile(array, 0.025)),
                float(np.quantile(array, 0.975)),
            ],
        }
    return {
        "resampling_unit": "chromosome",
        "replicates_requested": replicates,
        "replicates_valid": len(metrics["f1"]),
        "seed": seed,
        "fixed_threshold": threshold,
        "metrics": summary,
    }


def _write_json(path: Path, value: Any) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")


def train(args: argparse.Namespace) -> None:
    input_path = Path(args.labeled_features).resolve()
    output_path = Path(args.output_dir).resolve()
    partial_path = output_path.with_name(output_path.name + ".partial")
    if output_path.exists() or partial_path.exists():
        raise FileExistsError("Refusing to overwrite model output or partial directory")
    if not input_path.is_file() or input_path.stat().st_size == 0:
        raise ValueError("Missing or empty labeled feature input")
    rows = _read_labeled_rows(input_path)
    labels = np.asarray([int(row["label_exact_cds"]) for row in rows], dtype=np.uint8)
    groups = np.asarray([row["seqid"] for row in rows], dtype=object)
    if int(np.sum(labels)) < 100:
        raise ValueError("Fewer than 100 positive candidates is insufficient for this protocol")
    outer_splits = _validated_splits(
        labels, groups, n_splits=args.outer_folds, seed=args.seed
    )
    fold_by_row = np.full(len(rows), -1, dtype=int)
    for fold, (_, validation) in enumerate(outer_splits):
        fold_by_row[validation] = fold
    if np.any(fold_by_row < 0):
        raise AssertionError("Outer fold assignment is incomplete")

    partial_path.mkdir(parents=True)
    try:
        feature_results: dict[str, Any] = {}
        oof_results: dict[str, dict[str, Any]] = {}
        for feature_set in FEATURE_SETS:
            result = _run_nested_oof(
                rows,
                labels,
                groups,
                outer_splits,
                feature_set=feature_set,
                c_grid=args.c_grid,
                inner_folds=args.inner_folds,
                seed=args.seed,
            )
            probabilities = result["calibrated_probability"]
            review = _select_review_threshold(labels, probabilities)
            high = _select_high_confidence_threshold(
                labels,
                probabilities,
                minimum_precision_lower_bound=args.high_precision_lower,
                confidence=args.high_confidence,
                minimum_selected=args.high_minimum_selected,
            )
            feature_results[feature_set] = {
                "nested_grouped_oof_metrics": _probability_metrics(
                    labels, probabilities
                ),
                "raw_probability_metrics": _probability_metrics(
                    labels, result["raw_probability"]
                ),
                "descriptive_global_oof_thresholds": {
                    "review": review,
                    "high_confidence": high,
                    "warning": "thresholds selected on all OOF labels; use cross_fitted_policy for leakage-controlled performance",
                },
                "cross_fitted_policy": _cross_fitted_threshold_policy(
                    labels,
                    probabilities,
                    result["fold_ids"],
                    high_precision=args.high_precision_lower,
                    high_confidence=args.high_confidence,
                    high_minimum_selected=args.high_minimum_selected,
                ),
                "outer_folds": result["fold_records"],
                "cross_fitted_calibration": result["calibration_records"],
            }
            oof_results[feature_set] = result

        primary = oof_results["full"]
        primary_probabilities = primary["calibrated_probability"]
        review_threshold = feature_results["full"][
            "descriptive_global_oof_thresholds"
        ]["review"]
        high_threshold = feature_results["full"][
            "descriptive_global_oof_thresholds"
        ]["high_confidence"]
        bootstrap = _group_bootstrap(
            labels,
            primary_probabilities,
            groups,
            threshold=review_threshold["value"],
            replicates=args.bootstrap_replicates,
            seed=args.seed + 50000,
        )

        final_c, final_tuning = _choose_c(
            rows,
            labels,
            groups,
            feature_set="full",
            c_grid=args.c_grid,
            n_splits=args.outer_folds,
            seed=args.seed + 90000,
        )
        final_contract, final_estimator = _fit_estimator(
            rows,
            labels,
            feature_set="full",
            c_value=final_c,
        )
        final_slope, final_calibration_intercept = _fit_sigmoid(
            primary["logits"], labels
        )
        model = {
            "schema_version": COPY_MODEL_SCHEMA_VERSION,
            "generator": {
                "name": "PloidyPatch",
                "version": ploidypatch_version,
                "training_schema_version": TRAINING_SCHEMA_VERSION,
            },
            "development_scope": {
                "species": "Glycine max",
                "task": "annotation_copy_collapse_candidate_ranking",
                "candidate_universe": "exact_phased_CDS_union_of_miniprot_GeMoMa_LiftOn",
                "truth_use": "evaluator_only_during_model_development",
            },
            "training": {
                "labeled_features_file": input_path.name,
                "labeled_features_sha256": _sha256(input_path),
                "rows": len(rows),
                "positive_exact_cds": int(np.sum(labels)),
                "negative_candidates": int(len(labels) - np.sum(labels)),
                "chromosomes": sorted(np.unique(groups)),
                "chromosome_count": len(np.unique(groups)),
                "protocol": {
                    "outer_split": "StratifiedGroupKFold_by_seqid",
                    "outer_folds": args.outer_folds,
                    "inner_split": "StratifiedGroupKFold_by_seqid",
                    "inner_folds": args.inner_folds,
                    "hyperparameter_metric": "average_precision",
                    "C_grid": list(args.c_grid),
                    "random_seed": args.seed,
                    "class_weight": None,
                },
            },
            "feature_contract": final_contract,
            "estimator": {
                "family": "logistic_regression_l2",
                "solver": "lbfgs",
                "C": final_c,
                "intercept": float(final_estimator.intercept_[0]),
                "coefficients": [float(value) for value in final_estimator.coef_[0]],
                "coefficient_feature_order": final_contract[
                    "expanded_feature_names"
                ],
            },
            "calibration": {
                "method": "sigmoid_on_raw_logit",
                "slope": final_slope,
                "intercept": final_calibration_intercept,
                "fit_source": "all_nested_grouped_out_of_fold_logits",
            },
            "thresholds": {
                "review": review_threshold,
                "high_confidence": high_threshold,
            },
            "claim_boundary": {
                "automatic_approval": False,
                "output_role": "ranked_review_candidates",
                "copy_addition_requires": "human_or_orthogonal_evidence_review_before_patch_compilation",
                "external_holdout_status": "not_evaluated_when_model_was_frozen",
            },
        }
        report = {
            "schema_version": TRAINING_SCHEMA_VERSION,
            "input": {
                "path": str(input_path),
                "sha256": _sha256(input_path),
                "rows": len(rows),
                "positive_exact_cds": int(np.sum(labels)),
                "negative_candidates": int(len(labels) - np.sum(labels)),
                "chromosome_counts": dict(sorted(Counter(groups).items())),
                "chromosome_positive_counts": dict(
                    sorted(Counter(groups[labels == 1]).items())
                ),
            },
            "protocol": {
                "leakage_unit": "seqid_chromosome",
                "outer_folds": args.outer_folds,
                "inner_folds": args.inner_folds,
                "C_grid": list(args.c_grid),
                "seed": args.seed,
                "primary_feature_set": "full",
                "ablation_feature_sets": list(FEATURE_SETS[1:]),
                "high_confidence_precision_lower_bound": args.high_precision_lower,
                "high_confidence_one_sided_confidence": args.high_confidence,
                "high_confidence_minimum_selected": args.high_minimum_selected,
            },
            "feature_set_results": feature_results,
            "chromosome_bootstrap_primary_review": bootstrap,
            "final_fit": {
                "selected_C": final_c,
                "full_development_tuning": final_tuning,
                "calibration_slope": final_slope,
                "calibration_intercept": final_calibration_intercept,
            },
            "interpretation_guardrails": [
                "The global OOF threshold is a frozen deployment rule, not an unbiased estimate of thresholded performance.",
                "Cross-fitted threshold-policy metrics are the leakage-controlled development estimate.",
                "The chromosome bootstrap quantifies chromosome sampling variation, not species transfer uncertainty.",
                "Only a zero-retuning external species holdout can establish portability.",
            ],
        }

        _write_json(partial_path / "model.json", model)
        _write_json(partial_path / "training_report.json", report)
        environment = {
            "python": sys.version,
            "platform": platform.platform(),
            "ploidypatch": ploidypatch_version,
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "scikit_learn": sklearn.__version__,
            "executable": sys.executable,
            "script": str(Path(__file__).resolve()),
            "script_sha256": _sha256(Path(__file__).resolve()),
        }
        _write_json(partial_path / "environment.json", environment)

        with (partial_path / "fold_assignments.tsv").open(
            "x", encoding="utf-8", newline=""
        ) as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=("seqid", "outer_fold", "rows", "positive_exact_cds"),
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            for group in sorted(np.unique(groups)):
                indices = groups == group
                folds = set(fold_by_row[indices])
                if len(folds) != 1:
                    raise AssertionError("Chromosome was split across outer folds")
                writer.writerow(
                    {
                        "seqid": group,
                        "outer_fold": next(iter(folds)),
                        "rows": int(np.sum(indices)),
                        "positive_exact_cds": int(np.sum(labels[indices])),
                    }
                )

        with (partial_path / "oof_predictions.tsv").open(
            "x", encoding="utf-8", newline=""
        ) as handle:
            fields = (
                "candidate_digest",
                "seqid",
                "support_methods",
                "wgd_existing_partner",
                "label_exact_cds",
                "outer_fold",
                "full_raw_logit",
                "full_raw_probability",
                "full_calibrated_probability",
                "full_review_decision_global_threshold",
                "no_wgd_calibrated_probability",
                "no_method_quality_calibrated_probability",
                "method_support_only_calibrated_probability",
            )
            writer = csv.DictWriter(
                handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            for index, row in enumerate(rows):
                writer.writerow(
                    {
                        "candidate_digest": row["candidate_digest"],
                        "seqid": row["seqid"],
                        "support_methods": row["support_methods"],
                        "wgd_existing_partner": row["wgd_existing_partner"],
                        "label_exact_cds": row["label_exact_cds"],
                        "outer_fold": int(primary["fold_ids"][index]),
                        "full_raw_logit": format(primary["logits"][index], ".17g"),
                        "full_raw_probability": format(
                            primary["raw_probability"][index], ".17g"
                        ),
                        "full_calibrated_probability": format(
                            primary_probabilities[index], ".17g"
                        ),
                        "full_review_decision_global_threshold": int(
                            primary_probabilities[index] >= review_threshold["value"]
                        ),
                        "no_wgd_calibrated_probability": format(
                            oof_results["no_wgd"]["calibrated_probability"][index],
                            ".17g",
                        ),
                        "no_method_quality_calibrated_probability": format(
                            oof_results["no_method_quality"][
                                "calibrated_probability"
                            ][index],
                            ".17g",
                        ),
                        "method_support_only_calibrated_probability": format(
                            oof_results["method_support_only"][
                                "calibrated_probability"
                            ][index],
                            ".17g",
                        ),
                    }
                )

        artifact_names = sorted(
            path.name for path in partial_path.iterdir() if path.is_file()
        )
        with (partial_path / "SHA256SUMS").open(
            "x", encoding="utf-8", newline=""
        ) as handle:
            for name in artifact_names:
                handle.write(f"{_sha256(partial_path / name)}  {name}\n")
        os.replace(partial_path, output_path)
    except BaseException:
        raise


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Train the frozen PloidyPatch copy-candidate ranking model"
    )
    parser.add_argument("--labeled-features", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--seed", type=int, default=20260807)
    parser.add_argument("--outer-folds", type=int, default=5)
    parser.add_argument("--inner-folds", type=int, default=4)
    parser.add_argument(
        "--c-grid", type=float, nargs="+", default=list(DEFAULT_C_GRID)
    )
    parser.add_argument("--bootstrap-replicates", type=int, default=10000)
    parser.add_argument("--high-precision-lower", type=float, default=0.90)
    parser.add_argument("--high-confidence", type=float, default=0.95)
    parser.add_argument("--high-minimum-selected", type=int, default=30)
    args = parser.parse_args()
    if not 2 <= args.outer_folds or not 2 <= args.inner_folds:
        parser.error("outer and inner fold counts must be at least two")
    if any(value <= 0 for value in args.c_grid) or len(set(args.c_grid)) != len(
        args.c_grid
    ):
        parser.error("C grid values must be positive and unique")
    if not 0 < args.high_precision_lower <= 1:
        parser.error("high-precision lower bound must be within (0, 1]")
    if not 0.5 < args.high_confidence < 1:
        parser.error("high-confidence must be within (0.5, 1)")
    if args.high_minimum_selected < 1 or args.bootstrap_replicates < 100:
        parser.error("invalid minimum-selected or bootstrap replicate count")
    train(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
