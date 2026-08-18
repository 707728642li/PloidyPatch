#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np
import sklearn
from scipy.stats import rankdata
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import StratifiedGroupKFold

from ploidypatch.copy_model import (
    fit_copy_feature_contract,
    transform_copy_feature_row,
)


SCHEMA_VERSION = "ploidypatch.homeolog_topology_evaluation.v3"
TOPOLOGY_COMPONENT_FIELDS = (
    "cds_bp_ratio",
    "cds_segment_count_ratio",
    "phase_lcs_similarity",
    "junction_fraction_similarity",
    "coding_span_ratio",
)
OOF_SEEDS = (20260807, 20260817, 20260827, 20260837, 20260847)
BOOTSTRAP_REPLICATES = 2_000
ADDON_NAMES = (
    "topology_available",
    *TOPOLOGY_COMPONENT_FIELDS,
    "wgd_block_count_within_species_percentile",
    "wgd_block_length_within_species_percentile",
)
FEATURE_ADDON_INDEXES = {
    "baseline": (),
    "normalized_wgd_context": (0, 6, 7),
    "topology": (0, 1, 2, 3, 4, 5),
    "augmented": tuple(range(8)),
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def _read_index(path: Path) -> dict[str, dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or "candidate_digest" not in reader.fieldnames:
            raise ValueError(f"Missing candidate_digest column: {path}")
        rows = list(reader)
    indexed = {row["candidate_digest"]: row for row in rows}
    if len(indexed) != len(rows) or "" in indexed:
        raise ValueError(f"Empty or duplicate candidate digest: {path}")
    return indexed


def _parse_dataset(value: str) -> tuple[str, str, Path, Path, Path]:
    name, separator, payload = value.partition("=")
    parts = payload.split(",")
    if not separator or not name or len(parts) != 4:
        raise ValueError(
            "--dataset requires NAME=ROLE,TOPOLOGY,LABELS,SCORES"
        )
    role, topology, labels, scores = parts
    if role not in {"development", "post_holdout_diagnostic"}:
        raise ValueError(f"Unknown dataset role: {role}")
    return name, role, Path(topology), Path(labels), Path(scores)


def _percentiles(values: np.ndarray, available: np.ndarray) -> np.ndarray:
    output = np.zeros(len(values), dtype=float)
    selected = np.flatnonzero(available == 1)
    if len(selected) == 1:
        output[selected] = 0.5
    elif len(selected) > 1:
        output[selected] = (rankdata(values[selected], method="average") - 1) / (
            len(selected) - 1
        )
    return output


def _load_dataset(
    name: str, role: str, topology_path: Path, label_path: Path, score_path: Path
) -> dict[str, Any]:
    for path in (topology_path, label_path, score_path):
        if not path.is_file():
            raise FileNotFoundError(path)
    topology = _read_index(topology_path)
    labels = _read_index(label_path)
    scores = _read_index(score_path)
    if set(topology) != set(labels) or set(topology) != set(scores):
        raise ValueError(f"Candidate universes differ for {name}")
    digests = sorted(topology)
    y = np.asarray([int(labels[digest]["label_exact_cds"]) for digest in digests])
    groups = np.asarray([labels[digest]["seqid"] for digest in digests])
    probability = np.asarray(
        [float(scores[digest]["model_calibrated_probability"]) for digest in digests]
    )
    available = np.asarray(
        [int(topology[digest]["topology_available"]) for digest in digests]
    )
    topology_components = np.asarray(
        [
            [
                float(topology[digest][field]) if available[index] else 0.0
                for field in TOPOLOGY_COMPONENT_FIELDS
            ]
            for index, digest in enumerate(digests)
        ]
    )
    topology_score = np.asarray(
        [
            float(topology[digest]["topology_coherence_score"])
            if available[index]
            else 0.0
            for index, digest in enumerate(digests)
        ]
    )
    block_count = np.asarray(
        [
            float(topology[digest]["wgd_support_block_count"] or 0)
            for digest in digests
        ]
    )
    block_length = np.asarray(
        [
            float(topology[digest]["wgd_longest_block_pairs"] or 0)
            for digest in digests
        ]
    )
    normalized_wgd = np.column_stack(
        (
            _percentiles(block_count, available),
            _percentiles(block_length, available),
        )
    )
    addons = np.column_stack((available, topology_components, normalized_wgd))
    return {
        "name": name,
        "role": role,
        "paths": (topology_path, label_path, score_path),
        "digests": digests,
        "labels": y,
        "groups": groups,
        "probability": probability,
        "available": available,
        "topology_score": topology_score,
        "feature_rows": [labels[digest] for digest in digests],
        "addons": addons,
    }


def _fit_feature_set(
    train_rows: list[dict[str, str]],
    train_addons: np.ndarray,
    train_y: np.ndarray,
    test_rows: list[dict[str, str]],
    test_addons: np.ndarray,
    feature_set: str,
) -> tuple[np.ndarray, list[float], list[str]]:
    contract = fit_copy_feature_contract(train_rows, feature_set="full")
    train_base = np.asarray(
        [transform_copy_feature_row(row, contract) for row in train_rows],
        dtype=float,
    )
    test_base = np.asarray(
        [transform_copy_feature_row(row, contract) for row in test_rows],
        dtype=float,
    )
    addon_indexes = FEATURE_ADDON_INDEXES[feature_set]
    if addon_indexes:
        train_x = np.column_stack(
            (train_base, train_addons[:, addon_indexes])
        )
        test_x = np.column_stack((test_base, test_addons[:, addon_indexes]))
    else:
        train_x = train_base
        test_x = test_base
    feature_names = list(contract["expanded_feature_names"]) + [
        ADDON_NAMES[index] for index in addon_indexes
    ]
    model = LogisticRegression(
        C=1.0,
        solver="lbfgs",
        max_iter=5000,
        random_state=0,
    ).fit(train_x, train_y)
    probabilities = model.predict_proba(test_x)[:, 1]
    return probabilities, model.coef_[0].tolist(), feature_names


def _metrics(labels: np.ndarray, scores: np.ndarray) -> dict[str, float]:
    return {
        "average_precision": float(average_precision_score(labels, scores)),
        "roc_auc": float(roc_auc_score(labels, scores)),
    }


def _grouped_oof(
    dataset: dict[str, Any], matrix_name: str, *, seed: int
) -> np.ndarray:
    labels = dataset["labels"]
    predictions = np.full(len(labels), np.nan)
    splitter = StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=seed)
    placeholder = np.zeros((len(labels), 1), dtype=float)
    rows = dataset["feature_rows"]
    for train, test in splitter.split(placeholder, labels, dataset["groups"]):
        train_rows = [rows[index] for index in train]
        test_rows = [rows[index] for index in test]
        predictions[test], _, _ = _fit_feature_set(
            train_rows,
            dataset["addons"][train],
            labels[train],
            test_rows,
            dataset["addons"][test],
            matrix_name,
        )
    if np.isnan(predictions).any():
        raise AssertionError("Grouped OOF predictions are incomplete")
    return predictions


def _transfer(
    train_sets: list[dict[str, Any]], test_set: dict[str, Any], matrix_name: str
) -> tuple[np.ndarray, list[float], list[str]]:
    train_rows = [
        row for dataset in train_sets for row in dataset["feature_rows"]
    ]
    train_addons = np.concatenate(
        [dataset["addons"] for dataset in train_sets]
    )
    train_y = np.concatenate([dataset["labels"] for dataset in train_sets])
    return _fit_feature_set(
        train_rows,
        train_addons,
        train_y,
        test_set["feature_rows"],
        test_set["addons"],
        matrix_name,
    )


def _group_bootstrap_delta(
    labels: np.ndarray,
    left_scores: np.ndarray,
    right_scores: np.ndarray,
    groups: np.ndarray,
    *,
    seed: int = 20260807,
) -> dict[str, Any]:
    unique_groups, group_codes = np.unique(groups, return_inverse=True)
    observed = (
        average_precision_score(labels, left_scores)
        - average_precision_score(labels, right_scores)
    )
    rng = np.random.default_rng(seed)
    deltas: list[float] = []

    def weighted_ap(scores: np.ndarray, counts: np.ndarray) -> np.ndarray:
        order = np.argsort(scores, kind="mergesort")[::-1]
        ordered_scores = scores[order]
        ordered_labels = labels[order]
        ordered_group_codes = group_codes[order]
        weights = counts[:, ordered_group_codes]
        true_weights = weights * ordered_labels
        cumulative_true = np.cumsum(true_weights, axis=1)
        cumulative_total = np.cumsum(weights, axis=1)
        threshold_indices = np.concatenate(
            (np.flatnonzero(np.diff(ordered_scores)), [len(scores) - 1])
        )
        true_at_threshold = cumulative_true[:, threshold_indices]
        total_at_threshold = cumulative_total[:, threshold_indices]
        precision = np.divide(
            true_at_threshold,
            total_at_threshold,
            out=np.zeros_like(true_at_threshold, dtype=float),
            where=total_at_threshold != 0,
        )
        total_positive = cumulative_true[:, -1]
        recall = np.divide(
            true_at_threshold,
            total_positive[:, None],
            out=np.zeros_like(true_at_threshold, dtype=float),
            where=total_positive[:, None] != 0,
        )
        recall_increment = np.diff(
            np.column_stack((np.zeros(len(counts)), recall)), axis=1
        )
        return np.sum(recall_increment * precision, axis=1)

    batch_size = 16
    probabilities = np.full(len(unique_groups), 1 / len(unique_groups))
    positive_by_group = np.bincount(
        group_codes, weights=labels, minlength=len(unique_groups)
    )
    negative_by_group = np.bincount(
        group_codes, weights=1 - labels, minlength=len(unique_groups)
    )
    for start in range(0, BOOTSTRAP_REPLICATES, batch_size):
        rows = min(batch_size, BOOTSTRAP_REPLICATES - start)
        counts = rng.multinomial(len(unique_groups), probabilities, size=rows)
        left_ap = weighted_ap(left_scores, counts)
        right_ap = weighted_ap(right_scores, counts)
        valid = (
            np.isfinite(left_ap)
            & np.isfinite(right_ap)
            & (counts @ positive_by_group > 0)
            & (counts @ negative_by_group > 0)
        )
        deltas.extend((left_ap[valid] - right_ap[valid]).tolist())
    if not deltas:
        raise ValueError("Every group-bootstrap replicate had one label class")
    lower, upper = np.quantile(deltas, (0.025, 0.975))
    return {
        "observed_delta": float(observed),
        "ci_lower": float(lower),
        "ci_upper": float(upper),
        "replicates_requested": BOOTSTRAP_REPLICATES,
        "replicates_valid": len(deltas),
        "resampling_unit": "chromosome_group",
        "seed": seed,
    }


def evaluate(datasets: list[dict[str, Any]]) -> dict[str, Any]:
    development = [row for row in datasets if row["role"] == "development"]
    diagnostics = [
        row for row in datasets if row["role"] == "post_holdout_diagnostic"
    ]
    if len(development) < 2:
        raise ValueError("At least two development species are required")
    report: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "software": {"scikit_learn": sklearn.__version__},
        "policy": {
            "truth_blind_feature_generation": True,
            "development_species_only_for_fitting": True,
            "diagnostic_species_used_for_fitting": False,
            "model": "fold_local_copy_feature_contract_logistic_regression_l2_C_1",
            "stacking_on_in_sample_frozen_scores": False,
            "feature_contract_fit_scope": "training_fold_or_training_species_only",
            "grouped_validation": (
                "five_fold_stratified_by_chromosome_repeated_five_seeds"
            ),
            "uncertainty": "two_thousand_replicate_chromosome_group_bootstrap",
            "automatic_approval": False,
        },
        "feature_names": {
            "baseline": ["original_copy_feature_contract_full"],
            "normalized_wgd_context": [
                "original_copy_feature_contract_full",
                "topology_available",
                "wgd_block_count_within_species_percentile",
                "wgd_block_length_within_species_percentile",
            ],
            "topology": [
                "original_copy_feature_contract_full",
                "topology_available",
                *TOPOLOGY_COMPONENT_FIELDS,
            ],
            "augmented": [
                "original_copy_feature_contract_full",
                "topology_available",
                *TOPOLOGY_COMPONENT_FIELDS,
                "wgd_block_count_within_species_percentile",
                "wgd_block_length_within_species_percentile",
            ],
        },
        "datasets": {},
        "cross_species_transfer": [],
        "pooled_development_to_diagnostic": [],
    }
    for dataset in datasets:
        available = dataset["available"] == 1
        labels = dataset["labels"]
        entry: dict[str, Any] = {
            "role": dataset["role"],
            "candidates": len(labels),
            "positives": int(labels.sum()),
            "topology_available": int(available.sum()),
            "topology_available_positives": int(labels[available].sum()),
            "frozen_model": _metrics(labels, dataset["probability"]),
        }
        if available.any() and len(set(labels[available])) == 2:
            entry["topology_score_within_available"] = _metrics(
                labels[available], dataset["topology_score"][available]
            )
        if dataset["role"] == "development":
            predictions = {
                seed: {
                    matrix_name: _grouped_oof(
                        dataset, matrix_name, seed=seed
                    )
                    for matrix_name in report["feature_names"]
                }
                for seed in OOF_SEEDS
            }
            repeated: dict[str, Any] = {}
            for matrix_name in report["feature_names"]:
                metrics = [
                    _metrics(labels, predictions[seed][matrix_name])
                    for seed in OOF_SEEDS
                ]
                ap_values = [value["average_precision"] for value in metrics]
                repeated[matrix_name] = {
                    "per_seed": [
                        {"seed": seed, **value}
                        for seed, value in zip(OOF_SEEDS, metrics)
                    ],
                    "average_precision_mean": float(np.mean(ap_values)),
                    "average_precision_min": float(np.min(ap_values)),
                    "average_precision_max": float(np.max(ap_values)),
                }
                if matrix_name != "baseline":
                    deltas = [
                        metrics[index]["average_precision"]
                        - _metrics(labels, predictions[seed]["baseline"])[
                            "average_precision"
                        ]
                        for index, seed in enumerate(OOF_SEEDS)
                    ]
                    repeated[matrix_name]["average_precision_delta_vs_baseline"] = {
                        "mean": float(np.mean(deltas)),
                        "min": float(np.min(deltas)),
                        "max": float(np.max(deltas)),
                        "per_seed": [
                            {"seed": seed, "delta": float(delta)}
                            for seed, delta in zip(OOF_SEEDS, deltas)
                        ],
                        "primary_seed_group_bootstrap": _group_bootstrap_delta(
                            labels,
                            predictions[OOF_SEEDS[0]][matrix_name],
                            predictions[OOF_SEEDS[0]]["baseline"],
                            dataset["groups"],
                        ),
                    }
            entry["chromosome_grouped_oof_repeated"] = repeated
        report["datasets"][dataset["name"]] = entry

    for test_set in development:
        train_sets = [row for row in development if row is not test_set]
        transferred = {
            matrix_name: _transfer(train_sets, test_set, matrix_name)
            for matrix_name in report["feature_names"]
        }
        baseline_scores = transferred["baseline"][0]
        feature_results: dict[str, Any] = {}
        for matrix_name, (scores, coefficients, coefficient_names) in transferred.items():
            metrics = _metrics(test_set["labels"], scores)
            feature_results[matrix_name] = {
                **metrics,
                "coefficients": dict(zip(coefficient_names, coefficients)),
            }
            if matrix_name != "baseline":
                feature_results[matrix_name]["delta_vs_baseline"] = (
                    metrics["average_precision"]
                    - _metrics(test_set["labels"], baseline_scores)[
                        "average_precision"
                    ]
                )
                feature_results[matrix_name]["group_bootstrap"] = (
                    _group_bootstrap_delta(
                        test_set["labels"],
                        scores,
                        baseline_scores,
                        test_set["groups"],
                    )
                )
        report["cross_species_transfer"].append(
            {
                "train": [row["name"] for row in train_sets],
                "test": test_set["name"],
                "feature_sets": feature_results,
            }
        )

    for diagnostic in diagnostics:
        transferred = {
            matrix_name: _transfer(development, diagnostic, matrix_name)
            for matrix_name in report["feature_names"]
        }
        baseline_scores = transferred["baseline"][0]
        feature_results = {}
        for matrix_name, (scores, coefficients, coefficient_names) in transferred.items():
            metrics = _metrics(diagnostic["labels"], scores)
            feature_results[matrix_name] = {
                **metrics,
                "coefficients": dict(zip(coefficient_names, coefficients)),
            }
            if matrix_name != "baseline":
                feature_results[matrix_name]["delta_vs_baseline"] = (
                    metrics["average_precision"]
                    - _metrics(diagnostic["labels"], baseline_scores)[
                        "average_precision"
                    ]
                )
                feature_results[matrix_name]["group_bootstrap"] = (
                    _group_bootstrap_delta(
                        diagnostic["labels"],
                        scores,
                        baseline_scores,
                        diagnostic["groups"],
                    )
                )
        report["pooled_development_to_diagnostic"].append(
            {
                "train": [row["name"] for row in development],
                "test": diagnostic["name"],
                "claim_status": "post_holdout_diagnostic_not_new_formal_holdout",
                "feature_sets": feature_results,
            }
        )
    return report


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--dataset",
        action="append",
        required=True,
        help="NAME=ROLE,TOPOLOGY,LABELS,SCORES",
    )
    parser.add_argument("--output-json", required=True)
    args = parser.parse_args()
    parsed = [_parse_dataset(value) for value in args.dataset]
    names = [value[0] for value in parsed]
    if len(names) != len(set(names)):
        raise ValueError("Dataset names must be unique")
    datasets = [_load_dataset(*value) for value in parsed]
    report = evaluate(datasets)
    report["inputs"] = {
        dataset["name"]: {
            "topology": _sha256(dataset["paths"][0]),
            "labels": _sha256(dataset["paths"][1]),
            "scores": _sha256(dataset["paths"][2]),
        }
        for dataset in datasets
    }
    output = Path(args.output_json)
    if output.exists():
        raise FileExistsError(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="",
    )
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
