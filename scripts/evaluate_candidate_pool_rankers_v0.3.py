#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import platform
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import scipy
import sklearn
from scipy.optimize import minimize
from scipy.special import expit
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import GroupKFold

from ploidypatch import __version__ as ploidypatch_version
from ploidypatch.copy_features import FEATURE_COLUMNS
from ploidypatch.copy_model import fit_copy_feature_contract, transform_copy_feature_row
from ploidypatch.homeolog_ranker import TOPOLOGY_COMPONENT_FIELDS


SCHEMA_VERSION = "ploidypatch.candidate_pool_ranker_evaluation.v4"
OOF_SEEDS = (20260808, 20260818, 20260828, 20260838, 20260848)
SUPPORT_PATTERNS = (
    "gemoma",
    "lifton",
    "miniprot",
    "gemoma,lifton",
    "gemoma,miniprot",
    "lifton,miniprot",
    "gemoma,lifton,miniprot",
)
MODEL_NAMES = (
    "baseline",
    "global_refit",
    "offset_global",
    "offset_support_conditioned",
)
FIXED_C = 1.0


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def read_index(
    path: Path, key_fields: tuple[str, ...] = ("candidate_digest",)
) -> tuple[list[str], dict[str, dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Missing header: {path}")
        fields = list(reader.fieldnames)
        key_field = next((field for field in key_fields if field in fields), None)
        if key_field is None:
            raise ValueError(f"Missing candidate key {key_fields}: {path}")
        rows = list(reader)
    indexed = {row[key_field]: row for row in rows}
    if len(indexed) != len(rows) or "" in indexed:
        raise ValueError(f"Empty or duplicate digest: {path}")
    return fields, indexed


@dataclass
class Dataset:
    name: str
    role: str
    paths: tuple[Path, Path, Path]
    digests: list[str]
    rows: list[dict[str, str]]
    labels: np.ndarray
    groups: np.ndarray
    global_correction: np.ndarray
    support_correction: np.ndarray
    conflict_sets: np.ndarray
    topology_available: np.ndarray


def parse_dataset(value: str) -> Dataset:
    name, separator, payload = value.partition("=")
    parts = payload.split(",")
    if not separator or not name or len(parts) != 4:
        raise ValueError("--dataset requires NAME=ROLE,LABELS,TOPOLOGY,DECISIONS")
    role, label_value, topology_value, decision_value = parts
    if role not in {"development", "retrospective_diagnostic"}:
        raise ValueError(f"Invalid role for {name}: {role}")
    paths = (Path(label_value), Path(topology_value), Path(decision_value))
    if any(not path.is_file() for path in paths):
        raise FileNotFoundError(f"Missing dataset input for {name}")
    label_fields, labels = read_index(paths[0])
    _, topology = read_index(paths[1])
    _, decisions = read_index(paths[2], ("candidate_digest", "consensus_digest"))
    required_labels = {*FEATURE_COLUMNS, "label_exact_cds"}
    if not required_labels <= set(label_fields):
        raise ValueError(f"Labeled features lack contract columns for {name}")
    if set(labels) != set(topology) or set(labels) != set(decisions):
        raise ValueError(f"Candidate universes differ for {name}")
    digests = sorted(labels)
    rows = [labels[digest] for digest in digests]
    y = np.asarray([int(labels[digest]["label_exact_cds"]) for digest in digests], dtype=np.uint8)
    if not 0 < int(y.sum()) < len(y):
        raise ValueError(f"Dataset needs both label classes: {name}")
    available = np.asarray([int(topology[digest]["topology_available"]) for digest in digests], dtype=np.uint8)
    components = np.asarray(
        [
            [float(topology[digest][field]) if available[index] else 0.0 for field in TOPOLOGY_COMPONENT_FIELDS]
            for index, digest in enumerate(digests)
        ],
        dtype=float,
    )
    global_correction = np.column_stack((available, components))
    mean_coherence = components.mean(axis=1)
    support_columns: list[np.ndarray] = [components]
    observed_patterns = {row["support_methods"] for row in rows}
    if observed_patterns - set(SUPPORT_PATTERNS):
        raise ValueError(f"Unexpected support pattern for {name}: {sorted(observed_patterns - set(SUPPORT_PATTERNS))}")
    for pattern in SUPPORT_PATTERNS:
        indicator = np.asarray([row["support_methods"] == pattern for row in rows], dtype=float)
        support_columns.extend((available * indicator, available * mean_coherence * indicator))
    support_correction = np.column_stack(support_columns)
    return Dataset(
        name=name,
        role=role,
        paths=paths,
        digests=digests,
        rows=rows,
        labels=y,
        groups=np.asarray([row["seqid"] for row in rows]),
        global_correction=global_correction,
        support_correction=support_correction,
        conflict_sets=np.asarray([decisions[digest]["conflict_set_digest"] for digest in digests]),
        topology_available=available,
    )


def fit_offset(offset: np.ndarray, correction: np.ndarray, labels: np.ndarray) -> np.ndarray:
    def objective(beta: np.ndarray) -> tuple[float, np.ndarray]:
        logits = offset + correction @ beta
        loss = float(np.sum(np.logaddexp(0.0, logits) - labels * logits) + 0.5 / FIXED_C * np.dot(beta, beta))
        gradient = correction.T @ (expit(logits) - labels) + beta / FIXED_C
        return loss, gradient

    result = minimize(
        objective,
        np.zeros(correction.shape[1], dtype=float),
        method="L-BFGS-B",
        jac=True,
        options={"maxiter": 5000, "ftol": 1e-12, "gtol": 1e-8},
    )
    if not result.success or not np.isfinite(result.x).all():
        raise RuntimeError(f"Offset optimization failed: {result.message}")
    return np.asarray(result.x, dtype=float)


def fit_predict(train: list[Dataset], test: Dataset) -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    train_rows = [row for dataset in train for row in dataset.rows]
    train_y = np.concatenate([dataset.labels for dataset in train])
    contract = fit_copy_feature_contract(train_rows, feature_set="full")
    train_base = np.asarray([transform_copy_feature_row(row, contract) for row in train_rows], dtype=float)
    test_base = np.asarray([transform_copy_feature_row(row, contract) for row in test.rows], dtype=float)
    train_global = np.concatenate([dataset.global_correction for dataset in train])
    train_support = np.concatenate([dataset.support_correction for dataset in train])
    baseline = LogisticRegression(C=FIXED_C, solver="lbfgs", max_iter=5000, random_state=0).fit(train_base, train_y)
    global_refit = LogisticRegression(C=FIXED_C, solver="lbfgs", max_iter=5000, random_state=0).fit(
        np.column_stack((train_base, train_global)), train_y
    )
    train_offset = baseline.decision_function(train_base)
    beta_global = fit_offset(train_offset, train_global, train_y)
    beta_support = fit_offset(train_offset, train_support, train_y)
    test_offset = baseline.decision_function(test_base)
    scores = {
        "baseline": test_offset,
        "global_refit": global_refit.decision_function(np.column_stack((test_base, test.global_correction))),
        "offset_global": test_offset + test.global_correction @ beta_global,
        "offset_support_conditioned": test_offset + test.support_correction @ beta_support,
    }
    artifacts = {
        "base_feature_names": contract["expanded_feature_names"],
        "baseline_intercept": float(baseline.intercept_[0]),
        "baseline_coefficients": baseline.coef_[0].tolist(),
        "global_refit_intercept": float(global_refit.intercept_[0]),
        "global_refit_coefficients": global_refit.coef_[0].tolist(),
        "offset_global_coefficients": beta_global.tolist(),
        "offset_support_coefficients": beta_support.tolist(),
    }
    return scores, artifacts


def subset(dataset: Dataset, indexes: np.ndarray, suffix: str) -> Dataset:
    return Dataset(
        name=f"{dataset.name}_{suffix}",
        role=dataset.role,
        paths=dataset.paths,
        digests=[dataset.digests[index] for index in indexes],
        rows=[dataset.rows[index] for index in indexes],
        labels=dataset.labels[indexes],
        groups=dataset.groups[indexes],
        global_correction=dataset.global_correction[indexes],
        support_correction=dataset.support_correction[indexes],
        conflict_sets=dataset.conflict_sets[indexes],
        topology_available=dataset.topology_available[indexes],
    )


def metrics(labels: np.ndarray, scores: np.ndarray) -> dict[str, float]:
    return {
        "average_precision": float(average_precision_score(labels, scores)),
        "roc_auc": float(roc_auc_score(labels, scores)),
    }


def percentile_scores(values: np.ndarray) -> np.ndarray:
    if len(values) == 1:
        return np.asarray([0.5])
    order = np.argsort(values, kind="mergesort")
    ranks = np.empty(len(values), dtype=float)
    cursor = 0
    while cursor < len(order):
        end = cursor + 1
        while end < len(order) and values[order[end]] == values[order[cursor]]:
            end += 1
        ranks[order[cursor:end]] = (cursor + end - 1) / 2
        cursor = end
    return ranks / (len(values) - 1)


def conflict_metrics(dataset: Dataset, scores: np.ndarray) -> dict[str, Any]:
    groups: dict[str, list[int]] = {}
    for index, conflict in enumerate(dataset.conflict_sets):
        if conflict:
            groups.setdefault(str(conflict), []).append(index)
    evaluable = [indexes for indexes in groups.values() if len(indexes) >= 2 and int(dataset.labels[indexes].sum()) == 1]
    reciprocal_ranks: list[float] = []
    top1 = 0
    for indexes in evaluable:
        ordered = sorted(indexes, key=lambda index: (-scores[index], dataset.digests[index]))
        positive_rank = 1 + next(position for position, index in enumerate(ordered) if dataset.labels[index] == 1)
        reciprocal_ranks.append(1 / positive_rank)
        top1 += int(positive_rank == 1)
    return {
        "conflict_sets_total": len(groups),
        "evaluable_exactly_one_positive": len(evaluable),
        "top1_correct": top1,
        "top1_accuracy": top1 / len(evaluable) if evaluable else None,
        "mean_reciprocal_rank": float(np.mean(reciprocal_ranks)) if reciprocal_ranks else None,
    }


def review_metrics(dataset: Dataset, scores: np.ndarray) -> dict[str, Any]:
    output: dict[str, Any] = {}
    for fraction in (0.005, 0.01, 0.02):
        count = max(1, int(np.ceil(len(scores) * fraction)))
        selected = sorted(range(len(scores)), key=lambda index: (-scores[index], dataset.digests[index]))[:count]
        true_positive = int(dataset.labels[selected].sum())
        output[f"top_{fraction * 100:g}pct"] = {
            "reviewed": count,
            "true_positive": true_positive,
            "precision": true_positive / count,
            "recall_of_800_events": true_positive / 800,
        }
    return output


def grouped_oof(dataset: Dataset) -> tuple[dict[str, np.ndarray], list[dict[str, Any]]]:
    accumulators = {name: np.zeros(len(dataset.labels), dtype=float) for name in MODEL_NAMES}
    seed_reports: list[dict[str, Any]] = []
    placeholder = np.zeros((len(dataset.labels), 1), dtype=float)
    seen_partitions: set[tuple[tuple[str, ...], ...]] = set()
    for seed in OOF_SEEDS:
        predictions = {name: np.full(len(dataset.labels), np.nan) for name in MODEL_NAMES}
        fold_models: list[dict[str, Any]] = []
        splitter = GroupKFold(n_splits=5, shuffle=True, random_state=seed)
        splits = list(splitter.split(placeholder, dataset.labels, dataset.groups))
        partition = tuple(
            sorted(tuple(sorted(set(dataset.groups[test_indexes].tolist()))) for _, test_indexes in splits)
        )
        if partition in seen_partitions:
            raise ValueError(f"Repeated chromosome partition for {dataset.name}, seed {seed}")
        seen_partitions.add(partition)
        for fold, (train_indexes, test_indexes) in enumerate(splits, start=1):
            train_groups = set(dataset.groups[train_indexes].tolist())
            test_groups = set(dataset.groups[test_indexes].tolist())
            if train_groups & test_groups:
                raise AssertionError("Chromosome leakage between OOF train and test")
            if len(np.unique(dataset.labels[train_indexes])) != 2 or len(np.unique(dataset.labels[test_indexes])) != 2:
                raise ValueError(f"OOF fold lacks a label class for {dataset.name}, seed {seed}, fold {fold}")
            fold_scores, artifacts = fit_predict(
                [subset(dataset, train_indexes, f"train_{seed}_{fold}")],
                subset(dataset, test_indexes, f"test_{seed}_{fold}"),
            )
            for name in MODEL_NAMES:
                predictions[name][test_indexes] = fold_scores[name]
            fold_models.append(
                {
                    "fold": fold,
                    "train_groups": sorted(train_groups),
                    "test_groups": sorted(test_groups),
                    "artifacts": artifacts,
                }
            )
        if any(np.isnan(values).any() for values in predictions.values()):
            raise AssertionError("OOF predictions are incomplete")
        seed_reports.append(
            {
                "seed": seed,
                "metrics": {name: metrics(dataset.labels, predictions[name]) for name in MODEL_NAMES},
                "fold_models": fold_models,
            }
        )
        for name in MODEL_NAMES:
            accumulators[name] += predictions[name]
    return {name: values / len(OOF_SEEDS) for name, values in accumulators.items()}, seed_reports


def weighted_ap_contract(labels: np.ndarray, scores: np.ndarray, group_codes: np.ndarray, group_count: int) -> tuple[np.ndarray, np.ndarray]:
    order = np.argsort(scores, kind="mergesort")[::-1]
    ordered_scores = scores[order]
    ordered_labels = labels[order]
    ordered_groups = group_codes[order]
    thresholds = np.concatenate((np.flatnonzero(np.diff(ordered_scores)), [len(scores) - 1]))
    total = np.zeros((len(scores), group_count), dtype=np.int8)
    positive = np.zeros((len(scores), group_count), dtype=np.int8)
    row_indexes = np.arange(len(scores))
    total[row_indexes, ordered_groups] = 1
    positive[row_indexes, ordered_groups] = ordered_labels
    return np.cumsum(positive, axis=0)[thresholds], np.cumsum(total, axis=0)[thresholds]


def weighted_ap_batch(counts: np.ndarray, cumulative_positive: np.ndarray, cumulative_total: np.ndarray) -> np.ndarray:
    true_at = counts @ cumulative_positive.T
    total_at = counts @ cumulative_total.T
    precision = np.divide(true_at, total_at, out=np.zeros_like(true_at, dtype=float), where=total_at != 0)
    total_positive = true_at[:, -1]
    recall = np.divide(true_at, total_positive[:, None], out=np.zeros_like(true_at, dtype=float), where=total_positive[:, None] != 0)
    increments = np.diff(np.column_stack((np.zeros(len(counts)), recall)), axis=1)
    return np.sum(increments * precision, axis=1)


def group_bootstrap_delta(labels: np.ndarray, primary: np.ndarray, baseline: np.ndarray, groups: np.ndarray, *, replicates: int, seed: int) -> dict[str, Any]:
    unique_groups, group_codes = np.unique(groups, return_inverse=True)
    primary_contract = weighted_ap_contract(labels, primary, group_codes, len(unique_groups))
    baseline_contract = weighted_ap_contract(labels, baseline, group_codes, len(unique_groups))
    positive_by_group = np.bincount(group_codes, weights=labels, minlength=len(unique_groups))
    negative_by_group = np.bincount(group_codes, weights=1 - labels, minlength=len(unique_groups))
    rng = np.random.default_rng(seed)
    probabilities = np.full(len(unique_groups), 1 / len(unique_groups))
    deltas: list[float] = []
    for start in range(0, replicates, 64):
        rows = min(64, replicates - start)
        counts = rng.multinomial(len(unique_groups), probabilities, size=rows)
        primary_ap = weighted_ap_batch(counts, *primary_contract)
        baseline_ap = weighted_ap_batch(counts, *baseline_contract)
        valid = (counts @ positive_by_group > 0) & (counts @ negative_by_group > 0)
        deltas.extend((primary_ap[valid] - baseline_ap[valid]).tolist())
    lower, upper = np.quantile(deltas, (0.025, 0.975))
    return {
        "observed_delta": float(average_precision_score(labels, primary) - average_precision_score(labels, baseline)),
        "ci_lower": float(lower),
        "ci_upper": float(upper),
        "replicates_requested": replicates,
        "replicates_valid": len(deltas),
        "groups": unique_groups.tolist(),
        "seed": seed,
    }


def dataset_report(dataset: Dataset, scores: dict[str, np.ndarray]) -> dict[str, Any]:
    return {
        "counts": {
            "candidates": len(dataset.labels),
            "positives": int(dataset.labels.sum()),
            "topology_available": int(dataset.topology_available.sum()),
            "topology_positive": int(dataset.labels[dataset.topology_available == 1].sum()),
        },
        "metrics": {name: metrics(dataset.labels, scores[name]) for name in MODEL_NAMES},
        "conflicts": {name: conflict_metrics(dataset, scores[name]) for name in MODEL_NAMES},
        "review_budgets": {name: review_metrics(dataset, scores[name]) for name in MODEL_NAMES},
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Evaluate fixed conflict-aware rankers on v0.3 candidate pools")
    parser.add_argument("--dataset", action="append", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--protocol", required=True)
    parser.add_argument("--bootstrap-replicates", type=int, default=20_000)
    parser.add_argument("--bootstrap-seed", type=int, default=20260808)
    args = parser.parse_args()
    datasets = [parse_dataset(value) for value in args.dataset]
    protocol_path = Path(args.protocol)
    if not protocol_path.is_file() or protocol_path.stat().st_size == 0:
        raise FileNotFoundError("Missing v0.3 ranker protocol")
    names = [dataset.name for dataset in datasets]
    if len(names) != len(set(names)):
        raise ValueError("Dataset names must be unique")
    development = [dataset for dataset in datasets if dataset.role == "development"]
    diagnostics = [dataset for dataset in datasets if dataset.role == "retrospective_diagnostic"]
    if len(development) != 2 or not diagnostics:
        raise ValueError("Exactly two development species and at least one diagnostic are required")
    output = Path(args.output_dir)
    partial = Path(str(output) + ".partial")
    if output.exists() or partial.exists():
        raise FileExistsError("Refusing to overwrite v0.3 ranker evaluation")
    partial.mkdir(parents=True)

    oof_scores: dict[str, dict[str, np.ndarray]] = {}
    oof_reports: dict[str, Any] = {}
    for dataset in development:
        scores, seeds = grouped_oof(dataset)
        oof_scores[dataset.name] = scores
        oof_reports[dataset.name] = {
            **dataset_report(dataset, scores),
            "seed_metrics": seeds,
            "primary_ap_delta_bootstrap": group_bootstrap_delta(
                dataset.labels,
                scores["offset_support_conditioned"],
                scores["baseline"],
                dataset.groups,
                replicates=args.bootstrap_replicates,
                seed=args.bootstrap_seed,
            ),
        }

    pooled_labels = np.concatenate([dataset.labels for dataset in development])
    pooled_groups = np.concatenate(
        [np.asarray([f"{dataset.name}:{group}" for group in dataset.groups]) for dataset in development]
    )
    pooled_scores = {
        name: np.concatenate(
            [percentile_scores(oof_scores[dataset.name][name]) for dataset in development]
        )
        for name in MODEL_NAMES
    }
    pooled_report = {
        "score_harmonization": "within_species_percentile_before_pooling",
        "metrics": {name: metrics(pooled_labels, pooled_scores[name]) for name in MODEL_NAMES},
        "primary_ap_delta_bootstrap": group_bootstrap_delta(
            pooled_labels,
            pooled_scores["offset_support_conditioned"],
            pooled_scores["baseline"],
            pooled_groups,
            replicates=args.bootstrap_replicates,
            seed=args.bootstrap_seed,
        ),
    }

    transfer_reports: dict[str, Any] = {}
    transfer_predictions: dict[str, dict[str, np.ndarray]] = {}
    for train in development:
        test = next(dataset for dataset in development if dataset.name != train.name)
        scores, _ = fit_predict([train], test)
        transfer_predictions[f"{train.name}_to_{test.name}"] = scores
        transfer_reports[f"{train.name}_to_{test.name}"] = dataset_report(test, scores)

    diagnostic_reports: dict[str, Any] = {}
    diagnostic_predictions: dict[str, dict[str, np.ndarray]] = {}
    pooled_fit_artifacts: dict[str, Any] | None = None
    for diagnostic in diagnostics:
        scores, artifacts = fit_predict(development, diagnostic)
        diagnostic_predictions[diagnostic.name] = scores
        diagnostic_reports[diagnostic.name] = dataset_report(diagnostic, scores)
        pooled_fit_artifacts = artifacts

    criteria = {
        dataset.name: {
            "ap_delta": (
                oof_reports[dataset.name]["metrics"]["offset_support_conditioned"]["average_precision"]
                - oof_reports[dataset.name]["metrics"]["baseline"]["average_precision"]
            ),
            "conflict_top1_delta": (
                oof_reports[dataset.name]["conflicts"]["offset_support_conditioned"]["top1_accuracy"]
                - oof_reports[dataset.name]["conflicts"]["baseline"]["top1_accuracy"]
            ),
        }
        for dataset in development
    }
    retain = (
        all(value["ap_delta"] > 0 and value["conflict_top1_delta"] >= 0 for value in criteria.values())
        and pooled_report["metrics"]["offset_support_conditioned"]["average_precision"]
        > pooled_report["metrics"]["baseline"]["average_precision"]
    )
    report: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "code_commit": os.environ.get("PLOIDYPATCH_CODE_COMMIT", "unavailable_server_mirror"),
        "evaluation_contract": {
            "primary_estimator": "offset_support_conditioned",
            "fixed_C": FIXED_C,
            "oof_seeds": list(OOF_SEEDS),
            "folds": 5,
            "splitter": "GroupKFold(shuffle=True) with distinct-partition and label-class audits",
            "bootstrap_replicates": args.bootstrap_replicates,
            "bootstrap_seed": args.bootstrap_seed,
            "diagnostics_cannot_select_model": True,
            "automatic_approval": False,
        },
        "development_oof": oof_reports,
        "pooled_development_oof": pooled_report,
        "cross_species_transfer": transfer_reports,
        "retrospective_diagnostics": diagnostic_reports,
        "primary_retention_criteria": {"by_species": criteria, "retain_for_external_freeze": retain},
        "pooled_fit_artifacts": pooled_fit_artifacts,
        "inputs": [
            {
                "name": dataset.name,
                "role": dataset.role,
                "rows": len(dataset.labels),
                "positives": int(dataset.labels.sum()),
                "files": [
                    {"path": str(path), "bytes": path.stat().st_size, "sha256": sha256(path)}
                    for path in dataset.paths
                ],
            }
            for dataset in datasets
        ],
        "protocol_artifact": {
            "path": str(protocol_path),
            "bytes": protocol_path.stat().st_size,
            "sha256": sha256(protocol_path),
        },
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "scikit_learn": sklearn.__version__,
            "ploidypatch": ploidypatch_version,
        },
    }
    with (partial / "evaluation.json").open("x", encoding="utf-8", newline="") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
        handle.write("\n")
    with (partial / "predictions.tsv").open("x", encoding="utf-8", newline="") as handle:
        fields = ("dataset", "evaluation", "candidate_digest", "seqid", "label", "conflict_set", *MODEL_NAMES)
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for dataset in development:
            for index, digest in enumerate(dataset.digests):
                writer.writerow(
                    {
                        "dataset": dataset.name,
                        "evaluation": "mean_repeated_grouped_oof",
                        "candidate_digest": digest,
                        "seqid": dataset.groups[index],
                        "label": int(dataset.labels[index]),
                        "conflict_set": dataset.conflict_sets[index],
                        **{name: format(oof_scores[dataset.name][name][index], ".17g") for name in MODEL_NAMES},
                    }
                )
        for dataset in diagnostics:
            for index, digest in enumerate(dataset.digests):
                writer.writerow(
                    {
                        "dataset": dataset.name,
                        "evaluation": "pooled_development_to_retrospective_diagnostic",
                        "candidate_digest": digest,
                        "seqid": dataset.groups[index],
                        "label": int(dataset.labels[index]),
                        "conflict_set": dataset.conflict_sets[index],
                        **{name: format(diagnostic_predictions[dataset.name][name][index], ".17g") for name in MODEL_NAMES},
                    }
                )
    with (partial / "run_contract.tsv").open("x", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("field", "value"))
        writer.writerow(("code_commit", os.environ.get("PLOIDYPATCH_CODE_COMMIT", "unavailable_server_mirror")))
        writer.writerow(("primary_estimator", "offset_support_conditioned"))
        writer.writerow(("protocol_sha256", sha256(protocol_path)))
        writer.writerow(("development_species", ",".join(dataset.name for dataset in development)))
        writer.writerow(("diagnostic_species", ",".join(dataset.name for dataset in diagnostics)))
        writer.writerow(("automatic_approval", "false"))
    artifacts = sorted(path for path in partial.iterdir())
    with (partial / "SHA256SUMS").open("x", encoding="utf-8", newline="") as handle:
        for path in artifacts:
            handle.write(f"{sha256(path)}  {path.name}\n")
    os.replace(partial, output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
