#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import platform
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
import scipy
import sklearn
from scipy.optimize import minimize
from scipy.special import expit
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import GroupKFold

from ploidypatch import __version__ as ploidypatch_version
from ploidypatch.artifact_manifest import write_sha256sums, verify_sha256sums
from ploidypatch.baseline import _file_sha256
from ploidypatch.conflict_guard import compute_conflict_winner_guard
from ploidypatch.copy_features import FEATURE_COLUMNS
from ploidypatch.copy_model import transform_copy_feature_row
from ploidypatch.stable_ranker import (
    STABLE_REFERENCE_RANKER_SCHEMA_VERSION,
    fit_weighted_copy_feature_contract,
)
from ploidypatch.support_ranker import (
    CORRECTION_FEATURE_NAMES,
    support_conditioned_correction_vector,
)


EVALUATION_SCHEMA_VERSION = "ploidypatch.stable_reference_ranker_evaluation.v0.9"
OOF_SEEDS = (20260901, 20260911, 20260921, 20260931, 20260941)
BOOTSTRAP_SEEDS = {"actinidia": 20260951, "populus": 20260952}
BOOTSTRAP_REPLICATES = 20_000
FIXED_C = 1.0
FOLDS = 5
FRACTION_BUDGETS = (0.005, 0.01, 0.02)
FIXED_BUDGETS = (100, 250, 500)


def sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def _read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Missing TSV header: {path}")
        return list(reader.fieldnames), list(reader)


def _index(
    path: Path, key_fields: Sequence[str]
) -> tuple[list[str], dict[str, dict[str, str]]]:
    fields, rows = _read_tsv(path)
    key = next((field for field in key_fields if field in fields), None)
    if key is None:
        raise ValueError(f"Missing candidate key in {path}: {key_fields}")
    output = {row[key]: row for row in rows}
    if len(output) != len(rows) or "" in output:
        raise ValueError(f"Empty or duplicate candidate key: {path}")
    return fields, output


def _verify_truth_free_manifest(path: Path, artifact: Path) -> dict[str, Any]:
    if not path.is_file() or not artifact.is_file():
        raise FileNotFoundError(f"Missing truth-free artifact or manifest: {artifact}")
    with path.open("r", encoding="utf-8") as handle:
        manifest = json.load(handle)
    if not isinstance(manifest, dict) or manifest.get("truth_access") is not False:
        raise ValueError(f"Truth-free manifest assertion failed: {path}")
    expected = _file_sha256(artifact)
    observed: list[str] = []

    def visit(value: Any) -> None:
        if isinstance(value, Mapping):
            for key, item in value.items():
                if key == "sha256" and isinstance(item, str):
                    observed.append(item)
                else:
                    visit(item)
        elif isinstance(value, list):
            for item in value:
                visit(item)

    visit(manifest.get("outputs", {}))
    if observed.count(expected) != 1:
        raise ValueError(
            f"Manifest must bind the artifact SHA exactly once: {path} -> {artifact}"
        )
    return manifest


@dataclass(frozen=True)
class Dataset:
    name: str
    paths: tuple[Path, Path, Path, Path]
    manifest_paths: tuple[Path, Path]
    digests: tuple[str, ...]
    copy_rows: tuple[dict[str, str], ...]
    topology_rows: tuple[dict[str, str], ...]
    labels: np.ndarray
    groups: np.ndarray
    conflict_sets: np.ndarray
    correction: np.ndarray
    topology_available: np.ndarray


def parse_dataset(value: str) -> Dataset:
    name, separator, payload = value.partition("=")
    parts = payload.split(",")
    if not separator or name not in {"actinidia", "populus"} or len(parts) != 4:
        raise ValueError(
            "--dataset requires actinidia|populus=COPY,TOPOLOGY,LABELS,DECISIONS"
        )
    paths = tuple(Path(item) for item in parts)
    if any(not path.is_file() or path.stat().st_size == 0 for path in paths):
        raise FileNotFoundError(f"Missing or empty stable-ranker input for {name}")
    copy_path, topology_path, label_path, decisions_path = paths
    manifest_paths = (
        Path(str(copy_path) + ".manifest.json"),
        Path(str(topology_path) + ".manifest.json"),
    )
    _verify_truth_free_manifest(manifest_paths[0], copy_path)
    topology_manifest = _verify_truth_free_manifest(manifest_paths[1], topology_path)
    topology_inputs = json.dumps(topology_manifest.get("inputs", {}), sort_keys=True)
    if _file_sha256(copy_path) not in topology_inputs:
        raise ValueError("Topology manifest does not bind the stable copy features")

    copy_fields, copy_by_digest = _index(copy_path, ("candidate_digest",))
    topology_fields, topology_by_digest = _index(
        topology_path, ("candidate_digest",)
    )
    label_fields, label_by_digest = _index(
        label_path, ("candidate_digest", "consensus_digest")
    )
    decision_fields, decision_by_digest = _index(
        decisions_path, ("candidate_digest", "consensus_digest")
    )
    if not set(FEATURE_COLUMNS) <= set(copy_fields):
        raise ValueError(f"Stable copy features violate v1 feature contract: {name}")
    if not {
        "candidate_digest",
        "topology_available",
        "cds_bp_ratio",
        "cds_segment_count_ratio",
        "phase_lcs_similarity",
        "junction_fraction_similarity",
        "coding_span_ratio",
    } <= set(topology_fields):
        raise ValueError(f"Stable topology fields are incomplete: {name}")
    if "label_exact_cds" not in label_fields:
        raise ValueError(f"Missing exact-CDS label: {name}")
    if "status" in decision_fields:
        decision_by_digest = {
            digest: row
            for digest, row in decision_by_digest.items()
            if row["status"] == "accepted"
        }
    if "conflict_set_digest" not in decision_fields:
        raise ValueError(f"Missing conflict-set provenance: {name}")
    universes = tuple(
        set(mapping)
        for mapping in (
            copy_by_digest,
            topology_by_digest,
            label_by_digest,
            decision_by_digest,
        )
    )
    if any(universe != universes[0] for universe in universes[1:]):
        raise ValueError(f"Candidate universes differ: {name}")
    digests = tuple(sorted(universes[0]))
    copies = tuple(copy_by_digest[digest] for digest in digests)
    topologies = tuple(topology_by_digest[digest] for digest in digests)
    labels = np.asarray(
        [int(label_by_digest[digest]["label_exact_cds"]) for digest in digests],
        dtype=np.uint8,
    )
    if any(value not in {0, 1} for value in labels) or not 0 < int(labels.sum()) < len(labels):
        raise ValueError(f"Stable development dataset needs both classes: {name}")
    groups = np.asarray([row["seqid"] for row in copies])
    if "" in groups or len(np.unique(groups)) < FOLDS:
        raise ValueError(f"Insufficient target chromosome groups: {name}")
    correction = np.asarray(
        [
            support_conditioned_correction_vector(copy_row, topology_row)
            for copy_row, topology_row in zip(copies, topologies, strict=True)
        ],
        dtype=float,
    )
    available = np.asarray(
        [int(row["topology_available"]) for row in topologies], dtype=np.uint8
    )
    if np.any(correction[available == 0] != 0):
        raise AssertionError("Unavailable topology produced a nonzero correction")
    return Dataset(
        name=name,
        paths=paths,
        manifest_paths=manifest_paths,
        digests=digests,
        copy_rows=copies,
        topology_rows=topologies,
        labels=labels,
        groups=groups,
        conflict_sets=np.asarray(
            [decision_by_digest[digest]["conflict_set_digest"] for digest in digests]
        ),
        correction=correction,
        topology_available=available,
    )


def subset(dataset: Dataset, indexes: np.ndarray, suffix: str) -> Dataset:
    return Dataset(
        name=f"{dataset.name}:{suffix}",
        paths=dataset.paths,
        manifest_paths=dataset.manifest_paths,
        digests=tuple(dataset.digests[index] for index in indexes),
        copy_rows=tuple(dataset.copy_rows[index] for index in indexes),
        topology_rows=tuple(dataset.topology_rows[index] for index in indexes),
        labels=dataset.labels[indexes],
        groups=dataset.groups[indexes],
        conflict_sets=dataset.conflict_sets[indexes],
        correction=dataset.correction[indexes],
        topology_available=dataset.topology_available[indexes],
    )


def species_equal_weights(datasets: Sequence[Dataset]) -> np.ndarray:
    if not datasets or any(len(dataset.labels) == 0 for dataset in datasets):
        raise ValueError("Species balancing requires nonempty datasets")
    total = sum(len(dataset.labels) for dataset in datasets)
    count = len(datasets)
    weights = np.concatenate(
        [
            np.full(len(dataset.labels), total / (count * len(dataset.labels)))
            for dataset in datasets
        ]
    )
    if not math.isclose(float(weights.sum()), total, rel_tol=1e-12):
        raise AssertionError("Species-balanced weights must have mean one")
    totals = []
    cursor = 0
    for dataset in datasets:
        totals.append(float(weights[cursor : cursor + len(dataset.labels)].sum()))
        cursor += len(dataset.labels)
    if max(totals) - min(totals) > 1e-9 * total:
        raise AssertionError("Species-balanced total weights differ")
    return weights


def fit_offset(
    offset: np.ndarray,
    correction: np.ndarray,
    labels: np.ndarray,
    sample_weight: np.ndarray,
) -> np.ndarray:
    if not (
        len(offset) == len(correction) == len(labels) == len(sample_weight)
    ):
        raise ValueError("Offset fit inputs differ in length")

    def objective(beta: np.ndarray) -> tuple[float, np.ndarray]:
        logits = offset + correction @ beta
        residual = expit(logits) - labels
        loss = float(
            np.sum(sample_weight * (np.logaddexp(0.0, logits) - labels * logits))
            + 0.5 / FIXED_C * np.dot(beta, beta)
        )
        gradient = correction.T @ (sample_weight * residual) + beta / FIXED_C
        return loss, gradient

    result = minimize(
        objective,
        np.zeros(correction.shape[1], dtype=float),
        method="L-BFGS-B",
        jac=True,
        options={"maxiter": 5000, "ftol": 1e-12, "gtol": 1e-8},
    )
    if not result.success or not np.isfinite(result.x).all():
        raise RuntimeError(f"Stable topology offset optimization failed: {result.message}")
    return np.asarray(result.x, dtype=float)


def fit_predict(
    train: Sequence[Dataset], test: Dataset
) -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    train_rows = [row for dataset in train for row in dataset.copy_rows]
    labels = np.concatenate([dataset.labels for dataset in train])
    correction = np.concatenate([dataset.correction for dataset in train])
    weights = species_equal_weights(train)
    contract = fit_weighted_copy_feature_contract(
        train_rows, sample_weight=weights.tolist(), feature_set="full"
    )
    train_base = np.asarray(
        [transform_copy_feature_row(row, contract) for row in train_rows], dtype=float
    )
    test_base = np.asarray(
        [transform_copy_feature_row(row, contract) for row in test.copy_rows],
        dtype=float,
    )
    baseline = LogisticRegression(
        C=FIXED_C,
        solver="lbfgs",
        max_iter=5000,
        random_state=0,
    ).fit(train_base, labels, sample_weight=weights)
    train_offset = baseline.decision_function(train_base)
    beta = fit_offset(train_offset, correction, labels, weights)
    test_offset = baseline.decision_function(test_base)
    return (
        {
            "stable_copy_baseline": np.asarray(test_offset, dtype=float),
            "stable_reference_topology_raw": np.asarray(
                test_offset + test.correction @ beta, dtype=float
            ),
        },
        {
            "feature_contract": contract,
            "baseline_intercept": float(baseline.intercept_[0]),
            "baseline_coefficients": baseline.coef_[0].tolist(),
            "topology_offset_intercept": 0.0,
            "topology_offset_feature_order": list(CORRECTION_FEATURE_NAMES),
            "topology_offset_coefficients": beta.tolist(),
            "species_equal_weight_totals": {
                dataset.name.split(":", 1)[0]: float(
                    weights[
                        sum(len(prior.labels) for prior in train[:index]) : sum(
                            len(prior.labels) for prior in train[: index + 1]
                        )
                    ].sum()
                )
                for index, dataset in enumerate(train)
            },
            "sample_weight_mean": float(weights.mean()),
        },
    )


def apply_guard(
    dataset: Dataset, baseline: np.ndarray, primary: np.ndarray
) -> tuple[np.ndarray, dict[str, Any], np.ndarray]:
    computation = compute_conflict_winner_guard(
        digests=dataset.digests,
        conflict_by_digest={
            digest: str(dataset.conflict_sets[index])
            for index, digest in enumerate(dataset.digests)
        },
        baseline_scores={
            digest: float(baseline[index])
            for index, digest in enumerate(dataset.digests)
        },
        primary_scores={
            digest: float(primary[index])
            for index, digest in enumerate(dataset.digests)
        },
    )
    guarded = np.asarray(
        [computation["scores"][digest] for digest in dataset.digests], dtype=float
    )
    applied = np.asarray(
        [digest in computation["guarded_digests"] for digest in dataset.digests],
        dtype=np.uint8,
    )
    audit = {
        key: value
        for key, value in computation.items()
        if key
        in {
            "baseline_winner_mapping_sha256",
            "primary_winner_mapping_sha256",
            "guard_winner_mapping_sha256",
            "winner_mapping_sha256",
            "winner_mismatch_count",
        }
    }
    audit.update(
        {
            "conflict_sets": len(computation["conflicts"]),
            "guarded_sets": len(computation["guarded_sets"]),
            "guarded_candidates": len(computation["guarded_digests"]),
            "automatic_approvals": 0,
        }
    )
    return guarded, audit, applied


def grouped_oof(
    dataset: Dataset,
) -> tuple[dict[str, np.ndarray], list[dict[str, Any]], dict[str, Any], np.ndarray]:
    baseline_sum = np.zeros(len(dataset.labels), dtype=float)
    primary_sum = np.zeros(len(dataset.labels), dtype=float)
    partitions: list[dict[str, Any]] = []
    seen: set[tuple[tuple[str, ...], ...]] = set()
    placeholder = np.zeros((len(dataset.labels), 1), dtype=float)
    for seed in OOF_SEEDS:
        baseline = np.full(len(dataset.labels), np.nan)
        primary = np.full(len(dataset.labels), np.nan)
        splitter = GroupKFold(n_splits=FOLDS, shuffle=True, random_state=seed)
        splits = list(splitter.split(placeholder, dataset.labels, dataset.groups))
        signature = tuple(
            sorted(
                tuple(sorted(set(str(value) for value in dataset.groups[test])))
                for _, test in splits
            )
        )
        if signature in seen:
            raise ValueError(f"Repeated chromosome partition for {dataset.name}")
        seen.add(signature)
        fold_records: list[dict[str, Any]] = []
        for fold, (train_indexes, test_indexes) in enumerate(splits, start=1):
            train_groups = set(dataset.groups[train_indexes])
            test_groups = set(dataset.groups[test_indexes])
            if train_groups & test_groups:
                raise AssertionError("Target-chromosome leakage in stable OOF")
            if len(np.unique(dataset.labels[train_indexes])) != 2:
                raise ValueError("Stable OOF train fold lacks a label class")
            scores, artifact = fit_predict(
                [subset(dataset, train_indexes, f"train:{seed}:{fold}")],
                subset(dataset, test_indexes, f"test:{seed}:{fold}"),
            )
            baseline[test_indexes] = scores["stable_copy_baseline"]
            primary[test_indexes] = scores["stable_reference_topology_raw"]
            fold_records.append(
                {
                    "fold": fold,
                    "train_groups": sorted(str(value) for value in train_groups),
                    "test_groups": sorted(str(value) for value in test_groups),
                    "model_artifact": artifact,
                    "model_sha256": sha256_text(
                        json.dumps(artifact, sort_keys=True, separators=(",", ":"))
                    ),
                }
            )
        if np.isnan(baseline).any() or np.isnan(primary).any():
            raise AssertionError("Stable OOF predictions are incomplete")
        baseline_sum += baseline
        primary_sum += primary
        partitions.append({"seed": seed, "folds": fold_records})
    mean_baseline = baseline_sum / len(OOF_SEEDS)
    mean_primary = primary_sum / len(OOF_SEEDS)
    guarded, audit, applied = apply_guard(dataset, mean_baseline, mean_primary)
    return (
        {
            "stable_copy_baseline": mean_baseline,
            "stable_reference_topology_raw": mean_primary,
            "stable_reference_topology_guarded": guarded,
        },
        partitions,
        audit,
        applied,
    )


def metrics(labels: np.ndarray, scores: np.ndarray) -> dict[str, float]:
    return {
        "average_precision": float(average_precision_score(labels, scores)),
        "roc_auc": float(roc_auc_score(labels, scores)),
    }


def review_metrics(
    dataset: Dataset, scores: np.ndarray
) -> dict[str, dict[str, Any]]:
    budgets = [
        (f"top_{100 * fraction:g}pct", max(1, math.ceil(len(scores) * fraction)))
        for fraction in FRACTION_BUDGETS
    ]
    budgets.extend((f"top_{count}", min(len(scores), count)) for count in FIXED_BUDGETS)
    output: dict[str, dict[str, Any]] = {}
    ordering = sorted(
        range(len(scores)), key=lambda index: (-scores[index], dataset.digests[index])
    )
    for name, count in budgets:
        selected = ordering[:count]
        positives = int(dataset.labels[selected].sum())
        selected_text = "".join(f"{dataset.digests[index]}\n" for index in selected)
        output[name] = {
            "reviewed": count,
            "true_positive": positives,
            "precision": positives / count,
            "positive_candidate_recall": positives / int(dataset.labels.sum()),
            "selection_digest_sha256": sha256_text(selected_text),
        }
    return output


def _weighted_ap_contract(
    labels: np.ndarray, scores: np.ndarray, group_codes: np.ndarray, group_count: int
) -> tuple[np.ndarray, np.ndarray]:
    order = np.argsort(scores, kind="mergesort")[::-1]
    ordered_scores = scores[order]
    ordered_labels = labels[order]
    ordered_groups = group_codes[order]
    thresholds = np.concatenate(
        (np.flatnonzero(np.diff(ordered_scores)), [len(scores) - 1])
    )
    total = np.zeros((len(scores), group_count), dtype=np.int8)
    positive = np.zeros((len(scores), group_count), dtype=np.int8)
    indexes = np.arange(len(scores))
    total[indexes, ordered_groups] = 1
    positive[indexes, ordered_groups] = ordered_labels
    return np.cumsum(positive, axis=0)[thresholds], np.cumsum(total, axis=0)[thresholds]


def _weighted_ap_batch(
    counts: np.ndarray, cumulative_positive: np.ndarray, cumulative_total: np.ndarray
) -> np.ndarray:
    true_at = counts @ cumulative_positive.T
    total_at = counts @ cumulative_total.T
    precision = np.divide(
        true_at,
        total_at,
        out=np.zeros_like(true_at, dtype=float),
        where=total_at != 0,
    )
    total_positive = true_at[:, -1]
    recall = np.divide(
        true_at,
        total_positive[:, None],
        out=np.zeros_like(true_at, dtype=float),
        where=total_positive[:, None] != 0,
    )
    increments = np.diff(np.column_stack((np.zeros(len(counts)), recall)), axis=1)
    return np.sum(increments * precision, axis=1)


def group_bootstrap_delta(
    dataset: Dataset,
    primary: np.ndarray,
    baseline: np.ndarray,
    *,
    replicates: int,
    seed: int,
) -> dict[str, Any]:
    unique_groups, codes = np.unique(dataset.groups, return_inverse=True)
    primary_contract = _weighted_ap_contract(
        dataset.labels, primary, codes, len(unique_groups)
    )
    baseline_contract = _weighted_ap_contract(
        dataset.labels, baseline, codes, len(unique_groups)
    )
    positives = np.bincount(codes, weights=dataset.labels, minlength=len(unique_groups))
    negatives = np.bincount(
        codes, weights=1 - dataset.labels, minlength=len(unique_groups)
    )
    rng = np.random.default_rng(seed)
    deltas: list[float] = []
    probability = np.full(len(unique_groups), 1 / len(unique_groups))
    for start in range(0, replicates, 64):
        rows = min(64, replicates - start)
        counts = rng.multinomial(len(unique_groups), probability, size=rows)
        primary_ap = _weighted_ap_batch(counts, *primary_contract)
        baseline_ap = _weighted_ap_batch(counts, *baseline_contract)
        valid = (counts @ positives > 0) & (counts @ negatives > 0)
        deltas.extend((primary_ap[valid] - baseline_ap[valid]).tolist())
    lower, upper = np.quantile(np.asarray(deltas), (0.025, 0.975))
    return {
        "observed_delta": float(
            average_precision_score(dataset.labels, primary)
            - average_precision_score(dataset.labels, baseline)
        ),
        "ci_lower": float(lower),
        "ci_upper": float(upper),
        "replicates_requested": replicates,
        "replicates_valid": len(deltas),
        "seed": seed,
        "groups": len(unique_groups),
    }


def evaluation_record(
    dataset: Dataset,
    scores: dict[str, np.ndarray],
    guard_audit: dict[str, Any],
) -> dict[str, Any]:
    return {
        "counts": {
            "candidates": len(dataset.labels),
            "positives": int(dataset.labels.sum()),
            "negatives": int(len(dataset.labels) - dataset.labels.sum()),
            "topology_available": int(dataset.topology_available.sum()),
            "topology_positive": int(
                dataset.labels[dataset.topology_available == 1].sum()
            ),
            "topology_negative": int(
                (1 - dataset.labels[dataset.topology_available == 1]).sum()
            ),
        },
        "metrics": {name: metrics(dataset.labels, value) for name, value in scores.items()},
        "review_budgets": {
            name: review_metrics(dataset, value) for name, value in scores.items()
        },
        "guard_audit": guard_audit,
    }


def _write_predictions(
    path: Path,
    records: Sequence[
        tuple[str, str, Dataset, dict[str, np.ndarray], np.ndarray]
    ],
) -> None:
    fields = (
        "dataset",
        "evaluation",
        "candidate_digest",
        "seqid",
        "label_exact_cds",
        "conflict_set_digest",
        "topology_available",
        "stable_copy_baseline",
        "stable_reference_topology_raw",
        "stable_reference_topology_guarded",
        "conflict_guard_applied",
        "automatic_approval",
    )
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for dataset_name, evaluation, dataset, scores, applied in records:
            for index, digest in enumerate(dataset.digests):
                writer.writerow(
                    {
                        "dataset": dataset_name,
                        "evaluation": evaluation,
                        "candidate_digest": digest,
                        "seqid": dataset.groups[index],
                        "label_exact_cds": int(dataset.labels[index]),
                        "conflict_set_digest": dataset.conflict_sets[index],
                        "topology_available": int(dataset.topology_available[index]),
                        "stable_copy_baseline": format(
                            scores["stable_copy_baseline"][index], ".17g"
                        ),
                        "stable_reference_topology_raw": format(
                            scores["stable_reference_topology_raw"][index], ".17g"
                        ),
                        "stable_reference_topology_guarded": format(
                            scores["stable_reference_topology_guarded"][index], ".17g"
                        ),
                        "conflict_guard_applied": int(applied[index]),
                        "automatic_approval": 0,
                    }
                )


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Evaluate the frozen candidate-independent stable ranker v0.9"
    )
    parser.add_argument("--dataset", action="append", required=True)
    parser.add_argument("--protocol", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()
    protocol = Path(args.protocol)
    if not protocol.is_file() or protocol.stat().st_size == 0:
        raise FileNotFoundError("Missing stable v0.9 development protocol")
    datasets = [parse_dataset(value) for value in args.dataset]
    if sorted(dataset.name for dataset in datasets) != ["actinidia", "populus"]:
        raise ValueError("Stable v0.9 requires exactly Actinidia and Populus")
    by_name = {dataset.name: dataset for dataset in datasets}
    output = Path(args.output_dir)
    working = Path(str(output) + ".working")
    if output.exists() or working.exists():
        raise FileExistsError("Refusing to overwrite stable v0.9 evaluation")
    working.mkdir(parents=True)

    oof_predictions: dict[str, dict[str, np.ndarray]] = {}
    oof_partitions: dict[str, Any] = {}
    oof_audits: dict[str, Any] = {}
    oof_applied: dict[str, np.ndarray] = {}
    oof_reports: dict[str, Any] = {}
    for name in ("actinidia", "populus"):
        dataset = by_name[name]
        scores, partitions, guard, applied = grouped_oof(dataset)
        bootstrap = group_bootstrap_delta(
            dataset,
            scores["stable_reference_topology_guarded"],
            scores["stable_copy_baseline"],
            replicates=BOOTSTRAP_REPLICATES,
            seed=BOOTSTRAP_SEEDS[name],
        )
        oof_predictions[name] = scores
        oof_partitions[name] = partitions
        oof_audits[name] = guard
        oof_applied[name] = applied
        oof_reports[name] = {
            **evaluation_record(dataset, scores, guard),
            "ap_delta_bootstrap": bootstrap,
        }

    transfer_reports: dict[str, Any] = {}
    transfer_predictions: dict[str, dict[str, np.ndarray]] = {}
    transfer_applied: dict[str, np.ndarray] = {}
    transfer_model_audits: dict[str, Any] = {}
    for train_name, test_name in (("actinidia", "populus"), ("populus", "actinidia")):
        test = by_name[test_name]
        raw, artifact = fit_predict([by_name[train_name]], test)
        guarded, guard, applied = apply_guard(
            test,
            raw["stable_copy_baseline"],
            raw["stable_reference_topology_raw"],
        )
        scores = {**raw, "stable_reference_topology_guarded": guarded}
        key = f"{train_name}_to_{test_name}"
        transfer_predictions[key] = scores
        transfer_applied[key] = applied
        transfer_model_audits[key] = artifact
        transfer_reports[key] = evaluation_record(test, scores, guard)

    _, pooled_artifact = fit_predict(datasets, by_name["actinidia"])
    pooled_model = {
        "schema_version": STABLE_REFERENCE_RANKER_SCHEMA_VERSION,
        "model_version": "v0.9-development-candidate-independent",
        "truth_access": True,
        "training_species": ["actinidia", "populus"],
        "fixed_C": FIXED_C,
        "solver": "lbfgs",
        "max_iter": 5000,
        "automatic_approval": False,
        "calibrated_probability": False,
        "artifact": pooled_artifact,
    }

    top1_non_decrease = {}
    for name, report in oof_reports.items():
        budgets = report["review_budgets"]
        top1_non_decrease[f"oof:{name}"] = (
            budgets["stable_reference_topology_guarded"]["top_1pct"]["true_positive"]
            >= budgets["stable_copy_baseline"]["top_1pct"]["true_positive"]
        )
    for key, report in transfer_reports.items():
        budgets = report["review_budgets"]
        top1_non_decrease[f"transfer:{key}"] = (
            budgets["stable_reference_topology_guarded"]["top_1pct"]["true_positive"]
            >= budgets["stable_copy_baseline"]["top_1pct"]["true_positive"]
        )
    gates = {
        "oof_actinidia_ap_ci_lower_positive": (
            oof_reports["actinidia"]["ap_delta_bootstrap"]["observed_delta"] > 0
            and oof_reports["actinidia"]["ap_delta_bootstrap"]["ci_lower"] > 0
        ),
        "oof_populus_ap_ci_lower_positive": (
            oof_reports["populus"]["ap_delta_bootstrap"]["observed_delta"] > 0
            and oof_reports["populus"]["ap_delta_bootstrap"]["ci_lower"] > 0
        ),
        "transfer_actinidia_to_populus_ap_delta_positive": (
            transfer_reports["actinidia_to_populus"]["metrics"]["stable_reference_topology_guarded"]["average_precision"]
            > transfer_reports["actinidia_to_populus"]["metrics"]["stable_copy_baseline"]["average_precision"]
        ),
        "transfer_populus_to_actinidia_ap_delta_positive": (
            transfer_reports["populus_to_actinidia"]["metrics"]["stable_reference_topology_guarded"]["average_precision"]
            > transfer_reports["populus_to_actinidia"]["metrics"]["stable_copy_baseline"]["average_precision"]
        ),
        "top1pct_non_decrease_all_oof_and_transfer": all(top1_non_decrease.values()),
        "winner_mismatch_zero_all": all(
            report["guard_audit"]["winner_mismatch_count"] == 0
            for report in [*oof_reports.values(), *transfer_reports.values()]
        ),
        "automatic_approvals_zero_all": all(
            report["guard_audit"]["automatic_approvals"] == 0
            for report in [*oof_reports.values(), *transfer_reports.values()]
        ),
    }
    retain = all(gates.values())
    evaluation = {
        "schema_version": EVALUATION_SCHEMA_VERSION,
        "code_commit": os.environ.get("PLOIDYPATCH_CODE_COMMIT", "unavailable"),
        "status": (
            "retain_v09_for_untouched_walnut"
            if retain
            else "retire_v09_ranker_keep_chain_workflow"
        ),
        "retention_gates": gates,
        "top1pct_gate_components": top1_non_decrease,
        "oof": oof_reports,
        "cross_species_transfer": transfer_reports,
        "evaluation_contract": {
            "fixed_C": FIXED_C,
            "folds": FOLDS,
            "oof_seeds": list(OOF_SEEDS),
            "bootstrap_seeds": BOOTSTRAP_SEEDS,
            "bootstrap_replicates": BOOTSTRAP_REPLICATES,
            "species_weight_policy": (
                "equal_total_mass_per_species_rescaled_to_mean_row_weight_one"
            ),
            "winner_tie_break": "descending_score_then_candidate_digest",
            "automatic_approval": False,
            "calibration": False,
        },
        "inputs": {
            dataset.name: {
                "rows": len(dataset.labels),
                "positives": int(dataset.labels.sum()),
                "files": [
                    {
                        "path": str(path),
                        "bytes": path.stat().st_size,
                        "sha256": _file_sha256(path),
                    }
                    for path in (*dataset.paths, *dataset.manifest_paths)
                ],
            }
            for dataset in datasets
        },
        "protocol": {
            "path": str(protocol),
            "bytes": protocol.stat().st_size,
            "sha256": _file_sha256(protocol),
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
    with (working / "evaluation.json").open("x", encoding="utf-8") as handle:
        json.dump(evaluation, handle, indent=2, sort_keys=True)
        handle.write("\n")
    with (working / "partitions.json").open("x", encoding="utf-8") as handle:
        json.dump(oof_partitions, handle, indent=2, sort_keys=True)
        handle.write("\n")
    with (working / "pooled_model.json").open("x", encoding="utf-8") as handle:
        json.dump(pooled_model, handle, indent=2, sort_keys=True)
        handle.write("\n")
    with (working / "transfer_models.json").open("x", encoding="utf-8") as handle:
        json.dump(transfer_model_audits, handle, indent=2, sort_keys=True)
        handle.write("\n")
    prediction_records = []
    for name in ("actinidia", "populus"):
        prediction_records.append(
            (name, "mean_five_seed_grouped_oof", by_name[name], oof_predictions[name], oof_applied[name])
        )
    for key, scores in transfer_predictions.items():
        test_name = key.rsplit("_to_", 1)[1]
        prediction_records.append(
            (test_name, key, by_name[test_name], scores, transfer_applied[key])
        )
    _write_predictions(working / "predictions.tsv", prediction_records)
    with (working / "input_manifest.tsv").open("x", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("dataset", "role", "path", "bytes", "sha256"))
        for dataset in datasets:
            for role, path in zip(
                ("copy_features", "topology_features", "labels", "pool_decisions", "copy_manifest", "topology_manifest"),
                (*dataset.paths, *dataset.manifest_paths),
                strict=True,
            ):
                writer.writerow((dataset.name, role, str(path), path.stat().st_size, _file_sha256(path)))
    with (working / "run_contract.tsv").open("x", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("field", "value"))
        writer.writerow(("status", evaluation["status"]))
        writer.writerow(("protocol_sha256", _file_sha256(protocol)))
        writer.writerow(("automatic_approval", "false"))
        writer.writerow(("bootstrap_replicates", BOOTSTRAP_REPLICATES))
        writer.writerow(("all_retention_gates_pass", str(retain).lower()))
    write_sha256sums(working)
    verify_sha256sums(working, ignore_checksum_file=True)
    os.replace(working, output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
