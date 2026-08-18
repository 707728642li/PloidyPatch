#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import platform
import sys
from pathlib import Path
from typing import Any

import numpy as np
import sklearn
from sklearn.metrics import average_precision_score, roc_auc_score


SCHEMA_VERSION = "ploidypatch.frozen_homeolog_ranker_evaluation.v2"
REVIEW_BUDGETS = (0.005, 0.01, 0.02)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def _load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return value


def _read_index(path: Path) -> tuple[list[str], dict[str, dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or "candidate_digest" not in reader.fieldnames:
            raise ValueError(f"Missing candidate_digest: {path}")
        fields = list(reader.fieldnames)
        rows = list(reader)
    indexed = {row["candidate_digest"]: row for row in rows}
    if len(indexed) != len(rows) or "" in indexed:
        raise ValueError(f"Empty or duplicate candidate digest: {path}")
    return fields, indexed


def _metrics(labels: np.ndarray, scores: np.ndarray) -> dict[str, float]:
    return {
        "average_precision": float(average_precision_score(labels, scores)),
        "roc_auc": float(roc_auc_score(labels, scores)),
    }


def _weighted_ap_contract(
    labels: np.ndarray, scores: np.ndarray, group_codes: np.ndarray, groups: int
) -> tuple[np.ndarray, np.ndarray]:
    order = np.argsort(scores, kind="mergesort")[::-1]
    ordered_scores = scores[order]
    ordered_labels = labels[order]
    ordered_groups = group_codes[order]
    thresholds = np.concatenate(
        (np.flatnonzero(np.diff(ordered_scores)), [len(scores) - 1])
    )
    total_by_group = np.zeros((len(scores), groups), dtype=np.int32)
    positive_by_group = np.zeros((len(scores), groups), dtype=np.int32)
    row_indexes = np.arange(len(scores))
    total_by_group[row_indexes, ordered_groups] = 1
    positive_by_group[row_indexes, ordered_groups] = ordered_labels
    cumulative_total = np.cumsum(total_by_group, axis=0)[thresholds]
    cumulative_positive = np.cumsum(positive_by_group, axis=0)[thresholds]
    return cumulative_positive, cumulative_total


def _weighted_ap_batch(
    counts: np.ndarray,
    cumulative_positive: np.ndarray,
    cumulative_total: np.ndarray,
) -> np.ndarray:
    true_at_threshold = counts @ cumulative_positive.T
    total_at_threshold = counts @ cumulative_total.T
    precision = np.divide(
        true_at_threshold,
        total_at_threshold,
        out=np.zeros_like(true_at_threshold, dtype=float),
        where=total_at_threshold != 0,
    )
    total_positive = true_at_threshold[:, -1]
    recall = np.divide(
        true_at_threshold,
        total_positive[:, None],
        out=np.zeros_like(true_at_threshold, dtype=float),
        where=total_positive[:, None] != 0,
    )
    increments = np.diff(
        np.column_stack((np.zeros(len(counts)), recall)), axis=1
    )
    return np.sum(increments * precision, axis=1)


def _group_bootstrap_delta(
    labels: np.ndarray,
    topology_scores: np.ndarray,
    baseline_scores: np.ndarray,
    groups: np.ndarray,
    *,
    replicates: int,
    seed: int,
) -> dict[str, Any]:
    unique_groups, group_codes = np.unique(groups, return_inverse=True)
    if len(unique_groups) < 2:
        raise ValueError("Group bootstrap requires at least two groups")
    observed = float(
        average_precision_score(labels, topology_scores)
        - average_precision_score(labels, baseline_scores)
    )
    topology_contract = _weighted_ap_contract(
        labels, topology_scores, group_codes, len(unique_groups)
    )
    baseline_contract = _weighted_ap_contract(
        labels, baseline_scores, group_codes, len(unique_groups)
    )
    positive_by_group = np.bincount(
        group_codes, weights=labels, minlength=len(unique_groups)
    )
    negative_by_group = np.bincount(
        group_codes, weights=1 - labels, minlength=len(unique_groups)
    )
    probabilities = np.full(len(unique_groups), 1 / len(unique_groups))
    rng = np.random.default_rng(seed)
    deltas: list[float] = []
    batch_size = 64
    for start in range(0, replicates, batch_size):
        rows = min(batch_size, replicates - start)
        counts = rng.multinomial(len(unique_groups), probabilities, size=rows)
        topology_ap = _weighted_ap_batch(counts, *topology_contract)
        baseline_ap = _weighted_ap_batch(counts, *baseline_contract)
        valid = (
            np.isfinite(topology_ap)
            & np.isfinite(baseline_ap)
            & (counts @ positive_by_group > 0)
            & (counts @ negative_by_group > 0)
        )
        deltas.extend((topology_ap[valid] - baseline_ap[valid]).tolist())
    if not deltas:
        raise ValueError("Every group bootstrap replicate had one label class")
    lower, upper = np.quantile(deltas, (0.025, 0.975))
    return {
        "observed_delta": observed,
        "ci_lower": float(lower),
        "ci_upper": float(upper),
        "replicates_requested": replicates,
        "replicates_valid": len(deltas),
        "resampling_unit": "candidate_target_seqid_chromosome",
        "group_count": len(unique_groups),
        "groups": unique_groups.tolist(),
        "seed": seed,
    }


def _review_budget(
    labels: np.ndarray,
    scores: np.ndarray,
    digests: list[str],
    truth_event_count: int,
    fraction: float,
) -> dict[str, Any]:
    count = max(1, math.ceil(len(labels) * fraction))
    order = sorted(
        range(len(labels)), key=lambda index: (-scores[index], digests[index])
    )[:count]
    true_positive = int(np.sum(labels[order]))
    return {
        "fraction": fraction,
        "reviewed_candidates": count,
        "true_positive_exact_cds": true_positive,
        "precision": true_positive / count,
        "recall_of_all_hidden_events": true_positive / truth_event_count,
        "recall_of_candidate_universe_positives": true_positive / int(np.sum(labels)),
    }


def evaluate(args: argparse.Namespace) -> dict[str, Any]:
    score_path = Path(args.scores)
    label_path = Path(args.labels)
    score_manifest_path = Path(str(score_path) + ".manifest.json")
    label_manifest_path = Path(str(label_path) + ".manifest.json")
    for path in (score_path, label_path, score_manifest_path, label_manifest_path):
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f"Missing frozen evaluator input: {path}")
    score_manifest = _load_json(score_manifest_path)
    label_manifest = _load_json(label_manifest_path)
    if (
        score_manifest.get("schema_version")
        != "ploidypatch.homeolog_copy_rank_scores.v2"
        or score_manifest.get("truth_access") is not False
        or score_manifest.get("outputs", {}).get("scores", {}).get("sha256")
        != _sha256(score_path)
    ):
        raise ValueError("Frozen rank-score manifest failed integrity gate")
    if (
        label_manifest.get("schema_version")
        != "ploidypatch.copy_candidate_labels.v1"
        or label_manifest.get("evaluator_only") is not True
        or label_manifest.get("outputs", {}).get("labels", {}).get("sha256")
        != _sha256(label_path)
    ):
        raise ValueError("Evaluator label manifest failed integrity gate")

    score_fields, scores = _read_index(score_path)
    label_fields, labels = _read_index(label_path)
    required_score_fields = {
        "homeolog_baseline_rank_score",
        "homeolog_topology_rank_score",
        "homeolog_topology_available",
    }
    if not required_score_fields <= set(score_fields):
        raise ValueError("Rank-score table lacks frozen model outputs")
    if not {"label_exact_cds", "seqid"} <= set(label_fields):
        raise ValueError("Label table lacks exact-CDS labels or seqids")
    if set(scores) != set(labels):
        raise ValueError("Score and label candidate universes differ")
    digests = sorted(scores)
    y = np.asarray([int(labels[digest]["label_exact_cds"]) for digest in digests])
    if len(set(y)) != 2:
        raise ValueError("Frozen candidate universe does not contain both classes")
    groups = np.asarray([labels[digest]["seqid"] for digest in digests])
    baseline = np.asarray(
        [float(scores[digest]["homeolog_baseline_rank_score"]) for digest in digests]
    )
    topology = np.asarray(
        [float(scores[digest]["homeolog_topology_rank_score"]) for digest in digests]
    )
    available = np.asarray(
        [int(scores[digest]["homeolog_topology_available"]) for digest in digests]
    )
    if np.any((available != 0) & (available != 1)):
        raise ValueError("Topology availability is not binary")
    positives = int(np.sum(y))
    topology_positive = int(np.sum(y * available))
    if positives > args.truth_event_count:
        raise ValueError("Candidate positives exceed declared hidden events")
    bootstrap = _group_bootstrap_delta(
        y,
        topology,
        baseline,
        groups,
        replicates=args.bootstrap_replicates,
        seed=args.seed,
    )
    topology_coverage = topology_positive / positives
    primary_gate = (
        bootstrap["observed_delta"] > 0 and bootstrap["ci_lower"] > 0
    )
    coverage_gate = topology_coverage >= args.minimum_topology_coverage
    report: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "software": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "scikit_learn": sklearn.__version__,
        },
        "inputs": {
            "scores": {"path": str(score_path), "sha256": _sha256(score_path)},
            "labels": {"path": str(label_path), "sha256": _sha256(label_path)},
            "score_manifest_sha256": _sha256(score_manifest_path),
            "label_manifest_sha256": _sha256(label_manifest_path),
        },
        "policy": {
            "model_refit_on_external_species": False,
            "threshold_tuning": False,
            "score_interpretation": "uncalibrated_review_rank_not_probability",
            "automatic_approval": False,
            "primary_metric": "average_precision_topology_minus_baseline",
            "primary_success_gate": "observed_delta_positive_and_95pct_CI_lower_positive",
            "minimum_topology_coverage": args.minimum_topology_coverage,
        },
        "counts": {
            "candidates": len(y),
            "hidden_events": args.truth_event_count,
            "positive_exact_cds_candidates": positives,
            "candidate_recall_ceiling": positives / args.truth_event_count,
            "negative_candidates": int(len(y) - positives),
            "topology_available": int(np.sum(available)),
            "topology_available_positive_exact_cds": topology_positive,
            "topology_coverage_among_positive_candidates": topology_coverage,
        },
        "metrics": {
            "baseline": _metrics(y, baseline),
            "topology": _metrics(y, topology),
            "topology_minus_baseline_group_bootstrap": bootstrap,
        },
        "review_budgets": {
            f"top_{fraction * 100:g}pct": {
                "baseline": _review_budget(
                    y, baseline, digests, args.truth_event_count, fraction
                ),
                "topology": _review_budget(
                    y, topology, digests, args.truth_event_count, fraction
                ),
            }
            for fraction in REVIEW_BUDGETS
        },
        "gates": {
            "primary_delta_ap": "pass" if primary_gate else "fail",
            "topology_positive_coverage": "pass" if coverage_gate else "fail",
            "overall_predeclared_v2_portability": (
                "pass" if primary_gate and coverage_gate else "fail"
            ),
        },
    }
    return report


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Evaluate frozen homeolog rank scores on an external truth set"
    )
    parser.add_argument("--scores", required=True)
    parser.add_argument("--labels", required=True)
    parser.add_argument("--truth-event-count", type=int, required=True)
    parser.add_argument("--bootstrap-replicates", type=int, default=20_000)
    parser.add_argument("--seed", type=int, default=20260829)
    parser.add_argument("--minimum-topology-coverage", type=float, default=0.70)
    parser.add_argument("--output-json", required=True)
    args = parser.parse_args()
    if args.truth_event_count < 1 or args.bootstrap_replicates < 100:
        parser.error("truth events and bootstrap replicates must be positive")
    if not 0 <= args.minimum_topology_coverage <= 1:
        parser.error("topology coverage must be within [0, 1]")
    output_path = Path(args.output_json)
    if output_path.exists():
        raise FileExistsError(output_path)
    report = evaluate(args)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("x", encoding="utf-8", newline="") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
