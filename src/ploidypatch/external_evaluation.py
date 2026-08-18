"""Frozen reusable statistics and I/O primitives for external evaluation."""
from __future__ import annotations

import csv
from collections import defaultdict
import hashlib
import json
from pathlib import Path
from typing import Any, Sequence

import numpy as np
from sklearn.metrics import average_precision_score, roc_auc_score


FRACTION_BUDGETS = (0.005, 0.01, 0.02)
ABSOLUTE_BUDGETS = (100, 250, 500)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return value


def read_policy(path: Path) -> dict[str, str]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ["field", "value"]:
            raise ValueError("Policy must have field and value columns")
        rows = list(reader)
    policy = {row["field"]: row["value"] for row in rows}
    if len(policy) != len(rows) or "" in policy:
        raise ValueError("Policy contains empty or duplicate fields")
    return policy


def read_tsv(path: Path, key_fields: Sequence[str]) -> tuple[list[str], dict[str, dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Missing TSV header: {path}")
        fields = list(reader.fieldnames)
        key = next((field for field in key_fields if field in fields), None)
        if key is None:
            raise ValueError(f"Missing candidate key {key_fields}: {path}")
        rows = list(reader)
    indexed = {row[key]: row for row in rows}
    if len(indexed) != len(rows) or "" in indexed:
        raise ValueError(f"Empty or duplicate candidate key: {path}")
    return fields, indexed


def parse_named_paths(values: Sequence[str]) -> dict[str, Path]:
    output: dict[str, Path] = {}
    for value in values:
        name, separator, path = value.partition("=")
        if not separator or not name or not path or name in output:
            raise ValueError("--secondary-score requires unique NAME=PATH values")
        output[name] = Path(path)
    return output


def metrics(labels: np.ndarray, scores: np.ndarray) -> dict[str, float]:
    if not 0 < int(labels.sum()) < len(labels):
        raise ValueError("Ranking metrics require both label classes")
    return {
        "average_precision": float(average_precision_score(labels, scores)),
        "roc_auc": float(roc_auc_score(labels, scores)),
    }


def weighted_ap_contract(
    labels: np.ndarray,
    scores: np.ndarray,
    group_codes: np.ndarray,
    group_count: int,
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
    return np.cumsum(positive, axis=0)[thresholds], np.cumsum(total, axis=0)[
        thresholds
    ]


def weighted_ap_batch(
    counts: np.ndarray,
    cumulative_positive: np.ndarray,
    cumulative_total: np.ndarray,
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


def chromosome_bootstrap_delta(
    labels: np.ndarray,
    primary: np.ndarray,
    baseline: np.ndarray,
    groups: np.ndarray,
    *,
    replicates: int,
    seed: int,
) -> tuple[dict[str, Any], np.ndarray]:
    unique_groups, group_codes = np.unique(groups, return_inverse=True)
    primary_contract = weighted_ap_contract(
        labels, primary, group_codes, len(unique_groups)
    )
    baseline_contract = weighted_ap_contract(
        labels, baseline, group_codes, len(unique_groups)
    )
    positive_by_group = np.bincount(
        group_codes, weights=labels, minlength=len(unique_groups)
    )
    negative_by_group = np.bincount(
        group_codes, weights=1 - labels, minlength=len(unique_groups)
    )
    rng = np.random.default_rng(seed)
    probabilities = np.full(len(unique_groups), 1 / len(unique_groups))
    values: list[float] = []
    for start in range(0, replicates, 64):
        rows = min(64, replicates - start)
        counts = rng.multinomial(len(unique_groups), probabilities, size=rows)
        primary_ap = weighted_ap_batch(counts, *primary_contract)
        baseline_ap = weighted_ap_batch(counts, *baseline_contract)
        valid = (counts @ positive_by_group > 0) & (
            counts @ negative_by_group > 0
        )
        values.extend((primary_ap[valid] - baseline_ap[valid]).tolist())
    deltas = np.asarray(values, dtype=float)
    if not len(deltas):
        raise ValueError("No valid chromosome bootstrap replicate")
    lower, upper = np.quantile(deltas, (0.025, 0.975))
    observed = float(
        average_precision_score(labels, primary)
        - average_precision_score(labels, baseline)
    )
    return (
        {
            "observed_delta": observed,
            "ci_lower": float(lower),
            "ci_upper": float(upper),
            "replicates_requested": replicates,
            "replicates_valid": len(deltas),
            "groups": unique_groups.tolist(),
            "seed": seed,
        },
        deltas,
    )


def event_bootstrap_delta(
    primary_score: dict[str, Any],
    legacy_score: dict[str, Any],
    *,
    replicates: int,
    seed: int,
) -> tuple[dict[str, Any], np.ndarray]:
    primary_events = {
        row["event_id"]: int(bool(row["complete_cds_chain_recovery"]))
        for row in primary_score.get("event_details", [])
    }
    legacy_events = {
        row["event_id"]: int(bool(row["complete_cds_chain_recovery"]))
        for row in legacy_score.get("event_details", [])
    }
    if not primary_events or set(primary_events) != set(legacy_events):
        raise ValueError("H1 score files have empty or different event universes")
    event_ids = sorted(primary_events)
    paired = np.asarray(
        [primary_events[event] - legacy_events[event] for event in event_ids],
        dtype=float,
    )
    rng = np.random.default_rng(seed)
    deltas = np.empty(replicates, dtype=float)
    for start in range(0, replicates, 256):
        rows = min(256, replicates - start)
        samples = rng.integers(0, len(paired), size=(rows, len(paired)))
        deltas[start : start + rows] = paired[samples].mean(axis=1)
    lower, upper = np.quantile(deltas, (0.025, 0.975))
    return (
        {
            "events": len(event_ids),
            "primary_recovered": sum(primary_events.values()),
            "legacy_recovered": sum(legacy_events.values()),
            "primary_recall": sum(primary_events.values()) / len(event_ids),
            "legacy_recall": sum(legacy_events.values()) / len(event_ids),
            "observed_delta": float(paired.mean()),
            "ci_lower": float(lower),
            "ci_upper": float(upper),
            "replicates": replicates,
            "seed": seed,
        },
        deltas,
    )


def review_metrics(
    labels: np.ndarray,
    scores: np.ndarray,
    digests: Sequence[str],
) -> dict[str, Any]:
    ordered = sorted(
        range(len(scores)), key=lambda index: (-scores[index], digests[index])
    )
    output: dict[str, Any] = {}
    for fraction in FRACTION_BUDGETS:
        count = max(1, int(np.ceil(len(scores) * fraction)))
        selected = ordered[:count]
        positives = int(labels[selected].sum())
        output[f"top_{fraction * 100:g}pct"] = {
            "reviewed": count,
            "true_positive": positives,
            "precision": positives / count,
            "positive_candidate_recall": positives / int(labels.sum()),
        }
    for requested in ABSOLUTE_BUDGETS:
        count = min(requested, len(scores))
        selected = ordered[:count]
        positives = int(labels[selected].sum())
        output[f"top_{requested}"] = {
            "reviewed": count,
            "true_positive": positives,
            "precision": positives / count,
            "positive_candidate_recall": positives / int(labels.sum()),
        }
    return output


def conflict_metrics(
    labels: np.ndarray,
    scores: np.ndarray,
    digests: Sequence[str],
    conflict_sets: Sequence[str],
) -> dict[str, Any]:
    grouped: dict[str, list[int]] = defaultdict(list)
    for index, conflict in enumerate(conflict_sets):
        if conflict:
            grouped[conflict].append(index)
    evaluable = [
        indexes
        for indexes in grouped.values()
        if len(indexes) >= 2 and int(labels[indexes].sum()) == 1
    ]
    top1 = 0
    reciprocal_ranks: list[float] = []
    for indexes in evaluable:
        ordered = sorted(
            indexes, key=lambda index: (-scores[index], digests[index])
        )
        rank = 1 + next(
            position
            for position, index in enumerate(ordered)
            if labels[index] == 1
        )
        top1 += int(rank == 1)
        reciprocal_ranks.append(1 / rank)
    return {
        "conflict_sets_total": len(grouped),
        "evaluable_exactly_one_positive": len(evaluable),
        "top1_correct": top1,
        "top1_accuracy": top1 / len(evaluable) if evaluable else None,
        "mean_reciprocal_rank": (
            float(np.mean(reciprocal_ranks)) if reciprocal_ranks else None
        ),
    }


def optional_float(row: dict[str, str], field: str) -> float:
    value = row.get(field, "")
    if value == "":
        return 0.0
    parsed = float(value)
    if not np.isfinite(parsed):
        raise ValueError(f"Non-finite comparator field {field}")
    return parsed


def score_collateral_gate(score: dict[str, Any]) -> bool:
    return (
        score.get("quality_gate", {}).get("grade") == "pass"
        and score.get("collateral_changes", {}).get(
            "baseline_transcript_structures_missing_from_candidate"
        )
        == 0
    )


def summarize_pool_score(score: dict[str, Any]) -> dict[str, Any]:
    strict = score.get("strict_cds_chain", {})
    recovery = score.get("event_recovery", {})
    background = score.get("background_subtraction", {})
    collateral = score.get("collateral_changes", {})
    quality = score.get("quality_gate", {})
    return {
        "events": recovery.get("events"),
        "complete_cds_chain_recovery": recovery.get(
            "complete_cds_chain_recovery"
        ),
        "complete_cds_chain_recall": recovery.get("complete_cds_chain_recall"),
        "strict_cds_precision": strict.get("precision"),
        "strict_cds_recall": strict.get("recall"),
        "strict_cds_f1": strict.get("f1"),
        "differential_candidate_cds_chains": background.get(
            "differential_candidate_cds_chains"
        ),
        "baseline_transcript_structures_missing_from_candidate": collateral.get(
            "baseline_transcript_structures_missing_from_candidate"
        ),
        "quality_grade": quality.get("grade"),
    }
