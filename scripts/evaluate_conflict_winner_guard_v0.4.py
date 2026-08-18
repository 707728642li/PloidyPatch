#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from pathlib import Path
from typing import Any

import numpy as np
from sklearn.metrics import average_precision_score, roc_auc_score

from ploidypatch.conflict_guard import compute_conflict_winner_guard


SCHEMA_VERSION = "ploidypatch.conflict_winner_guard_evaluation.v1"
DEVELOPMENT_ROLES = {"development", "seen_external_development"}
DEFAULT_REPLICATES = 20_000
DEFAULT_SEED = 20260808


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def parse_input(value: str) -> tuple[str, str, Path]:
    name, separator, payload = value.partition("=")
    role, role_separator, path_value = payload.partition(",")
    if (
        not separator
        or not role_separator
        or not name
        or role
        not in {"development", "seen_external_development", "retrospective_diagnostic"}
    ):
        raise ValueError("--input requires NAME=ROLE,PATH")
    path = Path(path_value)
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"Missing input for {name}: {path}")
    return name, role, path


def read_rows(name: str, path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Missing header: {path}")
        fields = set(reader.fieldnames)
        if "dataset" in fields:
            required = {
                "candidate_digest",
                "seqid",
                "label",
                "conflict_set",
                "baseline",
                "offset_support_conditioned",
            }
            if not required <= fields:
                raise ValueError(f"Missing development prediction fields: {path}")
            selected = [row for row in reader if row["dataset"] == name]
            rows = [
                {
                    "candidate_digest": row["candidate_digest"],
                    "seqid": row["seqid"],
                    "label": row["label"],
                    "conflict_set": row["conflict_set"],
                    "baseline": row["baseline"],
                    "primary": row["offset_support_conditioned"],
                }
                for row in selected
            ]
        else:
            required = {
                "candidate_digest",
                "seqid",
                "label_exact_cds",
                "conflict_set_digest",
                "v03_baseline_logit",
                "v03_primary_rank_score",
            }
            if not required <= fields:
                raise ValueError(f"Missing external prediction fields: {path}")
            rows = [
                {
                    "candidate_digest": row["candidate_digest"],
                    "seqid": row["seqid"],
                    "label": row["label_exact_cds"],
                    "conflict_set": row["conflict_set_digest"],
                    "baseline": row["v03_baseline_logit"],
                    "primary": row["v03_primary_rank_score"],
                }
                for row in reader
            ]
    if not rows:
        raise ValueError(f"No rows selected for {name}")
    digests = [row["candidate_digest"] for row in rows]
    if "" in digests or len(digests) != len(set(digests)):
        raise ValueError(f"Empty or duplicate candidate digest for {name}")
    return rows


def apply_production_guard(
    rows: list[dict[str, str]],
) -> tuple[
    np.ndarray,
    np.ndarray,
    np.ndarray,
    dict[str, list[int]],
    set[str],
    dict[str, Any],
]:
    digests = [row["candidate_digest"] for row in rows]
    baseline_by_digest = {
        row["candidate_digest"]: float(row["baseline"]) for row in rows
    }
    primary_by_digest = {
        row["candidate_digest"]: float(row["primary"]) for row in rows
    }
    computation = compute_conflict_winner_guard(
        digests=digests,
        conflict_by_digest={
            row["candidate_digest"]: row["conflict_set"] for row in rows
        },
        baseline_scores=baseline_by_digest,
        primary_scores=primary_by_digest,
    )
    index_by_digest = {digest: index for index, digest in enumerate(digests)}
    conflicts = {
        conflict: [index_by_digest[digest] for digest in members]
        for conflict, members in computation["conflicts"].items()
    }
    return (
        np.asarray([baseline_by_digest[digest] for digest in digests]),
        np.asarray([primary_by_digest[digest] for digest in digests]),
        np.asarray([computation["scores"][digest] for digest in digests]),
        conflicts,
        computation["guarded_sets"],
        computation,
    )


def ranking_metrics(labels: np.ndarray, scores: np.ndarray) -> dict[str, float]:
    return {
        "average_precision": float(average_precision_score(labels, scores)),
        "roc_auc": float(roc_auc_score(labels, scores)),
    }


def review_metrics(
    labels: np.ndarray, scores: np.ndarray, digests: list[str]
) -> dict[str, dict[str, float | int | str]]:
    output: dict[str, dict[str, float | int | str]] = {}
    budgets = {
        "top_0.5pct": max(1, int(np.ceil(len(scores) * 0.005))),
        "top_1pct": max(1, int(np.ceil(len(scores) * 0.01))),
        "top_2pct": max(1, int(np.ceil(len(scores) * 0.02))),
        "top_100": min(100, len(scores)),
        "top_250": min(250, len(scores)),
        "top_500": min(500, len(scores)),
    }
    ranking = sorted(range(len(scores)), key=lambda index: (-scores[index], digests[index]))
    for name, count in budgets.items():
        positives = int(labels[ranking[:count]].sum())
        selected_digests = sorted(digests[index] for index in ranking[:count])
        output[name] = {
            "reviewed": count,
            "true_positive": positives,
            "precision": positives / count,
            "positive_candidate_recall": positives / int(labels.sum()),
            "selected_digest_sha256": hashlib.sha256(
                "".join(f"{digest}\n" for digest in selected_digests).encode("utf-8")
            ).hexdigest(),
        }
    return output


def conflict_metrics(
    labels: np.ndarray,
    scores: np.ndarray,
    digests: list[str],
    conflicts: dict[str, list[int]],
) -> dict[str, float | int | None]:
    evaluable = [
        indexes
        for indexes in conflicts.values()
        if int(labels[indexes].sum()) == 1
    ]
    top1 = 0
    reciprocal_ranks: list[float] = []
    for indexes in evaluable:
        ordered = sorted(indexes, key=lambda index: (-scores[index], digests[index]))
        positive_rank = 1 + next(
            rank for rank, index in enumerate(ordered) if labels[index] == 1
        )
        top1 += int(positive_rank == 1)
        reciprocal_ranks.append(1 / positive_rank)
    return {
        "conflict_sets_total": len(conflicts),
        "evaluable_exactly_one_positive": len(evaluable),
        "top1_correct": top1,
        "top1_accuracy": top1 / len(evaluable) if evaluable else None,
        "mean_reciprocal_rank": float(np.mean(reciprocal_ranks))
        if reciprocal_ranks
        else None,
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
    return (
        np.cumsum(positive, axis=0)[thresholds],
        np.cumsum(total, axis=0)[thresholds],
    )


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
    left: np.ndarray,
    right: np.ndarray,
    groups: np.ndarray,
    *,
    replicates: int,
    seed: int,
) -> dict[str, Any]:
    unique_groups, group_codes = np.unique(groups, return_inverse=True)
    left_contract = weighted_ap_contract(labels, left, group_codes, len(unique_groups))
    right_contract = weighted_ap_contract(
        labels, right, group_codes, len(unique_groups)
    )
    positive_by_group = np.bincount(
        group_codes, weights=labels, minlength=len(unique_groups)
    )
    negative_by_group = np.bincount(
        group_codes, weights=1 - labels, minlength=len(unique_groups)
    )
    probabilities = np.full(len(unique_groups), 1 / len(unique_groups))
    rng = np.random.default_rng(seed)
    batches = [
        rng.multinomial(
            len(unique_groups),
            probabilities,
            size=min(128, replicates - start),
        )
        for start in range(0, replicates, 128)
    ]

    def evaluate_batch(counts: np.ndarray) -> list[float]:
        left_ap = weighted_ap_batch(counts, *left_contract)
        right_ap = weighted_ap_batch(counts, *right_contract)
        valid = (counts @ positive_by_group > 0) & (
            counts @ negative_by_group > 0
        )
        return (left_ap[valid] - right_ap[valid]).tolist()

    with ThreadPoolExecutor(max_workers=min(8, len(batches))) as executor:
        evaluated = executor.map(evaluate_batch, batches)
        deltas = [delta for batch in evaluated for delta in batch]
    lower, upper = np.quantile(deltas, (0.025, 0.975))
    return {
        "observed_delta": float(
            average_precision_score(labels, left)
            - average_precision_score(labels, right)
        ),
        "ci_lower": float(lower),
        "ci_upper": float(upper),
        "replicates_requested": replicates,
        "replicates_valid": len(deltas),
        "groups": unique_groups.tolist(),
        "seed": seed,
    }


def evaluate_dataset(
    name: str,
    role: str,
    path: Path,
    *,
    replicates: int,
    seed: int,
) -> tuple[dict[str, Any], list[dict[str, str]]]:
    rows = read_rows(name, path)
    digests = [row["candidate_digest"] for row in rows]
    labels = np.asarray([int(row["label"]) for row in rows], dtype=np.uint8)
    groups = np.asarray([row["seqid"] for row in rows])
    if not 0 < int(labels.sum()) < len(labels):
        raise ValueError(f"Both label classes are required for {name}")
    (
        baseline,
        primary,
        guard,
        conflicts,
        guarded_sets,
        guard_computation,
    ) = apply_production_guard(rows)
    methods = {"baseline": baseline, "v03_primary": primary, "v04_guard": guard}
    metrics = {key: ranking_metrics(labels, value) for key, value in methods.items()}
    conflict = {
        key: conflict_metrics(labels, value, digests, conflicts)
        for key, value in methods.items()
    }
    reviews = {
        key: review_metrics(labels, value, digests) for key, value in methods.items()
    }
    guard_delta = metrics["v04_guard"]["average_precision"] - metrics["baseline"][
        "average_precision"
    ]
    primary_delta = metrics["v03_primary"]["average_precision"] - metrics[
        "baseline"
    ]["average_precision"]
    retention = guard_delta / primary_delta if primary_delta > 0 else None
    report = {
        "role": role,
        "counts": {
            "candidates": len(rows),
            "positives": int(labels.sum()),
            "chromosomes": len(set(groups)),
            "conflict_sets": len(conflicts),
            "winner_disagreement_sets": len(guarded_sets),
            "guarded_candidates": sum(
                len(conflicts[conflict_name]) for conflict_name in guarded_sets
            ),
        },
        "metrics": metrics,
        "conflicts": conflict,
        "review_budgets": reviews,
        "guard_minus_baseline_chromosome_bootstrap": chromosome_bootstrap_delta(
            labels, guard, baseline, groups, replicates=replicates, seed=seed
        ),
        "guard_minus_v03_primary_chromosome_bootstrap": chromosome_bootstrap_delta(
            labels, guard, primary, groups, replicates=replicates, seed=seed + 1
        ),
        "v03_ap_gain_retained_fraction": retention,
        "conflict_winner_identity": {
            "baseline_mapping_sha256": guard_computation[
                "baseline_winner_mapping_sha256"
            ],
            "v03_primary_mapping_sha256": guard_computation[
                "primary_winner_mapping_sha256"
            ],
            "v04_guard_mapping_sha256": guard_computation[
                "guard_winner_mapping_sha256"
            ],
            "combined_audit_sha256": guard_computation["winner_mapping_sha256"],
            "mismatch_count": guard_computation["winner_mismatch_count"],
            "identical_to_baseline": (
                guard_computation["baseline_winner_mapping_sha256"]
                == guard_computation["guard_winner_mapping_sha256"]
                and guard_computation["winner_mismatch_count"] == 0
            ),
        },
        "input": {
            "path": str(path),
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
        },
    }
    output_rows = []
    for index, row in enumerate(rows):
        conflict_name = row["conflict_set"]
        guarded = conflict_name in guarded_sets
        output_rows.append(
            {
                "dataset": name,
                "role": role,
                **row,
                "v04_guard": format(guard[index], ".17g"),
                "v04_guard_applied": int(guarded),
                "v04_topology_abstained": int(guarded),
                "v04_automatic_approval": 0,
            }
        )
    return report, output_rows


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Evaluate the label-independent v0.4 conflict winner guard"
    )
    parser.add_argument("--input", action="append", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--replicates", type=int, default=DEFAULT_REPLICATES)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    args = parser.parse_args()
    if args.replicates < 100:
        parser.error("--replicates must be at least 100")
    inputs = [parse_input(value) for value in args.input]
    names = [name for name, _, _ in inputs]
    if len(names) != len(set(names)):
        raise ValueError("Dataset names must be unique")
    output = Path(args.output_dir)
    partial = Path(str(output) + ".partial")
    if output.exists() or partial.exists():
        raise FileExistsError("Refusing to overwrite v0.4 evaluation")
    partial.mkdir(parents=True)

    completed: dict[str, tuple[dict[str, Any], list[dict[str, str]]]] = {}
    with ProcessPoolExecutor(max_workers=len(inputs)) as executor:
        futures = {
            name: executor.submit(
                evaluate_dataset,
                name,
                role,
                path,
                replicates=args.replicates,
                seed=args.seed + 100 * offset,
            )
            for offset, (name, role, path) in enumerate(inputs)
        }
        for name, future in futures.items():
            completed[name] = future.result()
    reports = {name: completed[name][0] for name, _, _ in inputs}
    predictions = [
        row
        for name, _, _ in inputs
        for row in completed[name][1]
    ]

    development = {
        name: report
        for name, report in reports.items()
        if report["role"] in DEVELOPMENT_ROLES
    }
    gates: dict[str, Any] = {}
    for name, report in development.items():
        gates[name] = {
            "guard_ap_above_baseline_with_positive_ci": (
                report["guard_minus_baseline_chromosome_bootstrap"]["observed_delta"]
                > 0
                and report["guard_minus_baseline_chromosome_bootstrap"]["ci_lower"]
                > 0
            ),
            "retain_at_least_90pct_v03_ap_gain": (
                report["v03_ap_gain_retained_fraction"] is not None
                and report["v03_ap_gain_retained_fraction"] >= 0.9
            ),
            "top_1pct_yield_not_below_baseline": (
                report["review_budgets"]["v04_guard"]["top_1pct"]["true_positive"]
                >= report["review_budgets"]["baseline"]["top_1pct"]["true_positive"]
            ),
            "conflict_top1_identical_to_baseline": report[
                "conflict_winner_identity"
            ]["identical_to_baseline"],
        }
    all_gates_pass = bool(gates) and all(
        all(dataset_gates.values()) for dataset_gates in gates.values()
    )
    report = {
        "schema_version": SCHEMA_VERSION,
        "policy": {
            "name": "baseline_fallback_on_conflict_winner_disagreement",
            "truth_access": False,
            "automatic_approval": False,
            "winner_tie_break": "descending_score_then_candidate_digest",
            "development_gate_policy": {
                "positive_guard_minus_baseline_ap_ci_lower": True,
                "minimum_v03_ap_gain_retained_fraction": 0.9,
                "top_1pct_yield_not_below_baseline": True,
                "every_conflict_top1_identical_to_baseline": True,
                "retrospective_diagnostics_cannot_select_policy": True,
                "apple_labels_seen_and_cannot_be_external_confirmation_for_v04": True,
            },
        },
        "bootstrap": {"replicates": args.replicates, "base_seed": args.seed},
        "datasets": reports,
        "development_gates": gates,
        "all_development_gates_pass": all_gates_pass,
        "next_required_evidence": (
            "freeze_v04_then_confirm_on_a_new_untouched_plant_species"
            if all_gates_pass
            else "revise_or_reject_v04_before_any_new_external_confirmation"
        ),
        "code_commit": os.environ.get("PLOIDYPATCH_CODE_COMMIT", "unavailable"),
    }
    report_path = partial / "evaluation.json"
    with report_path.open("x", encoding="utf-8", newline="") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
        handle.write("\n")
    prediction_path = partial / "predictions.tsv"
    fields = [
        "dataset",
        "role",
        "candidate_digest",
        "seqid",
        "label",
        "conflict_set",
        "baseline",
        "primary",
        "v04_guard",
        "v04_guard_applied",
        "v04_topology_abstained",
        "v04_automatic_approval",
    ]
    with prediction_path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(predictions)
    checksum_path = partial / "SHA256SUMS"
    with checksum_path.open("x", encoding="utf-8", newline="") as handle:
        for path in sorted(item for item in partial.iterdir() if item != checksum_path):
            handle.write(f"{sha256(path)}  {path.name}\n")
    os.replace(partial, output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
