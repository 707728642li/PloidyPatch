#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import platform
from collections import defaultdict
from pathlib import Path
from typing import Any, Sequence

import numpy as np
import scipy
import sklearn
from sklearn.metrics import average_precision_score, roc_auc_score

from ploidypatch import __version__ as ploidypatch_version


SCHEMA_VERSION = "ploidypatch.apple_external_ranker_evaluation.v1"
H1_BOOTSTRAP_SEED = 20260902
H2_BOOTSTRAP_SEED = 20260901
BOOTSTRAP_REPLICATES = 20_000
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


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Reveal frozen apple v0.3 endpoints exactly once"
    )
    parser.add_argument("--scores", required=True)
    parser.add_argument("--labels", required=True)
    parser.add_argument("--pool-decisions", required=True)
    parser.add_argument("--primary-pool-score", required=True)
    parser.add_argument("--legacy-pool-score", required=True)
    parser.add_argument("--evaluability", required=True)
    parser.add_argument("--policy", required=True)
    parser.add_argument("--execution-policy", required=True)
    parser.add_argument("--protocol", required=True)
    parser.add_argument("--secondary-score", action="append", default=[])
    parser.add_argument("--v02-scores")
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    paths = {
        name: Path(value)
        for name, value in {
            "scores": args.scores,
            "labels": args.labels,
            "pool_decisions": args.pool_decisions,
            "primary_pool_score": args.primary_pool_score,
            "legacy_pool_score": args.legacy_pool_score,
            "evaluability": args.evaluability,
            "policy": args.policy,
            "execution_policy": args.execution_policy,
            "protocol": args.protocol,
        }.items()
    }
    secondary_paths = parse_named_paths(args.secondary_score)
    for name, path in secondary_paths.items():
        paths[f"secondary_score:{name}"] = path
    if args.v02_scores:
        paths["v02_scores"] = Path(args.v02_scores)
    for role, path in paths.items():
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(f"Missing evaluator input {role}: {path}")

    policy = read_policy(paths["policy"])
    if policy.get("policy_id") != "ploidypatch_apple_external_validation_v0.3":
        raise ValueError("Wrong apple external validation policy")
    if int(policy["H2_bootstrap_seed"]) != H2_BOOTSTRAP_SEED:
        raise ValueError("H2 bootstrap seed differs from implementation")
    if policy.get("automatic_copy_addition_approval") != "false":
        raise ValueError("Apple policy must prohibit automatic approval")
    execution_policy = read_policy(paths["execution_policy"])
    if int(execution_policy["H1_bootstrap_seed"]) != H1_BOOTSTRAP_SEED:
        raise ValueError("H1 bootstrap seed differs from execution freeze")
    if int(execution_policy["H2_bootstrap_seed"]) != H2_BOOTSTRAP_SEED:
        raise ValueError("H2 bootstrap seed differs from execution freeze")
    if int(execution_policy["bootstrap_replicates"]) != BOOTSTRAP_REPLICATES:
        raise ValueError("Bootstrap count differs from execution freeze")

    score_fields, score_rows = read_tsv(paths["scores"], ("candidate_digest",))
    label_fields, label_rows = read_tsv(paths["labels"], ("candidate_digest",))
    _, decision_rows = read_tsv(
        paths["pool_decisions"], ("candidate_digest", "consensus_digest")
    )
    required_scores = {
        "candidate_digest",
        "seqid",
        "support_method_count",
        "v03_baseline_logit",
        "v03_primary_rank_score",
        "v03_topology_available",
        "v03_automatic_approval",
    }
    if not required_scores <= set(score_fields) or "label_exact_cds" not in label_fields:
        raise ValueError("Apple score or label table lacks frozen contract columns")
    if set(score_rows) != set(label_rows) or set(score_rows) != set(decision_rows):
        raise ValueError("Apple score, label and conflict candidate universes differ")

    score_manifest = load_json(Path(str(paths["scores"]) + ".manifest.json"))
    label_manifest = load_json(Path(str(paths["labels"]) + ".manifest.json"))
    if score_manifest.get("truth_access") is not False:
        raise ValueError("Score manifest does not prove truth-blind generation")
    if label_manifest.get("evaluator_only") is not True:
        raise ValueError("Label manifest is not evaluator-only")
    if score_manifest.get("outputs", {}).get("scores", {}).get("sha256") != sha256(
        paths["scores"]
    ):
        raise ValueError("Score checksum disagrees with manifest")
    if label_manifest.get("outputs", {}).get("labels", {}).get("sha256") != sha256(
        paths["labels"]
    ):
        raise ValueError("Label checksum disagrees with manifest")

    digests = sorted(score_rows)
    labels = np.asarray(
        [int(label_rows[digest]["label_exact_cds"]) for digest in digests],
        dtype=np.uint8,
    )
    if not 0 < int(labels.sum()) < len(labels):
        raise ValueError("Apple candidate universe lacks both classes")
    primary = np.asarray(
        [float(score_rows[digest]["v03_primary_rank_score"]) for digest in digests]
    )
    baseline = np.asarray(
        [float(score_rows[digest]["v03_baseline_logit"]) for digest in digests]
    )
    groups = np.asarray([score_rows[digest]["seqid"] for digest in digests])
    topology = np.asarray(
        [int(score_rows[digest]["v03_topology_available"]) for digest in digests],
        dtype=np.uint8,
    )
    if any(int(score_rows[digest]["v03_automatic_approval"]) for digest in digests):
        raise ValueError("Frozen apple scorer attempted automatic approval")
    conflicts = [decision_rows[digest].get("conflict_set_digest", "") for digest in digests]

    comparator_scores: dict[str, np.ndarray] = {
        "method_support_count": np.asarray(
            [float(score_rows[digest]["support_method_count"]) for digest in digests]
        ),
        "max_method_quality": np.asarray(
            [
                max(
                    optional_float(score_rows[digest], "miniprot_identity"),
                    optional_float(score_rows[digest], "gemoma_pAA"),
                    optional_float(score_rows[digest], "lifton_protein_identity"),
                )
                for digest in digests
            ]
        ),
    }
    if "v02_scores" in paths:
        _, v02_rows = read_tsv(paths["v02_scores"], ("candidate_digest",))
        if set(v02_rows) != set(digests):
            raise ValueError("v0.2 comparator candidate universe differs")
        comparator_scores["v02_topology"] = np.asarray(
            [float(v02_rows[digest]["homeolog_topology_rank_score"]) for digest in digests]
        )

    primary_pool_score = load_json(paths["primary_pool_score"])
    legacy_pool_score = load_json(paths["legacy_pool_score"])
    secondary_scores = {
        name: load_json(path) for name, path in secondary_paths.items()
    }
    evaluability = load_json(paths["evaluability"])
    h1, h1_deltas = event_bootstrap_delta(
        primary_pool_score,
        legacy_pool_score,
        replicates=BOOTSTRAP_REPLICATES,
        seed=H1_BOOTSTRAP_SEED,
    )
    if h1["events"] != evaluability.get("events"):
        raise ValueError("H1 event universe differs from evaluability report")
    h2, h2_deltas = chromosome_bootstrap_delta(
        labels,
        primary,
        baseline,
        groups,
        replicates=BOOTSTRAP_REPLICATES,
        seed=H2_BOOTSTRAP_SEED,
    )

    ranking_metrics = {
        "v03_primary": metrics(labels, primary),
        "v03_baseline": metrics(labels, baseline),
        **{
            name: metrics(labels, values)
            for name, values in comparator_scores.items()
        },
    }
    review = {
        "v03_primary": review_metrics(labels, primary, digests),
        "v03_baseline": review_metrics(labels, baseline, digests),
        **{
            name: review_metrics(labels, values, digests)
            for name, values in comparator_scores.items()
        },
    }
    conflict = {
        "v03_primary": conflict_metrics(labels, primary, digests, conflicts),
        "v03_baseline": conflict_metrics(labels, baseline, digests, conflicts),
    }
    topology_positive = int(labels[topology == 1].sum())
    topology_positive_coverage = topology_positive / int(labels.sum())

    h1_pass = h1["observed_delta"] > 0 and h1["ci_lower"] > 0
    h2_numerical_pass = h2["observed_delta"] > 0 and h2["ci_lower"] > 0
    h2["fixed_sequence_status"] = (
        "confirmatory_tested"
        if h1_pass
        else "descriptive_not_tested_due_to_H1_failure"
    )
    conflict_primary = conflict["v03_primary"]["top1_accuracy"]
    conflict_baseline = conflict["v03_baseline"]["top1_accuracy"]
    conflict_gate = (
        conflict_primary is not None
        and conflict_baseline is not None
        and conflict_primary >= conflict_baseline
    )
    review_gate = (
        review["v03_primary"]["top_1pct"]["true_positive"]
        >= review["v03_baseline"]["top_1pct"]["true_positive"]
    )
    topology_gate = topology_positive_coverage >= float(
        policy["minimum_topology_coverage_among_positive_candidates"]
    )
    all_pool_scores = {
        "primary_union": primary_pool_score,
        "legacy_union": legacy_pool_score,
        **secondary_scores,
    }
    collateral_by_arm = {
        name: score_collateral_gate(score)
        for name, score in all_pool_scores.items()
    }
    collateral_gate = all(collateral_by_arm.values())
    formal_evaluable = evaluability.get("formal_evaluable") is True
    gates = {
        "formal_evaluable": formal_evaluable,
        "H1_chain_preserving_ceiling": h1_pass,
        "H2_ranker_AP_numerical": h2_numerical_pass,
        "H2_tested_in_fixed_sequence": h1_pass,
        "conflict_top1_noninferiority": conflict_gate,
        "top_1pct_review_noninferiority": review_gate,
        "topology_positive_coverage": topology_gate,
        "zero_collateral_loss": collateral_gate,
        "automatic_approval_absent": True,
    }
    confirmatory_pass = all(gates.values())

    report: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "generator": {
            "ploidypatch": ploidypatch_version,
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "scikit_learn": sklearn.__version__,
        },
        "evaluation_role": "untouched_confirmatory_external_species",
        "inputs": {
            role: {"path": str(path), "bytes": path.stat().st_size, "sha256": sha256(path)}
            for role, path in paths.items()
        },
        "counts": {
            "candidates": len(labels),
            "positive_exact_cds_candidates": int(labels.sum()),
            "negative_candidates": len(labels) - int(labels.sum()),
            "target_chromosomes": len(set(groups.tolist())),
            "topology_available_candidates": int(topology.sum()),
            "topology_positive_candidates": topology_positive,
            "topology_positive_coverage": topology_positive_coverage,
            "automatic_approved": 0,
        },
        "H1_chain_preserving_candidate_ceiling": h1,
        "H2_primary_minus_baseline": h2,
        "ranking_metrics": ranking_metrics,
        "review_budgets": review,
        "conflict_sets": conflict,
        "candidate_pool_scores": {
            name: summarize_pool_score(score)
            for name, score in all_pool_scores.items()
        },
        "collateral_gate_by_arm": collateral_by_arm,
        "event_evaluability": evaluability,
        "gates": gates,
        "confirmatory_pass": confirmatory_pass,
        "claim_boundary": {
            "automatic_approval": False,
            "calibrated_probability": False,
            "interpretation": "review_priority_only",
            "failure_retained_without_apple_retuning": True,
        },
    }

    output = Path(args.output_dir)
    partial = Path(str(output) + ".partial")
    if output.exists() or partial.exists():
        raise FileExistsError("Refusing to overwrite apple external evaluation")
    partial.mkdir(parents=True)
    with (partial / "evaluation.json").open(
        "x", encoding="utf-8", newline=""
    ) as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
        handle.write("\n")
    candidate_fields = (
        "candidate_digest",
        "seqid",
        "label_exact_cds",
        "conflict_set_digest",
        "topology_available",
        "v03_primary_rank_score",
        "v03_baseline_logit",
        *comparator_scores,
    )
    with (partial / "candidates.tsv").open(
        "x", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(
            handle, fieldnames=candidate_fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for index, digest in enumerate(digests):
            writer.writerow(
                {
                    "candidate_digest": digest,
                    "seqid": groups[index],
                    "label_exact_cds": int(labels[index]),
                    "conflict_set_digest": conflicts[index],
                    "topology_available": int(topology[index]),
                    "v03_primary_rank_score": format(primary[index], ".17g"),
                    "v03_baseline_logit": format(baseline[index], ".17g"),
                    **{
                        name: format(values[index], ".17g")
                        for name, values in comparator_scores.items()
                    },
                }
            )
    with (partial / "bootstrap_deltas.tsv").open(
        "x", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("endpoint", "replicate", "delta"))
        for endpoint, values in (("H1", h1_deltas), ("H2", h2_deltas)):
            for index, value in enumerate(values, start=1):
                writer.writerow((endpoint, index, format(float(value), ".17g")))
    with (partial / "input_manifest.tsv").open(
        "x", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("role", "bytes", "sha256", "path"))
        for role, path in sorted(paths.items()):
            writer.writerow((role, path.stat().st_size, sha256(path), path))
    checksum_path = partial / "SHA256SUMS"
    with checksum_path.open("x", encoding="utf-8", newline="") as handle:
        for path in sorted(partial.iterdir()):
            if path == checksum_path:
                continue
            handle.write(f"{sha256(path)}  {path.name}\n")
    os.replace(partial, output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
