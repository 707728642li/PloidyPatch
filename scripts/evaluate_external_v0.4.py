#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import os
import platform
from pathlib import Path
from typing import Any, Sequence

import numpy as np
import scipy
import sklearn

from ploidypatch import __version__ as ploidypatch_version


COMMON_PATH = Path(__file__).with_name("evaluate_apple_external_v0.3.py")
SPEC = importlib.util.spec_from_file_location("ploidypatch_external_v03_common", COMMON_PATH)
if SPEC is None or SPEC.loader is None:
    raise ImportError(f"Cannot load external-evaluation primitives: {COMMON_PATH}")
COMMON = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(COMMON)

SCHEMA_VERSION = "ploidypatch.external_ranker_evaluation.v0.4"
POLICY_ID = "ploidypatch_populus_external_validation_v0.4"
GUARD_SCORE_SCHEMA = "ploidypatch.conflict_winner_guard_scores.v1"
POOL_SCHEMA = "ploidypatch.method_candidate_pool.v2"
LABEL_SCHEMA = "ploidypatch.external_candidate_labels.v1"
CUSTODY_SCHEMA = "ploidypatch.blind_run_custody.v1"
H1_BOOTSTRAP_SEED = 20261001
H2_BOOTSTRAP_SEED = 20261002
GUARD_V03_BOOTSTRAP_SEED = 20261003
BOOTSTRAP_REPLICATES = 20_000
MINIMUM_VALID_BOOTSTRAPS = 19_000
FRACTION_BUDGETS = (0.005, 0.01, 0.02)
ABSOLUTE_BUDGETS = (100, 250, 500)


def sha256(path: Path) -> str:
    return COMMON.sha256(path)


def load_json(path: Path) -> dict[str, Any]:
    return COMMON.load_json(path)


def verify_sha256sums(root: Path) -> dict[str, str]:
    checksum_path = root / "SHA256SUMS"
    if not checksum_path.is_file() or checksum_path.stat().st_size == 0:
        raise ValueError(f"Missing protocol SHA256SUMS: {root}")
    entries: dict[str, str] = {}
    for line_number, line in enumerate(
        checksum_path.read_text(encoding="utf-8").splitlines(), start=1
    ):
        digest, separator, name = line.partition("  ")
        path = root / name
        if (
            not separator
            or len(digest) != 64
            or not path.is_file()
            or sha256(path) != digest
        ):
            raise ValueError(f"Protocol checksum failure at line {line_number}: {name}")
        entries[name] = digest
    return entries


def selected_digest_sha256(digests: Sequence[str]) -> str:
    return hashlib.sha256(
        "".join(f"{digest}\n" for digest in digests).encode("utf-8")
    ).hexdigest()


def review_metrics(
    labels: np.ndarray,
    scores: np.ndarray,
    digests: Sequence[str],
) -> dict[str, Any]:
    ordered = sorted(
        range(len(scores)), key=lambda index: (-scores[index], digests[index])
    )
    output: dict[str, Any] = {}
    budgets = [
        (f"top_{fraction * 100:g}pct", max(1, int(np.ceil(len(scores) * fraction))))
        for fraction in FRACTION_BUDGETS
    ] + [
        (f"top_{requested}", min(requested, len(scores)))
        for requested in ABSOLUTE_BUDGETS
    ]
    for name, count in budgets:
        selected = ordered[:count]
        positives = int(labels[selected].sum())
        selected_digests = [digests[index] for index in selected]
        output[name] = {
            "reviewed": count,
            "true_positive": positives,
            "precision": positives / count,
            "positive_candidate_recall": positives / int(labels.sum()),
            "selected_digest_sha256": selected_digest_sha256(selected_digests),
        }
    return output


def read_inputs(
    values: dict[str, str], secondary: Sequence[str]
) -> tuple[dict[str, Path], dict[str, Path]]:
    paths = {name: Path(value) for name, value in values.items()}
    secondary_paths = COMMON.parse_named_paths(secondary)
    for name, path in secondary_paths.items():
        paths[f"secondary_score:{name}"] = path
    for role, path in paths.items():
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(f"Missing external evaluator input {role}: {path}")
    return paths, secondary_paths


def validate_custody(
    custody: dict[str, Any], scores: Path, score_manifest: Path
) -> None:
    forbidden_access = (
        "truth_mounted",
        "complete_target_annotation_mounted",
        "evaluator_references_mounted",
        "nas_data_mounted",
        "network_access",
    )
    if (
        custody.get("schema_version") != CUSTODY_SCHEMA
        or any(custody.get(field) is not False for field in forbidden_access)
        or custody.get("blind_outputs", {}).get("scores_sha256") != sha256(scores)
        or custody.get("blind_outputs", {}).get("score_manifest_sha256")
        != sha256(score_manifest)
        or not custody.get("runner_identity")
        or not custody.get("frozen_before_truth_reveal_at")
    ):
        raise ValueError("Blind-run custody contract is incomplete or violated")


def validate_evaluability(
    report: dict[str, Any], policy: dict[str, str]
) -> tuple[bool, dict[str, bool]]:
    sentinels = report.get("sentinels", {})
    events = int(report.get("events", 0))
    sentinel_gates = {
        "blind_noop_exact_recovery_zero": (
            sentinels.get("blind_noop_exact_recovery") == 0
        ),
        "complete_oracle_exact_recovery_one": (
            sentinels.get("complete_oracle_exact_recovery") == events
        ),
        "restoration_byte_identical": (
            sentinels.get("restoration_byte_identical") is True
        ),
        "blind_complete_genome_identical": (
            sentinels.get("blind_complete_genome_sha256_identical") is True
        ),
    }
    if not all(sentinel_gates.values()):
        raise ValueError(f"Sentinel violation invalidates run: {sentinel_gates}")
    complexity = report.get("complexity_bins", {})
    data_gates = {
        "minimum_events": events >= int(policy["minimum_formal_event_count"]),
        "minimum_target_chromosomes": int(report.get("target_chromosomes", 0))
        >= int(policy["minimum_target_chromosomes"]),
        "four_complexity_bins_present": len(complexity) == 4,
        "minimum_events_each_complexity_bin": bool(complexity)
        and min(int(value) for value in complexity.values())
        >= int(policy["minimum_events_per_complexity_bin"]),
    }
    return all(data_gates.values()), data_gates


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Reveal a preregistered untouched PloidyPatch v0.4 external test"
    )
    parser.add_argument("--scores", required=True)
    parser.add_argument("--labels", required=True)
    parser.add_argument("--pool-decisions", required=True)
    parser.add_argument("--pool-manifest", required=True)
    parser.add_argument("--primary-pool-score", required=True)
    parser.add_argument("--legacy-pool-score", required=True)
    parser.add_argument("--evaluability", required=True)
    parser.add_argument("--custody-manifest", required=True)
    parser.add_argument("--protocol-freeze", required=True)
    parser.add_argument("--composite-model-freeze", required=True)
    parser.add_argument("--policy", required=True)
    parser.add_argument("--secondary-score", action="append", default=[])
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    paths, secondary_paths = read_inputs(
        {
            "scores": args.scores,
            "labels": args.labels,
            "pool_decisions": args.pool_decisions,
            "pool_manifest": args.pool_manifest,
            "primary_pool_score": args.primary_pool_score,
            "legacy_pool_score": args.legacy_pool_score,
            "evaluability": args.evaluability,
            "custody_manifest": args.custody_manifest,
            "policy": args.policy,
        },
        args.secondary_score,
    )
    protocol_root = Path(args.protocol_freeze)
    composite_root = Path(args.composite_model_freeze)
    protocol_sums = verify_sha256sums(protocol_root)
    composite_sums = verify_sha256sums(composite_root)
    if paths["policy"].resolve().parent != protocol_root.resolve():
        raise ValueError("Policy must be the copy inside the protocol freeze")
    if protocol_sums.get(paths["policy"].name) != sha256(paths["policy"]):
        raise ValueError("Policy is not locked by protocol freeze")

    policy = COMMON.read_policy(paths["policy"])
    expected_policy = {
        "policy_id": POLICY_ID,
        "test_role": "untouched_confirmatory_external_species",
        "model_version": "PloidyPatch_ranker_v0.4",
        "automatic_copy_addition_approval": "false",
        "H1_bootstrap_seed": str(H1_BOOTSTRAP_SEED),
        "H2_bootstrap_seed": str(H2_BOOTSTRAP_SEED),
        "guard_v03_bootstrap_seed": str(GUARD_V03_BOOTSTRAP_SEED),
        "bootstrap_replicates": str(BOOTSTRAP_REPLICATES),
        "minimum_chromosome_bootstrap_valid_replicates": str(
            MINIMUM_VALID_BOOTSTRAPS
        ),
    }
    if any(policy.get(key) != value for key, value in expected_policy.items()):
        raise ValueError("External v0.4 policy identity or seeds differ")
    composite_manifest = load_json(composite_root / "composite_manifest.json")
    if (
        composite_manifest.get("schema_version")
        != "ploidypatch.composite_ranker.v0.4"
        or composite_manifest.get("automatic_approval") is not False
        or composite_sums.get("composite_manifest.json")
        != sha256(composite_root / "composite_manifest.json")
        or policy.get("composite_model_sha256sums_sha256")
        != sha256(composite_root / "SHA256SUMS")
    ):
        raise ValueError("Composite v0.4 model freeze differs from policy")

    score_manifest_path = Path(str(paths["scores"]) + ".manifest.json")
    label_manifest_path = Path(str(paths["labels"]) + ".manifest.json")
    if not score_manifest_path.is_file() or not label_manifest_path.is_file():
        raise ValueError("Score or label manifest is missing")
    score_manifest = load_json(score_manifest_path)
    label_manifest = load_json(label_manifest_path)
    paths["score_manifest"] = score_manifest_path
    paths["label_manifest"] = label_manifest_path
    pool_manifest = load_json(paths["pool_manifest"])
    if (
        score_manifest.get("schema_version") != GUARD_SCORE_SCHEMA
        or score_manifest.get("truth_access") is not False
        or score_manifest.get("outputs", {}).get("scores", {}).get("sha256")
        != sha256(paths["scores"])
        or score_manifest.get("winner_audit", {}).get("mismatch_count") != 0
        or score_manifest.get("winner_audit", {}).get("baseline_mapping_sha256")
        != score_manifest.get("winner_audit", {}).get("v04_guard_mapping_sha256")
        or score_manifest.get("inputs", {}).get("pool_decisions")
        != sha256(paths["pool_decisions"])
        or score_manifest.get("inputs", {}).get("pool_manifest")
        != sha256(paths["pool_manifest"])
    ):
        raise ValueError("v0.4 score manifest fails blind or winner-safety gate")
    if (
        pool_manifest.get("schema_version") != POOL_SCHEMA
        or pool_manifest.get("outputs", {}).get("decisions", {}).get("sha256")
        != sha256(paths["pool_decisions"])
    ):
        raise ValueError("Pool manifest fails schema or checksum gate")
    if (
        label_manifest.get("schema_version") != LABEL_SCHEMA
        or label_manifest.get("evaluator_only") is not True
        or label_manifest.get("outputs", {}).get("labels", {}).get("sha256")
        != sha256(paths["labels"])
        or label_manifest.get("blind_scores_sha256") != sha256(paths["scores"])
        or label_manifest.get("pool_manifest_sha256")
        != sha256(paths["pool_manifest"])
    ):
        raise ValueError("Evaluator labels do not bind to frozen blind artifacts")
    custody = load_json(paths["custody_manifest"])
    validate_custody(custody, paths["scores"], score_manifest_path)

    score_fields, score_rows = COMMON.read_tsv(paths["scores"], ("candidate_digest",))
    label_fields, label_rows = COMMON.read_tsv(paths["labels"], ("candidate_digest",))
    _, all_decision_rows = COMMON.read_tsv(
        paths["pool_decisions"], ("candidate_digest", "consensus_digest")
    )
    forbidden = [
        field for field in score_fields if "label" in field.lower() or "truth" in field.lower()
    ]
    required_scores = {
        "candidate_digest",
        "seqid",
        "support_method_count",
        "v03_baseline_logit",
        "v03_primary_rank_score",
        "v03_topology_available",
        "v04_primary_rank_score",
        "v04_conflict_guard_applied",
        "v04_topology_abstained",
        "v04_automatic_approval",
    }
    if forbidden or not required_scores <= set(score_fields) or "label_exact_cds" not in label_fields:
        raise ValueError("Score/label tables violate frozen evaluator columns")
    decision_rows = {
        digest: row
        for digest, row in all_decision_rows.items()
        if row.get("status") == "accepted"
    }
    if set(score_rows) != set(label_rows) or set(score_rows) != set(decision_rows):
        raise ValueError("Score, label and accepted-pool candidate universes differ")

    digests = sorted(score_rows)
    labels = np.asarray(
        [int(label_rows[digest]["label_exact_cds"]) for digest in digests],
        dtype=np.uint8,
    )
    if not 0 < int(labels.sum()) < len(labels):
        raise ValueError("External candidate labels require both classes")
    baseline = np.asarray(
        [float(score_rows[digest]["v03_baseline_logit"]) for digest in digests]
    )
    v03 = np.asarray(
        [float(score_rows[digest]["v03_primary_rank_score"]) for digest in digests]
    )
    v04 = np.asarray(
        [float(score_rows[digest]["v04_primary_rank_score"]) for digest in digests]
    )
    groups = np.asarray([score_rows[digest]["seqid"] for digest in digests])
    topology_available = np.asarray(
        [int(score_rows[digest]["v03_topology_available"]) for digest in digests],
        dtype=np.uint8,
    )
    topology_abstained = np.asarray(
        [int(score_rows[digest]["v04_topology_abstained"]) for digest in digests],
        dtype=np.uint8,
    )
    automatic = np.asarray(
        [int(score_rows[digest]["v04_automatic_approval"]) for digest in digests],
        dtype=np.uint8,
    )
    if not all(np.isfinite(values).all() for values in (baseline, v03, v04)):
        raise ValueError("Non-finite external score")
    if int(automatic.sum()):
        raise ValueError("Blind scorer attempted automatic approval")
    conflicts = [
        decision_rows[digest].get("conflict_set_digest", "") for digest in digests
    ]

    primary_pool_score = load_json(paths["primary_pool_score"])
    legacy_pool_score = load_json(paths["legacy_pool_score"])
    secondary_scores = {name: load_json(path) for name, path in secondary_paths.items()}
    evaluability = load_json(paths["evaluability"])
    formal_evaluable, data_gates = validate_evaluability(evaluability, policy)
    h1, h1_deltas = COMMON.event_bootstrap_delta(
        primary_pool_score,
        legacy_pool_score,
        replicates=BOOTSTRAP_REPLICATES,
        seed=H1_BOOTSTRAP_SEED,
    )
    if h1["events"] != int(evaluability.get("events", 0)):
        raise ValueError("H1 and evaluability event universes differ")
    h2, h2_deltas = COMMON.chromosome_bootstrap_delta(
        labels,
        v04,
        baseline,
        groups,
        replicates=BOOTSTRAP_REPLICATES,
        seed=H2_BOOTSTRAP_SEED,
    )
    guard_v03, guard_v03_deltas = COMMON.chromosome_bootstrap_delta(
        labels,
        v04,
        v03,
        groups,
        replicates=BOOTSTRAP_REPLICATES,
        seed=GUARD_V03_BOOTSTRAP_SEED,
    )

    support = np.asarray(
        [float(score_rows[digest]["support_method_count"]) for digest in digests]
    )
    quality = np.asarray(
        [
            max(
                COMMON.optional_float(score_rows[digest], "miniprot_identity"),
                COMMON.optional_float(score_rows[digest], "gemoma_pAA"),
                COMMON.optional_float(score_rows[digest], "lifton_protein_identity"),
            )
            for digest in digests
        ]
    )
    methods = {
        "baseline": baseline,
        "v03_primary": v03,
        "v04_guard": v04,
        "method_support_count": support,
        "max_method_quality": quality,
    }
    ranking_metrics = {
        name: COMMON.metrics(labels, scores) for name, scores in methods.items()
    }
    reviews = {name: review_metrics(labels, scores, digests) for name, scores in methods.items()}
    conflict = {
        name: COMMON.conflict_metrics(labels, scores, digests, conflicts)
        for name, scores in {"baseline": baseline, "v03_primary": v03, "v04_guard": v04}.items()
    }
    effective_topology = (topology_available == 1) & (topology_abstained == 0)
    topology_positive_coverage = int(labels[effective_topology].sum()) / int(labels.sum())
    v03_gain = (
        ranking_metrics["v03_primary"]["average_precision"]
        - ranking_metrics["baseline"]["average_precision"]
    )
    v04_gain = (
        ranking_metrics["v04_guard"]["average_precision"]
        - ranking_metrics["baseline"]["average_precision"]
    )
    retention = v04_gain / v03_gain if v03_gain > 0 else None
    all_pool_scores = {
        "primary_union": primary_pool_score,
        "legacy_union": legacy_pool_score,
        **secondary_scores,
    }
    collateral_by_arm = {
        name: COMMON.score_collateral_gate(score)
        for name, score in all_pool_scores.items()
    }

    h1_pass = h1["observed_delta"] > 0 and h1["ci_lower"] > 0
    h2_numerical_pass = h2["observed_delta"] > 0 and h2["ci_lower"] > 0
    h2_tested = formal_evaluable and h1_pass
    h2["fixed_sequence_status"] = (
        "confirmatory_tested" if h2_tested else "descriptive_not_tested"
    )
    bootstrap_gate = (
        int(h2["replicates_valid"]) >= MINIMUM_VALID_BOOTSTRAPS
        and int(guard_v03["replicates_valid"]) >= MINIMUM_VALID_BOOTSTRAPS
    )
    gates = {
        "formal_evaluable": formal_evaluable,
        "H1_chain_preserving_ceiling": h1_pass,
        "H2_guard_AP_numerical": h2_numerical_pass,
        "H2_tested_in_fixed_sequence": h2_tested,
        "conflict_winner_mapping_identical_to_baseline": True,
        "top_1pct_review_noninferiority": (
            reviews["v04_guard"]["top_1pct"]["true_positive"]
            >= reviews["baseline"]["top_1pct"]["true_positive"]
        ),
        "topology_positive_coverage": topology_positive_coverage
        >= float(policy["minimum_topology_coverage_among_positive_candidates"]),
        "retain_at_least_90pct_positive_v03_AP_gain": (
            retention is not None
            and retention >= float(policy["minimum_v03_AP_gain_retained_fraction"])
        ),
        "zero_collateral_loss_all_arms": all(collateral_by_arm.values()),
        "automatic_approval_absent": int(automatic.sum()) == 0,
        "bootstrap_valid_replicates": bootstrap_gate,
    }
    confirmatory_pass = all(gates.values())
    if not formal_evaluable:
        formal_outcome = "not_evaluable_without_rule_relaxation"
    elif confirmatory_pass:
        formal_outcome = "confirmatory_pass"
    else:
        formal_outcome = "formal_negative_external_result"

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
        "policy_id": POLICY_ID,
        "protocol_freeze_SHA256SUMS_sha256": sha256(protocol_root / "SHA256SUMS"),
        "composite_model_SHA256SUMS_sha256": sha256(composite_root / "SHA256SUMS"),
        "inputs": {
            role: {"path": str(path), "bytes": path.stat().st_size, "sha256": sha256(path)}
            for role, path in paths.items()
        },
        "counts": {
            "candidates": len(labels),
            "positive_exact_cds_candidates": int(labels.sum()),
            "negative_candidates": len(labels) - int(labels.sum()),
            "target_chromosomes": len(set(groups.tolist())),
            "effective_topology_positive_candidates": int(labels[effective_topology].sum()),
            "topology_positive_coverage": topology_positive_coverage,
            "automatic_approved": int(automatic.sum()),
        },
        "H1_chain_preserving_candidate_ceiling": h1,
        "H2_v04_guard_minus_baseline": h2,
        "descriptive_v04_guard_minus_v03": guard_v03,
        "v03_AP_gain_retained_fraction": retention,
        "ranking_metrics": ranking_metrics,
        "review_budgets": reviews,
        "conflict_sets": conflict,
        "winner_audit": score_manifest["winner_audit"],
        "candidate_pool_scores": {
            name: COMMON.summarize_pool_score(score)
            for name, score in all_pool_scores.items()
        },
        "collateral_gate_by_arm": collateral_by_arm,
        "event_evaluability": evaluability,
        "data_evaluability_gates": data_gates,
        "gates": gates,
        "confirmatory_pass": confirmatory_pass,
        "formal_outcome": formal_outcome,
        "claim_boundary": {
            "automatic_approval": False,
            "calibrated_probability": False,
            "interpretation": "review_priority_only",
            "no_same_species_v04_retry_after_reveal": True,
        },
    }

    output = Path(args.output_dir)
    partial = Path(str(output) + ".partial")
    if output.exists() or partial.exists():
        raise FileExistsError("Refusing to overwrite external v0.4 evaluation")
    partial.mkdir(parents=True)
    with (partial / "evaluation.json").open("x", encoding="utf-8", newline="") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
        handle.write("\n")
    with (partial / "candidates.tsv").open("x", encoding="utf-8", newline="") as handle:
        fields = (
            "candidate_digest",
            "seqid",
            "label_exact_cds",
            "conflict_set_digest",
            "v03_baseline_logit",
            "v03_primary_rank_score",
            "v04_primary_rank_score",
            "v04_topology_abstained",
        )
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for index, digest in enumerate(digests):
            writer.writerow(
                {
                    "candidate_digest": digest,
                    "seqid": groups[index],
                    "label_exact_cds": int(labels[index]),
                    "conflict_set_digest": conflicts[index],
                    "v03_baseline_logit": format(baseline[index], ".17g"),
                    "v03_primary_rank_score": format(v03[index], ".17g"),
                    "v04_primary_rank_score": format(v04[index], ".17g"),
                    "v04_topology_abstained": int(topology_abstained[index]),
                }
            )
    with (partial / "bootstrap_deltas.tsv").open("x", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("endpoint", "replicate", "delta"))
        for endpoint, values in (
            ("H1", h1_deltas),
            ("H2", h2_deltas),
            ("guard_v03_descriptive", guard_v03_deltas),
        ):
            for index, value in enumerate(values, start=1):
                writer.writerow((endpoint, index, format(float(value), ".17g")))
    checksum_path = partial / "SHA256SUMS"
    with checksum_path.open("x", encoding="utf-8", newline="") as handle:
        for path in sorted(partial.iterdir()):
            if path != checksum_path:
                handle.write(f"{sha256(path)}  {path.name}\n")
    os.replace(partial, output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
