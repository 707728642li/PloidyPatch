#!/usr/bin/env python3
"""Decompose the preregistered Actinidia topology-coverage gate.

This is a descriptive, post-reveal audit.  It verifies the frozen evaluation
and blind-score lineages, reproduces the published coverage numerator, and
partitions positive candidates into topology-availability/guard-abstention
states.  It never changes scores, labels, thresholds, or the formal outcome.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return value


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        return fields, list(reader)


def read_policy(path: Path) -> dict[str, str]:
    fields, rows = read_tsv(path)
    if fields != ["field", "value"]:
        raise ValueError("Actinidia policy must have exact field/value header")
    values: dict[str, str] = {}
    for row in rows:
        key = row["field"]
        if not key or key in values:
            raise ValueError(f"Duplicate or empty policy key: {key!r}")
        values[key] = row["value"]
    return values


def require_binary(row: dict[str, str], field: str) -> int:
    value = row.get(field, "")
    if value not in {"0", "1"}:
        raise ValueError(f"Non-binary {field}: {value!r}")
    return int(value)


def state_counts(rows: Iterable[dict[str, str]]) -> dict[str, int]:
    counts = {
        "positive_candidates": 0,
        "topology_available_positive": 0,
        "topology_unavailable_positive": 0,
        "topology_abstained_positive": 0,
        "effective_topology_positive": 0,
        "available_not_abstained_positive": 0,
        "available_abstained_positive": 0,
        "unavailable_not_abstained_positive": 0,
        "unavailable_abstained_positive": 0,
    }
    for row in rows:
        if require_binary(row, "label_exact_cds") != 1:
            continue
        available = require_binary(row, "v03_topology_available")
        abstained = require_binary(row, "v04_topology_abstained")
        counts["positive_candidates"] += 1
        counts[
            "topology_available_positive"
            if available
            else "topology_unavailable_positive"
        ] += 1
        if abstained:
            counts["topology_abstained_positive"] += 1
        if available and not abstained:
            counts["effective_topology_positive"] += 1
            counts["available_not_abstained_positive"] += 1
        elif available and abstained:
            counts["available_abstained_positive"] += 1
        elif not available and abstained:
            counts["unavailable_abstained_positive"] += 1
        else:
            counts["unavailable_not_abstained_positive"] += 1
    return counts


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, Any]]) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evaluation", type=Path, required=True)
    parser.add_argument("--candidates", type=Path, required=True)
    parser.add_argument("--evaluation-sums", type=Path, required=True)
    parser.add_argument("--scores", type=Path, required=True)
    parser.add_argument("--score-manifest", type=Path, required=True)
    parser.add_argument("--policy", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    evaluation = load_json(args.evaluation)
    score_manifest = load_json(args.score_manifest)
    policy = read_policy(args.policy)
    threshold_key = "minimum_topology_coverage_among_positive_candidates"
    threshold = float(policy[threshold_key])
    if not (0.0 <= threshold <= 1.0):
        raise ValueError("Topology coverage threshold is outside [0, 1]")

    expected_evaluation_files: dict[str, str] = {}
    with args.evaluation_sums.open(encoding="utf-8") as handle:
        for raw in handle:
            digest, relative = raw.rstrip("\n").split(maxsplit=1)
            expected_evaluation_files[relative.lstrip("*")] = digest
    for name, path in {
        "evaluation.json": args.evaluation,
        "candidates.tsv": args.candidates,
    }.items():
        if expected_evaluation_files.get(name) != sha256(path):
            raise ValueError(f"Evaluation checksum mismatch: {name}")

    score_sha = sha256(args.scores)
    if (
        evaluation["inputs"]["scores"]["sha256"] != score_sha
        or score_manifest["outputs"]["scores"]["sha256"] != score_sha
    ):
        raise ValueError("Blind scores disagree with evaluation or score manifest")
    if evaluation["inputs"]["policy"]["sha256"] != sha256(args.policy):
        raise ValueError("Policy differs from the frozen evaluation input")

    candidate_fields, candidates = read_tsv(args.candidates)
    score_fields, scores = read_tsv(args.scores)
    required_candidates = {"candidate_digest", "seqid", "label_exact_cds"}
    required_scores = {
        "candidate_digest",
        "seqid",
        "support_method_count",
        "support_methods",
        "wgd_existing_partner",
        "wgd_support_block_count",
        "v03_topology_available",
        "v04_conflict_guard_applied",
        "v04_topology_abstained",
    }
    if not required_candidates <= set(candidate_fields):
        raise ValueError("Candidate evaluation table lacks required fields")
    if not required_scores <= set(score_fields):
        raise ValueError("Blind score table lacks required fields")

    candidates_by_digest = {row["candidate_digest"]: row for row in candidates}
    scores_by_digest = {row["candidate_digest"]: row for row in scores}
    if (
        len(candidates_by_digest) != len(candidates)
        or len(scores_by_digest) != len(scores)
        or set(candidates_by_digest) != set(scores_by_digest)
    ):
        raise ValueError("Candidate and score universes differ or contain duplicates")

    joined: list[dict[str, str]] = []
    for digest in sorted(candidates_by_digest):
        candidate = candidates_by_digest[digest]
        score = scores_by_digest[digest]
        if candidate["seqid"] != score["seqid"]:
            raise ValueError(f"Seqid mismatch for {digest}")
        if require_binary(score, "v04_conflict_guard_applied") != require_binary(
            score, "v04_topology_abstained"
        ):
            raise ValueError("v0.4 abstention must equal conflict-guard application")
        joined.append({**score, "label_exact_cds": candidate["label_exact_cds"]})

    overall = state_counts(joined)
    total_candidates = len(joined)
    positive_total = overall["positive_candidates"]
    effective = overall["effective_topology_positive"]
    coverage = effective / positive_total
    minimum_required = math.ceil(threshold * positive_total)
    if total_candidates != int(evaluation["counts"]["candidates"]):
        raise ValueError("Candidate count does not reproduce formal evaluation")
    if positive_total != int(evaluation["counts"]["positive_exact_cds_candidates"]):
        raise ValueError("Positive count does not reproduce formal evaluation")
    if effective != int(evaluation["counts"]["effective_topology_positive_candidates"]):
        raise ValueError("Effective-topology numerator does not reproduce evaluation")
    if not math.isclose(
        coverage,
        float(evaluation["counts"]["topology_positive_coverage"]),
        rel_tol=0.0,
        abs_tol=1e-15,
    ):
        raise ValueError("Topology coverage does not reproduce formal evaluation")

    chromosome_groups: dict[str, list[dict[str, str]]] = defaultdict(list)
    support_groups: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in joined:
        chromosome_groups[row["seqid"]].append(row)
        support_groups[(row["support_method_count"], row["support_methods"])].append(row)

    def group_rows(groups: dict[Any, list[dict[str, str]]], keys: tuple[str, ...]):
        output: list[dict[str, Any]] = []
        for group, rows in sorted(groups.items(), key=lambda item: item[0]):
            values = group if isinstance(group, tuple) else (group,)
            counts = state_counts(rows)
            positives = counts["positive_candidates"]
            output.append(
                {
                    **dict(zip(keys, values, strict=True)),
                    "candidate_count": len(rows),
                    **counts,
                    "topology_positive_coverage": (
                        counts["effective_topology_positive"] / positives
                        if positives
                        else ""
                    ),
                }
            )
        return output

    chromosome_rows = group_rows(chromosome_groups, ("seqid",))
    support_rows = group_rows(
        support_groups, ("support_method_count", "support_methods")
    )

    output = args.output_dir
    working = Path(str(output) + ".working")
    if output.exists() or working.exists():
        raise FileExistsError("Refusing to overwrite topology audit")
    working.mkdir(parents=True)
    common_fields = [
        "candidate_count",
        *overall.keys(),
        "topology_positive_coverage",
    ]
    write_tsv(
        working / "by_chromosome.tsv",
        ["seqid", *common_fields],
        chromosome_rows,
    )
    write_tsv(
        working / "by_method_support.tsv",
        ["support_method_count", "support_methods", *common_fields],
        support_rows,
    )
    summary = {
        "schema_version": "ploidypatch.actinidia_topology_gate_audit.v0.5",
        "audit_role": "descriptive_post_reveal_no_model_or_threshold_change",
        "formal_outcome_unchanged": evaluation["formal_outcome"],
        "formal_confirmatory_pass_unchanged": evaluation["confirmatory_pass"],
        "threshold": threshold,
        "threshold_source_key": threshold_key,
        "coverage": coverage,
        "minimum_effective_positive_candidates_required": minimum_required,
        "effective_positive_candidate_shortfall": max(0, minimum_required - effective),
        "counts": {"candidates": total_candidates, **overall},
        "gate_reproduced": coverage >= threshold,
        "sources": {
            "evaluation_json_sha256": sha256(args.evaluation),
            "candidates_tsv_sha256": sha256(args.candidates),
            "blind_scores_sha256": score_sha,
            "blind_score_manifest_sha256": sha256(args.score_manifest),
            "policy_sha256": sha256(args.policy),
        },
    }
    with (working / "summary.json").open("x", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")
    output_files = sorted(path for path in working.iterdir() if path.is_file())
    with (working / "SHA256SUMS").open("x", encoding="utf-8", newline="") as handle:
        for path in output_files:
            handle.write(f"{sha256(path)}  {path.name}\n")
    os.replace(working, output)


if __name__ == "__main__":
    main()
