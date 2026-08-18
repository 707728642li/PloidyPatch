#!/usr/bin/env python3
"""Post-hoc, no-tuning diagnosis of the frozen maize v0.2 ranker result."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import Counter
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np
from sklearn.metrics import average_precision_score, roc_auc_score

from ploidypatch.copy_model import transform_copy_feature_row
from ploidypatch.homeolog_ranker import TOPOLOGY_ADDON_FIELDS, topology_addons


SCHEMA_VERSION = "ploidypatch.maize_homeolog_ranker_diagnostic.v1"
TOPOLOGY_COMPONENT_FIELDS = tuple(TOPOLOGY_ADDON_FIELDS[1:])


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"TSV has no header: {path}")
        return list(reader)


def keyed(rows: Iterable[dict[str, str]]) -> dict[str, dict[str, str]]:
    result: dict[str, dict[str, str]] = {}
    for row in rows:
        digest = row.get("candidate_digest", "")
        if not digest or digest in result:
            raise ValueError("Candidate digests must be nonempty and unique")
        result[digest] = row
    return result


def finite_float(row: Mapping[str, str], field: str) -> float:
    try:
        value = float(row[field])
    except (KeyError, ValueError) as exc:
        raise ValueError(f"Missing or nonnumeric {field}") from exc
    if not math.isfinite(value):
        raise ValueError(f"Nonfinite {field}")
    return value


def binary_metric(labels: Sequence[int], scores: Sequence[float]) -> dict[str, Any]:
    positives = int(sum(labels))
    result: dict[str, Any] = {
        "candidates": len(labels),
        "positives": positives,
        "positive_fraction": positives / len(labels) if labels else None,
    }
    if labels and 0 < positives < len(labels):
        result["average_precision"] = float(average_precision_score(labels, scores))
        result["roc_auc"] = float(roc_auc_score(labels, scores))
    else:
        result["average_precision"] = None
        result["roc_auc"] = None
    return result


def numeric_summary(values: Sequence[float]) -> dict[str, float | int | None]:
    if not values:
        return {"n": 0, "mean": None, "median": None, "q25": None, "q75": None}
    array = np.asarray(values, dtype=float)
    return {
        "n": len(values),
        "mean": float(np.mean(array)),
        "median": float(np.median(array)),
        "q25": float(np.quantile(array, 0.25)),
        "q75": float(np.quantile(array, 0.75)),
    }


def deterministic_ranks(
    rows: Sequence[dict[str, Any]], score_field: str
) -> dict[str, int]:
    ordered = sorted(
        rows,
        key=lambda row: (-float(row[score_field]), str(row["candidate_digest"])),
    )
    return {
        str(row["candidate_digest"]): rank
        for rank, row in enumerate(ordered, start=1)
    }


def top_k_summary(rows: Sequence[dict[str, Any]], k: int) -> dict[str, Any]:
    ordered: dict[str, list[dict[str, Any]]] = {}
    for estimator in ("baseline", "topology"):
        score_field = f"{estimator}_score"
        ordered[estimator] = sorted(
            rows,
            key=lambda row: (
                -float(row[score_field]),
                str(row["candidate_digest"]),
            ),
        )[:k]
    baseline_ids = {str(row["candidate_digest"]) for row in ordered["baseline"]}
    topology_ids = {str(row["candidate_digest"]) for row in ordered["topology"]}
    result: dict[str, Any] = {
        "requested_k": k,
        "actual_k": min(k, len(rows)),
        "set_overlap": len(baseline_ids & topology_ids),
        "baseline_only": len(baseline_ids - topology_ids),
        "topology_only": len(topology_ids - baseline_ids),
    }
    for estimator in ("baseline", "topology"):
        subset = ordered[estimator]
        positives = sum(int(row["label"]) for row in subset)
        result[estimator] = {
            "positives": positives,
            "precision": positives / len(subset) if subset else None,
            "topology_available": sum(int(row["topology_available"]) for row in subset),
        }
    result["topology_minus_baseline_positives"] = (
        result["topology"]["positives"] - result["baseline"]["positives"]
    )
    return result


def grouped_component_summary(
    rows: Sequence[dict[str, Any]], field: str
) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for label in (0, 1):
        values = [
            float(row[field])
            for row in rows
            if int(row["label"]) == label and int(row["topology_available"]) == 1
        ]
        result["positive" if label else "negative"] = numeric_summary(values)
    available = [row for row in rows if int(row["topology_available"]) == 1]
    labels = [int(row["label"]) for row in available]
    values = [float(row[field]) for row in available]
    result["available_subset_ranking"] = binary_metric(labels, values)
    return result


def diagnose(
    *,
    scores_path: Path,
    labels_path: Path,
    topology_path: Path,
    model_path: Path,
    output_path: Path,
) -> dict[str, Any]:
    if output_path.exists():
        raise FileExistsError(f"Refusing to overwrite diagnostic: {output_path}")
    score_rows = keyed(read_tsv(scores_path))
    label_rows = keyed(read_tsv(labels_path))
    topology_rows = keyed(read_tsv(topology_path))
    if set(score_rows) != set(label_rows) or set(score_rows) != set(topology_rows):
        raise ValueError("Score, label, and topology candidate universes differ")
    model = json.loads(model_path.read_text(encoding="utf-8"))
    if not isinstance(model, dict):
        raise ValueError("Model must be a JSON object")

    contract = model["base_feature_contract"]
    topology_estimator = model["estimators"]["topology"]
    base_count = len(contract["expanded_feature_names"])
    topology_coefficients = [float(value) for value in topology_estimator["coefficients"]]
    topology_base_coefficients = topology_coefficients[:base_count]
    addon_coefficients = topology_coefficients[base_count:]
    if len(addon_coefficients) != len(TOPOLOGY_ADDON_FIELDS):
        raise ValueError("Topology add-on coefficient contract changed")

    joined: list[dict[str, Any]] = []
    maximum_logit_error = 0.0
    for digest in sorted(score_rows):
        score = score_rows[digest]
        label = label_rows[digest]
        topology = topology_rows[digest]
        base_vector = transform_copy_feature_row(score, contract)
        addons = topology_addons(topology)
        topology_base_logit = float(topology_estimator["intercept"]) + sum(
            value * coefficient
            for value, coefficient in zip(
                base_vector, topology_base_coefficients, strict=True
            )
        )
        addon_logit = sum(
            value * coefficient
            for value, coefficient in zip(addons, addon_coefficients, strict=True)
        )
        observed_baseline_logit = finite_float(score, "homeolog_baseline_logit")
        observed_topology_logit = finite_float(score, "homeolog_topology_logit")
        maximum_logit_error = max(
            maximum_logit_error,
            abs(topology_base_logit + addon_logit - observed_topology_logit),
        )
        row: dict[str, Any] = {
            "candidate_digest": digest,
            "seqid": score["seqid"],
            "label": int(label["label_exact_cds"]),
            "support_methods": score.get("support_methods", ""),
            "support_method_count": int(score.get("support_method_count", "0")),
            "baseline_score": finite_float(score, "homeolog_baseline_rank_score"),
            "topology_score": finite_float(score, "homeolog_topology_rank_score"),
            "topology_available": int(topology["topology_available"]),
            "base_refit_logit_delta": topology_base_logit - observed_baseline_logit,
            "addon_logit_contribution": addon_logit,
            "total_logit_delta": observed_topology_logit - observed_baseline_logit,
        }
        for field in TOPOLOGY_COMPONENT_FIELDS:
            row[field] = float(topology[field]) if int(topology["topology_available"]) else 0.0
        joined.append(row)
    if maximum_logit_error > 1e-10:
        raise ValueError("Frozen topology score cannot be reconstructed from model")

    baseline_ranks = deterministic_ranks(joined, "baseline_score")
    topology_ranks = deterministic_ranks(joined, "topology_score")
    for row in joined:
        digest = str(row["candidate_digest"])
        row["rank_change_positive_is_improvement"] = (
            baseline_ranks[digest] - topology_ranks[digest]
        )

    result: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "analysis_role": "post_hoc_failure_diagnosis_not_model_selection",
        "policy": {
            "maize_refit": False,
            "maize_threshold_tuning": False,
            "changes_formal_gate": False,
            "automatic_approval": False,
        },
        "inputs": {
            "scores": {"path": str(scores_path), "sha256": sha256(scores_path)},
            "labels": {"path": str(labels_path), "sha256": sha256(labels_path)},
            "topology": {"path": str(topology_path), "sha256": sha256(topology_path)},
            "model": {"path": str(model_path), "sha256": sha256(model_path)},
        },
        "counts": {
            "candidates": len(joined),
            "positives": sum(int(row["label"]) for row in joined),
            "topology_available": sum(int(row["topology_available"]) for row in joined),
            "maximum_reconstructed_logit_error": maximum_logit_error,
        },
        "model_topology_addon_coefficients": dict(
            zip(TOPOLOGY_ADDON_FIELDS, addon_coefficients, strict=True)
        ),
    }

    labels = [int(row["label"]) for row in joined]
    result["overall_metrics"] = {
        "baseline": binary_metric(labels, [float(row["baseline_score"]) for row in joined]),
        "topology": binary_metric(labels, [float(row["topology_score"]) for row in joined]),
    }
    result["by_topology_availability"] = {}
    for available in (0, 1):
        subset = [row for row in joined if int(row["topology_available"]) == available]
        subset_labels = [int(row["label"]) for row in subset]
        result["by_topology_availability"][str(available)] = {
            "baseline": binary_metric(
                subset_labels, [float(row["baseline_score"]) for row in subset]
            ),
            "topology": binary_metric(
                subset_labels, [float(row["topology_score"]) for row in subset]
            ),
        }
    available_positive_fraction = result["by_topology_availability"]["1"]["baseline"][
        "positive_fraction"
    ]
    unavailable_positive_fraction = result["by_topology_availability"]["0"]["baseline"][
        "positive_fraction"
    ]
    result["topology_presence_enrichment"] = {
        "available_positive_fraction": available_positive_fraction,
        "unavailable_positive_fraction": unavailable_positive_fraction,
        "risk_ratio": (
            available_positive_fraction / unavailable_positive_fraction
            if unavailable_positive_fraction
            else None
        ),
    }

    result["review_budgets"] = {
        str(k): top_k_summary(joined, k) for k in (25, 50, 79, 100, 157, 200, 313)
    }
    result["topology_components_within_available"] = {
        field: grouped_component_summary(joined, field)
        for field in TOPOLOGY_COMPONENT_FIELDS
    }

    result["logit_decomposition"] = {}
    for group_name, predicate in {
        "all": lambda row: True,
        "positive": lambda row: int(row["label"]) == 1,
        "negative": lambda row: int(row["label"]) == 0,
        "available_positive": lambda row: int(row["label"]) == 1
        and int(row["topology_available"]) == 1,
        "available_negative": lambda row: int(row["label"]) == 0
        and int(row["topology_available"]) == 1,
        "unavailable_positive": lambda row: int(row["label"]) == 1
        and int(row["topology_available"]) == 0,
        "unavailable_negative": lambda row: int(row["label"]) == 0
        and int(row["topology_available"]) == 0,
    }.items():
        subset = [row for row in joined if predicate(row)]
        result["logit_decomposition"][group_name] = {
            field: numeric_summary([float(row[field]) for row in subset])
            for field in (
                "base_refit_logit_delta",
                "addon_logit_contribution",
                "total_logit_delta",
                "rank_change_positive_is_improvement",
            )
        }

    per_chromosome: dict[str, Any] = {}
    for seqid in sorted({str(row["seqid"]) for row in joined}, key=lambda value: int(value)):
        subset = [row for row in joined if str(row["seqid"]) == seqid]
        subset_labels = [int(row["label"]) for row in subset]
        baseline = binary_metric(
            subset_labels, [float(row["baseline_score"]) for row in subset]
        )
        topology = binary_metric(
            subset_labels, [float(row["topology_score"]) for row in subset]
        )
        per_chromosome[seqid] = {
            "baseline": baseline,
            "topology": topology,
            "delta_average_precision": (
                topology["average_precision"] - baseline["average_precision"]
                if topology["average_precision"] is not None
                and baseline["average_precision"] is not None
                else None
            ),
        }
    result["per_chromosome"] = per_chromosome

    by_method: dict[str, Any] = {}
    for support_methods in sorted({str(row["support_methods"]) for row in joined}):
        subset = [row for row in joined if str(row["support_methods"]) == support_methods]
        subset_labels = [int(row["label"]) for row in subset]
        baseline = binary_metric(
            subset_labels, [float(row["baseline_score"]) for row in subset]
        )
        topology = binary_metric(
            subset_labels, [float(row["topology_score"]) for row in subset]
        )
        by_method[support_methods] = {
            "baseline": baseline,
            "topology": topology,
            "delta_average_precision": (
                topology["average_precision"] - baseline["average_precision"]
                if topology["average_precision"] is not None
                and baseline["average_precision"] is not None
                else None
            ),
        }
    result["by_method_support_pattern"] = by_method
    result["topology_reason_counts"] = dict(
        sorted(Counter(row.get("topology_reason", "") for row in topology_rows.values()).items())
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("x", encoding="utf-8", newline="") as handle:
        json.dump(result, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return result


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scores", type=Path, required=True)
    parser.add_argument("--labels", type=Path, required=True)
    parser.add_argument("--topology", type=Path, required=True)
    parser.add_argument("--model", type=Path, required=True)
    parser.add_argument("--output-json", type=Path, required=True)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    report = diagnose(
        scores_path=args.scores,
        labels_path=args.labels,
        topology_path=args.topology,
        model_path=args.model,
        output_path=args.output_json,
    )
    print(json.dumps({"counts": report["counts"], "output": str(args.output_json)}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
