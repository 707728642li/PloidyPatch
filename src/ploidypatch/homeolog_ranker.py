from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any, Mapping, Sequence

from . import __version__
from .baseline import _file_sha256
from .copy_features import COPY_FEATURE_SCHEMA_VERSION
from .copy_model import _read_feature_rows, transform_copy_feature_row
from .homeolog_topology import HOMEOLOG_TOPOLOGY_SCHEMA_VERSION


HOMEOLOG_RANKER_SCHEMA_VERSION = "ploidypatch.homeolog_copy_ranker.v2"
HOMEOLOG_RANK_SCORE_SCHEMA_VERSION = "ploidypatch.homeolog_copy_rank_scores.v2"
HOMEOLOG_REVIEW_RANKING_SCHEMA_VERSION = (
    "ploidypatch.homeolog_review_rankings.v1"
)
TOPOLOGY_COMPONENT_FIELDS = (
    "cds_bp_ratio",
    "cds_segment_count_ratio",
    "phase_lcs_similarity",
    "junction_fraction_similarity",
    "coding_span_ratio",
)
TOPOLOGY_ADDON_FIELDS = ("topology_available", *TOPOLOGY_COMPONENT_FIELDS)


def _load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"Expected a JSON object: {path}")
    return value


def _read_topology_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"candidate_digest", *TOPOLOGY_ADDON_FIELDS}
        if reader.fieldnames is None or not required <= set(reader.fieldnames):
            raise ValueError("Homeolog-topology table lacks frozen input columns")
        fieldnames = list(reader.fieldnames)
        rows = list(reader)
    forbidden = [
        field
        for field in fieldnames
        if field.startswith("label_") or field.startswith("truth_")
    ]
    if forbidden:
        raise ValueError(
            "Refusing evaluator-only columns in truth-blind topology scoring: "
            + ", ".join(forbidden)
        )
    digests = [row["candidate_digest"] for row in rows]
    if any(not digest for digest in digests) or len(digests) != len(set(digests)):
        raise ValueError("Topology candidate digests must be nonempty and unique")
    return fieldnames, rows


def topology_addons(row: Mapping[str, str]) -> list[float]:
    try:
        available = int(row.get("topology_available", ""))
    except ValueError as exc:
        raise ValueError("topology_available must be zero or one") from exc
    if available not in {0, 1}:
        raise ValueError("topology_available must be zero or one")
    values = [float(available)]
    for field in TOPOLOGY_COMPONENT_FIELDS:
        raw = row.get(field, "")
        if available:
            try:
                value = float(raw)
            except ValueError as exc:
                raise ValueError(
                    f"Available topology row lacks numeric {field}"
                ) from exc
            if not math.isfinite(value) or not 0 <= value <= 1:
                raise ValueError(f"Topology component outside [0, 1]: {field}")
            values.append(value)
        else:
            if raw not in {"", "0", "0.0"}:
                raise ValueError(
                    f"Unavailable topology row unexpectedly contains {field}"
                )
            values.append(0.0)
    return values


def _validate_estimator(
    estimator: Mapping[str, Any], *, expected_names: Sequence[str]
) -> None:
    names = estimator.get("coefficient_feature_order")
    coefficients = estimator.get("coefficients")
    if names != list(expected_names):
        raise ValueError("Ranker coefficient feature order disagrees with contract")
    if not isinstance(coefficients, list) or len(coefficients) != len(expected_names):
        raise ValueError("Ranker coefficient count disagrees with contract")
    numeric = [estimator.get("intercept"), *coefficients]
    if any(
        not isinstance(value, (int, float)) or not math.isfinite(float(value))
        for value in numeric
    ):
        raise ValueError("Ranker estimator contains non-finite parameters")


def validate_homeolog_ranker(model: Mapping[str, Any]) -> None:
    if model.get("schema_version") != HOMEOLOG_RANKER_SCHEMA_VERSION:
        raise ValueError("Unsupported homeolog-ranker schema version")
    contract = model.get("base_feature_contract")
    estimators = model.get("estimators")
    if not isinstance(contract, Mapping) or not isinstance(estimators, Mapping):
        raise ValueError("Homeolog ranker lacks contract or estimators")
    base_names = contract.get("expanded_feature_names")
    if not isinstance(base_names, list) or not base_names:
        raise ValueError("Homeolog ranker lacks expanded base feature names")
    if model.get("topology_addon_fields") != list(TOPOLOGY_ADDON_FIELDS):
        raise ValueError("Homeolog ranker topology feature contract changed")
    baseline = estimators.get("baseline")
    topology = estimators.get("topology")
    if not isinstance(baseline, Mapping) or not isinstance(topology, Mapping):
        raise ValueError("Homeolog ranker requires baseline and topology estimators")
    _validate_estimator(baseline, expected_names=base_names)
    _validate_estimator(
        topology, expected_names=[*base_names, *TOPOLOGY_ADDON_FIELDS]
    )
    boundary = model.get("claim_boundary")
    if not isinstance(boundary, Mapping) or boundary.get("automatic_approval") is not False:
        raise ValueError("Homeolog ranker must disable automatic approval")


def _sigmoid(value: float) -> float:
    if value >= 0:
        exponential = math.exp(-value)
        return 1.0 / (1.0 + exponential)
    exponential = math.exp(value)
    return exponential / (1.0 + exponential)


def _score(vector: Sequence[float], estimator: Mapping[str, Any]) -> tuple[float, float]:
    coefficients = estimator["coefficients"]
    if len(vector) != len(coefficients):
        raise ValueError("Scoring vector length disagrees with estimator")
    logit = float(estimator["intercept"]) + sum(
        value * float(coefficient)
        for value, coefficient in zip(vector, coefficients, strict=True)
    )
    return logit, _sigmoid(logit)


def _percentile_ranks(values: Sequence[float]) -> list[float]:
    if not values:
        return []
    if len(values) == 1:
        return [0.5]
    ordered = sorted(range(len(values)), key=lambda index: (values[index], index))
    ranks = [0.0] * len(values)
    cursor = 0
    while cursor < len(ordered):
        end = cursor + 1
        while end < len(ordered) and values[ordered[end]] == values[ordered[cursor]]:
            end += 1
        average_zero_based_rank = (cursor + end - 1) / 2
        percentile = average_zero_based_rank / (len(values) - 1)
        for index in ordered[cursor:end]:
            ranks[index] = percentile
        cursor = end
    return ranks


def score_homeolog_copy_candidates(
    *,
    copy_feature_tsv_path: str | Path,
    topology_tsv_path: str | Path,
    model_json_path: str | Path,
    output_tsv_path: str | Path,
) -> dict[str, Any]:
    """Score candidates with a frozen, truth-blind plant homeolog ranker."""

    copy_path = Path(copy_feature_tsv_path)
    topology_path = Path(topology_tsv_path)
    model_path = Path(model_json_path)
    output_path = Path(output_tsv_path)
    copy_manifest_path = Path(str(copy_path) + ".manifest.json")
    topology_manifest_path = Path(str(topology_path) + ".manifest.json")
    output_manifest_path = Path(str(output_path) + ".manifest.json")
    required = (
        copy_path,
        topology_path,
        model_path,
        copy_manifest_path,
        topology_manifest_path,
    )
    for path in required:
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f"Missing or empty homeolog-ranker input: {path}")
    collisions = [path for path in (output_path, output_manifest_path) if path.exists()]
    if collisions:
        raise FileExistsError(
            "Refusing to overwrite homeolog-ranker artifact(s): "
            + ", ".join(str(path) for path in collisions)
        )

    model = _load_json(model_path)
    validate_homeolog_ranker(model)
    copy_manifest = _load_json(copy_manifest_path)
    topology_manifest = _load_json(topology_manifest_path)
    if (
        copy_manifest.get("schema_version") != COPY_FEATURE_SCHEMA_VERSION
        or copy_manifest.get("truth_access") is not False
        or copy_manifest.get("outputs", {}).get("features", {}).get("sha256")
        != _file_sha256(copy_path)
    ):
        raise ValueError("Copy-feature manifest fails schema, truth, or checksum gate")
    if (
        topology_manifest.get("schema_version") != HOMEOLOG_TOPOLOGY_SCHEMA_VERSION
        or topology_manifest.get("truth_access") is not False
        or topology_manifest.get("outputs", {}).get("features", {}).get("sha256")
        != _file_sha256(topology_path)
        or topology_manifest.get("inputs", {}).get("copy_features")
        != _file_sha256(copy_path)
    ):
        raise ValueError("Topology manifest fails schema, truth, lineage, or checksum gate")

    copy_fields, copy_rows = _read_feature_rows(copy_path)
    _, topology_rows = _read_topology_rows(topology_path)
    topology_by_digest = {row["candidate_digest"]: row for row in topology_rows}
    copy_digests = [row["candidate_digest"] for row in copy_rows]
    if set(copy_digests) != set(topology_by_digest):
        raise ValueError("Copy and topology candidate universes differ")

    contract = model["base_feature_contract"]
    baseline_estimator = model["estimators"]["baseline"]
    topology_estimator = model["estimators"]["topology"]
    scored: list[dict[str, Any]] = []
    primary_scores: list[float] = []
    available_count = 0
    for copy_row in copy_rows:
        topology_row = topology_by_digest[copy_row["candidate_digest"]]
        base_vector = transform_copy_feature_row(copy_row, contract)
        addons = topology_addons(topology_row)
        baseline_logit, baseline_score = _score(base_vector, baseline_estimator)
        primary_logit, primary_score = _score(
            [*base_vector, *addons], topology_estimator
        )
        available = int(addons[0])
        available_count += available
        primary_scores.append(primary_score)
        scored.append(
            {
                **copy_row,
                "homeolog_baseline_logit": format(baseline_logit, ".17g"),
                "homeolog_baseline_rank_score": format(baseline_score, ".17g"),
                "homeolog_topology_logit": format(primary_logit, ".17g"),
                "homeolog_topology_rank_score": format(primary_score, ".17g"),
                "homeolog_topology_available": available,
                "homeolog_evidence_tier": (
                    "existing_wgd_partner_topology"
                    if available
                    else "generic_candidate_topology_unavailable"
                ),
                "homeolog_automatic_approval": 0,
            }
        )
    percentiles = _percentile_ranks(primary_scores)
    output_fields = [
        *copy_fields,
        "homeolog_baseline_logit",
        "homeolog_baseline_rank_score",
        "homeolog_topology_logit",
        "homeolog_topology_rank_score",
        "homeolog_topology_rank_percentile",
        "homeolog_topology_available",
        "homeolog_evidence_tier",
        "homeolog_automatic_approval",
    ]
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=output_fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for row, percentile in zip(scored, percentiles, strict=True):
            row["homeolog_topology_rank_percentile"] = format(percentile, ".17g")
            writer.writerow(row)

    manifest: dict[str, Any] = {
        "schema_version": HOMEOLOG_RANK_SCORE_SCHEMA_VERSION,
        "generator": {"name": "PloidyPatch", "version": __version__},
        "truth_access": False,
        "inputs": {
            "copy_features": _file_sha256(copy_path),
            "topology_features": _file_sha256(topology_path),
            "model": _file_sha256(model_path),
        },
        "policy": {
            "primary_score": "homeolog_topology_rank_score",
            "score_interpretation": "uncalibrated_monotonic_review_rank_not_probability",
            "within_run_percentile": True,
            "automatic_approval": False,
            "topology_unavailable_policy": "retain_generic_rank_and_flag_unavailable",
        },
        "counts": {
            "candidates": len(scored),
            "topology_available": available_count,
            "topology_unavailable": len(scored) - available_count,
            "automatic_approved": 0,
        },
        "outputs": {
            "scores": {
                "file_name": output_path.name,
                "sha256": _file_sha256(output_path),
                "rows": len(scored),
            }
        },
    }
    with output_manifest_path.open("x", encoding="utf-8", newline="") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return manifest


def freeze_homeolog_review_rankings(
    *,
    score_tsv_path: str | Path,
    output_tsv_path: str | Path,
    review_budgets: Sequence[int] = (25, 50, 100, 200),
) -> dict[str, Any]:
    """Freeze deterministic baseline and topology review queues.

    The output deliberately contains ranks rather than decisions. It is used
    to establish top-K membership before orthogonal natural-validation
    evidence is opened.
    """

    score_path = Path(score_tsv_path)
    score_manifest_path = Path(str(score_path) + ".manifest.json")
    output_path = Path(output_tsv_path)
    output_manifest_path = Path(str(output_path) + ".manifest.json")
    for path in (score_path, score_manifest_path):
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f"Missing or empty review-ranking input: {path}")
    collisions = [path for path in (output_path, output_manifest_path) if path.exists()]
    if collisions:
        raise FileExistsError(
            "Refusing to overwrite review-ranking artifact(s): "
            + ", ".join(str(path) for path in collisions)
        )

    budgets = sorted(set(review_budgets))
    if not budgets or any(
        not isinstance(budget, int) or isinstance(budget, bool) or budget <= 0
        for budget in budgets
    ):
        raise ValueError("Review budgets must be positive integers")

    score_manifest = _load_json(score_manifest_path)
    if (
        score_manifest.get("schema_version") != HOMEOLOG_RANK_SCORE_SCHEMA_VERSION
        or score_manifest.get("truth_access") is not False
        or score_manifest.get("outputs", {}).get("scores", {}).get("sha256")
        != _file_sha256(score_path)
    ):
        raise ValueError("Rank-score manifest fails schema, truth, or checksum gate")

    with score_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "candidate_digest",
            "seqid",
            "start",
            "end",
            "strand",
            "homeolog_baseline_rank_score",
            "homeolog_topology_rank_score",
            "homeolog_topology_available",
        }
        if reader.fieldnames is None or not required <= set(reader.fieldnames):
            raise ValueError("Rank-score table lacks review-ranking columns")
        forbidden = [
            field
            for field in reader.fieldnames
            if field.startswith("label_") or field.startswith("truth_")
        ]
        if forbidden:
            raise ValueError(
                "Refusing evaluator-only columns in natural review freeze: "
                + ", ".join(forbidden)
            )
        rows = list(reader)
    digests = [row["candidate_digest"] for row in rows]
    if any(not digest for digest in digests) or len(digests) != len(set(digests)):
        raise ValueError("Review-ranking candidate digests must be nonempty and unique")

    estimators = {
        "baseline": "homeolog_baseline_rank_score",
        "topology": "homeolog_topology_rank_score",
    }
    output_fields = [
        "estimator",
        "review_rank",
        "candidate_digest",
        "seqid",
        "start",
        "end",
        "strand",
        "rank_score",
        "topology_available",
        "support_method_count",
        "support_methods",
        *[f"within_top_{budget}" for budget in budgets],
        "automatic_approval",
    ]
    ranked_rows: list[dict[str, Any]] = []
    topologies_by_estimator: dict[str, dict[str, int]] = {}
    for estimator, score_field in estimators.items():
        parsed: list[tuple[float, dict[str, str]]] = []
        for row in rows:
            try:
                score = float(row[score_field])
            except ValueError as exc:
                raise ValueError(f"Non-numeric review score: {score_field}") from exc
            if not math.isfinite(score) or not 0 <= score <= 1:
                raise ValueError(f"Review score outside [0, 1]: {score_field}")
            parsed.append((score, row))
        ordered = sorted(
            parsed,
            key=lambda item: (-item[0], item[1]["candidate_digest"]),
        )
        topologies_by_estimator[estimator] = {}
        for rank, (score, row) in enumerate(ordered, start=1):
            available = row["homeolog_topology_available"]
            if available not in {"0", "1"}:
                raise ValueError("homeolog_topology_available must be zero or one")
            ranked_rows.append(
                {
                    "estimator": estimator,
                    "review_rank": rank,
                    "candidate_digest": row["candidate_digest"],
                    "seqid": row["seqid"],
                    "start": row["start"],
                    "end": row["end"],
                    "strand": row["strand"],
                    "rank_score": format(score, ".17g"),
                    "topology_available": available,
                    "support_method_count": row.get("support_method_count", ""),
                    "support_methods": row.get("support_methods", ""),
                    **{
                        f"within_top_{budget}": int(rank <= budget)
                        for budget in budgets
                    },
                    "automatic_approval": 0,
                }
            )
        for budget in budgets:
            topologies_by_estimator[estimator][str(budget)] = sum(
                int(row["homeolog_topology_available"])
                for _, row in ordered[:budget]
            )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=output_fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(ranked_rows)

    manifest: dict[str, Any] = {
        "schema_version": HOMEOLOG_REVIEW_RANKING_SCHEMA_VERSION,
        "generator": {"name": "PloidyPatch", "version": __version__},
        "truth_access": False,
        "inputs": {"rank_scores": _file_sha256(score_path)},
        "policy": {
            "estimators": list(estimators),
            "review_budgets": budgets,
            "order": "descending_score_then_candidate_digest",
            "score_interpretation": "uncalibrated_review_rank_not_probability",
            "automatic_approval": False,
        },
        "counts": {
            "candidates": len(rows),
            "ranking_rows": len(ranked_rows),
            "topology_available_within_budget": topologies_by_estimator,
        },
        "outputs": {
            "rankings": {
                "file_name": output_path.name,
                "sha256": _file_sha256(output_path),
                "rows": len(ranked_rows),
            }
        },
    }
    with output_manifest_path.open("x", encoding="utf-8", newline="") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return manifest
