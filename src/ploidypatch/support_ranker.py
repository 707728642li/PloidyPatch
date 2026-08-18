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
from .homeolog_ranker import (
    TOPOLOGY_COMPONENT_FIELDS,
    _load_json,
    _percentile_ranks,
    _read_topology_rows,
    topology_addons,
)
from .homeolog_topology import HOMEOLOG_TOPOLOGY_SCHEMA_VERSION


SUPPORT_CONDITIONED_RANKER_SCHEMA_VERSION = (
    "ploidypatch.support_conditioned_homeolog_ranker.v1"
)
SUPPORT_CONDITIONED_SCORE_SCHEMA_VERSION = (
    "ploidypatch.support_conditioned_homeolog_scores.v1"
)
SUPPORT_PATTERNS = (
    "gemoma",
    "lifton",
    "miniprot",
    "gemoma,lifton",
    "gemoma,miniprot",
    "lifton,miniprot",
    "gemoma,lifton,miniprot",
)
CORRECTION_FEATURE_NAMES = (
    *(f"topology:{field}" for field in TOPOLOGY_COMPONENT_FIELDS),
    *(
        name
        for pattern in SUPPORT_PATTERNS
        for name in (
            f"topology_available:{pattern}",
            f"topology_mean_coherence:{pattern}",
        )
    ),
)


def support_conditioned_correction_vector(
    copy_row: Mapping[str, str], topology_row: Mapping[str, str]
) -> list[float]:
    """Build the frozen v0.3 topology-only offset design vector."""

    pattern = copy_row.get("support_methods", "")
    if pattern not in SUPPORT_PATTERNS:
        raise ValueError(f"Unsupported method support pattern: {pattern!r}")
    addons = topology_addons(topology_row)
    available = addons[0]
    components = addons[1:]
    mean_coherence = sum(components) / len(components)
    values = list(components)
    for level in SUPPORT_PATTERNS:
        indicator = float(pattern == level)
        values.extend((available * indicator, available * mean_coherence * indicator))
    if len(values) != len(CORRECTION_FEATURE_NAMES):
        raise AssertionError("Support-conditioned correction contract changed")
    if not available and any(values):
        raise AssertionError("Unavailable topology must produce exactly zero correction")
    return values


def _validate_estimator(
    estimator: Mapping[str, Any], expected_names: Sequence[str], *, allow_intercept: bool
) -> None:
    if estimator.get("coefficient_feature_order") != list(expected_names):
        raise ValueError("Estimator coefficient order disagrees with its contract")
    coefficients = estimator.get("coefficients")
    if not isinstance(coefficients, list) or len(coefficients) != len(expected_names):
        raise ValueError("Estimator coefficient count disagrees with its contract")
    numeric: list[Any] = list(coefficients)
    if allow_intercept:
        numeric.append(estimator.get("intercept"))
    elif estimator.get("intercept") not in {0, 0.0}:
        raise ValueError("Topology offset must not learn a second intercept")
    if any(
        not isinstance(value, (int, float)) or not math.isfinite(float(value))
        for value in numeric
    ):
        raise ValueError("Ranker contains non-finite numeric parameters")


def validate_support_conditioned_ranker(model: Mapping[str, Any]) -> None:
    if model.get("schema_version") != SUPPORT_CONDITIONED_RANKER_SCHEMA_VERSION:
        raise ValueError("Unsupported support-conditioned ranker schema")
    contract = model.get("base_feature_contract")
    estimators = model.get("estimators")
    if not isinstance(contract, Mapping) or not isinstance(estimators, Mapping):
        raise ValueError("Support-conditioned ranker lacks contract or estimators")
    base_names = contract.get("expanded_feature_names")
    if not isinstance(base_names, list) or not base_names:
        raise ValueError("Support-conditioned ranker lacks base feature names")
    baseline = estimators.get("baseline")
    correction = estimators.get("support_conditioned_topology_offset")
    if not isinstance(baseline, Mapping) or not isinstance(correction, Mapping):
        raise ValueError("Support-conditioned ranker lacks required estimators")
    _validate_estimator(baseline, base_names, allow_intercept=True)
    _validate_estimator(
        correction, CORRECTION_FEATURE_NAMES, allow_intercept=False
    )
    if model.get("support_patterns") != list(SUPPORT_PATTERNS):
        raise ValueError("Support-pattern contract changed")
    if model.get("topology_component_fields") != list(TOPOLOGY_COMPONENT_FIELDS):
        raise ValueError("Topology-component contract changed")
    boundary = model.get("claim_boundary")
    if (
        not isinstance(boundary, Mapping)
        or boundary.get("automatic_approval") is not False
        or boundary.get("calibrated_probability") is not False
    ):
        raise ValueError("v0.3 ranker must remain review-only and uncalibrated")


def _linear_score(vector: Sequence[float], estimator: Mapping[str, Any]) -> float:
    coefficients = estimator["coefficients"]
    if len(vector) != len(coefficients):
        raise ValueError("Scoring vector length disagrees with estimator")
    return float(estimator.get("intercept", 0.0)) + sum(
        float(value) * float(coefficient)
        for value, coefficient in zip(vector, coefficients, strict=True)
    )


def score_support_conditioned_candidates(
    *,
    copy_feature_tsv_path: str | Path,
    topology_tsv_path: str | Path,
    model_json_path: str | Path,
    output_tsv_path: str | Path,
) -> dict[str, Any]:
    """Apply the frozen v0.3 ranker without reading evaluator-only columns."""

    copy_path = Path(copy_feature_tsv_path)
    topology_path = Path(topology_tsv_path)
    model_path = Path(model_json_path)
    output_path = Path(output_tsv_path)
    copy_manifest_path = Path(str(copy_path) + ".manifest.json")
    topology_manifest_path = Path(str(topology_path) + ".manifest.json")
    output_manifest_path = Path(str(output_path) + ".manifest.json")
    for path in (
        copy_path,
        topology_path,
        model_path,
        copy_manifest_path,
        topology_manifest_path,
    ):
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f"Missing or empty v0.3 ranker input: {path}")
    collisions = [path for path in (output_path, output_manifest_path) if path.exists()]
    if collisions:
        raise FileExistsError(
            "Refusing to overwrite v0.3 score artifact(s): "
            + ", ".join(str(path) for path in collisions)
        )

    model = _load_json(model_path)
    validate_support_conditioned_ranker(model)
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
    correction_estimator = model["estimators"][
        "support_conditioned_topology_offset"
    ]
    scored: list[dict[str, Any]] = []
    primary_scores: list[float] = []
    available_count = 0
    for copy_row in copy_rows:
        topology_row = topology_by_digest[copy_row["candidate_digest"]]
        baseline_logit = _linear_score(
            transform_copy_feature_row(copy_row, contract), baseline_estimator
        )
        correction = _linear_score(
            support_conditioned_correction_vector(copy_row, topology_row),
            correction_estimator,
        )
        available = int(topology_row["topology_available"])
        if not available and correction != 0.0:
            raise AssertionError("Topology-unavailable candidate score changed")
        primary = baseline_logit + correction
        available_count += available
        primary_scores.append(primary)
        scored.append(
            {
                **copy_row,
                "v03_baseline_logit": format(baseline_logit, ".17g"),
                "v03_topology_correction": format(correction, ".17g"),
                "v03_primary_rank_score": format(primary, ".17g"),
                "v03_topology_available": available,
                "v03_automatic_approval": 0,
            }
        )

    percentiles = _percentile_ranks(primary_scores)
    output_fields = [
        *copy_fields,
        "v03_baseline_logit",
        "v03_topology_correction",
        "v03_primary_rank_score",
        "v03_primary_rank_percentile",
        "v03_topology_available",
        "v03_automatic_approval",
    ]
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=output_fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for row, percentile in zip(scored, percentiles, strict=True):
            row["v03_primary_rank_percentile"] = format(percentile, ".17g")
            writer.writerow(row)

    manifest: dict[str, Any] = {
        "schema_version": SUPPORT_CONDITIONED_SCORE_SCHEMA_VERSION,
        "generator": {"name": "PloidyPatch", "version": __version__},
        "truth_access": False,
        "inputs": {
            "copy_features": _file_sha256(copy_path),
            "topology_features": _file_sha256(topology_path),
            "model": _file_sha256(model_path),
        },
        "policy": {
            "primary_score": "v03_primary_rank_score",
            "interpretation": "uncalibrated_review_rank_not_probability",
            "automatic_approval": False,
            "topology_unavailable_policy": "exact_baseline_score",
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

