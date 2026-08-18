from __future__ import annotations

import math
from typing import Any, Mapping, Sequence

from .copy_features import FEATURE_COLUMNS
from .copy_model import (
    LOG1P_FIELDS,
    METHOD_PATTERNS,
    _feature_set_fields,
    _parse_float,
    _transform_numeric,
    fit_copy_feature_contract,
)


STABLE_REFERENCE_RANKER_SCHEMA_VERSION = (
    "ploidypatch.stable_reference_ranker.v0.9"
)


def _validated_weights(
    rows: Sequence[Mapping[str, str]], weights: Sequence[float]
) -> tuple[float, ...]:
    if not rows:
        raise ValueError("Cannot fit a stable feature contract on zero rows")
    if len(rows) != len(weights):
        raise ValueError("Feature rows and sample weights differ in length")
    parsed = tuple(float(value) for value in weights)
    if any(not math.isfinite(value) or value <= 0 for value in parsed):
        raise ValueError("Stable feature-contract weights must be finite and positive")
    return parsed


def _weighted_median(values: Sequence[float], weights: Sequence[float]) -> float:
    """Return a deterministic weighted median.

    When every weight is equal this definition is byte-for-byte compatible with
    ``statistics.median``: an exact half-weight boundary is represented by the
    mean of the adjacent observations.
    """

    if not values or len(values) != len(weights):
        raise ValueError("Weighted median requires equally sized nonempty inputs")
    ordered = sorted(zip(values, weights, strict=True), key=lambda item: item[0])
    half = math.fsum(weight for _, weight in ordered) / 2.0
    cumulative = 0.0
    for index, (value, weight) in enumerate(ordered):
        cumulative = math.fsum((cumulative, weight))
        if cumulative > half:
            return float(value)
        if cumulative == half:
            if index + 1 == len(ordered):
                return float(value)
            return (float(value) + float(ordered[index + 1][0])) / 2.0
    raise AssertionError("Positive weighted median did not cross half mass")


def fit_weighted_copy_feature_contract(
    rows: Sequence[Mapping[str, str]],
    *,
    sample_weight: Sequence[float],
    feature_set: str = "full",
) -> dict[str, Any]:
    """Fit the frozen copy-feature transform with explicit row weights.

    Stable v0.9 gives each development species equal total mass.  The weights
    affect imputation, centering, and scaling only; feature definitions and
    expansion order remain the existing copy-feature v1 contract.
    """

    weights = _validated_weights(rows, sample_weight)
    if all(weight == weights[0] for weight in weights):
        contract = fit_copy_feature_contract(rows, feature_set=feature_set)
        contract["missing_value_policy"] = (
            "weighted_training_median_plus_explicit_missing_indicator"
        )
        contract["scaling_policy"] = (
            "weighted_training_population_mean_and_standard_deviation"
        )
        contract["weight_policy"] = "caller_supplied_positive_weights"
        return contract
    total_weight = math.fsum(weights)
    numeric_fields, binary_fields, include_patterns, include_wgd_interactions = (
        _feature_set_fields(feature_set)
    )
    numeric: dict[str, dict[str, float | str]] = {}
    for field in numeric_fields:
        observed_values: list[float] = []
        observed_weights: list[float] = []
        parsed_values: list[float | None] = []
        for row, weight in zip(rows, weights, strict=True):
            parsed = _parse_float(row.get(field, ""), field)
            transformed = (
                None if parsed is None else _transform_numeric(parsed, field)
            )
            parsed_values.append(transformed)
            if transformed is not None:
                observed_values.append(transformed)
                observed_weights.append(weight)
        fill = (
            _weighted_median(observed_values, observed_weights)
            if observed_values
            else 0.0
        )
        filled = [fill if value is None else value for value in parsed_values]
        mean = math.fsum(
            value * weight
            for value, weight in zip(filled, weights, strict=True)
        ) / total_weight
        variance = math.fsum(
            weight * (value - mean) ** 2
            for value, weight in zip(filled, weights, strict=True)
        ) / total_weight
        scale = math.sqrt(variance)
        if scale == 0:
            scale = 1.0
        numeric[field] = {
            "transform": "log1p" if field in LOG1P_FIELDS else "identity",
            "impute_median": fill,
            "mean": mean,
            "scale": scale,
        }

    expanded: list[str] = []
    for field in numeric_fields:
        expanded.extend((f"numeric:{field}", f"missing:{field}"))
    expanded.extend(f"binary:{field}" for field in binary_fields)
    if include_patterns:
        expanded.extend(f"pattern:{pattern}" for pattern in METHOD_PATTERNS)
    if include_wgd_interactions:
        expanded.extend(f"wgd_pattern:{pattern}" for pattern in METHOD_PATTERNS)
    if len(expanded) != len(set(expanded)):
        raise AssertionError("Expanded stable feature names are not unique")
    return {
        "feature_set": feature_set,
        "required_input_columns": list(FEATURE_COLUMNS),
        "numeric_order": list(numeric_fields),
        "numeric": numeric,
        "binary": list(binary_fields),
        "method_patterns": list(METHOD_PATTERNS) if include_patterns else [],
        "wgd_pattern_interactions": include_wgd_interactions,
        "expanded_feature_names": expanded,
        "missing_value_policy": (
            "weighted_training_median_plus_explicit_missing_indicator"
        ),
        "scaling_policy": (
            "weighted_training_population_mean_and_standard_deviation"
        ),
        "weight_policy": "caller_supplied_positive_weights",
    }
