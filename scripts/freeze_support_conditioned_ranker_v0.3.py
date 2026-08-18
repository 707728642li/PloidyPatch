#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import platform
import sys
from pathlib import Path
from typing import Any

import numpy as np
import scipy
import sklearn
from scipy.optimize import minimize
from scipy.special import expit
from sklearn.linear_model import LogisticRegression

from ploidypatch import __version__ as ploidypatch_version
from ploidypatch.copy_features import FEATURE_COLUMNS
from ploidypatch.copy_model import (
    fit_copy_feature_contract,
    transform_copy_feature_row,
)
from ploidypatch.support_ranker import (
    CORRECTION_FEATURE_NAMES,
    SUPPORT_CONDITIONED_RANKER_SCHEMA_VERSION,
    SUPPORT_PATTERNS,
    support_conditioned_correction_vector,
    validate_support_conditioned_ranker,
)


FIXED_C = 1.0
CANONICAL_EVALUATION_SCHEMA = "ploidypatch.candidate_pool_ranker_evaluation.v4"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def read_index(path: Path) -> tuple[list[str], dict[str, dict[str, str]]]:
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


def parse_development(value: str) -> dict[str, Any]:
    name, separator, payload = value.partition("=")
    parts = payload.split(",")
    if not separator or not name or len(parts) != 2:
        raise ValueError("--development requires NAME=LABELED_FEATURES,TOPOLOGY")
    feature_path, topology_path = (Path(part) for part in parts)
    if not feature_path.is_file() or not topology_path.is_file():
        raise FileNotFoundError(f"Missing development input for {name}")
    feature_fields, features = read_index(feature_path)
    _, topology = read_index(topology_path)
    if not {*FEATURE_COLUMNS, "label_exact_cds"} <= set(feature_fields):
        raise ValueError(f"Development feature contract incomplete for {name}")
    if set(features) != set(topology):
        raise ValueError(f"Feature/topology candidate universes differ for {name}")
    digests = sorted(features)
    rows = [features[digest] for digest in digests]
    labels = np.asarray([int(row["label_exact_cds"]) for row in rows], dtype=float)
    if not 0 < int(labels.sum()) < len(labels):
        raise ValueError(f"Development dataset lacks both classes: {name}")
    patterns = {row["support_methods"] for row in rows}
    if patterns - set(SUPPORT_PATTERNS):
        raise ValueError(f"Unexpected support pattern for {name}")
    correction = np.asarray(
        [
            support_conditioned_correction_vector(
                features[digest], topology[digest]
            )
            for digest in digests
        ],
        dtype=float,
    )
    return {
        "name": name,
        "paths": (feature_path, topology_path),
        "rows": rows,
        "labels": labels,
        "correction": correction,
    }


def fit_offset(
    offset: np.ndarray, correction: np.ndarray, labels: np.ndarray
) -> np.ndarray:
    def objective(beta: np.ndarray) -> tuple[float, np.ndarray]:
        logits = offset + correction @ beta
        loss = float(
            np.sum(np.logaddexp(0.0, logits) - labels * logits)
            + 0.5 / FIXED_C * np.dot(beta, beta)
        )
        gradient = correction.T @ (expit(logits) - labels) + beta / FIXED_C
        return loss, gradient

    result = minimize(
        objective,
        np.zeros(correction.shape[1], dtype=float),
        method="L-BFGS-B",
        jac=True,
        options={"maxiter": 5000, "ftol": 1e-12, "gtol": 1e-8},
    )
    if not result.success or not np.isfinite(result.x).all():
        raise RuntimeError(f"Offset optimization failed: {result.message}")
    return np.asarray(result.x, dtype=float)


def maximum_absolute_delta(observed: np.ndarray, expected: list[float]) -> float:
    reference = np.asarray(expected, dtype=float)
    if observed.shape != reference.shape:
        raise ValueError("Canonical parameter vector length changed")
    return float(np.max(np.abs(observed - reference))) if len(observed) else 0.0


def offset_objective(
    offset: np.ndarray,
    correction: np.ndarray,
    labels: np.ndarray,
    beta: np.ndarray,
) -> float:
    logits = offset + correction @ beta
    return float(
        np.sum(np.logaddexp(0.0, logits) - labels * logits)
        + 0.5 / FIXED_C * np.dot(beta, beta)
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--development", action="append", required=True)
    parser.add_argument("--canonical-evaluation", required=True)
    parser.add_argument("--protocol", required=True)
    parser.add_argument("--output-json", required=True)
    args = parser.parse_args()

    if len(args.development) != 2:
        raise ValueError("The frozen v0.3 model requires exactly two development species")
    datasets = [parse_development(value) for value in args.development]
    if [dataset["name"] for dataset in datasets] != ["glycine", "brassica"]:
        raise ValueError("Development species order must remain glycine, brassica")
    evaluation_path = Path(args.canonical_evaluation)
    protocol_path = Path(args.protocol)
    output_path = Path(args.output_json)
    for path in (evaluation_path, protocol_path):
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(path)
    if output_path.exists():
        raise FileExistsError(output_path)
    with evaluation_path.open("r", encoding="utf-8") as handle:
        evaluation = json.load(handle)
    if (
        evaluation.get("schema_version") != CANONICAL_EVALUATION_SCHEMA
        or evaluation.get("primary_retention_criteria", {}).get(
            "retain_for_external_freeze"
        )
        is not True
        or evaluation.get("evaluation_contract", {}).get("primary_estimator")
        != "offset_support_conditioned"
    ):
        raise ValueError("Canonical evaluation did not authorize v0.3 freeze")

    rows = [row for dataset in datasets for row in dataset["rows"]]
    labels = np.concatenate([dataset["labels"] for dataset in datasets])
    correction = np.concatenate([dataset["correction"] for dataset in datasets])
    contract = fit_copy_feature_contract(rows, feature_set="full")
    base = np.asarray(
        [transform_copy_feature_row(row, contract) for row in rows], dtype=float
    )
    baseline = LogisticRegression(
        C=FIXED_C, solver="lbfgs", max_iter=5000, random_state=0
    ).fit(base, labels)
    beta = fit_offset(baseline.decision_function(base), correction, labels)

    canonical = evaluation.get("pooled_fit_artifacts")
    if not isinstance(canonical, dict):
        raise ValueError("Canonical evaluation lacks pooled fit artifacts")
    if canonical.get("base_feature_names") != contract["expanded_feature_names"]:
        raise ValueError("Refitted base-feature order differs from evaluation")
    intercept_delta = abs(
        float(baseline.intercept_[0]) - float(canonical["baseline_intercept"])
    )
    coefficient_delta = maximum_absolute_delta(
        baseline.coef_[0], canonical["baseline_coefficients"]
    )
    canonical_beta = np.asarray(
        canonical["offset_support_coefficients"], dtype=float
    )
    correction_delta = maximum_absolute_delta(
        beta, canonical["offset_support_coefficients"]
    )
    base_tolerance = 1e-12
    if max(intercept_delta, coefficient_delta) > base_tolerance:
        raise ValueError(
            "Refitted baseline does not numerically match canonical evaluation: "
            f"baseline_intercept={intercept_delta:.17g}, "
            f"baseline_coefficients={coefficient_delta:.17g}, "
            f"tolerance={base_tolerance:.17g}"
        )
    offset = baseline.decision_function(base)
    canonical_scores = offset + correction @ canonical_beta
    refit_scores = offset + correction @ beta
    score_delta = float(np.max(np.abs(refit_scores - canonical_scores)))
    canonical_order = np.argsort(canonical_scores, kind="mergesort")
    refit_order = np.argsort(refit_scores, kind="mergesort")
    rank_order_equal = bool(np.array_equal(canonical_order, refit_order))
    objective_delta = abs(
        offset_objective(offset, correction, labels, beta)
        - offset_objective(offset, correction, labels, canonical_beta)
    )
    refit_score_tolerance = 1e-6
    refit_objective_tolerance = 1e-8
    if (
        correction_delta > 1e-6
        or score_delta > refit_score_tolerance
        or objective_delta > refit_objective_tolerance
        or not rank_order_equal
    ):
        raise ValueError(
            "Independent topology-offset refit does not reproduce the canonical "
            "ranker: "
            f"coefficient_delta={correction_delta:.17g}, "
            f"score_delta={score_delta:.17g}, "
            f"objective_delta={objective_delta:.17g}, "
            f"rank_order_equal={rank_order_equal}"
        )

    model: dict[str, Any] = {
        "schema_version": SUPPORT_CONDITIONED_RANKER_SCHEMA_VERSION,
        "model_id": "PloidyPatch_support_conditioned_ranker_v0.3",
        "generator": {
            "name": "PloidyPatch",
            "version": ploidypatch_version,
            "code_commit": os.environ.get(
                "PLOIDYPATCH_CODE_COMMIT", "unavailable_server_mirror"
            ),
        },
        "selected_estimator": "baseline_logit_plus_support_conditioned_topology_offset",
        "fixed_regularization_C": FIXED_C,
        "base_feature_contract": contract,
        "support_patterns": list(SUPPORT_PATTERNS),
        "topology_component_fields": [
            name.removeprefix("topology:")
            for name in CORRECTION_FEATURE_NAMES[:5]
        ],
        "estimators": {
            "baseline": {
                "type": "logistic_regression_logit",
                "intercept": float(baseline.intercept_[0]),
                "coefficient_feature_order": contract[
                    "expanded_feature_names"
                ],
                "coefficients": baseline.coef_[0].tolist(),
            },
            "support_conditioned_topology_offset": {
                "type": "regularized_logistic_offset_correction",
                "intercept": 0.0,
                "coefficient_feature_order": list(CORRECTION_FEATURE_NAMES),
                "coefficients": canonical_beta.tolist(),
                "topology_unavailable_correction": 0.0,
                "parameter_source": "canonical_development_evaluation",
            },
        },
        "claim_boundary": {
            "score": "uncalibrated_review_rank",
            "calibrated_probability": False,
            "automatic_approval": False,
            "portable_threshold": False,
        },
        "development_training": {
            "species_order": [dataset["name"] for dataset in datasets],
            "candidates": int(len(labels)),
            "exact_positives": int(labels.sum()),
            "inputs": [
                {
                    "species": dataset["name"],
                    "candidates": len(dataset["rows"]),
                    "exact_positives": int(dataset["labels"].sum()),
                    "files": [
                        {
                            "path": str(path),
                            "bytes": path.stat().st_size,
                            "sha256": sha256(path),
                        }
                        for path in dataset["paths"]
                    ],
                }
                for dataset in datasets
            ],
        },
        "canonical_evaluation": {
            "path": str(evaluation_path),
            "sha256": sha256(evaluation_path),
            "schema_version": evaluation["schema_version"],
            "evaluation_code_commit": evaluation["code_commit"],
            "retain_for_external_freeze": True,
            "baseline_parameter_match_tolerance": base_tolerance,
            "independent_offset_refit_score_tolerance": refit_score_tolerance,
            "independent_offset_refit_objective_tolerance": (
                refit_objective_tolerance
            ),
            "maximum_absolute_deltas": {
                "baseline_intercept": intercept_delta,
                "baseline_coefficients": coefficient_delta,
                "independent_refit_correction_coefficients": correction_delta,
                "independent_refit_candidate_scores": score_delta,
                "independent_refit_objective": objective_delta,
            },
            "independent_refit_rank_order_equal": rank_order_equal,
        },
        "protocol": {
            "path": str(protocol_path),
            "sha256": sha256(protocol_path),
        },
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "scikit_learn": sklearn.__version__,
        },
    }
    validate_support_conditioned_ranker(model)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("x", encoding="utf-8", newline="") as handle:
        json.dump(model, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(json.dumps({"model": str(output_path), "sha256": sha256(output_path)}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
