from __future__ import annotations

import csv
import json
import math
import shutil
from collections import Counter
from pathlib import Path
from statistics import median
from typing import Any, Iterable, Mapping, Sequence

from . import __version__
from .baseline import _file_sha256
from .copy_features import COPY_FEATURE_SCHEMA_VERSION, FEATURE_COLUMNS
from .wgd_candidate_select import _candidate_blocks


COPY_MODEL_SCHEMA_VERSION = "ploidypatch.copy_candidate_model.v1"
COPY_SCORE_SCHEMA_VERSION = "ploidypatch.copy_candidate_scores.v1"
COPY_MODEL_SELECTION_SCHEMA_VERSION = "ploidypatch.copy_model_selection.v1"

METHOD_PATTERNS = (
    "gemoma",
    "lifton",
    "miniprot",
    "gemoma,lifton",
    "gemoma,miniprot",
    "lifton,miniprot",
    "gemoma,lifton,miniprot",
)

STRUCTURAL_NUMERIC = (
    "span_bp",
    "cds_segments",
    "cds_bp",
)
SUPPORT_NUMERIC = ("support_method_count",)
METHOD_NUMERIC = (
    "miniprot_identity",
    "miniprot_query_coverage",
    "miniprot_rank",
    "miniprot_score",
    "miniprot_positive",
    "miniprot_frameshifts",
    "miniprot_stop_codons",
    "gemoma_pAA",
    "gemoma_iAA",
    "gemoma_score_ratio",
    "gemoma_complete_exon_count",
    "gemoma_reference_exon_count",
    "gemoma_nps",
    "lifton_dna_identity",
    "lifton_protein_identity",
)
WGD_NUMERIC = (
    "wgd_support_block_count",
    "wgd_longest_block_pairs",
)
METHOD_BINARY = (
    "has_miniprot",
    "has_gemoma",
    "has_lifton",
    "gemoma_start_complete",
    "gemoma_stop_complete",
    "lifton_frameshift",
    "lifton_stop_missing",
    "lifton_start_missing",
    "lifton_inframe_indel",
)
WGD_BINARY = ("wgd_existing_partner",)
LOG1P_FIELDS = frozenset(
    {
        "span_bp",
        "cds_segments",
        "cds_bp",
        "miniprot_rank",
        "miniprot_score",
        "miniprot_frameshifts",
        "miniprot_stop_codons",
        "gemoma_complete_exon_count",
        "gemoma_reference_exon_count",
        "gemoma_nps",
        "wgd_support_block_count",
        "wgd_longest_block_pairs",
    }
)
FEATURE_SETS = (
    "full",
    "no_wgd",
    "no_method_quality",
    "method_support_only",
)
FEATURE_MASK_GROUPS = (
    "wgd_context",
    "method_quality",
)


def _feature_set_fields(
    feature_set: str,
) -> tuple[tuple[str, ...], tuple[str, ...], bool, bool]:
    if feature_set not in FEATURE_SETS:
        raise ValueError(
            f"Unknown feature set {feature_set!r}; expected one of {FEATURE_SETS}"
        )
    if feature_set == "full":
        return (
            STRUCTURAL_NUMERIC + SUPPORT_NUMERIC + METHOD_NUMERIC + WGD_NUMERIC,
            METHOD_BINARY + WGD_BINARY,
            True,
            True,
        )
    if feature_set == "no_wgd":
        return (
            STRUCTURAL_NUMERIC + SUPPORT_NUMERIC + METHOD_NUMERIC,
            METHOD_BINARY,
            True,
            False,
        )
    if feature_set == "no_method_quality":
        return (
            STRUCTURAL_NUMERIC + SUPPORT_NUMERIC + WGD_NUMERIC,
            ("has_miniprot", "has_gemoma", "has_lifton") + WGD_BINARY,
            True,
            True,
        )
    return (
        SUPPORT_NUMERIC,
        ("has_miniprot", "has_gemoma", "has_lifton"),
        True,
        False,
    )


def _parse_float(value: str, field: str) -> float | None:
    if value == "":
        return None
    try:
        parsed = float(value)
    except ValueError as exc:
        raise ValueError(f"Non-numeric value for {field}: {value!r}") from exc
    if not math.isfinite(parsed):
        raise ValueError(f"Non-finite value for {field}: {value!r}")
    return parsed


def _transform_numeric(value: float, field: str) -> float:
    if field in LOG1P_FIELDS:
        if value < 0:
            raise ValueError(f"Negative value cannot be log1p transformed: {field}")
        return math.log1p(value)
    return value


def fit_copy_feature_contract(
    rows: Sequence[Mapping[str, str]], *, feature_set: str = "full"
) -> dict[str, Any]:
    """Fit deterministic imputation/scaling without importing an ML runtime."""

    if not rows:
        raise ValueError("Cannot fit a copy-feature contract on zero rows")
    numeric_fields, binary_fields, include_patterns, include_wgd_interactions = (
        _feature_set_fields(feature_set)
    )
    numeric: dict[str, dict[str, float | str]] = {}
    for field in numeric_fields:
        observed = [
            _transform_numeric(value, field)
            for row in rows
            if (value := _parse_float(row.get(field, ""), field)) is not None
        ]
        fill = float(median(observed)) if observed else 0.0
        filled = [
            fill
            if (value := _parse_float(row.get(field, ""), field)) is None
            else _transform_numeric(value, field)
            for row in rows
        ]
        mean = sum(filled) / len(filled)
        variance = sum((value - mean) ** 2 for value in filled) / len(filled)
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
        raise AssertionError("Expanded copy-model feature names are not unique")
    return {
        "feature_set": feature_set,
        "required_input_columns": list(FEATURE_COLUMNS),
        "numeric_order": list(numeric_fields),
        "numeric": numeric,
        "binary": list(binary_fields),
        "method_patterns": list(METHOD_PATTERNS) if include_patterns else [],
        "wgd_pattern_interactions": include_wgd_interactions,
        "expanded_feature_names": expanded,
        "missing_value_policy": "training_median_plus_explicit_missing_indicator",
        "scaling_policy": "training_population_mean_and_standard_deviation",
    }


def transform_copy_feature_row(
    row: Mapping[str, str], contract: Mapping[str, Any]
) -> list[float]:
    pattern = row.get("support_methods", "")
    patterns = tuple(contract.get("method_patterns", ()))
    if patterns and pattern not in patterns:
        raise ValueError(f"Unsupported method support pattern: {pattern!r}")
    values: list[float] = []
    numeric = contract.get("numeric", {})
    numeric_order = contract.get("numeric_order", ())
    if set(numeric_order) != set(numeric) or len(numeric_order) != len(numeric):
        raise ValueError("Feature contract numeric order disagrees with specifications")
    for field in numeric_order:
        specification = numeric[field]
        raw = _parse_float(row.get(field, ""), field)
        missing = raw is None
        transformed = (
            float(specification["impute_median"])
            if missing
            else _transform_numeric(raw, field)
        )
        scale = float(specification["scale"])
        if not math.isfinite(scale) or scale <= 0:
            raise ValueError(f"Invalid model scale for {field}")
        values.extend(
            (
                (transformed - float(specification["mean"])) / scale,
                float(missing),
            )
        )
    for field in contract.get("binary", ()):
        raw = _parse_float(row.get(field, ""), field)
        if raw is None:
            raw = 0.0
        if raw not in {0.0, 1.0}:
            raise ValueError(f"Binary feature {field} must be zero or one")
        values.append(raw)
    values.extend(float(pattern == level) for level in patterns)
    if contract.get("wgd_pattern_interactions", False):
        wgd = _parse_float(row.get("wgd_existing_partner", ""), "wgd_existing_partner")
        if wgd not in {0.0, 1.0}:
            raise ValueError("wgd_existing_partner must be zero or one")
        values.extend(wgd * float(pattern == level) for level in patterns)
    expected = contract.get("expanded_feature_names", ())
    if len(values) != len(expected):
        raise ValueError("Feature contract and transformed vector lengths disagree")
    return values


def _sigmoid(value: float) -> float:
    if value >= 0:
        exponential = math.exp(-value)
        return 1.0 / (1.0 + exponential)
    exponential = math.exp(value)
    return exponential / (1.0 + exponential)


def validate_copy_model(model: Mapping[str, Any]) -> None:
    if model.get("schema_version") != COPY_MODEL_SCHEMA_VERSION:
        raise ValueError("Unsupported copy-model schema version")
    contract = model.get("feature_contract")
    estimator = model.get("estimator")
    calibration = model.get("calibration")
    thresholds = model.get("thresholds")
    if not all(isinstance(value, Mapping) for value in (contract, estimator, calibration, thresholds)):
        raise ValueError("Copy model lacks a feature contract, estimator, calibration, or thresholds")
    names = contract.get("expanded_feature_names", ())
    coefficients = estimator.get("coefficients", ())
    if not names or len(coefficients) != len(names):
        raise ValueError("Copy-model coefficient length does not match its feature contract")
    for value in (*coefficients, estimator.get("intercept"), calibration.get("slope"), calibration.get("intercept")):
        if not isinstance(value, (int, float)) or not math.isfinite(float(value)):
            raise ValueError("Copy model contains a non-finite numeric parameter")
    review = thresholds.get("review", {})
    if not isinstance(review, Mapping) or not isinstance(review.get("value"), (int, float)):
        raise ValueError("Copy model lacks a numeric review threshold")
    if not 0 <= float(review["value"]) <= 1:
        raise ValueError("Copy-model review threshold lies outside [0, 1]")
    high = thresholds.get("high_confidence", {})
    if not isinstance(high, Mapping):
        raise ValueError("Copy model lacks a high-confidence threshold record")
    high_value = high.get("value")
    if high_value is not None and (
        not isinstance(high_value, (int, float)) or not 0 <= float(high_value) <= 1
    ):
        raise ValueError("Invalid high-confidence copy-model threshold")


def predict_copy_candidate(
    row: Mapping[str, str], model: Mapping[str, Any]
) -> dict[str, float]:
    validate_copy_model(model)
    vector = transform_copy_feature_row(row, model["feature_contract"])
    estimator = model["estimator"]
    logit = float(estimator["intercept"]) + sum(
        float(coefficient) * value
        for coefficient, value in zip(estimator["coefficients"], vector, strict=True)
    )
    raw_probability = _sigmoid(logit)
    calibration = model["calibration"]
    calibrated_logit = float(calibration["slope"]) * logit + float(
        calibration["intercept"]
    )
    return {
        "raw_logit": logit,
        "raw_probability": raw_probability,
        "calibrated_probability": _sigmoid(calibrated_logit),
    }


def _read_feature_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or not set(FEATURE_COLUMNS) <= set(reader.fieldnames):
            raise ValueError("Copy-feature table lacks the frozen input columns")
        rows = list(reader)
        fieldnames = list(reader.fieldnames)
    evaluator_columns = [
        field
        for field in fieldnames
        if field.startswith("label_") or field.startswith("truth_")
    ]
    if evaluator_columns:
        raise ValueError(
            "Refusing evaluator-only columns in truth-blind copy scoring: "
            + ", ".join(evaluator_columns)
        )
    digests = [row["candidate_digest"] for row in rows]
    if any(not digest for digest in digests) or len(digests) != len(set(digests)):
        raise ValueError("Copy-feature candidate digests must be nonempty and unique")
    return fieldnames, rows


def _mask_feature_row(
    row: Mapping[str, str], mask_feature_groups: Sequence[str]
) -> dict[str, str]:
    """Remove declared evidence groups without refitting the frozen model."""

    unknown = set(mask_feature_groups) - set(FEATURE_MASK_GROUPS)
    if unknown:
        raise ValueError(
            "Unknown copy-model feature mask(s): " + ", ".join(sorted(unknown))
        )
    masked = dict(row)
    if "wgd_context" in mask_feature_groups:
        for field in WGD_NUMERIC:
            masked[field] = ""
        for field in WGD_BINARY:
            masked[field] = "0"
    if "method_quality" in mask_feature_groups:
        for field in METHOD_NUMERIC:
            masked[field] = ""
        for field in METHOD_BINARY:
            if not field.startswith("has_"):
                masked[field] = "0"
    return masked


def score_copy_candidate_features(
    *,
    feature_tsv_path: str | Path,
    model_json_path: str | Path,
    output_tsv_path: str | Path,
    mask_feature_groups: Iterable[str] = (),
) -> dict[str, Any]:
    """Score truth-blind candidate features using a frozen portable JSON model."""

    feature_path = Path(feature_tsv_path)
    masks = tuple(mask_feature_groups)
    if len(masks) != len(set(masks)):
        raise ValueError("Copy-model feature masks must be unique")
    unknown_masks = set(masks) - set(FEATURE_MASK_GROUPS)
    if unknown_masks:
        raise ValueError(
            "Unknown copy-model feature mask(s): "
            + ", ".join(sorted(unknown_masks))
        )
    feature_manifest_path = Path(str(feature_path) + ".manifest.json")
    model_path = Path(model_json_path)
    output_path = Path(output_tsv_path)
    manifest_path = Path(str(output_path) + ".manifest.json")
    for required in (feature_path, feature_manifest_path, model_path):
        if not required.is_file() or required.stat().st_size == 0:
            raise ValueError(f"Missing or empty copy-model input: {required}")
    collisions = [path for path in (output_path, manifest_path) if path.exists()]
    if collisions:
        raise FileExistsError(
            "Refusing to overwrite copy-score artifact(s): "
            + ", ".join(str(path) for path in collisions)
        )
    with model_path.open("r", encoding="utf-8") as handle:
        model = json.load(handle)
    validate_copy_model(model)
    with feature_manifest_path.open("r", encoding="utf-8") as handle:
        feature_manifest = json.load(handle)
    if (
        feature_manifest.get("schema_version") != COPY_FEATURE_SCHEMA_VERSION
        or feature_manifest.get("truth_access") is not False
        or feature_manifest.get("outputs", {}).get("features", {}).get("sha256")
        != _file_sha256(feature_path)
    ):
        raise ValueError("Copy-feature manifest fails schema, truth, or checksum gate")
    fieldnames, rows = _read_feature_rows(feature_path)
    score_fields = (
        "model_raw_logit",
        "model_raw_probability",
        "model_calibrated_probability",
        "model_review_decision",
        "model_high_confidence_decision",
    )
    if set(score_fields) & set(fieldnames):
        raise ValueError("Copy-feature table already contains model score columns")
    review_threshold = float(model["thresholds"]["review"]["value"])
    high_threshold = model["thresholds"]["high_confidence"].get("value")
    decisions: Counter[str] = Counter()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames + list(score_fields),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in rows:
            prediction = predict_copy_candidate(
                _mask_feature_row(row, masks), model
            )
            probability = prediction["calibrated_probability"]
            review = int(probability >= review_threshold)
            high = int(high_threshold is not None and probability >= float(high_threshold))
            decisions["review"] += review
            decisions["high_confidence"] += high
            writer.writerow(
                {
                    **row,
                    "model_raw_logit": format(prediction["raw_logit"], ".17g"),
                    "model_raw_probability": format(
                        prediction["raw_probability"], ".17g"
                    ),
                    "model_calibrated_probability": format(probability, ".17g"),
                    "model_review_decision": review,
                    "model_high_confidence_decision": high,
                }
            )
    manifest: dict[str, Any] = {
        "schema_version": COPY_SCORE_SCHEMA_VERSION,
        "generator": {"name": "PloidyPatch", "version": __version__},
        "truth_access": False,
        "inputs": {
            "features": {"file_name": feature_path.name, "sha256": _file_sha256(feature_path)},
            "feature_manifest": {
                "file_name": feature_manifest_path.name,
                "sha256": _file_sha256(feature_manifest_path),
            },
            "model": {"file_name": model_path.name, "sha256": _file_sha256(model_path)},
        },
        "model_schema_version": model["schema_version"],
        "counterfactual_feature_masks": list(masks),
        "interpretation": (
            "same_frozen_model_and_thresholds_with_declared_input_evidence_removed"
            if masks
            else "unmasked_frozen_model"
        ),
        "thresholds": model["thresholds"],
        "counts": {
            "candidates": len(rows),
            "review_selected": decisions["review"],
            "high_confidence_selected": decisions["high_confidence"],
        },
        "outputs": {
            "scores": {
                "file_name": output_path.name,
                "sha256": _file_sha256(output_path),
                "rows": len(rows),
            }
        },
    }
    with manifest_path.open("x", encoding="utf-8", newline="") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return manifest


def select_scored_copy_candidates(
    *,
    base_gff_path: str | Path,
    candidate_gff_path: str | Path,
    scored_tsv_path: str | Path,
    model_json_path: str | Path,
    output_gff_path: str | Path,
    selection_tsv_path: str | Path,
    policy: str = "review",
) -> dict[str, Any]:
    """Retain complete consensus hierarchies selected by a frozen model policy."""

    if policy not in {"review", "high_confidence"}:
        raise ValueError("Copy-model selection policy must be review or high_confidence")
    base_path = Path(base_gff_path)
    candidate_path = Path(candidate_gff_path)
    scored_path = Path(scored_tsv_path)
    scored_manifest_path = Path(str(scored_path) + ".manifest.json")
    model_path = Path(model_json_path)
    output_path = Path(output_gff_path)
    selection_path = Path(selection_tsv_path)
    manifest_path = Path(str(output_path) + ".manifest.json")
    for required in (
        base_path,
        candidate_path,
        scored_path,
        scored_manifest_path,
        model_path,
    ):
        if not required.is_file() or required.stat().st_size == 0:
            raise ValueError(f"Missing or empty copy-model selection input: {required}")
    collisions = [
        path for path in (output_path, selection_path, manifest_path) if path.exists()
    ]
    if collisions:
        raise FileExistsError(
            "Refusing to overwrite copy-model selection artifact(s): "
            + ", ".join(str(path) for path in collisions)
        )
    with model_path.open("r", encoding="utf-8") as handle:
        model = json.load(handle)
    validate_copy_model(model)
    with scored_manifest_path.open("r", encoding="utf-8") as handle:
        scored_manifest = json.load(handle)
    if (
        scored_manifest.get("schema_version") != COPY_SCORE_SCHEMA_VERSION
        or scored_manifest.get("truth_access") is not False
        or scored_manifest.get("inputs", {}).get("model", {}).get("sha256")
        != _file_sha256(model_path)
        or scored_manifest.get("outputs", {}).get("scores", {}).get("sha256")
        != _file_sha256(scored_path)
    ):
        raise ValueError("Copy-score manifest fails schema, model, truth, or checksum gate")
    threshold_record = model["thresholds"][policy]
    if threshold_record.get("value") is None:
        raise ValueError(f"Frozen model has no available {policy} threshold")

    blocks, gene_digests = _candidate_blocks(candidate_path, base_path)
    digest_to_gene = {digest: gene for gene, digest in gene_digests.items()}
    if len(digest_to_gene) != len(gene_digests):
        raise ValueError("Candidate GFF contains duplicate consensus digests")
    with scored_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        decision_field = f"model_{policy}_decision"
        required_fields = {
            "candidate_digest",
            "model_calibrated_probability",
            decision_field,
        }
        if reader.fieldnames is None or not required_fields <= set(reader.fieldnames):
            raise ValueError("Scored copy table lacks selection fields")
        scored_rows = list(reader)
    scored = {row["candidate_digest"]: row for row in scored_rows}
    if len(scored) != len(scored_rows) or any(not digest for digest in scored):
        raise ValueError("Scored copy table contains empty or duplicate digests")
    if set(scored) != set(digest_to_gene):
        raise ValueError("Scored table and candidate GFF digest universes differ")

    selected_digests: set[str] = set()
    decision_field = f"model_{policy}_decision"
    threshold = float(threshold_record["value"])
    for digest, row in scored.items():
        decision = row[decision_field]
        if decision not in {"0", "1"}:
            raise ValueError(f"Invalid scored decision for {digest}")
        probability = _parse_float(
            row["model_calibrated_probability"], "model_calibrated_probability"
        )
        if probability is None or not 0 <= probability <= 1:
            raise ValueError(f"Invalid scored probability for {digest}")
        if int(probability >= threshold) != int(decision):
            raise ValueError(f"Scored decision disagrees with frozen threshold for {digest}")
        if decision == "1":
            selected_digests.add(digest)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with base_path.open("rb") as source_handle, output_path.open("xb") as output_handle:
        shutil.copyfileobj(source_handle, output_handle, length=8 * 1024 * 1024)
        source_handle.seek(0, 2)
        if source_handle.tell():
            source_handle.seek(-1, 2)
            if source_handle.read(1) not in {b"\n", b"\r"}:
                output_handle.write(b"\n")
        output_handle.write(b"###\n")
        for digest in sorted(selected_digests):
            output_handle.write("".join(blocks[digest_to_gene[digest]]).encode("utf-8"))

    fields = (
        "gene_id",
        "candidate_digest",
        "calibrated_probability",
        "policy",
        "threshold",
        "status",
        "reason",
    )
    with selection_path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for digest in sorted(scored):
            accepted = digest in selected_digests
            writer.writerow(
                {
                    "gene_id": digest_to_gene[digest],
                    "candidate_digest": digest,
                    "calibrated_probability": scored[digest][
                        "model_calibrated_probability"
                    ],
                    "policy": policy,
                    "threshold": format(float(threshold_record["value"]), ".17g"),
                    "status": "accepted" if accepted else "rejected",
                    "reason": (
                        "frozen_model_threshold_pass"
                        if accepted
                        else "frozen_model_threshold_fail"
                    ),
                }
            )
    manifest: dict[str, Any] = {
        "schema_version": COPY_MODEL_SELECTION_SCHEMA_VERSION,
        "generator": {"name": "PloidyPatch", "version": __version__},
        "truth_access": False,
        "inputs": {
            "base_gff": {"file_name": base_path.name, "sha256": _file_sha256(base_path)},
            "candidate_gff": {
                "file_name": candidate_path.name,
                "sha256": _file_sha256(candidate_path),
            },
            "scores": {"file_name": scored_path.name, "sha256": _file_sha256(scored_path)},
            "score_manifest": {
                "file_name": scored_manifest_path.name,
                "sha256": _file_sha256(scored_manifest_path),
            },
            "model": {"file_name": model_path.name, "sha256": _file_sha256(model_path)},
        },
        "policy": {
            "name": policy,
            "threshold": threshold_record,
            "claim_boundary": "review_candidate_ranking_not_automatic_annotation_approval",
        },
        "counts": {
            "candidate_models": len(scored),
            "selected_models": len(selected_digests),
            "rejected_models": len(scored) - len(selected_digests),
        },
        "outputs": {
            "candidate_gff": {"file_name": output_path.name, "sha256": _file_sha256(output_path)},
            "selection": {
                "file_name": selection_path.name,
                "sha256": _file_sha256(selection_path),
                "rows": len(scored),
            },
        },
    }
    with manifest_path.open("x", encoding="utf-8", newline="") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return manifest
