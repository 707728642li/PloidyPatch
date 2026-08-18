#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import platform
import sys
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np
import sklearn
from sklearn.linear_model import LogisticRegression

from ploidypatch import __version__ as ploidypatch_version
from ploidypatch.copy_features import FEATURE_COLUMNS
from ploidypatch.copy_model import (
    fit_copy_feature_contract,
    transform_copy_feature_row,
)
from ploidypatch.homeolog_ranker import (
    HOMEOLOG_RANKER_SCHEMA_VERSION,
    TOPOLOGY_ADDON_FIELDS,
    _read_topology_rows,
    topology_addons,
    validate_homeolog_ranker,
)


TRAINING_SCHEMA_VERSION = "ploidypatch.homeolog_copy_ranker_training.v2"
FIXED_C = 1.0


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def _write_json(path: Path, value: Any) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")


def _parse_dataset(value: str) -> tuple[str, Path, Path]:
    name, separator, payload = value.partition("=")
    paths = payload.split(",")
    if not separator or not name or len(paths) != 2 or any(not item for item in paths):
        raise ValueError("--dataset requires NAME=LABELED_FEATURES,TOPOLOGY_FEATURES")
    if any(character not in "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-" for character in name):
        raise ValueError(f"Unsafe dataset name: {name!r}")
    return name, Path(paths[0]), Path(paths[1])


def _read_labeled_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {*FEATURE_COLUMNS, "label_exact_cds"}
        if reader.fieldnames is None or not required <= set(reader.fieldnames):
            raise ValueError(f"Labeled feature table lacks frozen columns: {path}")
        rows = list(reader)
    digests = [row["candidate_digest"] for row in rows]
    if any(not digest for digest in digests) or len(digests) != len(set(digests)):
        raise ValueError(f"Empty or duplicate labeled digest: {path}")
    if any(row["label_exact_cds"] not in {"0", "1"} for row in rows):
        raise ValueError(f"Non-binary exact-CDS label: {path}")
    return rows


def _estimator_record(
    model: LogisticRegression, feature_names: list[str]
) -> dict[str, Any]:
    return {
        "family": "logistic_regression_l2",
        "solver": "lbfgs",
        "C": FIXED_C,
        "random_state": 0,
        "intercept": float(model.intercept_[0]),
        "coefficients": [float(value) for value in model.coef_[0]],
        "coefficient_feature_order": feature_names,
    }


def train(dataset_values: list[str], output_dir: str) -> None:
    parsed = [_parse_dataset(value) for value in dataset_values]
    names = [name for name, _, _ in parsed]
    if len(names) < 2 or len(names) != len(set(names)):
        raise ValueError("Training requires at least two uniquely named species datasets")
    output_path = Path(output_dir)
    partial_path = Path(str(output_path) + ".partial")
    if output_path.exists() or partial_path.exists():
        raise FileExistsError("Refusing to overwrite model output or partial directory")

    all_rows: list[dict[str, str]] = []
    all_addons: list[list[float]] = []
    all_labels: list[int] = []
    dataset_records: list[dict[str, Any]] = []
    for name, labeled_path, topology_path in parsed:
        if not labeled_path.is_file() or not topology_path.is_file():
            raise FileNotFoundError(f"Missing training input for {name}")
        rows = _read_labeled_rows(labeled_path)
        _, topology_rows = _read_topology_rows(topology_path)
        topology_by_digest = {row["candidate_digest"]: row for row in topology_rows}
        digests = [row["candidate_digest"] for row in rows]
        if set(digests) != set(topology_by_digest):
            raise ValueError(f"Candidate universes differ for {name}")
        labels = [int(row["label_exact_cds"]) for row in rows]
        addons = [topology_addons(topology_by_digest[digest]) for digest in digests]
        available = sum(int(values[0]) for values in addons)
        all_rows.extend(rows)
        all_addons.extend(addons)
        all_labels.extend(labels)
        dataset_records.append(
            {
                "name": name,
                "labeled_features": {
                    "path": str(labeled_path),
                    "bytes": labeled_path.stat().st_size,
                    "sha256": _sha256(labeled_path),
                },
                "topology_features": {
                    "path": str(topology_path),
                    "bytes": topology_path.stat().st_size,
                    "sha256": _sha256(topology_path),
                },
                "rows": len(rows),
                "positive_exact_cds": sum(labels),
                "topology_available": available,
                "topology_positive_exact_cds": sum(
                    label for label, values in zip(labels, addons, strict=True) if values[0]
                ),
                "chromosomes": sorted({row["seqid"] for row in rows}),
            }
        )
    labels_array = np.asarray(all_labels, dtype=np.uint8)
    if int(np.sum(labels_array)) < 100 or int(np.sum(labels_array == 0)) < 100:
        raise ValueError("Training requires at least 100 positive and 100 negative rows")

    contract = fit_copy_feature_contract(all_rows, feature_set="full")
    base_matrix = np.asarray(
        [transform_copy_feature_row(row, contract) for row in all_rows], dtype=float
    )
    addon_matrix = np.asarray(all_addons, dtype=float)
    topology_matrix = np.column_stack((base_matrix, addon_matrix))
    baseline = LogisticRegression(
        C=FIXED_C, solver="lbfgs", max_iter=5000, random_state=0
    ).fit(base_matrix, labels_array)
    topology = LogisticRegression(
        C=FIXED_C, solver="lbfgs", max_iter=5000, random_state=0
    ).fit(topology_matrix, labels_array)
    base_names = list(contract["expanded_feature_names"])
    topology_names = [*base_names, *TOPOLOGY_ADDON_FIELDS]
    model: dict[str, Any] = {
        "schema_version": HOMEOLOG_RANKER_SCHEMA_VERSION,
        "generator": {
            "name": "PloidyPatch",
            "version": ploidypatch_version,
            "training_schema_version": TRAINING_SCHEMA_VERSION,
            "code_commit": os.environ.get(
                "PLOIDYPATCH_CODE_COMMIT", "unavailable_server_mirror"
            ),
        },
        "development_scope": {
            "species": names,
            "task": "plant_annotation_copy_collapse_candidate_ranking",
            "candidate_universe": "exact_phased_CDS_union_of_miniprot_GeMoMa_LiftOn",
            "truth_use": "training_only_after_truth_blind_features_were_frozen",
        },
        "base_feature_contract": contract,
        "topology_addon_fields": list(TOPOLOGY_ADDON_FIELDS),
        "estimators": {
            "baseline": _estimator_record(baseline, base_names),
            "topology": _estimator_record(topology, topology_names),
        },
        "primary_estimator": "topology",
        "score_contract": {
            "score": "sigmoid_of_raw_logit_for_monotonic_ranking_only",
            "calibrated_probability": False,
            "portable_threshold": None,
            "within_run_rank_percentile": True,
            "topology_missing_encoding": "availability_zero_and_components_zero",
        },
        "feature_selection_rationale": {
            "included": "existing_WGD_partner_CDS_topology",
            "excluded_from_primary": "within_species_normalized_WGD_block_context",
            "reason": "topology_has_supported_mechanistic_gain_in_both_development_species_whereas_normalized_WGD_context_is_inconsistent",
        },
        "claim_boundary": {
            "automatic_approval": False,
            "output_role": "ranked_review_candidates",
            "topology_unavailable": "retain_generic_rank_and_flag_unavailable",
            "external_v2_holdout_status": "not_evaluated_when_model_was_frozen",
        },
    }
    validate_homeolog_ranker(model)

    partial_path.mkdir(parents=True)
    try:
        _write_json(partial_path / "model.json", model)
        report = {
            "schema_version": TRAINING_SCHEMA_VERSION,
            "datasets": dataset_records,
            "pooled_counts": {
                "rows": len(all_rows),
                "positive_exact_cds": int(np.sum(labels_array)),
                "negative_candidates": int(np.sum(labels_array == 0)),
                "topology_available": int(np.sum(addon_matrix[:, 0])),
            },
            "fixed_training_protocol": {
                "estimator": "logistic_regression_l2",
                "C": FIXED_C,
                "solver": "lbfgs",
                "random_state": 0,
                "hyperparameter_selection_on_external_holdout": False,
                "calibration": None,
                "threshold": None,
            },
            "scientific_evaluation": {
                "artifact": "homeolog_topology_nested_evaluation_v0.2",
                "warning": "this final pooled fit has no in-sample performance claim",
            },
            "label_counts_by_species": {
                record["name"]: {
                    "rows": record["rows"],
                    "positive_exact_cds": record["positive_exact_cds"],
                }
                for record in dataset_records
            },
        }
        _write_json(partial_path / "training_report.json", report)
        environment = {
            "python": sys.version,
            "platform": platform.platform(),
            "ploidypatch": ploidypatch_version,
            "numpy": np.__version__,
            "scikit_learn": sklearn.__version__,
            "executable": sys.executable,
            "script": str(Path(__file__).resolve()),
            "script_sha256": _sha256(Path(__file__).resolve()),
        }
        _write_json(partial_path / "environment.json", environment)
        with (partial_path / "input_manifest.tsv").open(
            "x", encoding="utf-8", newline=""
        ) as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=("species", "role", "bytes", "sha256", "path"),
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            for record in dataset_records:
                for role in ("labeled_features", "topology_features"):
                    artifact = record[role]
                    writer.writerow(
                        {
                            "species": record["name"],
                            "role": role,
                            "bytes": artifact["bytes"],
                            "sha256": artifact["sha256"],
                            "path": artifact["path"],
                        }
                    )
        with (partial_path / "run_contract.tsv").open(
            "x", encoding="utf-8", newline=""
        ) as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(("field", "value"))
            writer.writerow(
                (
                    "code_commit",
                    os.environ.get(
                        "PLOIDYPATCH_CODE_COMMIT", "unavailable_server_mirror"
                    ),
                )
            )
            writer.writerow(("development_species", ",".join(names)))
            writer.writerow(("primary_estimator", "topology"))
            writer.writerow(("fixed_C", FIXED_C))
            writer.writerow(("calibrated_probability", "false"))
            writer.writerow(("portable_threshold", "none"))
            writer.writerow(("automatic_approval", "false"))
            writer.writerow(("external_v2_truth_access", "false"))
        artifact_names = sorted(path.name for path in partial_path.iterdir())
        with (partial_path / "SHA256SUMS").open(
            "x", encoding="utf-8", newline=""
        ) as handle:
            for name in artifact_names:
                handle.write(f"{_sha256(partial_path / name)}  {name}\n")
        os.replace(partial_path, output_path)
    except BaseException:
        raise


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Freeze the pooled plant homeolog copy-ranking model v0.2"
    )
    parser.add_argument("--dataset", action="append", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()
    train(args.dataset, args.output_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
