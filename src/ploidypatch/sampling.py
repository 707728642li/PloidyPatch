from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

from . import __version__
from .catalog import CATALOG_MANIFEST_SCHEMA_VERSION


SAMPLE_MANIFEST_SCHEMA_VERSION = "ploidypatch.candidate_sample_manifest.v1"
SAMPLE_COLUMNS = (
    "sampling_plan_row",
    "sampling_rank_sha256",
    "sampling_seed",
)


def _file_sha256(path: str | Path, block_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        while block := handle.read(block_size):
            digest.update(block)
    return digest.hexdigest()


def _read_tsv(path: str | Path) -> tuple[tuple[str, ...], list[dict[str, str]]]:
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"TSV has no header: {path}")
        return tuple(reader.fieldnames), list(reader)


def _rank(seed: int, candidate_id: str) -> str:
    return hashlib.sha256(f"{seed}\0{candidate_id}".encode("utf-8")).hexdigest()


def sample_candidate_catalog(
    *,
    catalog_path: str | Path,
    plan_path: str | Path,
    output_tsv_path: str | Path,
    seed: int,
    exclude_tsv_paths: tuple[str | Path, ...] = (),
) -> dict[str, Any]:
    """Take a deterministic evaluator-only sample from declared exact strata."""

    catalog_path = Path(catalog_path)
    plan_path = Path(plan_path)
    output_path = Path(output_tsv_path)
    output_manifest_path = Path(str(output_path) + ".manifest.json")
    collisions = [
        path for path in (output_path, output_manifest_path) if path.exists()
    ]
    if collisions:
        raise FileExistsError(
            "Refusing to overwrite sample artifact(s): "
            + ", ".join(str(path) for path in collisions)
        )

    catalog_manifest_path = Path(str(catalog_path) + ".manifest.json")
    if not catalog_manifest_path.is_file():
        raise FileNotFoundError(
            f"Candidate catalog manifest is required: {catalog_manifest_path}"
        )
    catalog_manifest = json.loads(
        catalog_manifest_path.read_text(encoding="utf-8")
    )
    if catalog_manifest.get("schema_version") != CATALOG_MANIFEST_SCHEMA_VERSION:
        raise ValueError("Unsupported candidate catalog manifest schema")
    catalog_sha256 = _file_sha256(catalog_path)
    if catalog_manifest.get("catalog", {}).get("sha256") != catalog_sha256:
        raise ValueError("Candidate catalog checksum does not match its manifest")
    source_text_sha256 = catalog_manifest.get("source", {}).get("text_sha256")
    if not source_text_sha256:
        raise ValueError("Candidate catalog manifest lacks source text SHA-256")

    catalog_columns, catalog_rows = _read_tsv(catalog_path)
    required_catalog_columns = {"candidate_id", "gene_id"}
    if not required_catalog_columns <= set(catalog_columns):
        raise ValueError("Candidate catalog lacks candidate_id or gene_id")
    if set(SAMPLE_COLUMNS) & set(catalog_columns):
        raise ValueError("Candidate catalog collides with sampling output columns")
    candidate_ids: set[str] = set()
    gene_ids: set[str] = set()
    for row in catalog_rows:
        if not row["candidate_id"] or not row["gene_id"]:
            raise ValueError("Candidate catalog contains an empty candidate_id or gene_id")
        if row["candidate_id"] in candidate_ids:
            raise ValueError(f"Duplicate candidate_id: {row['candidate_id']}")
        if row["gene_id"] in gene_ids:
            raise ValueError(f"Duplicate gene_id: {row['gene_id']}")
        candidate_ids.add(row["candidate_id"])
        gene_ids.add(row["gene_id"])

    excluded_candidate_ids: set[str] = set()
    exclusion_manifest = []
    for raw_exclusion_path in exclude_tsv_paths:
        exclusion_path = Path(raw_exclusion_path)
        exclusion_manifest_path = Path(str(exclusion_path) + ".manifest.json")
        if not exclusion_manifest_path.is_file():
            raise FileNotFoundError(
                f"Exclusion selection manifest is required: {exclusion_manifest_path}"
            )
        prior_manifest = json.loads(
            exclusion_manifest_path.read_text(encoding="utf-8")
        )
        exclusion_sha256 = _file_sha256(exclusion_path)
        if (
            prior_manifest.get("schema_version") != SAMPLE_MANIFEST_SCHEMA_VERSION
            or prior_manifest.get("output", {}).get("sha256") != exclusion_sha256
        ):
            raise ValueError(
                f"Exclusion selection does not match its manifest: {exclusion_path}"
            )
        exclusion_columns, exclusion_rows = _read_tsv(exclusion_path)
        if "candidate_id" not in exclusion_columns:
            raise ValueError(
                f"Exclusion TSV lacks candidate_id: {exclusion_path}"
            )
        file_ids: set[str] = set()
        for line_number, row in enumerate(exclusion_rows, start=2):
            candidate_id = row["candidate_id"]
            if not candidate_id or candidate_id in file_ids:
                raise ValueError(
                    f"Empty or duplicate exclusion candidate_id at "
                    f"{exclusion_path}:{line_number}"
                )
            if candidate_id not in candidate_ids:
                raise ValueError(
                    f"Exclusion candidate is absent from catalog: {candidate_id}"
                )
            file_ids.add(candidate_id)
        excluded_candidate_ids.update(file_ids)
        exclusion_manifest.append(
            {
                "file_name": exclusion_path.name,
                "sha256": exclusion_sha256,
                "manifest_file_name": exclusion_manifest_path.name,
                "manifest_sha256": _file_sha256(exclusion_manifest_path),
                "rows": len(exclusion_rows),
            }
        )

    plan_columns, plan_rows = _read_tsv(plan_path)
    if "sample_count" not in plan_columns or len(plan_columns) < 2:
        raise ValueError(
            "Sampling plan must contain at least one stratum column and sample_count"
        )
    if not plan_rows:
        raise ValueError("Sampling plan contains no strata")
    stratum_columns = tuple(
        column for column in plan_columns if column != "sample_count"
    )
    missing_stratum_columns = set(stratum_columns) - set(catalog_columns)
    if missing_stratum_columns:
        raise ValueError(
            "Sampling-plan column(s) absent from catalog: "
            + ", ".join(sorted(missing_stratum_columns))
        )

    candidates_by_stratum: dict[tuple[str, ...], list[dict[str, str]]] = defaultdict(
        list
    )
    original_candidates_by_stratum: Counter[tuple[str, ...]] = Counter()
    for row in catalog_rows:
        key = tuple(row[column] for column in stratum_columns)
        original_candidates_by_stratum[key] += 1
        if row["candidate_id"] in excluded_candidate_ids:
            continue
        candidates_by_stratum[key].append(row)

    seen_strata: set[tuple[str, ...]] = set()
    selections: list[dict[str, str]] = []
    stratum_reports: list[dict[str, Any]] = []
    for plan_row_number, plan_row in enumerate(plan_rows, start=2):
        key = tuple(plan_row[column] for column in stratum_columns)
        if key in seen_strata:
            raise ValueError(
                f"Duplicate sampling-plan stratum at line {plan_row_number}: {key}"
            )
        seen_strata.add(key)
        try:
            sample_count = int(plan_row["sample_count"])
        except (TypeError, ValueError) as exc:
            raise ValueError(
                f"Invalid sample_count at plan line {plan_row_number}"
            ) from exc
        if sample_count < 1:
            raise ValueError(
                f"sample_count must be positive at plan line {plan_row_number}"
            )
        available = candidates_by_stratum.get(key, [])
        original_available = original_candidates_by_stratum[key]
        if len(available) < sample_count:
            labels = ", ".join(
                f"{column}={value!r}"
                for column, value in zip(stratum_columns, key, strict=True)
            )
            raise ValueError(
                f"Sampling-plan stratum has {len(available)} candidates but requests "
                f"{sample_count}: {labels}"
            )
        ranked = sorted(
            available,
            key=lambda row: (_rank(seed, row["candidate_id"]), row["candidate_id"]),
        )
        selected = ranked[:sample_count]
        for row in selected:
            sampled = dict(row)
            sampled["sampling_plan_row"] = str(plan_row_number)
            sampled["sampling_rank_sha256"] = _rank(seed, row["candidate_id"])
            sampled["sampling_seed"] = str(seed)
            selections.append(sampled)
        stratum_reports.append(
            {
                "plan_row": plan_row_number,
                "stratum": {
                    column: value
                    for column, value in zip(stratum_columns, key, strict=True)
                },
                "available_candidates": len(available),
                "excluded_candidates": original_available - len(available),
                "requested_candidates": sample_count,
                "selected_candidates": len(selected),
            }
        )

    selected_ids = [row["candidate_id"] for row in selections]
    if len(selected_ids) != len(set(selected_ids)):
        raise AssertionError("A candidate was selected by more than one stratum")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=(*catalog_columns, *SAMPLE_COLUMNS),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(selections)

    manifest: dict[str, Any] = {
        "schema_version": SAMPLE_MANIFEST_SCHEMA_VERSION,
        "generator": {"name": "PloidyPatch", "version": __version__},
        "access": "evaluator_only",
        "source": {"text_sha256": source_text_sha256},
        "catalog": {
            "file_name": catalog_path.name,
            "sha256": catalog_sha256,
            "manifest_file_name": catalog_manifest_path.name,
            "manifest_sha256": _file_sha256(catalog_manifest_path),
            "rows": len(catalog_rows),
        },
        "plan": {
            "file_name": plan_path.name,
            "sha256": _file_sha256(plan_path),
            "stratum_columns": list(stratum_columns),
            "strata": len(plan_rows),
            "requested_candidates": sum(
                int(row["sample_count"]) for row in plan_rows
            ),
        },
        "selection": {
            "algorithm": "sha256_candidate_rank_v1",
            "seed": seed,
            "selected_candidates": len(selections),
            "strata": stratum_reports,
        },
        "output": {
            "file_name": output_path.name,
            "sha256": _file_sha256(output_path),
            "rows": len(selections),
        },
    }
    if exclusion_manifest:
        manifest["exclusions"] = {
            "inputs": exclusion_manifest,
            "unique_candidates": len(excluded_candidate_ids),
            "policy": "candidate_id_union_before_stratified_ranking",
        }
    with output_manifest_path.open("x", encoding="utf-8", newline="") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return manifest
