"""Truth-blind reproducibility gates for projected candidate models.

The primitives in this module compare complete adapter decision rows from two
independent executions.  A model is eligible only when the same model ID and
the entire decision row occur in both runs.  No biological label, performance
metric, or species-specific identifier is accepted by this API.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass
import json
import os
from pathlib import Path, PurePosixPath
from typing import Iterable

from .artifact_manifest import sha256_file
from .gff import parse_attributes


FORBIDDEN_DECISION_FIELDS = frozenset(
    {"label", "truth", "is_positive", "target_label", "event_label"}
)


@dataclass(frozen=True)
class DecisionComparison:
    """Exact model-level comparison of two adapter decision tables."""

    fields: tuple[str, ...]
    rows_a: dict[str, tuple[str, ...]]
    rows_b: dict[str, tuple[str, ...]]
    stable_models: frozenset[str]
    stable_accepted_models: frozenset[str]

    @property
    def all_models(self) -> tuple[str, ...]:
        return tuple(sorted(set(self.rows_a) | set(self.rows_b)))


def _read_decisions(path: Path) -> tuple[tuple[str, ...], dict[str, tuple[str, ...]]]:
    if not path.is_file() or path.is_symlink() or path.stat().st_size == 0:
        raise ValueError(f"Missing, empty or symlinked decision table: {path}")
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            raw_fields = next(reader)
        except StopIteration as error:
            raise ValueError(f"Empty decision table: {path}") from error
        fields = tuple(raw_fields)
        if not fields or len(fields) != len(set(fields)) or "model_id" not in fields:
            raise ValueError(f"Malformed decision header: {path}")
        if {field.casefold() for field in fields} & FORBIDDEN_DECISION_FIELDS:
            raise ValueError(f"Truth-bearing decision field is forbidden: {path}")
        model_index = fields.index("model_id")
        rows: dict[str, tuple[str, ...]] = {}
        for line_number, raw in enumerate(reader, start=2):
            row = tuple(raw)
            if len(row) != len(fields):
                raise ValueError(f"Malformed decision row {line_number}: {path}")
            model_id = row[model_index]
            if not model_id or model_id in rows:
                raise ValueError(f"Duplicate or empty model ID at row {line_number}: {path}")
            rows[model_id] = row
    if not rows:
        raise ValueError(f"Decision table has no models: {path}")
    return fields, rows


def compare_decision_tables(path_a: str | Path, path_b: str | Path) -> DecisionComparison:
    """Return models whose complete decision rows are byte-semantically identical."""

    source_a = Path(path_a)
    source_b = Path(path_b)
    fields_a, rows_a = _read_decisions(source_a)
    fields_b, rows_b = _read_decisions(source_b)
    if fields_a != fields_b:
        raise ValueError("Independent decision tables have different schemas")
    stable = frozenset(
        model_id
        for model_id in set(rows_a) & set(rows_b)
        if rows_a[model_id] == rows_b[model_id]
    )
    status_index = fields_a.index("status") if "status" in fields_a else None
    if status_index is None:
        raise ValueError("Decision tables lack a status field")
    stable_accepted = frozenset(
        model_id for model_id in stable if rows_a[model_id][status_index] == "accepted"
    )
    return DecisionComparison(
        fields=fields_a,
        rows_a=rows_a,
        rows_b=rows_b,
        stable_models=stable,
        stable_accepted_models=stable_accepted,
    )


def write_comparison_audit(
    *,
    comparison: DecisionComparison,
    decisions_a: str | Path,
    decisions_b: str | Path,
    output_tsv: str | Path,
    output_json: str | Path,
    method: str,
    reference: str,
) -> dict[str, object]:
    """Write a label-free, exact-universe audit for one method/reference arm."""

    destination = Path(output_tsv)
    manifest_path = Path(output_json)
    if any(path.exists() or path.is_symlink() for path in (destination, manifest_path)):
        raise FileExistsError("Refusing to overwrite reproducibility audit")
    destination.parent.mkdir(parents=True, exist_ok=True)
    status_index = comparison.fields.index("status")
    with destination.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=(
                "method",
                "reference",
                "model_id",
                "run_a_present",
                "run_b_present",
                "exact_decision_match",
                "run_a_status",
                "run_b_status",
                "stable_accepted",
            ),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for model_id in comparison.all_models:
            row_a = comparison.rows_a.get(model_id)
            row_b = comparison.rows_b.get(model_id)
            writer.writerow(
                {
                    "method": method,
                    "reference": reference,
                    "model_id": model_id,
                    "run_a_present": int(row_a is not None),
                    "run_b_present": int(row_b is not None),
                    "exact_decision_match": int(model_id in comparison.stable_models),
                    "run_a_status": "" if row_a is None else row_a[status_index],
                    "run_b_status": "" if row_b is None else row_b[status_index],
                    "stable_accepted": int(
                        model_id in comparison.stable_accepted_models
                    ),
                }
            )
    payload: dict[str, object] = {
        "schema_version": "ploidypatch.projection_reproducibility_audit.v1",
        "method": method,
        "reference": reference,
        "truth_access": False,
        "label_access": False,
        "selection_rule": "same_model_id_and_exact_complete_decision_row_in_both_runs",
        "decisions_a_sha256": sha256_file(decisions_a),
        "decisions_b_sha256": sha256_file(decisions_b),
        "decision_fields": list(comparison.fields),
        "models_a": len(comparison.rows_a),
        "models_b": len(comparison.rows_b),
        "models_union": len(comparison.all_models),
        "stable_models": len(comparison.stable_models),
        "stable_accepted_models": len(comparison.stable_accepted_models),
        "unstable_models": len(comparison.all_models)
        - len(comparison.stable_models),
        "audit_tsv_sha256": sha256_file(destination),
    }
    manifest_path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return payload


def _parse_feature(raw: str, *, line_number: int) -> tuple[list[str], dict[str, str]]:
    fields = raw.rstrip("\r\n").split("\t")
    if len(fields) != 9:
        raise ValueError(f"Malformed candidate GFF row {line_number}")
    attributes, malformed = parse_attributes(fields[8])
    if malformed:
        raise ValueError(f"Malformed candidate attributes at row {line_number}")
    return fields, attributes


def filter_candidate_gff_by_upstream_models(
    *,
    source: str | Path,
    allowed_models: Iterable[str],
    output: str | Path,
    model_attributes: tuple[str, ...] = ("upstream_model", "miniprot_model"),
) -> dict[str, int]:
    """Keep complete gene blocks whose exact upstream model passed both runs."""

    source_path = Path(source)
    destination = Path(output)
    if not source_path.is_file() or source_path.is_symlink():
        raise ValueError(f"Missing or symlinked candidate GFF: {source_path}")
    if destination.exists() or destination.is_symlink():
        raise FileExistsError(f"Refusing to overwrite stable candidate GFF: {destination}")
    allowed = frozenset(allowed_models)
    if not allowed:
        raise ValueError("Stable accepted model universe is empty")
    if not model_attributes or len(model_attributes) != len(set(model_attributes)):
        raise ValueError("Model-lineage attribute names are empty or duplicated")
    destination.parent.mkdir(parents=True, exist_ok=True)
    header: list[str] = []
    blocks: list[tuple[str, list[str]]] = []
    current_model: str | None = None
    current_rows: list[str] = []
    current_ids: set[str] = set()
    current_parents: set[str] = set()

    def finish_block() -> None:
        nonlocal current_model, current_rows, current_ids, current_parents
        if current_model is None:
            return
        if not current_rows or not current_ids or not current_parents <= current_ids:
            raise ValueError(f"Malformed candidate feature hierarchy: {current_model}")
        blocks.append((current_model, current_rows))
        current_model = None
        current_rows = []
        current_ids = set()
        current_parents = set()

    with source_path.open(encoding="utf-8", newline="") as handle:
        for line_number, raw in enumerate(handle, start=1):
            if not raw.strip() or raw.startswith("#"):
                if current_model is not None:
                    raise ValueError("Directive or blank row inside candidate gene block")
                header.append(raw)
                continue
            fields, attributes = _parse_feature(raw, line_number=line_number)
            feature_type = fields[2]
            lineage_values = {
                attributes[name]
                for name in model_attributes
                if attributes.get(name)
            }
            if len(lineage_values) > 1:
                raise ValueError(f"Conflicting candidate lineage at row {line_number}")
            lineage = next(iter(lineage_values), None)
            if feature_type == "gene":
                finish_block()
                model = lineage or ""
                feature_id = attributes.get("ID", "")
                if not model or not feature_id:
                    raise ValueError(f"Candidate gene lacks lineage at row {line_number}")
                current_model = model
            elif current_model is None:
                raise ValueError(f"Candidate child precedes a gene at row {line_number}")
            descendant_model = lineage
            if descendant_model is not None and descendant_model != current_model:
                raise ValueError(f"Mixed upstream models in candidate block: {current_model}")
            feature_id = attributes.get("ID")
            if feature_id:
                if feature_id in current_ids and feature_type not in {"CDS", "exon"}:
                    raise ValueError(f"Duplicate candidate feature ID: {feature_id}")
                current_ids.add(feature_id)
            for parent in attributes.get("Parent", "").split(","):
                if parent:
                    current_parents.add(parent)
            current_rows.append(raw)
    finish_block()
    models = [model for model, _rows in blocks]
    if len(models) != len(set(models)):
        raise ValueError("Candidate GFF repeats an upstream model")
    observed = frozenset(models)
    missing = allowed - observed
    if missing:
        raise ValueError(f"Stable accepted models are absent from candidate GFF: {len(missing)}")
    kept = [(model, rows) for model, rows in blocks if model in allowed]
    with destination.open("x", encoding="utf-8", newline="") as handle:
        if header:
            handle.writelines(header)
        else:
            handle.write("##gff-version 3\n")
        for _model, rows in kept:
            handle.writelines(rows)
    return {
        "source_models": len(blocks),
        "allowed_models": len(allowed),
        "kept_models": len(kept),
        "dropped_models": len(blocks) - len(kept),
    }


def verify_tree_manifest(*, root: str | Path, manifest: str | Path) -> dict[str, int]:
    """Verify an exact regular-file universe against a three-column TSV."""

    root_path = Path(root)
    manifest_path = Path(manifest)
    if not root_path.is_dir() or root_path.is_symlink():
        raise ValueError(f"Missing or symlinked reproducibility root: {root_path}")
    if not manifest_path.is_file() or manifest_path.is_symlink():
        raise ValueError(f"Missing or symlinked reproducibility manifest: {manifest_path}")
    with manifest_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != ("relative_path", "bytes", "sha256"):
            raise ValueError("Reproducibility tree manifest header differs")
        rows = list(reader)
    expected: dict[str, tuple[int, str]] = {}
    for row in rows:
        relative = row["relative_path"]
        pure = PurePosixPath(relative)
        if (
            not relative
            or pure.is_absolute()
            or "\\" in relative
            or any(part in {"", ".", ".."} for part in pure.parts)
            or relative in expected
        ):
            raise ValueError(f"Unsafe reproducibility manifest path: {relative}")
        try:
            size = int(row["bytes"])
        except ValueError as error:
            raise ValueError(f"Malformed reproducibility byte count: {relative}") from error
        digest = row["sha256"]
        if size < 0 or len(digest) != 64 or any(c not in "0123456789abcdef" for c in digest):
            raise ValueError(f"Malformed reproducibility binding: {relative}")
        expected[relative] = (size, digest)
    observed_paths: set[str] = set()
    byte_count = 0
    for path in sorted(root_path.rglob("*")):
        if path.is_symlink():
            raise ValueError(f"Symlink in reproducibility tree: {path}")
        if not path.is_file():
            continue
        relative = path.relative_to(root_path).as_posix()
        observed_paths.add(relative)
        expected_value = expected.get(relative)
        if expected_value is None:
            raise ValueError(f"Unexpected reproducibility file: {relative}")
        size = path.stat().st_size
        if expected_value != (size, sha256_file(path)):
            raise ValueError(f"Reproducibility file differs: {relative}")
        byte_count += size
    if observed_paths != set(expected):
        raise ValueError(
            f"Missing reproducibility files: {sorted(set(expected) - observed_paths)}"
        )
    return {"files": len(observed_paths), "bytes": byte_count}
