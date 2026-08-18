from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from . import __version__


CATALOG_SUMMARY_SCHEMA_VERSION = "ploidypatch.candidate_catalog_summary.v1"
MISSING_VALUE = "__MISSING__"


def _file_sha256(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _display_value(value: str | None) -> str:
    if value is None or value == "":
        return MISSING_VALUE
    return value


def _parse_crossings(crossings: Sequence[str]) -> tuple[tuple[str, str], ...]:
    parsed: list[tuple[str, str]] = []
    seen: set[tuple[str, str]] = set()
    for value in crossings:
        parts = tuple(part.strip() for part in value.split(","))
        if len(parts) != 2 or not all(parts):
            raise ValueError(
                f"Invalid cross specification {value!r}; expected COLUMN_A,COLUMN_B"
            )
        pair = (parts[0], parts[1])
        if pair in seen:
            raise ValueError(f"Duplicate cross specification: {value}")
        seen.add(pair)
        parsed.append(pair)
    return tuple(parsed)


def summarize_candidate_catalog(
    catalog_path: str | Path,
    columns: Sequence[str],
    crossings: Sequence[str] = (),
) -> dict[str, Any]:
    """Count explicit one-way and two-way strata in a candidate catalog."""

    requested_columns = tuple(columns)
    if not requested_columns:
        raise ValueError("At least one summary column is required")
    if len(set(requested_columns)) != len(requested_columns):
        raise ValueError("Summary columns must be unique")
    parsed_crossings = _parse_crossings(crossings)
    required_columns = {
        "candidate_id",
        "gene_id",
        *requested_columns,
        *(column for pair in parsed_crossings for column in pair),
    }

    one_way = {column: Counter() for column in requested_columns}
    joint = {pair: Counter() for pair in parsed_crossings}
    candidate_ids: set[str] = set()
    gene_ids: set[str] = set()
    row_count = 0
    path = Path(catalog_path)
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("Candidate catalog has no header")
        missing_columns = required_columns - set(reader.fieldnames)
        if missing_columns:
            raise ValueError(
                "Candidate catalog is missing required column(s): "
                + ", ".join(sorted(missing_columns))
            )
        for line_number, row in enumerate(reader, start=2):
            candidate_id = row["candidate_id"]
            gene_id = row["gene_id"]
            if not candidate_id or not gene_id:
                raise ValueError(
                    f"Missing candidate_id or gene_id at catalog line {line_number}"
                )
            if candidate_id in candidate_ids:
                raise ValueError(f"Duplicate candidate_id: {candidate_id}")
            if gene_id in gene_ids:
                raise ValueError(f"Duplicate gene_id: {gene_id}")
            candidate_ids.add(candidate_id)
            gene_ids.add(gene_id)
            row_count += 1
            for column, counter in one_way.items():
                counter[_display_value(row[column])] += 1
            for pair, counter in joint.items():
                counter[
                    (_display_value(row[pair[0]]), _display_value(row[pair[1]]))
                ] += 1

    report: dict[str, Any] = {
        "schema_version": CATALOG_SUMMARY_SCHEMA_VERSION,
        "generator": {"name": "PloidyPatch", "version": __version__},
        "catalog": {
            "file_name": path.name,
            "sha256": _file_sha256(path),
            "rows": row_count,
        },
        "missing_value_label": MISSING_VALUE,
        "one_way_counts": {
            column: dict(sorted(counter.items()))
            for column, counter in one_way.items()
        },
        "joint_counts": {
            f"{pair[0]}__x__{pair[1]}": [
                {pair[0]: values[0], pair[1]: values[1], "count": count}
                for values, count in sorted(counter.items())
            ]
            for pair, counter in joint.items()
        },
    }
    return report


def write_candidate_catalog_summary(
    catalog_path: str | Path,
    output_path: str | Path,
    columns: Sequence[str],
    crossings: Sequence[str] = (),
) -> dict[str, Any]:
    output = Path(output_path)
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite catalog summary: {output}")
    report = summarize_candidate_catalog(catalog_path, columns, crossings)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x", encoding="utf-8", newline="") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return report
