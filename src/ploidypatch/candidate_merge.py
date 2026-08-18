from __future__ import annotations

import csv
import json
import os
from collections import Counter
from pathlib import Path
from typing import Any, Iterable
from urllib.parse import quote

from . import __version__
from .baseline import SAFE_SOURCE, _file_sha256
from .gff import parse_attributes
from .io import open_text


CANDIDATE_MERGE_SCHEMA_VERSION = "ploidypatch.candidate_gff_merge.v1"
REFERENCE_ATTRIBUTES = ("ID", "Parent", "Derives_from")
ATTRIBUTE_SAFE = "._:-,+*()[]|/"


def _parse_inputs(values: Iterable[str]) -> tuple[tuple[str, Path], ...]:
    inputs: list[tuple[str, Path]] = []
    seen: set[str] = set()
    for value in values:
        label, separator, raw_path = value.partition("=")
        if not separator or not SAFE_SOURCE.fullmatch(label) or not raw_path:
            raise ValueError("Each candidate GFF input must use safe SOURCE=PATH syntax")
        if label in seen:
            raise ValueError(f"Duplicate candidate GFF source label: {label}")
        seen.add(label)
        inputs.append((label, Path(raw_path)))
    if not inputs:
        raise ValueError("At least one reference GFF input is required for namespacing")
    return tuple(inputs)


def _namespace_identifier(source: str, identifier: str) -> str:
    if not identifier:
        raise ValueError("Cannot namespace an empty GFF identifier")
    return f"{source}__{identifier}"


def _serialize_attributes(attributes: dict[str, str]) -> str:
    if not attributes:
        return "."
    return ";".join(
        f"{quote(key, safe=ATTRIBUTE_SAFE)}={quote(value, safe=ATTRIBUTE_SAFE)}"
        for key, value in attributes.items()
    )


def merge_candidate_gffs(
    *,
    candidate_inputs: Iterable[str],
    output_gff_path: str | Path,
    provenance_tsv_path: str | Path,
) -> dict[str, Any]:
    """Namespace and merge same-method reference GFFs without extra votes."""

    inputs = _parse_inputs(candidate_inputs)
    output_path = Path(output_gff_path)
    provenance_path = Path(provenance_tsv_path)
    manifest_path = Path(str(output_path) + ".manifest.json")
    partial_output = Path(str(output_path) + ".partial")
    partial_provenance = Path(str(provenance_path) + ".partial")
    collisions = [
        path
        for path in (
            output_path,
            provenance_path,
            manifest_path,
            partial_output,
            partial_provenance,
        )
        if path.exists()
    ]
    if collisions:
        raise FileExistsError(
            "Refusing to overwrite candidate-merge artifact(s): "
            + ", ".join(str(path) for path in collisions)
        )
    for _, path in inputs:
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f"Missing or empty candidate GFF input: {path}")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    provenance_path.parent.mkdir(parents=True, exist_ok=True)
    source_reports: dict[str, Any] = {}
    all_namespaced_ids: set[str] = set()
    total_features = 0
    total_ids = 0
    total_id_occurrences = 0
    total_repeated_id_lines = 0
    total_duplicate_feature_lines_skipped = 0
    try:
        with partial_output.open("x", encoding="utf-8", newline="") as output, (
            partial_provenance.open("x", encoding="utf-8", newline="")
        ) as provenance:
            output.write("##gff-version 3\n")
            writer = csv.DictWriter(
                provenance,
                fieldnames=(
                    "reference_source",
                    "input_line",
                    "feature_type",
                    "seqid",
                    "original_id",
                    "namespaced_id",
                ),
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            for source, path in inputs:
                source_ids: set[str] = set()
                source_id_signatures: dict[str, tuple[str, str, str, str, str]] = {}
                source_id_segment_lines: dict[
                    tuple[str, int, int, str], str
                ] = {}
                source_parent_refs: set[str] = set()
                source_features = 0
                source_ids_count = 0
                source_id_occurrences = 0
                source_repeated_id_lines = 0
                source_duplicate_feature_lines_skipped = 0
                type_counts: Counter[str] = Counter()
                stopped_at_fasta = False
                with open_text(path) as handle:
                    for line_number, raw_line in enumerate(handle, start=1):
                        stripped = raw_line.rstrip("\r\n")
                        if stripped == "##FASTA":
                            stopped_at_fasta = True
                            break
                        if not stripped or stripped.startswith("#"):
                            continue
                        fields = stripped.split("\t")
                        if len(fields) != 9:
                            raise ValueError(
                                f"Malformed candidate GFF line {line_number}: {path}"
                            )
                        try:
                            start = int(fields[3])
                            end = int(fields[4])
                        except ValueError as exc:
                            raise ValueError(
                                f"Non-integer candidate coordinate at line {line_number}: {path}"
                            ) from exc
                        if start < 1 or end < start:
                            raise ValueError(
                                f"Invalid candidate interval at line {line_number}: {path}"
                            )
                        attributes, malformed = parse_attributes(fields[8])
                        if malformed:
                            raise ValueError(
                                f"Malformed candidate attributes at line {line_number}: {path}"
                            )
                        original_id = attributes.get("ID", "")
                        namespaced_id = ""
                        if original_id:
                            namespaced_id = _namespace_identifier(source, original_id)
                            signature = (
                                fields[2],
                                fields[0],
                                fields[6],
                                attributes.get("Parent", ""),
                                attributes.get("Derives_from", ""),
                            )
                            segment = (original_id, start, end, fields[7])
                            if segment in source_id_segment_lines:
                                if source_id_segment_lines[segment] == stripped:
                                    source_duplicate_feature_lines_skipped += 1
                                    continue
                                raise ValueError(
                                    f"Non-identical duplicate ID segment "
                                    f"{original_id!r} in source {source} at "
                                    f"{fields[0]}:{start}-{end}"
                                )
                            source_id_segment_lines[segment] = stripped
                            if original_id in source_ids:
                                if source_id_signatures[original_id] != signature:
                                    raise ValueError(
                                        f"Conflicting repeated ID {original_id!r} in "
                                        f"source {source}"
                                    )
                                source_repeated_id_lines += 1
                            else:
                                if namespaced_id in all_namespaced_ids:
                                    raise AssertionError("Namespaced GFF ID collision")
                                source_ids.add(original_id)
                                source_id_signatures[original_id] = signature
                                all_namespaced_ids.add(namespaced_id)
                                source_ids_count += 1
                            attributes["ID"] = namespaced_id
                            source_id_occurrences += 1
                            writer.writerow(
                                {
                                    "reference_source": source,
                                    "input_line": line_number,
                                    "feature_type": fields[2],
                                    "seqid": fields[0],
                                    "original_id": original_id,
                                    "namespaced_id": namespaced_id,
                                }
                            )
                        for attribute in ("Parent", "Derives_from"):
                            if attribute not in attributes:
                                continue
                            identifiers = [
                                value for value in attributes[attribute].split(",") if value
                            ]
                            if not identifiers:
                                raise ValueError(
                                    f"Empty {attribute} at line {line_number}: {path}"
                                )
                            if attribute == "Parent":
                                source_parent_refs.update(identifiers)
                            attributes[attribute] = ",".join(
                                _namespace_identifier(source, value)
                                for value in identifiers
                            )
                        fields[8] = _serialize_attributes(attributes)
                        output.write("\t".join(fields) + "\n")
                        source_features += 1
                        type_counts[fields[2]] += 1
                orphan_parents = source_parent_refs - source_ids
                if orphan_parents:
                    raise ValueError(
                        f"Source {source} has Parent references without IDs: "
                        + ", ".join(sorted(orphan_parents)[:10])
                    )
                source_reports[source] = {
                    "file_name": path.name,
                    "bytes": path.stat().st_size,
                    "sha256": _file_sha256(path),
                    "features": source_features,
                    "ids": source_ids_count,
                    "id_occurrences": source_id_occurrences,
                    "repeated_id_lines": source_repeated_id_lines,
                    "duplicate_feature_lines_skipped": (
                        source_duplicate_feature_lines_skipped
                    ),
                    "feature_type_counts": dict(sorted(type_counts.items())),
                    "stopped_at_embedded_fasta": stopped_at_fasta,
                }
                total_features += source_features
                total_ids += source_ids_count
                total_id_occurrences += source_id_occurrences
                total_repeated_id_lines += source_repeated_id_lines
                total_duplicate_feature_lines_skipped += (
                    source_duplicate_feature_lines_skipped
                )
        if total_features == 0 or total_ids == 0:
            raise ValueError("Merged candidate GFF has no features or IDs")
        os.replace(partial_output, output_path)
        os.replace(partial_provenance, provenance_path)
    except BaseException:
        partial_output.unlink(missing_ok=True)
        partial_provenance.unlink(missing_ok=True)
        raise

    manifest: dict[str, Any] = {
        "schema_version": CANDIDATE_MERGE_SCHEMA_VERSION,
        "generator": {"name": "PloidyPatch", "version": __version__},
        "method_family_policy": {
            "reference_count": len(inputs),
            "vote_count_after_merge": 1,
            "rationale": "multiple references from one upstream method are not independent method votes",
            "identifier_namespace": "reference_source__original_id",
            "parent_integrity_required": True,
        },
        "inputs": source_reports,
        "counts": {
            "features": total_features,
            "ids": total_ids,
            "id_occurrences": total_id_occurrences,
            "repeated_id_lines": total_repeated_id_lines,
            "duplicate_feature_lines_skipped": (
                total_duplicate_feature_lines_skipped
            ),
            "reference_sources": len(inputs),
        },
        "outputs": {
            "gff": {
                "file_name": output_path.name,
                "bytes": output_path.stat().st_size,
                "sha256": _file_sha256(output_path),
            },
            "provenance": {
                "file_name": provenance_path.name,
                "bytes": provenance_path.stat().st_size,
                "sha256": _file_sha256(provenance_path),
                "rows": total_id_occurrences,
            },
        },
    }
    with manifest_path.open("x", encoding="utf-8", newline="") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return manifest
