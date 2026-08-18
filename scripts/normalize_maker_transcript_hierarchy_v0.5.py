#!/usr/bin/env python3
"""Normalize two narrowly defined MAKER transcript-hierarchy conventions.

The adapter never changes sequence IDs, coordinates, strands, phases or CDS
rows.  It only makes gene and transcript identifiers distinct and, for a
gene-with-direct-children provider, inserts one deterministic mRNA row whose
span is exactly the source gene span.
"""
from __future__ import annotations

import argparse
from collections import Counter
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
from typing import Iterable, TextIO


SCHEMA_VERSION = "ploidypatch.maker_transcript_hierarchy_compat.v0.5"
MODES = ("shared_gene_transcript_id", "gene_with_direct_children")
_ID = re.compile(r"(?:^|;)ID=([^;]+)(?=;|$)")
_PARENT = re.compile(r"(?:^|;)Parent=([^;]+)(?=;|$)")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(8 * 1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def _open_source_text(path: Path) -> TextIO:
    """Open a provider GFF3 without silently guessing any format but gzip."""
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def _replace_attribute(attributes: str, key: str, value: str) -> str:
    pattern = re.compile(rf"(^|;){re.escape(key)}=[^;]+(?=;|$)")
    replaced, count = pattern.subn(lambda match: f"{match.group(1)}{key}={value}", attributes)
    if count != 1:
        raise ValueError(f"Expected exactly one {key} attribute, observed {count}: {attributes}")
    return replaced


def _append_attribute(attributes: str, key: str, value: str) -> str:
    if re.search(rf"(?:^|;){re.escape(key)}=", attributes):
        raise ValueError(f"Refusing to overwrite existing {key} attribute: {attributes}")
    return f"{attributes};{key}={value}" if attributes else f"{key}={value}"


def _feature_id(attributes: str, *, line_number: int) -> str:
    matches = _ID.findall(attributes)
    if len(matches) != 1 or not matches[0]:
        raise ValueError(f"Feature line {line_number} requires exactly one non-empty ID")
    return matches[0]


def _parent_ids(attributes: str) -> tuple[str, ...]:
    matches = _PARENT.findall(attributes)
    if not matches:
        return ()
    if len(matches) != 1:
        raise ValueError("A feature may contain only one Parent attribute")
    values = tuple(value for value in matches[0].split(",") if value)
    if not values:
        raise ValueError("Parent attribute may not be empty")
    return values


def normalize_hierarchy(
    *, input_gff: str | Path, output_gff: str | Path, mode: str
) -> dict[str, object]:
    source = Path(input_gff)
    output = Path(output_gff)
    manifest_path = Path(str(output) + ".manifest.json")
    working = Path(str(output) + ".working")
    manifest_working = Path(str(manifest_path) + ".working")
    if mode not in MODES:
        raise ValueError(f"Unsupported hierarchy mode: {mode}")
    if not source.is_file() or source.is_symlink() or source.stat().st_size == 0:
        raise ValueError(f"Input GFF3 must be a non-empty regular file: {source}")
    if any(
        path.exists() or path.is_symlink()
        for path in (output, manifest_path, working, manifest_working)
    ):
        raise FileExistsError("Refusing to overwrite hierarchy adapter output")
    output.parent.mkdir(parents=True, exist_ok=True)

    counts: Counter[str] = Counter()
    source_gene_ids: set[str] = set()
    source_transcript_ids: set[str] = set()
    output_gene_ids: set[str] = set()
    output_transcript_ids: set[str] = set()
    source_cds_digest = hashlib.sha256()
    output_cds_digest = hashlib.sha256()
    gene_rows: dict[str, tuple[str, ...]] = {}
    observed_child_parents: set[str] = set()
    try:
        with _open_source_text(source) as incoming, working.open(
            "x", encoding="utf-8", newline=""
        ) as outgoing:
            for line_number, raw in enumerate(incoming, start=1):
                counts["input_lines"] += 1
                if not raw.strip() or raw.startswith("#"):
                    outgoing.write(raw)
                    counts["output_lines"] += 1
                    continue
                stripped = raw.rstrip("\r\n")
                fields = stripped.split("\t")
                if len(fields) != 9:
                    raise ValueError(f"Malformed GFF3 line {line_number}: expected 9 columns")
                try:
                    start, end = int(fields[3]), int(fields[4])
                except ValueError as error:
                    raise ValueError(f"Malformed coordinates on line {line_number}") from error
                if start < 1 or end < start or fields[6] not in {"+", "-", ".", "?"}:
                    raise ValueError(f"Invalid GFF3 interval or strand on line {line_number}")
                feature_type = fields[2]
                if feature_type == "CDS":
                    source_cds_digest.update((stripped + "\n").encode("utf-8"))
                if feature_type == "gene":
                    gene_id = _feature_id(fields[8], line_number=line_number)
                    if gene_id in source_gene_ids:
                        raise ValueError(f"Duplicate gene ID: {gene_id}")
                    source_gene_ids.add(gene_id)
                    normalized_gene_id = f"gene:{gene_id}"
                    if normalized_gene_id in output_gene_ids:
                        raise ValueError(f"Normalized gene ID collision: {normalized_gene_id}")
                    output_gene_ids.add(normalized_gene_id)
                    fields[8] = _replace_attribute(fields[8], "ID", normalized_gene_id)
                    fields[8] = _append_attribute(fields[8], "ploidypatch_original_gene_id", gene_id)
                    gene_rows[gene_id] = tuple(fields)
                    outgoing.write("\t".join(fields) + "\n")
                    counts["gene_ids_namespaced"] += 1
                    counts["output_lines"] += 1
                    if mode == "gene_with_direct_children":
                        transcript = list(fields)
                        transcript[2] = "mRNA"
                        transcript[7] = "."
                        transcript[8] = (
                            f"ID={gene_id};Parent={normalized_gene_id};"
                            "ploidypatch_synthetic_transcript=true"
                        )
                        outgoing.write("\t".join(transcript) + "\n")
                        output_transcript_ids.add(gene_id)
                        counts["synthetic_transcripts"] += 1
                        counts["output_lines"] += 1
                    continue
                if feature_type in {"mRNA", "transcript"}:
                    transcript_id = _feature_id(fields[8], line_number=line_number)
                    if transcript_id in source_transcript_ids:
                        raise ValueError(f"Duplicate transcript ID: {transcript_id}")
                    source_transcript_ids.add(transcript_id)
                    if transcript_id in output_transcript_ids:
                        raise ValueError(f"Output transcript ID collision: {transcript_id}")
                    output_transcript_ids.add(transcript_id)
                    parents = _parent_ids(fields[8])
                    if mode == "shared_gene_transcript_id":
                        if parents != (transcript_id,):
                            raise ValueError(
                                "Shared-ID mode requires each transcript Parent to equal its ID"
                            )
                        fields[8] = _replace_attribute(
                            fields[8], "Parent", f"gene:{transcript_id}"
                        )
                        counts["transcript_parents_namespaced"] += 1
                    else:
                        raise ValueError(
                            "Direct-child mode forbids source mRNA/transcript records"
                        )
                    outgoing.write("\t".join(fields) + "\n")
                    counts["output_lines"] += 1
                    continue

                observed_child_parents.update(_parent_ids(fields[8]))
                outgoing.write(stripped + "\n")
                counts["output_lines"] += 1
                if feature_type == "CDS":
                    output_cds_digest.update((stripped + "\n").encode("utf-8"))
                    counts["cds_rows"] += 1

        if not source_gene_ids:
            raise ValueError("Input GFF3 contains no gene records")
        if mode == "shared_gene_transcript_id":
            if source_gene_ids != source_transcript_ids:
                raise ValueError("Gene and transcript ID universes differ in shared-ID mode")
        elif source_transcript_ids:
            raise ValueError("Direct-child mode unexpectedly observed source transcripts")
        if output_transcript_ids != source_gene_ids:
            raise ValueError("Output transcript IDs do not exactly match source gene IDs")
        unknown_parents = observed_child_parents - output_transcript_ids - output_gene_ids
        if unknown_parents:
            raise ValueError(
                "Non-gene child Parent values are unresolved after normalization: "
                + ", ".join(sorted(unknown_parents)[:10])
            )
        if source_cds_digest.hexdigest() != output_cds_digest.hexdigest():
            raise RuntimeError("CDS rows changed during transcript hierarchy normalization")
        os.replace(working, output)
    except BaseException:
        if working.exists() and not working.is_symlink():
            working.unlink()
        raise

    manifest: dict[str, object] = {
        "schema_version": SCHEMA_VERSION,
        "mode": mode,
        "coordinate_or_cds_changes": False,
        "source": {
            "path": source.name,
            "bytes": source.stat().st_size,
            "sha256": _sha256(source),
        },
        "output": {
            "path": output.name,
            "bytes": output.stat().st_size,
            "sha256": _sha256(output),
        },
        "counts": dict(sorted(counts.items())),
        "source_genes": len(source_gene_ids),
        "output_transcripts": len(output_transcript_ids),
        "cds_rows_sha256": {
            "input": source_cds_digest.hexdigest(),
            "output": output_cds_digest.hexdigest(),
        },
    }
    try:
        with manifest_working.open("x", encoding="utf-8", newline="") as handle:
            json.dump(manifest, handle, indent=2, sort_keys=True)
            handle.write("\n")
        os.replace(manifest_working, manifest_path)
    except BaseException:
        if manifest_working.exists() and not manifest_working.is_symlink():
            manifest_working.unlink()
        if output.exists() and not output.is_symlink():
            output.unlink()
        raise
    return manifest


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-gff", required=True)
    parser.add_argument("--output-gff", required=True)
    parser.add_argument("--mode", required=True, choices=MODES)
    args = parser.parse_args(argv)
    report = normalize_hierarchy(
        input_gff=args.input_gff, output_gff=args.output_gff, mode=args.mode
    )
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
