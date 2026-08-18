from __future__ import annotations

import gzip
from collections.abc import Iterator
from pathlib import Path
from typing import TextIO


def open_text(path: str | Path) -> TextIO:
    """Open plain-text or gzip-compressed input with consistent decoding."""

    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("rt", encoding="utf-8", errors="replace")


def normalize_feature_id(value: str) -> str:
    """Remove common provider namespace prefixes without changing the ID body."""

    prefixes = {
        "cds",
        "exon",
        "gene",
        "mrna",
        "protein",
        "transcript",
    }
    if ":" in value:
        prefix, body = value.split(":", 1)
        if prefix.lower() in prefixes and body:
            return body
    return value


def parse_fasta_header_fields(header: str) -> dict[str, str]:
    """Extract whitespace-delimited ``key:value`` fields from a FASTA header.

    Ensembl sequence headers commonly use a protein identifier as their first
    token and carry the corresponding transcript identifier in a later
    ``transcript:<ID>`` field.  Keeping those namespaces separate prevents a
    provider naming convention from being mistaken for a missing sequence.
    The first occurrence of a field wins so malformed repeated metadata cannot
    silently change a relationship.
    """

    fields: dict[str, str] = {}
    for token in header.split()[1:]:
        key, separator, value = token.partition(":")
        if not separator or not key or not value:
            continue
        fields.setdefault(key.lower(), value)
    return fields


def fasta_relation_id(
    primary_id: str,
    header: str,
    relation: str = "transcript",
) -> tuple[str, str]:
    """Resolve a FASTA record to a related feature ID and report provenance.

    Explicit header metadata is preferred.  The normalized primary sequence
    identifier is used only as a documented fallback for providers whose FASTA
    headers omit the relationship field.
    """

    fields = parse_fasta_header_fields(header)
    aliases = (relation.lower(), f"{relation.lower()}_id")
    for alias in aliases:
        value = fields.get(alias)
        if value:
            return normalize_feature_id(value), f"header:{alias}"
    return normalize_feature_id(primary_id), "primary_fallback"


def iter_fasta(path: str | Path) -> Iterator[tuple[str, str, str]]:
    """Yield ``(normalized_primary_id, header, sequence)`` FASTA records."""

    header: str | None = None
    chunks: list[str] = []
    with open_text(path) as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    token = header.split(maxsplit=1)[0]
                    yield normalize_feature_id(token), header, "".join(chunks)
                header = line[1:].strip()
                chunks = []
            else:
                if header is None:
                    raise ValueError(f"Sequence encountered before FASTA header in {path}")
                chunks.append(line)
    if header is not None:
        token = header.split(maxsplit=1)[0]
        yield normalize_feature_id(token), header, "".join(chunks)


def read_fai(path: str | Path) -> dict[str, int]:
    """Read sequence lengths from a samtools FASTA index."""

    lengths: dict[str, int] = {}
    with open_text(path) as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 2:
                raise ValueError(f"Malformed FAI line {line_number} in {path}")
            seqid = fields[0]
            if seqid in lengths:
                raise ValueError(f"Duplicate FAI sequence ID {seqid!r} in {path}")
            try:
                length = int(fields[1])
            except ValueError as exc:
                raise ValueError(
                    f"Non-integer FAI length on line {line_number} in {path}"
                ) from exc
            if length < 0:
                raise ValueError(f"Negative FAI length on line {line_number} in {path}")
            lengths[seqid] = length
    return lengths
