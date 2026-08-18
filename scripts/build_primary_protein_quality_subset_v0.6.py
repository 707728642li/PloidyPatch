#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import os
from collections import Counter
from pathlib import Path

from ploidypatch.audit import PROTEIN_ALPHABET
from ploidypatch.gff import parse_attributes
from ploidypatch.io import fasta_relation_id, iter_fasta, normalize_feature_id, open_text


SCHEMA_VERSION = "ploidypatch.primary_protein_quality_subset.v0.6"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(8 * 1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def canonical_feature_id(value: str) -> str:
    """Strip only recognized provider namespaces until the ID is stable."""

    previous = value
    for _ in range(8):
        normalized = normalize_feature_id(previous)
        if normalized == previous:
            return normalized
        previous = normalized
    raise ValueError(f"Identifier has an excessive namespace depth: {value}")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Build an exact transcript-linked primary protein subset for quality audit"
    )
    parser.add_argument("--primary-gff", required=True, type=Path)
    parser.add_argument("--protein-fasta", required=True, type=Path)
    parser.add_argument("--output-fasta", required=True, type=Path)
    args = parser.parse_args()
    manifest_path = Path(str(args.output_fasta) + ".manifest.json")
    if args.output_fasta.exists() or manifest_path.exists():
        raise FileExistsError("Refusing to overwrite primary protein subset")
    for path in (args.primary_gff, args.protein_fasta):
        if not path.is_file() or path.is_symlink() or path.stat().st_size == 0:
            raise ValueError(f"Missing, empty or symlinked quality input: {path}")

    transcripts: set[str] = set()
    expected: set[str] = set()
    invalid_hierarchy_rows = 0
    with open_text(args.primary_gff) as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\r\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) != 9:
                invalid_hierarchy_rows += 1
                continue
            attributes, malformed = parse_attributes(fields[8])
            if malformed:
                invalid_hierarchy_rows += 1
                continue
            if fields[2] in {"mRNA", "transcript"}:
                transcript_id = attributes.get("ID", "")
                if not transcript_id:
                    invalid_hierarchy_rows += 1
                else:
                    transcripts.add(canonical_feature_id(transcript_id))
            elif fields[2] == "CDS":
                parents = {
                    canonical_feature_id(value)
                    for value in attributes.get("Parent", "").split(",")
                    if value
                }
                if len(parents) != 1:
                    invalid_hierarchy_rows += 1
                else:
                    expected.update(parents)
    if not expected:
        raise ValueError("Primary GFF contains no CDS-linked transcript IDs")

    records: dict[str, tuple[str, str]] = {}
    duplicate_relations: set[str] = set()
    invalid_protein_relations: set[str] = set()
    internal_stop_relations: set[str] = set()
    total_proteins = 0
    for record_id, header, sequence in iter_fasta(args.protein_fasta):
        total_proteins += 1
        relation, _ = fasta_relation_id(record_id, header)
        relation = canonical_feature_id(relation)
        if relation in records:
            duplicate_relations.add(relation)
            continue
        upper = sequence.upper()
        if not upper or set(upper) - PROTEIN_ALPHABET:
            invalid_protein_relations.add(relation)
        if "*" in upper.rstrip("*"):
            internal_stop_relations.add(relation)
        records[relation] = (header, sequence)

    mapped = expected & set(records)
    missing = expected - set(records)
    unexpected = set(records) - transcripts
    valid_mapped = mapped - duplicate_relations - invalid_protein_relations
    args.output_fasta.parent.mkdir(parents=True, exist_ok=True)
    with args.output_fasta.open("x", encoding="utf-8", newline="") as handle:
        for relation in sorted(valid_mapped):
            header, sequence = records[relation]
            handle.write(f">{header}\n")
            handle.writelines(
                sequence[index : index + 60] + "\n"
                for index in range(0, len(sequence), 60)
            )
        handle.flush()
        os.fsync(handle.fileno())

    manifest = {
        "schema_version": SCHEMA_VERSION,
        "role": "target_quality_only",
        "truth_access": False,
        "candidate_access": False,
        "mapping_policy": (
            "exact_transcript_relation_ID_after_repeated_recognized_namespace_stripping"
        ),
        "inputs": {
            "primary_gff_sha256": sha256(args.primary_gff),
            "protein_fasta_sha256": sha256(args.protein_fasta),
        },
        "counts": {
            "primary_transcripts": len(transcripts),
            "primary_CDS_transcripts": len(expected),
            "input_proteins": total_proteins,
            "unique_protein_relations": len(records),
            "mapped_primary_CDS_transcripts": len(mapped),
            "valid_mapped_primary_proteins": len(valid_mapped),
            "missing_primary_proteins": len(missing),
            "duplicate_protein_relations": len(duplicate_relations),
            "invalid_protein_relations": len(invalid_protein_relations),
            "internal_stop_relations": len(internal_stop_relations),
            "protein_relations_without_GFF_transcript": len(unexpected),
            "invalid_GFF_hierarchy_rows": invalid_hierarchy_rows,
        },
        "fractions": {
            "exact_unique_GFF_protein_mapping_fraction": len(valid_mapped) / len(expected),
            "valid_protein_sequence_fraction_among_mapped": (
                len(valid_mapped) / len(mapped) if mapped else 0.0
            ),
        },
        "examples": {
            "missing_primary_proteins": sorted(missing)[:10],
            "duplicate_protein_relations": sorted(duplicate_relations)[:10],
            "invalid_protein_relations": sorted(invalid_protein_relations)[:10],
            "internal_stop_relations": sorted(internal_stop_relations)[:10],
        },
        "output": {
            "file_name": args.output_fasta.name,
            "bytes": args.output_fasta.stat().st_size,
            "sha256": sha256(args.output_fasta),
            "records": len(valid_mapped),
        },
    }
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
