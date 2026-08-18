#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from ploidypatch.io import fasta_relation_id, iter_fasta, normalize_feature_id


SCHEMA_VERSION = "ploidypatch.primary_translation_quality.v0.6"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(8 * 1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def canonical_id(value: str) -> str:
    previous = value
    for _ in range(8):
        current = normalize_feature_id(previous)
        if current == previous:
            return current
        previous = current
    raise ValueError(f"Excessive feature namespace depth: {value}")


def read_proteins(path: Path, use_header_relation: bool) -> dict[str, str]:
    records: dict[str, str] = {}
    for primary, header, sequence in iter_fasta(path):
        relation = fasta_relation_id(primary, header)[0] if use_header_relation else primary
        relation = canonical_id(relation)
        if not relation or relation in records:
            raise ValueError(f"Empty or duplicate protein relation in {path}: {relation}")
        records[relation] = sequence.upper().rstrip("*")
    if not records:
        raise ValueError(f"Protein FASTA is empty: {path}")
    return records


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Compare provider proteins with genome/GFF-derived translations"
    )
    parser.add_argument("--provider-protein", required=True, type=Path)
    parser.add_argument("--translated-protein", required=True, type=Path)
    parser.add_argument("--output-json", required=True, type=Path)
    args = parser.parse_args()
    if args.output_json.exists():
        raise FileExistsError("Refusing to overwrite translation-quality report")
    for path in (args.provider_protein, args.translated_protein):
        if not path.is_file() or path.is_symlink() or path.stat().st_size == 0:
            raise ValueError(f"Missing, empty or symlinked translation input: {path}")

    provider = read_proteins(args.provider_protein, True)
    translated = read_proteins(args.translated_protein, False)
    shared = set(provider) & set(translated)
    exact = {identifier for identifier in shared if provider[identifier] == translated[identifier]}
    length_consistent = {
        identifier
        for identifier in shared
        if len(provider[identifier]) == len(translated[identifier])
    }
    missing_translation = set(provider) - set(translated)
    unexpected_translation = set(translated) - set(provider)
    mismatched = shared - exact
    report = {
        "schema_version": SCHEMA_VERSION,
        "role": "target_quality_only",
        "truth_access": False,
        "candidate_access": False,
        "comparison_policy": (
            "exact_amino_acid_identity_after_terminal_stop_removal_and_"
            "repeated_recognized_namespace_stripping"
        ),
        "inputs": {
            "provider_protein_sha256": sha256(args.provider_protein),
            "translated_protein_sha256": sha256(args.translated_protein),
        },
        "counts": {
            "provider_proteins": len(provider),
            "translated_proteins": len(translated),
            "shared_relations": len(shared),
            "exact_translations": len(exact),
            "length_consistent_translations": len(length_consistent),
            "sequence_mismatches": len(mismatched),
            "missing_translations": len(missing_translation),
            "unexpected_translations": len(unexpected_translation),
        },
        "fractions": {
            "valid_representative_translation_fraction": len(exact) / len(provider),
            "length_consistent_translation_fraction": len(length_consistent) / len(provider),
        },
        "examples": {
            "sequence_mismatches": sorted(mismatched)[:10],
            "missing_translations": sorted(missing_translation)[:10],
            "unexpected_translations": sorted(unexpected_translation)[:10],
        },
    }
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
