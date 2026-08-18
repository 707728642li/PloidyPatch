#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import re
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any, Iterator, TextIO


SCHEMA_VERSION = "ploidypatch.external_input_preflight.v0.4"
SOURCE_FIELDS = [
    "role",
    "species_id",
    "release",
    "artifact",
    "primary_seqid_regex",
    "source_path",
]
ARTIFACTS = {"genome", "gff3", "protein"}
ROLES = {"target", "candidate_reference", "evaluator_reference"}
GFF_ALIAS_KEYS = {
    "ID",
    "Name",
    "gene_id",
    "locus_tag",
    "protein_id",
    "orig_protein_id",
    "transcript_id",
}
PROTEIN_ALIAS_KEYS = {"ID", "locus", "gene", "transcript"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def fasta_records(path: Path) -> Iterator[tuple[str, str, set[str]]]:
    header = ""
    sequence: list[str] = []
    aliases: set[str] = set()
    with open_text(path) as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\r\n")
            if line.startswith(">"):
                if header:
                    yield header, "".join(sequence), aliases
                header = line[1:]
                if not header:
                    raise ValueError(f"Empty FASTA header at line {line_number}: {path}")
                tokens = header.split()
                aliases = {tokens[0]}
                for token in tokens[1:]:
                    key, separator, value = token.partition("=")
                    if separator and key in PROTEIN_ALIAS_KEYS and value:
                        aliases.add(value)
                    key, separator, value = token.partition(":")
                    if separator and key in PROTEIN_ALIAS_KEYS and value:
                        aliases.add(value)
                sequence = []
            elif not header:
                if line:
                    raise ValueError(f"FASTA sequence precedes header: {path}")
            else:
                sequence.append(line.strip())
    if header:
        yield header, "".join(sequence), aliases


def parse_attributes(value: str) -> dict[str, str]:
    attributes: dict[str, str] = {}
    for token in value.split(";"):
        key, separator, item = token.partition("=")
        if separator and key and item:
            attributes[key] = item
    return attributes


def inspect_genome(path: Path, pattern: str) -> dict[str, Any]:
    primary = re.compile(pattern) if pattern else None
    records = 0
    bases = 0
    primary_records = 0
    primary_bases = 0
    seqids: list[str] = []
    for header, sequence, _ in fasta_records(path):
        seqid = header.split()[0]
        if not sequence:
            raise ValueError(f"Empty genome FASTA sequence: {seqid}")
        records += 1
        bases += len(sequence)
        seqids.append(seqid)
        if primary is not None and primary.fullmatch(seqid):
            primary_records += 1
            primary_bases += len(sequence)
    if not records or (primary is not None and not primary_records):
        raise ValueError(f"Genome or primary-sequence selection is empty: {path}")
    return {
        "records": records,
        "bases": bases,
        "primary_records": primary_records,
        "primary_bases": primary_bases,
        "primary_bases_fraction": primary_bases / bases if bases else None,
        "seqid_sha256": hashlib.sha256(
            "".join(f"{seqid}\n" for seqid in seqids).encode("utf-8")
        ).hexdigest(),
    }


def inspect_gff(path: Path, pattern: str) -> tuple[dict[str, Any], set[str]]:
    primary = re.compile(pattern) if pattern else None
    feature_counts: Counter[str] = Counter()
    primary_feature_counts: Counter[str] = Counter()
    seqids: set[str] = set()
    aliases: set[str] = set()
    rows = 0
    with open_text(path) as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            if not raw_line or raw_line.startswith("#"):
                continue
            fields = raw_line.rstrip("\r\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"Malformed GFF line {line_number}: {path}")
            try:
                start, end = int(fields[3]), int(fields[4])
            except ValueError as exc:
                raise ValueError(f"Invalid GFF coordinate at line {line_number}") from exc
            if start < 1 or end < start:
                raise ValueError(f"Invalid GFF interval at line {line_number}")
            rows += 1
            seqids.add(fields[0])
            feature_counts[fields[2]] += 1
            if primary is not None and primary.fullmatch(fields[0]):
                primary_feature_counts[fields[2]] += 1
            attributes = parse_attributes(fields[8])
            for key in GFF_ALIAS_KEYS:
                if attributes.get(key):
                    aliases.update(
                        item for item in attributes[key].split(",") if item
                    )
    if not rows or not feature_counts.get("gene") or not feature_counts.get("CDS"):
        raise ValueError(f"GFF lacks rows, genes or CDS: {path}")
    return (
        {
            "rows": rows,
            "seqids": len(seqids),
            "feature_counts": dict(sorted(feature_counts.items())),
            "primary_feature_counts": dict(sorted(primary_feature_counts.items())),
            "alias_values": len(aliases),
            "seqid_sha256": hashlib.sha256(
                "".join(f"{seqid}\n" for seqid in sorted(seqids)).encode("utf-8")
            ).hexdigest(),
        },
        aliases,
    )


def inspect_proteins(path: Path, gff_aliases: set[str]) -> dict[str, Any]:
    records = 0
    residues = 0
    exact_alias_records = 0
    identifiers: set[str] = set()
    for header, sequence, aliases in fasta_records(path):
        identifier = header.split()[0]
        if identifier in identifiers or not sequence:
            raise ValueError(f"Duplicate or empty protein record: {identifier}")
        identifiers.add(identifier)
        records += 1
        residues += len(sequence)
        exact_alias_records += int(bool(aliases & gff_aliases))
    if not records:
        raise ValueError(f"Protein FASTA is empty: {path}")
    return {
        "records": records,
        "residues": residues,
        "records_with_any_exact_gff_alias": exact_alias_records,
        "exact_alias_fraction": exact_alias_records / records,
        "identifier_sha256": hashlib.sha256(
            "".join(f"{identifier}\n" for identifier in sorted(identifiers)).encode(
                "utf-8"
            )
        ).hexdigest(),
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Hash and syntax-check external inputs without enumerating WGD pairs"
    )
    parser.add_argument("--sources", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()
    output = Path(args.output_dir)
    partial = Path(str(output) + ".partial")
    if output.exists() or partial.exists():
        raise FileExistsError("Refusing to overwrite external input preflight")
    source_path = Path(args.sources)
    with source_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != SOURCE_FIELDS:
            raise ValueError("External source table fields differ from frozen schema")
        rows = list(reader)
    if not rows:
        raise ValueError("External source table is empty")
    seen = set()
    by_species: dict[str, dict[str, dict[str, str]]] = defaultdict(dict)
    for row in rows:
        key = (row["species_id"], row["artifact"])
        if (
            row["role"] not in ROLES
            or row["artifact"] not in ARTIFACTS
            or key in seen
        ):
            raise ValueError(f"Invalid or duplicate external source row: {key}")
        seen.add(key)
        by_species[row["species_id"]][row["artifact"]] = row
    if any(set(artifacts) != ARTIFACTS for artifacts in by_species.values()):
        raise ValueError("Every species requires genome, GFF3 and protein artifacts")
    if sum(row["role"] == "target" for row in rows if row["artifact"] == "genome") != 1:
        raise ValueError("External preflight requires exactly one target species")

    partial.mkdir(parents=True)
    def process_species(
        item: tuple[str, dict[str, dict[str, str]]]
    ) -> tuple[str, list[dict[str, Any]], dict[str, Any]]:
        species, artifacts = item
        species_manifest_rows = []
        for row in artifacts.values():
            path = Path(row["source_path"])
            if not path.is_file() or path.stat().st_size == 0:
                raise FileNotFoundError(f"Missing external source: {path}")
            species_manifest_rows.append(
                {
                    "role": row["role"],
                    "species_id": species,
                    "release": row["release"],
                    "artifact": row["artifact"],
                    "bytes": path.stat().st_size,
                    "sha256": sha256(path),
                    "source_path": str(path),
                }
            )
        genome_row = artifacts["genome"]
        gff_row = artifacts["gff3"]
        protein_row = artifacts["protein"]
        gff, aliases = inspect_gff(
            Path(gff_row["source_path"]), gff_row["primary_seqid_regex"]
        )
        species_report = {
            "role": genome_row["role"],
            "release": genome_row["release"],
            "primary_seqid_regex": genome_row["primary_seqid_regex"],
            "genome": inspect_genome(
                Path(genome_row["source_path"]), genome_row["primary_seqid_regex"]
            ),
            "gff3": gff,
            "protein": inspect_proteins(Path(protein_row["source_path"]), aliases),
        }
        return species, species_manifest_rows, species_report

    manifest_rows = []
    species_reports: dict[str, Any] = {}
    items = sorted(by_species.items())
    with ThreadPoolExecutor(max_workers=min(5, len(items))) as executor:
        results = list(executor.map(process_species, items))
    for species, species_manifest_rows, species_report in results:
        manifest_rows.extend(species_manifest_rows)
        species_reports[species] = species_report
    manifest_path = partial / "input_manifest.tsv"
    with manifest_path.open("x", encoding="utf-8", newline="") as handle:
        fields = [
            "role",
            "species_id",
            "release",
            "artifact",
            "bytes",
            "sha256",
            "source_path",
        ]
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(sorted(manifest_rows, key=lambda row: (row["role"], row["species_id"], row["artifact"])))
    report = {
        "schema_version": SCHEMA_VERSION,
        "scope": "checksum_syntax_counts_and_exact_alias_preflight_only",
        "wgd_pairs_enumerated": False,
        "candidate_counts_computed": False,
        "truth_labels_accessed": False,
        "selection_by_pair_yield_or_performance": False,
        "source_table": {
            "path": str(source_path),
            "sha256": sha256(source_path),
        },
        "species": species_reports,
    }
    with (partial / "metadata.json").open("x", encoding="utf-8", newline="") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
        handle.write("\n")
    checksum_path = partial / "SHA256SUMS"
    with checksum_path.open("x", encoding="utf-8", newline="") as handle:
        for path in sorted(partial.iterdir()):
            if path != checksum_path:
                handle.write(f"{sha256(path)}  {path.name}\n")
    os.replace(partial, output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
