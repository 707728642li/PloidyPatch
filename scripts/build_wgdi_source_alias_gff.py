#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from pathlib import Path
from typing import Any

from ploidypatch.gff import parse_attributes


SCHEMA_VERSION = "ploidypatch.wgdi_source_alias_gff.v1"
SAFE_ALIAS = re.compile(r"^[A-Za-z0-9_.:-]+$")
ALIAS_ATTRIBUTES = ("ID", "gene_id", "locus_tag", "Name")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def build_alias_gff(
    *, source_gff: Path, representatives_tsv: Path, output_gff: Path
) -> dict[str, Any]:
    manifest_path = Path(str(output_gff) + ".manifest.json")
    for required in (source_gff, representatives_tsv):
        if not required.is_file() or required.stat().st_size == 0:
            raise FileNotFoundError(f"Missing WGDI alias input: {required}")
    if output_gff.exists() or manifest_path.exists():
        raise FileExistsError("Refusing to overwrite WGDI source alias artifacts")

    with representatives_tsv.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or "gene_id" not in reader.fieldnames:
            raise ValueError("Representative table lacks gene_id")
        representative_rows = list(reader)
    aliases = [row["gene_id"] for row in representative_rows]
    alias_set = set(aliases)
    if len(alias_set) != len(aliases) or "" in alias_set:
        raise ValueError("Representative gene IDs are empty or duplicated")
    unsafe = sorted(alias for alias in alias_set if not SAFE_ALIAS.fullmatch(alias))
    if unsafe:
        raise ValueError(f"Unsafe representative gene alias: {unsafe[:5]}")

    rows_by_alias: dict[str, str] = {}
    feature_by_alias: dict[str, str] = {}
    source_gene_records = 0
    with source_gff.open("r", encoding="utf-8", newline="") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            if not raw_line.strip() or raw_line.startswith("#"):
                continue
            fields = raw_line.rstrip("\r\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"Malformed source GFF3 line {line_number}")
            if fields[2] != "gene":
                continue
            source_gene_records += 1
            attributes, malformed = parse_attributes(fields[8])
            if malformed:
                raise ValueError(
                    f"Malformed source GFF3 attributes at line {line_number}"
                )
            matched_aliases = {
                attributes.get(key, "") for key in ALIAS_ATTRIBUTES
            } & alias_set
            if not matched_aliases:
                continue
            if len(matched_aliases) != 1:
                raise ValueError(
                    f"Source gene matches multiple WGDI aliases at line {line_number}"
                )
            alias = next(iter(matched_aliases))
            feature_id = attributes.get("ID", "")
            if not feature_id:
                raise ValueError(f"Mapped gene lacks ID at line {line_number}")
            if alias in rows_by_alias:
                raise ValueError(f"WGDI alias maps to multiple genes: {alias}")
            existing_gene_id = attributes.get("gene_id", "")
            if existing_gene_id and existing_gene_id != alias:
                raise ValueError(
                    f"Conflicting pre-existing gene_id for alias {alias}"
                )
            if not existing_gene_id:
                fields[8] += f";gene_id={alias}"
            rows_by_alias[alias] = "\t".join(fields) + "\n"
            feature_by_alias[alias] = feature_id

    missing = sorted(alias_set - set(rows_by_alias))
    if missing:
        raise ValueError(
            "Representative IDs do not map uniquely through exact gene attributes: "
            + ", ".join(missing[:10])
        )
    output_gff.parent.mkdir(parents=True, exist_ok=True)
    with output_gff.open("x", encoding="utf-8", newline="") as handle:
        handle.write("##gff-version 3\n")
        for alias in aliases:
            handle.write(rows_by_alias[alias])
    manifest: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "role": "evaluator_only_identifier_alias_map",
        "truth_access": False,
        "coordinate_change": False,
        "mapping_policy": (
            "representative_gene_id_equals_unique_exact_"
            "source_gene_ID_gene_id_locus_tag_or_Name"
        ),
        "alias_attributes": list(ALIAS_ATTRIBUTES),
        "alias_attribute_added": "gene_id",
        "inputs": {
            "source_gff": {
                "bytes": source_gff.stat().st_size,
                "sha256": sha256(source_gff),
            },
            "representatives_tsv": {
                "bytes": representatives_tsv.stat().st_size,
                "sha256": sha256(representatives_tsv),
            },
        },
        "counts": {
            "source_gene_records": source_gene_records,
            "representative_gene_ids": len(aliases),
            "mapped_gene_ids": len(rows_by_alias),
            "unique_source_feature_ids": len(set(feature_by_alias.values())),
        },
        "output": {
            "bytes": output_gff.stat().st_size,
            "sha256": sha256(output_gff),
        },
    }
    with manifest_path.open("x", encoding="utf-8", newline="") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Build a gene-only GFF3 alias map for frozen WGDI IDs"
    )
    parser.add_argument("--source-gff", required=True, type=Path)
    parser.add_argument("--representatives", required=True, type=Path)
    parser.add_argument("--output-gff", required=True, type=Path)
    args = parser.parse_args()
    report = build_alias_gff(
        source_gff=args.source_gff,
        representatives_tsv=args.representatives,
        output_gff=args.output_gff,
    )
    print(json.dumps(report, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
