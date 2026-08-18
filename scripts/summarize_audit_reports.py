#!/usr/bin/env python3
"""Build a compact, deterministic TSV from PloidyPatch audit JSON reports."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any


FIELDS = (
    "dataset_id",
    "grade",
    "assembly_sequences",
    "assembly_bp",
    "genes",
    "transcripts",
    "expected_translation_transcripts",
    "protein_records",
    "protein_header_transcript_mappings",
    "protein_primary_id_fallbacks",
    "transcripts_missing_protein",
    "proteins_without_transcript",
    "transcripts_missing_cds",
    "cds_without_transcript",
    "protein_coding_transcripts_without_gff_cds",
    "translation_length_checks",
    "cds_protein_length_mismatches",
    "cds_not_multiple_of_three",
    "cds_frame_remainder_translation_consistent",
    "cds_frame_remainder_translation_mismatch",
    "nonzero_error_categories",
)


def nested(record: dict[str, Any], *keys: str, default: Any = 0) -> Any:
    value: Any = record
    for key in keys:
        if not isinstance(value, dict) or key not in value:
            return default
        value = value[key]
    return value


def summarize(path: Path) -> dict[str, Any]:
    report = json.loads(path.read_text(encoding="utf-8"))
    errors = nested(report, "quality_gate", "errors", default={})
    sources = nested(report, "protein", "relation_id_source_counts", default={})
    return {
        "dataset_id": path.stem,
        "grade": nested(report, "quality_gate", "grade", default="unknown"),
        "assembly_sequences": nested(report, "assembly", "sequence_count"),
        "assembly_bp": nested(report, "assembly", "total_bp"),
        "genes": nested(report, "gff3", "gene_ids"),
        "transcripts": nested(report, "gff3", "transcript_ids"),
        "expected_translation_transcripts": nested(
            report, "cross_checks", "expected_translation_transcript_ids"
        ),
        "protein_records": nested(report, "protein", "sequence_count"),
        "protein_header_transcript_mappings": sources.get("header:transcript", 0),
        "protein_primary_id_fallbacks": sources.get("primary_fallback", 0),
        "transcripts_missing_protein": nested(
            report, "cross_checks", "transcripts_missing_protein"
        ),
        "proteins_without_transcript": nested(
            report, "cross_checks", "proteins_without_transcript"
        ),
        "transcripts_missing_cds": nested(
            report, "cross_checks", "transcripts_missing_cds"
        ),
        "cds_without_transcript": nested(
            report, "cross_checks", "cds_without_transcript"
        ),
        "protein_coding_transcripts_without_gff_cds": nested(
            report, "gff3", "protein_coding_transcripts_without_cds"
        ),
        "translation_length_checks": nested(
            report, "cross_checks", "translation_length_checks"
        ),
        "cds_protein_length_mismatches": nested(
            report, "cross_checks", "cds_protein_length_mismatches"
        ),
        "cds_not_multiple_of_three": nested(report, "cds", "not_multiple_of_three"),
        "cds_frame_remainder_translation_consistent": nested(
            report, "cross_checks", "cds_frame_remainder_translation_consistent"
        ),
        "cds_frame_remainder_translation_mismatch": nested(
            report, "cross_checks", "cds_frame_remainder_translation_mismatch"
        ),
        "nonzero_error_categories": sum(bool(value) for value in errors.values()),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("report_dir", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    output = args.output or args.report_dir / "audit_summary.tsv"
    rows = [summarize(path) for path in sorted(args.report_dir.glob("*.json"))]
    if not rows:
        parser.error(f"no JSON reports found in {args.report_dir}")
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
