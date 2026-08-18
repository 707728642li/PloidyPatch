from __future__ import annotations

import hashlib
import json
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from statistics import median
from typing import Literal

from .gff import audit_gff
from .io import fasta_relation_id, iter_fasta, read_fai


DNA_ALPHABET = set("ACGTRYSWKMBDHVN")
PROTEIN_ALPHABET = set("ABCDEFGHIKLMNPQRSTVWXYZJUO*")
STOP_CODONS = {"TAA", "TAG", "TGA"}


@dataclass(frozen=True)
class FastaAuditIndex:
    """Per-relation metadata retained for cross-file consistency checks."""

    relation_ids: frozenset[str]
    sequence_lengths: dict[str, int]
    terminal_stop_ids: frozenset[str]


def _quantile(sorted_values: list[int], fraction: float) -> int | None:
    if not sorted_values:
        return None
    index = round((len(sorted_values) - 1) * fraction)
    return sorted_values[index]


def audit_fasta(
    path: str | Path,
    sequence_type: Literal["protein", "cds"],
) -> tuple[dict[str, object], FastaAuditIndex]:
    """Audit protein or CDS FASTA records."""

    path = Path(path)
    alphabet = PROTEIN_ALPHABET if sequence_type == "protein" else DNA_ALPHABET
    primary_ids: set[str] = set()
    duplicate_primary_ids: set[str] = set()
    relation_ids: set[str] = set()
    duplicate_relation_ids: set[str] = set()
    relation_source_counts: Counter[str] = Counter()
    relation_lengths: dict[str, int] = {}
    terminal_stop_ids: set[str] = set()
    lengths: list[int] = []
    invalid_character_records = 0
    invalid_character_counts: Counter[str] = Counter()
    invalid_character_examples: list[str] = []
    empty_records = 0
    internal_stop_records = 0
    terminal_stop_records = 0
    not_multiple_of_three = 0

    for record_id, header, sequence in iter_fasta(path):
        if record_id in primary_ids:
            duplicate_primary_ids.add(record_id)
        primary_ids.add(record_id)
        relation_id, relation_source = fasta_relation_id(record_id, header)
        relation_source_counts[relation_source] += 1
        if relation_id in relation_ids:
            duplicate_relation_ids.add(relation_id)
        relation_ids.add(relation_id)
        relation_lengths.setdefault(relation_id, len(sequence))
        sequence = sequence.upper()
        lengths.append(len(sequence))
        if not sequence:
            empty_records += 1
            continue
        invalid_characters = set(sequence) - alphabet
        if invalid_characters:
            invalid_character_records += 1
            invalid_character_counts.update(invalid_characters)
            if len(invalid_character_examples) < 10:
                invalid_character_examples.append(record_id)
        if sequence_type == "protein":
            has_terminal_stop = sequence.endswith("*")
            terminal_stop_records += int(has_terminal_stop)
            if has_terminal_stop:
                terminal_stop_ids.add(relation_id)
            internal_stop_records += int("*" in sequence[:-1])
        else:
            not_multiple_of_three += int(len(sequence) % 3 != 0)
            codons = [sequence[i : i + 3] for i in range(0, len(sequence) - 2, 3)]
            if codons:
                has_terminal_stop = codons[-1] in STOP_CODONS
                terminal_stop_records += int(has_terminal_stop)
                if has_terminal_stop:
                    terminal_stop_ids.add(relation_id)
                internal_stop_records += int(any(codon in STOP_CODONS for codon in codons[:-1]))

    sorted_lengths = sorted(lengths)
    report: dict[str, object] = {
        "path": str(path),
        "sequence_type": sequence_type,
        "sequence_count": len(lengths),
        "unique_ids": len(primary_ids),
        "duplicate_ids": len(duplicate_primary_ids),
        "relation_id_kind": "transcript",
        "unique_relation_ids": len(relation_ids),
        "duplicate_relation_ids": len(duplicate_relation_ids),
        "relation_id_source_counts": dict(sorted(relation_source_counts.items())),
        "empty_records": empty_records,
        "invalid_character_records": invalid_character_records,
        "invalid_character_counts": dict(sorted(invalid_character_counts.items())),
        "internal_stop_records": internal_stop_records,
        "terminal_stop_records": terminal_stop_records,
        "not_multiple_of_three": not_multiple_of_three if sequence_type == "cds" else None,
        "length": {
            "min": sorted_lengths[0] if sorted_lengths else None,
            "q25": _quantile(sorted_lengths, 0.25),
            "median": int(median(sorted_lengths)) if sorted_lengths else None,
            "q75": _quantile(sorted_lengths, 0.75),
            "max": sorted_lengths[-1] if sorted_lengths else None,
        },
        "examples": {
            "duplicate_ids": sorted(duplicate_primary_ids)[:10],
            "duplicate_relation_ids": sorted(duplicate_relation_ids)[:10],
            "invalid_character_ids": invalid_character_examples,
        },
    }
    index = FastaAuditIndex(
        relation_ids=frozenset(relation_ids),
        sequence_lengths=relation_lengths,
        terminal_stop_ids=frozenset(terminal_stop_ids),
    )
    return report, index


def _sha256(path: str | Path, block_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        while block := handle.read(block_size):
            digest.update(block)
    return digest.hexdigest()


def audit_bundle(
    *,
    gff_path: str | Path,
    protein_path: str | Path,
    cds_path: str | Path,
    fai_path: str | Path | None = None,
    checksums: bool = False,
) -> dict[str, object]:
    """Run the first PloidyPatch input-quality gate on one annotation bundle."""

    sequence_lengths = read_fai(fai_path) if fai_path else None
    gff_report, gff_ids = audit_gff(gff_path, sequence_lengths)
    protein_report, protein_index = audit_fasta(protein_path, "protein")
    cds_report, cds_index = audit_fasta(cds_path, "cds")

    transcript_ids = set(gff_ids.transcripts)
    expected_translation_ids = set(gff_ids.expected_translation_transcripts)
    protein_ids = set(protein_index.relation_ids)
    cds_ids = set(cds_index.relation_ids)
    protein_missing_for_transcript = expected_translation_ids - protein_ids
    protein_without_transcript = protein_ids - transcript_ids
    cds_missing_for_transcript = expected_translation_ids - cds_ids
    cds_without_transcript = cds_ids - transcript_ids
    comparable_translation_ids = expected_translation_ids & protein_ids & cds_ids
    translation_length_mismatches: list[dict[str, int | str]] = []
    frame_remainder_translation_consistent = 0
    frame_remainder_translation_mismatch = 0
    for transcript_id in sorted(comparable_translation_ids):
        protein_length = protein_index.sequence_lengths[transcript_id]
        cds_length = cds_index.sequence_lengths[transcript_id]
        translated_protein_length = protein_length - int(
            transcript_id in protein_index.terminal_stop_ids
        )
        expected_protein_length = cds_length // 3 - int(
            transcript_id in cds_index.terminal_stop_ids
        )
        lengths_match = translated_protein_length == expected_protein_length
        if cds_length % 3:
            if lengths_match:
                frame_remainder_translation_consistent += 1
            else:
                frame_remainder_translation_mismatch += 1
        if not lengths_match:
            translation_length_mismatches.append(
                {
                    "transcript_id": transcript_id,
                    "protein_length": protein_length,
                    "protein_length_without_terminal_stop": translated_protein_length,
                    "cds_length": cds_length,
                    "cds_floor_codons_without_terminal_stop": expected_protein_length,
                }
            )

    error_fields = {
        "gff_malformed_lines": int(gff_report["malformed_lines"]),
        "gff_invalid_coordinates": int(gff_report["invalid_coordinates"]),
        "gff_unknown_seqid_features": int(gff_report["unknown_seqid_features"]),
        "gff_out_of_bounds_features": int(gff_report["out_of_bounds_features"]),
        "gff_duplicate_gene_or_transcript_ids": int(
            gff_report["duplicate_gene_or_transcript_ids"]
        ),
        "gff_orphan_parent_ids": int(gff_report["orphan_parent_ids"]),
        "protein_duplicate_ids": int(protein_report["duplicate_ids"]),
        "protein_duplicate_relation_ids": int(
            protein_report["duplicate_relation_ids"]
        ),
        "protein_invalid_character_records": int(
            protein_report["invalid_character_records"]
        ),
        "cds_duplicate_ids": int(cds_report["duplicate_ids"]),
        "cds_duplicate_relation_ids": int(cds_report["duplicate_relation_ids"]),
        "cds_invalid_character_records": int(cds_report["invalid_character_records"]),
    }
    warning_fields = {
        "protein_internal_stop_records": int(protein_report["internal_stop_records"]),
        "cds_internal_stop_records": int(cds_report["internal_stop_records"]),
        "protein_coding_transcripts_without_gff_cds": int(
            gff_report["protein_coding_transcripts_without_cds"]
        ),
        "cds_protein_length_mismatches": len(translation_length_mismatches),
        "transcripts_missing_protein": len(protein_missing_for_transcript),
        "proteins_without_transcript": len(protein_without_transcript),
        "transcripts_missing_cds": len(cds_missing_for_transcript),
        "cds_without_transcript": len(cds_without_transcript),
    }
    observation_fields = {
        "cds_not_multiple_of_three": int(cds_report["not_multiple_of_three"] or 0),
        "cds_frame_remainder_translation_consistent": (
            frame_remainder_translation_consistent
        ),
        "cds_frame_remainder_translation_mismatch": (
            frame_remainder_translation_mismatch
        ),
    }
    if any(error_fields.values()):
        grade = "fail"
    elif any(warning_fields.values()):
        grade = "warn"
    else:
        grade = "pass"

    report: dict[str, object] = {
        "schema_version": "0.3.0",
        "assembly": {
            "fai_path": str(fai_path) if fai_path else None,
            "sequence_count": len(sequence_lengths) if sequence_lengths is not None else None,
            "total_bp": sum(sequence_lengths.values()) if sequence_lengths is not None else None,
            "longest_sequence_bp": max(sequence_lengths.values())
            if sequence_lengths
            else None,
        },
        "gff3": gff_report,
        "protein": protein_report,
        "cds": cds_report,
        "cross_checks": {
            "transcript_ids": len(transcript_ids),
            "expected_translation_transcript_ids": len(expected_translation_ids),
            "non_translation_transcript_ids": len(
                transcript_ids - expected_translation_ids
            ),
            "translation_length_checks": len(comparable_translation_ids),
            "cds_protein_length_mismatches": len(translation_length_mismatches),
            "cds_frame_remainder_translation_consistent": (
                frame_remainder_translation_consistent
            ),
            "cds_frame_remainder_translation_mismatch": (
                frame_remainder_translation_mismatch
            ),
            "transcripts_missing_protein": len(protein_missing_for_transcript),
            "proteins_without_transcript": len(protein_without_transcript),
            "transcripts_missing_cds": len(cds_missing_for_transcript),
            "cds_without_transcript": len(cds_without_transcript),
            "examples": {
                "transcripts_missing_protein": sorted(protein_missing_for_transcript)[:10],
                "proteins_without_transcript": sorted(protein_without_transcript)[:10],
                "transcripts_missing_cds": sorted(cds_missing_for_transcript)[:10],
                "cds_without_transcript": sorted(cds_without_transcript)[:10],
                "cds_protein_length_mismatches": translation_length_mismatches[:10],
            },
        },
        "quality_gate": {
            "grade": grade,
            "errors": error_fields,
            "warnings": warning_fields,
            "observations": observation_fields,
        },
    }
    if checksums:
        paths = [gff_path, protein_path, cds_path]
        if fai_path:
            paths.append(fai_path)
        report["sha256"] = {str(path): _sha256(path) for path in paths}
    return report


def write_json(report: dict[str, object], output: str | Path | None) -> None:
    payload = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if output:
        Path(output).write_text(payload, encoding="utf-8")
    else:
        print(payload, end="")
