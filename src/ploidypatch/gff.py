from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from urllib.parse import unquote

from .io import normalize_feature_id, open_text


TRANSCRIPT_TYPES = {"mRNA", "transcript"}
PROTEIN_CODING_BIOTYPES = {"protein_coding", "protein-coding"}


@dataclass(frozen=True)
class GffIdentifierSets:
    """Normalized identifier sets needed for cross-file consistency checks."""

    transcripts: frozenset[str]
    protein_coding_transcripts: frozenset[str]
    cds_parent_transcripts: frozenset[str]
    expected_translation_transcripts: frozenset[str]


def parse_attributes(text: str) -> tuple[dict[str, str], int]:
    """Parse GFF3 attributes and return the field map plus malformed count."""

    attributes: dict[str, str] = {}
    malformed = 0
    if text == "." or not text:
        return attributes, malformed
    for item in text.split(";"):
        item = item.strip()
        if not item:
            continue
        if "=" not in item:
            malformed += 1
            continue
        key, value = item.split("=", 1)
        key = unquote(key.strip())
        value = unquote(value.strip())
        if not key or key in attributes:
            malformed += 1
            continue
        attributes[key] = value
    return attributes, malformed


def _counter_dict(counter: Counter[str]) -> dict[str, int]:
    return dict(sorted(counter.items()))


def audit_gff(
    path: str | Path,
    sequence_lengths: dict[str, int] | None = None,
) -> tuple[dict[str, object], GffIdentifierSets]:
    """Audit GFF3 syntax, hierarchy, and coordinates.

    Returns the JSON-serializable report and normalized transcript identifiers.
    """

    path = Path(path)
    type_counts: Counter[str] = Counter()
    source_counts: Counter[str] = Counter()
    missing_id_by_type: Counter[str] = Counter()
    all_ids: set[str] = set()
    transcript_ids_raw: set[str] = set()
    transcript_ids_normalized: set[str] = set()
    protein_coding_transcript_ids: set[str] = set()
    transcript_biotype_counts: Counter[str] = Counter()
    gene_ids: set[str] = set()
    repeated_gene_or_transcript_ids: set[str] = set()
    parent_refs: set[str] = set()
    transcript_parents: dict[str, set[str]] = defaultdict(set)
    transcripts_with_cds: set[str] = set()
    transcripts_with_exon: set[str] = set()
    seqids_seen: set[str] = set()

    feature_lines = 0
    comment_lines = 0
    blank_lines = 0
    malformed_lines = 0
    malformed_attributes = 0
    invalid_coordinates = 0
    invalid_strand = 0
    invalid_phase = 0
    unknown_seqid_features = 0
    out_of_bounds_features = 0
    multi_parent_features = 0
    invalid_coordinate_examples: list[dict[str, object]] = []
    out_of_bounds_examples: list[dict[str, object]] = []

    with open_text(path) as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\n")
            if not line:
                blank_lines += 1
                continue
            if line.startswith("#"):
                comment_lines += 1
                continue
            fields = line.split("\t")
            if len(fields) != 9:
                malformed_lines += 1
                continue
            seqid, source, feature_type, start_text, end_text, _, strand, phase, attr = fields
            feature_lines += 1
            type_counts[feature_type] += 1
            source_counts[source] += 1
            seqids_seen.add(seqid)

            try:
                start = int(start_text)
                end = int(end_text)
                if start < 1 or end < start:
                    invalid_coordinates += 1
                    if len(invalid_coordinate_examples) < 10:
                        invalid_coordinate_examples.append(
                            {
                                "line_number": line_number,
                                "seqid": seqid,
                                "start": start,
                                "end": end,
                                "feature_type": feature_type,
                            }
                        )
            except ValueError:
                invalid_coordinates += 1
                start = end = -1
                if len(invalid_coordinate_examples) < 10:
                    invalid_coordinate_examples.append(
                        {
                            "line_number": line_number,
                            "seqid": seqid,
                            "start": start_text,
                            "end": end_text,
                            "feature_type": feature_type,
                        }
                    )

            if strand not in {"+", "-", ".", "?"}:
                invalid_strand += 1
            if phase not in {"0", "1", "2", "."}:
                invalid_phase += 1

            if sequence_lengths is not None:
                seq_length = sequence_lengths.get(seqid)
                if seq_length is None:
                    unknown_seqid_features += 1
                elif start >= 1 and end > seq_length:
                    out_of_bounds_features += 1
                    if len(out_of_bounds_examples) < 10:
                        out_of_bounds_examples.append(
                            {
                                "line_number": line_number,
                                "seqid": seqid,
                                "start": start,
                                "end": end,
                                "sequence_length": seq_length,
                                "feature_type": feature_type,
                            }
                        )

            attributes, malformed = parse_attributes(attr)
            malformed_attributes += malformed
            feature_id = attributes.get("ID")
            is_gene_or_tx = feature_type == "gene" or feature_type in TRANSCRIPT_TYPES
            if feature_id:
                if is_gene_or_tx and feature_id in all_ids:
                    repeated_gene_or_transcript_ids.add(feature_id)
                all_ids.add(feature_id)
                if feature_type == "gene":
                    gene_ids.add(feature_id)
                elif feature_type in TRANSCRIPT_TYPES:
                    transcript_ids_raw.add(feature_id)
                    normalized_id = normalize_feature_id(feature_id)
                    transcript_ids_normalized.add(normalized_id)
                    biotype = attributes.get("biotype") or attributes.get(
                        "transcript_biotype"
                    )
                    transcript_biotype_counts[biotype or "<missing>"] += 1
                    if biotype and biotype.lower() in PROTEIN_CODING_BIOTYPES:
                        protein_coding_transcript_ids.add(normalized_id)
            elif is_gene_or_tx:
                missing_id_by_type[feature_type] += 1

            parents = [p for p in attributes.get("Parent", "").split(",") if p]
            if len(parents) > 1:
                multi_parent_features += 1
            parent_refs.update(parents)
            if feature_type in TRANSCRIPT_TYPES and feature_id:
                transcript_parents[feature_id].update(parents)
            elif feature_type == "CDS":
                transcripts_with_cds.update(parents)
            elif feature_type == "exon":
                transcripts_with_exon.update(parents)

    orphan_parent_ids = parent_refs - all_ids
    genes_with_transcripts: Counter[str] = Counter()
    transcripts_without_gene_parent = 0
    for parents in transcript_parents.values():
        valid_gene_parents = parents & gene_ids
        if not valid_gene_parents:
            transcripts_without_gene_parent += 1
        for parent in valid_gene_parents:
            genes_with_transcripts[parent] += 1

    transcript_without_cds = transcript_ids_raw - transcripts_with_cds
    transcript_without_exon = transcript_ids_raw - transcripts_with_exon
    genes_without_transcript = gene_ids - set(genes_with_transcripts)
    cds_parent_transcript_ids = {
        normalize_feature_id(value) for value in transcripts_with_cds
    }
    expected_translation_transcript_ids = (
        protein_coding_transcript_ids | cds_parent_transcript_ids
    )
    protein_coding_without_cds = (
        protein_coding_transcript_ids - cds_parent_transcript_ids
    )

    report: dict[str, object] = {
        "path": str(path),
        "feature_lines": feature_lines,
        "comment_lines": comment_lines,
        "blank_lines": blank_lines,
        "type_counts": _counter_dict(type_counts),
        "source_counts": _counter_dict(source_counts),
        "sequence_ids_used": len(seqids_seen),
        "gene_ids": len(gene_ids),
        "transcript_ids": len(transcript_ids_raw),
        "transcript_biotype_counts": _counter_dict(transcript_biotype_counts),
        "protein_coding_transcripts": len(protein_coding_transcript_ids),
        "transcripts_with_cds": len(cds_parent_transcript_ids),
        "expected_translation_transcripts": len(
            expected_translation_transcript_ids
        ),
        "protein_coding_transcripts_without_cds": len(protein_coding_without_cds),
        "genes_with_multiple_transcripts": sum(
            count > 1 for count in genes_with_transcripts.values()
        ),
        "malformed_lines": malformed_lines,
        "malformed_attributes": malformed_attributes,
        "invalid_coordinates": invalid_coordinates,
        "invalid_strand": invalid_strand,
        "invalid_phase": invalid_phase,
        "unknown_seqid_features": unknown_seqid_features,
        "out_of_bounds_features": out_of_bounds_features,
        "duplicate_gene_or_transcript_ids": len(repeated_gene_or_transcript_ids),
        "missing_id_by_type": _counter_dict(missing_id_by_type),
        "orphan_parent_ids": len(orphan_parent_ids),
        "transcripts_without_gene_parent": transcripts_without_gene_parent,
        "transcripts_without_cds": len(transcript_without_cds),
        "transcripts_without_exon": len(transcript_without_exon),
        "genes_without_transcript": len(genes_without_transcript),
        "multi_parent_features": multi_parent_features,
        "examples": {
            "duplicate_gene_or_transcript_ids": sorted(repeated_gene_or_transcript_ids)[:10],
            "orphan_parent_ids": sorted(orphan_parent_ids)[:10],
            "transcripts_without_cds": sorted(transcript_without_cds)[:10],
            "protein_coding_transcripts_without_cds": sorted(
                protein_coding_without_cds
            )[:10],
            "unknown_sequence_ids": sorted(
                seqids_seen - set(sequence_lengths or {})
            )[:10]
            if sequence_lengths is not None
            else [],
            "invalid_coordinates": invalid_coordinate_examples,
            "out_of_bounds_features": out_of_bounds_examples,
        },
    }
    identifiers = GffIdentifierSets(
        transcripts=frozenset(transcript_ids_normalized),
        protein_coding_transcripts=frozenset(protein_coding_transcript_ids),
        cds_parent_transcripts=frozenset(cds_parent_transcript_ids),
        expected_translation_transcripts=frozenset(
            expected_translation_transcript_ids
        ),
    )
    return report, identifiers
