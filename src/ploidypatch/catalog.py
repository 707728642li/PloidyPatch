from __future__ import annotations

import bisect
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path
from typing import Any

from . import __version__
from .io import normalize_feature_id
from .perturb import (
    MISSING_GENE_EVENT,
    GffDocument,
    MissingGeneCandidate,
    candidate_id_for_gene,
    find_missing_gene_candidates,
    read_gff_document,
)
from .score import TranscriptSignature, build_annotation_index
from .synteny_io import discover_primary_chromosome_seqids


CATALOG_SCHEMA_VERSION = "ploidypatch.candidate_catalog.v1"
CATALOG_MANIFEST_SCHEMA_VERSION = "ploidypatch.candidate_catalog_manifest.v1"
CORE_COLUMNS = (
    "candidate_id",
    "event_type",
    "gene_id",
    "gene_id_normalized",
    "seqid",
    "start",
    "end",
    "strand",
    "gene_span_bp",
    "gene_span_bin",
    "gene_biotype",
    "transcript_count",
    "transcript_count_bin",
    "distinct_transcript_structures",
    "max_exons_per_transcript",
    "max_exons_bin",
    "max_cds_segments_per_transcript",
    "removed_record_count",
    "overlapping_other_genes",
    "genes_within_50kb",
    "structurally_unique",
)


def _file_sha256(path: str | Path, block_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        while block := handle.read(block_size):
            digest.update(block)
    return digest.hexdigest()


def _span_bin(span: int) -> str:
    if span < 1_000:
        return "lt_1kb"
    if span < 5_000:
        return "1_to_lt_5kb"
    if span < 20_000:
        return "5_to_lt_20kb"
    return "ge_20kb"


def _transcript_count_bin(count: int) -> str:
    if count == 1:
        return "1"
    if count == 2:
        return "2"
    return "ge_3"


def _exon_count_bin(count: int) -> str:
    if count == 1:
        return "1"
    if count <= 5:
        return "2_to_5"
    if count <= 10:
        return "6_to_10"
    return "ge_11"


def _interval_context(document: GffDocument) -> dict[int, tuple[int, int]]:
    genes_by_seqid: dict[str, list[Any]] = defaultdict(list)
    for record in document.records:
        if record.feature_type == "gene":
            genes_by_seqid[record.seqid].append(record)

    context: dict[int, tuple[int, int]] = {}
    for genes in genes_by_seqid.values():
        starts = sorted(record.start for record in genes)
        ends = sorted(record.end for record in genes)
        for gene in genes:
            overlapping = (
                bisect.bisect_right(starts, gene.end)
                - bisect.bisect_left(ends, gene.start)
                - 1
            )
            nearby = (
                bisect.bisect_right(starts, gene.end + 50_000)
                - bisect.bisect_left(ends, max(1, gene.start - 50_000))
                - 1
            )
            context[gene.line_number] = (max(0, overlapping), max(0, nearby))
    return context


def _signature_gene_ids(
    document: GffDocument,
) -> tuple[
    dict[str, TranscriptSignature],
    dict[TranscriptSignature, frozenset[str]],
]:
    index = build_annotation_index(document)
    signature_genes_mutable: dict[TranscriptSignature, set[str]] = defaultdict(set)
    transcript_signatures: dict[str, TranscriptSignature] = {}
    for transcript_id, model in index.transcripts.items():
        transcript_signatures[transcript_id] = model.signature
        signature_genes_mutable[model.signature].update(model.gene_ids)
    signature_genes = {
        signature: frozenset(gene_ids)
        for signature, gene_ids in signature_genes_mutable.items()
    }
    return transcript_signatures, signature_genes


def _candidate_row(
    candidate: MissingGeneCandidate,
    source_text_sha256: str,
    transcript_signatures: dict[str, TranscriptSignature],
    signature_gene_ids: dict[TranscriptSignature, frozenset[str]],
    interval_context: dict[int, tuple[int, int]],
) -> dict[str, str]:
    gene = candidate.gene
    gene_id = gene.feature_id
    if gene_id is None:
        raise AssertionError("Eligible candidate gene must have an ID")
    signatures = [transcript_signatures[value] for value in candidate.transcript_ids]
    max_exons = max(len(signature.exons) for signature in signatures)
    max_cds = max(len(signature.cds) for signature in signatures)
    structurally_unique = all(
        signature_gene_ids[signature] <= {gene_id} for signature in signatures
    )
    overlapping, nearby = interval_context[gene.line_number]
    span = gene.end - gene.start + 1
    biotype = gene.attributes.get("biotype") or gene.attributes.get("gene_biotype") or ""
    return {
        "candidate_id": candidate_id_for_gene(source_text_sha256, gene_id),
        "event_type": MISSING_GENE_EVENT,
        "gene_id": gene_id,
        "gene_id_normalized": normalize_feature_id(gene_id),
        "seqid": gene.seqid,
        "start": str(gene.start),
        "end": str(gene.end),
        "strand": gene.strand,
        "gene_span_bp": str(span),
        "gene_span_bin": _span_bin(span),
        "gene_biotype": biotype,
        "transcript_count": str(len(candidate.transcript_ids)),
        "transcript_count_bin": _transcript_count_bin(len(candidate.transcript_ids)),
        "distinct_transcript_structures": str(len(set(signatures))),
        "max_exons_per_transcript": str(max_exons),
        "max_exons_bin": _exon_count_bin(max_exons),
        "max_cds_segments_per_transcript": str(max_cds),
        "removed_record_count": str(len(candidate.removed_records)),
        "overlapping_other_genes": str(overlapping),
        "genes_within_50kb": str(nearby),
        "structurally_unique": "true" if structurally_unique else "false",
    }


def _read_external_strata(
    path: str | Path,
    column_prefix: str = "",
) -> tuple[tuple[str, ...], dict[str, dict[str, str]]]:
    path = Path(path)
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or "gene_id" not in reader.fieldnames:
            raise ValueError("External strata TSV must contain a gene_id column")
        source_columns = tuple(value for value in reader.fieldnames if value != "gene_id")
        extra_columns = tuple(f"{column_prefix}{value}" for value in source_columns)
        if len(set(extra_columns)) != len(extra_columns):
            raise ValueError("External strata columns are not unique after prefixing")
        collisions = set(extra_columns) & set(CORE_COLUMNS)
        if collisions:
            raise ValueError(
                "External strata columns collide with catalog columns: "
                + ", ".join(sorted(collisions))
            )
        rows: dict[str, dict[str, str]] = {}
        for line_number, row in enumerate(reader, start=2):
            gene_id = row.get("gene_id", "")
            if not gene_id:
                raise ValueError(f"Missing gene_id in external strata line {line_number}")
            if gene_id in rows:
                raise ValueError(f"Duplicate external strata gene_id: {gene_id}")
            rows[gene_id] = {
                output_column: row.get(source_column, "")
                for source_column, output_column in zip(
                    source_columns, extra_columns, strict=True
                )
            }
    return extra_columns, rows


def write_missing_gene_candidate_catalog(
    gff_path: str | Path,
    output_tsv_path: str | Path,
    external_strata_path: str | Path | None = None,
    external_strata_prefix: str = "",
    primary_chromosomes_only: bool = False,
) -> dict[str, Any]:
    """Write an evaluator-only candidate table for later stratified sampling."""

    source_path = Path(gff_path)
    output_path = Path(output_tsv_path)
    manifest_path = Path(str(output_path) + ".manifest.json")
    collisions = [path for path in (output_path, manifest_path) if path.exists()]
    if collisions:
        raise FileExistsError(
            "Refusing to overwrite candidate catalog artifact(s): "
            + ", ".join(str(path) for path in collisions)
        )

    chromosome_labels: dict[str, str] = {}
    selected_seqids: frozenset[str] | None = None
    if primary_chromosomes_only:
        selected_seqids, chromosome_labels = discover_primary_chromosome_seqids(
            source_path
        )
    document = read_gff_document(source_path, include_seqids=selected_seqids)
    source_text_sha256 = document.text_sha256
    candidates = find_missing_gene_candidates(document)
    transcript_signatures, signature_gene_ids = _signature_gene_ids(document)
    interval_context = _interval_context(document)
    rows = [
        _candidate_row(
            candidate,
            source_text_sha256,
            transcript_signatures,
            signature_gene_ids,
            interval_context,
        )
        for candidate in candidates
    ]

    extra_columns: tuple[str, ...] = ()
    external_rows: dict[str, dict[str, str]] = {}
    external_summary: dict[str, Any] | None = None
    if external_strata_path is not None:
        extra_columns, external_rows = _read_external_strata(
            external_strata_path, external_strata_prefix
        )
        catalog_gene_ids = {row["gene_id"] for row in rows}
        matched = catalog_gene_ids & set(external_rows)
        external_summary = {
            "file_name": Path(external_strata_path).name,
            "sha256": _file_sha256(external_strata_path),
            "rows": len(external_rows),
            "matched_candidates": len(matched),
            "candidates_without_external_strata": len(catalog_gene_ids - matched),
            "external_rows_not_eligible": len(set(external_rows) - catalog_gene_ids),
            "column_prefix": external_strata_prefix,
        }
        for row in rows:
            values = external_rows.get(row["gene_id"], {})
            row.update({column: values.get(column, "") for column in extra_columns})

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=(*CORE_COLUMNS, *extra_columns),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)

    structurally_ambiguous = sum(
        row["structurally_unique"] != "true" for row in rows
    )
    manifest: dict[str, Any] = {
        "schema_version": CATALOG_MANIFEST_SCHEMA_VERSION,
        "catalog_schema_version": CATALOG_SCHEMA_VERSION,
        "generator": {"name": "PloidyPatch", "version": __version__},
        "source": {
            "file_name": source_path.name,
            "file_sha256": _file_sha256(source_path),
            "text_sha256": source_text_sha256,
        },
        "event_type": MISSING_GENE_EVENT,
        "parameters": {
            "primary_chromosomes_only": primary_chromosomes_only,
            "chromosome_labels": chromosome_labels,
        },
        "catalog": {
            "file_name": output_path.name,
            "sha256": _file_sha256(output_path),
            "eligible_candidates": len(rows),
            "structurally_unique_candidates": len(rows) - structurally_ambiguous,
            "structurally_ambiguous_candidates": structurally_ambiguous,
        },
        "external_strata": external_summary,
        "access": "evaluator_only",
    }
    payload = json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    with manifest_path.open("x", encoding="utf-8", newline="") as handle:
        handle.write(payload)
    return manifest
