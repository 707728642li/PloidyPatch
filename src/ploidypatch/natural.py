from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

from . import __version__
from .baseline import (
    ProjectedModel,
    _file_sha256,
    _merge_intervals,
    _model_intervals,
    _overlap_bp,
    _parse_miniprot_models,
    _read_protein_map,
)
from .gff import TRANSCRIPT_TYPES
from .perturb import (
    BOUNDARY_SHIFT_EVENT,
    MISSING_GENE_EVENT,
    MISSING_INTERNAL_EXON_EVENT,
    SPLIT_GENE_EVENT,
    read_gff_document,
)
from .pav import _fetch_indexed_fasta, _read_fai_records
from .rna import JUNCTION_AGGREGATE_SCHEMA, JUNCTION_GROUP_AGGREGATE_SCHEMA


NATURAL_CANDIDATE_SCHEMA = "ploidypatch.natural_candidate_catalog.v2"
NATURAL_RNA_VALIDATION_SCHEMA = "ploidypatch.natural_rna_validation.v2"
NATURAL_SUMMARY_SCHEMA = "ploidypatch.natural_validation_summary.v2"
NATURAL_GRAPH_INPUT_SCHEMA = "ploidypatch.natural_graph_inputs.v1"
NATURAL_ASSEMBLY_CONTEXT_SCHEMA = "ploidypatch.natural_assembly_context.v1"
NATURAL_SECONDARY_RNA_SCHEMA = "ploidypatch.natural_secondary_rna_validation.v1"
ALTERNATIVE_STRUCTURE_EVENT = "annotation_alternative_cds_structure"
AMBIGUOUS_OPPOSITE_STRAND_EVENT = "ambiguous_opposite_strand_coding_overlap"


@dataclass(frozen=True)
class TargetGene:
    gene_id: str
    seqid: str
    strand: str
    start: int
    end: int
    cds: tuple[tuple[int, int], ...]
    transcript_chains: tuple[tuple[tuple[int, int], ...], ...]
    junctions: frozenset[tuple[str, int, int]]

    @property
    def cds_bp(self) -> int:
        return sum(end - start + 1 for start, end in self.cds)


@dataclass(frozen=True)
class GeneIntervalIndex:
    bins: dict[int, tuple[TargetGene, ...]]


GENE_BIN_SIZE = 100_000


def _junctions(
    seqid: str, intervals: Iterable[tuple[int, int]]
) -> frozenset[tuple[str, int, int]]:
    merged = _merge_intervals(list(intervals))
    return frozenset(
        (seqid, left_end, right_start)
        for (_, left_end), (right_start, _) in zip(merged, merged[1:])
        if right_start > left_end + 1
    )


def _read_target_genes(
    gff_path: str | Path,
) -> tuple[list[TargetGene], dict[tuple[str, str], GeneIntervalIndex]]:
    document = read_gff_document(gff_path)
    genes_by_id = {
        record.feature_id: record
        for record in document.records
        if record.feature_type == "gene" and record.feature_id
    }
    transcript_to_gene: dict[str, str] = {}
    for record in document.records:
        if record.feature_type not in TRANSCRIPT_TYPES or not record.feature_id:
            continue
        gene_parents = [parent for parent in record.parents if parent in genes_by_id]
        if len(gene_parents) == 1:
            transcript_to_gene[record.feature_id] = gene_parents[0]

    transcript_cds: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for record in document.records:
        if record.feature_type != "CDS":
            continue
        for parent in record.parents:
            if parent in transcript_to_gene:
                transcript_cds[parent].append((record.start, record.end))

    gene_transcripts: dict[str, list[tuple[tuple[int, int], ...]]] = defaultdict(list)
    for transcript_id, cds in transcript_cds.items():
        if not cds:
            continue
        gene_id = transcript_to_gene[transcript_id]
        gene_transcripts[gene_id].append(tuple(_merge_intervals(cds)))

    genes: list[TargetGene] = []
    for gene_id, transcript_chains in gene_transcripts.items():
        record = genes_by_id[gene_id]
        cds = tuple(
            _merge_intervals(
                [interval for chain in transcript_chains for interval in chain]
            )
        )
        genes.append(
            TargetGene(
                gene_id=gene_id,
                seqid=record.seqid,
                strand=record.strand,
                start=record.start,
                end=record.end,
                cds=cds,
                transcript_chains=tuple(sorted(set(transcript_chains))),
                junctions=frozenset(
                    junction
                    for chain in transcript_chains
                    for junction in _junctions(record.seqid, chain)
                ),
            )
        )

    grouped: dict[tuple[str, str], list[TargetGene]] = defaultdict(list)
    for gene in genes:
        grouped[(gene.seqid, gene.strand)].append(gene)
    indexes: dict[tuple[str, str], GeneIntervalIndex] = {}
    for key, values in grouped.items():
        raw_bins: dict[int, list[TargetGene]] = defaultdict(list)
        for gene in values:
            for bin_number in range(
                gene.start // GENE_BIN_SIZE,
                gene.end // GENE_BIN_SIZE + 1,
            ):
                raw_bins[bin_number].append(gene)
        indexes[key] = GeneIntervalIndex(
            bins={
                bin_number: tuple(
                    sorted(bin_genes, key=lambda gene: (gene.start, gene.end))
                )
                for bin_number, bin_genes in raw_bins.items()
            }
        )
    return genes, indexes


def _overlapping_genes(
    model: ProjectedModel,
    index: GeneIntervalIndex | None,
) -> list[tuple[TargetGene, int]]:
    if index is None:
        return []
    model_intervals = _model_intervals(model)
    nearby: dict[str, TargetGene] = {}
    for bin_number in range(
        model.start // GENE_BIN_SIZE,
        model.end // GENE_BIN_SIZE + 1,
    ):
        for gene in index.bins.get(bin_number, ()):
            nearby[gene.gene_id] = gene
    overlaps: list[tuple[TargetGene, int]] = []
    for gene in nearby.values():
        if gene.end < model.start:
            continue
        overlap = _overlap_bp(model_intervals, list(gene.cds))
        if overlap:
            overlaps.append((gene, overlap))
    return overlaps


def _classify_model(
    model: ProjectedModel,
    overlaps: list[tuple[TargetGene, int]],
    opposite_strand_overlaps: list[tuple[TargetGene, int]],
    *,
    max_existing_cds_overlap: float,
    min_boundary_extension_bp: int,
) -> dict[str, Any] | None:
    intervals = _model_intervals(model)
    model_bp = sum(end - start + 1 for start, end in intervals)
    total_overlap = _overlap_bp(
        intervals,
        _merge_intervals(
            [interval for gene, _ in overlaps for interval in gene.cds]
        ),
    )
    if total_overlap / model_bp <= max_existing_cds_overlap:
        opposite_overlap = _overlap_bp(
            intervals,
            _merge_intervals(
                [
                    interval
                    for gene, _ in opposite_strand_overlaps
                    for interval in gene.cds
                ]
            ),
        )
        if opposite_overlap / model_bp > max_existing_cds_overlap:
            event_type = AMBIGUOUS_OPPOSITE_STRAND_EVENT
            relevant = sorted(
                (
                    gene
                    for gene, overlap in opposite_strand_overlaps
                    if overlap / model_bp > max_existing_cds_overlap
                ),
                key=lambda gene: (gene.start, gene.gene_id),
            )
        else:
            event_type = MISSING_GENE_EVENT
            relevant = []
    else:
        significant = [
            gene
            for gene, overlap in overlaps
            if overlap / gene.cds_bp >= 0.2 and overlap / model_bp >= 0.1
        ]
        if len(significant) >= 2:
            event_type = SPLIT_GENE_EVENT
            relevant = sorted(significant, key=lambda gene: (gene.start, gene.gene_id))
        elif overlaps:
            best, _ = max(
                overlaps,
                key=lambda item: (
                    item[1] / model_bp,
                    item[1] / item[0].cds_bp,
                    item[1],
                ),
            )
            relevant = [best]
            if tuple(intervals) in best.transcript_chains:
                return None
            novel_segments = [
                interval
                for interval in intervals
                if _overlap_bp([interval], list(best.cds)) == 0
            ]
            internal = [
                interval
                for interval in novel_segments
                if interval[0] > best.cds[0][0] and interval[1] < best.cds[-1][1]
            ]
            extension = max(
                max(0, best.cds[0][0] - intervals[0][0]),
                max(0, intervals[-1][1] - best.cds[-1][1]),
            )
            if internal:
                event_type = MISSING_INTERNAL_EXON_EVENT
            elif extension >= min_boundary_extension_bp:
                event_type = BOUNDARY_SHIFT_EVENT
            else:
                event_type = ALTERNATIVE_STRUCTURE_EVENT
        else:
            event_type = MISSING_GENE_EVENT
            relevant = []

    target_cds = _merge_intervals(
        [interval for gene in relevant for interval in gene.cds]
    )
    target_junctions = frozenset(
        junction for gene in relevant for junction in gene.junctions
    )
    projected_junctions = _junctions(model.seqid, intervals)
    novel_segments = [
        interval
        for interval in intervals
        if not target_cds or _overlap_bp([interval], target_cds) == 0
    ]
    return {
        "event_type": event_type,
        "seqid": model.seqid,
        "start": model.start,
        "end": model.end,
        "strand": model.strand,
        "target_gene_ids": tuple(gene.gene_id for gene in relevant),
        "source": model.source,
        "query_id": model.query_id,
        "model_id": model.model_id,
        "identity": model.identity,
        "query_coverage": model.query_coverage,
        "frameshifts": model.frameshifts,
        "stop_codons": model.stop_codons,
        "projected_cds_segments": tuple(intervals),
        "projected_junctions": projected_junctions,
        "novel_junctions": projected_junctions - target_junctions,
        "novel_cds_segments": tuple(novel_segments),
    }


def _reciprocal_span_overlap(left: dict[str, Any], right: dict[str, Any]) -> float:
    overlap = max(
        0,
        min(left["end"], right["end"])
        - max(left["start"], right["start"])
        + 1,
    )
    denominator = min(
        left["end"] - left["start"] + 1,
        right["end"] - right["start"] + 1,
    )
    return overlap / denominator


def _cluster_candidates(candidates: list[dict[str, Any]]) -> list[list[dict[str, Any]]]:
    grouped: dict[tuple[Any, ...], list[dict[str, Any]]] = defaultdict(list)
    for candidate in candidates:
        grouped[
            (
                candidate["event_type"],
                candidate["seqid"],
                candidate["strand"],
                candidate["target_gene_ids"],
            )
        ].append(candidate)
    clusters: list[list[dict[str, Any]]] = []
    for key in sorted(grouped):
        current: list[dict[str, Any]] = []
        representative: dict[str, Any] | None = None
        for candidate in sorted(
            grouped[key],
            key=lambda item: (item["start"], item["end"], item["model_id"]),
        ):
            if (
                representative is None
                or _reciprocal_span_overlap(representative, candidate) < 0.5
            ):
                if current:
                    clusters.append(current)
                current = [candidate]
                representative = candidate
            else:
                current.append(candidate)
                representative = max(
                    current,
                    key=lambda item: (
                        item["identity"],
                        item["query_coverage"],
                        item["end"] - item["start"],
                    ),
                )
        if current:
            clusters.append(current)
    return clusters


def discover_natural_candidates(
    *,
    target_gff_path: str | Path,
    miniprot_gff_path: str | Path,
    protein_map_path: str | Path,
    output_tsv_path: str | Path,
    min_identity: float = 0.8,
    min_query_coverage: float = 0.8,
    max_existing_cds_overlap: float = 0.1,
    min_boundary_extension_bp: int = 30,
    near_best_score_fraction: float = 0.95,
) -> dict[str, Any]:
    """Discover RNA-blind natural structure hypotheses from haplotype proteins."""

    for name, value in (
        ("min_identity", min_identity),
        ("min_query_coverage", min_query_coverage),
        ("max_existing_cds_overlap", max_existing_cds_overlap),
        ("near_best_score_fraction", near_best_score_fraction),
    ):
        if not 0 <= value <= 1:
            raise ValueError(f"{name} must be between zero and one")
    if min_boundary_extension_bp < 1:
        raise ValueError("min_boundary_extension_bp must be positive")
    output = Path(output_tsv_path)
    manifest_path = Path(str(output) + ".manifest.json")
    if output.exists() or manifest_path.exists():
        raise FileExistsError(f"Refusing to overwrite natural candidates: {output}")

    protein_map = _read_protein_map(protein_map_path)
    models = _parse_miniprot_models(miniprot_gff_path, protein_map)
    _, gene_indexes = _read_target_genes(target_gff_path)
    query_scores: dict[str, list[float]] = defaultdict(list)
    for model in models:
        query_scores[model.query_id].append(model.score)
    query_mapping: dict[str, tuple[int, int]] = {}
    for query_id, scores in query_scores.items():
        best_score = max(scores)
        near_best = sum(
            score >= best_score * near_best_score_fraction for score in scores
        )
        query_mapping[query_id] = (len(scores), near_best)
    filter_counts: Counter[str] = Counter()
    preliminary: list[dict[str, Any]] = []
    for model in models:
        if model.rank != 1:
            filter_counts["non_primary_rank"] += 1
            continue
        if model.identity < min_identity:
            filter_counts["low_identity"] += 1
            continue
        if model.query_coverage < min_query_coverage:
            filter_counts["low_query_coverage"] += 1
            continue
        if model.frameshifts or model.stop_codons:
            filter_counts["disrupted_projection"] += 1
            continue
        if not model.cds:
            filter_counts["no_cds"] += 1
            continue
        overlaps = _overlapping_genes(
            model, gene_indexes.get((model.seqid, model.strand))
        )
        opposite_strand = "-" if model.strand == "+" else "+"
        opposite_strand_overlaps = _overlapping_genes(
            model, gene_indexes.get((model.seqid, opposite_strand))
        )
        candidate = _classify_model(
            model,
            overlaps,
            opposite_strand_overlaps,
            max_existing_cds_overlap=max_existing_cds_overlap,
            min_boundary_extension_bp=min_boundary_extension_bp,
        )
        if candidate is None:
            filter_counts["concordant_target_chain"] += 1
            continue
        (
            candidate["query_projection_count"],
            candidate["query_near_best_projection_count"],
        ) = query_mapping[model.query_id]
        preliminary.append(candidate)

    rows: list[dict[str, Any]] = []
    event_counts: Counter[str] = Counter()
    tier_counts: Counter[str] = Counter()
    for cluster in _cluster_candidates(preliminary):
        representative = max(
            cluster,
            key=lambda item: (
                item["identity"],
                item["query_coverage"],
                -item["frameshifts"],
                -item["stop_codons"],
                item["end"] - item["start"],
                item["model_id"],
            ),
        )
        sources = sorted({item["source"] for item in cluster})
        model_ids = sorted({item["model_id"] for item in cluster})
        query_ids = sorted({item["query_id"] for item in cluster})
        novel_junctions = sorted(
            {junction for item in cluster for junction in item["novel_junctions"]}
        )
        projected_junctions = sorted(
            {junction for item in cluster for junction in item["projected_junctions"]}
        )
        novel_segments = sorted(
            {segment for item in cluster for segment in item["novel_cds_segments"]}
        )
        event_key = "\t".join(
            (
                representative["event_type"],
                representative["seqid"],
                str(min(item["start"] for item in cluster)),
                str(max(item["end"] for item in cluster)),
                representative["strand"],
                ",".join(representative["target_gene_ids"]),
                ",".join(model_ids),
            )
        )
        candidate_id = "PPN-" + hashlib.sha256(event_key.encode()).hexdigest()[:20]
        tier = (
            "two_haplotype_sources"
            if len(sources) >= 2
            else "single_haplotype_source"
        )
        max_near_best = max(
            item["query_near_best_projection_count"] for item in cluster
        )
        mapping_tier = (
            "all_queries_unique_best_locus"
            if max_near_best == 1
            else "one_or_more_queries_multilocus"
        )
        event_counts[representative["event_type"]] += 1
        tier_counts[tier] += 1
        rows.append(
            {
                "candidate_id": candidate_id,
                "event_type": representative["event_type"],
                "seqid": representative["seqid"],
                "start": min(item["start"] for item in cluster),
                "end": max(item["end"] for item in cluster),
                "strand": representative["strand"],
                "target_gene_ids": ",".join(representative["target_gene_ids"]),
                "support_sources": ",".join(sources),
                "support_source_count": len(sources),
                "support_model_count": len(model_ids),
                "query_ids": ",".join(query_ids),
                "model_ids": ",".join(model_ids),
                "min_identity": min(item["identity"] for item in cluster),
                "min_query_coverage": min(item["query_coverage"] for item in cluster),
                "max_query_projection_count": max(
                    item["query_projection_count"] for item in cluster
                ),
                "max_query_near_best_projection_count": max_near_best,
                "mapping_specificity_tier": mapping_tier,
                "projected_cds_segments_json": json.dumps(
                    representative["projected_cds_segments"], separators=(",", ":")
                ),
                "projected_junction_count": len(projected_junctions),
                "novel_junction_count": len(novel_junctions),
                "novel_junctions_json": json.dumps(
                    novel_junctions, separators=(",", ":")
                ),
                "novel_cds_segments_json": json.dumps(
                    novel_segments, separators=(",", ":")
                ),
                "discovery_tier": tier,
                "rna_used_for_discovery": "false",
                "automatic_patch_policy": "review_required",
            }
        )

    fields = (
        "candidate_id",
        "event_type",
        "seqid",
        "start",
        "end",
        "strand",
        "target_gene_ids",
        "support_sources",
        "support_source_count",
        "support_model_count",
        "query_ids",
        "model_ids",
        "min_identity",
        "min_query_coverage",
        "max_query_projection_count",
        "max_query_near_best_projection_count",
        "mapping_specificity_tier",
        "projected_cds_segments_json",
        "projected_junction_count",
        "novel_junction_count",
        "novel_junctions_json",
        "novel_cds_segments_json",
        "discovery_tier",
        "rna_used_for_discovery",
        "automatic_patch_policy",
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(
            sorted(
                rows,
                key=lambda row: (row["seqid"], row["start"], row["candidate_id"]),
            )
        )

    manifest = {
        "schema_version": NATURAL_CANDIDATE_SCHEMA,
        "generator": {"name": "PloidyPatch", "version": __version__},
        "inputs": {
            "target_gff_sha256": _file_sha256(target_gff_path),
            "miniprot_gff_sha256": _file_sha256(miniprot_gff_path),
            "protein_map_sha256": _file_sha256(protein_map_path),
        },
        "parameters": {
            "min_identity": min_identity,
            "min_query_coverage": min_query_coverage,
            "max_existing_cds_overlap": max_existing_cds_overlap,
            "min_boundary_extension_bp": min_boundary_extension_bp,
            "near_best_score_fraction": near_best_score_fraction,
            "require_rank_1": True,
            "require_intact": True,
            "rna_used_for_discovery": False,
            "automatic_patch_policy": "review_required",
        },
        "counts": {
            "input_models": len(models),
            "preliminary_hypotheses": len(preliminary),
            "candidate_loci": len(rows),
            "filter_counts": dict(sorted(filter_counts.items())),
            "event_counts": dict(sorted(event_counts.items())),
            "discovery_tier_counts": dict(sorted(tier_counts.items())),
        },
        "output": {
            "file_name": output.name,
            "rows": len(rows),
            "sha256": _file_sha256(output),
        },
    }
    with manifest_path.open("x", encoding="utf-8", newline="") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return manifest


def validate_natural_candidates_with_rna(
    *,
    candidate_tsv_path: str | Path,
    junction_aggregate_tsv_path: str | Path,
    output_tsv_path: str | Path,
) -> dict[str, Any]:
    """Join a frozen RNA-blind catalog to independently aggregated junctions."""

    candidates_path = Path(candidate_tsv_path)
    junctions_path = Path(junction_aggregate_tsv_path)
    output = Path(output_tsv_path)
    manifest_path = Path(str(output) + ".manifest.json")
    if output.exists() or manifest_path.exists():
        raise FileExistsError(f"Refusing to overwrite RNA validation: {output}")

    candidate_manifest_path = Path(str(candidates_path) + ".manifest.json")
    junction_manifest_path = Path(str(junctions_path) + ".manifest.json")
    candidate_manifest = json.loads(candidate_manifest_path.read_text(encoding="utf-8"))
    junction_manifest = json.loads(junction_manifest_path.read_text(encoding="utf-8"))
    if candidate_manifest.get("schema_version") != NATURAL_CANDIDATE_SCHEMA:
        raise ValueError("Unsupported natural candidate manifest")
    if junction_manifest.get("schema_version") not in {
        JUNCTION_AGGREGATE_SCHEMA,
        JUNCTION_GROUP_AGGREGATE_SCHEMA,
    }:
        raise ValueError("Unsupported junction aggregate manifest")
    if candidate_manifest.get("output", {}).get("sha256") != _file_sha256(
        candidates_path
    ):
        raise ValueError("Natural candidate TSV checksum mismatch")
    if junction_manifest.get("output", {}).get("sha256") != _file_sha256(
        junctions_path
    ):
        raise ValueError("Junction aggregate TSV checksum mismatch")

    junction_support: dict[tuple[str, int, int], dict[str, str]] = {}
    with junctions_path.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            key = (
                row["seqid"],
                int(row["left_exon_end"]),
                int(row["right_exon_start"]),
            )
            junction_support[key] = row

    with candidates_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        candidate_fields = tuple(reader.fieldnames or ())
        candidates = list(reader)
    added_fields = (
        "rna_validation_state",
        "novel_junctions_tested",
        "primary_supported_novel_junctions",
        "primary_supported_fraction",
        "primary_supported_junction_read_sum",
        "rna_validation_tier",
        "rna_negative_evidence_policy",
    )
    state_counts: Counter[str] = Counter()
    tier_counts: Counter[str] = Counter()
    validated_rows: list[dict[str, Any]] = []
    for row in candidates:
        raw_junctions = json.loads(row["novel_junctions_json"])
        junctions = [
            (str(seqid), int(left), int(right))
            for seqid, left, right in raw_junctions
        ]
        supported = [
            junction_support[junction]
            for junction in junctions
            if junction in junction_support
            and junction_support[junction]["primary_support"] == "true"
        ]
        if not junctions:
            state = "not_assessable_no_novel_splice_junction"
        elif len(supported) == len(junctions):
            state = "all_novel_junctions_supported"
        elif supported:
            state = "some_novel_junctions_supported"
        else:
            state = "no_qualifying_junction_observed"
        source_count = int(row["support_source_count"])
        mapping_unique = row["mapping_specificity_tier"] == (
            "all_queries_unique_best_locus"
        )
        if (
            state == "all_novel_junctions_supported"
            and source_count >= 2
            and mapping_unique
        ):
            tier = "two_haplotype_unique_projection_plus_primary_rna"
        elif state == "all_novel_junctions_supported" and source_count >= 2:
            tier = "two_haplotype_multilocus_projection_plus_primary_rna"
        elif state == "all_novel_junctions_supported" and mapping_unique:
            tier = "single_haplotype_unique_projection_plus_primary_rna"
        elif state == "all_novel_junctions_supported":
            tier = "single_haplotype_multilocus_projection_plus_primary_rna"
        elif supported:
            tier = "partial_primary_rna_support"
        elif not junctions:
            tier = "homology_only_rna_not_assessable"
        else:
            tier = "homology_only_rna_missing"
        state_counts[state] += 1
        tier_counts[tier] += 1
        output_row: dict[str, Any] = dict(row)
        output_row.update(
            {
                "rna_validation_state": state,
                "novel_junctions_tested": len(junctions),
                "primary_supported_novel_junctions": len(supported),
                "primary_supported_fraction": (
                    len(supported) / len(junctions) if junctions else ""
                ),
                "primary_supported_junction_read_sum": sum(
                    int(value["primary_read_count"]) for value in supported
                ),
                "rna_validation_tier": tier,
                "rna_negative_evidence_policy": "absence_is_missing_not_contradiction",
            }
        )
        validated_rows.append(output_row)

    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=candidate_fields + added_fields, delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(validated_rows)
    manifest = {
        "schema_version": NATURAL_RNA_VALIDATION_SCHEMA,
        "generator": {"name": "PloidyPatch", "version": __version__},
        "inputs": {
            "candidate_tsv_sha256": _file_sha256(candidates_path),
            "candidate_manifest_sha256": _file_sha256(candidate_manifest_path),
            "junction_aggregate_tsv_sha256": _file_sha256(junctions_path),
            "junction_aggregate_manifest_sha256": _file_sha256(junction_manifest_path),
        },
        "policy": {
            "candidate_discovery_used_rna": False,
            "rna_is_held_out_validation": True,
            "strand_policy": "unstranded",
            "negative_evidence_policy": "absence_is_missing_not_contradiction",
            "automatic_patch_policy": "review_required",
        },
        "counts": {
            "candidates": len(validated_rows),
            "validation_state_counts": dict(sorted(state_counts.items())),
            "validation_tier_counts": dict(sorted(tier_counts.items())),
        },
        "output": {
            "file_name": output.name,
            "rows": len(validated_rows),
            "sha256": _file_sha256(output),
        },
    }
    with manifest_path.open("x", encoding="utf-8", newline="") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return manifest


def annotate_natural_assembly_context(
    *,
    validation_tsv_path: str | Path,
    genome_fasta_path: str | Path,
    genome_fai_path: str | Path,
    output_tsv_path: str | Path,
    flank_bp: int = 5000,
    max_ambiguous_fraction: float = 0.0,
) -> dict[str, Any]:
    """Add sequence-edge, ambiguity, and soft-mask context to frozen candidates."""

    if flank_bp < 0:
        raise ValueError("flank_bp must be non-negative")
    if not 0 <= max_ambiguous_fraction <= 1:
        raise ValueError("max_ambiguous_fraction must be within [0, 1]")
    validation = Path(validation_tsv_path)
    validation_manifest_path = Path(str(validation) + ".manifest.json")
    genome = Path(genome_fasta_path)
    genome_fai = Path(genome_fai_path)
    output = Path(output_tsv_path)
    manifest_path = Path(str(output) + ".manifest.json")
    if output.exists() or manifest_path.exists():
        raise FileExistsError(f"Refusing to overwrite assembly context: {output}")

    validation_manifest = json.loads(
        validation_manifest_path.read_text(encoding="utf-8")
    )
    if validation_manifest.get("schema_version") != NATURAL_RNA_VALIDATION_SCHEMA:
        raise ValueError("Unsupported natural RNA validation manifest")
    if validation_manifest.get("output", {}).get("sha256") != _file_sha256(
        validation
    ):
        raise ValueError("Natural RNA validation TSV checksum mismatch")

    fai_records = _read_fai_records(genome_fai)
    with validation.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        input_fields = tuple(reader.fieldnames or ())
        required = {"candidate_id", "seqid", "start", "end"}
        missing = required - set(input_fields)
        if missing:
            raise ValueError(
                "Natural validation lacks assembly-context column(s): "
                + ", ".join(sorted(missing))
            )
        rows = list(reader)

    added_fields = (
        "assembly_sequence_length",
        "distance_to_sequence_start_bp",
        "distance_to_sequence_end_bp",
        "left_flank_observed_bp",
        "right_flank_observed_bp",
        "locus_ambiguous_bp",
        "locus_ambiguous_fraction",
        "context_ambiguous_bp",
        "context_ambiguous_fraction",
        "locus_softmasked_bp",
        "locus_softmasked_fraction",
        "context_softmasked_bp",
        "context_softmasked_fraction",
        "assembly_context_state",
    )
    state_counts: Counter[str] = Counter()
    output_rows: list[dict[str, Any]] = []
    with genome.open("rb") as genome_handle:
        for line_number, row in enumerate(rows, start=2):
            seqid = row["seqid"]
            record = fai_records.get(seqid)
            if record is None:
                raise ValueError(
                    f"Candidate sequence {seqid!r} is absent from FAI at line "
                    f"{line_number}"
                )
            sequence_length = record[0]
            try:
                start = int(row["start"])
                end = int(row["end"])
            except ValueError as exc:
                raise ValueError(
                    f"Non-integer candidate interval at line {line_number}"
                ) from exc
            if start < 1 or end < start or end > sequence_length:
                raise ValueError(
                    f"Candidate interval outside {seqid} at line {line_number}: "
                    f"{start}-{end} / {sequence_length}"
                )

            locus = _fetch_indexed_fasta(
                genome_handle, record, start - 1, end
            )
            context_start = max(0, start - 1 - flank_bp)
            context_end = min(sequence_length, end + flank_bp)
            context = _fetch_indexed_fasta(
                genome_handle, record, context_start, context_end
            )
            locus_ambiguous = sum(base.upper() not in "ACGT" for base in locus)
            context_ambiguous = sum(
                base.upper() not in "ACGT" for base in context
            )
            locus_softmasked = sum(base in "acgt" for base in locus)
            context_softmasked = sum(base in "acgt" for base in context)
            locus_ambiguous_fraction = locus_ambiguous / len(locus)
            context_ambiguous_fraction = context_ambiguous / len(context)
            left_flank = start - 1 - context_start
            right_flank = context_end - end
            if locus_ambiguous_fraction > max_ambiguous_fraction:
                state = "locus_ambiguous_sequence"
            elif left_flank < flank_bp or right_flank < flank_bp:
                state = "truncated_at_sequence_edge"
            elif context_ambiguous_fraction > max_ambiguous_fraction:
                state = "flanking_ambiguous_sequence"
            else:
                state = "pass"
            state_counts[state] += 1
            output_row: dict[str, Any] = dict(row)
            output_row.update(
                {
                    "assembly_sequence_length": sequence_length,
                    "distance_to_sequence_start_bp": start - 1,
                    "distance_to_sequence_end_bp": sequence_length - end,
                    "left_flank_observed_bp": left_flank,
                    "right_flank_observed_bp": right_flank,
                    "locus_ambiguous_bp": locus_ambiguous,
                    "locus_ambiguous_fraction": f"{locus_ambiguous_fraction:.8f}",
                    "context_ambiguous_bp": context_ambiguous,
                    "context_ambiguous_fraction": f"{context_ambiguous_fraction:.8f}",
                    "locus_softmasked_bp": locus_softmasked,
                    "locus_softmasked_fraction": f"{locus_softmasked / len(locus):.8f}",
                    "context_softmasked_bp": context_softmasked,
                    "context_softmasked_fraction": f"{context_softmasked / len(context):.8f}",
                    "assembly_context_state": state,
                }
            )
            output_rows.append(output_row)

    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=input_fields + added_fields, delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(output_rows)
    manifest = {
        "schema_version": NATURAL_ASSEMBLY_CONTEXT_SCHEMA,
        "generator": {"name": "PloidyPatch", "version": __version__},
        "inputs": {
            "validation_tsv_sha256": _file_sha256(validation),
            "validation_manifest_sha256": _file_sha256(validation_manifest_path),
            "genome_fasta_sha256": _file_sha256(genome),
            "genome_fai_sha256": _file_sha256(genome_fai),
        },
        "parameters": {
            "flank_bp": flank_bp,
            "max_ambiguous_fraction": max_ambiguous_fraction,
            "ambiguous_definition": "base_not_ACGT_case_insensitive",
            "softmask_definition": "lowercase_acgt",
            "softmask_policy": "context_only_uncalibrated",
        },
        "counts": {
            "candidates": len(output_rows),
            "assembly_context_state_counts": dict(sorted(state_counts.items())),
        },
        "output": {
            "file_name": output.name,
            "rows": len(output_rows),
            "sha256": _file_sha256(output),
        },
    }
    with manifest_path.open("x", encoding="utf-8", newline="") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return manifest


def validate_natural_candidates_with_secondary_groups(
    *,
    primary_validation_tsv_path: str | Path,
    grouped_junction_tsv_path: str | Path,
    output_tsv_path: str | Path,
) -> dict[str, Any]:
    """Join frozen candidates to conservative non-primary filename-stem groups."""

    validation = Path(primary_validation_tsv_path)
    validation_manifest_path = Path(str(validation) + ".manifest.json")
    junctions = Path(grouped_junction_tsv_path)
    junction_manifest_path = Path(str(junctions) + ".manifest.json")
    output = Path(output_tsv_path)
    manifest_path = Path(str(output) + ".manifest.json")
    if output.exists() or manifest_path.exists():
        raise FileExistsError(
            f"Refusing to overwrite secondary RNA validation: {output}"
        )
    validation_manifest = json.loads(
        validation_manifest_path.read_text(encoding="utf-8")
    )
    junction_manifest = json.loads(junction_manifest_path.read_text(encoding="utf-8"))
    input_validation_schema = validation_manifest.get("schema_version")
    if input_validation_schema not in {
        NATURAL_RNA_VALIDATION_SCHEMA,
        NATURAL_ASSEMBLY_CONTEXT_SCHEMA,
    }:
        raise ValueError("Unsupported primary natural validation manifest")
    if junction_manifest.get("schema_version") != JUNCTION_GROUP_AGGREGATE_SCHEMA:
        raise ValueError("Secondary validation requires a grouped junction aggregate")
    if validation_manifest.get("output", {}).get("sha256") != _file_sha256(
        validation
    ):
        raise ValueError("Primary natural RNA validation TSV checksum mismatch")
    if junction_manifest.get("output", {}).get("sha256") != _file_sha256(
        junctions
    ):
        raise ValueError("Grouped junction aggregate TSV checksum mismatch")

    junction_support: dict[tuple[str, int, int], dict[str, str]] = {}
    with junctions.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required_junction_fields = {
            "seqid",
            "left_exon_end",
            "right_exon_start",
            "secondary_read_count",
            "secondary_groups_ge_threshold",
            "secondary_group_support",
        }
        missing = required_junction_fields - set(reader.fieldnames or ())
        if missing:
            raise ValueError(
                "Grouped junction aggregate lacks column(s): "
                + ", ".join(sorted(missing))
            )
        for row in reader:
            key = (
                row["seqid"],
                int(row["left_exon_end"]),
                int(row["right_exon_start"]),
            )
            junction_support[key] = row

    with validation.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        input_fields = tuple(reader.fieldnames or ())
        required = {
            "candidate_id",
            "event_type",
            "novel_junctions_json",
            "rna_validation_state",
        }
        missing = required - set(input_fields)
        if missing:
            raise ValueError(
                "Primary natural validation lacks column(s): "
                + ", ".join(sorted(missing))
            )
        rows = list(reader)

    added_fields = (
        "secondary_group_validation_state",
        "secondary_group_supported_novel_junctions",
        "secondary_group_supported_fraction",
        "secondary_group_min_supporting_groups",
        "secondary_group_junction_read_sum",
        "secondary_group_validation_tier",
        "secondary_group_negative_evidence_policy",
    )
    state_counts: Counter[str] = Counter()
    tier_counts: Counter[str] = Counter()
    event_state_counts: Counter[tuple[str, str]] = Counter()
    output_rows: list[dict[str, Any]] = []
    for row in rows:
        raw_junctions = json.loads(row["novel_junctions_json"])
        candidate_junctions = [
            (str(seqid), int(left), int(right))
            for seqid, left, right in raw_junctions
        ]
        supported = [
            junction_support[junction]
            for junction in candidate_junctions
            if junction in junction_support
            and junction_support[junction]["secondary_group_support"] == "true"
        ]
        if not candidate_junctions:
            state = "not_assessable_no_novel_splice_junction"
        elif len(supported) == len(candidate_junctions):
            state = "all_novel_junctions_group_recurrent"
        elif supported:
            state = "some_novel_junctions_group_recurrent"
        else:
            state = "no_group_recurrent_junction_observed"
        if not candidate_junctions:
            tier = "secondary_groups_not_assessable"
        elif (
            state == "all_novel_junctions_group_recurrent"
            and row["rna_validation_state"] == "all_novel_junctions_supported"
        ):
            tier = "primary_and_secondary_groups_all_junctions"
        elif state == "all_novel_junctions_group_recurrent":
            tier = "secondary_groups_all_primary_not_all"
        elif supported:
            tier = "secondary_groups_partial"
        else:
            tier = "secondary_groups_missing"
        state_counts[state] += 1
        tier_counts[tier] += 1
        event_state_counts[(row["event_type"], state)] += 1
        output_row: dict[str, Any] = dict(row)
        output_row.update(
            {
                "secondary_group_validation_state": state,
                "secondary_group_supported_novel_junctions": len(supported),
                "secondary_group_supported_fraction": (
                    len(supported) / len(candidate_junctions)
                    if candidate_junctions
                    else ""
                ),
                "secondary_group_min_supporting_groups": (
                    min(int(value["secondary_groups_ge_threshold"]) for value in supported)
                    if supported
                    else 0
                ),
                "secondary_group_junction_read_sum": sum(
                    int(value["secondary_read_count"]) for value in supported
                ),
                "secondary_group_validation_tier": tier,
                "secondary_group_negative_evidence_policy": (
                    "absence_is_missing_not_contradiction"
                ),
            }
        )
        output_rows.append(output_row)

    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=input_fields + added_fields, delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(output_rows)
    manifest = {
        "schema_version": NATURAL_SECONDARY_RNA_SCHEMA,
        "generator": {"name": "PloidyPatch", "version": __version__},
        "inputs": {
            "primary_validation_tsv_sha256": _file_sha256(validation),
            "primary_validation_manifest_sha256": _file_sha256(
                validation_manifest_path
            ),
            "grouped_junction_tsv_sha256": _file_sha256(junctions),
            "grouped_junction_manifest_sha256": _file_sha256(
                junction_manifest_path
            ),
        },
        "policy": {
            "input_validation_schema": input_validation_schema,
            "candidate_discovery_used_secondary_rna": False,
            "sample_group_interpretation": (
                "filename_stem_only_not_biological_metadata"
            ),
            "secondary_group_support_is_context_not_truth": True,
            "negative_evidence_policy": "absence_is_missing_not_contradiction",
            "automatic_patch_policy": "review_required",
        },
        "counts": {
            "candidates": len(output_rows),
            "validation_state_counts": dict(sorted(state_counts.items())),
            "validation_tier_counts": dict(sorted(tier_counts.items())),
            "event_by_state": [
                {"event_type": event, "state": state, "count": count}
                for (event, state), count in sorted(event_state_counts.items())
            ],
        },
        "output": {
            "file_name": output.name,
            "rows": len(output_rows),
            "sha256": _file_sha256(output),
        },
    }
    with manifest_path.open("x", encoding="utf-8", newline="") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return manifest


def summarize_natural_validation(
    *,
    validation_tsv_path: str | Path,
    output_json_path: str | Path,
) -> dict[str, Any]:
    """Write auditable event-by-evidence cross counts for a validation table."""

    validation = Path(validation_tsv_path)
    validation_manifest_path = Path(str(validation) + ".manifest.json")
    output = Path(output_json_path)
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite natural summary: {output}")
    validation_manifest = json.loads(
        validation_manifest_path.read_text(encoding="utf-8")
    )
    if validation_manifest.get("schema_version") != NATURAL_RNA_VALIDATION_SCHEMA:
        raise ValueError("Unsupported natural RNA validation manifest")
    if validation_manifest.get("output", {}).get("sha256") != _file_sha256(
        validation
    ):
        raise ValueError("Natural RNA validation TSV checksum mismatch")

    columns = (
        "event_type",
        "discovery_tier",
        "rna_validation_state",
        "rna_validation_tier",
    )
    crosses = (
        ("event_type", "rna_validation_state"),
        ("event_type", "rna_validation_tier"),
        ("discovery_tier", "rna_validation_state"),
    )
    one_way = {column: Counter() for column in columns}
    joint = {cross: Counter() for cross in crosses}
    candidate_ids: set[str] = set()
    with validation.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"candidate_id", *columns}
        missing = required - set(reader.fieldnames or ())
        if missing:
            raise ValueError(
                "Natural validation lacks summary column(s): "
                + ", ".join(sorted(missing))
            )
        for line_number, row in enumerate(reader, start=2):
            candidate_id = row["candidate_id"]
            if not candidate_id or candidate_id in candidate_ids:
                raise ValueError(
                    f"Empty or duplicate candidate ID at line {line_number}"
                )
            candidate_ids.add(candidate_id)
            for column in columns:
                one_way[column][row[column] or "__MISSING__"] += 1
            for cross in crosses:
                joint[cross][(row[cross[0]], row[cross[1]])] += 1

    report = {
        "schema_version": NATURAL_SUMMARY_SCHEMA,
        "generator": {"name": "PloidyPatch", "version": __version__},
        "input": {
            "file_name": validation.name,
            "rows": len(candidate_ids),
            "sha256": _file_sha256(validation),
            "manifest_sha256": _file_sha256(validation_manifest_path),
        },
        "one_way_counts": {
            column: dict(sorted(counter.items()))
            for column, counter in one_way.items()
        },
        "joint_counts": {
            f"{left}__x__{right}": [
                {left: values[0], right: values[1], "count": count}
                for values, count in sorted(counter.items())
            ]
            for (left, right), counter in joint.items()
        },
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x", encoding="utf-8", newline="") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return report


def prepare_natural_graph_inputs(
    *,
    validation_tsv_path: str | Path,
    output_candidate_tsv_path: str | Path,
    output_evidence_tsv_path: str | Path,
) -> dict[str, Any]:
    """Adapt held-out natural validation into correlation-aware graph inputs."""

    validation = Path(validation_tsv_path)
    validation_manifest_path = Path(str(validation) + ".manifest.json")
    candidate_output = Path(output_candidate_tsv_path)
    evidence_output = Path(output_evidence_tsv_path)
    manifest_path = Path(str(candidate_output) + ".manifest.json")
    collisions = [
        path
        for path in (candidate_output, evidence_output, manifest_path)
        if path.exists()
    ]
    if collisions:
        raise FileExistsError(
            "Refusing to overwrite natural graph input(s): "
            + ", ".join(str(path) for path in collisions)
        )
    validation_manifest = json.loads(
        validation_manifest_path.read_text(encoding="utf-8")
    )
    input_schema = validation_manifest.get("schema_version")
    if input_schema not in {
        NATURAL_RNA_VALIDATION_SCHEMA,
        NATURAL_ASSEMBLY_CONTEXT_SCHEMA,
        NATURAL_SECONDARY_RNA_SCHEMA,
    }:
        raise ValueError("Unsupported natural validation manifest")
    if validation_manifest.get("output", {}).get("sha256") != _file_sha256(
        validation
    ):
        raise ValueError("Natural RNA validation TSV checksum mismatch")

    event_state = {
        MISSING_GENE_EVENT: "missing_annotation",
        MISSING_INTERNAL_EXON_EVENT: "missing_exon",
        BOUNDARY_SHIFT_EVENT: "boundary_error",
        SPLIT_GENE_EVENT: "split_error",
        ALTERNATIVE_STRUCTURE_EVENT: "uncertain",
        AMBIGUOUS_OPPOSITE_STRAND_EVENT: "uncertain",
    }
    with validation.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        required = {
            "candidate_id",
            "event_type",
            "target_gene_ids",
            "query_ids",
            "seqid",
            "start",
            "end",
            "strand",
            "support_sources",
            "min_identity",
            "min_query_coverage",
            "mapping_specificity_tier",
            "novel_junctions_tested",
            "primary_supported_novel_junctions",
            "rna_validation_state",
        }
        missing = required - set(reader.fieldnames or ())
        if missing:
            raise ValueError(
                "Natural validation lacks graph column(s): "
                + ", ".join(sorted(missing))
            )

    candidate_fields = (
        "event_id",
        "event_type",
        "candidate_id",
        "gene_ids",
        "seqid",
        "start",
        "end",
        "strand",
        "relationship_state",
        "wgd_event",
        "subgenome",
        "proposal_uri",
    )
    evidence_fields = (
        "candidate_id",
        "evidence_id",
        "evidence_type",
        "direction",
        "scope",
        "strength",
        "reliability",
        "independent_group",
        "source",
        "details",
    )
    candidate_rows: list[dict[str, Any]] = []
    evidence_rows: list[dict[str, Any]] = []
    graph_state_counts: Counter[str] = Counter()
    evidence_type_counts: Counter[str] = Counter()
    input_sha = _file_sha256(validation)
    for row in rows:
        candidate_id = row["candidate_id"]
        state = event_state.get(row["event_type"])
        if state is None:
            raise ValueError(f"Unsupported natural event type: {row['event_type']}")
        gene_ids = row["target_gene_ids"] or row["query_ids"]
        candidate_rows.append(
            {
                "event_id": candidate_id,
                "event_type": state,
                "candidate_id": candidate_id,
                "gene_ids": gene_ids,
                "seqid": row["seqid"],
                "start": row["start"],
                "end": row["end"],
                "strand": row["strand"],
                "relationship_state": "",
                "wgd_event": "",
                "subgenome": "",
                "proposal_uri": f"sha256:{input_sha}#{candidate_id}",
            }
        )
        graph_state_counts[state] += 1
        neutral = state == "uncertain"
        protein_strength = min(
            float(row["min_identity"]), float(row["min_query_coverage"])
        )
        evidence_rows.append(
            {
                "candidate_id": candidate_id,
                "evidence_id": f"{candidate_id}:protein_projection",
                "evidence_type": "protein_projection",
                "direction": "context" if neutral else "support",
                "scope": "coding",
                "strength": protein_strength,
                "reliability": 0.69,
                "independent_group": "black_haplotype_projection_annotation",
                "source": "Black_Hap1_and_Black_Hap2_miniprot",
                "details": json.dumps(
                    {
                        "support_sources": row["support_sources"],
                        "correlation": "same_material_related_annotation_pipelines",
                    },
                    separators=(",", ":"),
                ),
            }
        )
        evidence_type_counts["protein_projection"] += 1
        tested = int(row["novel_junctions_tested"])
        supported = int(row["primary_supported_novel_junctions"])
        if tested and supported:
            evidence_rows.append(
                {
                    "candidate_id": candidate_id,
                    "evidence_id": f"{candidate_id}:primary_rna_junction",
                    "evidence_type": "rna_junction",
                    "direction": "context" if neutral else "support",
                    "scope": "both",
                    "strength": supported / tested,
                    "reliability": 0.8,
                    "independent_group": "black_primary_short_read_rna",
                    "source": "Black-1_Black-2_Black-3_junction_aggregate",
                    "details": row["rna_validation_state"],
                }
            )
            evidence_type_counts["rna_junction"] += 1
        if "secondary_group_supported_novel_junctions" in row:
            secondary_supported = int(
                row["secondary_group_supported_novel_junctions"]
            )
            if tested and secondary_supported:
                evidence_rows.append(
                    {
                        "candidate_id": candidate_id,
                        "evidence_id": (
                            f"{candidate_id}:secondary_group_rna_junction"
                        ),
                        "evidence_type": "rna_junction",
                        "direction": "context",
                        "scope": "both",
                        "strength": secondary_supported / tested,
                        "reliability": 0.6,
                        "independent_group": (
                            "nonblack_filename_stem_short_read_rna"
                        ),
                        "source": "57_nonprimary_BAM_group_aggregate",
                        "details": json.dumps(
                            {
                                "state": row[
                                    "secondary_group_validation_state"
                                ],
                                "minimum_supporting_groups": row[
                                    "secondary_group_min_supporting_groups"
                                ],
                                "junction_read_sum": row[
                                    "secondary_group_junction_read_sum"
                                ],
                                "interpretation": (
                                    "context_only_filename_stem_groups"
                                ),
                            },
                            separators=(",", ":"),
                        ),
                    }
                )
                evidence_type_counts["rna_junction"] += 1
        if row["mapping_specificity_tier"] != "all_queries_unique_best_locus":
            evidence_rows.append(
                {
                    "candidate_id": candidate_id,
                    "evidence_id": f"{candidate_id}:projection_specificity",
                    "evidence_type": "mappability",
                    "direction": "context" if neutral else "contradict",
                    "scope": "coding",
                    "strength": 1.0,
                    "reliability": 0.8,
                    "independent_group": "protein_projection_specificity",
                    "source": "miniprot_near_best_projection_count",
                    "details": row["mapping_specificity_tier"],
                }
            )
            evidence_type_counts["mappability"] += 1
        if "assembly_context_state" in row:
            assembly_state = row["assembly_context_state"]
            if assembly_state != "pass":
                evidence_rows.append(
                    {
                        "candidate_id": candidate_id,
                        "evidence_id": f"{candidate_id}:assembly_context",
                        "evidence_type": "assembly_gap",
                        "direction": "context" if neutral else "contradict",
                        "scope": "both",
                        "strength": 1.0,
                        "reliability": 0.95,
                        "independent_group": "black_primary_assembly_sequence",
                        "source": "Black_Primary_FASTA_indexed_context",
                        "details": json.dumps(
                            {
                                "state": assembly_state,
                                "locus_ambiguous_fraction": row[
                                    "locus_ambiguous_fraction"
                                ],
                                "context_ambiguous_fraction": row[
                                    "context_ambiguous_fraction"
                                ],
                                "left_flank_observed_bp": row[
                                    "left_flank_observed_bp"
                                ],
                                "right_flank_observed_bp": row[
                                    "right_flank_observed_bp"
                                ],
                            },
                            separators=(",", ":"),
                        ),
                    }
                )
                evidence_type_counts["assembly_gap"] += 1
            repeat_fraction = max(
                float(row["locus_softmasked_fraction"]),
                float(row["context_softmasked_fraction"]),
            )
            if repeat_fraction > 0:
                evidence_rows.append(
                    {
                        "candidate_id": candidate_id,
                        "evidence_id": f"{candidate_id}:softmask_context",
                        "evidence_type": "repeat_overlap",
                        "direction": "context",
                        "scope": "both",
                        "strength": repeat_fraction,
                        "reliability": 0.5,
                        "independent_group": "black_primary_assembly_softmask",
                        "source": "Black_Primary_FASTA_lowercase_mask",
                        "details": "uncalibrated_context_only",
                    }
                )
                evidence_type_counts["repeat_overlap"] += 1

    candidate_output.parent.mkdir(parents=True, exist_ok=True)
    with candidate_output.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=candidate_fields, delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(candidate_rows)
    evidence_output.parent.mkdir(parents=True, exist_ok=True)
    with evidence_output.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=evidence_fields, delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(evidence_rows)

    manifest = {
        "schema_version": NATURAL_GRAPH_INPUT_SCHEMA,
        "generator": {"name": "PloidyPatch", "version": __version__},
        "input": {
            "validation_tsv_sha256": input_sha,
            "validation_manifest_sha256": _file_sha256(validation_manifest_path),
        },
        "policy": {
            "haplotype_projection_reliability": 0.69,
            "primary_short_read_rna_reliability": 0.8,
            "secondary_group_short_read_rna_reliability": 0.6,
            "secondary_group_rna_context_only": True,
            "multilocus_projection_reliability": 0.8,
            "haplotype_sources_share_one_independent_group": True,
            "rna_absence_generates_contradiction": False,
            "neutral_events_emit_context_only": True,
            "assembly_uncertainty_reliability": 0.95,
            "softmask_policy": "context_only_uncalibrated",
            "confidence_status": "transparent_uncalibrated",
        },
        "counts": {
            "candidates": len(candidate_rows),
            "evidence_edges": len(evidence_rows),
            "graph_state_counts": dict(sorted(graph_state_counts.items())),
            "evidence_type_counts": dict(sorted(evidence_type_counts.items())),
        },
        "outputs": {
            "candidates": {
                "file_name": candidate_output.name,
                "sha256": _file_sha256(candidate_output),
            },
            "evidence": {
                "file_name": evidence_output.name,
                "sha256": _file_sha256(evidence_output),
            },
        },
    }
    with manifest_path.open("x", encoding="utf-8", newline="") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return manifest
