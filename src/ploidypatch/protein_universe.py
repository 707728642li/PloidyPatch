"""Build exact protein-supported annotation universes for synteny analyses.

The source annotation is never rewritten in place.  This module creates a
separate, checksummed subset containing only genes with an exact provider
protein relation.  Provider-wide exact structural duplicate rows may be
collapsed, while multi-locus gene/transcript identities and fuzzy mappings
remain hard errors.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
import json
import os
from pathlib import Path
import shutil
import tempfile
from typing import Any, Iterable

from .artifact_manifest import sha256_file, verify_sha256sums, write_sha256sums
from .audit import PROTEIN_ALPHABET
from .io import iter_fasta, normalize_feature_id, parse_fasta_header_fields
from .perturb import GffRecord, read_gff_document
from .synteny_io import _normalized_gene_records


PROTEIN_UNIVERSE_SCHEMA = "ploidypatch.protein_supported_gene_universe.v1"


@dataclass(frozen=True)
class ProteinRepresentative:
    gene_id: str
    protein_id: str
    relation_source: str
    sequence: str


def _normalized_parents(record: GffRecord) -> tuple[str, ...]:
    return tuple(sorted({normalize_feature_id(value) for value in record.parents}))


def _structural_key(record: GffRecord) -> tuple[Any, ...]:
    """Return the exact biological identity used for duplicate collapse."""

    feature_id = normalize_feature_id(record.feature_id) if record.feature_id else ""
    return (
        record.feature_type,
        feature_id,
        record.seqid,
        record.start,
        record.end,
        record.strand,
        record.phase,
        _normalized_parents(record),
    )


def collapse_exact_structural_duplicates(
    records: Iterable[GffRecord],
) -> tuple[tuple[GffRecord, ...], dict[str, Any]]:
    """Collapse exact-coordinate hierarchy duplicates and reject identity drift.

    Attribute variants that do not alter ID, Parent, coordinates, strand or
    phase are recorded but cannot create another WGDI gene.  Gene and
    transcript IDs are required to identify one structural locus.  Repeated
    exon/CDS IDs at distinct coordinates remain valid multipart features.
    """

    groups: dict[tuple[Any, ...], list[GffRecord]] = defaultdict(list)
    identity_structures: dict[tuple[str, str], set[tuple[Any, ...]]] = defaultdict(set)
    input_count = 0
    for record in records:
        input_count += 1
        key = _structural_key(record)
        groups[key].append(record)
        if record.feature_type in {"gene", "mRNA", "transcript"} and record.feature_id:
            identity = (record.feature_type, normalize_feature_id(record.feature_id))
            identity_structures[identity].add(key)

    collisions = {
        identity: values
        for identity, values in identity_structures.items()
        if len(values) != 1
    }
    if collisions:
        examples = ",".join(
            f"{feature_type}:{feature_id}"
            for feature_type, feature_id in sorted(collisions)[:10]
        )
        raise ValueError(
            "Gene/transcript identity occurs at multiple structures after exact "
            f"normalization: {examples}"
        )

    selected: list[GffRecord] = []
    duplicate_rows = 0
    duplicate_groups = 0
    attribute_variant_groups = 0
    feature_duplicate_rows: Counter[str] = Counter()
    for values in groups.values():
        if len(values) > 1:
            duplicate_groups += 1
            duplicate_rows += len(values) - 1
            feature_duplicate_rows[values[0].feature_type] += len(values) - 1
            if len({value.raw_line.rstrip("\r\n").split("\t", 8)[8] for value in values}) > 1:
                attribute_variant_groups += 1
        # Lexical selection makes volatile provider attributes deterministic.
        selected.append(min(values, key=lambda value: value.raw_line))

    selected.sort(key=lambda value: value.line_number)
    audit: dict[str, Any] = {
        "input_feature_rows": input_count,
        "retained_feature_rows_after_exact_collapse": len(selected),
        "collapsed_duplicate_rows": duplicate_rows,
        "duplicate_structural_groups": duplicate_groups,
        "duplicate_groups_with_attribute_variants": attribute_variant_groups,
        "collapsed_rows_by_feature_type": dict(sorted(feature_duplicate_rows.items())),
        "gene_or_transcript_multilocus_collisions": 0,
    }
    return tuple(selected), audit


def _coding_gene_relations(
    records: tuple[GffRecord, ...],
    genes: dict[str, GffRecord],
    transcript_to_gene: dict[str, str],
) -> tuple[set[str], dict[str, str], int, int]:
    coding_genes: set[str] = set()
    protein_to_genes: dict[str, set[str]] = defaultdict(set)
    orphan_without_parent = 0
    unresolved_with_parent = 0
    for record in records:
        if record.feature_type != "CDS":
            continue
        parents = {normalize_feature_id(parent) for parent in record.parents}
        if not parents:
            # A parentless provider fragment is not a gene model and cannot be
            # rescued by its ID/coordinates.  It is excluded and counted.
            orphan_without_parent += 1
            continue
        mapped = {
            transcript_to_gene[parent]
            for parent in parents
            if parent in transcript_to_gene
        } | {parent for parent in parents if parent in genes}
        if len(mapped) != 1:
            unresolved_with_parent += 1
            continue
        gene_id = next(iter(mapped))
        coding_genes.add(gene_id)
        protein_id = record.attributes.get("protein_id")
        if protein_id:
            protein_to_genes[normalize_feature_id(protein_id)].add(gene_id)
    ambiguous = {
        protein_id: mapped
        for protein_id, mapped in protein_to_genes.items()
        if len(mapped) != 1
    }
    if ambiguous:
        raise ValueError("A GFF CDS protein_id maps to multiple normalized genes")
    return (
        coding_genes,
        {protein_id: next(iter(mapped)) for protein_id, mapped in protein_to_genes.items()},
        orphan_without_parent,
        unresolved_with_parent,
    )


def _provider_representatives(
    protein_path: str | Path,
    *,
    genes: dict[str, GffRecord],
    transcript_to_gene: dict[str, str],
    protein_to_gene: dict[str, str],
    coding_genes: set[str],
) -> tuple[dict[str, ProteinRepresentative], dict[str, Any]]:
    representatives: dict[str, ProteinRepresentative] = {}
    provider_ids: set[str] = set()
    provider_records = 0
    unmapped = 0
    mapped_non_coding = 0
    alternatives = 0
    relation_sources: Counter[str] = Counter()
    for protein_id, header, raw_sequence in iter_fasta(protein_path):
        provider_records += 1
        if protein_id in provider_ids:
            raise ValueError(f"Duplicate provider protein identifier: {protein_id}")
        provider_ids.add(protein_id)
        fields = parse_fasta_header_fields(header)
        candidates: dict[str, str] = {}
        explicit_gene = fields.get("gene") or fields.get("gene_id")
        if explicit_gene:
            gene_id = normalize_feature_id(explicit_gene)
            if gene_id in genes:
                candidates[gene_id] = "header:gene"
        explicit_transcript = fields.get("transcript") or fields.get("transcript_id")
        if explicit_transcript:
            transcript_id = normalize_feature_id(explicit_transcript)
            if transcript_id in transcript_to_gene:
                candidates[transcript_to_gene[transcript_id]] = (
                    "header:transcript->gff:Parent"
                )
        if protein_id in protein_to_gene:
            candidates[protein_to_gene[protein_id]] = (
                "FASTA_first_token->GFF_CDS:protein_id"
            )
        if protein_id in transcript_to_gene:
            candidates[transcript_to_gene[protein_id]] = (
                "FASTA_first_token->GFF_transcript"
            )
        if len(candidates) > 1:
            raise ValueError(
                f"Provider protein has discordant exact gene mappings: {protein_id}"
            )
        if not candidates:
            unmapped += 1
            continue
        gene_id, relation_source = next(iter(candidates.items()))
        if gene_id not in coding_genes:
            mapped_non_coding += 1
            continue
        sequence = raw_sequence.upper().rstrip("*")
        if not sequence or set(sequence) - PROTEIN_ALPHABET:
            raise ValueError(f"Invalid provider protein sequence: {protein_id}")
        relation_sources[relation_source] += 1
        candidate = ProteinRepresentative(
            gene_id=gene_id,
            protein_id=protein_id,
            relation_source=relation_source,
            sequence=sequence,
        )
        current = representatives.get(gene_id)
        if current is None or (-len(sequence), protein_id) < (
            -len(current.sequence),
            current.protein_id,
        ):
            if current is not None:
                alternatives += 1
            representatives[gene_id] = candidate
        else:
            alternatives += 1
    selected_ids = [item.protein_id for item in representatives.values()]
    if len(selected_ids) != len(set(selected_ids)):
        raise ValueError("One provider protein maps to multiple genes")
    return representatives, {
        "provider_protein_records": provider_records,
        "unmapped_provider_proteins": unmapped,
        "provider_proteins_mapped_to_non_coding_gene": mapped_non_coding,
        "mapped_alternative_isoforms_not_selected": alternatives,
        "relation_source_counts": dict(sorted(relation_sources.items())),
    }


def _record_gene_id(
    record: GffRecord,
    *,
    genes: dict[str, GffRecord],
    transcript_to_gene: dict[str, str],
) -> str | None:
    if record.feature_type == "gene" and record.feature_id:
        return normalize_feature_id(record.feature_id)
    if record.feature_type in {"mRNA", "transcript"} and record.feature_id:
        return transcript_to_gene.get(normalize_feature_id(record.feature_id))
    mapped = {
        transcript_to_gene[parent]
        for parent in _normalized_parents(record)
        if parent in transcript_to_gene
    } | {parent for parent in _normalized_parents(record) if parent in genes}
    if len(mapped) > 1:
        raise ValueError(
            f"Feature has parents in multiple normalized genes on line {record.line_number}"
        )
    return next(iter(mapped)) if mapped else None


def _write_tsv(path: Path, header: tuple[str, ...], rows: Iterable[tuple[Any, ...]]) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        handle.write("\t".join(header) + "\n")
        for row in rows:
            handle.write("\t".join(str(value) for value in row) + "\n")


def build_protein_supported_gene_universe(
    *,
    primary_gff_path: str | Path,
    provider_protein_path: str | Path,
    output_dir: str | Path,
    species_id: str,
    holdout_id: str,
    policy_id: str,
) -> dict[str, Any]:
    """Create an atomic, exact protein-supported GFF/protein subset."""

    source_gff = Path(primary_gff_path)
    source_protein = Path(provider_protein_path)
    output = Path(output_dir)
    for path in (source_gff, source_protein):
        if not path.is_file() or path.is_symlink() or path.stat().st_size == 0:
            raise ValueError(f"Missing, empty or symlinked protein-universe input: {path}")
    if output.exists() or output.is_symlink():
        raise FileExistsError(f"Refusing to overwrite protein universe: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    working = Path(tempfile.mkdtemp(prefix=f".{output.name}.working.", dir=output.parent))
    try:
        document = read_gff_document(source_gff)
        records, collapse_audit = collapse_exact_structural_duplicates(document.records)
        genes, transcript_to_gene, hierarchy = _normalized_gene_records(records)
        (
            coding_genes,
            protein_to_gene,
            orphan_cds_without_parent,
            unresolved_cds_with_parent,
        ) = _coding_gene_relations(
            records, genes, transcript_to_gene
        )
        if not coding_genes or unresolved_cds_with_parent:
            raise ValueError(
                "Protein-universe CDS hierarchy is not uniquely gene-resolved: "
                f"coding_genes={len(coding_genes)}, "
                f"unresolved_CDS_with_parent={unresolved_cds_with_parent}"
            )
        representatives, provider_audit = _provider_representatives(
            source_protein,
            genes=genes,
            transcript_to_gene=transcript_to_gene,
            protein_to_gene=protein_to_gene,
            coding_genes=coding_genes,
        )
        selected_genes = frozenset(representatives)
        if not selected_genes:
            raise ValueError("No coding gene has an exact provider protein")

        retained_records = [
            record
            for record in records
            if _record_gene_id(
                record, genes=genes, transcript_to_gene=transcript_to_gene
            )
            in selected_genes
        ]
        retained_record_ids = {
            normalize_feature_id(record.feature_id)
            for record in retained_records
            if record.feature_id
        }
        for record in retained_records:
            missing_parents = {
                parent
                for parent in _normalized_parents(record)
                if parent not in retained_record_ids
            }
            if missing_parents:
                raise ValueError(
                    f"Retained feature has missing parents on line {record.line_number}"
                )

        gff_output = working / "protein_supported.gff3"
        with gff_output.open("x", encoding="utf-8", newline="") as handle:
            handle.write("##gff-version 3\n")
            for record in retained_records:
                handle.write(record.raw_line.rstrip("\r\n") + "\n")

        protein_output = working / "representative.protein.fa"
        with protein_output.open("x", encoding="utf-8", newline="") as handle:
            for gene_id in sorted(selected_genes):
                representative = representatives[gene_id]
                handle.write(f">{representative.protein_id} gene:{gene_id}\n")
                for offset in range(0, len(representative.sequence), 60):
                    handle.write(representative.sequence[offset : offset + 60] + "\n")

        representatives_output = working / "representatives.tsv"
        _write_tsv(
            representatives_output,
            ("gene_id", "protein_id", "relation_source", "protein_length"),
            (
                (
                    gene_id,
                    representatives[gene_id].protein_id,
                    representatives[gene_id].relation_source,
                    len(representatives[gene_id].sequence),
                )
                for gene_id in sorted(selected_genes)
            ),
        )
        excluded_output = working / "excluded_genes.tsv"
        _write_tsv(
            excluded_output,
            ("gene_id", "reason"),
            (
                (gene_id, "no_exact_provider_protein")
                for gene_id in sorted(coding_genes - selected_genes)
            ),
        )

        output_proteins = list(iter_fasta(protein_output))
        if len(output_proteins) != len(selected_genes):
            raise ValueError("Representative protein output is not one record per gene")
        manifest: dict[str, Any] = {
            "schema_version": PROTEIN_UNIVERSE_SCHEMA,
            "holdout_id": holdout_id,
            "policy_id": policy_id,
            "species_id": species_id,
            "truth_access": False,
            "candidate_access": False,
            "fuzzy_mapping_used": False,
            "source_annotation_preserved": True,
            "usage": "protein_dependent_WGDI_only",
            "mapping_policy": (
                "exact_header_gene_or_transcript_or_GFF_CDS_protein_id; "
                "longest_sequence_then_protein_id"
            ),
            "duplicate_policy": (
                "collapse_only_same_normalized_ID_type_seqid_coordinates_strand_"
                "phase_and_parent_set;reject_gene_or_transcript_multilocus"
            ),
            "hierarchy": hierarchy,
            "counts": {
                **collapse_audit,
                **provider_audit,
                "normalized_genes": len(genes),
                "coding_genes": len(coding_genes),
                "genes_with_exact_provider_protein": len(selected_genes),
                "genes_excluded_without_exact_provider_protein": len(
                    coding_genes - selected_genes
                ),
                "retained_feature_rows": len(retained_records),
                "orphan_CDS_rows_without_Parent_excluded": orphan_cds_without_parent,
                "unresolved_CDS_rows_with_Parent": unresolved_cds_with_parent,
            },
            "inputs": {
                "primary_gff_sha256": sha256_file(source_gff),
                "provider_protein_sha256": sha256_file(source_protein),
            },
            "outputs": {
                "protein_supported_gff3_sha256": sha256_file(gff_output),
                "representative_protein_fasta_sha256": sha256_file(protein_output),
                "representatives_tsv_sha256": sha256_file(representatives_output),
                "excluded_genes_tsv_sha256": sha256_file(excluded_output),
            },
        }
        (working / "manifest.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        write_sha256sums(working)
        verify_sha256sums(working, ignore_checksum_file=True)
        os.replace(working, output)
        return manifest
    except BaseException:
        shutil.rmtree(working, ignore_errors=True)
        raise
