from __future__ import annotations

import csv
import hashlib
import json
import re
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from . import __version__
from .audit import PROTEIN_ALPHABET
from .io import (
    fasta_relation_id,
    iter_fasta,
    normalize_feature_id,
    open_text,
    parse_fasta_header_fields,
    read_fai,
)
from .perturb import GffRecord, read_gff_document


WGDI_INPUT_MANIFEST_SCHEMA = "ploidypatch.wgdi_input_manifest.v1"
SAFE_PREFIX = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")
CHROMOSOME_ALIAS = re.compile(
    r"^chr(?:omosome)?[-_]?([A-Za-z0-9]+)$", re.IGNORECASE
)


@dataclass(frozen=True)
class RepresentativeProtein:
    gene_id: str
    protein_id: str
    relation_source: str
    sequence: str


def _file_sha256(path: str | Path, block_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        while block := handle.read(block_size):
            digest.update(block)
    return digest.hexdigest()


def _write_exclusive(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8", newline="") as handle:
        handle.write(text)


def _normalized_gene_records(
    records: tuple[GffRecord, ...],
) -> tuple[dict[str, GffRecord], dict[str, str], str]:
    genes: dict[str, GffRecord] = {}
    transcript_to_gene: dict[str, str] = {}
    if not genes:
        for record in records:
            if (
                record.feature_type in {"mRNA", "transcript"}
                and record.feature_id
                and not record.parents
            ):
                transcript_id = normalize_feature_id(record.feature_id)
                if transcript_id in genes:
                    raise ValueError(
                        "Root transcript ID collision after normalization: "
                        f"{transcript_id}"
                    )
                genes[transcript_id] = record
                transcript_to_gene[transcript_id] = transcript_id
        if genes:
            return genes, transcript_to_gene, "root_transcript_as_gene"
    for record in records:
        if record.feature_type == "gene" and record.feature_id:
            gene_id = normalize_feature_id(record.feature_id)
            if gene_id in genes:
                raise ValueError(f"Gene ID collision after normalization: {gene_id}")
            genes[gene_id] = record
    for record in records:
        if record.feature_type not in {"mRNA", "transcript"} or not record.feature_id:
            continue
        gene_parents = {
            normalize_feature_id(value)
            for value in record.parents
            if normalize_feature_id(value) in genes
        }
        if len(gene_parents) != 1:
            continue
        transcript_id = normalize_feature_id(record.feature_id)
        gene_id = next(iter(gene_parents))
        previous = transcript_to_gene.get(transcript_id)
        if previous is not None and previous != gene_id:
            raise ValueError(
                f"Transcript maps to multiple normalized genes: {transcript_id}"
            )
        transcript_to_gene[transcript_id] = gene_id
    return genes, transcript_to_gene, "explicit_gene_parent"


def _select_representative_proteins(
    protein_path: str | Path,
    genes: dict[str, GffRecord],
    transcript_to_gene: dict[str, str],
) -> tuple[dict[str, RepresentativeProtein], dict[str, int]]:
    representatives: dict[str, RepresentativeProtein] = {}
    relation_sources: Counter[str] = Counter()
    protein_records = 0
    unmapped_records = 0
    empty_records = 0
    invalid_records = 0
    invalid_character_counts: Counter[str] = Counter()
    for protein_id, header, raw_sequence in iter_fasta(protein_path):
        protein_records += 1
        fields = parse_fasta_header_fields(header)
        explicit_gene = fields.get("gene") or fields.get("gene_id")
        if explicit_gene:
            gene_id = normalize_feature_id(explicit_gene)
            relation_source = "header:gene"
        else:
            transcript_id, transcript_source = fasta_relation_id(
                protein_id, header, relation="transcript"
            )
            gene_id = transcript_to_gene.get(transcript_id, "")
            relation_source = f"{transcript_source}->gff:Parent"
        if not gene_id or gene_id not in genes:
            unmapped_records += 1
            continue
        sequence = raw_sequence.upper().rstrip("*")
        if not sequence:
            empty_records += 1
            continue
        invalid_characters = set(sequence) - PROTEIN_ALPHABET
        if invalid_characters:
            invalid_records += 1
            invalid_character_counts.update(invalid_characters)
            continue
        relation_sources[relation_source] += 1
        candidate = RepresentativeProtein(
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
            representatives[gene_id] = candidate
    summary = {
        "protein_records": protein_records,
        "unmapped_protein_records": unmapped_records,
        "empty_protein_records": empty_records,
        "invalid_protein_records": invalid_records,
        "invalid_protein_character_kinds": len(invalid_character_counts),
        "genes_with_representative": len(representatives),
        "relation_source_kinds": len(relation_sources),
    }
    for source, count in sorted(relation_sources.items()):
        summary[f"relation_source:{source}"] = count
    for character, count in sorted(invalid_character_counts.items()):
        summary[f"invalid_protein_character:{character}"] = count
    return representatives, summary


def discover_primary_chromosome_seqids(
    gff_path: str | Path,
) -> tuple[frozenset[str], dict[str, str]]:
    labels: dict[str, str] = {}
    with open_text(gff_path) as handle:
        for raw_line in handle:
            if not raw_line or raw_line.startswith("#"):
                continue
            fields = raw_line.rstrip("\r\n").split("\t")
            if len(fields) != 9 or fields[2] not in {"region", "chromosome"}:
                continue
            attributes = {}
            for item in fields[8].split(";"):
                key, separator, value = item.partition("=")
                if separator and key and key not in attributes:
                    attributes[key] = value
            if fields[2] == "region":
                if attributes.get("genome") == "chromosome":
                    chromosome = attributes.get("chromosome") or attributes.get(
                        "Name"
                    )
                else:
                    aliases = [
                        alias.strip()
                        for alias in attributes.get("Alias", "").split(",")
                        if alias.strip()
                    ]
                    explicit_aliases = [
                        match.group(1)
                        for alias in aliases
                        if (match := CHROMOSOME_ALIAS.fullmatch(alias))
                    ]
                    if not explicit_aliases:
                        continue
                    normalized_explicit = {
                        alias.lstrip("0") or "0" for alias in explicit_aliases
                    }
                    chromosome = next(
                        (
                            alias
                            for alias in aliases
                            if alias.isdigit()
                            and (alias.lstrip("0") or "0") in normalized_explicit
                        ),
                        explicit_aliases[0],
                    )
            else:
                chromosome = attributes.get("Name") or attributes.get("chromosome")
                if not chromosome:
                    chromosome = attributes.get("ID", "").removeprefix(
                        "chromosome:"
                    )
                chromosome = chromosome or fields[0]
            if chromosome:
                labels[fields[0]] = chromosome
    if not labels:
        raise ValueError(
            "No primary chromosomes were identified from NCBI region records "
            "or Ensembl chromosome/explicit chromosome-alias records"
        )
    return frozenset(labels), labels


def prepare_wgdi_inputs(
    gff_path: str | Path,
    protein_path: str | Path,
    fai_path: str | Path,
    output_dir: str | Path,
    prefix: str,
    min_genes_per_seqid: int = 1,
    primary_chromosomes_only: bool = False,
) -> dict[str, Any]:
    """Create normalized gene-level WGDI inputs with auditable representatives."""

    if not SAFE_PREFIX.fullmatch(prefix):
        raise ValueError(
            "prefix must start with an alphanumeric character and contain only "
            "letters, digits, dots, underscores, or hyphens"
        )
    if min_genes_per_seqid < 1:
        raise ValueError("min_genes_per_seqid must be at least 1")
    output_dir = Path(output_dir)
    gff_output = output_dir / f"{prefix}.wgdi.gff"
    lens_output = output_dir / f"{prefix}.wgdi.lens"
    protein_output = output_dir / f"{prefix}.wgdi.pep.fa"
    representatives_output = output_dir / f"{prefix}.representatives.tsv"
    manifest_output = output_dir / f"{prefix}.wgdi_inputs.manifest.json"
    output_paths = (
        gff_output,
        lens_output,
        protein_output,
        representatives_output,
        manifest_output,
    )
    collisions = [path for path in output_paths if path.exists()]
    if collisions:
        raise FileExistsError(
            "Refusing to overwrite WGDI input artifact(s): "
            + ", ".join(str(path) for path in collisions)
        )

    chromosome_labels: dict[str, str] = {}
    selected_seqids: frozenset[str] | None = None
    if primary_chromosomes_only:
        selected_seqids, chromosome_labels = discover_primary_chromosome_seqids(
            gff_path
        )
    document = read_gff_document(gff_path, include_seqids=selected_seqids)
    sequence_lengths = read_fai(fai_path)
    selected_records = document.records
    genes, transcript_to_gene, gene_record_policy = _normalized_gene_records(
        selected_records
    )
    representatives, protein_summary = _select_representative_proteins(
        protein_path, genes, transcript_to_gene
    )

    genes_by_seqid: dict[str, list[tuple[str, GffRecord]]] = defaultdict(list)
    for gene_id, record in genes.items():
        if gene_id in representatives:
            genes_by_seqid[record.seqid].append((gene_id, record))
    retained_seqids = {
        seqid
        for seqid, values in genes_by_seqid.items()
        if len(values) >= min_genes_per_seqid
    }
    missing_lengths = retained_seqids - set(sequence_lengths)
    if missing_lengths:
        examples = ", ".join(sorted(missing_lengths)[:10])
        raise ValueError(f"Retained GFF3 sequence IDs absent from FAI: {examples}")

    retained_genes: list[tuple[str, GffRecord, int]] = []
    gff_lines: list[str] = []
    lens_lines: list[str] = []
    for seqid in sorted(retained_seqids):
        ordered = sorted(
            genes_by_seqid[seqid],
            key=lambda item: (item[1].start, item[1].end, item[0]),
        )
        lens_lines.append(f"{seqid}\t{sequence_lengths[seqid]}\t{len(ordered)}\n")
        for order, (gene_id, record) in enumerate(ordered, start=1):
            retained_genes.append((gene_id, record, order))
            gff_lines.append(
                f"{seqid}\t{gene_id}\t{record.start}\t{record.end}\t"
                f"{record.strand}\t{order}\n"
            )

    protein_lines: list[str] = []
    representative_lines = [
        "gene_id\tprotein_id\trelation_source\tprotein_length\tseqid\tgene_order\n"
    ]
    for gene_id, record, order in retained_genes:
        representative = representatives[gene_id]
        protein_lines.append(f">{gene_id}\n")
        sequence = representative.sequence
        protein_lines.extend(
            sequence[index : index + 60] + "\n"
            for index in range(0, len(sequence), 60)
        )
        representative_lines.append(
            f"{gene_id}\t{representative.protein_id}\t"
            f"{representative.relation_source}\t{len(sequence)}\t"
            f"{record.seqid}\t{order}\n"
        )

    _write_exclusive(gff_output, "".join(gff_lines))
    _write_exclusive(lens_output, "".join(lens_lines))
    _write_exclusive(protein_output, "".join(protein_lines))
    _write_exclusive(representatives_output, "".join(representative_lines))
    retained_gene_ids = {gene_id for gene_id, _, _ in retained_genes}
    manifest: dict[str, Any] = {
        "schema_version": WGDI_INPUT_MANIFEST_SCHEMA,
        "generator": {"name": "PloidyPatch", "version": __version__},
        "source": {
            "gff": {"file_name": Path(gff_path).name, "sha256": _file_sha256(gff_path)},
            "protein": {
                "file_name": Path(protein_path).name,
                "sha256": _file_sha256(protein_path),
            },
            "fai": {"file_name": Path(fai_path).name, "sha256": _file_sha256(fai_path)},
        },
        "parameters": {
            "prefix": prefix,
            "min_genes_per_seqid": min_genes_per_seqid,
            "primary_chromosomes_only": primary_chromosomes_only,
            "gene_record_policy": gene_record_policy,
        },
        "counts": {
            "input_genes": len(genes),
            **protein_summary,
            "retained_seqids": len(retained_seqids),
            "excluded_seqids": len(genes_by_seqid) - len(retained_seqids),
            "retained_genes": len(retained_gene_ids),
            "genes_excluded_by_seqid_threshold": len(representatives)
            - len(retained_gene_ids),
            "genes_without_representative": len(genes) - len(representatives),
        },
        "chromosome_labels": chromosome_labels,
        "outputs": {
            "gff": {"file_name": gff_output.name, "sha256": _file_sha256(gff_output)},
            "lens": {
                "file_name": lens_output.name,
                "sha256": _file_sha256(lens_output),
            },
            "protein": {
                "file_name": protein_output.name,
                "sha256": _file_sha256(protein_output),
            },
            "representatives": {
                "file_name": representatives_output.name,
                "sha256": _file_sha256(representatives_output),
            },
        },
    }
    _write_exclusive(
        manifest_output,
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
    )
    return manifest
