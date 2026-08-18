from __future__ import annotations

import csv
import hashlib
import json
import re
from collections import Counter
from contextlib import ExitStack
from pathlib import Path
from typing import Any, BinaryIO

from . import __version__
from .baseline import _parse_miniprot_models
from .io import iter_fasta, open_text
from .perturb import _file_sha256, read_gff_document
from .score import build_annotation_index


PAF_DELETION_SCHEMA_VERSION = "ploidypatch.paf_target_deletions.v2"
PAV_GENE_CATALOG_SCHEMA_VERSION = "ploidypatch.assembly_pav_gene_catalog.v2"
PAV_RECONCILED_CATALOG_SCHEMA_VERSION = (
    "ploidypatch.reconciled_assembly_pav_gene_catalog.v1"
)
COMBINED_PAF_SCHEMA_VERSION = "ploidypatch.combined_chromosome_paf.v1"
PAV_PROTEIN_SUBSET_SCHEMA_VERSION = "ploidypatch.pav_protein_subset.v1"
PAV_PROTEIN_SCREEN_SCHEMA_VERSION = "ploidypatch.pav_protein_screen.v3"
PAV_WGDI_ANNOTATION_SCHEMA_VERSION = "ploidypatch.pav_wgdi_annotation.v1"
PAV_EVENT_CATALOG_SCHEMA_VERSION = "ploidypatch.pav_event_catalog.v1"
CIGAR_PATTERN = re.compile(r"(\d+)([MIDNSHP=X])")
QUERY_CONSUMING = {"M", "I", "=", "X"}
TARGET_CONSUMING = {"M", "D", "N", "=", "X"}
ALIGNED = {"M", "=", "X"}


def _read_fai_records(
    path: str | Path,
) -> dict[str, tuple[int, int, int, int]]:
    """Read the fields required for random access to an uncompressed FASTA."""

    records: dict[str, tuple[int, int, int, int]] = {}
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            fields = raw_line.rstrip("\r\n").split("\t")
            if len(fields) < 5:
                raise ValueError(f"Malformed FAI line {line_number} in {path}")
            seqid = fields[0]
            if seqid in records:
                raise ValueError(f"Duplicate FAI sequence ID {seqid!r} in {path}")
            try:
                length, offset, line_bases, line_width = map(int, fields[1:5])
            except ValueError as exc:
                raise ValueError(
                    f"Non-integer FAI field on line {line_number} in {path}"
                ) from exc
            if (
                length < 0
                or offset < 0
                or line_bases < 1
                or line_width < line_bases
            ):
                raise ValueError(f"Invalid FAI geometry on line {line_number} in {path}")
            records[seqid] = (length, offset, line_bases, line_width)
    if not records:
        raise ValueError(f"FAI contains no records: {path}")
    return records


def _fetch_indexed_fasta(
    handle: BinaryIO,
    record: tuple[int, int, int, int],
    start: int,
    end: int,
) -> str:
    """Fetch a zero-based, half-open interval from an indexed FASTA."""

    length, offset, line_bases, line_width = record
    if start < 0 or end < start or end > length:
        raise ValueError(f"FASTA interval outside sequence: {start}-{end} / {length}")

    def byte_position(position: int) -> int:
        line, column = divmod(position, line_bases)
        return offset + line * line_width + column

    byte_start = byte_position(start)
    byte_end = byte_position(end)
    handle.seek(byte_start)
    raw = handle.read(byte_end - byte_start)
    sequence = raw.replace(b"\n", b"").replace(b"\r", b"").decode("ascii")
    if len(sequence) != end - start:
        raise ValueError(
            "FASTA/FAI random-access length mismatch: "
            f"expected {end - start}, observed {len(sequence)}"
        )
    return sequence


def _query_flank_metrics(
    *,
    handle: BinaryIO,
    fai_records: dict[str, tuple[int, int, int, int]],
    seqid: str,
    breakpoint: int,
    flank_bp: int,
    max_ambiguous_fraction: float,
) -> dict[str, str | int]:
    record = fai_records.get(seqid)
    if record is None:
        raise ValueError(f"Query sequence {seqid!r} is absent from the supplied FAI")
    sequence_length = record[0]
    if breakpoint < 0 or breakpoint > sequence_length:
        raise ValueError(
            f"Query breakpoint outside {seqid}: {breakpoint} / {sequence_length}"
        )
    start = max(0, breakpoint - flank_bp)
    end = min(sequence_length, breakpoint + flank_bp)
    sequence = _fetch_indexed_fasta(handle, record, start, end).upper()
    expected_bp = 2 * flank_bp
    observed_bp = len(sequence)
    ambiguous_bp = sum(base not in "ACGT" for base in sequence)
    ambiguous_fraction = ambiguous_bp / observed_bp if observed_bp else 0.0
    if observed_bp < expected_bp:
        status = "truncated_at_sequence_end"
    elif ambiguous_fraction > max_ambiguous_fraction:
        status = "ambiguous_sequence"
    else:
        status = "pass"
    return {
        "query_flank_start": start + 1 if observed_bp else "",
        "query_flank_end": end if observed_bp else "",
        "query_flank_expected_bp": expected_bp,
        "query_flank_observed_bp": observed_bp,
        "query_flank_ambiguous_bp": ambiguous_bp,
        "query_flank_ambiguous_fraction": f"{ambiguous_fraction:.8f}",
        "query_flank_status": status,
    }


def _parse_tags(fields: list[str]) -> dict[str, str]:
    tags = {}
    for field in fields:
        parts = field.split(":", 2)
        if len(parts) == 3:
            tags[parts[0]] = parts[2]
    return tags


def _parse_cigar(cigar: str) -> list[tuple[int, str]]:
    operations = [
        (int(length), operation)
        for length, operation in CIGAR_PATTERN.findall(cigar)
    ]
    if (
        not operations
        or "".join(
            f"{length}{operation}" for length, operation in operations
        )
        != cigar
    ):
        raise ValueError(f"Malformed PAF cg CIGAR: {cigar}")
    return operations


def combine_chromosome_pafs(
    *,
    input_dir_path: str | Path,
    chromosome_pair_tsv_path: str | Path,
    output_paf_path: str | Path,
) -> dict[str, Any]:
    """Combine per-query PAFs in declared chromosome order with validation."""

    input_dir = Path(input_dir_path)
    output_paf = Path(output_paf_path)
    manifest_path = Path(str(output_paf) + ".manifest.json")
    collisions = [path for path in (output_paf, manifest_path) if path.exists()]
    if collisions:
        raise FileExistsError(
            "Refusing to overwrite combined PAF artifact(s): "
            + ", ".join(str(path) for path in collisions)
        )
    pairs: list[tuple[str, str]] = []
    seen_queries: set[str] = set()
    with Path(chromosome_pair_tsv_path).open(
        "r", encoding="utf-8", newline=""
    ) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"query_seqid", "target_seqid"}
        if reader.fieldnames is None or not required <= set(reader.fieldnames):
            raise ValueError(
                "Chromosome-pair TSV requires query_seqid and target_seqid"
            )
        for line_number, row in enumerate(reader, start=2):
            query_seqid = row["query_seqid"]
            target_seqid = row["target_seqid"]
            if not query_seqid or not target_seqid:
                raise ValueError(
                    f"Incomplete chromosome pair at line {line_number}"
                )
            if query_seqid in seen_queries:
                raise ValueError(f"Duplicate query_seqid: {query_seqid}")
            seen_queries.add(query_seqid)
            pairs.append((query_seqid, target_seqid))
    if not pairs:
        raise ValueError("Chromosome-pair TSV contains no pairs")

    output_paf.parent.mkdir(parents=True, exist_ok=True)
    inputs: list[dict[str, Any]] = []
    counts: Counter[str] = Counter()
    with output_paf.open("xb") as destination:
        for query_seqid, expected_target in pairs:
            source_path = input_dir / f"{query_seqid}.paf"
            if not source_path.is_file() or source_path.stat().st_size == 0:
                raise ValueError(f"Missing or empty chromosome PAF: {source_path}")
            source_rows = 0
            with source_path.open("rb") as source:
                for line_number, raw_line in enumerate(source, start=1):
                    if not raw_line.strip() or raw_line.startswith(b"#"):
                        destination.write(raw_line)
                        continue
                    fields = raw_line.rstrip(b"\r\n").split(b"\t")
                    if len(fields) < 12:
                        raise ValueError(
                            f"Malformed PAF line {line_number} in {source_path}"
                        )
                    observed_query = fields[0].decode("utf-8")
                    observed_target = fields[5].decode("utf-8")
                    if observed_query != query_seqid:
                        raise ValueError(
                            f"Unexpected query {observed_query!r} in {source_path}"
                        )
                    source_rows += 1
                    counts["alignments_total"] += 1
                    if observed_target == expected_target:
                        counts["alignments_expected_target"] += 1
                    else:
                        counts["alignments_other_target"] += 1
                    destination.write(raw_line)
            if source_rows == 0:
                raise ValueError(f"Chromosome PAF contains no alignments: {source_path}")
            inputs.append(
                {
                    "query_seqid": query_seqid,
                    "expected_target_seqid": expected_target,
                    "file_name": source_path.name,
                    "rows": source_rows,
                    "sha256": _file_sha256(source_path),
                }
            )
    manifest = {
        "schema_version": COMBINED_PAF_SCHEMA_VERSION,
        "generator": {"name": "PloidyPatch", "version": __version__},
        "input": {
            "directory_name": input_dir.name,
            "chromosome_pairs": {
                "file_name": Path(chromosome_pair_tsv_path).name,
                "sha256": _file_sha256(chromosome_pair_tsv_path),
                "pairs": len(pairs),
            },
            "pafs": inputs,
        },
        "counts": dict(sorted(counts.items())),
        "output": {
            "file_name": output_paf.name,
            "rows": counts["alignments_total"],
            "sha256": _file_sha256(output_paf),
        },
    }
    with manifest_path.open("x", encoding="utf-8", newline="") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return manifest


def subset_catalog_proteins(
    *,
    catalog_tsv_path: str | Path,
    protein_fasta_path: str | Path,
    output_fasta_path: str | Path,
) -> dict[str, Any]:
    """Extract one gene-keyed representative protein per PAV catalog row."""

    output_fasta = Path(output_fasta_path)
    manifest_path = Path(str(output_fasta) + ".manifest.json")
    collisions = [path for path in (output_fasta, manifest_path) if path.exists()]
    if collisions:
        raise FileExistsError(
            "Refusing to overwrite PAV protein artifact(s): "
            + ", ".join(str(path) for path in collisions)
        )
    gene_ids: list[str] = []
    seen_catalog_genes: set[str] = set()
    with Path(catalog_tsv_path).open(
        "r", encoding="utf-8", newline=""
    ) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or "gene_id" not in reader.fieldnames:
            raise ValueError("PAV catalog lacks gene_id")
        for line_number, row in enumerate(reader, start=2):
            gene_id = row["gene_id"]
            if not gene_id or gene_id in seen_catalog_genes:
                raise ValueError(
                    f"Empty or duplicate catalog gene_id at line {line_number}"
                )
            seen_catalog_genes.add(gene_id)
            gene_ids.append(gene_id)
    if not gene_ids:
        raise ValueError("PAV catalog contains no genes")

    wanted = set(gene_ids)
    proteins: dict[str, str] = {}
    for record_id, _header, sequence in iter_fasta(protein_fasta_path):
        if record_id not in wanted:
            continue
        if record_id in proteins:
            raise ValueError(f"Duplicate candidate protein ID: {record_id}")
        protein = sequence.upper().rstrip("*")
        if not protein:
            raise ValueError(f"Empty candidate protein: {record_id}")
        proteins[record_id] = protein
    missing = wanted - set(proteins)
    if missing:
        raise ValueError(
            "Catalog gene(s) absent from representative protein FASTA: "
            + ", ".join(sorted(missing)[:10])
        )

    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    with output_fasta.open("x", encoding="utf-8", newline="") as handle:
        for gene_id in gene_ids:
            handle.write(f">{gene_id}\n")
            sequence = proteins[gene_id]
            for index in range(0, len(sequence), 70):
                handle.write(sequence[index : index + 70] + "\n")
    lengths = [len(proteins[gene_id]) for gene_id in gene_ids]
    manifest = {
        "schema_version": PAV_PROTEIN_SUBSET_SCHEMA_VERSION,
        "generator": {"name": "PloidyPatch", "version": __version__},
        "inputs": {
            "catalog": {
                "file_name": Path(catalog_tsv_path).name,
                "sha256": _file_sha256(catalog_tsv_path),
                "genes": len(gene_ids),
            },
            "representative_proteins": {
                "file_name": Path(protein_fasta_path).name,
                "sha256": _file_sha256(protein_fasta_path),
            },
        },
        "protein_lengths_aa": {
            "minimum": min(lengths),
            "maximum": max(lengths),
            "total": sum(lengths),
        },
        "output": {
            "file_name": output_fasta.name,
            "records": len(gene_ids),
            "sha256": _file_sha256(output_fasta),
        },
    }
    with manifest_path.open("x", encoding="utf-8", newline="") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return manifest


def screen_pav_candidate_proteins(
    *,
    catalog_tsv_path: str | Path,
    miniprot_gff_path: str | Path,
    protein_map_path: str | Path,
    chromosome_pair_tsv_path: str | Path,
    output_tsv_path: str | Path,
    min_identity: float = 0.5,
    min_query_coverage: float = 0.5,
    min_partial_query_coverage: float = 0.2,
    max_local_distance_bp: int = 100000,
    max_breakpoint_spread_bp: int = 1000,
) -> dict[str, Any]:
    """Classify protein evidence without conflating plant homeologs with loci."""

    if not 0 <= min_identity <= 1:
        raise ValueError("Minimum identity must be within [0, 1]")
    if not 0 <= min_partial_query_coverage <= min_query_coverage <= 1:
        raise ValueError(
            "Query coverage thresholds must satisfy 0 <= partial <= full <= 1"
        )
    if max_local_distance_bp < 0 or max_breakpoint_spread_bp < 0:
        raise ValueError("Breakpoint distance thresholds must be non-negative")
    output_tsv = Path(output_tsv_path)
    manifest_path = Path(str(output_tsv) + ".manifest.json")
    collisions = [path for path in (output_tsv, manifest_path) if path.exists()]
    if collisions:
        raise FileExistsError(
            "Refusing to overwrite PAV protein-screen artifact(s): "
            + ", ".join(str(path) for path in collisions)
        )

    catalog_fields = [
        "candidate_id",
        "gene_id",
        "seqid",
        "gene_start",
        "gene_end",
        "strand",
        "gene_span_bp",
        "transcript_structures",
        "cds_union_bp",
        "support_count",
        "support_sources",
    ]
    catalog_rows: list[dict[str, str]] = []
    catalog_input_fields: list[str] = []
    genes: set[str] = set()
    with Path(catalog_tsv_path).open(
        "r", encoding="utf-8", newline=""
    ) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        catalog_input_fields = list(reader.fieldnames or [])
        missing = set(catalog_fields) - set(catalog_input_fields)
        if reader.fieldnames is None or missing:
            raise ValueError(
                "PAV catalog is missing required column(s): "
                + ", ".join(sorted(missing))
            )
        for line_number, row in enumerate(reader, start=2):
            gene_id = row["gene_id"]
            if not gene_id or gene_id in genes:
                raise ValueError(
                    f"Empty or duplicate PAV gene_id at line {line_number}"
                )
            genes.add(gene_id)
            support_sources = [
                source for source in row["support_sources"].split(",") if source
            ]
            if not support_sources:
                raise ValueError(f"No support sources at catalog line {line_number}")
            breakpoint_fields = {
                f"{source}_{suffix}"
                for source in support_sources
                for suffix in ("query_seqid", "query_breakpoint")
            }
            missing_breakpoint_fields = breakpoint_fields - set(row)
            if missing_breakpoint_fields:
                raise ValueError(
                    f"PAV catalog line {line_number} lacks breakpoint field(s): "
                    + ", ".join(sorted(missing_breakpoint_fields))
                )
            catalog_rows.append(dict(row))
    if not catalog_rows:
        raise ValueError("PAV catalog contains no candidates")

    primary_query_seqids: set[str] = set()
    with Path(chromosome_pair_tsv_path).open(
        "r", encoding="utf-8", newline=""
    ) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"query_seqid", "target_seqid"}
        if reader.fieldnames is None or not required <= set(reader.fieldnames):
            raise ValueError(
                "Chromosome-pair TSV requires query_seqid and target_seqid"
            )
        for line_number, row in enumerate(reader, start=2):
            query_seqid = row["query_seqid"]
            if not query_seqid or query_seqid in primary_query_seqids:
                raise ValueError(
                    f"Empty or duplicate query_seqid at pair line {line_number}"
                )
            primary_query_seqids.add(query_seqid)
    if not primary_query_seqids:
        raise ValueError("Chromosome-pair TSV contains no primary query sequences")

    query_metadata: dict[str, tuple[str, int]] = {}
    query_by_gene: dict[str, str] = {}
    with Path(protein_map_path).open(
        "r", encoding="utf-8", newline=""
    ) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"query_id", "source", "source_record_id", "length_aa"}
        missing = required - set(reader.fieldnames or [])
        if reader.fieldnames is None or missing:
            raise ValueError(
                "Protein map is missing required column(s): "
                + ", ".join(sorted(missing))
            )
        for line_number, row in enumerate(reader, start=2):
            query_id = row["query_id"]
            gene_id = row["source_record_id"]
            try:
                length = int(row["length_aa"])
            except ValueError as exc:
                raise ValueError(
                    f"Invalid protein length in map line {line_number}"
                ) from exc
            if not query_id or length < 1 or query_id in query_metadata:
                raise ValueError(
                    f"Invalid or duplicate query_id in map line {line_number}"
                )
            query_metadata[query_id] = (row["source"], length)
            if gene_id in genes:
                if gene_id in query_by_gene:
                    raise ValueError(f"Multiple protein-map queries for {gene_id}")
                query_by_gene[gene_id] = query_id
    missing_queries = genes - set(query_by_gene)
    if missing_queries:
        raise ValueError(
            "PAV gene(s) absent from protein map: "
            + ", ".join(sorted(missing_queries)[:10])
        )

    models = _parse_miniprot_models(miniprot_gff_path, query_metadata)
    models_by_query: dict[str, list[Any]] = {}
    for model in models:
        models_by_query.setdefault(model.query_id, []).append(model)

    result_fields = [
        *catalog_input_fields,
        "query_id",
        "expected_query_seqid",
        "breakpoint_min",
        "breakpoint_max",
        "breakpoint_spread_bp",
        "models_total",
        "models_local_total",
        "models_local_substantial",
        "models_unlocalized_total",
        "models_unlocalized_substantial",
        "models_off_locus_total",
        "models_off_locus_substantial",
        "screen_class",
        "screen_decision",
        "best_evidence_scope",
        "best_model_id",
        "best_seqid",
        "best_start",
        "best_end",
        "best_strand",
        "best_rank",
        "best_score",
        "best_identity",
        "best_query_coverage",
        "best_frameshifts",
        "best_stop_codons",
    ]
    counts: Counter[str] = Counter()
    result_rows: list[dict[str, str | int]] = []
    for catalog_row in catalog_rows:
        query_id = query_by_gene[catalog_row["gene_id"]]
        query_models = models_by_query.get(query_id, [])
        support_sources = catalog_row["support_sources"].split(",")
        breakpoint_seqids = [
            catalog_row[f"{source}_query_seqid"] for source in support_sources
        ]
        try:
            breakpoints = [
                int(catalog_row[f"{source}_query_breakpoint"])
                for source in support_sources
            ]
        except ValueError as exc:
            raise ValueError(
                f"Non-integer query breakpoint for {catalog_row['gene_id']}"
            ) from exc
        expected_query_seqid = breakpoint_seqids[0]
        breakpoint_min = min(breakpoints)
        breakpoint_max = max(breakpoints)
        breakpoint_spread = breakpoint_max - breakpoint_min
        breakpoint = (breakpoint_min + breakpoint_max) // 2
        breakpoint_consistent = (
            len(set(breakpoint_seqids)) == 1
            and expected_query_seqid in primary_query_seqids
            and breakpoint_spread <= max_breakpoint_spread_bp
        )

        def tier(model: Any) -> int:
            if model.identity < min_identity:
                return 0
            if model.query_coverage >= min_query_coverage:
                return 3 if model.frameshifts == 0 and model.stop_codons == 0 else 2
            if model.query_coverage >= min_partial_query_coverage:
                return 1
            return 0

        def model_rank(model: Any) -> tuple[int, float, float, float, int]:
            return (
                tier(model),
                model.query_coverage,
                model.identity,
                model.score,
                -model.rank,
            )

        def best_model(candidates: list[Any]) -> Any | None:
            return max(candidates, key=model_rank, default=None)

        local_models: list[Any] = []
        unlocalized_models: list[Any] = []
        off_locus_models: list[Any] = []
        if breakpoint_consistent:
            for model in query_models:
                if model.seqid == expected_query_seqid:
                    distance = max(model.start - breakpoint, breakpoint - model.end, 0)
                    if distance <= max_local_distance_bp:
                        local_models.append(model)
                    else:
                        off_locus_models.append(model)
                elif model.seqid in primary_query_seqids:
                    off_locus_models.append(model)
                else:
                    unlocalized_models.append(model)
        local_substantial = [model for model in local_models if tier(model) > 0]
        unlocalized_substantial = [
            model for model in unlocalized_models if tier(model) > 0
        ]
        off_locus_substantial = [
            model for model in off_locus_models if tier(model) > 0
        ]

        if not breakpoint_consistent:
            screen_class = "breakpoint_disagreement"
            decision = "abstain_breakpoint_disagreement"
            evidence_scope = "undetermined"
            best = best_model(query_models)
        elif local_substantial:
            best = best_model(local_substantial)
            evidence_scope = "local"
            if tier(best) == 3:
                screen_class = "local_intact_homolog"
                decision = "exclude_local_sequence_present"
            elif tier(best) == 2:
                screen_class = "local_disrupted_homolog"
                decision = "abstain_local_disrupted_sequence"
            else:
                screen_class = "local_partial_homolog"
                decision = "abstain_local_partial_sequence"
        elif unlocalized_substantial:
            best = best_model(unlocalized_substantial)
            evidence_scope = "unlocalized"
            tier_label = {3: "intact", 2: "disrupted", 1: "partial"}[tier(best)]
            screen_class = f"unlocalized_{tier_label}_homolog"
            decision = "abstain_unlocalized_sequence_evidence"
        elif off_locus_substantial:
            best = best_model(off_locus_substantial)
            evidence_scope = "off_locus"
            tier_label = {3: "intact", 2: "disrupted", 1: "partial"}[tier(best)]
            screen_class = f"off_locus_{tier_label}_homolog_only"
            decision = "retain_sequence_absence_candidate"
        else:
            screen_class = "no_substantial_homolog"
            decision = "retain_sequence_absence_candidate"
            evidence_scope = "none"
            best = best_model(query_models)
        counts[f"screen_class_{screen_class}"] += 1
        counts[f"screen_decision_{decision}"] += 1
        result: dict[str, str | int] = {
            **{field: catalog_row[field] for field in catalog_input_fields},
            "query_id": query_id,
            "expected_query_seqid": expected_query_seqid,
            "breakpoint_min": breakpoint_min,
            "breakpoint_max": breakpoint_max,
            "breakpoint_spread_bp": breakpoint_spread,
            "models_total": len(query_models),
            "models_local_total": len(local_models),
            "models_local_substantial": len(local_substantial),
            "models_unlocalized_total": len(unlocalized_models),
            "models_unlocalized_substantial": len(unlocalized_substantial),
            "models_off_locus_total": len(off_locus_models),
            "models_off_locus_substantial": len(off_locus_substantial),
            "screen_class": screen_class,
            "screen_decision": decision,
            "best_evidence_scope": evidence_scope,
            "best_model_id": "",
            "best_seqid": "",
            "best_start": "",
            "best_end": "",
            "best_strand": "",
            "best_rank": "",
            "best_score": "",
            "best_identity": "",
            "best_query_coverage": "",
            "best_frameshifts": "",
            "best_stop_codons": "",
        }
        if best is not None:
            result.update(
                {
                    "best_model_id": best.model_id,
                    "best_seqid": best.seqid,
                    "best_start": best.start,
                    "best_end": best.end,
                    "best_strand": best.strand,
                    "best_rank": best.rank,
                    "best_score": f"{best.score:.8f}",
                    "best_identity": f"{best.identity:.8f}",
                    "best_query_coverage": f"{best.query_coverage:.8f}",
                    "best_frameshifts": best.frameshifts,
                    "best_stop_codons": best.stop_codons,
                }
            )
        result_rows.append(result)

    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    with output_tsv.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=result_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(result_rows)
    manifest = {
        "schema_version": PAV_PROTEIN_SCREEN_SCHEMA_VERSION,
        "generator": {"name": "PloidyPatch", "version": __version__},
        "inputs": {
            "catalog": {
                "file_name": Path(catalog_tsv_path).name,
                "sha256": _file_sha256(catalog_tsv_path),
                "candidates": len(catalog_rows),
            },
            "miniprot_gff": {
                "file_name": Path(miniprot_gff_path).name,
                "sha256": _file_sha256(miniprot_gff_path),
                "models": len(models),
            },
            "protein_map": {
                "file_name": Path(protein_map_path).name,
                "sha256": _file_sha256(protein_map_path),
            },
            "chromosome_pairs": {
                "file_name": Path(chromosome_pair_tsv_path).name,
                "sha256": _file_sha256(chromosome_pair_tsv_path),
                "primary_query_sequences": len(primary_query_seqids),
            },
        },
        "parameters": {
            "min_identity": min_identity,
            "min_query_coverage": min_query_coverage,
            "min_partial_query_coverage": min_partial_query_coverage,
            "intact_requires_zero_frameshifts_and_stop_codons": True,
            "max_local_distance_bp": max_local_distance_bp,
            "max_breakpoint_spread_bp": max_breakpoint_spread_bp,
            "off_locus_homologs_do_not_negate_copy_specific_absence": True,
            "unlocalized_homologs_force_abstention": True,
        },
        "counts": dict(sorted(counts.items())),
        "output": {
            "file_name": output_tsv.name,
            "rows": len(result_rows),
            "sha256": _file_sha256(output_tsv),
        },
    }
    with manifest_path.open("x", encoding="utf-8", newline="") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return manifest


def annotate_pav_with_wgdi(
    *,
    protein_screen_tsv_path: str | Path,
    wgdi_gene_tsv_path: str | Path,
    output_tsv_path: str | Path,
) -> dict[str, Any]:
    """Combine locus-aware absence screening with plant WGD evidence."""

    output_tsv = Path(output_tsv_path)
    manifest_path = Path(str(output_tsv) + ".manifest.json")
    collisions = [path for path in (output_tsv, manifest_path) if path.exists()]
    if collisions:
        raise FileExistsError(
            "Refusing to overwrite PAV/WGDI artifact(s): "
            + ", ".join(str(path) for path in collisions)
        )

    screen_rows: list[dict[str, str]] = []
    screen_genes: set[str] = set()
    with Path(protein_screen_tsv_path).open(
        "r", encoding="utf-8", newline=""
    ) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        screen_fields = list(reader.fieldnames or [])
        required = {"gene_id", "screen_class", "screen_decision"}
        missing = required - set(screen_fields)
        if not screen_fields or missing:
            raise ValueError(
                "Protein screen is missing required column(s): "
                + ", ".join(sorted(missing))
            )
        for line_number, row in enumerate(reader, start=2):
            gene_id = row["gene_id"]
            if not gene_id or gene_id in screen_genes:
                raise ValueError(
                    f"Empty or duplicate screen gene_id at line {line_number}"
                )
            screen_genes.add(gene_id)
            screen_rows.append(row)
    if not screen_rows:
        raise ValueError("Protein screen contains no candidates")

    wgdi_by_gene: dict[str, dict[str, str]] = {}
    with Path(wgdi_gene_tsv_path).open(
        "r", encoding="utf-8", newline=""
    ) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        wgdi_fields = list(reader.fieldnames or [])
        required = {
            "gene_id",
            "subgenome",
            "expected_progenitor",
            "expected_progenitor_supported",
            "cross_progenitor_supported",
            "synteny_stratum",
        }
        missing = required - set(wgdi_fields)
        if not wgdi_fields or missing:
            raise ValueError(
                "WGDI evidence is missing required column(s): "
                + ", ".join(sorted(missing))
            )
        for line_number, row in enumerate(reader, start=2):
            gene_id = row["gene_id"]
            if not gene_id or gene_id in wgdi_by_gene:
                raise ValueError(
                    f"Empty or duplicate WGDI gene_id at line {line_number}"
                )
            wgdi_by_gene[gene_id] = row
    missing_wgdi = screen_genes - set(wgdi_by_gene)
    if missing_wgdi:
        raise ValueError(
            "Protein-screen gene(s) absent from WGDI evidence: "
            + ", ".join(sorted(missing_wgdi)[:10])
        )

    appended_wgdi_fields = [field for field in wgdi_fields if field != "gene_id"]
    fieldnames = [
        *screen_fields,
        *(f"wgdi_{field}" for field in appended_wgdi_fields),
        "plant_pav_evidence_class",
        "plant_pav_decision",
    ]
    output_rows: list[dict[str, str]] = []
    counts: Counter[str] = Counter()
    for screen_row in screen_rows:
        wgdi_row = wgdi_by_gene[screen_row["gene_id"]]
        screen_decision = screen_row["screen_decision"]
        expected_supported = wgdi_row["expected_progenitor_supported"].lower()
        cross_supported = wgdi_row["cross_progenitor_supported"].lower()
        if expected_supported not in {"true", "false"} or cross_supported not in {
            "true",
            "false",
        }:
            raise ValueError(
                f"Invalid WGDI boolean for {screen_row['gene_id']}"
            )
        if screen_decision == "retain_sequence_absence_candidate":
            if expected_supported == "true":
                evidence_class = "wgd_anchored_sequence_absence_candidate"
                final_decision = "retain_high_confidence_candidate"
            elif cross_supported == "true":
                evidence_class = "cross_only_synteny_uncertain"
                final_decision = "abstain_cross_only_synteny"
            else:
                evidence_class = "no_synteny_support_uncertain"
                final_decision = "abstain_no_synteny_support"
        elif screen_decision.startswith("exclude_"):
            evidence_class = "sequence_presence_contradiction"
            final_decision = screen_decision
        else:
            evidence_class = "sequence_evidence_uncertain"
            final_decision = screen_decision
        counts[f"plant_pav_evidence_class_{evidence_class}"] += 1
        counts[f"plant_pav_decision_{final_decision}"] += 1
        output_rows.append(
            {
                **screen_row,
                **{
                    f"wgdi_{field}": wgdi_row[field]
                    for field in appended_wgdi_fields
                },
                "plant_pav_evidence_class": evidence_class,
                "plant_pav_decision": final_decision,
            }
        )

    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    with output_tsv.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(output_rows)
    manifest = {
        "schema_version": PAV_WGDI_ANNOTATION_SCHEMA_VERSION,
        "generator": {"name": "PloidyPatch", "version": __version__},
        "inputs": {
            "protein_screen": {
                "file_name": Path(protein_screen_tsv_path).name,
                "sha256": _file_sha256(protein_screen_tsv_path),
                "rows": len(screen_rows),
            },
            "wgdi_gene_evidence": {
                "file_name": Path(wgdi_gene_tsv_path).name,
                "sha256": _file_sha256(wgdi_gene_tsv_path),
                "rows": len(wgdi_by_gene),
            },
        },
        "decision_policy": {
            "high_confidence_requires": [
                "locus-aware protein screen retains sequence-absence candidate",
                "expected progenitor has WGDI collinearity support",
            ],
            "high_confidence_is_not_biological_truth": True,
            "cross_only_or_no_synteny_support_forces_abstention": True,
        },
        "counts": dict(sorted(counts.items())),
        "output": {
            "file_name": output_tsv.name,
            "rows": len(output_rows),
            "sha256": _file_sha256(output_tsv),
        },
    }
    with manifest_path.open("x", encoding="utf-8", newline="") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return manifest


def collapse_high_confidence_pav_events(
    *,
    annotated_gene_tsv_path: str | Path,
    output_tsv_path: str | Path,
) -> dict[str, Any]:
    """Collapse high-confidence gene candidates into structural absence events."""

    output_tsv = Path(output_tsv_path)
    manifest_path = Path(str(output_tsv) + ".manifest.json")
    collisions = [path for path in (output_tsv, manifest_path) if path.exists()]
    if collisions:
        raise FileExistsError(
            "Refusing to overwrite PAV event artifact(s): "
            + ", ".join(str(path) for path in collisions)
        )
    required = {
        "gene_id",
        "gene_start",
        "seqid",
        "support_sources",
        "screen_class",
        "plant_pav_decision",
        "wgdi_subgenome",
        "wgdi_chromosome_label",
        "wgdi_synteny_stratum",
    }
    groups: dict[tuple[tuple[str, ...], ...], list[dict[str, str]]] = {}
    input_rows = 0
    retained_gene_ids: set[str] = set()
    with Path(annotated_gene_tsv_path).open(
        "r", encoding="utf-8", newline=""
    ) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        field_set = set(reader.fieldnames or [])
        missing = required - field_set
        if reader.fieldnames is None or missing:
            raise ValueError(
                "Annotated PAV catalog is missing required column(s): "
                + ", ".join(sorted(missing))
            )
        for line_number, row in enumerate(reader, start=2):
            input_rows += 1
            if row["plant_pav_decision"] != "retain_high_confidence_candidate":
                continue
            gene_id = row["gene_id"]
            if not gene_id or gene_id in retained_gene_ids:
                raise ValueError(
                    f"Empty or duplicate retained gene_id at line {line_number}"
                )
            retained_gene_ids.add(gene_id)
            sources = [source for source in row["support_sources"].split(",") if source]
            if not sources:
                raise ValueError(f"No structural sources at line {line_number}")
            source_fields = {
                f"{source}_{suffix}"
                for source in sources
                for suffix in (
                    "deletion_id",
                    "deletion_start",
                    "deletion_end",
                    "deletion_bp",
                    "query_seqid",
                    "query_breakpoint",
                )
            }
            missing_source_fields = source_fields - field_set
            if missing_source_fields:
                raise ValueError(
                    f"Annotated PAV line {line_number} lacks structural field(s): "
                    + ", ".join(sorted(missing_source_fields))
                )
            key = tuple(
                (
                    source,
                    row["seqid"],
                    row[f"{source}_deletion_start"],
                    row[f"{source}_deletion_end"],
                    row[f"{source}_query_seqid"],
                    row[f"{source}_query_breakpoint"],
                )
                for source in sources
            )
            groups.setdefault(key, []).append(row)
    if not retained_gene_ids:
        raise ValueError("Annotated PAV catalog has no high-confidence candidates")

    fieldnames = [
        "event_id",
        "event_status",
        "seqid",
        "chromosome_label",
        "subgenome",
        "target_start_min",
        "target_start_max",
        "target_end_min",
        "target_end_max",
        "target_boundary_spread_bp",
        "deletion_bp_min",
        "deletion_bp_max",
        "query_seqid",
        "query_breakpoint_min",
        "query_breakpoint_max",
        "query_breakpoint_spread_bp",
        "source_count",
        "support_sources",
        "source_deletion_ids",
        "gene_count",
        "gene_ids",
        "screen_classes",
        "wgdi_synteny_strata",
    ]
    event_rows: list[dict[str, str | int]] = []
    counts: Counter[str] = Counter()
    for key, rows in groups.items():
        sources = [entry[0] for entry in key]
        seqids = {row["seqid"] for row in rows}
        subgenomes = {row["wgdi_subgenome"] for row in rows}
        chromosome_labels = {row["wgdi_chromosome_label"] for row in rows}
        query_seqids = {
            row[f"{source}_query_seqid"] for row in rows for source in sources
        }
        if (
            len(seqids) != 1
            or len(subgenomes) != 1
            or len(chromosome_labels) != 1
            or len(query_seqids) != 1
        ):
            raise ValueError("A structural PAV event crosses incompatible contexts")
        starts = [int(row[f"{source}_deletion_start"]) for source in sources for row in rows]
        ends = [int(row[f"{source}_deletion_end"]) for source in sources for row in rows]
        deletion_sizes = [
            int(row[f"{source}_deletion_bp"]) for source in sources for row in rows
        ]
        breakpoints = [
            int(row[f"{source}_query_breakpoint"])
            for source in sources
            for row in rows
        ]
        gene_ids = sorted(
            (row["gene_id"] for row in rows),
            key=lambda gene_id: next(
                int(row["gene_start"]) for row in rows if row["gene_id"] == gene_id
            ),
        )
        source_deletion_ids = [
            f"{source}:{rows[0][f'{source}_deletion_id']}" for source in sources
        ]
        digest_material = "\t".join(
            [next(iter(seqids)), *["|".join(entry) for entry in key]]
        )
        event_id = "PPPAEV-" + hashlib.sha256(
            digest_material.encode("utf-8")
        ).hexdigest()[:20]
        event_rows.append(
            {
                "event_id": event_id,
                "event_status": "high_confidence_sequence_absence_candidate",
                "seqid": next(iter(seqids)),
                "chromosome_label": next(iter(chromosome_labels)),
                "subgenome": next(iter(subgenomes)),
                "target_start_min": min(starts),
                "target_start_max": max(starts),
                "target_end_min": min(ends),
                "target_end_max": max(ends),
                "target_boundary_spread_bp": (max(starts) - min(starts))
                + (max(ends) - min(ends)),
                "deletion_bp_min": min(deletion_sizes),
                "deletion_bp_max": max(deletion_sizes),
                "query_seqid": next(iter(query_seqids)),
                "query_breakpoint_min": min(breakpoints),
                "query_breakpoint_max": max(breakpoints),
                "query_breakpoint_spread_bp": max(breakpoints) - min(breakpoints),
                "source_count": len(sources),
                "support_sources": ",".join(sources),
                "source_deletion_ids": ",".join(source_deletion_ids),
                "gene_count": len(gene_ids),
                "gene_ids": ",".join(gene_ids),
                "screen_classes": ",".join(
                    sorted({row["screen_class"] for row in rows})
                ),
                "wgdi_synteny_strata": ",".join(
                    sorted({row["wgdi_synteny_stratum"] for row in rows})
                ),
            }
        )
        counts["events_total"] += 1
        counts["genes_total"] += len(gene_ids)
        counts[f"events_with_{len(gene_ids)}_genes"] += 1
        counts[f"events_subgenome_{next(iter(subgenomes))}"] += 1
    event_rows.sort(
        key=lambda row: (row["seqid"], int(row["target_start_min"]), row["event_id"])
    )
    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    with output_tsv.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(event_rows)
    manifest = {
        "schema_version": PAV_EVENT_CATALOG_SCHEMA_VERSION,
        "generator": {"name": "PloidyPatch", "version": __version__},
        "input": {
            "file_name": Path(annotated_gene_tsv_path).name,
            "sha256": _file_sha256(annotated_gene_tsv_path),
            "rows": input_rows,
            "high_confidence_genes": len(retained_gene_ids),
        },
        "policy": {
            "input_decision_required": "retain_high_confidence_candidate",
            "grouping_key": (
                "ordered source-specific deletion boundaries and query breakpoints"
            ),
            "event_status_is_not_biological_truth": True,
        },
        "counts": dict(sorted(counts.items())),
        "output": {
            "file_name": output_tsv.name,
            "rows": len(event_rows),
            "sha256": _file_sha256(output_tsv),
        },
    }
    with manifest_path.open("x", encoding="utf-8", newline="") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return manifest


def extract_paf_target_deletions(
    *,
    paf_path: str | Path,
    output_tsv_path: str | Path,
    min_deletion_bp: int = 100,
    min_flanking_anchor_bp: int = 1000,
    min_mapq: int = 20,
    require_primary: bool = True,
    chromosome_pair_tsv_path: str | Path | None = None,
    query_genome_path: str | Path | None = None,
    query_fai_path: str | Path | None = None,
    query_flank_bp: int = 5000,
    max_query_flank_ambiguous_fraction: float = 0.0,
    require_clean_query_flanks: bool = True,
) -> dict[str, Any]:
    """Extract target-only intervals from assembly-alignment CIGAR strings."""

    if min_deletion_bp < 1 or min_flanking_anchor_bp < 0 or min_mapq < 0:
        raise ValueError("PAF deletion thresholds must be non-negative")
    if query_flank_bp < 0:
        raise ValueError("Query flank size must be non-negative")
    if not 0 <= max_query_flank_ambiguous_fraction <= 1:
        raise ValueError("Query-flank ambiguous fraction must be within [0, 1]")
    if (query_genome_path is None) != (query_fai_path is None):
        raise ValueError("Query genome and FAI must be supplied together")
    output_tsv = Path(output_tsv_path)
    manifest_path = Path(str(output_tsv) + ".manifest.json")
    collisions = [path for path in (output_tsv, manifest_path) if path.exists()]
    if collisions:
        raise FileExistsError(
            "Refusing to overwrite PAF-deletion artifact(s): "
            + ", ".join(str(path) for path in collisions)
        )
    output_tsv.parent.mkdir(parents=True, exist_ok=True)

    expected_targets: dict[str, str] | None = None
    if chromosome_pair_tsv_path is not None:
        expected_targets = {}
        with Path(chromosome_pair_tsv_path).open(
            "r", encoding="utf-8", newline=""
        ) as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            required = {"query_seqid", "target_seqid"}
            if reader.fieldnames is None or not required <= set(reader.fieldnames):
                raise ValueError(
                    "Chromosome-pair TSV requires query_seqid and target_seqid"
                )
            for line_number, row in enumerate(reader, start=2):
                query_seqid = row["query_seqid"]
                if query_seqid in expected_targets:
                    raise ValueError(
                        "Duplicate query_seqid in chromosome-pair TSV at line "
                        f"{line_number}: {query_seqid}"
                    )
                expected_targets[query_seqid] = row["target_seqid"]
        if not expected_targets or any(
            not query_seqid or not target_seqid
            for query_seqid, target_seqid in expected_targets.items()
        ):
            raise ValueError("Chromosome-pair TSV contains no complete pairs")

    fieldnames = [
        "deletion_id",
        "target_seqid",
        "target_start",
        "target_end",
        "deletion_bp",
        "query_seqid",
        "query_breakpoint",
        "strand",
        "mapq",
        "alignment_identity",
        "left_anchor_bp",
        "right_anchor_bp",
        "query_flank_start",
        "query_flank_end",
        "query_flank_expected_bp",
        "query_flank_observed_bp",
        "query_flank_ambiguous_bp",
        "query_flank_ambiguous_fraction",
        "query_flank_status",
        "paf_line",
    ]
    counts: Counter[str] = Counter()
    output_rows = 0
    fai_records = _read_fai_records(query_fai_path) if query_fai_path else None
    with ExitStack() as stack:
        source = stack.enter_context(open_text(paf_path))
        destination = stack.enter_context(
            output_tsv.open("x", encoding="utf-8", newline="")
        )
        query_handle = (
            stack.enter_context(Path(query_genome_path).open("rb"))
            if query_genome_path
            else None
        )
        writer = csv.DictWriter(destination, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for line_number, raw_line in enumerate(source, start=1):
            if not raw_line.strip() or raw_line.startswith("#"):
                continue
            counts["alignments_total"] += 1
            fields = raw_line.rstrip("\n").split("\t")
            if len(fields) < 12:
                raise ValueError(f"Malformed PAF line {line_number}: fewer than 12 fields")
            query_name = fields[0]
            query_start = int(fields[2])
            query_end = int(fields[3])
            strand = fields[4]
            target_name = fields[5]
            target_start = int(fields[7])
            target_end = int(fields[8])
            matching_bases = int(fields[9])
            alignment_block_length = int(fields[10])
            mapq = int(fields[11])
            tags = _parse_tags(fields[12:])
            if (
                expected_targets is not None
                and expected_targets.get(query_name) != target_name
            ):
                counts["alignments_outside_expected_chromosome_pair"] += 1
                continue
            if mapq < min_mapq:
                counts["alignments_mapq_below_threshold"] += 1
                continue
            if require_primary and tags.get("tp") != "P":
                counts["alignments_non_primary"] += 1
                continue
            cigar = tags.get("cg")
            if cigar is None:
                counts["alignments_missing_cigar"] += 1
                continue
            operations = _parse_cigar(cigar)
            query_consumed = sum(
                length for length, operation in operations if operation in QUERY_CONSUMING
            )
            target_consumed = sum(
                length for length, operation in operations if operation in TARGET_CONSUMING
            )
            if query_consumed != query_end - query_start:
                raise ValueError(
                    f"PAF query span/CIGAR mismatch at line {line_number}: "
                    f"{query_end - query_start} != {query_consumed}"
                )
            if target_consumed != target_end - target_start:
                raise ValueError(
                    f"PAF target span/CIGAR mismatch at line {line_number}: "
                    f"{target_end - target_start} != {target_consumed}"
                )
            counts["alignments_passing"] += 1
            query_offset = 0
            target_offset = 0
            alignment_digest = hashlib.sha256(
                raw_line.rstrip("\n").encode("utf-8")
            ).hexdigest()[:20]
            identity = (
                matching_bases / alignment_block_length
                if alignment_block_length
                else 0.0
            )
            for operation_index, (length, operation) in enumerate(operations):
                if operation == "D":
                    counts["target_deletions_total"] += 1
                    left_anchor = (
                        operations[operation_index - 1][0]
                        if operation_index > 0
                        and operations[operation_index - 1][1] in ALIGNED
                        else 0
                    )
                    right_anchor = (
                        operations[operation_index + 1][0]
                        if operation_index + 1 < len(operations)
                        and operations[operation_index + 1][1] in ALIGNED
                        else 0
                    )
                    if length < min_deletion_bp:
                        counts["target_deletions_below_size"] += 1
                    elif min(left_anchor, right_anchor) < min_flanking_anchor_bp:
                        counts["target_deletions_below_anchor"] += 1
                    else:
                        genomic_query_breakpoint = (
                            query_start + query_offset
                            if strand == "+"
                            else query_end - query_offset
                        )
                        flank_metrics: dict[str, str | int] = {
                            "query_flank_start": "",
                            "query_flank_end": "",
                            "query_flank_expected_bp": "",
                            "query_flank_observed_bp": "",
                            "query_flank_ambiguous_bp": "",
                            "query_flank_ambiguous_fraction": "",
                            "query_flank_status": "not_evaluated",
                        }
                        if query_handle is not None and fai_records is not None:
                            flank_metrics = _query_flank_metrics(
                                handle=query_handle,
                                fai_records=fai_records,
                                seqid=query_name,
                                breakpoint=genomic_query_breakpoint,
                                flank_bp=query_flank_bp,
                                max_ambiguous_fraction=(
                                    max_query_flank_ambiguous_fraction
                                ),
                            )
                            status = str(flank_metrics["query_flank_status"])
                            counts[f"target_deletions_query_flank_{status}"] += 1
                            if require_clean_query_flanks and status != "pass":
                                counts["target_deletions_excluded_by_query_flank"] += 1
                                if operation in QUERY_CONSUMING:
                                    query_offset += length
                                if operation in TARGET_CONSUMING:
                                    target_offset += length
                                continue
                        target_interval_start = target_start + target_offset + 1
                        target_interval_end = target_start + target_offset + length
                        deletion_id = (
                            f"PPD-{alignment_digest}-{operation_index:05d}"
                        )
                        writer.writerow(
                            {
                                "deletion_id": deletion_id,
                                "target_seqid": target_name,
                                "target_start": target_interval_start,
                                "target_end": target_interval_end,
                                "deletion_bp": length,
                                "query_seqid": query_name,
                                "query_breakpoint": genomic_query_breakpoint,
                                "strand": strand,
                                "mapq": mapq,
                                "alignment_identity": f"{identity:.8f}",
                                "left_anchor_bp": left_anchor,
                                "right_anchor_bp": right_anchor,
                                **flank_metrics,
                                "paf_line": line_number,
                            }
                        )
                        output_rows += 1
                if operation in QUERY_CONSUMING:
                    query_offset += length
                if operation in TARGET_CONSUMING:
                    target_offset += length

    manifest = {
        "schema_version": PAF_DELETION_SCHEMA_VERSION,
        "generator": {"name": "PloidyPatch", "version": __version__},
        "input": {
            "file_name": Path(paf_path).name,
            "sha256": _file_sha256(paf_path),
        },
        "parameters": {
            "min_deletion_bp": min_deletion_bp,
            "min_flanking_anchor_bp": min_flanking_anchor_bp,
            "min_mapq": min_mapq,
            "require_primary": require_primary,
            "query_flank_bp": query_flank_bp,
            "max_query_flank_ambiguous_fraction": (
                max_query_flank_ambiguous_fraction
            ),
            "require_clean_query_flanks": require_clean_query_flanks,
        },
        "counts": dict(sorted(counts.items())),
        "output": {
            "file_name": output_tsv.name,
            "rows": output_rows,
            "sha256": _file_sha256(output_tsv),
        },
    }
    if chromosome_pair_tsv_path is not None:
        manifest["input"]["chromosome_pairs"] = {
            "file_name": Path(chromosome_pair_tsv_path).name,
            "sha256": _file_sha256(chromosome_pair_tsv_path),
            "pairs": len(expected_targets or {}),
        }
    if query_genome_path is not None and query_fai_path is not None:
        manifest["input"]["query_genome"] = {
            "file_name": Path(query_genome_path).name,
            "sha256": _file_sha256(query_genome_path),
        }
        manifest["input"]["query_fai"] = {
            "file_name": Path(query_fai_path).name,
            "sha256": _file_sha256(query_fai_path),
            "sequences": len(fai_records or {}),
        }
    with manifest_path.open("x", encoding="utf-8", newline="") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return manifest


def _merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    merged: list[list[int]] = []
    for start, end in sorted(intervals):
        if merged and start <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return [(start, end) for start, end in merged]


def _interval_overlap_bp(
    intervals: list[tuple[int, int]], target_start: int, target_end: int
) -> int:
    return sum(
        max(0, min(end, target_end) - max(start, target_start) + 1)
        for start, end in intervals
    )


def catalog_genes_in_target_deletions(
    *,
    gff_path: str | Path,
    deletion_tsv_path: str | Path,
    output_tsv_path: str | Path,
    min_gene_coverage: float = 0.95,
    min_cds_coverage: float = 1.0,
) -> dict[str, Any]:
    """Catalog coding genes contained in one well-anchored target deletion."""

    if not 0 <= min_gene_coverage <= 1 or not 0 <= min_cds_coverage <= 1:
        raise ValueError("Gene and CDS coverage thresholds must be within [0, 1]")
    output_tsv = Path(output_tsv_path)
    manifest_path = Path(str(output_tsv) + ".manifest.json")
    collisions = [path for path in (output_tsv, manifest_path) if path.exists()]
    if collisions:
        raise FileExistsError(
            "Refusing to overwrite PAV-catalog artifact(s): "
            + ", ".join(str(path) for path in collisions)
        )

    deletions_by_seqid: dict[str, list[dict[str, str]]] = {}
    with Path(deletion_tsv_path).open(
        "r", encoding="utf-8", newline=""
    ) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "deletion_id",
            "target_seqid",
            "target_start",
            "target_end",
            "deletion_bp",
            "query_seqid",
            "query_breakpoint",
            "strand",
            "mapq",
            "alignment_identity",
            "left_anchor_bp",
            "right_anchor_bp",
        }
        if reader.fieldnames is None or not required <= set(reader.fieldnames):
            raise ValueError(
                "Deletion TSV is missing required column(s): "
                + ", ".join(sorted(required - set(reader.fieldnames or [])))
            )
        for row in reader:
            deletions_by_seqid.setdefault(row["target_seqid"], []).append(row)
    for rows in deletions_by_seqid.values():
        rows.sort(key=lambda row: (int(row["target_start"]), int(row["target_end"])))

    document = read_gff_document(gff_path)
    annotation = build_annotation_index(document)
    gene_records = {
        record.feature_id: record
        for record in document.records
        if record.feature_type == "gene" and record.feature_id
    }
    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "candidate_id",
        "gene_id",
        "seqid",
        "gene_start",
        "gene_end",
        "strand",
        "gene_span_bp",
        "transcript_structures",
        "cds_union_bp",
        "deletion_id",
        "deletion_start",
        "deletion_end",
        "deletion_bp",
        "gene_coverage",
        "cds_coverage",
        "query_seqid",
        "query_breakpoint",
        "alignment_strand",
        "mapq",
        "alignment_identity",
        "left_anchor_bp",
        "right_anchor_bp",
        "query_flank_start",
        "query_flank_end",
        "query_flank_expected_bp",
        "query_flank_observed_bp",
        "query_flank_ambiguous_bp",
        "query_flank_ambiguous_fraction",
        "query_flank_status",
    ]
    counts: Counter[str] = Counter()
    candidate_rows = []
    for gene_id, gene in gene_records.items():
        counts["genes_total"] += 1
        signatures = annotation.gene_signatures.get(gene_id, frozenset())
        cds_intervals = _merge_intervals(
            [
                (start, end)
                for signature in signatures
                for start, end, _ in signature.cds
            ]
        )
        if not cds_intervals:
            counts["genes_without_cds"] += 1
            continue
        counts["coding_genes"] += 1
        gene_span_bp = gene.end - gene.start + 1
        cds_union_bp = sum(end - start + 1 for start, end in cds_intervals)
        best = None
        for deletion in deletions_by_seqid.get(gene.seqid, []):
            deletion_start = int(deletion["target_start"])
            deletion_end = int(deletion["target_end"])
            if deletion_start > gene.end:
                break
            if deletion_end < gene.start:
                continue
            gene_overlap = max(
                0, min(gene.end, deletion_end) - max(gene.start, deletion_start) + 1
            )
            gene_coverage = gene_overlap / gene_span_bp
            cds_coverage = (
                _interval_overlap_bp(cds_intervals, deletion_start, deletion_end)
                / cds_union_bp
            )
            rank = (cds_coverage, gene_coverage, int(deletion["deletion_bp"]))
            if best is None or rank > best[0]:
                best = (rank, deletion, gene_coverage, cds_coverage)
        if best is None:
            counts["coding_genes_without_deletion_overlap"] += 1
            continue
        counts["coding_genes_overlapping_deletion"] += 1
        _, deletion, gene_coverage, cds_coverage = best
        if gene_coverage < min_gene_coverage:
            counts["genes_below_gene_coverage"] += 1
            continue
        if cds_coverage < min_cds_coverage:
            counts["genes_below_cds_coverage"] += 1
            continue
        candidate_digest = hashlib.sha256(
            f"{gene_id}\t{deletion['deletion_id']}".encode("utf-8")
        ).hexdigest()[:20]
        candidate_rows.append(
            {
                "candidate_id": f"PPPAV-{candidate_digest}",
                "gene_id": gene_id,
                "seqid": gene.seqid,
                "gene_start": gene.start,
                "gene_end": gene.end,
                "strand": gene.strand,
                "gene_span_bp": gene_span_bp,
                "transcript_structures": len(signatures),
                "cds_union_bp": cds_union_bp,
                "deletion_id": deletion["deletion_id"],
                "deletion_start": deletion["target_start"],
                "deletion_end": deletion["target_end"],
                "deletion_bp": deletion["deletion_bp"],
                "gene_coverage": f"{gene_coverage:.8f}",
                "cds_coverage": f"{cds_coverage:.8f}",
                "query_seqid": deletion["query_seqid"],
                "query_breakpoint": deletion["query_breakpoint"],
                "alignment_strand": deletion["strand"],
                "mapq": deletion["mapq"],
                "alignment_identity": deletion["alignment_identity"],
                "left_anchor_bp": deletion["left_anchor_bp"],
                "right_anchor_bp": deletion["right_anchor_bp"],
                "query_flank_start": deletion.get("query_flank_start", ""),
                "query_flank_end": deletion.get("query_flank_end", ""),
                "query_flank_expected_bp": deletion.get(
                    "query_flank_expected_bp", ""
                ),
                "query_flank_observed_bp": deletion.get(
                    "query_flank_observed_bp", ""
                ),
                "query_flank_ambiguous_bp": deletion.get(
                    "query_flank_ambiguous_bp", ""
                ),
                "query_flank_ambiguous_fraction": deletion.get(
                    "query_flank_ambiguous_fraction", ""
                ),
                "query_flank_status": deletion.get(
                    "query_flank_status", "not_evaluated"
                ),
            }
        )
        counts["genes_accepted"] += 1

    candidate_rows.sort(
        key=lambda row: (row["seqid"], int(row["gene_start"]), row["gene_id"])
    )
    with output_tsv.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(candidate_rows)
    manifest = {
        "schema_version": PAV_GENE_CATALOG_SCHEMA_VERSION,
        "generator": {"name": "PloidyPatch", "version": __version__},
        "inputs": {
            "gff": {"file_name": Path(gff_path).name, "sha256": _file_sha256(gff_path)},
            "deletions": {
                "file_name": Path(deletion_tsv_path).name,
                "sha256": _file_sha256(deletion_tsv_path),
            },
        },
        "parameters": {
            "min_gene_coverage": min_gene_coverage,
            "min_cds_coverage": min_cds_coverage,
            "single_deletion_required": True,
        },
        "counts": dict(sorted(counts.items())),
        "output": {
            "file_name": output_tsv.name,
            "rows": len(candidate_rows),
            "sha256": _file_sha256(output_tsv),
        },
    }
    with manifest_path.open("x", encoding="utf-8", newline="") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return manifest


def _parse_labeled_paths(values: list[str]) -> list[tuple[str, Path]]:
    parsed: list[tuple[str, Path]] = []
    seen_labels: set[str] = set()
    for value in values:
        label, separator, path_text = value.partition("=")
        if (
            not separator
            or not path_text
            or re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*", label) is None
        ):
            raise ValueError(
                "Catalog input must use a safe LABEL=PATH value: " + value
            )
        if label in seen_labels:
            raise ValueError(f"Duplicate catalog label: {label}")
        seen_labels.add(label)
        parsed.append((label, Path(path_text)))
    if not parsed:
        raise ValueError("At least one labeled catalog is required")
    return parsed


def reconcile_pav_gene_catalogs(
    *,
    catalog_inputs: list[str],
    output_tsv_path: str | Path,
    min_support: int = 2,
) -> dict[str, Any]:
    """Retain gene-level PAV candidates supported by multiple evidence runs."""

    inputs = _parse_labeled_paths(catalog_inputs)
    if min_support < 1 or min_support > len(inputs):
        raise ValueError("Minimum support must be within the number of catalogs")
    output_tsv = Path(output_tsv_path)
    manifest_path = Path(str(output_tsv) + ".manifest.json")
    collisions = [path for path in (output_tsv, manifest_path) if path.exists()]
    if collisions:
        raise FileExistsError(
            "Refusing to overwrite reconciled PAV artifact(s): "
            + ", ".join(str(path) for path in collisions)
        )

    gene_fields = [
        "gene_id",
        "seqid",
        "gene_start",
        "gene_end",
        "strand",
        "gene_span_bp",
        "transcript_structures",
        "cds_union_bp",
    ]
    evidence_fields = [
        "candidate_id",
        "deletion_id",
        "deletion_start",
        "deletion_end",
        "deletion_bp",
        "gene_coverage",
        "cds_coverage",
        "query_seqid",
        "query_breakpoint",
        "alignment_strand",
        "mapq",
        "alignment_identity",
        "left_anchor_bp",
        "right_anchor_bp",
        "query_flank_ambiguous_fraction",
        "query_flank_status",
    ]
    required = set(gene_fields + evidence_fields[:-2])
    canonical_by_gene: dict[str, dict[str, str]] = {}
    evidence_by_gene: dict[str, dict[str, dict[str, str]]] = {}
    input_manifests: list[dict[str, Any]] = []
    counts: Counter[str] = Counter()
    for label, path in inputs:
        rows_read = 0
        seen_genes: set[str] = set()
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            missing = required - set(reader.fieldnames or [])
            if reader.fieldnames is None or missing:
                raise ValueError(
                    f"Catalog {label} is missing required column(s): "
                    + ", ".join(sorted(missing))
                )
            for line_number, row in enumerate(reader, start=2):
                gene_id = row["gene_id"]
                if not gene_id:
                    raise ValueError(f"Empty gene_id in {label} at line {line_number}")
                if gene_id in seen_genes:
                    raise ValueError(
                        f"Duplicate gene_id {gene_id!r} in {label} at line {line_number}"
                    )
                seen_genes.add(gene_id)
                canonical = {field: row[field] for field in gene_fields}
                previous = canonical_by_gene.setdefault(gene_id, canonical)
                if previous != canonical:
                    raise ValueError(
                        f"Gene metadata disagree across catalogs for {gene_id}"
                    )
                evidence_by_gene.setdefault(gene_id, {})[label] = {
                    field: row.get(field, "") for field in evidence_fields
                }
                rows_read += 1
        counts[f"input_rows_{label}"] = rows_read
        input_manifests.append(
            {
                "label": label,
                "file_name": path.name,
                "sha256": _file_sha256(path),
                "rows": rows_read,
            }
        )

    labels = [label for label, _ in inputs]
    fieldnames = [
        "candidate_id",
        *gene_fields,
        "support_count",
        "support_sources",
    ]
    for label in labels:
        fieldnames.extend(
            f"{label}_{field}" for field in evidence_fields
        )
    output_rows: list[dict[str, str | int]] = []
    for gene_id, source_evidence in evidence_by_gene.items():
        support_labels = [label for label in labels if label in source_evidence]
        support = len(support_labels)
        counts[f"genes_with_{support}_source_support"] += 1
        if support < min_support:
            counts["genes_below_minimum_support"] += 1
            continue
        digest = hashlib.sha256(
            (gene_id + "\t" + "\t".join(support_labels)).encode("utf-8")
        ).hexdigest()[:20]
        row: dict[str, str | int] = {
            "candidate_id": f"PPRPAV-{digest}",
            **canonical_by_gene[gene_id],
            "support_count": support,
            "support_sources": ",".join(support_labels),
        }
        for label in labels:
            evidence = source_evidence.get(label, {})
            for field in evidence_fields:
                row[f"{label}_{field}"] = evidence.get(field, "")
        output_rows.append(row)
        counts["genes_accepted"] += 1
    output_rows.sort(
        key=lambda row: (row["seqid"], int(row["gene_start"]), row["gene_id"])
    )
    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    with output_tsv.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(output_rows)
    manifest = {
        "schema_version": PAV_RECONCILED_CATALOG_SCHEMA_VERSION,
        "generator": {"name": "PloidyPatch", "version": __version__},
        "inputs": input_manifests,
        "parameters": {"min_support": min_support, "labels": labels},
        "counts": dict(sorted(counts.items())),
        "output": {
            "file_name": output_tsv.name,
            "rows": len(output_rows),
            "sha256": _file_sha256(output_tsv),
        },
    }
    with manifest_path.open("x", encoding="utf-8", newline="") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return manifest
