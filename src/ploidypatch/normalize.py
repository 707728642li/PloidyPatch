from __future__ import annotations

import csv
import hashlib
import json
import os
import shutil
import tempfile
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

from . import __version__
from .audit import DNA_ALPHABET, PROTEIN_ALPHABET, STOP_CODONS
from .gff import TRANSCRIPT_TYPES, parse_attributes
from .io import fasta_relation_id, iter_fasta, normalize_feature_id, open_text
from .perturb import GffDocument, GffRecord, read_gff_document
from .synteny_io import discover_primary_chromosome_seqids


PRIMARY_BUNDLE_SCHEMA_VERSION = "ploidypatch.primary_bundle_manifest.v1"
PRIMARY_ANNOTATION_BUNDLE_SCHEMA_VERSION = (
    "ploidypatch.primary_annotation_bundle_manifest.v1"
)
PROVIDER_GFF3_COMPATIBILITY_SCHEMA_VERSION = (
    "ploidypatch.provider_gff3_compatibility_manifest.v1"
)
PRIMARY_SEQID_TABLE_COLUMNS = ("seqid", "chromosome_label")
OUTPUT_NAMES = {
    "gff3": "primary_chromosomes.gff3",
    "genome": "primary_chromosomes.genome.fa",
    "fai": "primary_chromosomes.genome.fa.fai",
    "protein": "primary_chromosomes.protein.fa",
    "cds": "primary_chromosomes.cds.fa",
    "exclusions": "excluded_translations.tsv",
    "manifest": "manifest.json",
}
PRIMARY_ANNOTATION_OUTPUT_NAMES = {
    "gff3": "primary_chromosomes.gff3",
    "genome": "primary_chromosomes.genome.fa",
    "fai": "primary_chromosomes.genome.fa.fai",
    "manifest": "manifest.json",
}
PROVIDER_GFF3_OUTPUT_NAMES = {
    "gff3": "sanitized.gff3",
    "manifest": "manifest.json",
}


def _file_sha256(path: str | Path, block_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        while block := handle.read(block_size):
            digest.update(block)
    return digest.hexdigest()


def read_primary_seqid_table(
    path: str | Path,
) -> tuple[frozenset[str], dict[str, str]]:
    """Read an explicit, auditable primary-chromosome selection table.

    Some plant community annotations do not emit NCBI ``region`` or Ensembl
    ``chromosome`` records.  A strict two-column table avoids silently
    guessing primary chromosomes from provider-specific sequence names.
    """

    table_path = Path(path)
    with table_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != list(PRIMARY_SEQID_TABLE_COLUMNS):
            raise ValueError(
                "Primary seqid table must have exactly these tab-separated "
                "columns: " + ", ".join(PRIMARY_SEQID_TABLE_COLUMNS)
            )
        labels: dict[str, str] = {}
        for row_number, row in enumerate(reader, start=2):
            seqid = (row.get("seqid") or "").strip()
            chromosome_label = (row.get("chromosome_label") or "").strip()
            if not seqid or not chromosome_label:
                raise ValueError(
                    "Primary seqid table contains an empty value on line "
                    f"{row_number}"
                )
            if seqid in labels:
                raise ValueError(
                    f"Primary seqid table contains duplicate seqid {seqid!r}"
                )
            labels[seqid] = chromosome_label
    if not labels:
        raise ValueError("Primary seqid table contains no chromosome rows")
    return frozenset(labels), labels


def _repair_unescaped_note_semicolons(text: str) -> tuple[str, int] | None:
    """Repair only bare semicolon continuations of a GFF3 ``Note`` value.

    Community annotations occasionally contain ``Note=first;second`` where the
    semicolon is intended as text rather than an attribute separator.  This
    helper deliberately refuses every broader malformed-attribute pattern.
    """

    items = text.split(";")
    repaired: list[str] = []
    previous_key: str | None = None
    repair_count = 0
    for item in items:
        stripped = item.strip()
        if not stripped:
            continue
        if "=" in stripped:
            key = stripped.split("=", 1)[0].strip()
            if not key:
                return None
            repaired.append(stripped)
            previous_key = key
            continue
        if previous_key != "Note" or not repaired:
            return None
        repaired[-1] += "%3B" + stripped
        repair_count += 1
    if repair_count == 0:
        return None
    repaired_text = ";".join(repaired)
    _, malformed = parse_attributes(repaired_text)
    if malformed:
        return None
    return repaired_text, repair_count


def normalize_provider_gff3(
    *,
    gff_path: str | Path,
    output_dir: str | Path,
    repair_unescaped_note_semicolons: bool = False,
    drop_invalid_intron_intervals: bool = False,
    strip_embedded_fasta: bool = False,
) -> dict[str, Any]:
    """Create a strictly parseable GFF3 with narrow, audited compatibility fixes.

    All compatibility actions are disabled by default.  Enabled actions are
    intentionally limited to provider defects that do not alter gene, exon or
    CDS coordinates: unescaped semicolons inside ``Note`` values, redundant
    negative-length ``intron`` records, and an embedded FASTA section when a
    separate genome FASTA is used downstream.
    """

    source = Path(gff_path)
    output = Path(output_dir)
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite output directory: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(
        tempfile.mkdtemp(prefix=f".{output.name}.partial.", dir=output.parent)
    )
    try:
        paths = {
            role: staging / name for role, name in PROVIDER_GFF3_OUTPUT_NAMES.items()
        }
        counts: Counter[str] = Counter()
        repaired_examples: list[dict[str, Any]] = []
        dropped_examples: list[dict[str, Any]] = []
        embedded_fasta_line: int | None = None
        in_embedded_fasta = False
        with open_text(source) as source_handle, paths["gff3"].open(
            "x", encoding="utf-8", newline=""
        ) as output_handle:
            for line_number, raw_line in enumerate(source_handle, start=1):
                counts["input_lines"] += 1
                stripped = raw_line.rstrip("\r\n")
                if in_embedded_fasta:
                    counts["stripped_embedded_fasta_lines"] += 1
                    continue
                if stripped.startswith("##FASTA"):
                    if not strip_embedded_fasta:
                        raise ValueError(
                            "Embedded FASTA section encountered on line "
                            f"{line_number}; enable strip_embedded_fasta explicitly"
                        )
                    embedded_fasta_line = line_number
                    in_embedded_fasta = True
                    counts["stripped_embedded_fasta_directives"] += 1
                    continue
                if not stripped or stripped.startswith("#"):
                    output_handle.write(raw_line)
                    counts["output_nonfeature_lines"] += 1
                    continue
                fields = stripped.split("\t")
                if len(fields) != 9:
                    raise ValueError(
                        f"Malformed GFF3 feature line {line_number}: expected 9 columns"
                    )
                try:
                    start = int(fields[3])
                    end = int(fields[4])
                except ValueError as exc:
                    raise ValueError(
                        f"Malformed GFF3 coordinates on line {line_number}"
                    ) from exc
                if start < 1 or end < start:
                    can_drop = (
                        drop_invalid_intron_intervals
                        and fields[2] == "intron"
                        and start >= 1
                        and end < start
                    )
                    if not can_drop:
                        raise ValueError(
                            f"Invalid GFF3 interval on line {line_number}"
                        )
                    counts["dropped_invalid_intron_intervals"] += 1
                    if len(dropped_examples) < 20:
                        dropped_examples.append(
                            {
                                "line_number": line_number,
                                "seqid": fields[0],
                                "start": start,
                                "end": end,
                            }
                        )
                    continue
                _, malformed = parse_attributes(fields[8])
                if malformed:
                    repaired = (
                        _repair_unescaped_note_semicolons(fields[8])
                        if repair_unescaped_note_semicolons
                        else None
                    )
                    if repaired is None:
                        raise ValueError(
                            "Malformed GFF3 attributes on line "
                            f"{line_number}: {malformed} invalid field(s)"
                        )
                    fields[8], repair_count = repaired
                    counts["repaired_note_records"] += 1
                    counts["repaired_note_semicolons"] += repair_count
                    if len(repaired_examples) < 20:
                        repaired_examples.append(
                            {
                                "line_number": line_number,
                                "seqid": fields[0],
                                "feature_type": fields[2],
                                "repaired_semicolons": repair_count,
                            }
                        )
                    raw_line = "\t".join(fields) + "\n"
                output_handle.write(raw_line)
                counts["output_feature_lines"] += 1

        validated = read_gff_document(paths["gff3"])
        if len(validated.records) != counts["output_feature_lines"]:
            raise RuntimeError("Strict validation record count does not match output")
        manifest: dict[str, Any] = {
            "schema_version": PROVIDER_GFF3_COMPATIBILITY_SCHEMA_VERSION,
            "generator": {"name": "PloidyPatch", "version": __version__},
            "source": {
                "file_name": source.name,
                "bytes": source.stat().st_size,
                "sha256": _file_sha256(source),
            },
            "policy": {
                "repair_unescaped_note_semicolons": (
                    repair_unescaped_note_semicolons
                ),
                "drop_invalid_intron_intervals": drop_invalid_intron_intervals,
                "strip_embedded_fasta": strip_embedded_fasta,
            },
            "observed": {
                **dict(sorted(counts.items())),
                "embedded_fasta_line": embedded_fasta_line,
                "repaired_examples": repaired_examples,
                "dropped_examples": dropped_examples,
                "strict_validation_records": len(validated.records),
            },
            "output": {
                "file_name": paths["gff3"].name,
                "bytes": paths["gff3"].stat().st_size,
                "sha256": _file_sha256(paths["gff3"]),
            },
        }
        with paths["manifest"].open("x", encoding="utf-8", newline="") as handle:
            json.dump(manifest, handle, indent=2, sort_keys=True)
            handle.write("\n")
        os.replace(staging, output)
        return manifest
    except BaseException:
        shutil.rmtree(staging, ignore_errors=True)
        raise


def _read_relevant_fasta(
    path: str | Path,
    expected_relation_ids: set[str],
) -> dict[str, tuple[str, str, str]]:
    records: dict[str, tuple[str, str, str]] = {}
    for record_id, header, sequence in iter_fasta(path):
        relation_id, _ = fasta_relation_id(record_id, header)
        if relation_id not in expected_relation_ids:
            continue
        if relation_id in records:
            raise ValueError(f"Duplicate FASTA relation ID {relation_id!r} in {path}")
        records[relation_id] = (record_id, header, sequence)
    return records


def _translation_exclusion_reasons(
    transcript_id: str,
    protein_records: dict[str, tuple[str, str, str]],
    cds_records: dict[str, tuple[str, str, str]],
) -> tuple[str, ...]:
    reasons: list[str] = []
    protein_record = protein_records.get(transcript_id)
    cds_record = cds_records.get(transcript_id)
    if protein_record is None:
        reasons.append("missing_protein")
    if cds_record is None:
        reasons.append("missing_cds")
    if protein_record is None or cds_record is None:
        return tuple(reasons)

    protein = protein_record[2].upper()
    cds = cds_record[2].upper()
    if not protein:
        reasons.append("empty_protein")
    if not cds:
        reasons.append("empty_cds")
    if set(protein) - PROTEIN_ALPHABET:
        reasons.append("invalid_protein_character")
    if set(cds) - DNA_ALPHABET:
        reasons.append("invalid_cds_character")
    if "*" in protein[:-1]:
        reasons.append("protein_internal_stop")
    codons = [cds[index : index + 3] for index in range(0, len(cds) - 2, 3)]
    if any(codon in STOP_CODONS for codon in codons[:-1]):
        reasons.append("cds_internal_stop")
    protein_length = len(protein) - int(protein.endswith("*"))
    cds_length = len(cds) // 3 - int(bool(codons) and codons[-1] in STOP_CODONS)
    if protein_length != cds_length:
        reasons.append("translation_length_mismatch")
    return tuple(reasons)


def _translation_context(
    document: GffDocument,
) -> tuple[set[str], dict[str, str], dict[str, set[str]], dict[str, str]]:
    transcript_raw_by_normalized: dict[str, str] = {}
    transcript_gene: dict[str, str] = {}
    gene_transcripts: dict[str, set[str]] = defaultdict(set)
    expected: set[str] = set()
    for record in document.records:
        if record.feature_type in TRANSCRIPT_TYPES and record.feature_id:
            normalized = normalize_feature_id(record.feature_id)
            existing = transcript_raw_by_normalized.get(normalized)
            if existing is not None and existing != record.feature_id:
                raise ValueError(
                    f"Transcript IDs collide after normalization: {existing}, "
                    f"{record.feature_id}"
                )
            transcript_raw_by_normalized[normalized] = record.feature_id
            if len(record.parents) == 1:
                transcript_gene[normalized] = record.parents[0]
                gene_transcripts[record.parents[0]].add(normalized)
            biotype = record.attributes.get("biotype") or record.attributes.get(
                "transcript_biotype"
            )
            if biotype and biotype.lower() in {"protein_coding", "protein-coding"}:
                expected.add(normalized)
        if record.feature_type == "CDS":
            expected.update(normalize_feature_id(parent) for parent in record.parents)
    unknown = expected - set(transcript_raw_by_normalized)
    if unknown:
        raise ValueError(
            "CDS parent(s) are not recognized transcripts: "
            + ", ".join(sorted(unknown)[:10])
        )
    return expected, transcript_gene, gene_transcripts, transcript_raw_by_normalized


def _retained_records(
    document: GffDocument,
    excluded_transcripts: set[str],
    gene_transcripts: dict[str, set[str]],
    transcript_raw_by_normalized: dict[str, str],
) -> tuple[list[GffRecord], set[str]]:
    excluded_raw = {
        transcript_raw_by_normalized[transcript_id]
        for transcript_id in excluded_transcripts
    }
    excluded_genes = {
        gene_id
        for gene_id, transcripts in gene_transcripts.items()
        if transcripts and transcripts <= excluded_transcripts
    }
    blocked_ids = set(excluded_raw) | excluded_genes
    retained: list[GffRecord] = []
    for record in document.records:
        feature_id = record.feature_id
        if feature_id and feature_id in blocked_ids:
            continue
        blocked_parents = set(record.parents) & blocked_ids
        if blocked_parents:
            if blocked_parents != set(record.parents):
                raise ValueError(
                    "Cannot subset a feature with mixed retained/excluded parents: "
                    f"line {record.line_number}"
                )
            if feature_id:
                blocked_ids.add(feature_id)
            continue
        retained.append(record)
    return retained, excluded_genes


def _write_filtered_gff(
    document: GffDocument,
    retained_records: list[GffRecord],
    primary_seqids: frozenset[str],
    output_path: Path,
) -> None:
    retained_lines = {record.line_number for record in retained_records}
    with output_path.open("x", encoding="utf-8", newline="") as handle:
        for line_number, raw_line in enumerate(document.lines, start=1):
            if raw_line.startswith("##FASTA"):
                raise ValueError("Embedded FASTA sections are not supported")
            if raw_line.startswith("##sequence-region"):
                fields = raw_line.rstrip("\r\n").split()
                if len(fields) >= 2 and fields[1] in primary_seqids:
                    handle.write(raw_line)
                continue
            if raw_line.startswith("#") or not raw_line.strip():
                handle.write(raw_line)
            elif line_number in retained_lines:
                handle.write(raw_line)


def _write_fasta_records(
    records: dict[str, tuple[str, str, str]],
    retained_ids: set[str],
    output_path: Path,
    line_bases: int = 70,
) -> int:
    missing = retained_ids - set(records)
    if missing:
        raise ValueError(
            f"Cannot write {output_path.name}; missing relation IDs: "
            + ", ".join(sorted(missing)[:10])
        )
    written = 0
    with output_path.open("x", encoding="utf-8", newline="") as handle:
        for relation_id, (_, header, sequence) in records.items():
            if relation_id not in retained_ids:
                continue
            handle.write(f">{header}\n")
            for index in range(0, len(sequence), line_bases):
                handle.write(sequence[index : index + line_bases] + "\n")
            written += 1
    return written


def _write_primary_genome(
    source_path: str | Path,
    primary_seqids: frozenset[str],
    output_path: Path,
    fai_path: Path,
    line_bases: int = 60,
    canonical_headers: bool = False,
) -> dict[str, int]:
    lengths: dict[str, int] = {}
    fai_rows: list[str] = []
    with output_path.open("xb") as handle:
        for record_id, header, sequence in iter_fasta(source_path):
            if record_id not in primary_seqids:
                continue
            if record_id in lengths:
                raise ValueError(f"Duplicate genome sequence ID: {record_id}")
            output_header = record_id if canonical_headers else header
            header_bytes = f">{output_header}\n".encode("utf-8")
            handle.write(header_bytes)
            offset = handle.tell()
            sequence_bytes = sequence.encode("ascii")
            if not sequence_bytes:
                raise ValueError(f"Empty primary genome sequence: {record_id}")
            for index in range(0, len(sequence_bytes), line_bases):
                handle.write(sequence_bytes[index : index + line_bases] + b"\n")
            lengths[record_id] = len(sequence_bytes)
            first_line_bases = min(line_bases, len(sequence_bytes))
            fai_rows.append(
                f"{record_id}\t{len(sequence_bytes)}\t{offset}\t"
                f"{first_line_bases}\t{first_line_bases + 1}\n"
            )
    missing = set(primary_seqids) - set(lengths)
    if missing:
        raise ValueError(
            "Primary chromosomes are missing from genome FASTA: "
            + ", ".join(sorted(missing))
        )
    with fai_path.open("x", encoding="utf-8", newline="") as handle:
        handle.writelines(fai_rows)
    return lengths


def prepare_primary_annotation_bundle(
    *,
    gff_path: str | Path,
    genome_path: str | Path,
    output_dir: str | Path,
    primary_seqid_table_path: str | Path | None = None,
    canonical_fasta_headers: bool = False,
) -> dict[str, Any]:
    """Subset a genome and GFF to declared primary chromosomes.

    Unlike :func:`prepare_ncbi_primary_bundle`, this operation is deliberately
    structural: it preserves every GFF record on retained chromosomes and does
    not inspect proteins, CDS FASTA records, or translation quality.  It is
    therefore suitable for fair annotation-transfer baseline inputs.
    """

    output = Path(output_dir)
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite output directory: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(
        tempfile.mkdtemp(prefix=f".{output.name}.partial.", dir=output.parent)
    )
    try:
        paths = {
            role: staging / name
            for role, name in PRIMARY_ANNOTATION_OUTPUT_NAMES.items()
        }
        if primary_seqid_table_path is None:
            primary_seqids, chromosome_labels = discover_primary_chromosome_seqids(
                gff_path
            )
            selection_method = "GFF region/chromosome feature discovery"
        else:
            primary_seqids, chromosome_labels = read_primary_seqid_table(
                primary_seqid_table_path
            )
            selection_method = "explicit primary seqid table"
        document = read_gff_document(gff_path, include_seqids=primary_seqids)
        retained_records = document.records
        annotated_seqids = {record.seqid for record in retained_records}
        missing_annotation_seqids = set(primary_seqids) - annotated_seqids
        if missing_annotation_seqids:
            raise ValueError(
                "Explicit primary seqid(s) have no GFF records: "
                + ", ".join(sorted(missing_annotation_seqids))
            )
        _write_filtered_gff(
            document, retained_records, primary_seqids, paths["gff3"]
        )
        genome_lengths = _write_primary_genome(
            genome_path,
            primary_seqids,
            paths["genome"],
            paths["fai"],
            canonical_headers=canonical_fasta_headers,
        )

        feature_counts = Counter(record.feature_type for record in retained_records)
        source_paths = {"gff3": Path(gff_path), "genome": Path(genome_path)}
        if primary_seqid_table_path is not None:
            source_paths["primary_seqid_table"] = Path(primary_seqid_table_path)
        manifest: dict[str, Any] = {
            "schema_version": PRIMARY_ANNOTATION_BUNDLE_SCHEMA_VERSION,
            "generator": {"name": "PloidyPatch", "version": __version__},
            "sources": {
                role: {
                    "file_name": path.name,
                    "bytes": path.stat().st_size,
                    "sha256": _file_sha256(path),
                }
                for role, path in source_paths.items()
            },
            "selection": {
                "method": selection_method,
                "translation_filtering": False,
                "primary_chromosomes": len(primary_seqids),
                "chromosome_labels": dict(sorted(chromosome_labels.items())),
                "retained_gff_records": len(retained_records),
                "retained_feature_counts": dict(sorted(feature_counts.items())),
                "genome_bp": sum(genome_lengths.values()),
            },
            "outputs": {},
        }
        if canonical_fasta_headers:
            manifest["selection"]["genome_fasta_headers"] = "exact_seqid"
        for role in ("gff3", "genome", "fai"):
            path = paths[role]
            manifest["outputs"][role] = {
                "file_name": path.name,
                "bytes": path.stat().st_size,
                "sha256": _file_sha256(path),
            }
        with paths["manifest"].open("x", encoding="utf-8", newline="") as handle:
            json.dump(manifest, handle, indent=2, sort_keys=True)
            handle.write("\n")
        os.replace(staging, output)
        return manifest
    except BaseException:
        shutil.rmtree(staging, ignore_errors=True)
        raise


def prepare_ncbi_primary_bundle(
    *,
    gff_path: str | Path,
    protein_path: str | Path,
    cds_path: str | Path,
    genome_path: str | Path,
    output_dir: str | Path,
) -> dict[str, Any]:
    """Create a clean, relation-consistent NCBI primary-chromosome bundle."""

    output = Path(output_dir)
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite output directory: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(
        tempfile.mkdtemp(prefix=f".{output.name}.partial.", dir=output.parent)
    )
    try:
        paths = {role: staging / name for role, name in OUTPUT_NAMES.items()}
        primary_seqids, chromosome_labels = discover_primary_chromosome_seqids(gff_path)
        document = read_gff_document(gff_path, include_seqids=primary_seqids)
        (
            expected_transcripts,
            transcript_gene,
            gene_transcripts,
            transcript_raw_by_normalized,
        ) = _translation_context(document)
        protein_records = _read_relevant_fasta(protein_path, expected_transcripts)
        cds_records = _read_relevant_fasta(cds_path, expected_transcripts)
        exclusion_reasons = {
            transcript_id: reasons
            for transcript_id in sorted(expected_transcripts)
            if (
                reasons := _translation_exclusion_reasons(
                    transcript_id, protein_records, cds_records
                )
            )
        }
        excluded_transcripts = set(exclusion_reasons)
        retained_translations = expected_transcripts - excluded_transcripts
        retained_records, excluded_genes = _retained_records(
            document,
            excluded_transcripts,
            gene_transcripts,
            transcript_raw_by_normalized,
        )

        _write_filtered_gff(document, retained_records, primary_seqids, paths["gff3"])
        protein_count = _write_fasta_records(
            protein_records, retained_translations, paths["protein"]
        )
        cds_count = _write_fasta_records(
            cds_records, retained_translations, paths["cds"]
        )
        genome_lengths = _write_primary_genome(
            genome_path, primary_seqids, paths["genome"], paths["fai"]
        )

        reason_counts: Counter[str] = Counter()
        with paths["exclusions"].open("x", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(("transcript_id", "gene_id", "reasons"))
            for transcript_id, reasons in exclusion_reasons.items():
                reason_counts.update(reasons)
                writer.writerow(
                    (
                        transcript_id,
                        transcript_gene.get(transcript_id, ""),
                        ";".join(reasons),
                    )
                )

        source_paths = {
            "gff3": Path(gff_path),
            "protein": Path(protein_path),
            "cds": Path(cds_path),
            "genome": Path(genome_path),
        }
        manifest: dict[str, Any] = {
            "schema_version": PRIMARY_BUNDLE_SCHEMA_VERSION,
            "generator": {"name": "PloidyPatch", "version": __version__},
            "sources": {
                role: {
                    "file_name": path.name,
                    "bytes": path.stat().st_size,
                    "sha256": _file_sha256(path),
                }
                for role, path in source_paths.items()
            },
            "selection": {
                "method": "NCBI region feature with genome=chromosome",
                "primary_chromosomes": len(primary_seqids),
                "chromosome_labels": chromosome_labels,
                "source_translation_transcripts": len(expected_transcripts),
                "retained_translation_transcripts": len(retained_translations),
                "excluded_translation_transcripts": len(excluded_transcripts),
                "excluded_genes_without_retained_transcript": len(excluded_genes),
                "exclusion_reason_counts": dict(sorted(reason_counts.items())),
                "retained_gff_records": len(retained_records),
                "genome_bp": sum(genome_lengths.values()),
            },
            "outputs": {},
        }
        for role in ("gff3", "genome", "fai", "protein", "cds", "exclusions"):
            path = paths[role]
            manifest["outputs"][role] = {
                "file_name": path.name,
                "bytes": path.stat().st_size,
                "sha256": _file_sha256(path),
            }
        manifest["outputs"]["protein"]["records"] = protein_count
        manifest["outputs"]["cds"]["records"] = cds_count
        with paths["manifest"].open("x", encoding="utf-8", newline="") as handle:
            json.dump(manifest, handle, indent=2, sort_keys=True)
            handle.write("\n")
        os.replace(staging, output)
        return manifest
    except BaseException:
        shutil.rmtree(staging, ignore_errors=True)
        raise
