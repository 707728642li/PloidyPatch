from __future__ import annotations

import csv
import hashlib
import json
import random
import re
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

from . import __version__
from .baseline import _file_sha256
from .consensus import _merge_intervals
from .gff import parse_attributes


CSV_FIELD_SIZE_LIMIT = 2**31 - 1
csv.field_size_limit(CSV_FIELD_SIZE_LIMIT)


ISOSEQ_PREPARE_SCHEMA = "ploidypatch.isoseq_b73_prepare.v1"
ISOSEQ_VALIDATION_SCHEMA = "ploidypatch.isoseq_candidate_validation.v2"
ISOSEQ_QUERY_FILTER_SCHEMA = "ploidypatch.candidate_query_paf_filter.v1"
ALIGNMENT_STRAND_SOURCES = ("query_orientation", "minimap2_ts")
CIGAR_TOKEN = re.compile(r"(\d+)([MIDNSHP=X])")
CANONICAL_PLUS = frozenset({"GTAG", "GCAG", "ATAC"})
CANONICAL_MINUS = frozenset({"CTAC", "CTGC", "GTAT"})


def _refuse_outputs(paths: Iterable[Path]) -> None:
    collisions = [str(path) for path in paths if path.exists()]
    if collisions:
        raise FileExistsError("Refusing to overwrite output(s): " + ", ".join(collisions))


def prepare_b73_isoseq_transcripts(
    *,
    fasta_path: str | Path,
    count_csv_path: str | Path,
    output_fasta_path: str | Path,
    output_count_tsv_path: str | Path,
    minimum_b73_full_length_reads: int = 2,
) -> dict[str, Any]:
    """Select transcript sequences by a predeclared pure-B73 read threshold."""

    if minimum_b73_full_length_reads < 1:
        raise ValueError("minimum_b73_full_length_reads must be positive")
    fasta = Path(fasta_path)
    counts = Path(count_csv_path)
    output_fasta = Path(output_fasta_path)
    output_counts = Path(output_count_tsv_path)
    manifest_path = Path(str(output_fasta) + ".manifest.json")
    for required in (fasta, counts):
        if not required.is_file() or required.stat().st_size == 0:
            raise ValueError(f"Missing or empty Iso-Seq input: {required}")
    _refuse_outputs((output_fasta, output_counts, manifest_path))

    b73_columns = ("EM1", "R1", "END1")
    b73_counts: dict[str, tuple[int, int, int, int]] = {}
    with counts.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames or "id" not in reader.fieldnames:
            raise ValueError("Iso-Seq count table lacks id")
        missing = set(b73_columns) - set(reader.fieldnames)
        if missing:
            raise ValueError("Iso-Seq count table lacks B73 columns: " + ",".join(sorted(missing)))
        for line_number, row in enumerate(reader, start=2):
            transcript_id = (row.get("id") or "").strip()
            if not transcript_id or transcript_id in b73_counts:
                raise ValueError(f"Empty or duplicate Iso-Seq id at line {line_number}")
            try:
                values = tuple(int(row[column]) for column in b73_columns)
            except (TypeError, ValueError) as exc:
                raise ValueError(f"Invalid B73 count at line {line_number}") from exc
            if any(value < 0 for value in values):
                raise ValueError(f"Negative B73 count at line {line_number}")
            b73_counts[transcript_id] = (*values, sum(values))

    selected = {
        transcript_id
        for transcript_id, values in b73_counts.items()
        if values[-1] >= minimum_b73_full_length_reads
    }
    observed_ids: set[str] = set()
    selected_bases = 0
    selected_records = 0
    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    keep = False
    with fasta.open("r", encoding="ascii", newline="") as source, output_fasta.open(
        "x", encoding="ascii", newline=""
    ) as destination:
        for line_number, raw_line in enumerate(source, start=1):
            if raw_line.startswith(">"):
                transcript_id = raw_line[1:].strip().split()[0]
                if not transcript_id or transcript_id in observed_ids:
                    raise ValueError(f"Empty or duplicate FASTA id at line {line_number}")
                observed_ids.add(transcript_id)
                keep = transcript_id in selected
                if keep:
                    selected_records += 1
                    destination.write(f">{transcript_id}\n")
            elif keep:
                sequence = raw_line.strip()
                if not sequence:
                    continue
                if set(sequence.upper()) - set("ACGTUN"):
                    raise ValueError(f"Invalid FASTA sequence at line {line_number}")
                selected_bases += len(sequence)
                destination.write(sequence.upper().replace("U", "T") + "\n")
    missing_selected = selected - observed_ids
    if missing_selected:
        raise ValueError(f"Selected count IDs absent from FASTA: {len(missing_selected)}")
    unknown_fasta = observed_ids - set(b73_counts)
    if unknown_fasta:
        raise ValueError(f"FASTA IDs absent from count table: {len(unknown_fasta)}")

    output_counts.parent.mkdir(parents=True, exist_ok=True)
    with output_counts.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("transcript_id", "EM1", "R1", "END1", "b73_full_length_reads"))
        for transcript_id in sorted(selected):
            writer.writerow((transcript_id, *b73_counts[transcript_id]))

    manifest: dict[str, Any] = {
        "schema_version": ISOSEQ_PREPARE_SCHEMA,
        "generator": {"name": "PloidyPatch", "version": __version__},
        "inputs": {
            "fasta": {"file_name": fasta.name, "sha256": _file_sha256(fasta)},
            "counts": {"file_name": counts.name, "sha256": _file_sha256(counts)},
        },
        "parameters": {
            "b73_columns": list(b73_columns),
            "minimum_b73_full_length_reads": minimum_b73_full_length_reads,
            "other_genotypes_excluded": True,
        },
        "counts": {
            "count_table_transcripts": len(b73_counts),
            "fasta_transcripts": len(observed_ids),
            "selected_transcripts": selected_records,
            "selected_bases": selected_bases,
        },
        "outputs": {
            "fasta": {"file_name": output_fasta.name, "sha256": _file_sha256(output_fasta)},
            "counts": {"file_name": output_counts.name, "sha256": _file_sha256(output_counts)},
        },
    }
    with manifest_path.open("x", encoding="utf-8", newline="") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return manifest


@dataclass(frozen=True)
class Candidate:
    digest: str
    seqid: str
    start: int
    end: int
    strand: str
    cds: tuple[tuple[int, int], ...]

    @property
    def cds_bp(self) -> int:
        return sum(end - start + 1 for start, end in self.cds)

    @property
    def junctions(self) -> tuple[tuple[int, int], ...]:
        return tuple((left[1], right[0]) for left, right in zip(self.cds, self.cds[1:]))


@dataclass(frozen=True)
class Alignment:
    query: str
    query_length: int
    seqid: str
    start: int
    end: int
    strand: str | None
    mapq: int
    query_coverage: float
    identity: float
    alignment_score: int
    primary: bool
    exonic_blocks: tuple[tuple[int, int], ...]
    junctions: tuple[tuple[int, int], ...]


def _parse_paf_tags(fields: list[str]) -> dict[str, str]:
    return {
        tag[:2]: tag[5:]
        for tag in fields[12:]
        if len(tag) >= 6 and tag[2:5] in {":i:", ":A:", ":Z:", ":f:"}
    }


def _reference_transcript_strand(
    query_target_strand: str,
    transcript_relation: str | None,
    alignment_strand_source: str,
) -> str | None:
    """Return transcript strand relative to the reference target.

    In minimap2 splice PAF, column 5 says whether query and target have the
    same orientation.  The ``ts`` tag says whether the read/query strand is
    the same as the inferred transcript strand.  Their sign product is the
    biological transcript strand relative to the target.  Unspliced reads do
    not carry ``ts`` and therefore cannot provide strand-specific evidence for
    an unstranded cDNA library.
    """

    if alignment_strand_source not in ALIGNMENT_STRAND_SOURCES:
        raise ValueError(
            "alignment_strand_source must be one of: "
            + ",".join(ALIGNMENT_STRAND_SOURCES)
        )
    if query_target_strand not in {"+", "-"}:
        raise ValueError("PAF query/target strand must be + or -")
    if alignment_strand_source == "query_orientation":
        return query_target_strand
    if transcript_relation is None:
        return None
    if transcript_relation not in {"+", "-"}:
        raise ValueError("minimap2 ts tag must be + or -")
    return "+" if query_target_strand == transcript_relation else "-"


class FastaIndex:
    def __init__(self, fasta_path: Path, fai_path: Path):
        self.handle = fasta_path.open("rb")
        self.entries: dict[str, tuple[int, int, int, int]] = {}
        with fai_path.open("r", encoding="utf-8", newline="") as handle:
            for line_number, raw_line in enumerate(handle, start=1):
                fields = raw_line.rstrip("\n").split("\t")
                if len(fields) < 5:
                    raise ValueError(f"Malformed FAI line {line_number}")
                self.entries[fields[0]] = tuple(int(value) for value in fields[1:5])

    def close(self) -> None:
        self.handle.close()

    def fetch(self, seqid: str, start: int, end: int) -> str:
        if seqid not in self.entries:
            raise ValueError(f"Sequence missing from FAI: {seqid}")
        length, offset, line_bases, line_width = self.entries[seqid]
        start = max(1, start)
        end = min(length, end)
        if start > end:
            return ""
        zero_start = start - 1
        zero_end = end - 1
        file_start = offset + (zero_start // line_bases) * line_width + zero_start % line_bases
        file_end = offset + (zero_end // line_bases) * line_width + zero_end % line_bases
        self.handle.seek(file_start)
        raw = self.handle.read(file_end - file_start + 1)
        return raw.replace(b"\n", b"").replace(b"\r", b"").decode("ascii").upper()


def _read_candidates(path: Path) -> list[Candidate]:
    transcripts: dict[str, tuple[str, str, int, int, str]] = {}
    cds_by_transcript: dict[str, list[tuple[int, int]]] = defaultdict(list)
    with path.open("r", encoding="utf-8", newline="") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            if not raw_line or raw_line.startswith("#"):
                continue
            fields = raw_line.rstrip("\r\n").split("\t")
            if len(fields) != 9 or fields[1] != "PloidyPatchConsensus":
                continue
            attributes, malformed = parse_attributes(fields[8])
            if malformed:
                raise ValueError(f"Malformed candidate attributes at line {line_number}")
            if fields[2] in {"mRNA", "transcript"}:
                transcript_id = attributes.get("ID", "")
                digest = attributes.get("consensus_digest", "")
                if not transcript_id or len(digest) != 64 or transcript_id in transcripts:
                    raise ValueError(f"Invalid candidate transcript at line {line_number}")
                transcripts[transcript_id] = (
                    fields[0], fields[6], int(fields[3]), int(fields[4]), digest
                )
            elif fields[2] == "CDS":
                parents = [value for value in attributes.get("Parent", "").split(",") if value]
                if len(parents) != 1:
                    raise ValueError(f"Candidate CDS lacks one parent at line {line_number}")
                cds_by_transcript[parents[0]].append((int(fields[3]), int(fields[4])))
    orphan = set(cds_by_transcript) - set(transcripts)
    if orphan:
        raise ValueError("Candidate CDS references unknown transcript")
    candidates: list[Candidate] = []
    digests: set[str] = set()
    for transcript_id, (seqid, strand, start, end, digest) in transcripts.items():
        cds = _merge_intervals(cds_by_transcript.get(transcript_id, ()))
        if not cds or strand not in {"+", "-"} or digest in digests:
            raise ValueError(f"Invalid candidate hierarchy: {transcript_id}")
        digests.add(digest)
        candidates.append(Candidate(digest, seqid, start, end, strand, cds))
    return candidates


def filter_candidate_query_paf(
    *,
    candidate_gff_path: str | Path,
    paf_inputs: Iterable[tuple[str, str | Path]],
    output_paf_path: str | Path,
    output_count_tsv_path: str | Path,
    output_summary_tsv_path: str | Path,
    output_manifest_json_path: str | Path,
    alignment_strand_source: str = "query_orientation",
) -> dict[str, Any]:
    """Retain all alignments for queries that can overlap a candidate.

    Selection is made using candidate strand and span, but once a query is
    selected every one of its alignments is retained so downstream uniqueness
    and near-best checks remain lossless.  Prefixing query names by the input
    accession makes independently generated read identifiers unambiguous.
    """

    if alignment_strand_source not in ALIGNMENT_STRAND_SOURCES:
        raise ValueError(
            "alignment_strand_source must be one of: "
            + ",".join(ALIGNMENT_STRAND_SOURCES)
        )
    candidate_gff = Path(candidate_gff_path)
    output_paf = Path(output_paf_path)
    output_counts = Path(output_count_tsv_path)
    output_summary = Path(output_summary_tsv_path)
    output_manifest = Path(output_manifest_json_path)
    inputs: list[tuple[str, Path]] = []
    accessions: set[str] = set()
    for accession, raw_path in paf_inputs:
        path = Path(raw_path)
        if not re.fullmatch(r"[A-Za-z0-9_.-]+", accession) or "|" in accession:
            raise ValueError(f"Unsafe PAF accession: {accession!r}")
        if accession in accessions:
            raise ValueError(f"Duplicate PAF accession: {accession}")
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f"Missing or empty PAF input: {path}")
        accessions.add(accession)
        inputs.append((accession, path))
    if not candidate_gff.is_file() or candidate_gff.stat().st_size == 0:
        raise ValueError(f"Missing or empty candidate GFF: {candidate_gff}")
    if not inputs:
        raise ValueError("At least one PAF input is required")
    _refuse_outputs((output_paf, output_counts, output_summary, output_manifest))

    bin_bp = 1_000_000
    candidate_bins: dict[tuple[str, str, int], list[tuple[int, int]]] = defaultdict(list)
    candidates = _read_candidates(candidate_gff)
    for candidate in candidates:
        for bin_index in range(
            (candidate.start - 1) // bin_bp,
            (candidate.end - 1) // bin_bp + 1,
        ):
            candidate_bins[(candidate.seqid, candidate.strand, bin_index)].append(
                (candidate.start, candidate.end)
            )

    def can_overlap_candidate(fields: list[str], strand: str) -> bool:
        start, end = int(fields[7]) + 1, int(fields[8])
        for bin_index in range((start - 1) // bin_bp, (end - 1) // bin_bp + 1):
            for candidate_start, candidate_end in candidate_bins.get(
                (fields[5], strand, bin_index), ()
            ):
                if start <= candidate_end and end >= candidate_start:
                    return True
        return False

    output_paf.parent.mkdir(parents=True, exist_ok=True)
    output_counts.parent.mkdir(parents=True, exist_ok=True)
    output_summary.parent.mkdir(parents=True, exist_ok=True)
    output_manifest.parent.mkdir(parents=True, exist_ok=True)
    count_rows: list[tuple[str, int]] = []
    summary_rows: list[dict[str, Any]] = []
    input_manifests: list[dict[str, Any]] = []
    with output_paf.open("x", encoding="utf-8", newline="") as destination:
        for accession, paf in inputs:
            mapped_queries: set[str] = set()
            retained_queries: set[str] = set()
            raw_alignments = 0
            strand_available = 0
            strand_unavailable = 0
            source_sha256 = hashlib.sha256()
            with paf.open("r", encoding="utf-8", newline="") as source:
                for line_number, raw_line in enumerate(source, start=1):
                    source_sha256.update(raw_line.encode("utf-8"))
                    fields = raw_line.rstrip("\r\n").split("\t")
                    if len(fields) < 12 or not fields[0]:
                        raise ValueError(
                            f"Malformed {accession} PAF line {line_number}"
                        )
                    raw_alignments += 1
                    mapped_queries.add(fields[0])
                    tags = _parse_paf_tags(fields)
                    strand = _reference_transcript_strand(
                        fields[4], tags.get("ts"), alignment_strand_source
                    )
                    if strand is None:
                        strand_unavailable += 1
                        continue
                    strand_available += 1
                    if can_overlap_candidate(fields, strand):
                        retained_queries.add(fields[0])
            retained_alignments = 0
            with paf.open("r", encoding="utf-8", newline="") as source:
                for line_number, raw_line in enumerate(source, start=1):
                    fields = raw_line.rstrip("\r\n").split("\t")
                    if len(fields) < 12:
                        raise ValueError(
                            f"Malformed {accession} PAF line {line_number}"
                        )
                    if fields[0] not in retained_queries:
                        continue
                    fields[0] = f"{accession}|{fields[0]}"
                    destination.write("\t".join(fields) + "\n")
                    retained_alignments += 1
            count_rows.extend(
                (f"{accession}|{query}", 1) for query in retained_queries
            )
            summary_rows.append(
                {
                    "accession": accession,
                    "mapped_queries": len(mapped_queries),
                    "retained_queries": len(retained_queries),
                    "raw_alignments": raw_alignments,
                    "retained_alignments": retained_alignments,
                    "strand_available_alignments": strand_available,
                    "strand_unavailable_alignments": strand_unavailable,
                }
            )
            input_manifests.append(
                {
                    "accession": accession,
                    "file_name": paf.name,
                    "bytes": paf.stat().st_size,
                    "sha256": source_sha256.hexdigest(),
                }
            )
    if not count_rows or output_paf.stat().st_size == 0:
        raise ValueError("Candidate query filter retained no mapped reads")

    with output_counts.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("transcript_id", "full_length_reads"))
        writer.writerows(sorted(count_rows))
    with output_summary.open("x", encoding="utf-8", newline="") as handle:
        fieldnames = tuple(summary_rows[0])
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(summary_rows)

    manifest: dict[str, Any] = {
        "schema_version": ISOSEQ_QUERY_FILTER_SCHEMA,
        "generator": {"name": "PloidyPatch", "version": __version__},
        "inputs": {
            "candidate_gff": {
                "file_name": candidate_gff.name,
                "sha256": _file_sha256(candidate_gff),
            },
            "paf": input_manifests,
        },
        "parameters": {
            "alignment_strand_source": alignment_strand_source,
            "minimap2_ts_derivation": (
                "reference_transcript_strand=paf_query_target_strand*"
                "ts_query_transcript_relation"
                if alignment_strand_source == "minimap2_ts"
                else None
            ),
            "lossless_all_alignments_for_retained_query": True,
            "query_prefix": "accession|query_id",
            "candidate_bin_bp": bin_bp,
        },
        "counts": {
            "candidate_models": len(candidates),
            "mapped_queries": sum(row["mapped_queries"] for row in summary_rows),
            "retained_queries": len(count_rows),
            "raw_alignments": sum(row["raw_alignments"] for row in summary_rows),
            "retained_alignments": sum(
                row["retained_alignments"] for row in summary_rows
            ),
            "strand_available_alignments": sum(
                row["strand_available_alignments"] for row in summary_rows
            ),
            "strand_unavailable_alignments": sum(
                row["strand_unavailable_alignments"] for row in summary_rows
            ),
        },
        "outputs": {
            "paf": {
                "file_name": output_paf.name,
                "bytes": output_paf.stat().st_size,
                "sha256": _file_sha256(output_paf),
            },
            "counts": {
                "file_name": output_counts.name,
                "sha256": _file_sha256(output_counts),
            },
            "summary": {
                "file_name": output_summary.name,
                "sha256": _file_sha256(output_summary),
            },
        },
    }
    with output_manifest.open("x", encoding="utf-8", newline="") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return manifest


def _parse_cigar(cigar: str, target_start_zero: int) -> tuple[tuple[tuple[int, int], ...], tuple[tuple[int, int], ...]]:
    tokens = [(int(length), operation) for length, operation in CIGAR_TOKEN.findall(cigar)]
    if not tokens or "".join(f"{length}{operation}" for length, operation in tokens) != cigar:
        raise ValueError(f"Malformed minimap2 cg CIGAR: {cigar}")
    reference_cursor = target_start_zero
    exon_start = reference_cursor + 1
    blocks: list[tuple[int, int]] = []
    junctions: list[tuple[int, int]] = []
    for length, operation in tokens:
        if operation in {"M", "D", "=", "X"}:
            reference_cursor += length
        elif operation == "N":
            if reference_cursor >= exon_start:
                blocks.append((exon_start, reference_cursor))
            junctions.append((reference_cursor, reference_cursor + length + 1))
            reference_cursor += length
            exon_start = reference_cursor + 1
    if reference_cursor >= exon_start:
        blocks.append((exon_start, reference_cursor))
    return _merge_intervals(blocks), tuple(junctions)


def _read_paf(
    path: Path, alignment_strand_source: str
) -> tuple[dict[str, list[Alignment]], Counter[str]]:
    grouped: dict[str, list[Alignment]] = defaultdict(list)
    strand_counts: Counter[str] = Counter()
    with path.open("r", encoding="utf-8", newline="") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            fields = raw_line.rstrip("\r\n").split("\t")
            if len(fields) < 12:
                raise ValueError(f"Malformed PAF line {line_number}")
            tags = _parse_paf_tags(fields)
            if "cg" not in tags or "AS" not in tags:
                raise ValueError(f"PAF line {line_number} lacks cg or AS tag")
            strand = _reference_transcript_strand(
                fields[4], tags.get("ts"), alignment_strand_source
            )
            strand_counts[
                "available" if strand is not None else "unavailable"
            ] += 1
            query_length = int(fields[1])
            aligned_query = int(fields[3]) - int(fields[2])
            block_length = int(fields[10])
            blocks, junctions = _parse_cigar(tags["cg"], int(fields[7]))
            grouped[fields[0]].append(
                Alignment(
                    query=fields[0],
                    query_length=query_length,
                    seqid=fields[5],
                    start=int(fields[7]) + 1,
                    end=int(fields[8]),
                    strand=strand,
                    mapq=int(fields[11]),
                    query_coverage=aligned_query / query_length if query_length else 0.0,
                    identity=int(fields[9]) / block_length if block_length else 0.0,
                    alignment_score=int(tags["AS"]),
                    primary=tags.get("tp", "") == "P",
                    exonic_blocks=blocks,
                    junctions=junctions,
                )
            )
    return grouped, strand_counts


def _overlap_bp(left: Iterable[tuple[int, int]], right: Iterable[tuple[int, int]]) -> int:
    left_intervals = tuple(left)
    right_intervals = tuple(right)
    i = j = total = 0
    while i < len(left_intervals) and j < len(right_intervals):
        left_start, left_end = left_intervals[i]
        right_start, right_end = right_intervals[j]
        total += max(0, min(left_end, right_end) - max(left_start, right_start) + 1)
        if left_end <= right_end:
            i += 1
        else:
            j += 1
    return total


def _junction_is_canonical(index: FastaIndex, seqid: str, strand: str, junction: tuple[int, int]) -> bool:
    donor_end, acceptor_start = junction
    motif = index.fetch(seqid, donor_end + 1, donor_end + 2) + index.fetch(seqid, acceptor_start - 2, acceptor_start - 1)
    return motif in (CANONICAL_PLUS if strand == "+" else CANONICAL_MINUS)


def validate_isoseq_candidate_chains(
    *,
    candidate_gff_path: str | Path,
    paf_path: str | Path,
    selected_count_tsv_path: str | Path,
    genome_fasta_path: str | Path,
    output_evidence_tsv_path: str | Path,
    minimum_query_coverage: float = 0.90,
    minimum_identity: float = 0.98,
    minimum_mapq: int = 20,
    maximum_secondary_score_fraction: float = 0.95,
    minimum_candidate_cds_coverage: float = 0.90,
    flank_bp: int = 5000,
    alignment_strand_source: str = "query_orientation",
) -> dict[str, Any]:
    """Validate frozen coding candidates with independently mapped full-length RNA."""

    if not 0 < minimum_query_coverage <= 1 or not 0 < minimum_identity <= 1:
        raise ValueError("Alignment coverage and identity thresholds must be in (0,1]")
    if not 0 < maximum_secondary_score_fraction <= 1:
        raise ValueError("maximum_secondary_score_fraction must be in (0,1]")
    if not 0 < minimum_candidate_cds_coverage <= 1 or flank_bp < 0:
        raise ValueError("Invalid candidate coverage or flank threshold")
    if alignment_strand_source not in ALIGNMENT_STRAND_SOURCES:
        raise ValueError(
            "alignment_strand_source must be one of: "
            + ",".join(ALIGNMENT_STRAND_SOURCES)
        )
    candidate_gff = Path(candidate_gff_path)
    paf = Path(paf_path)
    selected_counts = Path(selected_count_tsv_path)
    genome = Path(genome_fasta_path)
    fai = Path(str(genome) + ".fai")
    output = Path(output_evidence_tsv_path)
    manifest_path = Path(str(output) + ".manifest.json")
    for required in (candidate_gff, paf, selected_counts, genome, fai):
        if not required.is_file() or required.stat().st_size == 0:
            raise ValueError(f"Missing or empty validation input: {required}")
    _refuse_outputs((output, manifest_path))

    b73_counts: dict[str, int] = {}
    with selected_counts.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or "transcript_id" not in reader.fieldnames:
            raise ValueError("Selected transcript counts lack transcript_id")
        if "full_length_reads" in reader.fieldnames:
            count_field = "full_length_reads"
        elif "b73_full_length_reads" in reader.fieldnames:
            count_field = "b73_full_length_reads"
        else:
            raise ValueError(
                "Selected transcript counts lack full_length_reads or "
                "b73_full_length_reads"
            )
        for row in reader:
            transcript_id = row["transcript_id"]
            if transcript_id in b73_counts:
                raise ValueError("Duplicate selected transcript count")
            count = int(row[count_field])
            if count < 1:
                raise ValueError("Selected full-length read counts must be positive")
            b73_counts[transcript_id] = count
    candidates = _read_candidates(candidate_gff)
    paf_groups, paf_strand_counts = _read_paf(paf, alignment_strand_source)
    unknown_queries = set(paf_groups) - set(b73_counts)
    if unknown_queries:
        raise ValueError(f"PAF contains unselected transcript queries: {len(unknown_queries)}")

    alignment_reasons: Counter[str] = Counter()
    preliminary: list[Alignment] = []
    for query in sorted(b73_counts):
        alignments = paf_groups.get(query, [])
        primaries = [alignment for alignment in alignments if alignment.primary]
        if not primaries:
            alignment_reasons["no_primary_alignment"] += 1
            continue
        primary = max(primaries, key=lambda alignment: (alignment.alignment_score, alignment.mapq))
        other_scores = [alignment.alignment_score for alignment in alignments if alignment is not primary]
        if primary.strand is None:
            alignment_reasons["transcript_strand_unavailable"] += 1
        elif primary.query_coverage < minimum_query_coverage:
            alignment_reasons["query_coverage_below_threshold"] += 1
        elif primary.identity < minimum_identity:
            alignment_reasons["identity_below_threshold"] += 1
        elif primary.mapq < minimum_mapq:
            alignment_reasons["mapq_below_threshold"] += 1
        elif other_scores and max(other_scores) >= maximum_secondary_score_fraction * primary.alignment_score:
            alignment_reasons["near_best_secondary_alignment"] += 1
        else:
            preliminary.append(primary)

    index = FastaIndex(genome, fai)
    try:
        junction_read_support: Counter[tuple[str, str, int, int]] = Counter()
        for alignment in preliminary:
            assert alignment.strand is not None
            for donor_end, acceptor_start in set(alignment.junctions):
                junction_read_support[(alignment.seqid, alignment.strand, donor_end, acceptor_start)] += b73_counts[alignment.query]
        usable: list[Alignment] = []
        for alignment in preliminary:
            assert alignment.strand is not None
            unsupported_noncanonical = [
                junction
                for junction in alignment.junctions
                if not _junction_is_canonical(index, alignment.seqid, alignment.strand, junction)
                and junction_read_support[(alignment.seqid, alignment.strand, *junction)] < 2
            ]
            if unsupported_noncanonical:
                alignment_reasons["unsupported_noncanonical_junction"] += 1
            else:
                usable.append(alignment)
        alignment_reasons["usable"] = len(usable)

        bin_bp = 1_000_000
        alignment_bins: dict[tuple[str, str, int], list[int]] = defaultdict(list)
        for alignment_index, alignment in enumerate(usable):
            for bin_index in range((alignment.start - 1) // bin_bp, (alignment.end - 1) // bin_bp + 1):
                alignment_bins[(alignment.seqid, alignment.strand, bin_index)].append(alignment_index)

        rows: list[dict[str, Any]] = []
        state_counts: Counter[str] = Counter()
        for candidate in sorted(candidates, key=lambda item: item.digest):
            flank_sequence = index.fetch(candidate.seqid, candidate.start - flank_bp, candidate.end + flank_bp)
            candidate_junctions = set(candidate.junctions)
            noncanonical_candidate_junctions = sum(
                not _junction_is_canonical(index, candidate.seqid, candidate.strand, junction)
                for junction in candidate_junctions
            )
            nearby_indices: set[int] = set()
            for bin_index in range((candidate.start - 1) // bin_bp, (candidate.end - 1) // bin_bp + 1):
                nearby_indices.update(alignment_bins.get((candidate.seqid, candidate.strand, bin_index), ()))
            overlaps: list[tuple[Alignment, float, set[tuple[int, int]]]] = []
            for alignment_index in sorted(nearby_indices):
                alignment = usable[alignment_index]
                if alignment.end < candidate.start or alignment.start > candidate.end:
                    continue
                coverage = _overlap_bp(candidate.cds, alignment.exonic_blocks) / candidate.cds_bp
                matched = candidate_junctions.intersection(alignment.junctions)
                overlaps.append((alignment, coverage, matched))
            best_coverage = max((coverage for _, coverage, _ in overlaps), default=0.0)
            matching_union: set[tuple[int, int]] = set()
            for _, _, matched in overlaps:
                matching_union.update(matched)

            supporting: list[Alignment] = []
            if "N" in flank_sequence:
                state = "not_assessable"
            elif not candidate_junctions:
                supporting = [alignment for alignment, coverage, _ in overlaps if coverage >= minimum_candidate_cds_coverage]
                state = "single_exon_span_support" if supporting else "no_qualifying_observation"
            else:
                full = [
                    alignment
                    for alignment, coverage, matched in overlaps
                    if coverage >= minimum_candidate_cds_coverage and matched == candidate_junctions
                ]
                if full:
                    supporting = full
                    state = "full_chain_supported"
                elif matching_union == candidate_junctions:
                    supporting = [alignment for alignment, _, matched in overlaps if matched]
                    state = "all_junctions_supported"
                elif matching_union:
                    supporting = [alignment for alignment, _, matched in overlaps if matched]
                    state = "partial_junction_support"
                else:
                    state = "no_qualifying_observation"
            state_counts[state] += 1
            supporting_ids = sorted({alignment.query for alignment in supporting})
            supporting_reads = sum(b73_counts[query] for query in supporting_ids)
            supporting_groups = sorted(
                {
                    query.split("|", 1)[0]
                    for query in supporting_ids
                    if "|" in query and query.split("|", 1)[0]
                }
            )
            rows.append(
                {
                    "candidate_digest": candidate.digest,
                    "seqid": candidate.seqid,
                    "start": candidate.start,
                    "end": candidate.end,
                    "strand": candidate.strand,
                    "cds_segments": len(candidate.cds),
                    "cds_bp": candidate.cds_bp,
                    "candidate_junctions": len(candidate_junctions),
                    "noncanonical_candidate_junctions": noncanonical_candidate_junctions,
                    "qualifying_overlapping_transcripts": len(overlaps),
                    "best_candidate_cds_coverage": f"{best_coverage:.12g}",
                    "matched_candidate_junctions": len(matching_union),
                    "evidence_state": state,
                    "supporting_transcripts": ",".join(supporting_ids),
                    "supporting_full_length_reads": supporting_reads,
                    "supporting_b73_full_length_reads": supporting_reads,
                    "supporting_query_groups": ",".join(supporting_groups),
                }
            )
    finally:
        index.close()

    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=tuple(rows[0]) if rows else (), delimiter="\t", lineterminator="\n")
        if rows:
            writer.writeheader()
            writer.writerows(rows)
    manifest: dict[str, Any] = {
        "schema_version": ISOSEQ_VALIDATION_SCHEMA,
        "generator": {"name": "PloidyPatch", "version": __version__},
        "inputs": {
            "candidate_gff": {"file_name": candidate_gff.name, "sha256": _file_sha256(candidate_gff)},
            "paf": {"file_name": paf.name, "sha256": _file_sha256(paf)},
            "selected_counts": {"file_name": selected_counts.name, "sha256": _file_sha256(selected_counts)},
            "genome_fasta": {"file_name": genome.name, "sha256": _file_sha256(genome)},
            "genome_fai": {"file_name": fai.name, "sha256": _file_sha256(fai)},
        },
        "parameters": {
            "minimum_query_coverage": minimum_query_coverage,
            "minimum_identity": minimum_identity,
            "minimum_mapq": minimum_mapq,
            "maximum_secondary_score_fraction": maximum_secondary_score_fraction,
            "minimum_candidate_cds_coverage": minimum_candidate_cds_coverage,
            "flank_bp": flank_bp,
            "single_exon_is_not_full_chain_positive": True,
            "selected_count_field": count_field,
            "alignment_strand_source": alignment_strand_source,
            "minimap2_ts_derivation": (
                "reference_transcript_strand=paf_query_target_strand*"
                "ts_query_transcript_relation"
                if alignment_strand_source == "minimap2_ts"
                else None
            ),
        },
        "counts": {
            "selected_transcripts": len(b73_counts),
            "paf_queries": len(paf_groups),
            "candidate_models": len(candidates),
            "alignment_decisions": dict(sorted(alignment_reasons.items())),
            "paf_strand_availability": dict(sorted(paf_strand_counts.items())),
            "evidence_states": dict(sorted(state_counts.items())),
        },
        "output": {"file_name": output.name, "sha256": _file_sha256(output), "rows": len(rows)},
    }
    with manifest_path.open("x", encoding="utf-8", newline="") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return manifest


def join_isoseq_review_rankings(
    *,
    evidence_tsv_path: str | Path,
    review_rankings_tsv_path: str | Path,
    output_tsv_path: str | Path,
    output_summary_json_path: str | Path,
    review_budgets: Iterable[int] = (25, 50, 100, 200),
    comparator_estimator: str = "baseline",
    primary_estimator: str = "topology",
) -> dict[str, Any]:
    evidence_path = Path(evidence_tsv_path)
    rankings_path = Path(review_rankings_tsv_path)
    output = Path(output_tsv_path)
    summary_path = Path(output_summary_json_path)
    budgets = tuple(sorted(set(review_budgets)))
    if not budgets or any(value < 1 for value in budgets):
        raise ValueError("review_budgets must contain positive values")
    for required in (evidence_path, rankings_path):
        if not required.is_file() or required.stat().st_size == 0:
            raise ValueError(f"Missing review evidence input: {required}")
    _refuse_outputs((output, summary_path))
    with evidence_path.open("r", encoding="utf-8", newline="") as handle:
        evidence_rows = list(csv.DictReader(handle, delimiter="\t"))
    evidence = {row["candidate_digest"]: row for row in evidence_rows}
    if len(evidence) != len(evidence_rows):
        raise ValueError("Evidence table contains duplicate candidate digests")
    with rankings_path.open("r", encoding="utf-8", newline="") as handle:
        ranking_rows = list(csv.DictReader(handle, delimiter="\t"))
    missing = {row["candidate_digest"] for row in ranking_rows} - set(evidence)
    if missing:
        raise ValueError(f"Rankings contain candidates absent from evidence: {len(missing)}")
    joined_rows = [{**row, **evidence[row["candidate_digest"]]} for row in ranking_rows]
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=tuple(joined_rows[0]) if joined_rows else (), delimiter="\t", lineterminator="\n")
        if joined_rows:
            writer.writeheader()
            writer.writerows(joined_rows)

    estimators = sorted({row["estimator"] for row in ranking_rows})
    if (
        not comparator_estimator
        or not primary_estimator
        or comparator_estimator == primary_estimator
    ):
        raise ValueError("Review comparison must name two distinct estimators")
    missing_estimators = {comparator_estimator, primary_estimator} - set(estimators)
    if missing_estimators:
        raise ValueError(
            "Review rankings lack requested estimator(s): "
            + ",".join(sorted(missing_estimators))
        )
    states = (
        "full_chain_supported", "all_junctions_supported", "partial_junction_support",
        "single_exon_span_support", "no_qualifying_observation", "not_assessable",
    )
    yield_rows: list[dict[str, Any]] = []
    for estimator in estimators:
        estimator_rows = [row for row in joined_rows if row["estimator"] == estimator]
        for budget in budgets:
            selected = [row for row in estimator_rows if int(row["review_rank"]) <= budget]
            counts = Counter(row["evidence_state"] for row in selected)
            yield_rows.append({"estimator": estimator, "review_budget": budget, **{state: counts[state] for state in states}})
    primary_budget = 100
    top_sets = {
        estimator: {
            row["candidate_digest"]
            for row in ranking_rows
            if row["estimator"] == estimator and int(row["review_rank"]) <= primary_budget
        }
        for estimator in estimators
    }
    primary_yield = {
        row["estimator"]: row["full_chain_supported"]
        for row in yield_rows if row["review_budget"] == primary_budget
    }
    estimator_delta = {
        "primary_estimator": primary_estimator,
        "comparator_estimator": comparator_estimator,
        "delta_full_chain_supported": (
            primary_yield[primary_estimator] - primary_yield[comparator_estimator]
        ),
    }
    primary_report: dict[str, Any] = {
        "full_chain_supported": primary_yield,
        "top_set_overlap": len(set.intersection(*top_sets.values())) if top_sets else 0,
        "estimator_delta": estimator_delta,
    }
    if primary_estimator == "topology" and comparator_estimator == "baseline":
        primary_report["topology_minus_baseline"] = estimator_delta[
            "delta_full_chain_supported"
        ]
    summary: dict[str, Any] = {
        "schema_version": "ploidypatch.isoseq_review_yield.v1",
        "inputs": {
            "evidence_sha256": _file_sha256(evidence_path),
            "review_rankings_sha256": _file_sha256(rankings_path),
        },
        "parameters": {
            "review_budgets": list(budgets),
            "primary_budget": primary_budget,
            "primary_state": "full_chain_supported",
            "comparator_estimator": comparator_estimator,
            "primary_estimator": primary_estimator,
        },
        "yield": yield_rows,
        "primary": primary_report,
        "output": {"file_name": output.name, "sha256": _file_sha256(output), "rows": len(joined_rows)},
    }
    with summary_path.open("x", encoding="utf-8", newline="") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return summary


def _quantile(values: list[float], probability: float) -> float:
    ordered = sorted(values)
    position = probability * (len(ordered) - 1)
    lower = int(position)
    upper = min(lower + 1, len(ordered) - 1)
    fraction = position - lower
    return ordered[lower] * (1 - fraction) + ordered[upper] * fraction


def bootstrap_isoseq_review_yield(
    *,
    evidence_tsv_path: str | Path,
    review_rankings_tsv_path: str | Path,
    output_json_path: str | Path,
    review_budgets: Iterable[int] = (25, 50, 100, 200),
    replicates: int = 20_000,
    seed: int = 20260808,
    alpha: float = 0.05,
    comparator_estimator: str = "baseline",
    primary_estimator: str = "topology",
) -> dict[str, Any]:
    """Evaluate top-K evidence yield with chromosome-aware resampling."""

    if replicates < 100 or not 0 < alpha < 1:
        raise ValueError("Invalid bootstrap replicate count or alpha")
    budgets = tuple(sorted(set(review_budgets)))
    if not budgets or any(value < 1 for value in budgets):
        raise ValueError("review_budgets must contain positive values")
    evidence_path = Path(evidence_tsv_path)
    rankings_path = Path(review_rankings_tsv_path)
    output = Path(output_json_path)
    for required in (evidence_path, rankings_path):
        if not required.is_file() or required.stat().st_size == 0:
            raise ValueError(f"Missing Iso-Seq bootstrap input: {required}")
    _refuse_outputs((output,))
    with evidence_path.open("r", encoding="utf-8", newline="") as handle:
        evidence_rows = list(csv.DictReader(handle, delimiter="\t"))
    evidence = {row["candidate_digest"]: row for row in evidence_rows}
    if len(evidence) != len(evidence_rows):
        raise ValueError("Evidence table contains duplicate candidate digests")
    with rankings_path.open("r", encoding="utf-8", newline="") as handle:
        ranking_rows = list(csv.DictReader(handle, delimiter="\t"))
    if {row["candidate_digest"] for row in ranking_rows} - set(evidence):
        raise ValueError("Review rankings contain candidates absent from evidence")
    estimators = sorted({row["estimator"] for row in ranking_rows})
    if (
        not comparator_estimator
        or not primary_estimator
        or comparator_estimator == primary_estimator
    ):
        raise ValueError("Bootstrap comparison must name two distinct estimators")
    missing_estimators = {comparator_estimator, primary_estimator} - set(estimators)
    if missing_estimators:
        raise ValueError(
            "Review rankings lack requested estimator(s): "
            + ",".join(sorted(missing_estimators))
        )
    seqids = sorted({row["seqid"] for row in evidence_rows})
    primary = {
        digest: row["evidence_state"] == "full_chain_supported"
        for digest, row in evidence.items()
    }
    candidates_by_seqid: dict[str, list[str]] = defaultdict(list)
    for digest, row in evidence.items():
        candidates_by_seqid[row["seqid"]].append(digest)
    for digests in candidates_by_seqid.values():
        digests.sort()

    top_sets: dict[tuple[str, int], set[str]] = {}
    for estimator in estimators:
        estimator_rows = [row for row in ranking_rows if row["estimator"] == estimator]
        for budget in budgets:
            selected = {
                row["candidate_digest"]
                for row in estimator_rows
                if int(row["review_rank"]) <= budget
            }
            if len(selected) != min(budget, len(estimator_rows)):
                raise ValueError(f"Non-unique or incomplete top-{budget} for {estimator}")
            top_sets[(estimator, budget)] = selected

    rng = random.Random(seed)
    random_reports: list[dict[str, Any]] = []
    lower_probability = alpha / 2
    upper_probability = 1 - alpha / 2
    for estimator in estimators:
        for budget in budgets:
            selected = top_sets[(estimator, budget)]
            observed = sum(primary[digest] for digest in selected)
            allocation = Counter(evidence[digest]["seqid"] for digest in selected)
            samples: list[int] = []
            for _ in range(replicates):
                support = 0
                for seqid, count in sorted(allocation.items()):
                    support += sum(
                        primary[digest]
                        for digest in rng.sample(candidates_by_seqid[seqid], count)
                    )
                samples.append(support)
            expected = sum(samples) / replicates
            random_reports.append(
                {
                    "estimator": estimator,
                    "review_budget": budget,
                    "observed_full_chain_supported": observed,
                    "chromosome_allocation": dict(sorted(allocation.items())),
                    "random_mean": expected,
                    "random_ci_lower": _quantile(samples, lower_probability),
                    "random_ci_upper": _quantile(samples, upper_probability),
                    "fold_enrichment_over_random_mean": observed / expected if expected else None,
                    "empirical_p_random_ge_observed": (1 + sum(value >= observed for value in samples)) / (replicates + 1),
                }
            )

    primary_budget = 100
    if primary_budget not in budgets:
        raise ValueError("The predeclared primary review budget 100 is required")
    support_by_estimator_seqid: dict[str, Counter[str]] = {}
    for estimator in (comparator_estimator, primary_estimator):
        counter: Counter[str] = Counter()
        for digest in top_sets[(estimator, primary_budget)]:
            if primary[digest]:
                counter[evidence[digest]["seqid"]] += 1
        support_by_estimator_seqid[estimator] = counter
    observed_delta = sum(support_by_estimator_seqid[primary_estimator].values()) - sum(
        support_by_estimator_seqid[comparator_estimator].values()
    )
    delta_samples = [
        sum(
            support_by_estimator_seqid[primary_estimator][sampled_seqid]
            - support_by_estimator_seqid[comparator_estimator][sampled_seqid]
            for sampled_seqid in (rng.choice(seqids) for _ in seqids)
        )
        for _ in range(replicates)
    ]
    report: dict[str, Any] = {
        "schema_version": "ploidypatch.isoseq_review_bootstrap.v2",
        "inputs": {
            "evidence_sha256": _file_sha256(evidence_path),
            "review_rankings_sha256": _file_sha256(rankings_path),
        },
        "parameters": {
            "review_budgets": list(budgets),
            "primary_budget": primary_budget,
            "primary_state": "full_chain_supported",
            "replicates": replicates,
            "seed": seed,
            "alpha": alpha,
            "random_null": "sample_without_replacement_within_observed_top_k_chromosome_allocation",
            "paired_delta_resampling_unit": "target_chromosome",
            "comparator_estimator": comparator_estimator,
            "primary_estimator": primary_estimator,
        },
        "counts": {
            "candidates": len(evidence),
            "full_chain_supported": sum(primary.values()),
            "target_chromosomes": len(seqids),
        },
        "random_ranking_enrichment": random_reports,
        "primary_estimator_delta": {
            "primary_estimator": primary_estimator,
            "comparator_estimator": comparator_estimator,
            "observed_delta_supported_per_100": observed_delta,
            "ci_lower": _quantile(delta_samples, lower_probability),
            "ci_upper": _quantile(delta_samples, upper_probability),
            "probability_delta_gt_zero": sum(value > 0 for value in delta_samples) / replicates,
            "probability_delta_eq_zero": sum(value == 0 for value in delta_samples) / replicates,
            "comparator_support_by_chromosome": dict(
                sorted(support_by_estimator_seqid[comparator_estimator].items())
            ),
            "primary_support_by_chromosome": dict(
                sorted(support_by_estimator_seqid[primary_estimator].items())
            ),
        },
    }
    if primary_estimator == "topology" and comparator_estimator == "baseline":
        report["primary_topology_minus_baseline"] = {
            **report["primary_estimator_delta"],
            "baseline_support_by_chromosome": report["primary_estimator_delta"][
                "comparator_support_by_chromosome"
            ],
            "topology_support_by_chromosome": report["primary_estimator_delta"][
                "primary_support_by_chromosome"
            ],
        }
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x", encoding="utf-8", newline="") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return report
