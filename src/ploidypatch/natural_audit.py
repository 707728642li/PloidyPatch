from __future__ import annotations

import csv
import hashlib
import json
import math
from bisect import bisect_right
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

from . import __version__
from .baseline import _file_sha256
from .consensus import _merge_intervals
from .gff import parse_attributes
from .io import open_text
from .isoseq_validation import Candidate, FastaIndex, _overlap_bp, _read_candidates


CSV_FIELD_SIZE_LIMIT = 2**31 - 1
csv.field_size_limit(CSV_FIELD_SIZE_LIMIT)


NATURAL_CDS_EXPORT_SCHEMA = "ploidypatch.natural_candidate_cds_export.v1"
NATURAL_AUDIT_SCHEMA = "ploidypatch.natural_candidate_biological_audit.v1"
STOP_CODONS = frozenset({"TAA", "TAG", "TGA"})
DNA_COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


@dataclass(frozen=True)
class FeatureInterval:
    seqid: str
    start: int
    end: int
    strand: str
    identifier: str


@dataclass(frozen=True)
class SelfAlignment:
    query: str
    query_length: int
    query_start: int
    query_end: int
    strand: str
    seqid: str
    start: int
    end: int
    matches: int
    block_length: int
    mapq: int
    alignment_score: int

    @property
    def query_coverage(self) -> float:
        return (
            (self.query_end - self.query_start) / self.query_length
            if self.query_length
            else 0.0
        )

    @property
    def identity(self) -> float:
        return self.matches / self.block_length if self.block_length else 0.0


class FeatureIntervalIndex:
    def __init__(self, intervals: Iterable[FeatureInterval]):
        grouped: dict[str, list[FeatureInterval]] = defaultdict(list)
        for interval in intervals:
            grouped[interval.seqid].append(interval)
        self.intervals: dict[str, tuple[FeatureInterval, ...]] = {}
        self.starts: dict[str, tuple[int, ...]] = {}
        self.prefix_max_end: dict[str, tuple[int, ...]] = {}
        for seqid, values in grouped.items():
            ordered = tuple(
                sorted(values, key=lambda item: (item.start, item.end, item.identifier))
            )
            maximum = 0
            prefix: list[int] = []
            for interval in ordered:
                maximum = max(maximum, interval.end)
                prefix.append(maximum)
            self.intervals[seqid] = ordered
            self.starts[seqid] = tuple(item.start for item in ordered)
            self.prefix_max_end[seqid] = tuple(prefix)

    def overlaps(self, seqid: str, start: int, end: int) -> list[FeatureInterval]:
        intervals = self.intervals.get(seqid, ())
        if not intervals:
            return []
        starts = self.starts[seqid]
        prefix = self.prefix_max_end[seqid]
        index = bisect_right(starts, end) - 1
        output: list[FeatureInterval] = []
        while index >= 0 and prefix[index] >= start:
            interval = intervals[index]
            if interval.end >= start:
                output.append(interval)
            index -= 1
        return output

    def nearest_distance(self, seqid: str, start: int, end: int) -> int | None:
        intervals = self.intervals.get(seqid, ())
        if not intervals:
            return None
        starts = self.starts[seqid]
        prefix = self.prefix_max_end[seqid]
        right = bisect_right(starts, end)
        left_distance = None
        if right:
            maximum_end = prefix[right - 1]
            if maximum_end >= start:
                return 0
            left_distance = start - maximum_end - 1
        right_distance = (
            intervals[right].start - end - 1 if right < len(intervals) else None
        )
        distances = [
            distance
            for distance in (left_distance, right_distance)
            if distance is not None
        ]
        return min(distances) if distances else None


def _refuse_outputs(paths: Iterable[Path]) -> None:
    collisions = [str(path) for path in paths if path.exists()]
    if collisions:
        raise FileExistsError("Refusing to overwrite output(s): " + ", ".join(collisions))


def _reverse_complement(sequence: str) -> str:
    return sequence.translate(DNA_COMPLEMENT)[::-1]


def _candidate_cds(index: FastaIndex, candidate: Candidate) -> str:
    sequence = "".join(
        index.fetch(candidate.seqid, start, end) for start, end in candidate.cds
    )
    return sequence if candidate.strand == "+" else _reverse_complement(sequence)


def export_natural_candidate_cds(
    *,
    candidate_gff_path: str | Path,
    genome_fasta_path: str | Path,
    output_fasta_path: str | Path,
) -> dict[str, Any]:
    candidate_gff = Path(candidate_gff_path)
    genome = Path(genome_fasta_path)
    fai = Path(str(genome) + ".fai")
    output = Path(output_fasta_path)
    manifest_path = Path(str(output) + ".manifest.json")
    for required in (candidate_gff, genome, fai):
        if not required.is_file() or required.stat().st_size == 0:
            raise ValueError(f"Missing or empty candidate CDS input: {required}")
    _refuse_outputs((output, manifest_path))
    candidates = sorted(_read_candidates(candidate_gff), key=lambda item: item.digest)
    if not candidates:
        raise ValueError("No PloidyPatchConsensus candidates in GFF")
    index = FastaIndex(genome, fai)
    output.parent.mkdir(parents=True, exist_ok=True)
    total_bases = 0
    try:
        with output.open("x", encoding="ascii", newline="") as handle:
            for candidate in candidates:
                sequence = _candidate_cds(index, candidate)
                if len(sequence) != candidate.cds_bp:
                    raise AssertionError("Extracted CDS length differs from GFF")
                total_bases += len(sequence)
                handle.write(
                    f">{candidate.digest} {candidate.seqid}:{candidate.start}-"
                    f"{candidate.end}({candidate.strand})\n"
                )
                for start in range(0, len(sequence), 60):
                    handle.write(sequence[start : start + 60] + "\n")
    finally:
        index.close()
    manifest: dict[str, Any] = {
        "schema_version": NATURAL_CDS_EXPORT_SCHEMA,
        "generator": {"name": "PloidyPatch", "version": __version__},
        "inputs": {
            "candidate_gff": _file_sha256(candidate_gff),
            "genome_fasta": _file_sha256(genome),
            "genome_fai": _file_sha256(fai),
        },
        "counts": {"candidates": len(candidates), "cds_bases": total_bases},
        "output": {
            "file_name": output.name,
            "sha256": _file_sha256(output),
        },
        "truth_access": False,
        "automatic_approval": False,
    }
    with manifest_path.open("x", encoding="utf-8", newline="") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return manifest


def _read_base_features(
    path: Path,
) -> tuple[FeatureIntervalIndex, FeatureIntervalIndex, int, int]:
    genes: list[FeatureInterval] = []
    cds: list[FeatureInterval] = []
    with open_text(path) as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            if not raw_line or raw_line.startswith("#"):
                continue
            fields = raw_line.rstrip("\r\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"Malformed base GFF line {line_number}")
            if fields[2] not in {"gene", "CDS"}:
                continue
            try:
                start, end = int(fields[3]), int(fields[4])
            except ValueError as exc:
                raise ValueError(f"Invalid base GFF coordinate at line {line_number}") from exc
            attributes, malformed = parse_attributes(fields[8])
            if malformed:
                raise ValueError(f"Malformed base GFF attributes at line {line_number}")
            if fields[2] == "gene":
                identifier = attributes.get("ID") or f"gene_line_{line_number}"
                genes.append(
                    FeatureInterval(fields[0], start, end, fields[6], identifier)
                )
            else:
                parents = attributes.get("Parent", f"cds_line_{line_number}")
                cds.append(FeatureInterval(fields[0], start, end, fields[6], parents))
    if not genes or not cds:
        raise ValueError("Base GFF must contain gene and CDS features")
    return FeatureIntervalIndex(genes), FeatureIntervalIndex(cds), len(genes), len(cds)


def _read_seqid_map(path: Path | None) -> dict[str, str]:
    if path is None:
        return {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or not {
            "source_seqid",
            "target_seqid",
        } <= set(reader.fieldnames):
            raise ValueError("Repeat seqid map lacks source_seqid and target_seqid")
        rows = list(reader)
    mapping = {row["source_seqid"]: row["target_seqid"] for row in rows}
    if (
        len(mapping) != len(rows)
        or "" in mapping
        or any(not value for value in mapping.values())
        or len(set(mapping.values())) != len(mapping)
    ):
        raise ValueError("Repeat seqid map must be a nonempty one-to-one mapping")
    return mapping


def _read_repeat_features(
    path: Path | None, seqid_map_path: Path | None
) -> tuple[FeatureIntervalIndex | None, int, set[str], int, set[str]]:
    if path is None:
        if seqid_map_path is not None:
            raise ValueError("Repeat seqid map was supplied without a repeat GFF")
        return None, 0, set(), 0, set()
    mapping = _read_seqid_map(seqid_map_path)
    intervals: list[FeatureInterval] = []
    observed_source_seqids: set[str] = set()
    excluded_source_seqids: set[str] = set()
    excluded_feature_count = 0
    with open_text(path) as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            if not raw_line or raw_line.startswith("#"):
                continue
            fields = raw_line.rstrip("\r\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"Malformed repeat GFF line {line_number}")
            source_seqid = fields[0]
            if mapping and source_seqid not in mapping:
                excluded_source_seqids.add(source_seqid)
                excluded_feature_count += 1
                continue
            observed_source_seqids.add(source_seqid)
            seqid = mapping.get(source_seqid, source_seqid)
            try:
                start, end = int(fields[3]), int(fields[4])
            except ValueError as exc:
                raise ValueError(f"Invalid repeat coordinate at line {line_number}") from exc
            intervals.append(
                FeatureInterval(
                    seqid,
                    start,
                    end,
                    fields[6],
                    f"repeat_line_{line_number}",
                )
            )
    if not intervals:
        raise ValueError("Repeat GFF has no feature intervals")
    return (
        FeatureIntervalIndex(intervals),
        len(intervals),
        observed_source_seqids,
        excluded_feature_count,
        excluded_source_seqids,
    )


def _read_rankings(
    path: Path, candidate_digests: set[str]
) -> tuple[list[dict[str, str]], tuple[str, ...]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or not {
            "estimator",
            "review_rank",
            "candidate_digest",
        } <= set(reader.fieldnames):
            raise ValueError("Review rankings lack required fields")
        rows = list(reader)
    estimators = tuple(sorted({row["estimator"] for row in rows}))
    if not estimators:
        raise ValueError("Review rankings are empty")
    for estimator in estimators:
        selected = [row for row in rows if row["estimator"] == estimator]
        digests = [row["candidate_digest"] for row in selected]
        ranks = [int(row["review_rank"]) for row in selected]
        if set(digests) != candidate_digests or len(digests) != len(set(digests)):
            raise ValueError(f"Ranking universe differs for estimator {estimator}")
        if sorted(ranks) != list(range(1, len(ranks) + 1)):
            raise ValueError(f"Review ranks are not a permutation for {estimator}")
    return rows, estimators


def _read_isoseq(path: Path | None, candidate_digests: set[str]) -> dict[str, dict[str, str]]:
    if path is None:
        return {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or not {
            "candidate_digest",
            "evidence_state",
        } <= set(reader.fieldnames):
            raise ValueError("Iso-Seq evidence lacks required fields")
        rows = list(reader)
    indexed = {row["candidate_digest"]: row for row in rows}
    if len(indexed) != len(rows) or set(indexed) != candidate_digests:
        raise ValueError("Iso-Seq evidence universe differs from candidates")
    return indexed


def _read_self_alignments(path: Path | None) -> dict[str, list[SelfAlignment]]:
    if path is None:
        return {}
    grouped: dict[str, list[SelfAlignment]] = defaultdict(list)
    with path.open("r", encoding="utf-8", newline="") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            fields = raw_line.rstrip("\r\n").split("\t")
            if len(fields) < 12:
                raise ValueError(f"Malformed self-map PAF line {line_number}")
            tags = {
                field[:2]: field[5:]
                for field in fields[12:]
                if len(field) >= 6 and field[2:5] in {":i:", ":A:", ":Z:", ":f:"}
            }
            score = int(tags.get("AS", fields[9]))
            grouped[fields[0]].append(
                SelfAlignment(
                    query=fields[0],
                    query_length=int(fields[1]),
                    query_start=int(fields[2]),
                    query_end=int(fields[3]),
                    strand=fields[4],
                    seqid=fields[5],
                    start=int(fields[7]) + 1,
                    end=int(fields[8]),
                    matches=int(fields[9]),
                    block_length=int(fields[10]),
                    mapq=int(fields[11]),
                    alignment_score=score,
                )
            )
    return grouped


def _sequence_metrics(sequence: str) -> dict[str, Any]:
    codons = [sequence[index : index + 3] for index in range(0, len(sequence) - 2, 3)]
    terminal_stop = bool(codons and codons[-1] in STOP_CODONS)
    internal = codons[:-1] if terminal_stop else codons
    internal_stops = sum(codon in STOP_CODONS for codon in internal)
    ambiguous = sum(base not in "ACGT" for base in sequence)
    valid = len(sequence) - ambiguous
    counts = Counter(base for base in sequence if base in "ACGT")
    entropy = 0.0
    if valid:
        entropy = -sum(
            (count / valid) * math.log2(count / valid) for count in counts.values()
        ) / 2.0
    max_homopolymer = 0
    current = 0
    previous = ""
    for base in sequence:
        if base == previous:
            current += 1
        else:
            previous = base
            current = 1
        max_homopolymer = max(max_homopolymer, current)
    complete_orf = (
        len(sequence) % 3 == 0
        and sequence.startswith("ATG")
        and terminal_stop
        and internal_stops == 0
        and ambiguous == 0
    )
    return {
        "cds_sequence_sha256": hashlib.sha256(sequence.encode("ascii")).hexdigest(),
        "cds_length": len(sequence),
        "cds_modulo_three": len(sequence) % 3,
        "start_codon": sequence[:3],
        "terminal_codon": sequence[-3:] if len(sequence) >= 3 else sequence,
        "has_atg_start": int(sequence.startswith("ATG")),
        "has_terminal_stop": int(terminal_stop),
        "internal_stop_codons": internal_stops,
        "ambiguous_bases": ambiguous,
        "gc_fraction": (
            (counts.get("G", 0) + counts.get("C", 0)) / valid if valid else 0.0
        ),
        "mononucleotide_entropy": entropy,
        "maximum_homopolymer": max_homopolymer,
        "complete_orf": int(complete_orf),
    }


def _merged_overlap_bp(
    candidate_intervals: Iterable[tuple[int, int]],
    features: Iterable[FeatureInterval],
) -> int:
    feature_intervals = _merge_intervals(
        (feature.start, feature.end) for feature in features
    )
    return _overlap_bp(candidate_intervals, feature_intervals)


def _is_expected_alignment(candidate: Candidate, alignment: SelfAlignment) -> bool:
    if alignment.seqid != candidate.seqid or alignment.strand != candidate.strand:
        return False
    overlap = max(0, min(candidate.end, alignment.end) - max(candidate.start, alignment.start) + 1)
    candidate_span = candidate.end - candidate.start + 1
    alignment_span = alignment.end - alignment.start + 1
    return overlap / candidate_span >= 0.8 and overlap / alignment_span >= 0.8


def _distinct_locus_count(alignments: Iterable[SelfAlignment]) -> int:
    grouped: dict[tuple[str, str], list[tuple[int, int]]] = defaultdict(list)
    for alignment in alignments:
        grouped[(alignment.seqid, alignment.strand)].append(
            (alignment.start, alignment.end)
        )
    count = 0
    for intervals in grouped.values():
        merged: list[list[int]] = []
        for start, end in sorted(intervals):
            if merged and start <= merged[-1][1] + 1000:
                merged[-1][1] = max(merged[-1][1], end)
            else:
                merged.append([start, end])
        count += len(merged)
    return count


def _mappability_metrics(
    candidate: Candidate,
    alignments: list[SelfAlignment] | None,
    *,
    minimum_query_coverage: float,
    minimum_identity: float,
    near_equal_score_fraction: float,
) -> dict[str, Any]:
    if alignments is None:
        return {
            "self_alignments": "",
            "high_quality_loci": "",
            "expected_locus_found": "",
            "expected_locus_query_coverage": "",
            "expected_locus_identity": "",
            "expected_locus_mapq": "",
            "near_equal_nonexpected_loci": "",
            "best_nonexpected_score_fraction": "",
            "mappability_class": "not_evaluated",
        }
    high_quality = [
        alignment
        for alignment in alignments
        if alignment.query_coverage >= minimum_query_coverage
        and alignment.identity >= minimum_identity
    ]
    expected = [
        alignment
        for alignment in high_quality
        if _is_expected_alignment(candidate, alignment)
    ]
    if not expected:
        return {
            "self_alignments": len(alignments),
            "high_quality_loci": _distinct_locus_count(high_quality),
            "expected_locus_found": 0,
            "expected_locus_query_coverage": "",
            "expected_locus_identity": "",
            "expected_locus_mapq": "",
            "near_equal_nonexpected_loci": "",
            "best_nonexpected_score_fraction": "",
            "mappability_class": "expected_locus_missing",
        }
    best_expected = max(expected, key=lambda item: item.alignment_score)
    alternatives = [
        alignment
        for alignment in high_quality
        if not _is_expected_alignment(candidate, alignment)
    ]
    near_equal = [
        alignment
        for alignment in alternatives
        if alignment.alignment_score
        >= near_equal_score_fraction * best_expected.alignment_score
    ]
    near_equal_loci = _distinct_locus_count(near_equal)
    best_alternative = max(
        (alignment.alignment_score for alignment in alternatives), default=0
    )
    if near_equal_loci >= 3:
        category = "high_copy_or_repeat"
    elif near_equal_loci:
        category = "duplicated_locus"
    else:
        category = "unique_locus"
    return {
        "self_alignments": len(alignments),
        "high_quality_loci": _distinct_locus_count(high_quality),
        "expected_locus_found": 1,
        "expected_locus_query_coverage": best_expected.query_coverage,
        "expected_locus_identity": best_expected.identity,
        "expected_locus_mapq": best_expected.mapq,
        "near_equal_nonexpected_loci": near_equal_loci,
        "best_nonexpected_score_fraction": (
            best_alternative / best_expected.alignment_score
            if best_expected.alignment_score
            else 0.0
        ),
        "mappability_class": category,
    }


def _count_rows(rows: Iterable[dict[str, Any]], field: str) -> dict[str, int]:
    return dict(sorted(Counter(str(row[field]) for row in rows).items()))


def audit_natural_candidates(
    *,
    candidate_gff_path: str | Path,
    base_gff_path: str | Path,
    genome_fasta_path: str | Path,
    review_rankings_tsv_path: str | Path,
    output_tsv_path: str | Path,
    output_summary_json_path: str | Path,
    isoseq_evidence_tsv_path: str | Path | None = None,
    self_map_paf_path: str | Path | None = None,
    repeat_gff_path: str | Path | None = None,
    repeat_seqid_map_tsv_path: str | Path | None = None,
    repeat_flank_bp: int = 2000,
    minimum_full_length_read_support: int = 1,
    review_budgets: Iterable[int] = (25, 50, 100, 200),
    minimum_query_coverage: float = 0.90,
    minimum_identity: float = 0.98,
    near_equal_score_fraction: float = 0.95,
) -> dict[str, Any]:
    if not 0 < minimum_query_coverage <= 1 or not 0 < minimum_identity <= 1:
        raise ValueError("Mappability thresholds must be in (0,1]")
    if not 0 < near_equal_score_fraction <= 1:
        raise ValueError("Near-equal score fraction must be in (0,1]")
    budgets = tuple(sorted(set(int(value) for value in review_budgets)))
    if not budgets or budgets[0] < 1:
        raise ValueError("Review budgets must be positive")
    if repeat_flank_bp < 0:
        raise ValueError("Repeat flank length must be non-negative")
    if minimum_full_length_read_support < 1:
        raise ValueError("Minimum full-length read support must be positive")

    candidate_gff = Path(candidate_gff_path)
    base_gff = Path(base_gff_path)
    genome = Path(genome_fasta_path)
    fai = Path(str(genome) + ".fai")
    rankings = Path(review_rankings_tsv_path)
    isoseq_path = Path(isoseq_evidence_tsv_path) if isoseq_evidence_tsv_path else None
    self_map_path = Path(self_map_paf_path) if self_map_paf_path else None
    repeat_path = Path(repeat_gff_path) if repeat_gff_path else None
    repeat_map_path = (
        Path(repeat_seqid_map_tsv_path) if repeat_seqid_map_tsv_path else None
    )
    output = Path(output_tsv_path)
    summary_path = Path(output_summary_json_path)
    manifest_path = Path(str(output) + ".manifest.json")
    required = [candidate_gff, base_gff, genome, fai, rankings]
    required.extend(
        path
        for path in (isoseq_path, self_map_path, repeat_path, repeat_map_path)
        if path is not None
    )
    for path in required:
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f"Missing or empty natural-audit input: {path}")
    _refuse_outputs((output, summary_path, manifest_path))

    candidates = sorted(_read_candidates(candidate_gff), key=lambda item: item.digest)
    candidate_by_digest = {candidate.digest: candidate for candidate in candidates}
    if len(candidate_by_digest) != len(candidates) or not candidates:
        raise ValueError("Natural candidate universe is empty or duplicated")
    candidate_digests = set(candidate_by_digest)
    ranking_rows, estimators = _read_rankings(rankings, candidate_digests)
    isoseq = _read_isoseq(isoseq_path, candidate_digests)
    self_alignments = _read_self_alignments(self_map_path)
    unknown_queries = set(self_alignments) - candidate_digests
    if unknown_queries:
        raise ValueError(f"Self-map PAF contains unknown candidate queries: {len(unknown_queries)}")
    gene_index, cds_index, base_gene_count, base_cds_count = _read_base_features(base_gff)
    (
        repeat_index,
        repeat_feature_count,
        repeat_source_seqids,
        repeat_excluded_feature_count,
        repeat_excluded_source_seqids,
    ) = _read_repeat_features(repeat_path, repeat_map_path)

    candidate_metrics: dict[str, dict[str, Any]] = {}
    fasta_index = FastaIndex(genome, fai)
    try:
        for candidate in candidates:
            sequence = _candidate_cds(fasta_index, candidate)
            if len(sequence) != candidate.cds_bp:
                raise AssertionError("Candidate CDS extraction length mismatch")
            genes = gene_index.overlaps(candidate.seqid, candidate.start, candidate.end)
            cds = cds_index.overlaps(candidate.seqid, candidate.start, candidate.end)
            same_genes = [item for item in genes if item.strand == candidate.strand]
            opposite_genes = [item for item in genes if item.strand != candidate.strand]
            same_cds = [item for item in cds if item.strand == candidate.strand]
            opposite_cds = [item for item in cds if item.strand != candidate.strand]
            gene_overlap = _merged_overlap_bp(
                ((candidate.start, candidate.end),), genes
            )
            cds_overlap = _merged_overlap_bp(candidate.cds, cds)
            if same_genes and opposite_genes:
                collision = "mixed_strand_gene_overlap"
            elif same_genes:
                collision = "same_strand_gene_overlap"
            elif opposite_genes:
                collision = "opposite_strand_gene_overlap"
            else:
                collision = "intergenic"
            evidence = isoseq.get(candidate.digest, {})
            evidence_state = evidence.get("evidence_state", "not_evaluated")
            supporting_read_text = evidence.get(
                "supporting_full_length_reads",
                evidence.get("supporting_b73_full_length_reads", "0"),
            )
            try:
                supporting_full_length_reads = int(supporting_read_text or 0)
            except ValueError as exc:
                raise ValueError(
                    "Natural evidence has a non-integer full-length read count"
                ) from exc
            if supporting_full_length_reads < 0:
                raise ValueError("Natural evidence has a negative read count")
            if repeat_index is None:
                repeat_span: list[FeatureInterval] = []
                repeat_cds: list[FeatureInterval] = []
                repeat_flank: list[FeatureInterval] = []
                repeat_span_bp: int | str = ""
                repeat_cds_bp: int | str = ""
                repeat_flank_overlap_bp: int | str = ""
                repeat_context = "not_evaluated"
            else:
                repeat_span = repeat_index.overlaps(
                    candidate.seqid, candidate.start, candidate.end
                )
                repeat_cds = [
                    interval
                    for start, end in candidate.cds
                    for interval in repeat_index.overlaps(candidate.seqid, start, end)
                ]
                flank_start = max(1, candidate.start - repeat_flank_bp)
                flank_end = candidate.end + repeat_flank_bp
                repeat_flank = repeat_index.overlaps(
                    candidate.seqid, flank_start, flank_end
                )
                repeat_span_bp = _merged_overlap_bp(
                    ((candidate.start, candidate.end),), repeat_span
                )
                repeat_cds_bp = _merged_overlap_bp(candidate.cds, repeat_cds)
                repeat_flank_overlap_bp = _merged_overlap_bp(
                    ((flank_start, flank_end),), repeat_flank
                )
                if repeat_cds_bp:
                    repeat_context = "cds_repeat_overlap"
                elif repeat_span_bp:
                    repeat_context = "intronic_or_span_repeat_overlap"
                elif repeat_flank:
                    repeat_context = "flanking_repeat_only"
                else:
                    repeat_context = "repeat_distant"
            metrics = {
                "candidate_digest": candidate.digest,
                "seqid": candidate.seqid,
                "start": candidate.start,
                "end": candidate.end,
                "strand": candidate.strand,
                "span_bp": candidate.end - candidate.start + 1,
                "cds_segments": len(candidate.cds),
                **_sequence_metrics(sequence),
                "overlapping_base_genes": len({item.identifier for item in genes}),
                "same_strand_base_genes": len(
                    {item.identifier for item in same_genes}
                ),
                "opposite_strand_base_genes": len(
                    {item.identifier for item in opposite_genes}
                ),
                "base_gene_overlap_bp": gene_overlap,
                "base_gene_overlap_fraction": gene_overlap
                / (candidate.end - candidate.start + 1),
                "base_cds_overlap_bp": cds_overlap,
                "base_cds_overlap_fraction": cds_overlap / candidate.cds_bp,
                "same_strand_base_cds_overlap_bp": _merged_overlap_bp(
                    candidate.cds, same_cds
                ),
                "opposite_strand_base_cds_overlap_bp": _merged_overlap_bp(
                    candidate.cds, opposite_cds
                ),
                "nearest_base_gene_distance_bp": gene_index.nearest_distance(
                    candidate.seqid, candidate.start, candidate.end
                ),
                "collision_class": collision,
                "repeat_features_overlapping_span": (
                    len(repeat_span) if repeat_index is not None else ""
                ),
                "repeat_span_overlap_bp": repeat_span_bp,
                "repeat_span_overlap_fraction": (
                    repeat_span_bp / (candidate.end - candidate.start + 1)
                    if isinstance(repeat_span_bp, int)
                    else ""
                ),
                "repeat_cds_overlap_bp": repeat_cds_bp,
                "repeat_cds_overlap_fraction": (
                    repeat_cds_bp / candidate.cds_bp
                    if isinstance(repeat_cds_bp, int)
                    else ""
                ),
                "repeat_flank_overlap_bp": repeat_flank_overlap_bp,
                "repeat_context_class": repeat_context,
                **_mappability_metrics(
                    candidate,
                    self_alignments.get(candidate.digest)
                    if self_map_path is not None
                    else None,
                    minimum_query_coverage=minimum_query_coverage,
                    minimum_identity=minimum_identity,
                    near_equal_score_fraction=near_equal_score_fraction,
                ),
                "isoseq_evidence_state": evidence_state,
                "isoseq_supporting_transcripts": evidence.get(
                    "supporting_transcripts", ""
                ),
                "isoseq_supporting_full_length_reads": evidence.get(
                    "supporting_full_length_reads",
                    evidence.get("supporting_b73_full_length_reads", ""),
                ),
                "isoseq_supporting_query_groups": evidence.get(
                    "supporting_query_groups", ""
                ),
            }
            metrics["strict_collision_free"] = int(
                metrics["base_gene_overlap_bp"] == 0
                and metrics["base_cds_overlap_bp"] == 0
            )
            metrics["transcript_chain_supported"] = int(
                evidence_state == "full_chain_supported" and len(candidate.cds) >= 2
                and supporting_full_length_reads
                >= minimum_full_length_read_support
            )
            metrics["single_exon_span_supported"] = int(
                evidence_state == "single_exon_span_support"
            )
            metrics["repeat_risk_flag"] = int(
                metrics["mappability_class"] == "high_copy_or_repeat"
                or (
                    isinstance(metrics["repeat_cds_overlap_bp"], int)
                    and metrics["repeat_cds_overlap_bp"] > 0
                )
                or metrics["mononucleotide_entropy"] < 0.75
                or metrics["maximum_homopolymer"] >= 20
            )
            metrics["case_study_ready"] = int(
                metrics["complete_orf"] == 1
                and metrics["strict_collision_free"] == 1
                and metrics["mappability_class"]
                in {"unique_locus", "duplicated_locus"}
                and metrics["repeat_risk_flag"] == 0
                and metrics["transcript_chain_supported"] == 1
            )
            candidate_metrics[candidate.digest] = metrics
    finally:
        fasta_index.close()

    output_rows: list[dict[str, Any]] = []
    for ranking in sorted(
        ranking_rows,
        key=lambda row: (row["estimator"], int(row["review_rank"])),
    ):
        digest = ranking["candidate_digest"]
        output_rows.append(
            {
                "estimator": ranking["estimator"],
                "review_rank": int(ranking["review_rank"]),
                **candidate_metrics[digest],
            }
        )
    output_fields = list(output_rows[0])
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=output_fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for row in output_rows:
            writer.writerow(
                {
                    key: format(value, ".17g") if isinstance(value, float) else value
                    for key, value in row.items()
                }
            )

    unique_rows = list(candidate_metrics.values())
    estimator_budgets: list[dict[str, Any]] = []
    for estimator in estimators:
        estimator_rows = [row for row in output_rows if row["estimator"] == estimator]
        for budget in budgets:
            selected = [row for row in estimator_rows if row["review_rank"] <= budget]
            estimator_budgets.append(
                {
                    "estimator": estimator,
                    "review_budget": budget,
                    "candidates": len(selected),
                    "complete_orf": sum(row["complete_orf"] for row in selected),
                    "strict_collision_free": sum(
                        row["strict_collision_free"] for row in selected
                    ),
                    "repeat_risk_flag": sum(
                        row["repeat_risk_flag"] for row in selected
                    ),
                    "transcript_chain_supported": sum(
                        row["transcript_chain_supported"] for row in selected
                    ),
                    "single_exon_span_supported": sum(
                        row["single_exon_span_supported"] for row in selected
                    ),
                    "case_study_ready": sum(
                        row["case_study_ready"] for row in selected
                    ),
                }
            )
    summary: dict[str, Any] = {
        "schema_version": NATURAL_AUDIT_SCHEMA,
        "generator": {"name": "PloidyPatch", "version": __version__},
        "interpretation": {
            "review_only": True,
            "automatic_approval": False,
            "case_study_ready_is_not_automatic_acceptance": True,
            "repeat_risk_definition": (
                "official_repeat_CDS_overlap_or_at_least_three_near_equal_"
                "nonexpected_loci_or_"
                "mononucleotide_entropy_below_0.75_or_homopolymer_at_least_20"
            ),
            "strict_collision_free_definition": (
                "zero_overlap_with_any_existing_gene_span_and_CDS"
            ),
        },
        "parameters": {
            "review_budgets": list(budgets),
            "minimum_self_map_query_coverage": minimum_query_coverage,
            "minimum_self_map_identity": minimum_identity,
            "near_equal_score_fraction": near_equal_score_fraction,
            "near_equal_loci_merge_distance_bp": 1000,
            "high_copy_nonexpected_loci_threshold": 3,
            "repeat_flank_bp": repeat_flank_bp,
            "minimum_full_length_read_support": minimum_full_length_read_support,
            "repeat_seqid_map_policy": (
                "exact_allowlist_with_unmapped_source_seqids_counted_and_excluded"
                if repeat_map_path is not None
                else "identity"
            ),
        },
        "counts": {
            "candidates": len(unique_rows),
            "ranking_rows": len(output_rows),
            "estimators": list(estimators),
            "base_genes": base_gene_count,
            "base_cds_features": base_cds_count,
            "repeat_features": repeat_feature_count,
            "repeat_source_seqids": sorted(repeat_source_seqids),
            "repeat_features_excluded_by_seqid_map": (
                repeat_excluded_feature_count
            ),
            "repeat_source_seqids_excluded_by_seqid_map": sorted(
                repeat_excluded_source_seqids
            ),
            "complete_orf": sum(row["complete_orf"] for row in unique_rows),
            "strict_collision_free": sum(
                row["strict_collision_free"] for row in unique_rows
            ),
            "repeat_risk_flag": sum(
                row["repeat_risk_flag"] for row in unique_rows
            ),
            "transcript_chain_supported": sum(
                row["transcript_chain_supported"] for row in unique_rows
            ),
            "single_exon_span_supported": sum(
                row["single_exon_span_supported"] for row in unique_rows
            ),
            "case_study_ready": sum(
                row["case_study_ready"] for row in unique_rows
            ),
            "orf_state": _count_rows(unique_rows, "complete_orf"),
            "collision_class": _count_rows(unique_rows, "collision_class"),
            "mappability_class": _count_rows(unique_rows, "mappability_class"),
            "isoseq_evidence_state": _count_rows(
                unique_rows, "isoseq_evidence_state"
            ),
            "repeat_context_class": _count_rows(
                unique_rows, "repeat_context_class"
            ),
        },
        "review_yield": estimator_budgets,
        "inputs": {
            "candidate_gff": _file_sha256(candidate_gff),
            "base_gff": _file_sha256(base_gff),
            "genome_fasta": _file_sha256(genome),
            "genome_fai": _file_sha256(fai),
            "review_rankings": _file_sha256(rankings),
            "isoseq_evidence": _file_sha256(isoseq_path) if isoseq_path else None,
            "self_map_paf": _file_sha256(self_map_path) if self_map_path else None,
            "repeat_gff": _file_sha256(repeat_path) if repeat_path else None,
            "repeat_seqid_map": (
                _file_sha256(repeat_map_path) if repeat_map_path else None
            ),
        },
        "output": {
            "file_name": output.name,
            "rows": len(output_rows),
            "sha256": _file_sha256(output),
        },
    }
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    with summary_path.open("x", encoding="utf-8", newline="") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")
    manifest = {
        "schema_version": NATURAL_AUDIT_SCHEMA,
        "truth_access": isoseq_path is not None,
        "automatic_approval": False,
        "inputs": summary["inputs"],
        "outputs": {
            "audit_tsv": summary["output"],
            "summary_json": {
                "file_name": summary_path.name,
                "sha256": _file_sha256(summary_path),
            },
        },
    }
    with manifest_path.open("x", encoding="utf-8", newline="") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return summary
