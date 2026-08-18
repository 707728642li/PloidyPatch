from __future__ import annotations

import hashlib
import json
import csv
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from . import __version__
from .perturb import (
    SUPPORTED_ANNOTATION_EVENTS,
    TRUTH_SCHEMA_VERSION,
    GffDocument,
    read_gff_document,
)


SCORE_SCHEMA_VERSION = "ploidypatch.annotation_repair_score.v5"
TRANSCRIPT_TYPES = {"mRNA", "transcript"}


@dataclass(frozen=True)
class TranscriptSignature:
    """Identifier-independent transcript structure used for strict scoring."""

    seqid: str
    strand: str
    start: int
    end: int
    exons: tuple[tuple[int, int], ...]
    cds: tuple[tuple[int, int, str], ...]

    @property
    def digest(self) -> str:
        payload = json.dumps(
            {
                "seqid": self.seqid,
                "strand": self.strand,
                "start": self.start,
                "end": self.end,
                "exons": self.exons,
                "cds": self.cds,
            },
            sort_keys=True,
            separators=(",", ":"),
        )
        return hashlib.sha256(payload.encode("utf-8")).hexdigest()


@dataclass(frozen=True)
class TranscriptModel:
    transcript_id: str
    gene_ids: tuple[str, ...]
    signature: TranscriptSignature


@dataclass(frozen=True)
class CdsChainSignature:
    """Identifier- and UTR-independent phased CDS chain."""

    seqid: str
    strand: str
    cds: tuple[tuple[int, int, str], ...]


@dataclass(frozen=True)
class AnnotationIndex:
    transcripts: dict[str, TranscriptModel]
    signatures: frozenset[TranscriptSignature]
    gene_signatures: dict[str, frozenset[TranscriptSignature]]


def build_annotation_index(document: GffDocument) -> AnnotationIndex:
    """Build identifier and structural indexes from a parsed GFF3 document."""

    transcript_records = {}
    exon_by_parent: dict[str, list[tuple[int, int]]] = defaultdict(list)
    cds_by_parent: dict[str, list[tuple[int, int, str]]] = defaultdict(list)
    for record in document.records:
        if record.feature_type in TRANSCRIPT_TYPES and record.feature_id:
            if record.feature_id in transcript_records:
                raise ValueError(f"Duplicate transcript ID: {record.feature_id}")
            transcript_records[record.feature_id] = record
        elif record.feature_type == "exon":
            for parent in record.parents:
                exon_by_parent[parent].append((record.start, record.end))
        elif record.feature_type == "CDS":
            for parent in record.parents:
                cds_by_parent[parent].append((record.start, record.end, record.phase))

    transcripts: dict[str, TranscriptModel] = {}
    gene_signatures_mutable: dict[str, set[TranscriptSignature]] = defaultdict(set)
    for transcript_id, record in transcript_records.items():
        signature = TranscriptSignature(
            seqid=record.seqid,
            strand=record.strand,
            start=record.start,
            end=record.end,
            exons=tuple(sorted(set(exon_by_parent.get(transcript_id, [])))),
            cds=tuple(sorted(set(cds_by_parent.get(transcript_id, [])))),
        )
        model = TranscriptModel(
            transcript_id=transcript_id,
            gene_ids=tuple(sorted(set(record.parents))),
            signature=signature,
        )
        transcripts[transcript_id] = model
        for gene_id in model.gene_ids:
            gene_signatures_mutable[gene_id].add(signature)

    gene_signatures = {
        gene_id: frozenset(signatures)
        for gene_id, signatures in gene_signatures_mutable.items()
    }
    return AnnotationIndex(
        transcripts=transcripts,
        signatures=frozenset(model.signature for model in transcripts.values()),
        gene_signatures=gene_signatures,
    )


def _file_sha256(path: str | Path, block_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        while block := handle.read(block_size):
            digest.update(block)
    return digest.hexdigest()


def _rate(numerator: int, denominator: int) -> float | None:
    return numerator / denominator if denominator else None


def _f1(precision: float | None, recall: float | None) -> float | None:
    if precision is None or recall is None or precision + recall == 0:
        return None
    return 2 * precision * recall / (precision + recall)


def _set_metrics(
    truth_features: set[tuple[Any, ...]],
    candidate_features: set[tuple[Any, ...]],
) -> dict[str, int | float | None]:
    true_positive = len(truth_features & candidate_features)
    false_positive = len(candidate_features - truth_features)
    false_negative = len(truth_features - candidate_features)
    precision = _rate(true_positive, true_positive + false_positive)
    recall = _rate(true_positive, true_positive + false_negative)
    return {
        "true_positive": true_positive,
        "false_positive": false_positive,
        "false_negative": false_negative,
        "precision": precision,
        "recall": recall,
        "f1": _f1(precision, recall),
        "truth_features": len(truth_features),
        "candidate_novel_features": len(candidate_features),
    }


def _partial_feature_sets(
    signatures: set[TranscriptSignature],
) -> tuple[
    set[tuple[Any, ...]],
    set[tuple[Any, ...]],
    set[tuple[Any, ...]],
]:
    exon_segments: set[tuple[Any, ...]] = set()
    splice_junctions: set[tuple[Any, ...]] = set()
    cds_segments: set[tuple[Any, ...]] = set()
    for signature in signatures:
        for start, end in signature.exons:
            exon_segments.add((signature.seqid, signature.strand, start, end))
        sorted_exons = sorted(signature.exons)
        for left, right in zip(sorted_exons, sorted_exons[1:]):
            splice_junctions.add(
                (signature.seqid, signature.strand, left[1], right[0])
            )
        for start, end, phase in signature.cds:
            cds_segments.add(
                (signature.seqid, signature.strand, start, end, phase)
            )
    return exon_segments, splice_junctions, cds_segments


def _cds_chain_signatures(
    signatures: set[TranscriptSignature] | frozenset[TranscriptSignature],
) -> set[CdsChainSignature]:
    return {
        CdsChainSignature(signature.seqid, signature.strand, signature.cds)
        for signature in signatures
        if signature.cds
    }


def _merged_cds_intervals(
    signatures: set[TranscriptSignature],
) -> dict[tuple[str, str], list[tuple[int, int]]]:
    raw: dict[tuple[str, str], list[tuple[int, int]]] = defaultdict(list)
    for signature in signatures:
        key = (signature.seqid, signature.strand)
        raw[key].extend((start, end) for start, end, _ in signature.cds)
    merged: dict[tuple[str, str], list[tuple[int, int]]] = {}
    for key, intervals in raw.items():
        combined: list[list[int]] = []
        for start, end in sorted(intervals):
            if combined and start <= combined[-1][1] + 1:
                combined[-1][1] = max(combined[-1][1], end)
            else:
                combined.append([start, end])
        merged[key] = [(start, end) for start, end in combined]
    return merged


def _interval_bp(intervals: dict[tuple[str, str], list[tuple[int, int]]]) -> int:
    return sum(
        end - start + 1
        for values in intervals.values()
        for start, end in values
    )


def _intersection_bp(
    left: dict[tuple[str, str], list[tuple[int, int]]],
    right: dict[tuple[str, str], list[tuple[int, int]]],
) -> int:
    total = 0
    for key in set(left) & set(right):
        left_values = left[key]
        right_values = right[key]
        left_index = 0
        right_index = 0
        while left_index < len(left_values) and right_index < len(right_values):
            left_start, left_end = left_values[left_index]
            right_start, right_end = right_values[right_index]
            start = max(left_start, right_start)
            end = min(left_end, right_end)
            if start <= end:
                total += end - start + 1
            if left_end <= right_end:
                left_index += 1
            else:
                right_index += 1
    return total


def _cds_coverage_metrics(
    truth_signatures: set[TranscriptSignature],
    candidate_signatures: set[TranscriptSignature],
    background_signatures: set[TranscriptSignature] | None = None,
) -> dict[str, int | float | None]:
    truth_intervals = _merged_cds_intervals(truth_signatures)
    candidate_intervals = _merged_cds_intervals(candidate_signatures)
    if background_signatures:
        candidate_intervals = _subtract_merged_intervals(
            candidate_intervals,
            _merged_cds_intervals(background_signatures),
        )
    truth_bp = _interval_bp(truth_intervals)
    candidate_bp = _interval_bp(candidate_intervals)
    shared_bp = _intersection_bp(truth_intervals, candidate_intervals)
    precision = _rate(shared_bp, candidate_bp)
    recall = _rate(shared_bp, truth_bp)
    return {
        "intersection_bp": shared_bp,
        "truth_bp": truth_bp,
        "candidate_novel_bp": candidate_bp,
        "precision": precision,
        "recall": recall,
        "f1": _f1(precision, recall),
    }


def _subtract_merged_intervals(
    left: dict[tuple[str, str], list[tuple[int, int]]],
    right: dict[tuple[str, str], list[tuple[int, int]]],
) -> dict[tuple[str, str], list[tuple[int, int]]]:
    result: dict[tuple[str, str], list[tuple[int, int]]] = {}
    for key, left_intervals in left.items():
        right_intervals = right.get(key, [])
        fragments: list[tuple[int, int]] = []
        for left_start, left_end in left_intervals:
            cursor = left_start
            for right_start, right_end in right_intervals:
                if right_end < cursor:
                    continue
                if right_start > left_end:
                    break
                if right_start > cursor:
                    fragments.append((cursor, min(left_end, right_start - 1)))
                cursor = max(cursor, right_end + 1)
                if cursor > left_end:
                    break
            if cursor <= left_end:
                fragments.append((cursor, left_end))
        if fragments:
            result[key] = fragments
    return result


def _read_event_strata(
    path: str | Path,
    columns: tuple[str, ...],
) -> dict[str, dict[str, str]]:
    rows: dict[str, dict[str, str]] = {}
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"gene_id", *columns}
        if reader.fieldnames is None or not required <= set(reader.fieldnames):
            raise ValueError(
                "Event strata TSV is missing required column(s): "
                + ", ".join(sorted(required - set(reader.fieldnames or [])))
            )
        for line_number, row in enumerate(reader, start=2):
            gene_id = row["gene_id"]
            if not gene_id or gene_id in rows:
                raise ValueError(
                    f"Empty or duplicate gene_id in event strata line {line_number}"
                )
            rows[gene_id] = {column: row[column] for column in columns}
    return rows


def _stratum_recall_row(
    labels: dict[str, str],
    evaluations: list[dict[str, Any]],
    novel_candidate: set[TranscriptSignature],
) -> dict[str, Any]:
    truth_signatures: set[TranscriptSignature] = set()
    for evaluation in evaluations:
        truth_signatures.update(evaluation["signatures"])
    recovered_signatures = truth_signatures & novel_candidate
    truth_cds_chains = _cds_chain_signatures(truth_signatures)
    candidate_cds_chains = _cds_chain_signatures(novel_candidate)
    recovered_cds_chains = truth_cds_chains & candidate_cds_chains
    hidden_exons, hidden_junctions, hidden_cds = _partial_feature_sets(
        truth_signatures
    )
    candidate_exons, candidate_junctions, candidate_cds = _partial_feature_sets(
        novel_candidate
    )
    cds_coverage = _cds_coverage_metrics(truth_signatures, novel_candidate)
    return {
        **labels,
        "events": len(evaluations),
        "strict_transcript_structures": {
            "truth": len(truth_signatures),
            "recovered": len(recovered_signatures),
            "recall": _rate(len(recovered_signatures), len(truth_signatures)),
        },
        "strict_cds_chains": {
            "truth": len(truth_cds_chains),
            "recovered": len(recovered_cds_chains),
            "recall": _rate(len(recovered_cds_chains), len(truth_cds_chains)),
        },
        "complete_events": {
            "recovered": sum(
                int(evaluation["complete"]) for evaluation in evaluations
            ),
            "recall": _rate(
                sum(int(evaluation["complete"]) for evaluation in evaluations),
                len(evaluations),
            ),
        },
        "exact_gene_groups": {
            "recovered": sum(
                int(evaluation["exact_gene"]) for evaluation in evaluations
            ),
            "recall": _rate(
                sum(int(evaluation["exact_gene"]) for evaluation in evaluations),
                len(evaluations),
            ),
        },
        "partial_recall": {
            "exon_segments": {
                "truth": len(hidden_exons),
                "recovered": len(hidden_exons & candidate_exons),
                "recall": _rate(
                    len(hidden_exons & candidate_exons), len(hidden_exons)
                ),
            },
            "splice_junctions": {
                "truth": len(hidden_junctions),
                "recovered": len(hidden_junctions & candidate_junctions),
                "recall": _rate(
                    len(hidden_junctions & candidate_junctions),
                    len(hidden_junctions),
                ),
            },
            "cds_segments_with_phase": {
                "truth": len(hidden_cds),
                "recovered": len(hidden_cds & candidate_cds),
                "recall": _rate(len(hidden_cds & candidate_cds), len(hidden_cds)),
            },
            "cds_nucleotide_coverage": {
                "truth_bp": cds_coverage["truth_bp"],
                "recovered_bp": cds_coverage["intersection_bp"],
                "recall": cds_coverage["recall"],
            },
        },
    }


def _stratified_recall(
    *,
    evaluations: list[dict[str, Any]],
    novel_candidate: set[TranscriptSignature],
    strata_path: str | Path,
    columns: tuple[str, ...],
) -> dict[str, Any]:
    strata = _read_event_strata(strata_path, columns)
    evaluation_rows: list[tuple[dict[str, Any], dict[str, str]]] = []
    for evaluation in evaluations:
        gene_id = evaluation["event"]["target"]["gene_id"]
        values = strata.get(gene_id)
        if values is None:
            raise ValueError(f"Hidden event is absent from event strata TSV: {gene_id}")
        evaluation_rows.append((evaluation, values))

    marginal: dict[str, list[dict[str, Any]]] = {}
    for column in columns:
        grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
        for evaluation, values in evaluation_rows:
            grouped[values[column]].append(evaluation)
        marginal[column] = [
            _stratum_recall_row({column: value}, grouped[value], novel_candidate)
            for value in sorted(grouped)
        ]

    joint_groups: dict[tuple[str, ...], list[dict[str, Any]]] = defaultdict(list)
    for evaluation, values in evaluation_rows:
        key = tuple(values[column] for column in columns)
        joint_groups[key].append(evaluation)
    joint = [
        _stratum_recall_row(
            {
                column: value
                for column, value in zip(columns, key, strict=True)
            },
            joint_groups[key],
            novel_candidate,
        )
        for key in sorted(joint_groups)
    ]
    return {
        "source": {
            "file_name": Path(strata_path).name,
            "sha256": _file_sha256(strata_path),
        },
        "columns": list(columns),
        "marginal": marginal,
        "joint": joint,
        "note": "Recall is stratified by hidden events; precision remains global.",
    }


def score_annotation_repair(
    source_gff_path: str | Path,
    perturbed_gff_path: str | Path,
    candidate_gff_path: str | Path,
    truth_path: str | Path,
    include_event_details: bool = False,
    event_strata_path: str | Path | None = None,
    stratum_columns: tuple[str, ...] = (),
    control_candidate_gff_path: str | Path | None = None,
) -> dict[str, Any]:
    """Score strict identifier-independent recovery of hidden gene structures.

    The source annotation and hidden truth are evaluator-only. Candidate
    methods receive only the perturbed inputs and declared public evidence.
    """

    truth = json.loads(Path(truth_path).read_text(encoding="utf-8"))
    if truth.get("schema_version") != TRUTH_SCHEMA_VERSION:
        raise ValueError("Unsupported or missing hidden-truth schema version")
    unsupported = {
        event.get("event_type")
        for event in truth.get("events", [])
        if event.get("event_type") not in SUPPORTED_ANNOTATION_EVENTS
    }
    if unsupported:
        raise ValueError(f"Unsupported event type(s) for this scorer: {unsupported}")

    source_document = read_gff_document(source_gff_path)
    if source_document.text_sha256 != truth["source"]["text_sha256"]:
        raise ValueError("Source GFF3 text checksum does not match hidden truth")
    source = build_annotation_index(source_document)
    perturbed_document = read_gff_document(perturbed_gff_path)
    expected_perturbed_sha = truth.get("perturbation", {}).get(
        "perturbed_text_sha256"
    )
    if (
        expected_perturbed_sha is not None
        and perturbed_document.text_sha256 != expected_perturbed_sha
    ):
        raise ValueError("Perturbed GFF3 text checksum does not match hidden truth")
    perturbed = build_annotation_index(perturbed_document)
    source_resolved = Path(source_gff_path).resolve()
    perturbed_resolved = Path(perturbed_gff_path).resolve()
    candidate_resolved = Path(candidate_gff_path).resolve()
    if candidate_resolved == source_resolved:
        candidate = source
    elif candidate_resolved == perturbed_resolved:
        candidate = perturbed
    else:
        candidate = build_annotation_index(read_gff_document(candidate_gff_path))

    hidden_signatures: set[TranscriptSignature] = set()
    event_truth: list[
        tuple[
            dict[str, Any],
            frozenset[TranscriptSignature],
            frozenset[frozenset[TranscriptSignature]],
            frozenset[TranscriptSignature],
        ]
    ] = []
    introduced_error_signatures: set[TranscriptSignature] = set()
    for event in truth.get("events", []):
        signatures: set[TranscriptSignature] = set()
        for transcript_id in event["target"]["transcript_ids"]:
            model = source.transcripts.get(transcript_id)
            if model is None:
                raise ValueError(
                    f"Hidden-truth transcript is absent from source GFF3: {transcript_id}"
                )
            signatures.add(model.signature)
        if not signatures:
            raise ValueError(f"Event has no scorable transcripts: {event['event_id']}")
        source_gene_ids = event["target"].get(
            "gene_ids", [event["target"]["gene_id"]]
        )
        truth_gene_groups: set[frozenset[TranscriptSignature]] = set()
        for gene_id in source_gene_ids:
            group = source.gene_signatures.get(gene_id)
            if not group:
                raise ValueError(
                    f"Hidden-truth gene is absent from source GFF3: {gene_id}"
                )
            truth_gene_groups.add(group)
        event_errors: set[TranscriptSignature] = set()
        for gene_id in event["target"].get("perturbed_gene_ids", []):
            event_errors.update(perturbed.gene_signatures.get(gene_id, ()))
        frozen_signatures = frozenset(signatures)
        hidden_signatures.update(frozen_signatures)
        introduced_error_signatures.update(event_errors)
        event_truth.append(
            (
                event,
                frozen_signatures,
                frozenset(truth_gene_groups),
                frozenset(event_errors),
            )
        )

    ambiguous_hidden = hidden_signatures & set(perturbed.signatures)
    if ambiguous_hidden:
        examples = sorted(signature.digest for signature in ambiguous_hidden)[:5]
        raise ValueError(
            "Hidden transcript structures remain in the perturbed annotation; "
            f"the benchmark is structurally ambiguous (examples: {examples})"
        )

    candidate_signatures = set(candidate.signatures)
    perturbed_signatures = set(perturbed.signatures)
    raw_novel_candidate = candidate_signatures - perturbed_signatures
    control_novel: set[TranscriptSignature] = set()
    control_retains_source: bool | None = None
    if control_candidate_gff_path is not None:
        control = build_annotation_index(read_gff_document(control_candidate_gff_path))
        control_signatures = set(control.signatures)
        control_novel = control_signatures - set(source.signatures)
        control_retains_source = set(source.signatures) <= control_signatures
        if not control_retains_source:
            raise ValueError(
                "Complete-annotation control candidate removes source structures"
            )
    novel_candidate = raw_novel_candidate - control_novel
    recovered_hidden = candidate_signatures & hidden_signatures
    true_novel = novel_candidate & hidden_signatures
    if recovered_hidden != true_novel:
        raise AssertionError("Recovered hidden structures must be novel to the blind input")

    false_novel = novel_candidate - hidden_signatures
    missed_hidden = hidden_signatures - candidate_signatures
    collateral_baseline = perturbed_signatures - introduced_error_signatures
    collateral_missing = collateral_baseline - candidate_signatures
    exact_candidate_gene_sets = set(candidate.gene_signatures.values())
    hidden_cds_chains = _cds_chain_signatures(hidden_signatures)
    candidate_cds_chains = _cds_chain_signatures(candidate_signatures)
    perturbed_cds_chains = _cds_chain_signatures(perturbed_signatures)
    raw_novel_candidate_cds_chains = candidate_cds_chains - perturbed_cds_chains
    control_novel_cds_chains = _cds_chain_signatures(control_novel)
    novel_candidate_cds_chains = (
        raw_novel_candidate_cds_chains - control_novel_cds_chains
    )
    ambiguous_hidden_cds_chains = hidden_cds_chains & perturbed_cds_chains
    if ambiguous_hidden_cds_chains:
        raise ValueError(
            "Hidden CDS chains remain in the perturbed annotation; "
            "the CDS benchmark is structurally ambiguous"
        )
    exact_candidate_gene_cds_sets = {
        frozen
        for signatures in candidate.gene_signatures.values()
        if (frozen := frozenset(_cds_chain_signatures(signatures)))
    }

    if (event_strata_path is None) != (not stratum_columns):
        raise ValueError(
            "event_strata_path and at least one stratum column must be used together"
        )
    if len(set(stratum_columns)) != len(stratum_columns):
        raise ValueError("Stratum columns must be unique")

    complete_events = 0
    exact_gene_events = 0
    cds_scorable_events = 0
    complete_cds_events = 0
    exact_cds_gene_events = 0
    events_with_introduced_errors = 0
    complete_error_removal_events = 0
    event_details: list[dict[str, Any]] = []
    event_evaluations: list[dict[str, Any]] = []
    for event, signatures, truth_gene_groups, event_errors in event_truth:
        event_cds_chains = frozenset(_cds_chain_signatures(signatures))
        recovered_count = len(signatures & candidate_signatures)
        errors_remaining = event_errors & candidate_signatures
        is_error_removed = not errors_remaining
        is_complete = signatures <= candidate_signatures and is_error_removed
        is_exact_gene = (
            all(group in exact_candidate_gene_sets for group in truth_gene_groups)
            and is_error_removed
        )
        complete_events += int(is_complete)
        exact_gene_events += int(is_exact_gene)
        events_with_introduced_errors += int(bool(event_errors))
        complete_error_removal_events += int(bool(event_errors) and is_error_removed)
        cds_scorable_events += int(bool(event_cds_chains))
        is_complete_cds = (
            bool(event_cds_chains)
            and event_cds_chains <= candidate_cds_chains
            and is_error_removed
        )
        is_exact_cds_gene = bool(event_cds_chains) and (
            all(
                frozenset(_cds_chain_signatures(group))
                in exact_candidate_gene_cds_sets
                for group in truth_gene_groups
            )
            and is_error_removed
        )
        complete_cds_events += int(is_complete_cds)
        exact_cds_gene_events += int(is_exact_cds_gene)
        event_evaluations.append(
            {
                "event": event,
                "signatures": signatures,
                "complete": is_complete,
                "exact_gene": is_exact_gene,
                "complete_cds": is_complete_cds,
                "exact_cds_gene": is_exact_cds_gene,
            }
        )
        if include_event_details:
            event_details.append(
                {
                    "event_id": event["event_id"],
                    "event_type": event["event_type"],
                    "target_gene_id": event["target"]["gene_id"],
                    "target_gene_ids": event["target"].get(
                        "gene_ids", [event["target"]["gene_id"]]
                    ),
                    "truth_transcript_structures": len(signatures),
                    "recovered_transcript_structures": recovered_count,
                    "complete_transcript_recovery": is_complete,
                    "exact_gene_grouping": is_exact_gene,
                    "truth_cds_chains": len(event_cds_chains),
                    "recovered_cds_chains": len(
                        event_cds_chains & candidate_cds_chains
                    ),
                    "complete_cds_chain_recovery": is_complete_cds,
                    "exact_cds_gene_grouping": is_exact_cds_gene,
                    "introduced_error_structures": len(event_errors),
                    "introduced_error_structures_remaining": len(errors_remaining),
                    "complete_error_removal": is_error_removed,
                }
            )

    event_count = len(event_truth)
    tp = len(true_novel)
    fp = len(false_novel)
    fn = len(missed_hidden)
    precision = _rate(tp, tp + fp)
    recall = _rate(tp, tp + fn)
    f1 = _f1(precision, recall)
    cds_chain_metrics = _set_metrics(
        set(hidden_cds_chains), set(novel_candidate_cds_chains)
    )
    (
        hidden_exons,
        hidden_junctions,
        hidden_cds_segments,
    ) = _partial_feature_sets(hidden_signatures)
    (
        candidate_exons,
        candidate_junctions,
        candidate_cds_segments,
    ) = _partial_feature_sets(novel_candidate)
    if control_novel:
        control_exons, control_junctions, control_cds_segments = _partial_feature_sets(
            control_novel
        )
        candidate_exons -= control_exons
        candidate_junctions -= control_junctions
        candidate_cds_segments -= control_cds_segments
    report: dict[str, Any] = {
        "schema_version": SCORE_SCHEMA_VERSION,
        "evaluator": {"name": "PloidyPatch", "version": __version__},
        "inputs": {
            "source_gff_sha256": _file_sha256(source_gff_path),
            "perturbed_gff_sha256": _file_sha256(perturbed_gff_path),
            "candidate_gff_sha256": _file_sha256(candidate_gff_path),
            "hidden_truth_sha256": _file_sha256(truth_path),
        },
        "evaluation_mode": (
            "paired_complete_annotation_difference"
            if control_candidate_gff_path is not None
            else "single_candidate"
        ),
        "strict_transcript_structure": {
            "true_positive": tp,
            "false_positive": fp,
            "false_negative": fn,
            "precision": precision,
            "recall": recall,
            "f1": f1,
            "hidden_truth_structures": len(hidden_signatures),
            "candidate_novel_structures": len(novel_candidate),
        },
        "strict_cds_chain": cds_chain_metrics,
        "partial_structure": {
            "exon_segments": _set_metrics(hidden_exons, candidate_exons),
            "splice_junctions": _set_metrics(
                hidden_junctions, candidate_junctions
            ),
            "cds_segments_with_phase": _set_metrics(
                hidden_cds_segments, candidate_cds_segments
            ),
            "cds_nucleotide_coverage": _cds_coverage_metrics(
                hidden_signatures,
                raw_novel_candidate,
                background_signatures=control_novel,
            ),
        },
        "event_recovery": {
            "events": event_count,
            "complete_transcript_recovery": complete_events,
            "complete_transcript_recall": _rate(complete_events, event_count),
            "exact_gene_grouping": exact_gene_events,
            "exact_gene_recall": _rate(exact_gene_events, event_count),
            "cds_scorable_events": cds_scorable_events,
            "complete_cds_chain_recovery": complete_cds_events,
            "complete_cds_chain_recall": _rate(
                complete_cds_events, cds_scorable_events
            ),
            "exact_cds_gene_grouping": exact_cds_gene_events,
            "exact_cds_gene_recall": _rate(
                exact_cds_gene_events, cds_scorable_events
            ),
            "events_with_introduced_errors": events_with_introduced_errors,
            "complete_error_removal": complete_error_removal_events,
            "complete_error_removal_recall": _rate(
                complete_error_removal_events, events_with_introduced_errors
            ),
        },
        "collateral_changes": {
            "baseline_transcript_structures_missing_from_candidate": len(
                collateral_missing
            )
        },
        "quality_gate": {
            "source_truth_checksum_match": True,
            "perturbed_truth_checksum_match": (
                True if expected_perturbed_sha is not None else None
            ),
            "hidden_structures_absent_from_perturbed": True,
            "introduced_error_structures_identified": len(
                introduced_error_signatures
            ),
            "grade": "pass",
        },
    }
    if control_candidate_gff_path is not None:
        report["inputs"]["control_candidate_gff_sha256"] = _file_sha256(
            control_candidate_gff_path
        )
        report["background_subtraction"] = {
            "candidate_novel_transcript_structures_before_subtraction": len(
                raw_novel_candidate
            ),
            "control_novel_transcript_structures": len(control_novel),
            "shared_background_transcript_structures": len(
                raw_novel_candidate & control_novel
            ),
            "differential_candidate_transcript_structures": len(novel_candidate),
            "candidate_novel_cds_chains_before_subtraction": len(
                raw_novel_candidate_cds_chains
            ),
            "control_novel_cds_chains": len(control_novel_cds_chains),
            "shared_background_cds_chains": len(
                raw_novel_candidate_cds_chains & control_novel_cds_chains
            ),
            "differential_candidate_cds_chains": len(
                novel_candidate_cds_chains
            ),
        }
        report["quality_gate"]["control_retains_all_source_structures"] = (
            control_retains_source
        )
    if include_event_details:
        report["event_details"] = event_details
    if event_strata_path is not None:
        report["stratified_recall"] = _stratified_recall(
            evaluations=event_evaluations,
            novel_candidate=novel_candidate,
            strata_path=event_strata_path,
            columns=stratum_columns,
        )
    return report
