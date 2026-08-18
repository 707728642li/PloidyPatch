from __future__ import annotations

import hashlib
import json
import csv
from collections import Counter, defaultdict
from itertools import zip_longest
from pathlib import Path
from typing import Any

from . import __version__
from .io import iter_fasta
from .perturb import (
    MISSING_GENE_EVENT,
    TRUTH_SCHEMA_VERSION,
    _file_sha256,
    read_gff_document,
)
from .score import TranscriptSignature, build_annotation_index


MASK_TRUTH_SCHEMA_VERSION = "ploidypatch.masked_gap_truth.v1"
MASK_MANIFEST_SCHEMA_VERSION = "ploidypatch.masked_gap_manifest.v1"
ABSTENTION_SCORE_SCHEMA_VERSION = "ploidypatch.masked_gap_abstention_score.v2"
MASK_AUDIT_SCHEMA_VERSION = "ploidypatch.masked_gap_genome_audit.v1"
MASK_SELECTION_SUMMARY_SCHEMA_VERSION = (
    "ploidypatch.masked_gap_selection_summary.v1"
)


def _merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    merged: list[list[int]] = []
    for start, end in sorted(intervals):
        if merged and start <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return [(start, end) for start, end in merged]


def _write_masked_genome(
    source_genome_path: str | Path,
    intervals_by_seqid: dict[str, list[tuple[int, int]]],
    output_genome_path: Path,
    output_fai_path: Path,
    line_bases: int = 60,
) -> tuple[dict[str, int], dict[str, str]]:
    lengths: dict[str, int] = {}
    sequence_hashes: dict[str, str] = {}
    seen_mask_seqids: set[str] = set()
    fai_rows: list[str] = []
    with output_genome_path.open("xb") as handle:
        for record_id, header, sequence in iter_fasta(source_genome_path):
            if record_id in lengths:
                raise ValueError(f"Duplicate genome sequence ID: {record_id}")
            sequence_bytes = bytearray(sequence.encode("ascii"))
            for start, end in intervals_by_seqid.get(record_id, []):
                if start < 1 or end > len(sequence_bytes) or end < start:
                    raise ValueError(
                        f"Mask interval is outside {record_id}: {start}-{end}"
                    )
                sequence_bytes[start - 1 : end] = b"N" * (end - start + 1)
                seen_mask_seqids.add(record_id)
            handle.write(f">{header}\n".encode("utf-8"))
            offset = handle.tell()
            for index in range(0, len(sequence_bytes), line_bases):
                handle.write(bytes(sequence_bytes[index : index + line_bases]) + b"\n")
            lengths[record_id] = len(sequence_bytes)
            sequence_hashes[record_id] = hashlib.sha256(sequence_bytes).hexdigest()
            if not sequence_bytes:
                raise ValueError(f"Empty genome sequence: {record_id}")
            first_line_bases = min(line_bases, len(sequence_bytes))
            fai_rows.append(
                f"{record_id}\t{len(sequence_bytes)}\t{offset}\t"
                f"{first_line_bases}\t{first_line_bases + 1}\n"
            )
    missing_seqids = set(intervals_by_seqid) - seen_mask_seqids
    if missing_seqids:
        raise ValueError(
            "Masked sequence ID(s) are absent from genome FASTA: "
            + ", ".join(sorted(missing_seqids))
        )
    with output_fai_path.open("x", encoding="utf-8", newline="") as handle:
        handle.writelines(fai_rows)
    return lengths, sequence_hashes


def create_masked_gap_control(
    *,
    source_genome_path: str | Path,
    hidden_truth_path: str | Path,
    output_genome_path: str | Path,
    output_mask_truth_path: str | Path,
    background_gff_path: str | Path | None = None,
) -> dict[str, Any]:
    """Mask selected annotation-missing loci without changing coordinates."""

    output_genome = Path(output_genome_path)
    output_fai = Path(str(output_genome) + ".fai")
    output_mask_truth = Path(output_mask_truth_path)
    manifest_path = Path(str(output_genome) + ".manifest.json")
    collisions = [
        path
        for path in (output_genome, output_fai, output_mask_truth, manifest_path)
        if path.exists()
    ]
    if collisions:
        raise FileExistsError(
            "Refusing to overwrite masked-gap artifact(s): "
            + ", ".join(str(path) for path in collisions)
        )
    truth = json.loads(Path(hidden_truth_path).read_text(encoding="utf-8"))
    if truth.get("schema_version") != TRUTH_SCHEMA_VERSION:
        raise ValueError("Unsupported hidden-truth schema for genome masking")

    background_by_seqid: dict[str, list[TranscriptSignature]] = defaultdict(list)
    if background_gff_path is not None:
        background = build_annotation_index(read_gff_document(background_gff_path))
        for signature in background.signatures:
            background_by_seqid[signature.seqid].append(signature)

    intervals_by_seqid: dict[str, list[tuple[int, int]]] = defaultdict(list)
    mask_events: list[dict[str, Any]] = []
    excluded_events: list[dict[str, Any]] = []
    for event in truth.get("events", []):
        if event.get("event_type") != MISSING_GENE_EVENT:
            raise ValueError("Masked-gap control supports missing-gene events only")
        target = event["target"]
        start = int(target["start"])
        end = int(target["end"])
        seqid = target["seqid"]
        background_claims = set()
        if background_gff_path is not None:
            for signature in background_by_seqid.get(seqid, []):
                if signature.start > end or signature.end < start:
                    continue
                feature_intervals = set(signature.exons)
                feature_intervals.update(
                    (feature_start, feature_end)
                    for feature_start, feature_end, _ in signature.cds
                )
                if any(
                    feature_start <= end and feature_end >= start
                    for feature_start, feature_end in feature_intervals
                ):
                    background_claims.add(signature)
        if background_claims:
            excluded_events.append(
                {
                    "event_id": event["event_id"],
                    "reason": "preexisting_background_feature_overlap",
                    "overlapping_transcript_structures": len(background_claims),
                }
            )
            continue
        intervals_by_seqid[seqid].append((start, end))
        mask_events.append(
            {
                "event_id": event["event_id"],
                "seqid": seqid,
                "start": start,
                "end": end,
                "span_bp": end - start + 1,
            }
        )
    if not mask_events:
        raise ValueError("Hidden truth contains no eligible maskable events")
    merged_by_seqid = {
        seqid: _merge_intervals(intervals)
        for seqid, intervals in intervals_by_seqid.items()
    }

    output_genome.parent.mkdir(parents=True, exist_ok=True)
    output_mask_truth.parent.mkdir(parents=True, exist_ok=True)
    lengths, sequence_hashes = _write_masked_genome(
        source_genome_path, merged_by_seqid, output_genome, output_fai
    )
    masked_union_bp = sum(
        end - start + 1
        for intervals in merged_by_seqid.values()
        for start, end in intervals
    )
    mask_truth: dict[str, Any] = {
        "schema_version": MASK_TRUTH_SCHEMA_VERSION,
        "generator": {"name": "PloidyPatch", "version": __version__},
        "source": {
            "genome_file_name": Path(source_genome_path).name,
            "genome_sha256": _file_sha256(source_genome_path),
            "hidden_truth_file_name": Path(hidden_truth_path).name,
            "hidden_truth_sha256": _file_sha256(hidden_truth_path),
        },
        "selection": {
            "requested_events": len(truth.get("events", [])),
            "eligible_events": len(mask_events),
            "excluded_events": len(excluded_events),
            "exclusion_rule": (
                "exclude_target_span_overlapping_background_exon_or_cds"
                if background_gff_path is not None
                else None
            ),
            "excluded_event_details": excluded_events,
        },
        "mask": {
            "character": "N",
            "events": len(mask_events),
            "union_bp": masked_union_bp,
            "event_intervals": mask_events,
        },
        "masked_sequence_sha256": sequence_hashes,
    }
    background_source = None
    if background_gff_path is not None:
        background_source = {
            "file_name": Path(background_gff_path).name,
            "sha256": _file_sha256(background_gff_path),
            "role": "exclude_masks_overlapping_retained_exon_or_cds",
        }
        mask_truth["source"]["background_gff"] = background_source
    with output_mask_truth.open("x", encoding="utf-8", newline="") as handle:
        json.dump(mask_truth, handle, indent=2, sort_keys=True)
        handle.write("\n")

    manifest: dict[str, Any] = {
        "schema_version": MASK_MANIFEST_SCHEMA_VERSION,
        "generator": {"name": "PloidyPatch", "version": __version__},
        "control_type": "assembly_gap_masked_gene",
        "source": {
            "genome_file_name": Path(source_genome_path).name,
            "genome_sha256": _file_sha256(source_genome_path),
            "hidden_truth_sha256": _file_sha256(hidden_truth_path),
        },
        "mask": {
            "character": "N",
            "requested_events": len(truth.get("events", [])),
            "events": len(mask_events),
            "excluded_background_overlap_events": len(excluded_events),
            "union_bp": masked_union_bp,
            "sequence_ids": len(merged_by_seqid),
        },
        "assembly": {
            "sequences": len(lengths),
            "total_bp": sum(lengths.values()),
        },
        "outputs": {
            "genome": {
                "file_name": output_genome.name,
                "sha256": _file_sha256(output_genome),
            },
            "fai": {
                "file_name": output_fai.name,
                "sha256": _file_sha256(output_fai),
            },
            "evaluator_mask_truth": {
                "file_name": output_mask_truth.name,
                "sha256": _file_sha256(output_mask_truth),
                "storage": "separate_evaluator_directory",
            },
        },
    }
    if background_source is not None:
        manifest["source"]["background_gff"] = background_source
    with manifest_path.open("x", encoding="utf-8", newline="") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return manifest


def score_masked_gap_abstention(
    *,
    perturbed_gff_path: str | Path,
    candidate_gff_path: str | Path,
    mask_truth_path: str | Path,
    include_event_details: bool = False,
) -> dict[str, Any]:
    """Count novel candidate structures that should abstain at masked loci."""

    mask_truth = json.loads(Path(mask_truth_path).read_text(encoding="utf-8"))
    if mask_truth.get("schema_version") != MASK_TRUTH_SCHEMA_VERSION:
        raise ValueError("Unsupported masked-gap truth schema")
    perturbed = build_annotation_index(read_gff_document(perturbed_gff_path))
    if Path(candidate_gff_path).resolve() == Path(perturbed_gff_path).resolve():
        candidate = perturbed
    else:
        candidate = build_annotation_index(read_gff_document(candidate_gff_path))
    novel = set(candidate.signatures) - set(perturbed.signatures)
    novel_by_seqid: dict[str, list[TranscriptSignature]] = defaultdict(list)
    for signature in novel:
        novel_by_seqid[signature.seqid].append(signature)
    baseline_by_seqid: dict[str, list[TranscriptSignature]] = defaultdict(list)
    for signature in perturbed.signatures:
        baseline_by_seqid[signature.seqid].append(signature)

    intervals = [
        (
            event["event_id"],
            event["seqid"],
            int(event["start"]),
            int(event["end"]),
        )
        for event in mask_truth["mask"]["event_intervals"]
    ]
    spanning_signatures = set()
    claiming_signatures = set()
    claimed_intervals: dict[str, list[tuple[int, int]]] = defaultdict(list)
    baseline_claimed_intervals: dict[str, list[tuple[int, int]]] = defaultdict(list)
    span_crossing_events = 0
    false_repair_events = 0
    preexisting_overlap_events = 0
    clean_false_repair_events = 0
    event_details: list[dict[str, Any]] = []
    for event_id, seqid, start, end in intervals:
        spanning = {
            signature
            for signature in novel_by_seqid.get(seqid, [])
            if signature.start <= end and signature.end >= start
        }
        claiming = set()
        event_claimed_intervals: list[tuple[int, int]] = []
        for signature in spanning:
            feature_intervals = set(signature.exons)
            feature_intervals.update(
                (feature_start, feature_end)
                for feature_start, feature_end, _ in signature.cds
            )
            signature_claims = False
            for feature_start, feature_end in feature_intervals:
                overlap_start = max(start, feature_start)
                overlap_end = min(end, feature_end)
                if overlap_start <= overlap_end:
                    claimed_intervals[seqid].append((overlap_start, overlap_end))
                    event_claimed_intervals.append((overlap_start, overlap_end))
                    signature_claims = True
            if signature_claims:
                claiming.add(signature)
        baseline_claiming = set()
        for signature in baseline_by_seqid.get(seqid, []):
            if signature.start > end or signature.end < start:
                continue
            feature_intervals = set(signature.exons)
            feature_intervals.update(
                (feature_start, feature_end)
                for feature_start, feature_end, _ in signature.cds
            )
            for feature_start, feature_end in feature_intervals:
                overlap_start = max(start, feature_start)
                overlap_end = min(end, feature_end)
                if overlap_start <= overlap_end:
                    baseline_claimed_intervals[seqid].append(
                        (overlap_start, overlap_end)
                    )
                    baseline_claiming.add(signature)
        spanning_signatures.update(spanning)
        claiming_signatures.update(claiming)
        span_crossing_events += int(bool(spanning))
        false_repair_events += int(bool(claiming))
        preexisting_overlap_events += int(bool(baseline_claiming))
        if not baseline_claiming:
            clean_false_repair_events += int(bool(claiming))
        if include_event_details:
            event_details.append(
                {
                    "event_id": event_id,
                    "seqid": seqid,
                    "start": start,
                    "end": end,
                    "candidate_structures_spanning": len(spanning),
                    "candidate_structures_claiming_features": len(claiming),
                    "candidate_feature_bp_inside_mask": sum(
                        interval_end - interval_start + 1
                        for interval_start, interval_end in _merge_intervals(
                            event_claimed_intervals
                        )
                    ),
                    "preexisting_structures_claiming_features": len(
                        baseline_claiming
                    ),
                    "candidate_signature_digests": sorted(
                        signature.digest for signature in claiming
                    ),
                }
            )
    event_count = len(intervals)
    clean_event_count = event_count - preexisting_overlap_events
    unique_claimed_bp = sum(
        interval_end - interval_start + 1
        for values in claimed_intervals.values()
        for interval_start, interval_end in _merge_intervals(values)
    )
    report = {
        "schema_version": ABSTENTION_SCORE_SCHEMA_VERSION,
        "evaluator": {"name": "PloidyPatch", "version": __version__},
        "inputs": {
            "perturbed_gff_sha256": _file_sha256(perturbed_gff_path),
            "candidate_gff_sha256": _file_sha256(candidate_gff_path),
            "mask_truth_sha256": _file_sha256(mask_truth_path),
        },
        "abstention": {
            "masked_events": event_count,
            "events_with_false_repair": false_repair_events,
            "event_false_repair_rate": (
                false_repair_events / event_count if event_count else None
            ),
            "events_abstained": event_count - false_repair_events,
            "event_abstention_rate": (
                (event_count - false_repair_events) / event_count
                if event_count
                else None
            ),
            "novel_structures_with_feature_claim_in_masked_loci": len(
                claiming_signatures
            ),
            "unique_candidate_feature_bp_inside_mask": unique_claimed_bp,
            "events_with_transcript_span_crossing_masked_locus": (
                span_crossing_events
            ),
            "event_transcript_span_crossing_rate": (
                span_crossing_events / event_count if event_count else None
            ),
            "novel_structures_spanning_masked_loci": len(spanning_signatures),
            "all_novel_candidate_structures": len(novel),
            "evaluable_events_without_preexisting_feature_overlap": (
                clean_event_count
            ),
            "evaluable_events_with_false_repair": clean_false_repair_events,
            "evaluable_event_false_repair_rate": (
                clean_false_repair_events / clean_event_count
                if clean_event_count
                else None
            ),
            "evaluable_event_abstention_rate": (
                (clean_event_count - clean_false_repair_events) / clean_event_count
                if clean_event_count
                else None
            ),
        },
        "control_diagnostics": {
            "events_with_preexisting_feature_overlap": preexisting_overlap_events,
            "preexisting_feature_overlap_rate": (
                preexisting_overlap_events / event_count if event_count else None
            ),
            "unique_preexisting_feature_bp_inside_mask": sum(
                interval_end - interval_start + 1
                for values in baseline_claimed_intervals.values()
                for interval_start, interval_end in _merge_intervals(values)
            ),
        },
        "quality_gate": {
            "masked_events_present": event_count > 0,
            "grade": "pass" if event_count > 0 else "fail",
        },
    }
    if include_event_details:
        report["event_details"] = event_details
    return report


def audit_masked_gap_genome(
    *,
    source_genome_path: str | Path,
    masked_genome_path: str | Path,
    mask_truth_path: str | Path,
) -> dict[str, Any]:
    """Verify that a coordinate-stable control differs only by declared N masks."""

    mask_truth = json.loads(Path(mask_truth_path).read_text(encoding="utf-8"))
    if mask_truth.get("schema_version") != MASK_TRUTH_SCHEMA_VERSION:
        raise ValueError("Unsupported masked-gap truth schema")
    intervals_by_seqid: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for event in mask_truth["mask"]["event_intervals"]:
        intervals_by_seqid[event["seqid"]].append(
            (int(event["start"]), int(event["end"]))
        )
    merged_by_seqid = {
        seqid: _merge_intervals(values)
        for seqid, values in intervals_by_seqid.items()
    }

    sequence_ids_match = True
    headers_match = True
    lengths_match = True
    outside_mask_identical = True
    inside_mask_all_n = True
    sequences = 0
    total_bp = 0
    masked_union_bp = 0
    changed_bp = 0
    original_n_inside_mask = 0
    seen_seqids: set[str] = set()
    source_records = iter_fasta(source_genome_path)
    masked_records = iter_fasta(masked_genome_path)
    for source_record, masked_record in zip_longest(source_records, masked_records):
        if source_record is None or masked_record is None:
            sequence_ids_match = False
            lengths_match = False
            continue
        source_id, source_header, source_sequence = source_record
        masked_id, masked_header, masked_sequence = masked_record
        sequences += 1
        total_bp += len(source_sequence)
        if source_id != masked_id:
            sequence_ids_match = False
            continue
        seen_seqids.add(source_id)
        headers_match &= source_header == masked_header
        lengths_match &= len(source_sequence) == len(masked_sequence)
        if len(source_sequence) != len(masked_sequence):
            continue
        cursor = 0
        for start, end in merged_by_seqid.get(source_id, []):
            if start < 1 or end > len(source_sequence) or end < start:
                raise ValueError(
                    f"Mask interval is outside {source_id}: {start}-{end}"
                )
            mask_start = start - 1
            if source_sequence[cursor:mask_start] != masked_sequence[cursor:mask_start]:
                outside_mask_identical = False
            source_slice = source_sequence[mask_start:end]
            masked_slice = masked_sequence[mask_start:end]
            masked_union_bp += end - start + 1
            changed_bp += sum(
                source_base != masked_base
                for source_base, masked_base in zip(source_slice, masked_slice)
            )
            original_n_inside_mask += sum(
                source_base.upper() == "N" for source_base in source_slice
            )
            if any(masked_base.upper() != "N" for masked_base in masked_slice):
                inside_mask_all_n = False
            cursor = end
        if source_sequence[cursor:] != masked_sequence[cursor:]:
            outside_mask_identical = False
    declared_seqids_present = set(merged_by_seqid) <= seen_seqids
    declared_union_bp = int(mask_truth["mask"]["union_bp"])
    union_bp_matches_truth = masked_union_bp == declared_union_bp
    checks = {
        "sequence_ids_match": sequence_ids_match,
        "headers_match": headers_match,
        "lengths_match": lengths_match,
        "declared_sequence_ids_present": declared_seqids_present,
        "outside_mask_identical": outside_mask_identical,
        "inside_mask_all_n": inside_mask_all_n,
        "union_bp_matches_truth": union_bp_matches_truth,
    }
    return {
        "schema_version": MASK_AUDIT_SCHEMA_VERSION,
        "auditor": {"name": "PloidyPatch", "version": __version__},
        "inputs": {
            "source_genome_sha256": _file_sha256(source_genome_path),
            "masked_genome_sha256": _file_sha256(masked_genome_path),
            "mask_truth_sha256": _file_sha256(mask_truth_path),
        },
        "assembly": {"sequences": sequences, "total_bp": total_bp},
        "mask": {
            "declared_union_bp": declared_union_bp,
            "observed_union_bp": masked_union_bp,
            "changed_bp": changed_bp,
            "original_n_inside_mask": original_n_inside_mask,
        },
        "checks": checks,
        "quality_gate": {
            "grade": "pass" if all(checks.values()) else "fail",
            "all_checks_pass": all(checks.values()),
        },
    }


def summarize_masked_gap_selection(
    *,
    mask_truth_path: str | Path,
    hidden_truth_path: str | Path,
    strata_tsv_path: str | Path,
    columns: tuple[str, ...],
) -> dict[str, Any]:
    """Summarize clean-control retention without exposing event identifiers."""

    if not columns or len(set(columns)) != len(columns):
        raise ValueError("At least one unique stratum column is required")
    mask_truth = json.loads(Path(mask_truth_path).read_text(encoding="utf-8"))
    hidden_truth = json.loads(Path(hidden_truth_path).read_text(encoding="utf-8"))
    if mask_truth.get("schema_version") != MASK_TRUTH_SCHEMA_VERSION:
        raise ValueError("Unsupported masked-gap truth schema")
    if hidden_truth.get("schema_version") != TRUTH_SCHEMA_VERSION:
        raise ValueError("Unsupported hidden-truth schema")

    event_to_gene = {
        event["event_id"]: event["target"]["gene_id"]
        for event in hidden_truth.get("events", [])
    }
    eligible_event_ids = {
        event["event_id"] for event in mask_truth["mask"]["event_intervals"]
    }
    unknown_eligible = eligible_event_ids - set(event_to_gene)
    if unknown_eligible:
        raise ValueError("Mask truth contains events absent from hidden truth")

    strata_by_gene: dict[str, tuple[str, ...]] = {}
    with Path(strata_tsv_path).open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = set(reader.fieldnames or [])
        required = {"gene_id", *columns}
        if not required <= fieldnames:
            raise ValueError(
                "Strata TSV is missing required column(s): "
                + ", ".join(sorted(required - fieldnames))
            )
        for line_number, row in enumerate(reader, start=2):
            gene_id = row["gene_id"]
            if gene_id in strata_by_gene:
                raise ValueError(
                    f"Duplicate gene_id in strata TSV at line {line_number}: {gene_id}"
                )
            strata_by_gene[gene_id] = tuple(row[column] for column in columns)

    requested_counts: Counter[tuple[str, ...]] = Counter()
    eligible_counts: Counter[tuple[str, ...]] = Counter()
    missing_genes = []
    for event_id, gene_id in event_to_gene.items():
        values = strata_by_gene.get(gene_id)
        if values is None:
            missing_genes.append(gene_id)
            continue
        requested_counts[values] += 1
        if event_id in eligible_event_ids:
            eligible_counts[values] += 1
    if missing_genes:
        raise ValueError(
            f"Strata TSV is missing {len(missing_genes)} hidden-truth gene(s)"
        )

    groups = []
    for values in sorted(requested_counts):
        requested = requested_counts[values]
        eligible = eligible_counts[values]
        groups.append(
            {
                "values": dict(zip(columns, values)),
                "requested_events": requested,
                "eligible_events": eligible,
                "excluded_events": requested - eligible,
                "retention_rate": eligible / requested,
            }
        )
    return {
        "schema_version": MASK_SELECTION_SUMMARY_SCHEMA_VERSION,
        "summarizer": {"name": "PloidyPatch", "version": __version__},
        "inputs": {
            "mask_truth_sha256": _file_sha256(mask_truth_path),
            "hidden_truth_sha256": _file_sha256(hidden_truth_path),
            "strata_tsv_sha256": _file_sha256(strata_tsv_path),
        },
        "columns": list(columns),
        "totals": {
            "requested_events": sum(requested_counts.values()),
            "eligible_events": sum(eligible_counts.values()),
            "excluded_events": (
                sum(requested_counts.values()) - sum(eligible_counts.values())
            ),
        },
        "groups": groups,
    }
