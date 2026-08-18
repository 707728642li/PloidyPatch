from __future__ import annotations

import csv
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

from . import __version__
from .baseline import _file_sha256
from .perturb import TRUTH_SCHEMA_VERSION, read_gff_document
from .structure_hypothesis import CodingChain, _coding_chains


STRUCTURE_HYPOTHESIS_SCORE_SCHEMA_VERSION = (
    "ploidypatch.structure_hypothesis_score.v1"
)
SUPPORTED_EVENT_TYPES = (
    "annotation_boundary_shift",
    "annotation_fused_gene",
    "annotation_missing_internal_exon",
    "annotation_split_gene",
)


def _rate(numerator: int, denominator: int) -> float | None:
    return numerator / denominator if denominator else None


def _chain_signature(chain: CodingChain) -> tuple[Any, ...]:
    return (chain.seqid, chain.strand, chain.cds)


def _read_hypotheses(path: str | Path) -> list[dict[str, str]]:
    required = {
        "hypothesis_id",
        "event_type",
        "candidate_group_ids",
        "annotation_gene_ids",
    }
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or not required <= set(reader.fieldnames):
            raise ValueError(
                "Structure-hypothesis TSV is missing required columns: "
                + ", ".join(sorted(required))
            )
        rows = list(reader)
    ids = [row["hypothesis_id"] for row in rows]
    if any(not hypothesis_id for hypothesis_id in ids):
        raise ValueError("Structure-hypothesis TSV contains an empty ID")
    if len(ids) != len(set(ids)):
        raise ValueError("Structure-hypothesis TSV contains duplicate IDs")
    return rows


def score_structure_hypotheses(
    *,
    source_gff_path: str | Path,
    perturbed_gff_path: str | Path,
    candidate_gff_path: str | Path,
    hypotheses_tsv_path: str | Path,
    truth_path: str | Path,
    control_hypotheses_tsv_path: str | Path | None = None,
    include_event_details: bool = False,
) -> dict[str, Any]:
    """Score typed audit hypotheses without requiring erroneous-model removal."""

    with Path(truth_path).open("r", encoding="utf-8") as handle:
        truth = json.load(handle)
    if truth.get("schema_version") != TRUTH_SCHEMA_VERSION:
        raise ValueError("Unsupported hidden-truth schema")

    source_document = read_gff_document(source_gff_path)
    perturbed_document = read_gff_document(perturbed_gff_path)
    source_checksum_match = (
        source_document.text_sha256 == truth["source"]["text_sha256"]
    )
    perturbed_checksum_match = (
        perturbed_document.text_sha256
        == truth["perturbation"]["perturbed_text_sha256"]
    )

    source_chains = _coding_chains(source_gff_path, candidates_only=False)
    source_by_transcript = {
        chain.transcript_id: chain for chain in source_chains
    }
    if len(source_by_transcript) != len(source_chains):
        raise ValueError("Source annotation contains duplicate coding transcript IDs")
    candidate_chains = _coding_chains(candidate_gff_path, candidates_only=True)
    candidate_by_group = {
        chain.candidate_group: chain for chain in candidate_chains
    }
    if len(candidate_by_group) != len(candidate_chains):
        raise ValueError("Candidate GFF3 contains duplicate candidate groups")

    hypothesis_rows = _read_hypotheses(hypotheses_tsv_path)
    control_ids: set[str] = set()
    if control_hypotheses_tsv_path is not None:
        control_ids = {
            row["hypothesis_id"]
            for row in _read_hypotheses(control_hypotheses_tsv_path)
        }

    event_catalog: list[dict[str, Any]] = []
    for event in truth["events"]:
        event_type = event["event_type"]
        if event_type not in SUPPORTED_EVENT_TYPES:
            continue
        transcript_ids = tuple(event["target"].get("transcript_ids", ()))
        missing_transcripts = [
            transcript_id
            for transcript_id in transcript_ids
            if transcript_id not in source_by_transcript
        ]
        if missing_transcripts:
            raise ValueError(
                f"Truth references absent source transcript(s): {missing_transcripts}"
            )
        expected_signatures = tuple(
            sorted(
                _chain_signature(source_by_transcript[transcript_id])
                for transcript_id in transcript_ids
            )
        )
        event_catalog.append(
            {
                "event_id": event["event_id"],
                "event_type": event_type,
                "annotation_gene_ids": tuple(
                    sorted(event["target"].get("perturbed_gene_ids", ()))
                ),
                "expected_signatures": expected_signatures,
                "expected_chain_count": len(expected_signatures),
            }
        )

    event_by_id = {event["event_id"]: event for event in event_catalog}
    true_event_by_hypothesis: dict[str, str] = {}
    unresolved_candidate_groups: set[str] = set()
    for row in hypothesis_rows:
        candidate_groups = tuple(
            sorted(
                value
                for value in row["candidate_group_ids"].split(",")
                if value
            )
        )
        missing_groups = [
            group for group in candidate_groups if group not in candidate_by_group
        ]
        unresolved_candidate_groups.update(missing_groups)
        if missing_groups:
            continue
        candidate_signatures = tuple(
            sorted(
                _chain_signature(candidate_by_group[group])
                for group in candidate_groups
            )
        )
        annotation_gene_ids = tuple(
            sorted(
                value
                for value in row["annotation_gene_ids"].split(",")
                if value
            )
        )
        matches = [
            event
            for event in event_catalog
            if event["event_type"] == row["event_type"]
            and event["annotation_gene_ids"] == annotation_gene_ids
            and event["expected_signatures"] == candidate_signatures
        ]
        if len(matches) > 1:
            raise AssertionError(
                f"Hypothesis matches multiple hidden events: {row['hypothesis_id']}"
            )
        if matches:
            true_event_by_hypothesis[row["hypothesis_id"]] = matches[0][
                "event_id"
            ]

    all_ids = {row["hypothesis_id"] for row in hypothesis_rows}
    differential_ids = all_ids - control_ids
    true_ids = set(true_event_by_hypothesis)
    differential_true_ids = differential_ids & true_ids
    recovered_event_ids = {
        true_event_by_hypothesis[hypothesis_id]
        for hypothesis_id in differential_true_ids
    }

    raw_true = len(true_ids)
    differential_true = len(differential_true_ids)
    differential_false = len(differential_ids - true_ids)
    raw_metrics = {
        "hypotheses": len(all_ids),
        "true_positive": raw_true,
        "false_positive": len(all_ids) - raw_true,
        "precision": _rate(raw_true, len(all_ids)),
    }
    differential_metrics = {
        "hypotheses": len(differential_ids),
        "true_positive": differential_true,
        "false_positive": differential_false,
        "precision": _rate(differential_true, len(differential_ids)),
        "shared_control_hypotheses": len(all_ids & control_ids),
        "control_hypotheses": len(control_ids),
    }

    event_type_metrics: dict[str, dict[str, Any]] = {}
    rows_by_id = {row["hypothesis_id"]: row for row in hypothesis_rows}
    for event_type in SUPPORTED_EVENT_TYPES:
        event_ids = {
            event["event_id"]
            for event in event_catalog
            if event["event_type"] == event_type
        }
        hypothesis_ids = {
            hypothesis_id
            for hypothesis_id in differential_ids
            if rows_by_id[hypothesis_id]["event_type"] == event_type
        }
        event_true_ids = hypothesis_ids & true_ids
        event_recovered = {
            true_event_by_hypothesis[hypothesis_id]
            for hypothesis_id in event_true_ids
        }
        expected_chains = sum(
            event_by_id[event_id]["expected_chain_count"] for event_id in event_ids
        )
        recovered_chains = sum(
            event_by_id[event_id]["expected_chain_count"]
            for event_id in event_recovered
        )
        event_type_metrics[event_type] = {
            "events": len(event_ids),
            "recovered_events": len(event_recovered),
            "event_recall": _rate(len(event_recovered), len(event_ids)),
            "expected_cds_chains": expected_chains,
            "recovered_cds_chains": recovered_chains,
            "cds_chain_recall": _rate(recovered_chains, expected_chains),
            "differential_hypotheses": len(hypothesis_ids),
            "true_positive": len(event_true_ids),
            "false_positive": len(hypothesis_ids - true_ids),
            "precision": _rate(len(event_true_ids), len(hypothesis_ids)),
        }

    total_events = len(event_catalog)
    total_expected_chains = sum(
        event["expected_chain_count"] for event in event_catalog
    )
    recovered_chains = sum(
        event_by_id[event_id]["expected_chain_count"]
        for event_id in recovered_event_ids
    )
    quality_pass = (
        source_checksum_match
        and perturbed_checksum_match
        and not unresolved_candidate_groups
    )
    report: dict[str, Any] = {
        "schema_version": STRUCTURE_HYPOTHESIS_SCORE_SCHEMA_VERSION,
        "evaluator": {"name": "PloidyPatch", "version": __version__},
        "evaluation_mode": (
            "paired_complete_annotation_difference"
            if control_hypotheses_tsv_path is not None
            else "raw_hypotheses"
        ),
        "inputs": {
            "source_gff_sha256": _file_sha256(source_gff_path),
            "perturbed_gff_sha256": _file_sha256(perturbed_gff_path),
            "candidate_gff_sha256": _file_sha256(candidate_gff_path),
            "hypotheses_tsv_sha256": _file_sha256(hypotheses_tsv_path),
            "hidden_truth_sha256": _file_sha256(truth_path),
            "control_hypotheses_tsv_sha256": (
                _file_sha256(control_hypotheses_tsv_path)
                if control_hypotheses_tsv_path is not None
                else None
            ),
        },
        "quality_gate": {
            "grade": "pass" if quality_pass else "fail",
            "source_truth_checksum_match": source_checksum_match,
            "perturbed_truth_checksum_match": perturbed_checksum_match,
            "unresolved_candidate_group_count": len(
                unresolved_candidate_groups
            ),
        },
        "raw_proposals": raw_metrics,
        "paired_difference": differential_metrics,
        "event_recovery": {
            "events": total_events,
            "recovered_events": len(recovered_event_ids),
            "event_recall": _rate(len(recovered_event_ids), total_events),
            "expected_cds_chains": total_expected_chains,
            "recovered_cds_chains": recovered_chains,
            "cds_chain_recall": _rate(recovered_chains, total_expected_chains),
        },
        "event_type_metrics": event_type_metrics,
    }
    if include_event_details:
        hypotheses_by_event: dict[str, list[str]] = defaultdict(list)
        for hypothesis_id in differential_true_ids:
            hypotheses_by_event[
                true_event_by_hypothesis[hypothesis_id]
            ].append(hypothesis_id)
        report["event_details"] = [
            {
                "event_id": event["event_id"],
                "event_type": event["event_type"],
                "recovered": event["event_id"] in recovered_event_ids,
                "expected_cds_chains": event["expected_chain_count"],
                "true_hypothesis_ids": sorted(
                    hypotheses_by_event[event["event_id"]]
                ),
            }
            for event in event_catalog
        ]
    return report
