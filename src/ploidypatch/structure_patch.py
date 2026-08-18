from __future__ import annotations

import csv
import json
from collections import defaultdict, deque
from pathlib import Path
from typing import Any, Iterable

from . import __version__
from .baseline import _file_sha256, _merge_intervals
from .perturb import GffDocument, GffRecord, read_gff_document
from .structure_hypothesis import (
    TRANSCRIPT_TYPES,
    CodingChain,
    _coding_chains,
    _is_missing_internal_segment,
    _is_terminal_boundary_change,
)


STRUCTURE_EDIT_SPEC_SCHEMA_VERSION = "ploidypatch.structure_patch_edits.v1"
SUPPORTED_PATCH_EVENTS = frozenset(
    {
        "annotation_boundary_shift",
        "annotation_fused_gene",
        "annotation_missing_internal_exon",
        "annotation_split_gene",
    }
)


def _read_hypotheses(path: str | Path) -> list[dict[str, str]]:
    required = {
        "hypothesis_id",
        "event_type",
        "candidate_group_ids",
        "annotation_gene_ids",
        "annotation_transcript_ids",
    }
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or not required <= set(reader.fieldnames):
            raise ValueError(
                "Structure hypotheses are missing required columns: "
                + ", ".join(sorted(required))
            )
        rows = list(reader)
    ids = [row["hypothesis_id"] for row in rows]
    if any(not value for value in ids) or len(ids) != len(set(ids)):
        raise ValueError("Structure hypothesis IDs must be non-empty and unique")
    return rows


def _split_values(text: str) -> tuple[str, ...]:
    return tuple(sorted(value for value in text.split(",") if value))


def _source_hierarchy(
    document: GffDocument,
    gene_ids: tuple[str, ...],
    expected_transcript_ids: tuple[str, ...],
) -> tuple[GffRecord, ...]:
    records_by_id: dict[str, list[GffRecord]] = defaultdict(list)
    children_by_parent: dict[str, list[GffRecord]] = defaultdict(list)
    for record in document.records:
        if record.feature_id:
            records_by_id[record.feature_id].append(record)
        for parent in record.parents:
            children_by_parent[parent].append(record)

    selected: dict[int, GffRecord] = {}
    selected_ids: set[str] = set()
    transcript_ids: set[str] = set()
    queue: deque[str] = deque(gene_ids)
    while queue:
        feature_id = queue.popleft()
        if feature_id in selected_ids:
            continue
        records = records_by_id.get(feature_id, [])
        if len(records) != 1:
            raise ValueError(
                f"Patch target ID must resolve exactly once: {feature_id}"
            )
        record = records[0]
        selected_ids.add(feature_id)
        selected[record.line_number] = record
        for child in children_by_parent.get(feature_id, ()):  # includes no-ID CDS
            selected[child.line_number] = child
            if child.feature_type in TRANSCRIPT_TYPES and child.feature_id:
                transcript_ids.add(child.feature_id)
            # Leaf features in real plant GFF3 files often reuse one CDS ID
            # across all segments of a transcript.  Only resolve an ID through
            # records_by_id when it is itself a parent of another feature;
            # otherwise the already selected child line is the complete leaf.
            if child.feature_id and child.feature_id in children_by_parent:
                queue.append(child.feature_id)

    if tuple(sorted(transcript_ids)) != expected_transcript_ids:
        raise ValueError(
            "Patch target must contain exactly the hypothesis transcripts; "
            f"observed={sorted(transcript_ids)}, expected={expected_transcript_ids}"
        )
    for gene_id in gene_ids:
        gene = records_by_id[gene_id][0]
        if gene.feature_type != "gene":
            raise ValueError(f"Patch target is not a gene: {gene_id}")
        direct = children_by_parent.get(gene_id, ())
        if not direct or any(
            child.feature_type not in TRANSCRIPT_TYPES for child in direct
        ):
            raise ValueError(
                f"Patch target gene has non-transcript direct children: {gene_id}"
            )
    for record in selected.values():
        if record.parents and not set(record.parents) <= selected_ids:
            raise ValueError(
                f"Patch target contains a shared external parent at line {record.line_number}"
            )
    return tuple(selected[line] for line in sorted(selected))


def _validate_topology(
    event_type: str,
    candidates: tuple[CodingChain, ...],
    annotation: tuple[CodingChain, ...],
) -> None:
    valid = False
    if event_type == "annotation_missing_internal_exon":
        valid = (
            len(candidates) == 1
            and len(annotation) == 1
            and _is_missing_internal_segment(candidates[0], annotation[0])
        )
    elif event_type == "annotation_boundary_shift":
        valid = (
            len(candidates) == 1
            and len(annotation) == 1
            and _is_terminal_boundary_change(candidates[0], annotation[0])
        )
    elif event_type == "annotation_split_gene":
        valid = (
            len(candidates) == 1
            and len(annotation) == 2
            and annotation[0].cds_set.isdisjoint(annotation[1].cds_set)
            and annotation[0].cds_set | annotation[1].cds_set
            == candidates[0].cds_set
        )
    elif event_type == "annotation_fused_gene":
        valid = (
            len(candidates) == 2
            and len(annotation) == 1
            and candidates[0].cds_set.isdisjoint(candidates[1].cds_set)
            and candidates[0].cds_set | candidates[1].cds_set
            == annotation[0].cds_set
        )
    if not valid:
        raise ValueError(f"Hypothesis no longer satisfies exact topology: {event_type}")


def _replacement_lines(
    hypothesis_id: str,
    event_type: str,
    candidates: tuple[CodingChain, ...],
) -> tuple[str, ...]:
    lines: list[str] = []
    token = hypothesis_id.removeprefix("PPH-")
    ordered = sorted(
        candidates,
        key=lambda item: (item.seqid, item.start, item.end, item.candidate_group),
    )
    for part, candidate in enumerate(ordered, start=1):
        suffix = f"{token}_{part}" if len(ordered) > 1 else token
        gene_id = f"PPR_gene_{suffix}"
        transcript_id = f"PPR_tx_{suffix}"
        provenance = (
            f"ploidypatch_hypothesis={hypothesis_id};"
            f"repair_event_type={event_type};"
            f"candidate_group={candidate.candidate_group};"
            f"support_group_count={candidate.support_group_count};"
            f"support_groups={','.join(candidate.support_groups)}"
        )
        prefix = f"{candidate.seqid}\tPloidyPatchRepair"
        lines.append(
            f"{prefix}\tgene\t{candidate.start}\t{candidate.end}\t.\t"
            f"{candidate.strand}\t.\tID={gene_id};{provenance}\n"
        )
        lines.append(
            f"{prefix}\tmRNA\t{candidate.start}\t{candidate.end}\t.\t"
            f"{candidate.strand}\t.\tID={transcript_id};Parent={gene_id};"
            f"{provenance}\n"
        )
        intervals = _merge_intervals(
            [(start, end) for start, end, _ in candidate.cds]
        )
        for exon_number, (start, end) in enumerate(intervals, start=1):
            lines.append(
                f"{prefix}\texon\t{start}\t{end}\t.\t{candidate.strand}\t.\t"
                f"ID=PPR_exon_{suffix}_{exon_number};Parent={transcript_id}\n"
            )
        for start, end, phase in candidate.cds:
            lines.append(
                f"{prefix}\tCDS\t{start}\t{end}\t.\t{candidate.strand}\t"
                f"{phase}\tParent={transcript_id}\n"
            )
    return tuple(lines)


def compile_structure_patch_edits(
    *,
    annotation_gff_path: str | Path,
    candidate_gff_path: str | Path,
    hypotheses_tsv_path: str | Path,
    output_edits_json_path: str | Path,
    allowed_event_types: Iterable[str],
    min_support_group_count: int = 2,
) -> dict[str, Any]:
    """Compile explicitly allowed, exact hypotheses into immutable line edits."""

    output = Path(output_edits_json_path)
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite structure edits: {output}")
    allowed = frozenset(allowed_event_types)
    if not allowed or not allowed <= SUPPORTED_PATCH_EVENTS:
        raise ValueError(
            "At least one supported --event-type must be explicitly allowed"
        )
    if min_support_group_count < 1:
        raise ValueError("min_support_group_count must be positive")

    document = read_gff_document(annotation_gff_path)
    annotation_chains = _coding_chains(
        annotation_gff_path, candidates_only=False
    )
    annotation_by_id = {
        chain.transcript_id: chain for chain in annotation_chains
    }
    candidates = _coding_chains(candidate_gff_path, candidates_only=True)
    candidate_by_group = {chain.candidate_group: chain for chain in candidates}
    rows = _read_hypotheses(hypotheses_tsv_path)

    compiled: list[dict[str, Any]] = []
    operations_by_line: dict[int, dict[str, Any]] = {}
    skipped_counts: dict[str, int] = defaultdict(int)
    for row in rows:
        event_type = row["event_type"]
        if event_type not in allowed:
            skipped_counts["event_type_not_allowed"] += 1
            continue
        support_text = row.get("support_group_count_min") or row.get(
            "support_source_count_min", "0"
        )
        try:
            support_count = int(support_text)
        except ValueError as exc:
            raise ValueError(
                f"Invalid support count for {row['hypothesis_id']}"
            ) from exc
        if support_count < min_support_group_count:
            skipped_counts["support_below_threshold"] += 1
            continue

        candidate_groups = _split_values(row["candidate_group_ids"])
        transcript_ids = _split_values(row["annotation_transcript_ids"])
        gene_ids = _split_values(row["annotation_gene_ids"])
        try:
            selected_candidates = tuple(
                candidate_by_group[group] for group in candidate_groups
            )
            selected_annotation = tuple(
                annotation_by_id[transcript_id]
                for transcript_id in transcript_ids
            )
        except KeyError as exc:
            raise ValueError(
                f"Hypothesis references an absent candidate/transcript: {exc}"
            ) from exc
        if tuple(
            sorted(
                {gene_id for chain in selected_annotation for gene_id in chain.gene_ids}
            )
        ) != gene_ids:
            raise ValueError(
                f"Hypothesis gene/transcript mapping changed: {row['hypothesis_id']}"
            )
        _validate_topology(event_type, selected_candidates, selected_annotation)
        hierarchy = _source_hierarchy(document, gene_ids, transcript_ids)
        target_lines = {record.line_number for record in hierarchy}
        overlap = target_lines & set(operations_by_line)
        if overlap:
            raise ValueError(
                f"Selected hypotheses overlap source lines: {sorted(overlap)}"
            )
        insertion_line = min(target_lines)
        replacements = _replacement_lines(
            row["hypothesis_id"], event_type, selected_candidates
        )
        for record in hierarchy:
            operations_by_line[record.line_number] = {
                "source_line_number": record.line_number,
                "source_raw_line": record.raw_line,
                "replacement_lines": (
                    list(replacements)
                    if record.line_number == insertion_line
                    else []
                ),
            }
        compiled.append(
            {
                "hypothesis_id": row["hypothesis_id"],
                "event_type": event_type,
                "candidate_group_ids": list(candidate_groups),
                "annotation_gene_ids": list(gene_ids),
                "annotation_transcript_ids": list(transcript_ids),
                "support_group_count_min": support_count,
                "source_line_count": len(target_lines),
                "replacement_line_count": len(replacements),
            }
        )

    if not compiled:
        raise ValueError("No structure hypotheses passed the explicit patch policy")
    payload = {
        "schema_version": STRUCTURE_EDIT_SPEC_SCHEMA_VERSION,
        "generator": {"name": "PloidyPatch", "version": __version__},
        "source": {
            "file_name": Path(annotation_gff_path).name,
            "sha256": _file_sha256(annotation_gff_path),
        },
        "inputs": {
            "candidate_gff": {
                "file_name": Path(candidate_gff_path).name,
                "sha256": _file_sha256(candidate_gff_path),
            },
            "hypotheses_tsv": {
                "file_name": Path(hypotheses_tsv_path).name,
                "sha256": _file_sha256(hypotheses_tsv_path),
            },
        },
        "parameters": {
            "allowed_event_types": sorted(allowed),
            "min_support_group_count": min_support_group_count,
            "target_hierarchy_policy": "complete_single_transcript_gene_only",
            "topology_revalidated": True,
        },
        "event_ids": sorted(row["hypothesis_id"] for row in compiled),
        "compiled_hypotheses": compiled,
        "counts": {
            "input_hypotheses": len(rows),
            "compiled_hypotheses": len(compiled),
            "operations": len(operations_by_line),
            "skipped": dict(sorted(skipped_counts.items())),
        },
        "operations": [
            operations_by_line[line_number]
            for line_number in sorted(operations_by_line)
        ],
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x", encoding="utf-8", newline="") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return payload
