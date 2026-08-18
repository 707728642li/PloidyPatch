from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter, defaultdict
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
from typing import Any

from . import __version__
from .baseline import _file_sha256
from .perturb import read_gff_document


STRUCTURE_HYPOTHESIS_SCHEMA_VERSION = "ploidypatch.structure_hypotheses.v1"
TRANSCRIPT_TYPES = frozenset({"mRNA", "transcript"})
BIN_SIZE = 100_000


@dataclass(frozen=True)
class CodingChain:
    transcript_id: str
    gene_ids: tuple[str, ...]
    seqid: str
    strand: str
    start: int
    end: int
    cds: tuple[tuple[int, int, str], ...]
    candidate_group: str = ""
    candidate_category: str = ""
    support_source_count: int = 0
    support_sources: tuple[str, ...] = ()
    support_group_count: int = 0
    support_groups: tuple[str, ...] = ()

    @property
    def cds_set(self) -> frozenset[tuple[int, int, str]]:
        return frozenset(self.cds)

    @property
    def cds_bp(self) -> int:
        return sum(end - start + 1 for start, end, _ in self.cds)


def _coding_chains(
    path: str | Path,
    *,
    candidates_only: bool,
) -> list[CodingChain]:
    document = read_gff_document(path)
    transcript_records = {
        record.feature_id: record
        for record in document.records
        if record.feature_type in TRANSCRIPT_TYPES and record.feature_id
    }
    cds_by_parent: dict[str, list[tuple[int, int, str]]] = defaultdict(list)
    for record in document.records:
        if record.feature_type != "CDS":
            continue
        for parent in record.parents:
            cds_by_parent[parent].append((record.start, record.end, record.phase))

    chains: list[CodingChain] = []
    seen_candidate_groups: set[str] = set()
    for transcript_id, record in transcript_records.items():
        candidate_group = record.attributes.get("candidate_group", "")
        if candidates_only != bool(candidate_group):
            continue
        cds = tuple(sorted(set(cds_by_parent.get(transcript_id, ()))))
        if not cds:
            continue
        if candidate_group in seen_candidate_groups:
            raise ValueError(f"Duplicate candidate group in GFF3: {candidate_group}")
        if candidate_group:
            seen_candidate_groups.add(candidate_group)
        support_text = record.attributes.get("support_source_count", "0")
        try:
            support_source_count = int(support_text)
        except ValueError as exc:
            raise ValueError(
                f"Invalid support_source_count for {transcript_id}: {support_text}"
            ) from exc
        chains.append(
            CodingChain(
                transcript_id=transcript_id,
                gene_ids=record.parents,
                seqid=record.seqid,
                strand=record.strand,
                start=record.start,
                end=record.end,
                cds=cds,
                candidate_group=candidate_group,
                candidate_category=record.attributes.get(
                    "candidate_category", ""
                ),
                support_source_count=support_source_count,
                support_sources=tuple(
                    sorted(
                        source
                        for source in record.attributes.get(
                            "support_sources", ""
                        ).split(",")
                        if source
                    )
                ),
                support_group_count=int(
                    record.attributes.get(
                        "support_group_count", support_text
                    )
                ),
                support_groups=tuple(
                    sorted(
                        group
                        for group in record.attributes.get(
                            "support_groups",
                            record.attributes.get("support_sources", ""),
                        ).split(",")
                        if group
                    )
                ),
            )
        )
    return chains


def _annotation_bins(
    chains: list[CodingChain],
) -> dict[tuple[str, str, int], tuple[int, ...]]:
    mutable: dict[tuple[str, str, int], list[int]] = defaultdict(list)
    for index, chain in enumerate(chains):
        for bin_number in range(
            chain.start // BIN_SIZE, chain.end // BIN_SIZE + 1
        ):
            mutable[(chain.seqid, chain.strand, bin_number)].append(index)
    return {key: tuple(values) for key, values in mutable.items()}


def _nearby_annotation_chains(
    candidate: CodingChain,
    annotation: list[CodingChain],
    bins: dict[tuple[str, str, int], tuple[int, ...]],
) -> list[CodingChain]:
    indices: set[int] = set()
    for bin_number in range(
        candidate.start // BIN_SIZE, candidate.end // BIN_SIZE + 1
    ):
        indices.update(bins.get((candidate.seqid, candidate.strand, bin_number), ()))
    return [
        annotation[index]
        for index in sorted(indices)
        if annotation[index].start <= candidate.end
        and annotation[index].end >= candidate.start
    ]


def _is_missing_internal_segment(
    candidate: CodingChain,
    annotation: CodingChain,
) -> bool:
    candidate_set = candidate.cds_set
    annotation_set = annotation.cds_set
    if len(annotation_set) < 2 or not annotation_set < candidate_set:
        return False
    added = candidate_set - annotation_set
    if len(added) != 1:
        return False
    start, end, _ = next(iter(added))
    return candidate.start < start and end < candidate.end


def _is_terminal_boundary_change(
    candidate: CodingChain,
    annotation: CodingChain,
) -> bool:
    candidate_set = candidate.cds_set
    annotation_set = annotation.cds_set
    if len(candidate_set) != len(annotation_set):
        return False
    candidate_only = candidate_set - annotation_set
    annotation_only = annotation_set - candidate_set
    if len(candidate_only) != 1 or len(annotation_only) != 1:
        return False
    candidate_segment = next(iter(candidate_only))
    annotation_segment = next(iter(annotation_only))
    c_start, c_end, c_phase = candidate_segment
    a_start, a_end, a_phase = annotation_segment
    if c_phase != a_phase or (c_start == a_start) == (c_end == a_end):
        return False
    ordered = sorted(candidate_set)
    return candidate_segment in {ordered[0], ordered[-1]}


def _split_matches(
    candidate: CodingChain,
    nearby: list[CodingChain],
) -> list[tuple[CodingChain, CodingChain]]:
    matches: list[tuple[CodingChain, CodingChain]] = []
    for first, second in combinations(nearby, 2):
        if not first.gene_ids or not second.gene_ids:
            continue
        if set(first.gene_ids).intersection(second.gene_ids):
            continue
        first_set = first.cds_set
        second_set = second.cds_set
        if first_set.isdisjoint(second_set) and first_set | second_set == candidate.cds_set:
            matches.append((first, second))
    return matches


def _hypothesis_id(
    event_type: str,
    candidate_groups: tuple[str, ...],
    annotation_transcripts: tuple[str, ...],
) -> str:
    payload = json.dumps(
        (event_type, candidate_groups, annotation_transcripts),
        separators=(",", ":"),
    )
    digest = hashlib.sha256(payload.encode("utf-8")).hexdigest()[:20]
    return f"PPH-{digest}"


def _hypothesis_row(
    *,
    event_type: str,
    candidates: tuple[CodingChain, ...],
    annotation: tuple[CodingChain, ...],
    topology_rule: str,
    proposed_action: str,
) -> dict[str, Any]:
    candidate_groups = tuple(sorted(item.candidate_group for item in candidates))
    annotation_transcripts = tuple(
        sorted(item.transcript_id for item in annotation)
    )
    annotation_genes = tuple(
        sorted({gene_id for item in annotation for gene_id in item.gene_ids})
    )
    source_sets = [set(item.support_sources) for item in candidates]
    group_sets = [set(item.support_groups) for item in candidates]
    source_intersection = set.intersection(*source_sets) if source_sets else set()
    source_union = set.union(*source_sets) if source_sets else set()
    group_intersection = set.intersection(*group_sets) if group_sets else set()
    group_union = set.union(*group_sets) if group_sets else set()
    return {
        "hypothesis_id": _hypothesis_id(
            event_type, candidate_groups, annotation_transcripts
        ),
        "event_type": event_type,
        "candidate_group_ids": ",".join(candidate_groups),
        "annotation_gene_ids": ",".join(annotation_genes),
        "annotation_transcript_ids": ",".join(annotation_transcripts),
        "candidate_chain_count": len(candidates),
        "support_source_count_min": min(
            item.support_source_count for item in candidates
        ),
        "support_sources_intersection": ",".join(sorted(source_intersection)),
        "support_sources_union": ",".join(sorted(source_union)),
        "support_group_count_min": min(
            item.support_group_count for item in candidates
        ),
        "support_groups_intersection": ",".join(sorted(group_intersection)),
        "support_groups_union": ",".join(sorted(group_union)),
        "topology_rule": topology_rule,
        "proposed_action": proposed_action,
    }


def infer_structure_hypotheses(
    *,
    annotation_gff_path: str | Path,
    candidate_gff_path: str | Path,
    output_tsv_path: str | Path,
    candidate_topology_tsv_path: str | Path,
) -> dict[str, Any]:
    """Infer conservative structure-event hypotheses from exact CDS topology."""

    output = Path(output_tsv_path)
    topology_output = Path(candidate_topology_tsv_path)
    manifest_path = Path(str(output) + ".manifest.json")
    if len({output.resolve(), topology_output.resolve(), manifest_path.resolve()}) != 3:
        raise ValueError("Structure-hypothesis output paths must be distinct")
    collisions = [
        path for path in (output, topology_output, manifest_path) if path.exists()
    ]
    if collisions:
        raise FileExistsError(
            "Refusing to overwrite structure-hypothesis artifact(s): "
            + ", ".join(str(path) for path in collisions)
        )

    annotation = _coding_chains(annotation_gff_path, candidates_only=False)
    candidates = _coding_chains(candidate_gff_path, candidates_only=True)
    bins = _annotation_bins(annotation)
    nearby_by_group: dict[str, list[CodingChain]] = {}
    boundary_matches: dict[str, list[CodingChain]] = defaultdict(list)
    missing_matches: dict[str, list[CodingChain]] = defaultdict(list)
    split_matches: dict[
        str, list[tuple[CodingChain, CodingChain]]
    ] = defaultdict(list)
    fusion_components: dict[str, list[CodingChain]] = defaultdict(list)
    annotation_by_id = {item.transcript_id: item for item in annotation}

    for candidate in candidates:
        nearby = _nearby_annotation_chains(candidate, annotation, bins)
        nearby_by_group[candidate.candidate_group] = nearby
        if candidate.candidate_category == "cds_extension_or_missing_segment":
            for current in nearby:
                if _is_missing_internal_segment(candidate, current):
                    missing_matches[candidate.candidate_group].append(current)
                if _is_terminal_boundary_change(candidate, current):
                    boundary_matches[candidate.candidate_group].append(current)
        elif candidate.candidate_category == "spans_multiple_annotated_genes":
            split_matches[candidate.candidate_group].extend(
                _split_matches(candidate, nearby)
            )
        elif candidate.candidate_category == "alternative_cds_chain_within_gene":
            for current in nearby:
                if candidate.cds_set < current.cds_set:
                    fusion_components[current.transcript_id].append(candidate)

    hypotheses: list[dict[str, Any]] = []
    candidate_hypothesis_counts: Counter[str] = Counter()
    for candidate in candidates:
        group = candidate.candidate_group
        if len(missing_matches[group]) == 1 and not boundary_matches[group]:
            row = _hypothesis_row(
                event_type="annotation_missing_internal_exon",
                candidates=(candidate,),
                annotation=(missing_matches[group][0],),
                topology_rule="annotation_chain_plus_one_internal_cds_segment",
                proposed_action="replace_transcript_cds_chain",
            )
            hypotheses.append(row)
            candidate_hypothesis_counts[group] += 1
        if len(boundary_matches[group]) == 1 and not missing_matches[group]:
            row = _hypothesis_row(
                event_type="annotation_boundary_shift",
                candidates=(candidate,),
                annotation=(boundary_matches[group][0],),
                topology_rule="one_terminal_cds_segment_differs_at_one_boundary",
                proposed_action="replace_transcript_cds_chain",
            )
            hypotheses.append(row)
            candidate_hypothesis_counts[group] += 1
        if len(split_matches[group]) == 1:
            row = _hypothesis_row(
                event_type="annotation_split_gene",
                candidates=(candidate,),
                annotation=split_matches[group][0],
                topology_rule="candidate_chain_equals_union_of_two_gene_chains",
                proposed_action="merge_split_gene_models",
            )
            hypotheses.append(row)
            candidate_hypothesis_counts[group] += 1

    for transcript_id, components in fusion_components.items():
        current = annotation_by_id[transcript_id]
        pairs = [
            (first, second)
            for first, second in combinations(
                sorted(components, key=lambda item: item.candidate_group), 2
            )
            if first.cds_set.isdisjoint(second.cds_set)
            and first.cds_set | second.cds_set == current.cds_set
        ]
        if len(pairs) != 1:
            continue
        first, second = pairs[0]
        row = _hypothesis_row(
            event_type="annotation_fused_gene",
            candidates=(first, second),
            annotation=(current,),
            topology_rule="two_disjoint_candidate_chains_partition_annotation_chain",
            proposed_action="split_fused_gene_model",
        )
        hypotheses.append(row)
        candidate_hypothesis_counts[first.candidate_group] += 1
        candidate_hypothesis_counts[second.candidate_group] += 1

    fusion_component_counts: Counter[str] = Counter()
    for components in fusion_components.values():
        fusion_component_counts.update(
            candidate.candidate_group for candidate in components
        )

    unique_hypotheses = {row["hypothesis_id"]: row for row in hypotheses}
    if len(unique_hypotheses) != len(hypotheses):
        raise AssertionError("Duplicate structure-hypothesis ID generated")
    hypotheses = sorted(
        hypotheses,
        key=lambda row: (
            row["event_type"],
            row["annotation_gene_ids"],
            row["candidate_group_ids"],
        ),
    )

    hypothesis_fields = (
        "hypothesis_id",
        "event_type",
        "candidate_group_ids",
        "annotation_gene_ids",
        "annotation_transcript_ids",
        "candidate_chain_count",
        "support_source_count_min",
        "support_sources_intersection",
        "support_sources_union",
        "support_group_count_min",
        "support_groups_intersection",
        "support_groups_union",
        "topology_rule",
        "proposed_action",
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=hypothesis_fields,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(hypotheses)

    topology_fields = (
        "candidate_group_id",
        "candidate_category",
        "support_source_count",
        "support_sources",
        "support_group_count",
        "support_groups",
        "cds_segments",
        "cds_bp",
        "nearby_annotation_transcripts",
        "missing_internal_match_count",
        "terminal_boundary_match_count",
        "split_pair_match_count",
        "fusion_component_match_count",
        "emitted_hypothesis_count",
        "status",
    )
    topology_output.parent.mkdir(parents=True, exist_ok=True)
    status_counts: Counter[str] = Counter()
    with topology_output.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=topology_fields,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for candidate in sorted(candidates, key=lambda item: item.candidate_group):
            group = candidate.candidate_group
            emitted = candidate_hypothesis_counts[group]
            if emitted:
                status = "hypothesis_emitted"
            elif any(
                (
                    missing_matches[group],
                    boundary_matches[group],
                    split_matches[group],
                )
            ):
                status = "ambiguous_topology"
            elif fusion_component_counts[group]:
                status = "unpaired_fusion_component"
            else:
                status = "no_exact_topology"
            status_counts[status] += 1
            fusion_count = fusion_component_counts[group]
            writer.writerow(
                {
                    "candidate_group_id": group,
                    "candidate_category": candidate.candidate_category,
                    "support_source_count": candidate.support_source_count,
                    "support_sources": ",".join(candidate.support_sources),
                    "support_group_count": candidate.support_group_count,
                    "support_groups": ",".join(candidate.support_groups),
                    "cds_segments": len(candidate.cds),
                    "cds_bp": candidate.cds_bp,
                    "nearby_annotation_transcripts": len(
                        nearby_by_group[group]
                    ),
                    "missing_internal_match_count": len(
                        missing_matches[group]
                    ),
                    "terminal_boundary_match_count": len(
                        boundary_matches[group]
                    ),
                    "split_pair_match_count": len(split_matches[group]),
                    "fusion_component_match_count": fusion_count,
                    "emitted_hypothesis_count": emitted,
                    "status": status,
                }
            )

    event_counts = Counter(row["event_type"] for row in hypotheses)
    manifest = {
        "schema_version": STRUCTURE_HYPOTHESIS_SCHEMA_VERSION,
        "generator": {"name": "PloidyPatch", "version": __version__},
        "inputs": {
            "annotation_gff": {
                "file_name": Path(annotation_gff_path).name,
                "sha256": _file_sha256(annotation_gff_path),
            },
            "candidate_gff": {
                "file_name": Path(candidate_gff_path).name,
                "sha256": _file_sha256(candidate_gff_path),
            },
        },
        "parameters": {
            "ambiguity_policy": "emit_only_unique_exact_topology",
            "truth_access": False,
        },
        "counts": {
            "annotation_coding_transcripts": len(annotation),
            "candidate_chains": len(candidates),
            "hypotheses": len(hypotheses),
            "event_type_counts": dict(sorted(event_counts.items())),
            "candidate_status_counts": dict(sorted(status_counts.items())),
        },
        "outputs": {
            "hypotheses": {
                "file_name": output.name,
                "rows": len(hypotheses),
                "sha256": _file_sha256(output),
            },
            "candidate_topology": {
                "file_name": topology_output.name,
                "rows": len(candidates),
                "sha256": _file_sha256(topology_output),
            },
        },
    }
    with manifest_path.open("x", encoding="utf-8", newline="") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return manifest
