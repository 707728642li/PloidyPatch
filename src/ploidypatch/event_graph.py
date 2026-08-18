from __future__ import annotations

import csv
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path
from typing import Any

from . import __version__


EVENT_GRAPH_SCHEMA_VERSION = "ploidypatch.event_graph.v1"
EVENT_DECISION_SCHEMA_VERSION = "ploidypatch.event_decisions.v1"
EVENT_COHERENCE_SCHEMA_VERSION = "ploidypatch.event_coherence.v0.2"

LOCUS_STATES = frozenset(
    {
        "accepted",
        "split_error",
        "fusion_error",
        "missing_exon",
        "boundary_error",
        "missing_annotation",
        "collapsed_copies",
        "haplotig_duplicate",
        "pseudogene",
        "true_absence",
        "assembly_uncertain",
        "uncertain",
    }
)
RELATIONSHIP_STATES = frozenset(
    {
        "",
        "allele",
        "speciation_ortholog",
        "wgd_homeolog",
        "tandem_duplicate",
        "proximal_duplicate",
        "transposed_duplicate",
        "dispersed_duplicate",
        "lineage_specific",
        "uncertain_relationship",
    }
)
EVIDENCE_DIRECTIONS = frozenset({"support", "contradict", "context"})
EVIDENCE_SCOPES = frozenset({"coding", "isoform", "both", "relationship"})
EVIDENCE_PHASES = frozenset({"before_edit", "after_edit", "shared"})

EVIDENCE_WEIGHTS = {
    "protein_projection": 2.0,
    "cds_chain_concordance": 2.0,
    "domain_concordance": 1.0,
    "synteny": 1.5,
    "wgd_quota": 1.5,
    "subgenome_match": 1.0,
    "rna_junction": 2.0,
    "rna_coverage": 1.0,
    "long_read_full_length": 3.0,
    "assembly_gap": 2.0,
    "read_depth": 1.5,
    "copy_number": 2.0,
    "mappability": 1.0,
    "repeat_overlap": 0.75,
    "splice_site": 1.5,
    "gene_tree": 1.5,
    "pseudogene_lesion": 2.0,
    "whole_genome_alignment_absence": 2.5,
    "pav_deletion": 2.5,
    "haplotype_alignment": 2.0,
}

STRUCTURE_STATES = frozenset(
    {"split_error", "fusion_error", "missing_exon", "boundary_error"}
)


def _file_sha256(path: str | Path, block_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        while block := handle.read(block_size):
            digest.update(block)
    return digest.hexdigest()


def _read_tsv(path: str | Path) -> tuple[list[str], list[dict[str, str]]]:
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"TSV has no header: {path}")
        return list(reader.fieldnames), list(reader)


def _parse_probability(value: str, field: str, line_number: int) -> float:
    try:
        parsed = float(value)
    except ValueError as exc:
        raise ValueError(f"Non-numeric {field} at line {line_number}") from exc
    if not 0.0 <= parsed <= 1.0:
        raise ValueError(f"{field} is outside [0,1] at line {line_number}")
    return parsed


def _candidate_rows(path: str | Path) -> list[dict[str, Any]]:
    fields, rows = _read_tsv(path)
    required = {"event_id", "event_type", "candidate_id", "gene_ids"}
    if not required <= set(fields):
        raise ValueError(
            "Candidate TSV is missing required columns: "
            + ", ".join(sorted(required - set(fields)))
        )
    candidates: list[dict[str, Any]] = []
    seen_candidates: set[str] = set()
    for line_number, row in enumerate(rows, start=2):
        event_id = row["event_id"]
        candidate_id = row["candidate_id"]
        event_type = row["event_type"]
        if not event_id:
            raise ValueError(f"Empty event_id at line {line_number}")
        if not candidate_id or candidate_id in seen_candidates:
            raise ValueError(f"Empty or duplicate candidate_id at line {line_number}")
        if event_type not in LOCUS_STATES:
            raise ValueError(f"Unknown locus state at line {line_number}: {event_type}")
        relationship_state = row.get("relationship_state", "")
        if relationship_state not in RELATIONSHIP_STATES:
            raise ValueError(
                f"Unknown relationship state at line {line_number}: "
                f"{relationship_state}"
            )
        gene_ids = tuple(value for value in row["gene_ids"].split(",") if value)
        if not gene_ids:
            raise ValueError(f"Candidate has no gene_ids at line {line_number}")
        seen_candidates.add(candidate_id)
        candidates.append(
            {
                "event_id": event_id,
                "event_type": event_type,
                "candidate_id": candidate_id,
                "gene_ids": gene_ids,
                "seqid": row.get("seqid", ""),
                "start": row.get("start", ""),
                "end": row.get("end", ""),
                "strand": row.get("strand", ""),
                "relationship_state": relationship_state,
                "wgd_event": row.get("wgd_event", ""),
                "subgenome": row.get("subgenome", ""),
                "proposal_uri": row.get("proposal_uri", ""),
            }
        )
    if not candidates:
        raise ValueError("Candidate TSV contains no events")
    return candidates


def _evidence_rows(
    path: str | Path, candidate_ids: frozenset[str]
) -> list[dict[str, Any]]:
    fields, rows = _read_tsv(path)
    required = {
        "candidate_id",
        "evidence_id",
        "evidence_type",
        "direction",
        "scope",
        "strength",
        "reliability",
        "independent_group",
        "source",
    }
    if not required <= set(fields):
        raise ValueError(
            "Evidence TSV is missing required columns: "
            + ", ".join(sorted(required - set(fields)))
        )
    evidence: list[dict[str, Any]] = []
    seen_ids: set[str] = set()
    for line_number, row in enumerate(rows, start=2):
        candidate_id = row["candidate_id"]
        evidence_id = row["evidence_id"]
        evidence_type = row["evidence_type"]
        direction = row["direction"]
        scope = row["scope"]
        if candidate_id not in candidate_ids:
            raise ValueError(
                f"Evidence references unknown candidate at line {line_number}: "
                f"{candidate_id}"
            )
        if not evidence_id or evidence_id in seen_ids:
            raise ValueError(f"Empty or duplicate evidence_id at line {line_number}")
        if evidence_type not in EVIDENCE_WEIGHTS:
            raise ValueError(
                f"Unknown evidence_type at line {line_number}: {evidence_type}"
            )
        if direction not in EVIDENCE_DIRECTIONS:
            raise ValueError(f"Unknown evidence direction at line {line_number}")
        if scope not in EVIDENCE_SCOPES:
            raise ValueError(f"Unknown evidence scope at line {line_number}")
        if not row["independent_group"] or not row["source"]:
            raise ValueError(
                f"Evidence lacks independent_group or source at line {line_number}"
            )
        evidence_phase = row.get("evidence_phase", "").strip()
        if "evidence_phase" in fields and evidence_phase not in EVIDENCE_PHASES:
            raise ValueError(
                f"Unknown or empty evidence_phase at line {line_number}: "
                f"{evidence_phase!r}"
            )
        strength = _parse_probability(row["strength"], "strength", line_number)
        reliability = _parse_probability(
            row["reliability"], "reliability", line_number
        )
        sign = 1.0 if direction == "support" else -1.0
        contribution = (
            0.0
            if direction == "context"
            else sign * EVIDENCE_WEIGHTS[evidence_type] * strength * reliability
        )
        seen_ids.add(evidence_id)
        evidence.append(
            {
                "candidate_id": candidate_id,
                "evidence_id": evidence_id,
                "evidence_type": evidence_type,
                "direction": direction,
                "scope": scope,
                "strength": strength,
                "reliability": reliability,
                "independent_group": row["independent_group"],
                "source": row["source"],
                "algorithm_family": (
                    row.get("algorithm_family", "").strip() or row["source"]
                ),
                "biological_source_group": (
                    row.get("biological_source_group", "").strip()
                    or row["independent_group"]
                ),
                "evidence_phase": evidence_phase,
                "details": row.get("details", ""),
                "weight": EVIDENCE_WEIGHTS[evidence_type],
                "raw_contribution": contribution,
            }
        )
    return evidence


def _dependency_components(evidence: list[dict[str, Any]]) -> list[list[int]]:
    """Connect evidence that shares any declared dependency axis."""

    parents = list(range(len(evidence)))

    def find(index: int) -> int:
        while parents[index] != index:
            parents[index] = parents[parents[index]]
            index = parents[index]
        return index

    def union(left: int, right: int) -> None:
        left_root = find(left)
        right_root = find(right)
        if left_root != right_root:
            parents[right_root] = left_root

    owners: dict[tuple[str, str], int] = {}
    axes = ("independent_group", "algorithm_family", "biological_source_group")
    for index, edge in enumerate(evidence):
        for axis in axes:
            key = (axis, edge[axis])
            if key in owners:
                union(index, owners[key])
            else:
                owners[key] = index
    components: dict[int, list[int]] = defaultdict(list)
    for index in range(len(evidence)):
        components[find(index)].append(index)
    return sorted(components.values(), key=lambda values: values[0])


def _dependency_scope_score(
    evidence: list[dict[str, Any]], allowed_scopes: frozenset[str]
) -> dict[str, Any]:
    scoped = [
        edge
        for edge in evidence
        if edge["scope"] in allowed_scopes and edge["direction"] != "context"
    ]
    component_contributions: dict[str, float] = {}
    for component_number, indices in enumerate(
        _dependency_components(scoped), start=1
    ):
        values = [scoped[index]["raw_contribution"] for index in indices]
        positive = max((value for value in values if value > 0), default=0.0)
        negative = min((value for value in values if value < 0), default=0.0)
        component_contributions[f"component_{component_number}"] = positive + negative
    adjusted = sum(component_contributions.values())
    raw = sum(edge["raw_contribution"] for edge in scoped)
    return {
        "logit": adjusted if scoped else None,
        "raw_edge_logit": raw if scoped else None,
        "dependency_components": len(component_contributions),
        "contributing_edges": len(scoped),
        "component_contributions": component_contributions,
        "dependency_policy": (
            "connected_components_sharing_independent_group_algorithm_family_"
            "or_biological_source_group_then_max_support_plus_max_contradiction"
        ),
    }


def _coherence_scope(
    evidence: list[dict[str, Any]], allowed_scopes: frozenset[str]
) -> dict[str, Any]:
    before = _dependency_scope_score(
        [
            edge
            for edge in evidence
            if edge["evidence_phase"] in {"before_edit", "shared"}
        ],
        allowed_scopes,
    )
    after = _dependency_scope_score(
        [
            edge
            for edge in evidence
            if edge["evidence_phase"] in {"after_edit", "shared"}
        ],
        allowed_scopes,
    )
    if before["logit"] is None or after["logit"] is None:
        gain = None
        raw_gain = None
        dependency_penalty = None
    else:
        gain = after["logit"] - before["logit"]
        raw_gain = after["raw_edge_logit"] - before["raw_edge_logit"]
        dependency_penalty = raw_gain - gain
    return {
        "before_edit": before,
        "after_edit": after,
        "gain": gain,
        "raw_edge_gain": raw_gain,
        "dependency_adjustment": dependency_penalty,
    }


def _coherence_summary(
    candidate: dict[str, Any], evidence: list[dict[str, Any]]
) -> dict[str, Any]:
    scopes = {
        "coding": _coherence_scope(evidence, frozenset({"coding", "both"})),
        "isoform": _coherence_scope(evidence, frozenset({"isoform", "both"})),
        "relationship": _coherence_scope(evidence, frozenset({"relationship"})),
    }
    required_names = ["coding"]
    if candidate["event_type"] in STRUCTURE_STATES:
        required_names.append("isoform")
    if candidate["relationship_state"]:
        required_names.append("relationship")
    required_gains = [scopes[name]["gain"] for name in required_names]
    available = all(value is not None for value in required_gains)
    overall_gain = (
        sum(value for value in required_gains if value is not None)
        / len(required_gains)
        if available
        else None
    )
    return {
        "schema_version": EVENT_COHERENCE_SCHEMA_VERSION,
        "available": available,
        "required_scopes": required_names,
        "overall_gain": overall_gain,
        "scopes": scopes,
        "claim_boundary": "uncalibrated_review_ranking_not_automatic_approval",
    }


def _scope_score(
    evidence: list[dict[str, Any]], allowed_scopes: frozenset[str]
) -> dict[str, Any]:
    by_group: dict[str, list[float]] = defaultdict(list)
    contributing_edges = 0
    for edge in evidence:
        if edge["scope"] not in allowed_scopes or edge["direction"] == "context":
            continue
        by_group[edge["independent_group"]].append(edge["raw_contribution"])
        contributing_edges += 1
    group_contributions: dict[str, float] = {}
    for group, values in by_group.items():
        positive = max((value for value in values if value > 0), default=0.0)
        negative = min((value for value in values if value < 0), default=0.0)
        group_contributions[group] = positive + negative
    score = sum(group_contributions.values())
    if not by_group:
        confidence = None
    elif score >= 0:
        confidence = 1.0 / (1.0 + math.exp(-score))
    else:
        exp_score = math.exp(score)
        confidence = exp_score / (1.0 + exp_score)
    return {
        "logit": score if by_group else None,
        "confidence": confidence,
        "independent_groups": len(by_group),
        "contributing_edges": contributing_edges,
        "group_contributions": dict(sorted(group_contributions.items())),
        "correlation_policy": "max_support_plus_max_contradiction_per_group",
    }


def _strong(
    evidence: list[dict[str, Any]], evidence_type: str, direction: str
) -> bool:
    return any(
        edge["evidence_type"] == evidence_type
        and edge["direction"] == direction
        and edge["strength"] * edge["reliability"] >= 0.7
        for edge in evidence
    )


def _constraints(
    candidate: dict[str, Any], evidence: list[dict[str, Any]]
) -> list[dict[str, str]]:
    violations: list[dict[str, str]] = []
    state = candidate["event_type"]
    relationship = candidate["relationship_state"]

    if state == "true_absence":
        for evidence_type in ("protein_projection", "long_read_full_length"):
            if _strong(evidence, evidence_type, "contradict"):
                violations.append(
                    {
                        "constraint": "absence_vs_intact_locus_evidence",
                        "severity": "hard",
                        "message": (
                            f"true_absence conflicts with strong {evidence_type} "
                            "support"
                        ),
                    }
                )
        if _strong(evidence, "assembly_gap", "contradict"):
            violations.append(
                {
                    "constraint": "absence_vs_assembly_gap",
                    "severity": "hard",
                    "message": "assembly gap requires assembly_uncertain, not true_absence",
                }
            )
        if not any(
            _strong(evidence, evidence_type, "support")
            for evidence_type in ("whole_genome_alignment_absence", "pav_deletion")
        ):
            violations.append(
                {
                    "constraint": "absence_requires_sequence_evidence",
                    "severity": "hard",
                    "message": "true_absence lacks strong sequence-level absence evidence",
                }
            )

    if relationship == "wgd_homeolog":
        if not candidate["wgd_event"]:
            violations.append(
                {
                    "constraint": "wgd_homeolog_requires_named_event",
                    "severity": "hard",
                    "message": "wgd_homeolog is not attached to a named WGD event",
                }
            )
        if not _strong(evidence, "synteny", "support"):
            violations.append(
                {
                    "constraint": "wgd_homeolog_requires_synteny",
                    "severity": "hard",
                    "message": "wgd_homeolog lacks strong synteny support",
                }
            )

    if state == "collapsed_copies" and not any(
        _strong(evidence, evidence_type, "support")
        for evidence_type in ("copy_number", "read_depth")
    ):
        violations.append(
            {
                "constraint": "collapse_requires_copy_evidence",
                "severity": "hard",
                "message": "collapsed_copies lacks strong copy-number/read-depth support",
            }
        )

    if state == "haplotig_duplicate" and not _strong(
        evidence, "haplotype_alignment", "support"
    ):
        violations.append(
            {
                "constraint": "haplotig_requires_allelic_evidence",
                "severity": "hard",
                "message": "haplotig_duplicate lacks strong haplotype alignment",
            }
        )

    if state in STRUCTURE_STATES | {"missing_annotation"} and not any(
        _strong(evidence, evidence_type, "support")
        for evidence_type in (
            "protein_projection",
            "cds_chain_concordance",
            "rna_junction",
            "long_read_full_length",
        )
    ):
        violations.append(
            {
                "constraint": "repair_requires_structure_evidence",
                "severity": "hard",
                "message": "repair hypothesis lacks strong independent structure evidence",
            }
        )
    if (
        state in STRUCTURE_STATES | {"missing_annotation"}
        and _strong(evidence, "mappability", "contradict")
        and not _strong(evidence, "synteny", "support")
    ):
        violations.append(
            {
                "constraint": "multilocus_repair_requires_synteny",
                "severity": "hard",
                "message": (
                    "multi-locus projection cannot select a plant repair "
                    "without locus-disambiguating synteny"
                ),
            }
        )
    if (
        state in STRUCTURE_STATES | {"missing_annotation"}
        and _strong(evidence, "assembly_gap", "contradict")
    ):
        violations.append(
            {
                "constraint": "repair_requires_resolved_assembly_context",
                "severity": "hard",
                "message": (
                    "repair overlaps or abuts unresolved assembly sequence "
                    "context"
                ),
            }
        )
    return violations


def _decision(
    candidate: dict[str, Any], coding: dict[str, Any], isoform: dict[str, Any],
    relationship: dict[str, Any], violations: list[dict[str, str]]
) -> tuple[str, str]:
    if violations:
        return "abstain", "hard_constraint_violation"
    state = candidate["event_type"]
    required = [coding["confidence"]]
    if state in STRUCTURE_STATES:
        required.append(isoform["confidence"])
    if candidate["relationship_state"]:
        required.append(relationship["confidence"])
    if any(value is None for value in required):
        return "review", "required_evidence_scope_missing"
    values = [value for value in required if value is not None]
    if all(value >= 0.8 for value in values):
        return "accept_high_confidence", "all_required_confidences_at_least_0.8"
    if any(value <= 0.2 for value in values):
        return "reject", "required_confidence_at_most_0.2"
    return "review", "intermediate_or_conflicting_evidence"


def _selection_score(row: dict[str, Any]) -> float | None:
    values = [row["coding_model"]["confidence"]]
    if row["event_type"] in STRUCTURE_STATES:
        values.append(row["isoform_model"]["confidence"])
    if row["relationship_state"]:
        values.append(row["relationship_model"]["confidence"])
    if any(value is None for value in values):
        return None
    defined = [value for value in values if value is not None]
    return sum(defined) / len(defined)


def _select_compatible_alternatives(decisions: list[dict[str, Any]]) -> None:
    by_event: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in decisions:
        row["selected"] = False
        row["selection_score"] = _selection_score(row)
        by_event[row["event_id"]].append(row)
    for alternatives in by_event.values():
        accepted = [
            row for row in alternatives if row["decision"] == "accept_high_confidence"
        ]
        if len(accepted) == 1:
            accepted[0]["selected"] = True
            continue
        if len(accepted) < 2:
            continue
        accepted.sort(
            key=lambda row: (
                -(row["selection_score"] or 0.0),
                row["candidate_id"],
            )
        )
        top_score = accepted[0]["selection_score"] or 0.0
        second_score = accepted[1]["selection_score"] or 0.0
        if top_score - second_score >= 0.1:
            accepted[0]["selected"] = True
            for row in accepted[1:]:
                row["decision"] = "review"
                row["decision_reason"] = "alternative_lost_joint_selection"
        else:
            for row in accepted:
                row["decision"] = "review"
                row["decision_reason"] = "competing_high_confidence_alternatives"


def infer_event_graph(
    candidate_tsv_path: str | Path,
    evidence_tsv_path: str | Path,
    *,
    output_json_path: str | Path | None = None,
    decisions_tsv_path: str | Path | None = None,
) -> dict[str, Any]:
    """Build an auditable evidence graph and constrained transparent scores."""

    candidates = _candidate_rows(candidate_tsv_path)
    candidate_ids = frozenset(row["candidate_id"] for row in candidates)
    evidence = _evidence_rows(evidence_tsv_path, candidate_ids)
    evidence_by_candidate: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for edge in evidence:
        evidence_by_candidate[edge["candidate_id"]].append(edge)

    decisions: list[dict[str, Any]] = []
    for candidate in candidates:
        edges = sorted(
            evidence_by_candidate.get(candidate["candidate_id"], []),
            key=lambda edge: edge["evidence_id"],
        )
        after_edit_edges = [
            edge
            for edge in edges
            if not edge["evidence_phase"]
            or edge["evidence_phase"] in {"after_edit", "shared"}
        ]
        coding = _scope_score(after_edit_edges, frozenset({"coding", "both"}))
        isoform = _scope_score(after_edit_edges, frozenset({"isoform", "both"}))
        relationship = _scope_score(
            after_edit_edges, frozenset({"relationship"})
        )
        violations = _constraints(candidate, after_edit_edges)
        decision, reason = _decision(
            candidate, coding, isoform, relationship, violations
        )
        coherence = _coherence_summary(candidate, edges)
        if violations:
            risk_tier = "abstain"
            risk_reason = "hard_constraint_violation"
        elif not coherence["available"]:
            risk_tier = "abstain"
            risk_reason = "paired_before_after_evidence_missing"
        elif coherence["overall_gain"] is not None and coherence["overall_gain"] > 0:
            risk_tier = "review_ranked"
            risk_reason = "positive_dependency_adjusted_coherence_gain"
        else:
            risk_tier = "abstain"
            risk_reason = "nonpositive_dependency_adjusted_coherence_gain"
        decisions.append(
            {
                **candidate,
                "coding_model": coding,
                "isoform_model": isoform,
                "relationship_model": relationship,
                "event_coherence": coherence,
                "risk_output": {
                    "tier": risk_tier,
                    "reason": risk_reason,
                    "automatic_approval": False,
                },
                "decision": decision,
                "decision_reason": reason,
                "constraint_violations": violations,
                "evidence_edges": edges,
            }
        )

    _select_compatible_alternatives(decisions)

    graph = {
        "schema_version": EVENT_GRAPH_SCHEMA_VERSION,
        "decision_schema_version": EVENT_DECISION_SCHEMA_VERSION,
        "generator": {"name": "PloidyPatch", "version": __version__},
        "inputs": {
            "candidates": {
                "file_name": Path(candidate_tsv_path).name,
                "sha256": _file_sha256(candidate_tsv_path),
            },
            "evidence": {
                "file_name": Path(evidence_tsv_path).name,
                "sha256": _file_sha256(evidence_tsv_path),
            },
        },
        "model": {
            "evidence_weights": dict(sorted(EVIDENCE_WEIGHTS.items())),
            "correlation_policy": "max support plus strongest contradiction per independent group",
            "accept_threshold": 0.8,
            "reject_threshold": 0.2,
            "confidence_status": "transparent_uncalibrated_v1",
            "event_coherence": {
                "schema_version": EVENT_COHERENCE_SCHEMA_VERSION,
                "comparison": "after_edit_minus_before_edit",
                "dependency_axes": [
                    "independent_group",
                    "algorithm_family",
                    "biological_source_group",
                ],
                "automatic_approval": False,
            },
        },
        "counts": {
            "events": len({row["event_id"] for row in candidates}),
            "candidates": len(candidates),
            "evidence_edges": len(evidence),
            "accept_high_confidence": sum(
                row["decision"] == "accept_high_confidence" for row in decisions
            ),
            "review": sum(row["decision"] == "review" for row in decisions),
            "reject": sum(row["decision"] == "reject" for row in decisions),
            "abstain": sum(row["decision"] == "abstain" for row in decisions),
            "selected": sum(row["selected"] for row in decisions),
            "risk_review_ranked": sum(
                row["risk_output"]["tier"] == "review_ranked"
                for row in decisions
            ),
            "risk_abstain": sum(
                row["risk_output"]["tier"] == "abstain" for row in decisions
            ),
        },
        "decisions": decisions,
    }

    if output_json_path is not None:
        output = Path(output_json_path)
        if output.exists():
            raise FileExistsError(f"Refusing to overwrite event graph: {output}")
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(
            json.dumps(graph, indent=2, sort_keys=True, ensure_ascii=False) + "\n",
            encoding="utf-8",
            newline="",
        )
    if decisions_tsv_path is not None:
        output = Path(decisions_tsv_path)
        if output.exists():
            raise FileExistsError(f"Refusing to overwrite decisions TSV: {output}")
        output.parent.mkdir(parents=True, exist_ok=True)
        fields = [
            "event_id",
            "candidate_id",
            "event_type",
            "gene_ids",
            "relationship_state",
            "wgd_event",
            "subgenome",
            "coding_confidence",
            "isoform_confidence",
            "relationship_confidence",
            "decision",
            "decision_reason",
            "selection_score",
            "selected",
            "coherence_available",
            "coherence_gain",
            "risk_tier",
            "risk_reason",
            "automatic_approval",
            "constraint_violation_count",
        ]
        with output.open("x", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            for row in decisions:
                writer.writerow(
                    {
                        "event_id": row["event_id"],
                        "candidate_id": row["candidate_id"],
                        "event_type": row["event_type"],
                        "gene_ids": ",".join(row["gene_ids"]),
                        "relationship_state": row["relationship_state"],
                        "wgd_event": row["wgd_event"],
                        "subgenome": row["subgenome"],
                        "coding_confidence": row["coding_model"]["confidence"],
                        "isoform_confidence": row["isoform_model"]["confidence"],
                        "relationship_confidence": row["relationship_model"][
                            "confidence"
                        ],
                        "decision": row["decision"],
                        "decision_reason": row["decision_reason"],
                        "selection_score": row["selection_score"],
                        "selected": int(row["selected"]),
                        "coherence_available": int(
                            row["event_coherence"]["available"]
                        ),
                        "coherence_gain": row["event_coherence"]["overall_gain"],
                        "risk_tier": row["risk_output"]["tier"],
                        "risk_reason": row["risk_output"]["reason"],
                        "automatic_approval": int(
                            row["risk_output"]["automatic_approval"]
                        ),
                        "constraint_violation_count": len(
                            row["constraint_violations"]
                        ),
                    }
                )
    return graph
