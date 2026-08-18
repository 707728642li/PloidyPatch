from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter, defaultdict, deque
from pathlib import Path
from typing import Any

from . import __version__
from .baseline import _file_sha256
from .perturb import read_gff_document
from .score import build_annotation_index
from .structure_perturb import (
    _simple_genes,
    copy_collapse_partner_ids,
)


COPY_PAIR_SAMPLE_SCHEMA_VERSION = "ploidypatch.copy_pair_sample.v1"
COPY_PAIR_SAMPLE_MANIFEST_SCHEMA_VERSION = (
    "ploidypatch.copy_pair_sample_manifest.v1"
)
SAMPLE_METADATA_FIELDS = (
    "collapsed_gene_id",
    "retained_gene_id",
    "target_seqid",
    "coding_cds_segments",
    "coding_cds_bp",
    "coding_complexity_bin",
    "balanced_selection_rank",
)
COMPLEXITY_BINS = ("one", "two_to_three", "four_to_six", "seven_plus")


def _complexity_bin(segments: int) -> str:
    if segments == 1:
        return "one"
    if segments <= 3:
        return "two_to_three"
    if segments <= 6:
        return "four_to_six"
    return "seven_plus"


def _rank(seed: int, namespace: str, key: str) -> str:
    return hashlib.sha256(f"{seed}\0{namespace}\0{key}".encode()).hexdigest()


def sample_balanced_copy_pairs(
    *,
    source_gff_path: str | Path,
    pair_tsv_path: str | Path,
    output_pair_tsv_path: str | Path,
    decisions_tsv_path: str | Path,
    count: int,
    seed: int,
    balance_group_field: str | None = None,
) -> dict[str, Any]:
    """Select evaluator copy-collapse pairs across chromosomes and complexity."""

    if count < 1:
        raise ValueError("Balanced copy-pair count must be positive")
    source_path = Path(source_gff_path)
    pair_path = Path(pair_tsv_path)
    output_path = Path(output_pair_tsv_path)
    decisions_path = Path(decisions_tsv_path)
    manifest_path = Path(str(output_path) + ".manifest.json")
    for path in (source_path, pair_path):
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f"Missing or empty copy-pair sampler input: {path}")
    collisions = [
        path for path in (output_path, decisions_path, manifest_path) if path.exists()
    ]
    if collisions:
        raise FileExistsError(
            "Refusing to overwrite copy-pair sampling artifact(s): "
            + ", ".join(str(path) for path in collisions)
        )

    with pair_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or not {"gene_id_a", "gene_id_b"} <= set(
            reader.fieldnames
        ):
            raise ValueError("Copy-pair table lacks gene_id_a and gene_id_b")
        input_fields = list(reader.fieldnames)
        input_rows = list(reader)
    if balance_group_field is not None:
        if balance_group_field not in input_fields:
            raise ValueError(
                f"Copy-pair table lacks balance group field: {balance_group_field}"
            )
        if any(not row[balance_group_field] for row in input_rows):
            raise ValueError("Copy-pair balance group values must be nonempty")
    collisions = set(input_fields) & set(SAMPLE_METADATA_FIELDS)
    if collisions:
        raise ValueError(
            "Copy-pair input already contains sampler metadata: "
            + ", ".join(sorted(collisions))
        )
    pair_keys: set[tuple[str, str]] = set()
    for line_number, row in enumerate(input_rows, start=2):
        gene_a, gene_b = row["gene_id_a"], row["gene_id_b"]
        if not gene_a or not gene_b or gene_a == gene_b:
            raise ValueError(f"Invalid copy pair at line {line_number}")
        pair = tuple(sorted((gene_a, gene_b)))
        if pair in pair_keys:
            raise ValueError(f"Duplicate copy pair at line {line_number}")
        pair_keys.add(pair)

    document = read_gff_document(source_path)
    simple_genes = _simple_genes(document)
    annotation_index = build_annotation_index(document)
    phased_cds_chain_counts = Counter(
        (
            model.signature.seqid,
            model.signature.strand,
            model.signature.cds,
        )
        for model in annotation_index.transcripts.values()
        if model.signature.cds
    )
    candidates: list[dict[str, Any]] = []
    decisions: list[dict[str, str]] = []
    decision_by_pair: dict[tuple[str, str], dict[str, str]] = {}
    for row in input_rows:
        gene_a, gene_b = row["gene_id_a"], row["gene_id_b"]
        pair = tuple(sorted((gene_a, gene_b)))
        base_decision = {
            **row,
            **{field: "" for field in SAMPLE_METADATA_FIELDS},
            "status": "rejected",
            "reason": "",
        }
        decision_by_pair[pair] = base_decision
        if gene_a not in simple_genes or gene_b not in simple_genes:
            base_decision["reason"] = "unscorable_single_isoform_cds_model_missing"
            decisions.append(base_decision)
            continue
        collapsed_id, retained_id = copy_collapse_partner_ids(
            gene_a, gene_b, seed=seed
        )
        collapsed = simple_genes[collapsed_id]
        retained = simple_genes[retained_id]
        if collapsed.gene.seqid == retained.gene.seqid:
            base_decision["reason"] = "same_seqid_pair"
            decisions.append(base_decision)
            continue
        collapsed_phased_cds_chain = (
            collapsed.gene.seqid,
            collapsed.gene.strand,
            tuple(
                sorted(
                    set(
                        (record.start, record.end, record.phase)
                        for record in collapsed.cds
                    )
                )
            ),
        )
        if phased_cds_chain_counts[collapsed_phased_cds_chain] != 1:
            base_decision["reason"] = (
                "collapsed_phased_cds_chain_not_unique_in_source"
            )
            decisions.append(base_decision)
            continue
        segments = len(collapsed.cds)
        cds_bp = sum(record.end - record.start + 1 for record in collapsed.cds)
        complexity = _complexity_bin(segments)
        pair_identifier = row.get("pair_id") or f"{pair[0]}|{pair[1]}"
        rank = _rank(seed, "balanced_copy_pair", pair_identifier)
        metadata = {
            "collapsed_gene_id": collapsed_id,
            "retained_gene_id": retained_id,
            "target_seqid": collapsed.gene.seqid,
            "coding_cds_segments": str(segments),
            "coding_cds_bp": str(cds_bp),
            "coding_complexity_bin": complexity,
            "balanced_selection_rank": rank,
        }
        base_decision.update(metadata)
        candidates.append(
            {
                "pair": pair,
                "genes": frozenset(pair),
                "row": row,
                "metadata": metadata,
                "seqid": collapsed.gene.seqid,
                "balance_group": (
                    row[balance_group_field]
                    if balance_group_field is not None
                    else collapsed.gene.seqid
                ),
                "complexity": complexity,
                "rank": rank,
            }
        )
        decisions.append(base_decision)

    queues: dict[str, dict[str, deque[dict[str, Any]]]] = defaultdict(dict)
    by_group: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for candidate in candidates:
        by_group[(candidate["balance_group"], candidate["complexity"])].append(
            candidate
        )
    for (balance_group, complexity), rows in by_group.items():
        rows.sort(key=lambda item: (item["rank"], item["pair"]))
        queues[balance_group][complexity] = deque(rows)
    balance_namespace = (
        "target_seqid" if balance_group_field is None else balance_group_field
    )
    balance_groups = sorted(
        queues,
        key=lambda value: (_rank(seed, balance_namespace, value), value),
    )
    cursors = {group: 0 for group in balance_groups}
    selected: list[dict[str, Any]] = []
    occupied: set[str] = set()

    def take_one(balance_group: str) -> dict[str, Any] | None:
        for offset in range(len(COMPLEXITY_BINS)):
            index = (cursors[balance_group] + offset) % len(COMPLEXITY_BINS)
            complexity = COMPLEXITY_BINS[index]
            queue = queues[balance_group].get(complexity)
            if not queue:
                continue
            while queue:
                candidate = queue.popleft()
                if occupied & candidate["genes"]:
                    decision_by_pair[candidate["pair"]]["reason"] = (
                        "shared_gene_competition"
                    )
                    continue
                cursors[balance_group] = (index + 1) % len(COMPLEXITY_BINS)
                return candidate
        return None

    while len(selected) < count:
        progress = False
        for balance_group in balance_groups:
            candidate = take_one(balance_group)
            if candidate is None:
                continue
            progress = True
            selected.append(candidate)
            occupied.update(candidate["genes"])
            decision = decision_by_pair[candidate["pair"]]
            decision["status"] = "selected"
            decision["reason"] = (
                "balanced_chromosome_complexity_hash_rank"
                if balance_group_field is None
                else "balanced_declared_group_complexity_hash_rank"
            )
            if len(selected) == count:
                break
        if not progress:
            break

    selected_pairs = {candidate["pair"] for candidate in selected}
    for candidate in candidates:
        decision = decision_by_pair[candidate["pair"]]
        if candidate["pair"] not in selected_pairs and not decision["reason"]:
            decision["reason"] = "not_selected_after_balanced_hash_rank"

    output_fields = [*input_fields, *SAMPLE_METADATA_FIELDS]
    decision_fields = [*output_fields, "status", "reason"]
    output_path.parent.mkdir(parents=True, exist_ok=True)
    decisions_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=output_fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for candidate in selected:
            writer.writerow({**candidate["row"], **candidate["metadata"]})
    with decisions_path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=decision_fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(decisions)

    selected_by_seqid = Counter(candidate["seqid"] for candidate in selected)
    selected_by_complexity = Counter(
        candidate["complexity"] for candidate in selected
    )
    selected_by_balance_group = Counter(
        candidate["balance_group"] for candidate in selected
    )
    reason_counts = Counter(row["reason"] for row in decisions)
    manifest: dict[str, Any] = {
        "schema_version": COPY_PAIR_SAMPLE_MANIFEST_SCHEMA_VERSION,
        "pair_schema_version": COPY_PAIR_SAMPLE_SCHEMA_VERSION,
        "generator": {"name": "PloidyPatch", "version": __version__},
        "access": "evaluator_only",
        "inputs": {
            "source_gff": {
                "file_name": source_path.name,
                "sha256": _file_sha256(source_path),
            },
            "pairs": {
                "file_name": pair_path.name,
                "sha256": _file_sha256(pair_path),
                "rows": len(input_rows),
            },
        },
        "parameters": {
            "requested_count": count,
            "seed": seed,
            "collapsed_partner_policy": "same_as_structure_perturbation_sha256_rank",
            "target_chromosome_policy": (
                "global_round_robin"
                if balance_group_field is None
                else "reported_not_primary_balance_axis"
            ),
            "balance_group_field": balance_group_field,
            "balance_group_policy": (
                "target_seqid_global_round_robin"
                if balance_group_field is None
                else "declared_group_global_round_robin"
            ),
            "within_chromosome_policy": "round_robin_fixed_CDS_segment_bins_then_sha256_rank",
            "complexity_bins": list(COMPLEXITY_BINS),
            "require_cross_seqid": True,
            "require_operable_single_isoform_cds_model": True,
            "require_unique_collapsed_phased_cds_chain_in_source": True,
        },
        "counts": {
            "input_pairs": len(input_rows),
            "operable_cross_seqid_pairs": len(candidates),
            "selected_pairs": len(selected),
            "selection_shortfall": max(0, count - len(selected)),
            "selected_by_target_seqid": dict(sorted(selected_by_seqid.items())),
            "selected_by_balance_group": dict(
                sorted(selected_by_balance_group.items())
            ),
            "selected_by_complexity": dict(sorted(selected_by_complexity.items())),
            "decision_reason_counts": dict(sorted(reason_counts.items())),
        },
        "outputs": {
            "selected_pairs": {
                "file_name": output_path.name,
                "sha256": _file_sha256(output_path),
                "rows": len(selected),
            },
            "decisions": {
                "file_name": decisions_path.name,
                "sha256": _file_sha256(decisions_path),
                "rows": len(decisions),
            },
        },
    }
    with manifest_path.open("x", encoding="utf-8", newline="") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return manifest
