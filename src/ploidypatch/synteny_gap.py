from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from . import __version__
from .gff import parse_attributes
from .perturb import _file_sha256
from .wgdi_summary import BLOCK_HEADER, SOURCE_LABEL


SYNTENY_GAP_SCHEMA_VERSION = "ploidypatch.wgdi_synteny_gap_candidates.v1"
SYNTENY_MODEL_SELECTION_SCHEMA_VERSION = (
    "ploidypatch.synteny_gap_model_selection.v1"
)


@dataclass(frozen=True)
class OrderedGene:
    seqid: str
    gene_id: str
    start: int
    end: int
    strand: str
    order: int


@dataclass(frozen=True)
class OrderedPair:
    query_gene: str
    query_order: int
    target_gene: str
    target_order: int


@dataclass(frozen=True)
class OrderedBlock:
    block_id: str
    score: float
    pvalue: float
    declared_pairs: int
    query_seqid: str
    target_seqid: str
    orientation: str
    pairs: tuple[OrderedPair, ...]


def _read_ordered_gff(
    path: str | Path,
) -> tuple[dict[str, OrderedGene], dict[tuple[str, int], OrderedGene]]:
    by_id: dict[str, OrderedGene] = {}
    by_order: dict[tuple[str, int], OrderedGene] = {}
    with Path(path).open("r", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            fields = raw_line.rstrip("\r\n").split("\t")
            if len(fields) != 6:
                raise ValueError(f"Malformed WGDI GFF line {line_number} in {path}")
            seqid, gene_id, start_text, end_text, strand, order_text = fields
            try:
                start, end, order = int(start_text), int(end_text), int(order_text)
            except ValueError as exc:
                raise ValueError(
                    f"Non-integer WGDI GFF field at line {line_number} in {path}"
                ) from exc
            if start < 1 or end < start or order < 1 or strand not in {"+", "-", "."}:
                raise ValueError(f"Invalid WGDI GFF row at line {line_number} in {path}")
            if gene_id in by_id:
                raise ValueError(f"Duplicate WGDI gene ID: {gene_id}")
            if (seqid, order) in by_order:
                raise ValueError(f"Duplicate WGDI gene order: {seqid}:{order}")
            gene = OrderedGene(seqid, gene_id, start, end, strand, order)
            by_id[gene_id] = gene
            by_order[(seqid, order)] = gene
    if not by_id:
        raise ValueError(f"WGDI GFF contains no genes: {path}")
    return by_id, by_order


def _parse_ordered_blocks(path: str | Path) -> list[OrderedBlock]:
    blocks: list[OrderedBlock] = []
    header: dict[str, str] | None = None
    pairs: list[OrderedPair] = []

    def finish() -> None:
        nonlocal header, pairs
        if header is None:
            return
        block = OrderedBlock(
            block_id=header["block_id"],
            score=float(header["score"]),
            pvalue=float(header["pvalue"]),
            declared_pairs=int(header["pair_count"]),
            query_seqid=header["query_seqid"],
            target_seqid=header["target_seqid"],
            orientation=header["orientation"],
            pairs=tuple(pairs),
        )
        if block.declared_pairs != len(block.pairs):
            raise ValueError(
                f"WGDI block {block.block_id} declares {block.declared_pairs} "
                f"pairs but contains {len(block.pairs)}"
            )
        blocks.append(block)
        header = None
        pairs = []

    with Path(path).open("r", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith("#"):
                finish()
                match = BLOCK_HEADER.fullmatch(line)
                if match is None:
                    raise ValueError(
                        f"Malformed WGDI block header on line {line_number}: {line}"
                    )
                header = match.groupdict()
                continue
            if header is None:
                raise ValueError(
                    f"WGDI pair before block header on line {line_number}"
                )
            fields = line.split()
            if len(fields) != 5:
                raise ValueError(
                    f"Malformed WGDI pair on line {line_number}: expected 5 fields"
                )
            try:
                query_order, target_order = int(fields[1]), int(fields[3])
            except ValueError as exc:
                raise ValueError(
                    f"Non-integer WGDI order on line {line_number}"
                ) from exc
            pairs.append(OrderedPair(fields[0], query_order, fields[2], target_order))
    finish()
    if not blocks:
        raise ValueError(f"WGDI collinearity file contains no blocks: {path}")
    return blocks


def infer_wgdi_synteny_gaps(
    *,
    query_wgdi_gff_path: str | Path,
    target_wgdi_gff_path: str | Path,
    collinearity_path: str | Path,
    source_label: str,
    output_tsv_path: str | Path,
    expected_chromosome_pair_tsv_path: str | Path | None = None,
    max_query_intervening_genes: int = 0,
    min_target_excess_genes: int = 1,
    max_target_gap_genes: int = 5,
    max_query_locus_bp: int = 500000,
) -> dict[str, Any]:
    """Infer target-gene hypotheses from gaps between adjacent syntenic anchors."""

    if SOURCE_LABEL.fullmatch(source_label) is None:
        raise ValueError(f"Unsafe source label: {source_label}")
    if (
        max_query_intervening_genes < 0
        or min_target_excess_genes < 1
        or max_target_gap_genes < 1
        or max_query_locus_bp < 1
    ):
        raise ValueError("Synteny-gap thresholds must be positive")
    output_tsv = Path(output_tsv_path)
    manifest_path = Path(str(output_tsv) + ".manifest.json")
    collisions = [path for path in (output_tsv, manifest_path) if path.exists()]
    if collisions:
        raise FileExistsError(
            "Refusing to overwrite synteny-gap artifact(s): "
            + ", ".join(str(path) for path in collisions)
        )

    query_by_id, _ = _read_ordered_gff(query_wgdi_gff_path)
    target_by_id, target_by_order = _read_ordered_gff(target_wgdi_gff_path)
    blocks = _parse_ordered_blocks(collinearity_path)
    expected_pairs: dict[str, tuple[str, str]] | None = None
    if expected_chromosome_pair_tsv_path is not None:
        expected_pairs = {}
        with Path(expected_chromosome_pair_tsv_path).open(
            "r", encoding="utf-8", newline=""
        ) as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            required = {"query_seqid", "target_seqid", "source_label"}
            missing = required - set(reader.fieldnames or [])
            if reader.fieldnames is None or missing:
                raise ValueError(
                    "Expected chromosome-pair TSV is missing column(s): "
                    + ", ".join(sorted(missing))
                )
            for line_number, row in enumerate(reader, start=2):
                query_seqid = row["query_seqid"]
                if not query_seqid or query_seqid in expected_pairs:
                    raise ValueError(
                        f"Empty or duplicate expected query at line {line_number}"
                    )
                expected_pairs[query_seqid] = (
                    row["source_label"],
                    row["target_seqid"],
                )
        if not expected_pairs:
            raise ValueError("Expected chromosome-pair TSV contains no pairs")
    counts: Counter[str] = Counter()
    candidates: dict[tuple[str, int, int, str], dict[str, Any]] = {}

    for block in blocks:
        counts["blocks_total"] += 1
        for pair in block.pairs:
            query_gene = query_by_id.get(pair.query_gene)
            target_gene = target_by_id.get(pair.target_gene)
            if query_gene is None or target_gene is None:
                raise ValueError(
                    f"WGDI block {block.block_id} references a gene absent from GFF"
                )
            if (
                query_gene.seqid != block.query_seqid
                or query_gene.order != pair.query_order
                or target_gene.seqid != block.target_seqid
                or target_gene.order != pair.target_order
            ):
                raise ValueError(
                    f"WGDI block/GFF disagreement in block {block.block_id}"
                )
        if expected_pairs is not None and expected_pairs.get(block.query_seqid) != (
            source_label,
            block.target_seqid,
        ):
            counts["blocks_outside_expected_chromosome_pair"] += 1
            continue
        counts["blocks_expected_chromosome_pair"] += 1
        for left_pair, right_pair in zip(block.pairs, block.pairs[1:]):
            counts["adjacent_anchor_pairs_total"] += 1
            query_gap = abs(right_pair.query_order - left_pair.query_order) - 1
            target_gap = abs(right_pair.target_order - left_pair.target_order) - 1
            if query_gap < 0 or target_gap < 0:
                raise ValueError(f"Repeated anchor order in block {block.block_id}")
            if query_gap > max_query_intervening_genes:
                counts["anchor_pairs_query_gap_too_large"] += 1
                continue
            if target_gap > max_target_gap_genes:
                counts["anchor_pairs_target_gap_too_large"] += 1
                continue
            target_excess = target_gap - query_gap
            if target_excess < min_target_excess_genes:
                counts["anchor_pairs_without_target_excess"] += 1
                continue

            left_query = query_by_id[left_pair.query_gene]
            right_query = query_by_id[right_pair.query_gene]
            physical_left, physical_right = sorted(
                (left_query, right_query), key=lambda gene: (gene.start, gene.end)
            )
            locus_start = physical_left.end + 1
            locus_end = physical_right.start - 1
            if locus_end < locus_start:
                counts["anchor_pairs_without_physical_gap"] += 1
                continue
            locus_span = locus_end - locus_start + 1
            if locus_span > max_query_locus_bp:
                counts["anchor_pairs_locus_too_large"] += 1
                continue

            low_target_order = min(left_pair.target_order, right_pair.target_order)
            high_target_order = max(left_pair.target_order, right_pair.target_order)
            gap_genes = []
            for target_order in range(low_target_order + 1, high_target_order):
                gene = target_by_order.get((block.target_seqid, target_order))
                if gene is None:
                    raise ValueError(
                        f"Missing target gene order {block.target_seqid}:{target_order}"
                    )
                gap_genes.append(gene)
            if len(gap_genes) != target_gap:
                raise AssertionError("Target gap expansion is inconsistent")
            gap_material = (
                f"{source_label}\t{block.query_seqid}\t{locus_start}\t{locus_end}\t"
                f"{left_pair.query_gene}\t{right_pair.query_gene}"
            )
            gap_id = "PPSG-" + hashlib.sha256(
                gap_material.encode("utf-8")
            ).hexdigest()[:20]
            for target_gene in gap_genes:
                counts["raw_gene_hypotheses"] += 1
                key = (block.query_seqid, locus_start, locus_end, target_gene.gene_id)
                existing = candidates.get(key)
                block_key = f"{source_label}:{block.block_id}"
                if existing is None:
                    candidate_material = gap_material + "\t" + target_gene.gene_id
                    candidates[key] = {
                        "candidate_id": "PPSGC-"
                        + hashlib.sha256(candidate_material.encode("utf-8")).hexdigest()[:20],
                        "gap_id": gap_id,
                        "source_label": source_label,
                        "query_seqid": block.query_seqid,
                        "locus_start": locus_start,
                        "locus_end": locus_end,
                        "locus_span_bp": locus_span,
                        "query_left_anchor_gene": physical_left.gene_id,
                        "query_left_anchor_order": physical_left.order,
                        "query_right_anchor_gene": physical_right.gene_id,
                        "query_right_anchor_order": physical_right.order,
                        "query_intervening_genes": query_gap,
                        "target_seqid": block.target_seqid,
                        "target_gene_id": target_gene.gene_id,
                        "target_gene_order": target_gene.order,
                        "target_gene_start": target_gene.start,
                        "target_gene_end": target_gene.end,
                        "target_gene_strand": target_gene.strand,
                        "target_gap_genes": target_gap,
                        "target_excess_genes": target_excess,
                        "best_block_id": block_key,
                        "best_block_score": block.score,
                        "best_block_pvalue": block.pvalue,
                        "best_block_pairs": len(block.pairs),
                        "best_block_orientation": block.orientation,
                        "supporting_block_ids": {block_key},
                    }
                else:
                    existing["supporting_block_ids"].add(block_key)
                    current_rank = (
                        existing["best_block_score"],
                        -existing["best_block_pvalue"],
                        existing["best_block_pairs"],
                    )
                    new_rank = (block.score, -block.pvalue, len(block.pairs))
                    if new_rank > current_rank:
                        existing.update(
                            {
                                "best_block_id": block_key,
                                "best_block_score": block.score,
                                "best_block_pvalue": block.pvalue,
                                "best_block_pairs": len(block.pairs),
                                "best_block_orientation": block.orientation,
                            }
                        )

    fieldnames = [
        "candidate_id",
        "gap_id",
        "source_label",
        "query_seqid",
        "locus_start",
        "locus_end",
        "locus_span_bp",
        "query_left_anchor_gene",
        "query_left_anchor_order",
        "query_right_anchor_gene",
        "query_right_anchor_order",
        "query_intervening_genes",
        "target_seqid",
        "target_gene_id",
        "target_gene_order",
        "target_gene_start",
        "target_gene_end",
        "target_gene_strand",
        "target_gap_genes",
        "target_excess_genes",
        "best_block_id",
        "best_block_score",
        "best_block_pvalue",
        "best_block_pairs",
        "best_block_orientation",
        "supporting_block_count",
        "supporting_block_ids",
    ]
    rows = []
    for candidate in candidates.values():
        support_ids = sorted(candidate.pop("supporting_block_ids"))
        candidate["supporting_block_count"] = len(support_ids)
        candidate["supporting_block_ids"] = ",".join(support_ids)
        candidate["best_block_score"] = f"{candidate['best_block_score']:.8f}"
        candidate["best_block_pvalue"] = f"{candidate['best_block_pvalue']:.8g}"
        rows.append(candidate)
    rows.sort(
        key=lambda row: (
            row["query_seqid"],
            int(row["locus_start"]),
            row["target_seqid"],
            int(row["target_gene_order"]),
        )
    )
    counts["unique_gene_hypotheses"] = len(rows)
    counts["unique_gap_loci"] = len({row["gap_id"] for row in rows})
    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    with output_tsv.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    manifest = {
        "schema_version": SYNTENY_GAP_SCHEMA_VERSION,
        "generator": {"name": "PloidyPatch", "version": __version__},
        "inputs": {
            "query_wgdi_gff": {
                "file_name": Path(query_wgdi_gff_path).name,
                "sha256": _file_sha256(query_wgdi_gff_path),
                "genes": len(query_by_id),
            },
            "target_wgdi_gff": {
                "file_name": Path(target_wgdi_gff_path).name,
                "sha256": _file_sha256(target_wgdi_gff_path),
                "genes": len(target_by_id),
            },
            "collinearity": {
                "file_name": Path(collinearity_path).name,
                "sha256": _file_sha256(collinearity_path),
                "blocks": len(blocks),
            },
        },
        "parameters": {
            "source_label": source_label,
            "max_query_intervening_genes": max_query_intervening_genes,
            "min_target_excess_genes": min_target_excess_genes,
            "max_target_gap_genes": max_target_gap_genes,
            "max_query_locus_bp": max_query_locus_bp,
        },
        "counts": dict(sorted(counts.items())),
        "output": {
            "file_name": output_tsv.name,
            "rows": len(rows),
            "sha256": _file_sha256(output_tsv),
        },
    }
    if expected_chromosome_pair_tsv_path is not None:
        manifest["inputs"]["expected_chromosome_pairs"] = {
            "file_name": Path(expected_chromosome_pair_tsv_path).name,
            "sha256": _file_sha256(expected_chromosome_pair_tsv_path),
            "pairs": len(expected_pairs or {}),
        }
    with manifest_path.open("x", encoding="utf-8", newline="") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return manifest


def select_synteny_gap_models(
    *,
    gap_tsv_paths: list[str | Path],
    baseline_decisions_tsv_path: str | Path,
    adapted_candidate_gff_path: str | Path,
    output_selection_tsv_path: str | Path,
    output_candidate_gff_path: str | Path,
) -> dict[str, Any]:
    """Select accepted protein models fully contained in blind synteny gaps."""

    if not gap_tsv_paths:
        raise ValueError("At least one synteny-gap TSV is required")
    output_selection = Path(output_selection_tsv_path)
    output_candidate = Path(output_candidate_gff_path)
    manifest_path = Path(str(output_selection) + ".manifest.json")
    collisions = [
        path
        for path in (output_selection, output_candidate, manifest_path)
        if path.exists()
    ]
    if collisions:
        raise FileExistsError(
            "Refusing to overwrite synteny-model artifact(s): "
            + ", ".join(str(path) for path in collisions)
        )

    gap_fields = [
        "candidate_id",
        "gap_id",
        "source_label",
        "query_seqid",
        "locus_start",
        "locus_end",
        "locus_span_bp",
        "query_left_anchor_gene",
        "query_right_anchor_gene",
        "target_seqid",
        "target_gene_id",
        "target_gene_order",
        "target_gap_genes",
        "target_excess_genes",
        "best_block_id",
        "best_block_score",
        "best_block_pvalue",
        "best_block_pairs",
        "supporting_block_count",
        "supporting_block_ids",
    ]
    gaps: list[dict[str, str]] = []
    gap_input_manifest = []
    gap_candidate_ids: set[str] = set()
    for raw_path in gap_tsv_paths:
        path = Path(raw_path)
        rows_read = 0
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            missing = set(gap_fields) - set(reader.fieldnames or [])
            if reader.fieldnames is None or missing:
                raise ValueError(
                    f"Gap TSV {path} is missing column(s): "
                    + ", ".join(sorted(missing))
                )
            for line_number, row in enumerate(reader, start=2):
                candidate_id = row["candidate_id"]
                if not candidate_id or candidate_id in gap_candidate_ids:
                    raise ValueError(
                        f"Empty or duplicate gap candidate at {path}:{line_number}"
                    )
                gap_candidate_ids.add(candidate_id)
                gaps.append(row)
                rows_read += 1
        gap_input_manifest.append(
            {
                "file_name": path.name,
                "sha256": _file_sha256(path),
                "rows": rows_read,
            }
        )

    decision_required = {
        "model_id",
        "query_id",
        "source",
        "seqid",
        "start",
        "end",
        "strand",
        "score",
        "rank",
        "identity",
        "query_coverage",
        "frameshifts",
        "stop_codons",
        "status",
        "reason",
    }
    accepted_by_query: dict[str, list[dict[str, str]]] = {}
    accepted_model_ids: set[str] = set()
    decision_rows = 0
    with Path(baseline_decisions_tsv_path).open(
        "r", encoding="utf-8", newline=""
    ) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = decision_required - set(reader.fieldnames or [])
        if reader.fieldnames is None or missing:
            raise ValueError(
                "Baseline decisions are missing column(s): "
                + ", ".join(sorted(missing))
            )
        for line_number, row in enumerate(reader, start=2):
            decision_rows += 1
            model_id = row["model_id"]
            if not model_id:
                raise ValueError(f"Empty baseline model ID at line {line_number}")
            if row["status"] != "accepted":
                continue
            if model_id in accepted_model_ids:
                raise ValueError(f"Duplicate accepted model ID: {model_id}")
            accepted_model_ids.add(model_id)
            accepted_by_query.setdefault(row["query_id"], []).append(row)

    selected_by_model: dict[str, dict[str, str]] = {}
    counts: Counter[str] = Counter()
    for gap in gaps:
        counts["gap_gene_hypotheses_total"] += 1
        query_id = f"{gap['source_label']}__{gap['target_gene_id']}"
        local_models = []
        locus_start = int(gap["locus_start"])
        locus_end = int(gap["locus_end"])
        for model in accepted_by_query.get(query_id, []):
            if model["seqid"] != gap["query_seqid"]:
                continue
            if int(model["start"]) < locus_start or int(model["end"]) > locus_end:
                continue
            local_models.append(model)
        if not local_models:
            counts["gap_gene_hypotheses_without_local_accepted_model"] += 1
            continue
        counts["gap_gene_hypotheses_with_local_accepted_model"] += 1
        best = max(
            local_models,
            key=lambda model: (
                -int(model["rank"]),
                float(model["identity"]),
                float(model["query_coverage"]),
                float(model["score"]),
            ),
        )
        selection = {
            **{field: gap[field] for field in gap_fields},
            "query_id": query_id,
            "model_id": best["model_id"],
            "model_seqid": best["seqid"],
            "model_start": best["start"],
            "model_end": best["end"],
            "model_strand": best["strand"],
            "model_score": best["score"],
            "model_rank": best["rank"],
            "model_identity": best["identity"],
            "model_query_coverage": best["query_coverage"],
            "model_frameshifts": best["frameshifts"],
            "model_stop_codons": best["stop_codons"],
        }
        previous = selected_by_model.get(best["model_id"])
        if previous is None:
            selected_by_model[best["model_id"]] = selection
        else:
            previous_rank = (
                int(previous["supporting_block_count"]),
                float(previous["best_block_score"]),
                -float(previous["best_block_pvalue"]),
            )
            selection_rank = (
                int(selection["supporting_block_count"]),
                float(selection["best_block_score"]),
                -float(selection["best_block_pvalue"]),
            )
            if selection_rank > previous_rank:
                selected_by_model[best["model_id"]] = selection
            counts["duplicate_model_hypotheses"] += 1

    selection_fields = [
        *gap_fields,
        "query_id",
        "model_id",
        "model_seqid",
        "model_start",
        "model_end",
        "model_strand",
        "model_score",
        "model_rank",
        "model_identity",
        "model_query_coverage",
        "model_frameshifts",
        "model_stop_codons",
    ]
    selections = sorted(
        selected_by_model.values(),
        key=lambda row: (
            row["model_seqid"],
            int(row["model_start"]),
            row["model_id"],
        ),
    )
    counts["selected_models"] = len(selections)
    output_selection.parent.mkdir(parents=True, exist_ok=True)
    output_candidate.parent.mkdir(parents=True, exist_ok=True)
    with output_selection.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=selection_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(selections)

    selected_ids = set(selected_by_model)
    baseline_features_total = 0
    baseline_features_selected = 0
    selected_gene_features = 0
    current_selected: bool | None = None
    with Path(adapted_candidate_gff_path).open(
        "r", encoding="utf-8", newline=""
    ) as source, output_candidate.open(
        "x", encoding="utf-8", newline=""
    ) as destination:
        for line_number, raw_line in enumerate(source, start=1):
            line = raw_line.rstrip("\r\n")
            if not line or line.startswith("#"):
                destination.write(raw_line)
                continue
            fields = line.split("\t")
            if len(fields) != 9:
                raise ValueError(
                    f"Malformed adapted candidate GFF line {line_number}"
                )
            if fields[1] != "PloidyPatchBaseline":
                destination.write(raw_line)
                continue
            baseline_features_total += 1
            if fields[2] == "gene":
                attributes, malformed = parse_attributes(fields[8])
                if malformed or "miniprot_model" not in attributes:
                    raise ValueError(
                        f"Baseline gene lacks miniprot_model at line {line_number}"
                    )
                current_selected = attributes["miniprot_model"] in selected_ids
                selected_gene_features += int(current_selected)
            elif current_selected is None:
                raise ValueError(
                    f"Baseline child feature precedes gene at line {line_number}"
                )
            if current_selected:
                destination.write(raw_line)
                baseline_features_selected += 1
    if selected_gene_features != len(selected_ids):
        raise ValueError(
            "Selected model/GFF gene count mismatch: "
            f"{len(selected_ids)} != {selected_gene_features}"
        )

    manifest = {
        "schema_version": SYNTENY_MODEL_SELECTION_SCHEMA_VERSION,
        "generator": {"name": "PloidyPatch", "version": __version__},
        "inputs": {
            "synteny_gaps": gap_input_manifest,
            "baseline_decisions": {
                "file_name": Path(baseline_decisions_tsv_path).name,
                "sha256": _file_sha256(baseline_decisions_tsv_path),
                "rows": decision_rows,
                "accepted_models": len(accepted_model_ids),
            },
            "adapted_candidate_gff": {
                "file_name": Path(adapted_candidate_gff_path).name,
                "sha256": _file_sha256(adapted_candidate_gff_path),
            },
        },
        "policy": {
            "baseline_status_required": "accepted",
            "model_must_be_fully_contained_between_syntenic_anchors": True,
            "one_best_model_per_gap_gene_hypothesis": True,
            "deduplicate_by_miniprot_model_id": True,
        },
        "counts": {
            **dict(sorted(counts.items())),
            "adapted_baseline_features_total": baseline_features_total,
            "adapted_baseline_features_selected": baseline_features_selected,
        },
        "outputs": {
            "selection": {
                "file_name": output_selection.name,
                "rows": len(selections),
                "sha256": _file_sha256(output_selection),
            },
            "candidate_gff": {
                "file_name": output_candidate.name,
                "selected_gene_models": selected_gene_features,
                "sha256": _file_sha256(output_candidate),
            },
        },
    }
    with manifest_path.open("x", encoding="utf-8", newline="") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return manifest
