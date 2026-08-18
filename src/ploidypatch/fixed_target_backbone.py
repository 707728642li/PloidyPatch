from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import shutil
from bisect import bisect_right
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

from . import __version__
from .artifact_manifest import verify_sha256sums, write_sha256sums
from .gff import TRANSCRIPT_TYPES, parse_attributes
from .io import normalize_feature_id, open_text
from .wgdi_summary import WgdiBlock, _file_sha256, parse_wgdi_collinearity


FIXED_TARGET_BACKBONE_SCHEMA_VERSION = "ploidypatch.fixed_target_backbone.v0.6"
FIXED_TARGET_BACKBONE_GENE_SCHEMA_VERSION = (
    "ploidypatch.fixed_target_backbone_genes.v0.6"
)
FIXED_TARGET_BACKBONE_EDGE_SCHEMA_VERSION = (
    "ploidypatch.fixed_target_backbone_edges.v0.6"
)
FIXED_TARGET_BACKBONE_CELL_SCHEMA_VERSION = (
    "ploidypatch.fixed_target_backbone_cells.v0.6"
)
FIXED_TARGET_BACKBONE_SEMANTIC_INVARIANTS = {
    "input_permutation_deterministic": True,
    "reverse_duplicate_statistics_audit_only": True,
    "ambiguous_endpoint_block_quarantine": True,
}


@dataclass(frozen=True)
class FixedTargetBackbonePolicy:
    primary_seqids: tuple[str, ...]
    min_block_pairs: int = 20
    require_cross_seqid: bool = True
    require_reciprocal_unique: bool = True

    def normalized(self) -> dict[str, Any]:
        if self.min_block_pairs < 2:
            raise ValueError("min_block_pairs must be at least two")
        if self.require_cross_seqid is not True:
            raise ValueError("v0.6 fixed target backbone requires cross-seqid pairs")
        if self.require_reciprocal_unique is not True:
            raise ValueError("v0.6 fixed target backbone requires reciprocal uniqueness")
        if not self.primary_seqids:
            raise ValueError("primary_seqids must not be empty")
        if any(
            not isinstance(seqid, str)
            or not seqid
            or any(character in seqid for character in "\t\r\n")
            for seqid in self.primary_seqids
        ):
            raise ValueError("primary_seqids contain an empty or unsafe value")
        if len(set(self.primary_seqids)) != len(self.primary_seqids):
            raise ValueError("primary_seqids contain duplicates")
        return {
            "primary_seqids": sorted(self.primary_seqids),
            "min_block_pairs": self.min_block_pairs,
            "require_cross_seqid": True,
            "require_reciprocal_unique": True,
            "candidate_pair_evidence": "forbidden",
            "relationship_scope": "truth_blind_fixed_target_self_synteny",
        }


@dataclass(frozen=True)
class _BaseGene:
    gene_id: str
    wgdi_gene_id: str
    seqid: str
    start: int
    end: int
    strand: str
    coding_start: int
    coding_end: int
    gene_order: int

    @property
    def midpoint2(self) -> int:
        return self.start + self.end


@dataclass(frozen=True)
class _CanonicalBlock:
    canonical_block_id: str
    seqid_a: str
    seqid_b: str
    orientation: str
    pairs: tuple[tuple[str, str], ...]


def _canonical_json(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True)


def _digest_id(prefix: str, value: Any) -> str:
    digest = hashlib.sha256(_canonical_json(value).encode("utf-8")).hexdigest()[:24]
    return f"{prefix}{digest}"


def _reject_symlink(path: Path, *, label: str) -> None:
    absolute = path.absolute()
    for component in (absolute, *absolute.parents):
        if component.exists() and component.is_symlink():
            raise ValueError(f"{label} traverses a symbolic link: {component}")


def _required_regular_file(path: str | Path, *, label: str) -> Path:
    value = Path(path)
    _reject_symlink(value, label=label)
    if not value.is_file() or value.stat().st_size == 0:
        raise ValueError(f"Missing or empty {label}: {value}")
    return value


def _parse_query_gff(
    path: Path, primary_seqids: set[str]
) -> dict[str, tuple[str, int, int, str, int]]:
    genes: dict[str, tuple[str, int, int, str, int]] = {}
    orders_by_seqid: dict[str, set[int]] = defaultdict(set)
    with path.open("r", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            fields = raw_line.rstrip("\r\n").split("\t")
            if len(fields) != 6:
                raise ValueError(f"Malformed WGDI query GFF line {line_number}")
            seqid, gene_id, start_text, end_text, strand, order_text = fields
            try:
                start, end, order = int(start_text), int(end_text), int(order_text)
            except ValueError as exc:
                raise ValueError(
                    f"Non-integer WGDI query coordinate/order on line {line_number}"
                ) from exc
            if (
                not gene_id
                or gene_id in genes
                or seqid not in primary_seqids
                or start < 1
                or end < start
                or strand not in {"+", "-"}
                or order < 1
                or order in orders_by_seqid[seqid]
            ):
                raise ValueError(f"Invalid or duplicate WGDI query gene on line {line_number}")
            genes[gene_id] = (seqid, start, end, strand, order)
            orders_by_seqid[seqid].add(order)
    if not genes:
        raise ValueError("WGDI query GFF contains no primary-seqid genes")
    for seqid, orders in orders_by_seqid.items():
        if orders != set(range(1, len(orders) + 1)):
            raise ValueError(f"WGDI query gene orders are not contiguous on {seqid}")
    return genes


def _map_wgdi_ids_to_feature_ids(
    path: Path, wgdi_gene_ids: set[str]
) -> tuple[dict[str, str], dict[str, Any]]:
    exact_aliases: dict[str, set[str]] = defaultdict(set)
    normalized_ids: dict[str, set[str]] = defaultdict(set)
    gene_records = 0
    with open_text(path) as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\r\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) != 9:
                raise ValueError(f"Malformed base GFF3 line {line_number}")
            if fields[2] != "gene":
                continue
            attributes, malformed = parse_attributes(fields[8])
            feature_id = attributes.get("ID", "")
            if malformed or not feature_id:
                raise ValueError(f"Invalid base gene attributes on line {line_number}")
            gene_records += 1
            for key in ("ID", "gene_id", "locus_tag", "Name"):
                alias = attributes.get(key, "")
                if alias:
                    exact_aliases[alias].add(feature_id)
            normalized_ids[normalize_feature_id(feature_id)].add(feature_id)

    mapping: dict[str, str] = {}
    modes: Counter[str] = Counter()
    problems: dict[str, Any] = {"missing": [], "ambiguous": {}}
    for wgdi_gene_id in sorted(wgdi_gene_ids):
        matches = exact_aliases.get(wgdi_gene_id, set())
        mode = "exact_attribute"
        if not matches:
            matches = normalized_ids.get(wgdi_gene_id, set())
            mode = "normalized_feature_ID_fallback"
        if len(matches) == 1:
            mapping[wgdi_gene_id] = next(iter(matches))
            modes[mode] += 1
        elif not matches:
            problems["missing"].append(wgdi_gene_id)
        else:
            problems["ambiguous"][wgdi_gene_id] = sorted(matches)
    if problems["missing"] or problems["ambiguous"]:
        raise ValueError(
            "WGDI gene IDs do not map uniquely to base gene features: "
            + json.dumps(
                {
                    "missing": problems["missing"][:10],
                    "ambiguous": dict(list(problems["ambiguous"].items())[:10]),
                },
                sort_keys=True,
            )
        )
    reverse: dict[str, list[str]] = defaultdict(list)
    for wgdi_gene_id, feature_id in mapping.items():
        reverse[feature_id].append(wgdi_gene_id)
    collisions = {
        feature_id: sorted(values)
        for feature_id, values in reverse.items()
        if len(values) > 1
    }
    if collisions:
        raise ValueError(
            "Multiple WGDI gene IDs map to one base gene feature: "
            + json.dumps(dict(list(sorted(collisions.items()))[:10]), sort_keys=True)
        )
    return mapping, {
        "file_name": path.name,
        "sha256": _file_sha256(path),
        "gene_records": gene_records,
        "mapped_wgdi_genes": len(mapping),
        "alias_keys": ["ID", "gene_id", "locus_tag", "Name"],
        "mapping_policy": (
            "unique_exact_attribute_then_unique_normalized_feature_ID_fallback"
        ),
        "mapping_mode_counts": dict(sorted(modes.items())),
        "bidirectional_unique": True,
    }


def _parse_base_genes(
    path: Path,
    query_genes: dict[str, tuple[str, int, int, str, int]],
    feature_id_by_wgdi_id: dict[str, str],
) -> dict[str, _BaseGene]:
    gene_meta: dict[str, tuple[str, int, int, str]] = {}
    transcript_meta: dict[str, tuple[str, str, tuple[str, ...]]] = {}
    cds_rows: list[tuple[int, str, int, int, str, tuple[str, ...]]] = []
    with path.open("r", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\r\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) != 9:
                raise ValueError(f"Malformed base GFF3 line {line_number}")
            seqid, _, feature_type, start_text, end_text, _, strand, phase, raw = fields
            try:
                start, end = int(start_text), int(end_text)
            except ValueError as exc:
                raise ValueError(f"Non-integer base GFF3 coordinate on line {line_number}") from exc
            attributes, malformed = parse_attributes(raw)
            if malformed or start < 1 or end < start or strand not in {"+", "-", "."}:
                raise ValueError(f"Invalid base GFF3 row on line {line_number}")
            if feature_type == "gene":
                gene_id = attributes.get("ID", "")
                if not gene_id or gene_id in gene_meta:
                    raise ValueError(f"Empty or duplicate base gene ID on line {line_number}")
                gene_meta[gene_id] = (seqid, start, end, strand)
            elif feature_type in TRANSCRIPT_TYPES:
                transcript_id = attributes.get("ID", "")
                parents = tuple(value for value in attributes.get("Parent", "").split(",") if value)
                if not transcript_id or transcript_id in transcript_meta or len(parents) != 1:
                    raise ValueError(f"Invalid base transcript on line {line_number}")
                transcript_meta[transcript_id] = (seqid, strand, parents)
            elif feature_type == "CDS":
                parents = tuple(value for value in attributes.get("Parent", "").split(",") if value)
                if not parents or phase not in {"0", "1", "2", "."}:
                    raise ValueError(f"Invalid base CDS on line {line_number}")
                cds_rows.append((line_number, seqid, start, end, strand, parents))

    for transcript_id, (seqid, strand, parents) in transcript_meta.items():
        parent = parents[0]
        if parent not in gene_meta:
            raise ValueError(f"Transcript {transcript_id} references an unknown gene")
        gene_seqid, gene_start, gene_end, gene_strand = gene_meta[parent]
        if seqid != gene_seqid or strand != gene_strand:
            raise ValueError(f"Transcript {transcript_id} disagrees with its gene locus")

    coding_bounds: dict[str, list[int]] = {}
    for line_number, seqid, start, end, strand, parents in cds_rows:
        resolved_genes: set[str] = set()
        for parent in parents:
            if parent in transcript_meta:
                transcript_seqid, transcript_strand, gene_parents = transcript_meta[parent]
                if seqid != transcript_seqid or strand != transcript_strand:
                    raise ValueError(f"CDS disagrees with transcript on line {line_number}")
                resolved_genes.update(gene_parents)
            elif parent in gene_meta:
                resolved_genes.add(parent)
            else:
                raise ValueError(f"CDS references an unknown parent on line {line_number}")
        if len(resolved_genes) != 1:
            raise ValueError(f"CDS maps to multiple genes on line {line_number}")
        gene_id = next(iter(resolved_genes))
        gene_seqid, gene_start, gene_end, gene_strand = gene_meta[gene_id]
        if (
            seqid != gene_seqid
            or strand != gene_strand
            or start < gene_start
            or end > gene_end
        ):
            raise ValueError(f"CDS lies outside its gene locus on line {line_number}")
        bounds = coding_bounds.setdefault(gene_id, [start, end])
        bounds[0] = min(bounds[0], start)
        bounds[1] = max(bounds[1], end)

    output: dict[str, _BaseGene] = {}
    for wgdi_gene_id, query in query_genes.items():
        gene_id = feature_id_by_wgdi_id[wgdi_gene_id]
        if gene_id not in coding_bounds:
            raise ValueError(
                f"WGDI query gene lacks CDS in base GFF3: {wgdi_gene_id} -> {gene_id}"
            )
        seqid, start, end, strand, order = query
        if gene_meta[gene_id] != (seqid, start, end, strand):
            raise ValueError(f"WGDI query/base locus mismatch for gene {gene_id}")
        coding_start, coding_end = coding_bounds[gene_id]
        output[gene_id] = _BaseGene(
            gene_id=gene_id,
            wgdi_gene_id=wgdi_gene_id,
            seqid=seqid,
            start=start,
            end=end,
            strand=strand,
            coding_start=coding_start,
            coding_end=coding_end,
            gene_order=order,
        )
    return output


def _canonicalize_blocks(
    blocks: Iterable[WgdiBlock],
    genes: dict[str, _BaseGene],
    feature_id_by_wgdi_id: dict[str, str],
    policy: dict[str, Any],
) -> tuple[list[_CanonicalBlock], dict[str, int], dict[str, Any]]:
    seen_source_ids: set[str] = set()
    merged: dict[str, dict[str, Any]] = {}
    counts: Counter[str] = Counter()
    statistic_observations: list[tuple[float, float]] = []
    allowlist = set(policy["primary_seqids"])
    for block in blocks:
        if block.block_id in seen_source_ids:
            raise ValueError(f"Duplicate WGDI source block ID: {block.block_id}")
        seen_source_ids.add(block.block_id)
        if block.query_seqid not in allowlist or block.target_seqid not in allowlist:
            raise ValueError(f"WGDI block uses a seqid outside the primary allowlist: {block.block_id}")
        if block.query_seqid == block.target_seqid:
            counts["same_seqid_blocks_skipped"] += 1
            continue
        if len(block.pairs) < policy["min_block_pairs"]:
            counts["short_blocks_skipped"] += 1
            continue
        if not math.isfinite(block.score) or not math.isfinite(block.pvalue):
            raise ValueError(f"WGDI block has a non-finite score/pvalue: {block.block_id}")
        missing = {
            gene for pair in block.pairs for gene in pair
        } - set(feature_id_by_wgdi_id)
        if missing:
            raise ValueError(
                f"WGDI block {block.block_id} references unknown genes: "
                + ", ".join(sorted(missing)[:10])
            )
        mapped_pairs = tuple(
            (feature_id_by_wgdi_id[query_gene], feature_id_by_wgdi_id[target_gene])
            for query_gene, target_gene in block.pairs
        )
        query_counts = Counter(query_gene for query_gene, _ in mapped_pairs)
        target_counts = Counter(target_gene for _, target_gene in mapped_pairs)
        if any(value > 1 for value in query_counts.values()) or any(
            value > 1 for value in target_counts.values()
        ):
            counts["ambiguous_endpoint_blocks_quarantined"] += 1
            counts["ambiguous_endpoint_rows_quarantined"] += len(mapped_pairs)
            continue
        statistic_observations.append((block.score, block.pvalue))
        for query_gene, target_gene in mapped_pairs:
            if (
                genes[query_gene].seqid != block.query_seqid
                or genes[target_gene].seqid != block.target_seqid
            ):
                raise ValueError(f"WGDI block header/anchor seqids disagree: {block.block_id}")

        if block.query_seqid < block.target_seqid:
            seqid_a, seqid_b = block.query_seqid, block.target_seqid
            pairs = mapped_pairs
        else:
            seqid_a, seqid_b = block.target_seqid, block.query_seqid
            pairs = tuple((target, query) for query, target in mapped_pairs)
        ordered_pairs = tuple(
            sorted(
                pairs,
                key=lambda pair: (
                    genes[pair[0]].gene_order,
                    genes[pair[1]].gene_order,
                    pair[0],
                    pair[1],
                ),
            )
        )
        semantic = {
            "seqid_a": seqid_a,
            "seqid_b": seqid_b,
            "orientation": block.orientation,
            "pairs": ordered_pairs,
        }
        canonical_id = _digest_id("PPTB-B-", semantic)
        existing = merged.get(canonical_id)
        if existing is None:
            merged[canonical_id] = {
                **semantic,
                "statistics": {(block.score, block.pvalue)},
            }
        else:
            if (
                existing["seqid_a"] != seqid_a
                or existing["seqid_b"] != seqid_b
                or existing["orientation"] != block.orientation
                or existing["pairs"] != ordered_pairs
            ):
                raise AssertionError("Canonical WGDI block digest collision")
            existing["statistics"].add((block.score, block.pvalue))
            if len(existing["statistics"]) > 1:
                counts["duplicate_geometry_statistic_discordances"] += 1
            counts["reverse_or_exact_duplicate_blocks_collapsed"] += 1
    canonical = [
        _CanonicalBlock(
            canonical_block_id=block_id,
            seqid_a=value["seqid_a"],
            seqid_b=value["seqid_b"],
            orientation=value["orientation"],
            pairs=value["pairs"],
        )
        for block_id, value in sorted(merged.items())
    ]
    counts["source_blocks"] = len(seen_source_ids)
    counts["canonical_qualifying_blocks"] = len(canonical)
    discordant_groups = [
        {
            "canonical_block_id": block_id,
            "distinct_statistics": [
                {
                    "score": format(score, ".17g"),
                    "pvalue": format(pvalue, ".17g"),
                }
                for score, pvalue in sorted(value["statistics"])
            ],
        }
        for block_id, value in sorted(merged.items())
        if len(value["statistics"]) > 1
    ]
    scores = [score for score, _ in statistic_observations]
    pvalues = [pvalue for _, pvalue in statistic_observations]
    statistic_audit: dict[str, Any] = {
        "semantic_use": False,
        "source_qualifying_block_observations": len(statistic_observations),
        "score_min": format(min(scores), ".17g") if scores else None,
        "score_max": format(max(scores), ".17g") if scores else None,
        "pvalue_min": format(min(pvalues), ".17g") if pvalues else None,
        "pvalue_max": format(max(pvalues), ".17g") if pvalues else None,
        "duplicate_geometry_discordant_group_count": len(discordant_groups),
        "duplicate_geometry_discordant_groups": discordant_groups,
    }
    return canonical, dict(sorted(counts.items())), statistic_audit


def _backbone_id(
    genes: dict[str, _BaseGene],
    blocks: list[_CanonicalBlock],
    policy: dict[str, Any],
) -> str:
    value = {
        "policy": policy,
        "genes": [
            {
                "gene_id": gene.gene_id,
                "wgdi_gene_id": gene.wgdi_gene_id,
                "seqid": gene.seqid,
                "start": gene.start,
                "end": gene.end,
                "strand": gene.strand,
                "coding_start": gene.coding_start,
                "coding_end": gene.coding_end,
                "gene_order": gene.gene_order,
            }
            for gene in sorted(genes.values(), key=lambda item: item.gene_id)
        ],
        "blocks": [
            {
                "id": block.canonical_block_id,
                "seqid_a": block.seqid_a,
                "seqid_b": block.seqid_b,
                "orientation": block.orientation,
                "pairs": block.pairs,
            }
            for block in blocks
        ],
    }
    return _digest_id("PPTB-", value)


def _gene_rows(backbone_id: str, genes: dict[str, _BaseGene]) -> list[dict[str, Any]]:
    return [
        {
            "backbone_id": backbone_id,
            "gene_id": gene.gene_id,
            "wgdi_gene_id": gene.wgdi_gene_id,
            "seqid": gene.seqid,
            "start": gene.start,
            "end": gene.end,
            "strand": gene.strand,
            "gene_order": gene.gene_order,
            "midpoint2": gene.midpoint2,
            "coding_start": gene.coding_start,
            "coding_end": gene.coding_end,
        }
        for gene in sorted(
            genes.values(), key=lambda item: (item.seqid, item.gene_order, item.gene_id)
        )
    ]


def _edge_rows(
    backbone_id: str,
    genes: dict[str, _BaseGene],
    blocks: list[_CanonicalBlock],
) -> list[dict[str, Any]]:
    evidence: dict[tuple[str, str], dict[str, Any]] = {}
    for block in blocks:
        for left, right in block.pairs:
            pair = tuple(sorted((left, right)))
            item = evidence.setdefault(
                pair,
                {
                    "blocks": set(),
                    "longest": len(block.pairs),
                },
            )
            item["blocks"].add(block.canonical_block_id)
            item["longest"] = max(item["longest"], len(block.pairs))
    partners: dict[str, set[str]] = defaultdict(set)
    for left, right in evidence:
        partners[left].add(right)
        partners[right].add(left)
    rows: list[dict[str, Any]] = []
    for (left, right), item in sorted(evidence.items()):
        accepted = len(partners[left]) == 1 and len(partners[right]) == 1
        rows.append(
            {
                "backbone_id": backbone_id,
                "edge_id": _digest_id("PPTB-E-", (left, right)),
                "gene_id_a": left,
                "gene_id_b": right,
                "seqid_a": genes[left].seqid,
                "seqid_b": genes[right].seqid,
                "support_block_count": len(item["blocks"]),
                "support_block_ids": ",".join(sorted(item["blocks"])),
                "longest_block_pairs": item["longest"],
                "status": "accepted" if accepted else "rejected",
                "reason": (
                    "base_only_cross_seqid_reciprocal_self_synteny"
                    if accepted
                    else "base_only_nonreciprocal_multiple_partners"
                ),
                "relationship_scope": "truth_blind_fixed_target_self_synteny",
            }
        )
    return rows


def _cell_rows(
    backbone_id: str,
    genes: dict[str, _BaseGene],
    blocks: list[_CanonicalBlock],
    accepted_pairs: set[tuple[str, str]],
) -> tuple[list[dict[str, Any]], int]:
    rows: list[dict[str, Any]] = []
    inconsistent = 0

    def emit(
        block: _CanonicalBlock,
        first_source: str,
        first_partner: str,
        second_source: str,
        second_partner: str,
    ) -> None:
        nonlocal inconsistent
        endpoints = sorted(
            ((first_source, first_partner), (second_source, second_partner)),
            key=lambda pair: (genes[pair[0]].gene_order, pair[0], pair[1]),
        )
        (source_left, partner_at_left), (source_right, partner_at_right) = endpoints
        partner_delta = (
            genes[partner_at_right].gene_order - genes[partner_at_left].gene_order
        )
        expected_positive = block.orientation == "plus"
        if partner_delta == 0 or (partner_delta > 0) != expected_positive:
            inconsistent += 1
            return
        semantic = {
            "block": block.canonical_block_id,
            "source_seqid": genes[source_left].seqid,
            "source_left": source_left,
            "source_right": source_right,
            "partner_seqid": genes[partner_at_left].seqid,
            "partner_at_source_left": partner_at_left,
            "partner_at_source_right": partner_at_right,
        }
        rows.append(
            {
                "backbone_id": backbone_id,
                "cell_id": _digest_id("PPTB-C-", semantic),
                "canonical_block_id": block.canonical_block_id,
                "source_seqid": genes[source_left].seqid,
                "source_left_gene_id": source_left,
                "source_left_order": genes[source_left].gene_order,
                "source_left_midpoint2": genes[source_left].midpoint2,
                "source_right_gene_id": source_right,
                "source_right_order": genes[source_right].gene_order,
                "source_right_midpoint2": genes[source_right].midpoint2,
                "partner_seqid": genes[partner_at_left].seqid,
                "partner_at_source_left_gene_id": partner_at_left,
                "partner_at_source_left_order": genes[partner_at_left].gene_order,
                "partner_at_source_left_midpoint2": genes[partner_at_left].midpoint2,
                "partner_at_source_right_gene_id": partner_at_right,
                "partner_at_source_right_order": genes[partner_at_right].gene_order,
                "partner_at_source_right_midpoint2": genes[partner_at_right].midpoint2,
                "orientation": block.orientation,
                "block_pair_count": len(block.pairs),
            }
        )

    for block in blocks:
        for first, second in zip(block.pairs, block.pairs[1:], strict=False):
            if (
                tuple(sorted(first)) not in accepted_pairs
                or tuple(sorted(second)) not in accepted_pairs
            ):
                continue
            emit(block, first[0], first[1], second[0], second[1])
            emit(block, first[1], first[0], second[1], second[0])
    return sorted(rows, key=lambda row: row["cell_id"]), inconsistent


def _merge_intervals(intervals: Iterable[tuple[int, int]]) -> tuple[tuple[int, int], ...]:
    merged: list[list[int]] = []
    for start, end in sorted(intervals):
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return tuple((start, end) for start, end in merged)


def _in_intervals(value: int, intervals: tuple[tuple[int, int], ...]) -> bool:
    starts = [start for start, _ in intervals]
    index = bisect_right(starts, value) - 1
    return index >= 0 and value <= intervals[index][1]


def _backbone_applicability_metrics(
    genes: dict[str, _BaseGene],
    cell_rows: list[dict[str, Any]],
    primary_seqids: Iterable[str],
) -> dict[str, Any]:
    intervals_by_source: dict[str, list[tuple[int, int]]] = defaultdict(list)
    intervals_by_source_partner: dict[tuple[str, str], list[tuple[int, int]]] = defaultdict(list)
    cell_counts: Counter[str] = Counter()
    for row in cell_rows:
        seqid = str(row["source_seqid"])
        partner = str(row["partner_seqid"])
        interval = (
            int(row["source_left_midpoint2"]),
            int(row["source_right_midpoint2"]),
        )
        intervals_by_source[seqid].append(interval)
        intervals_by_source_partner[(seqid, partner)].append(interval)
        cell_counts[seqid] += 1
    merged_source = {
        seqid: _merge_intervals(intervals)
        for seqid, intervals in intervals_by_source.items()
    }
    merged_partner = {
        key: _merge_intervals(intervals)
        for key, intervals in intervals_by_source_partner.items()
    }
    partners_by_source: dict[str, tuple[str, ...]] = defaultdict(tuple)
    for seqid, partner in sorted(merged_partner):
        partners_by_source[seqid] = (*partners_by_source[seqid], partner)

    covered_genes = 0
    uniquely_partnered_covered_genes = 0
    for gene in genes.values():
        if not _in_intervals(gene.midpoint2, merged_source.get(gene.seqid, ())):
            continue
        covered_genes += 1
        compatible_partners = {
            partner
            for partner in partners_by_source.get(gene.seqid, ())
            if _in_intervals(
                gene.midpoint2,
                merged_partner[(gene.seqid, partner)],
            )
        }
        uniquely_partnered_covered_genes += len(compatible_partners) == 1

    primary = tuple(primary_seqids)
    chromosomes_with_cells = [seqid for seqid in primary if cell_counts[seqid] > 0]
    return {
        "primary_gene_midpoint_backbone_coverage": covered_genes / len(genes),
        "primary_chromosome_cell_coverage_fraction": (
            len(chromosomes_with_cells) / len(primary)
        ),
        "minimum_cells_per_covered_chromosome_observed": (
            min(cell_counts[seqid] for seqid in chromosomes_with_cells)
            if chromosomes_with_cells
            else 0
        ),
        "unique_partner_chromosome_fraction_among_covered_genes": (
            uniquely_partnered_covered_genes / covered_genes if covered_genes else 0.0
        ),
        "covered_primary_genes": covered_genes,
        "uniquely_partnered_covered_primary_genes": uniquely_partnered_covered_genes,
        "primary_chromosomes_with_cells": len(chromosomes_with_cells),
        "directed_cells_by_primary_chromosome": {
            seqid: cell_counts[seqid] for seqid in sorted(primary)
        },
    }


def _write_tsv(path: Path, fieldnames: tuple[str, ...], rows: list[dict[str, Any]]) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
        handle.flush()
        os.fsync(handle.fileno())


def _write_json(path: Path, value: dict[str, Any]) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())


def build_fixed_target_backbone(
    *,
    base_gff_path: str | Path,
    query_wgdi_gff_path: str | Path,
    base_only_collinearity_path: str | Path,
    policy: FixedTargetBackbonePolicy,
    output_dir_path: str | Path,
) -> dict[str, Any]:
    """Build an immutable target self-synteny backbone without candidate input."""

    normalized_policy = policy.normalized()
    base_path = _required_regular_file(base_gff_path, label="base GFF3")
    query_path = _required_regular_file(query_wgdi_gff_path, label="WGDI query GFF")
    collinearity_path = _required_regular_file(
        base_only_collinearity_path, label="base-only WGDI collinearity"
    )
    output = Path(output_dir_path)
    working = Path(f"{output}.working")
    _reject_symlink(output, label="output directory")
    _reject_symlink(working, label="working directory")
    if output.exists() or working.exists():
        raise FileExistsError(f"Refusing to overwrite fixed target backbone: {output}")

    query_genes = _parse_query_gff(query_path, set(normalized_policy["primary_seqids"]))
    feature_id_by_wgdi_id, alias_manifest = _map_wgdi_ids_to_feature_ids(
        base_path, set(query_genes)
    )
    genes = _parse_base_genes(base_path, query_genes, feature_id_by_wgdi_id)
    raw_blocks = parse_wgdi_collinearity(collinearity_path, "base_only")
    blocks, block_counts, block_statistic_audit = _canonicalize_blocks(
        raw_blocks, genes, feature_id_by_wgdi_id, normalized_policy
    )
    backbone_id = _backbone_id(genes, blocks, normalized_policy)
    gene_rows = _gene_rows(backbone_id, genes)
    edge_rows = _edge_rows(backbone_id, genes, blocks)
    accepted_pairs = {
        tuple(sorted((str(row["gene_id_a"]), str(row["gene_id_b"]))))
        for row in edge_rows
        if row["status"] == "accepted"
    }
    cell_rows, inconsistent_cells = _cell_rows(
        backbone_id, genes, blocks, accepted_pairs
    )
    applicability_metrics = _backbone_applicability_metrics(
        genes, cell_rows, normalized_policy["primary_seqids"]
    )

    try:
        output.parent.mkdir(parents=True, exist_ok=True)
        _reject_symlink(output.parent, label="output parent")
        working.mkdir()
        genes_path = working / "genes.tsv"
        edges_path = working / "edges.tsv"
        cells_path = working / "cells.tsv"
        manifest_path = working / "manifest.json"
        _write_tsv(genes_path, tuple(gene_rows[0]), gene_rows)
        edge_fields = (
            "backbone_id", "edge_id", "gene_id_a", "gene_id_b", "seqid_a",
            "seqid_b", "support_block_count", "support_block_ids",
            "longest_block_pairs", "status", "reason", "relationship_scope",
        )
        cell_fields = (
            "backbone_id", "cell_id", "canonical_block_id", "source_seqid",
            "source_left_gene_id", "source_left_order", "source_left_midpoint2",
            "source_right_gene_id", "source_right_order", "source_right_midpoint2",
            "partner_seqid", "partner_at_source_left_gene_id",
            "partner_at_source_left_order", "partner_at_source_left_midpoint2",
            "partner_at_source_right_gene_id", "partner_at_source_right_order",
            "partner_at_source_right_midpoint2", "orientation", "block_pair_count",
        )
        _write_tsv(
            edges_path,
            edge_fields,
            edge_rows,
        )
        _write_tsv(cells_path, cell_fields, cell_rows)
        manifest: dict[str, Any] = {
            "schema_version": FIXED_TARGET_BACKBONE_SCHEMA_VERSION,
            "backbone_id": backbone_id,
            "generator": {"name": "PloidyPatch", "version": __version__},
            "truth_access": False,
            "candidate_access": False,
            "semantic_invariants": dict(FIXED_TARGET_BACKBONE_SEMANTIC_INVARIANTS),
            "applicable": bool(blocks and cell_rows),
            "applicability_reason": (
                "unambiguous_base_only_backbone_available"
                if blocks and cell_rows
                else "no_unambiguous_orientation_consistent_backbone_cells"
            ),
            "inputs": {
                "base_gff": {"file_name": base_path.name, "sha256": _file_sha256(base_path)},
                "query_wgdi_gff": {
                    "file_name": query_path.name,
                    "sha256": _file_sha256(query_path),
                },
                "base_only_collinearity": {
                    "file_name": collinearity_path.name,
                    "sha256": _file_sha256(collinearity_path),
                },
            },
            "policy": normalized_policy,
            "schemas": {
                "genes": FIXED_TARGET_BACKBONE_GENE_SCHEMA_VERSION,
                "edges": FIXED_TARGET_BACKBONE_EDGE_SCHEMA_VERSION,
                "cells": FIXED_TARGET_BACKBONE_CELL_SCHEMA_VERSION,
            },
            "counts": {
                "genes": len(gene_rows),
                "edges": len(edge_rows),
                "accepted_edges": sum(row["status"] == "accepted" for row in edge_rows),
                "directed_cells": len(cell_rows),
                "orientation_inconsistent_cells_skipped": inconsistent_cells,
                **block_counts,
            },
            "identifier_mapping": alias_manifest,
            "block_statistic_audit": block_statistic_audit,
            "applicability_metrics": applicability_metrics,
            "outputs": {
                "genes": {"file_name": genes_path.name, "sha256": _file_sha256(genes_path), "rows": len(gene_rows)},
                "edges": {"file_name": edges_path.name, "sha256": _file_sha256(edges_path), "rows": len(edge_rows)},
                "cells": {"file_name": cells_path.name, "sha256": _file_sha256(cells_path), "rows": len(cell_rows)},
            },
            "claim_boundary": (
                "truth-blind fixed target self-synteny backbone; not candidate evidence, "
                "gene-tree proof, or event-confirmed homeology"
            ),
        }
        _write_json(manifest_path, manifest)
        write_sha256sums(working)
        verify_sha256sums(working, ignore_checksum_file=True)
        os.replace(working, output)
        return manifest
    except Exception:
        if working.exists():
            shutil.rmtree(working)
        raise
