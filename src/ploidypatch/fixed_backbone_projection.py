from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import shutil
from bisect import bisect_left, bisect_right
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

from . import __version__
from .artifact_manifest import sha256_file, verify_sha256sums, write_sha256sums


SCHEMA_VERSION = "ploidypatch.fixed_backbone_projection.v0.6"
CATALOG_FIELDS = (
    "candidate_digest",
    "candidate_gene_id",
    "query_protein_id",
    "seqid",
    "start",
    "end",
    "strand",
)
HIT_FIELDS = (
    "query_protein_id",
    "subject_gene_id",
    "pident",
    "alignment_length",
    "query_length",
    "evalue",
    "bitscore",
)
OUTPUT_FIELDS = (
    "gene_id",
    "consensus_digest",
    "pair_id",
    "partner_gene_id",
    "support_block_count",
    "longest_block_pairs",
    "status",
    "reason",
    "evidence_type",
    "source_base_gene_id",
    "anchor_cell_count",
    "compatible_partner_count",
    "qualifying_base_hit_count",
    "backbone_id",
)


@dataclass(frozen=True)
class FixedBackboneProjectionPolicy:
    minimum_query_coverage: float = 0.50
    maximum_evalue: float = 1e-5

    def normalized(self) -> dict[str, Any]:
        if not 0 < self.minimum_query_coverage <= 1:
            raise ValueError("minimum_query_coverage must be in (0, 1]")
        if not math.isfinite(self.maximum_evalue) or self.maximum_evalue <= 0:
            raise ValueError("maximum_evalue must be finite and positive")
        return {
            "minimum_query_coverage": self.minimum_query_coverage,
            "maximum_evalue": self.maximum_evalue,
            "candidate_candidate_evidence": "forbidden",
            "conflict_peer_inheritance": "forbidden",
            "gap_extrapolation": "forbidden",
            "multiple_compatible_partners": "abstain",
            "direct_overlap_failure_fallback": "forbidden",
        }


@dataclass(frozen=True)
class _Gene:
    gene_id: str
    seqid: str
    coding_start: int
    coding_end: int
    midpoint2: int


@dataclass(frozen=True)
class _Candidate:
    digest: str
    gene_id: str
    protein_id: str
    seqid: str
    start: int
    end: int
    strand: str

    @property
    def midpoint2(self) -> int:
        return self.start + self.end


def _read_exact_tsv(path: Path, fields: tuple[str, ...]) -> list[dict[str, str]]:
    if not path.is_file() or path.is_symlink() or path.stat().st_size == 0:
        raise ValueError(f"Missing, empty or symlinked TSV: {path}")
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != fields:
            raise ValueError(f"TSV columns differ for {path}: expected {fields}")
        return list(reader)


def _load_backbone(root: Path) -> tuple[str, dict[str, _Gene], list[dict[str, str]], list[dict[str, str]], dict[str, Any]]:
    verify_sha256sums(root, ignore_checksum_file=True)
    manifest_path = root / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("schema_version") != "ploidypatch.fixed_target_backbone.v0.6":
        raise ValueError("Fixed target backbone schema differs")
    if manifest.get("candidate_access") is not False or manifest.get("truth_access") is not False:
        raise ValueError("Fixed target backbone access declaration differs")
    backbone_id = manifest.get("backbone_id")
    if not isinstance(backbone_id, str) or not backbone_id.startswith("PPTB-"):
        raise ValueError("Fixed target backbone ID is invalid")
    gene_fields = (
        "backbone_id", "gene_id", "wgdi_gene_id", "seqid", "start", "end",
        "strand", "gene_order", "midpoint2", "coding_start", "coding_end",
    )
    genes: dict[str, _Gene] = {}
    for row in _read_exact_tsv(root / "genes.tsv", gene_fields):
        gene_id = row["gene_id"]
        if row["backbone_id"] != backbone_id or not gene_id or gene_id in genes:
            raise ValueError("Backbone genes contain an invalid ID or binding")
        genes[gene_id] = _Gene(
            gene_id=gene_id,
            seqid=row["seqid"],
            coding_start=int(row["coding_start"]),
            coding_end=int(row["coding_end"]),
            midpoint2=int(row["midpoint2"]),
        )
    edges = _read_exact_tsv(
        root / "edges.tsv",
        (
            "backbone_id", "edge_id", "gene_id_a", "gene_id_b", "seqid_a",
            "seqid_b", "support_block_count", "support_block_ids",
            "longest_block_pairs", "status", "reason", "relationship_scope",
        ),
    )
    cells = _read_exact_tsv(
        root / "cells.tsv",
        (
            "backbone_id", "cell_id", "canonical_block_id", "source_seqid",
            "source_left_gene_id", "source_left_order", "source_left_midpoint2",
            "source_right_gene_id", "source_right_order", "source_right_midpoint2",
            "partner_seqid", "partner_at_source_left_gene_id",
            "partner_at_source_left_order", "partner_at_source_left_midpoint2",
            "partner_at_source_right_gene_id", "partner_at_source_right_order",
            "partner_at_source_right_midpoint2", "orientation", "block_pair_count",
        ),
    )
    if any(row["backbone_id"] != backbone_id for row in (*edges, *cells)):
        raise ValueError("Backbone table ID binding differs")
    return backbone_id, genes, edges, cells, manifest


def _parse_candidates(path: Path, known_seqids: set[str]) -> list[_Candidate]:
    rows = _read_exact_tsv(path, CATALOG_FIELDS)
    candidates: list[_Candidate] = []
    digests: set[str] = set()
    genes: set[str] = set()
    proteins: set[str] = set()
    for line_number, row in enumerate(rows, start=2):
        try:
            start, end = int(row["start"]), int(row["end"])
        except ValueError as exc:
            raise ValueError(f"Candidate coordinates are non-integer at line {line_number}") from exc
        digest, gene_id, protein_id = (
            row["candidate_digest"], row["candidate_gene_id"], row["query_protein_id"]
        )
        if (
            len(digest) != 64
            or any(character not in "0123456789abcdef" for character in digest)
            or not gene_id
            or not protein_id
            or digest in digests
            or gene_id in genes
            or protein_id in proteins
            or row["seqid"] not in known_seqids
            or start < 1
            or end < start
            or row["strand"] not in {"+", "-"}
        ):
            raise ValueError(f"Invalid or duplicate candidate at line {line_number}")
        digests.add(digest)
        genes.add(gene_id)
        proteins.add(protein_id)
        candidates.append(
            _Candidate(digest, gene_id, protein_id, row["seqid"], start, end, row["strand"])
        )
    if not candidates:
        raise ValueError("Candidate catalog is empty")
    return candidates


def _parse_hits(
    path: Path,
    candidates: Iterable[_Candidate],
    base_gene_ids: set[str],
    policy: dict[str, Any],
) -> tuple[dict[str, set[str]], int]:
    protein_ids = {candidate.protein_id for candidate in candidates}
    qualifying: dict[str, set[str]] = defaultdict(set)
    seen_rows: set[tuple[str, str, str, str, str, str, str]] = set()
    rows = _read_exact_tsv(path, HIT_FIELDS)
    for line_number, row in enumerate(rows, start=2):
        key = tuple(row[field] for field in HIT_FIELDS)
        if key in seen_rows:
            raise ValueError(f"Duplicate fixed-base hit row at line {line_number}")
        seen_rows.add(key)
        query, subject = row["query_protein_id"], row["subject_gene_id"]
        if query not in protein_ids or subject not in base_gene_ids:
            raise ValueError(f"Hit references an unknown candidate/base ID at line {line_number}")
        try:
            alignment_length = int(row["alignment_length"])
            query_length = int(row["query_length"])
            pident = float(row["pident"])
            evalue = float(row["evalue"])
            bitscore = float(row["bitscore"])
        except ValueError as exc:
            raise ValueError(f"Malformed fixed-base hit at line {line_number}") from exc
        if (
            alignment_length < 1
            or query_length < alignment_length
            or not 0 <= pident <= 100
            or not math.isfinite(evalue)
            or evalue < 0
            or not math.isfinite(bitscore)
            or bitscore < 0
        ):
            raise ValueError(f"Invalid fixed-base hit values at line {line_number}")
        if (
            alignment_length / query_length >= policy["minimum_query_coverage"]
            and evalue <= policy["maximum_evalue"]
        ):
            qualifying[query].add(subject)
    return qualifying, len(rows)


def _load_hits_manifest(
    path: Path,
    *,
    backbone_id: str,
    catalog_path: Path,
    hits_path: Path,
    base_gene_count: int,
    policy: dict[str, Any],
) -> dict[str, Any]:
    if not path.is_file() or path.is_symlink() or path.stat().st_size == 0:
        raise ValueError(f"Missing, empty or symlinked fixed-base hit manifest: {path}")
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict) or value.get("schema_version") != (
        "ploidypatch.fixed_base_protein_hits.v0.6"
    ):
        raise ValueError("Fixed-base hit manifest schema differs")
    if (
        value.get("truth_access") is not False
        or value.get("candidate_candidate_hits") is not False
        or value.get("backbone_id") != backbone_id
        or value.get("database_scope")
        != "fixed_backbone_representative_proteins_only"
        or value.get("exhaustive_target_search") is not True
        or value.get("candidate_catalog_sha256") != sha256_file(catalog_path)
        or value.get("base_gene_count") != base_gene_count
    ):
        raise ValueError("Fixed-base hit manifest role or input binding differs")
    output = value.get("output")
    parameters = value.get("parameters")
    if (
        not isinstance(output, dict)
        or output.get("hits_sha256") != sha256_file(hits_path)
        or not isinstance(output.get("rows"), int)
        or output["rows"] < 0
        or not isinstance(parameters, dict)
        or parameters.get("maximum_target_sequences") != 0
        or parameters.get("more_sensitive") is not True
        or parameters.get("outfmt_fields") != list(HIT_FIELDS)
    ):
        raise ValueError("Fixed-base hit manifest output or exhaustive-search parameters differ")
    try:
        upstream_evalue = float(parameters["maximum_evalue"])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("Fixed-base hit manifest maximum_evalue is invalid") from exc
    if (
        not math.isfinite(upstream_evalue)
        or upstream_evalue < policy["maximum_evalue"]
    ):
        raise ValueError("Upstream hit search is narrower than the projection policy")
    for key in ("query_protein_fasta_sha256", "base_protein_fasta_sha256"):
        digest = value.get(key)
        if not isinstance(digest, str) or len(digest) != 64 or any(
            character not in "0123456789abcdef" for character in digest
        ):
            raise ValueError(f"Fixed-base hit manifest lacks a valid {key}")
    return value


def _pair_id(backbone_id: str, digest: str, partner: str) -> str:
    value = f"{backbone_id}\0{digest}\0{partner}".encode()
    return "PPTP-" + hashlib.sha256(value).hexdigest()[:24]


def _interval_index(
    values: Iterable[Any], start: Any, end: Any
) -> tuple[list[Any], list[int], list[int]]:
    ordered = sorted(values, key=lambda value: (start(value), end(value), str(value)))
    starts: list[int] = []
    prefix_max_end: list[int] = []
    running = -1
    for value in ordered:
        starts.append(start(value))
        running = max(running, end(value))
        prefix_max_end.append(running)
    return ordered, starts, prefix_max_end


def _overlapping_intervals(
    index: tuple[list[Any], list[int], list[int]],
    start: int,
    end: int,
    value_end: Any,
) -> list[Any]:
    values, starts, prefix_max_end = index
    position = bisect_right(starts, end) - 1
    output: list[Any] = []
    while position >= 0 and prefix_max_end[position] >= start:
        value = values[position]
        if value_end(value) >= start:
            output.append(value)
        position -= 1
    return output


def _strict_point_intervals(
    index: tuple[list[Any], list[int], list[int]], point: int, value_end: Any
) -> list[Any]:
    values, starts, prefix_max_end = index
    position = bisect_left(starts, point) - 1
    output: list[Any] = []
    while position >= 0 and prefix_max_end[position] > point:
        value = values[position]
        if value_end(value) > point:
            output.append(value)
        position -= 1
    return output


def project_candidates_to_fixed_backbone(
    *,
    backbone_dir_path: str | Path,
    candidate_catalog_tsv_path: str | Path,
    fixed_base_hits_tsv_path: str | Path,
    fixed_base_hits_manifest_path: str | Path,
    policy: FixedBackboneProjectionPolicy,
    output_dir_path: str | Path,
) -> dict[str, Any]:
    """Project each candidate independently onto a sealed target-only backbone."""

    backbone_root = Path(backbone_dir_path)
    catalog_path = Path(candidate_catalog_tsv_path)
    hits_path = Path(fixed_base_hits_tsv_path)
    hits_manifest_path = Path(fixed_base_hits_manifest_path)
    output = Path(output_dir_path)
    working = Path(f"{output}.working")
    if output.exists() or working.exists():
        raise FileExistsError("Refusing to overwrite fixed-backbone projection")
    normalized_policy = policy.normalized()
    backbone_id, genes, edge_rows, cell_rows, backbone_manifest = _load_backbone(backbone_root)
    candidates = _parse_candidates(catalog_path, {gene.seqid for gene in genes.values()})
    hits_manifest = _load_hits_manifest(
        hits_manifest_path,
        backbone_id=backbone_id,
        catalog_path=catalog_path,
        hits_path=hits_path,
        base_gene_count=len(genes),
        policy=normalized_policy,
    )
    qualifying_hits, hit_row_count = _parse_hits(
        hits_path, candidates, set(genes), normalized_policy
    )
    if hits_manifest["output"]["rows"] != hit_row_count:
        raise ValueError("Fixed-base hit manifest row count differs")

    accepted_partner: dict[str, tuple[str, dict[str, str]]] = {}
    for row in edge_rows:
        if row["status"] != "accepted":
            continue
        for source, partner in (
            (row["gene_id_a"], row["gene_id_b"]),
            (row["gene_id_b"], row["gene_id_a"]),
        ):
            if source in accepted_partner and accepted_partner[source][0] != partner:
                raise ValueError("Backbone accepted edges are not reciprocal unique")
            accepted_partner[source] = (partner, row)
    base_by_seqid: dict[str, list[_Gene]] = defaultdict(list)
    for gene in genes.values():
        base_by_seqid[gene.seqid].append(gene)
    cells_by_seqid: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in cell_rows:
        cells_by_seqid[row["source_seqid"]].append(row)
    base_indexes = {
        seqid: _interval_index(
            values,
            lambda gene: gene.coding_start,
            lambda gene: gene.coding_end,
        )
        for seqid, values in base_by_seqid.items()
    }
    cell_indexes = {
        seqid: _interval_index(
            values,
            lambda row: int(row["source_left_midpoint2"]),
            lambda row: int(row["source_right_midpoint2"]),
        )
        for seqid, values in cells_by_seqid.items()
    }

    output_rows: list[dict[str, Any]] = []
    reasons: Counter[str] = Counter()
    for candidate in sorted(candidates, key=lambda item: item.digest):
        hits = qualifying_hits.get(candidate.protein_id, set())
        overlaps = _overlapping_intervals(
            base_indexes[candidate.seqid],
            candidate.start,
            candidate.end,
            lambda gene: gene.coding_end,
        )
        partner = ""
        support_blocks = 0
        longest_block = 0
        anchor_count = 0
        compatible: set[str] = set()
        source_gene = ""
        evidence_type = ""
        if len(overlaps) > 1:
            reason = "multiple_overlapping_base_coding_loci"
        elif len(overlaps) == 1:
            evidence_type = "direct_overlap_fixed_edge"
            source_gene = overlaps[0].gene_id
            edge = accepted_partner.get(source_gene)
            if edge is None:
                reason = "overlapping_base_locus_lacks_unique_fixed_edge"
            else:
                expected, edge_row = edge
                compatible = {expected} & hits
                if compatible == {expected}:
                    partner = expected
                    support_blocks = int(edge_row["support_block_count"])
                    longest_block = int(edge_row["longest_block_pairs"])
                    reason = "accepted_direct_overlap_unique_fixed_partner"
                else:
                    reason = "direct_fixed_partner_lacks_qualifying_base_hit"
        else:
            evidence_type = "strict_gap_anchor_projection"
            active_cells = _strict_point_intervals(
                cell_indexes.get(candidate.seqid, ([], [], [])),
                candidate.midpoint2,
                lambda row: int(row["source_right_midpoint2"]),
            )
            anchor_count = len(active_cells)
            for row in active_cells:
                low, high = sorted(
                    (
                        int(row["partner_at_source_left_midpoint2"]),
                        int(row["partner_at_source_right_midpoint2"]),
                    )
                )
                for hit_gene_id in hits:
                    hit_gene = genes[hit_gene_id]
                    if (
                        hit_gene.seqid == row["partner_seqid"]
                        and low < hit_gene.midpoint2 < high
                    ):
                        compatible.add(hit_gene_id)
            if not active_cells:
                reason = "no_strict_fixed_anchor_bracket"
            elif len(compatible) == 0:
                reason = "anchor_bracket_has_no_qualifying_base_hit"
            elif len(compatible) > 1:
                reason = "multiple_compatible_fixed_partners"
            else:
                partner = next(iter(compatible))
                supporting_cells = []
                for row in active_cells:
                    low, high = sorted(
                        (
                            int(row["partner_at_source_left_midpoint2"]),
                            int(row["partner_at_source_right_midpoint2"]),
                        )
                    )
                    gene = genes[partner]
                    if gene.seqid == row["partner_seqid"] and low < gene.midpoint2 < high:
                        supporting_cells.append(row)
                support_blocks = len({row["canonical_block_id"] for row in supporting_cells})
                longest_block = max(int(row["block_pair_count"]) for row in supporting_cells)
                reason = "accepted_strict_gap_unique_fixed_partner"
        status = "accepted" if partner else "rejected"
        reasons[reason] += 1
        output_rows.append(
            {
                "gene_id": candidate.gene_id,
                "consensus_digest": candidate.digest,
                "pair_id": _pair_id(backbone_id, candidate.digest, partner) if partner else "",
                "partner_gene_id": partner,
                "support_block_count": support_blocks if partner else "",
                "longest_block_pairs": longest_block if partner else "",
                "status": status,
                "reason": reason,
                "evidence_type": evidence_type,
                "source_base_gene_id": source_gene,
                "anchor_cell_count": anchor_count,
                "compatible_partner_count": len(compatible),
                "qualifying_base_hit_count": len(hits),
                "backbone_id": backbone_id,
            }
        )

    try:
        output.parent.mkdir(parents=True, exist_ok=True)
        working.mkdir()
        selection_path = working / "selection.tsv"
        with selection_path.open("x", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=OUTPUT_FIELDS, delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            writer.writerows(output_rows)
            handle.flush()
            os.fsync(handle.fileno())
        manifest: dict[str, Any] = {
            "schema_version": SCHEMA_VERSION,
            "generator": {"name": "PloidyPatch", "version": __version__},
            "backbone_id": backbone_id,
            "truth_access": False,
            "candidate_universe_invariant_by_construction": True,
            "policy": normalized_policy,
            "inputs": {
                "backbone_manifest_sha256": sha256_file(backbone_root / "manifest.json"),
                "candidate_catalog_sha256": sha256_file(catalog_path),
                "fixed_base_hits_sha256": sha256_file(hits_path),
                "fixed_base_hits_manifest_sha256": sha256_file(hits_manifest_path),
            },
            "upstream_backbone_claim_boundary": backbone_manifest["claim_boundary"],
            "counts": {
                "candidates": len(output_rows),
                "fixed_base_hit_rows": hit_row_count,
                "accepted": sum(row["status"] == "accepted" for row in output_rows),
                "rejected": sum(row["status"] == "rejected" for row in output_rows),
                "reasons": dict(sorted(reasons.items())),
            },
            "output": {
                "file_name": "selection.tsv",
                "sha256": sha256_file(selection_path),
                "rows": len(output_rows),
            },
            "claim_boundary": (
                "candidate-independent projection onto truth-blind fixed target self-synteny; "
                "not event-confirmed homeology and never automatic annotation approval"
            ),
        }
        manifest_path = working / "manifest.json"
        manifest_path.write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
            newline="",
        )
        write_sha256sums(working)
        verify_sha256sums(working, ignore_checksum_file=True)
        os.replace(working, output)
        return manifest
    except Exception:
        if working.exists():
            shutil.rmtree(working)
        raise
