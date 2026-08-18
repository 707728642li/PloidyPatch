from __future__ import annotations

import csv
import hashlib
import json
import re
from collections import defaultdict
from pathlib import Path
from typing import Any

from . import __version__
from .gff import parse_attributes
from .io import open_text
from .wgdi_summary import _file_sha256, _read_query_gff, parse_wgdi_collinearity


SELF_WGD_PAIR_SCHEMA_VERSION = "ploidypatch.self_wgd_pairs.v2"
SELF_WGD_PAIR_MANIFEST_SCHEMA_VERSION = "ploidypatch.self_wgd_pair_manifest.v2"
WGD_EVENT_LABEL = re.compile(r"^[A-Za-z][A-Za-z0-9_.-]*$")


def _pair_id(wgd_event: str, gene_a: str, gene_b: str) -> str:
    digest = hashlib.sha256(
        f"{wgd_event}\0{gene_a}\0{gene_b}".encode("utf-8")
    ).hexdigest()[:20]
    return f"PPSW-{digest}"


def _map_wgdi_to_feature_ids(
    source_gff_path: str | Path | None,
    wgdi_gene_ids: set[str],
) -> tuple[dict[str, str], dict[str, Any] | None]:
    if source_gff_path is None:
        return {gene_id: gene_id for gene_id in wgdi_gene_ids}, None

    aliases: dict[str, set[str]] = defaultdict(set)
    source_path = Path(source_gff_path)
    gene_records = 0
    with open_text(source_path) as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\r\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) != 9:
                raise ValueError(
                    f"Malformed source GFF3 line {line_number}: expected 9 fields"
                )
            if fields[2] != "gene":
                continue
            attributes, malformed = parse_attributes(fields[8])
            if malformed:
                raise ValueError(
                    f"Malformed source GFF3 attributes at line {line_number}"
                )
            feature_id = attributes.get("ID", "")
            if not feature_id:
                continue
            gene_records += 1
            for key in ("ID", "gene_id", "locus_tag"):
                alias = attributes.get(key, "")
                if alias:
                    aliases[alias].add(feature_id)

    mapping: dict[str, str] = {}
    ambiguous: dict[str, set[str]] = {}
    for wgdi_gene_id in wgdi_gene_ids:
        feature_ids = aliases.get(wgdi_gene_id, set())
        if len(feature_ids) == 1:
            mapping[wgdi_gene_id] = next(iter(feature_ids))
        elif len(feature_ids) > 1:
            ambiguous[wgdi_gene_id] = feature_ids
    missing = wgdi_gene_ids - set(mapping) - set(ambiguous)
    if missing or ambiguous:
        examples = {
            "missing": sorted(missing)[:10],
            "ambiguous": {
                key: sorted(ambiguous[key])
                for key in sorted(ambiguous)[:10]
            },
        }
        raise ValueError(
            "WGDI gene IDs do not map uniquely to source gene feature IDs: "
            + json.dumps(examples, sort_keys=True)
        )
    return mapping, {
        "file_name": source_path.name,
        "sha256": _file_sha256(source_path),
        "gene_records": gene_records,
        "mapped_wgdi_genes": len(mapping),
        "alias_keys": ["ID", "gene_id", "locus_tag"],
        "mapping_policy": "unique_exact_attribute_value",
    }


def infer_self_wgdi_pairs(
    *,
    query_wgdi_gff_path: str | Path,
    collinearity_path: str | Path,
    source_gff_path: str | Path | None = None,
    output_pair_tsv_path: str | Path,
    decisions_tsv_path: str | Path,
    wgd_event: str,
    min_block_pairs: int = 20,
    require_different_seqids: bool = True,
    require_reciprocal_unique: bool = True,
) -> dict[str, Any]:
    """Infer conservative paleopolyploid duplicate pairs from self-WGDI.

    The output is syntenic pair evidence tied to a named candidate WGD event;
    it is not a gene-tree proof of homeology.
    """

    if not WGD_EVENT_LABEL.fullmatch(wgd_event):
        raise ValueError(f"Unsafe WGD event label: {wgd_event}")
    if min_block_pairs < 2:
        raise ValueError("min_block_pairs must be at least two")
    output_pair = Path(output_pair_tsv_path)
    decisions = Path(decisions_tsv_path)
    manifest_path = Path(str(output_pair) + ".manifest.json")
    existing = [path for path in (output_pair, decisions, manifest_path) if path.exists()]
    if existing:
        raise FileExistsError(
            "Refusing to overwrite self-WGD pair artifact(s): "
            + ", ".join(str(path) for path in existing)
        )

    genes = _read_query_gff(query_wgdi_gff_path)
    blocks = parse_wgdi_collinearity(collinearity_path, "self_wgdi")
    pair_blocks: dict[tuple[str, str], set[str]] = defaultdict(set)
    pair_best_score: dict[tuple[str, str], float] = {}
    pair_best_pvalue: dict[tuple[str, str], float] = {}
    pair_longest_block: dict[tuple[str, str], int] = defaultdict(int)
    skipped_counts: dict[str, int] = defaultdict(int)
    for block in blocks:
        if len(block.pairs) < min_block_pairs:
            skipped_counts["pairs_in_short_blocks"] += len(block.pairs)
            continue
        for query_gene, target_gene in block.pairs:
            if query_gene not in genes or target_gene not in genes:
                raise ValueError(
                    "Self-WGDI block references a gene absent from WGDI GFF: "
                    f"{query_gene}, {target_gene}"
                )
            if query_gene == target_gene:
                skipped_counts["self_pairs"] += 1
                continue
            gene_a, gene_b = sorted((query_gene, target_gene))
            if require_different_seqids and genes[gene_a][0] == genes[gene_b][0]:
                skipped_counts["same_seqid_pairs"] += 1
                continue
            pair = (gene_a, gene_b)
            pair_blocks[pair].add(block.block_id)
            pair_best_score[pair] = max(
                block.score, pair_best_score.get(pair, block.score)
            )
            pair_best_pvalue[pair] = min(
                block.pvalue, pair_best_pvalue.get(pair, block.pvalue)
            )
            pair_longest_block[pair] = max(
                pair_longest_block[pair], len(block.pairs)
            )

    partners: dict[str, set[str]] = defaultdict(set)
    for gene_a, gene_b in pair_blocks:
        partners[gene_a].add(gene_b)
        partners[gene_b].add(gene_a)

    feature_id_by_wgdi_id, source_gff_manifest = _map_wgdi_to_feature_ids(
        source_gff_path, set(genes)
    )

    fields = (
        "pair_id",
        "gene_id_a",
        "gene_id_b",
        "wgdi_gene_id_a",
        "wgdi_gene_id_b",
        "wgd_event",
        "seqid_a",
        "seqid_b",
        "support_block_count",
        "support_block_ids",
        "best_block_score",
        "best_block_pvalue",
        "longest_block_pairs",
        "relationship_scope",
    )
    decision_fields = (*fields, "status", "reason")
    accepted_rows: list[dict[str, str]] = []
    decision_rows: list[dict[str, str]] = []
    for wgdi_gene_a, wgdi_gene_b in sorted(pair_blocks):
        gene_a = feature_id_by_wgdi_id[wgdi_gene_a]
        gene_b = feature_id_by_wgdi_id[wgdi_gene_b]
        row = {
            "pair_id": _pair_id(wgd_event, gene_a, gene_b),
            "gene_id_a": gene_a,
            "gene_id_b": gene_b,
            "wgdi_gene_id_a": wgdi_gene_a,
            "wgdi_gene_id_b": wgdi_gene_b,
            "wgd_event": wgd_event,
            "seqid_a": genes[wgdi_gene_a][0],
            "seqid_b": genes[wgdi_gene_b][0],
            "support_block_count": str(
                len(pair_blocks[(wgdi_gene_a, wgdi_gene_b)])
            ),
            "support_block_ids": ",".join(
                sorted(pair_blocks[(wgdi_gene_a, wgdi_gene_b)])
            ),
            "best_block_score": str(
                pair_best_score[(wgdi_gene_a, wgdi_gene_b)]
            ),
            "best_block_pvalue": str(
                pair_best_pvalue[(wgdi_gene_a, wgdi_gene_b)]
            ),
            "longest_block_pairs": str(
                pair_longest_block[(wgdi_gene_a, wgdi_gene_b)]
            ),
            "relationship_scope": "candidate_named_wgd_syntenic_duplicate",
        }
        if require_reciprocal_unique and (
            len(partners[wgdi_gene_a]) != 1
            or len(partners[wgdi_gene_b]) != 1
        ):
            status, reason = "rejected", "nonreciprocal_multiple_partners"
        else:
            status, reason = "accepted", "cross_seqid_reciprocal_self_synteny"
            accepted_rows.append(row)
        decision_rows.append({**row, "status": status, "reason": reason})

    output_pair.parent.mkdir(parents=True, exist_ok=True)
    decisions.parent.mkdir(parents=True, exist_ok=True)
    with output_pair.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(accepted_rows)
    with decisions.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=decision_fields,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(decision_rows)

    reason_counts: dict[str, int] = defaultdict(int)
    for row in decision_rows:
        reason_counts[row["reason"]] += 1
    manifest: dict[str, Any] = {
        "schema_version": SELF_WGD_PAIR_MANIFEST_SCHEMA_VERSION,
        "pair_schema_version": SELF_WGD_PAIR_SCHEMA_VERSION,
        "generator": {"name": "PloidyPatch", "version": __version__},
        "inputs": {
            "query_wgdi_gff": {
                "file_name": Path(query_wgdi_gff_path).name,
                "sha256": _file_sha256(query_wgdi_gff_path),
                "genes": len(genes),
            },
            "collinearity": {
                "file_name": Path(collinearity_path).name,
                "sha256": _file_sha256(collinearity_path),
                "blocks": len(blocks),
                "pair_rows": sum(len(block.pairs) for block in blocks),
            },
            "source_gff": source_gff_manifest,
        },
        "parameters": {
            "wgd_event": wgd_event,
            "min_block_pairs": min_block_pairs,
            "require_different_seqids": require_different_seqids,
            "require_reciprocal_unique": require_reciprocal_unique,
            "relationship_scope": "candidate_named_wgd_syntenic_duplicate",
            "gene_identifier_policy": (
                "unique_source_gff_alias_mapping"
                if source_gff_path is not None
                else "wgdi_id_is_feature_id"
            ),
        },
        "counts": {
            "candidate_pairs": len(pair_blocks),
            "accepted_pairs": len(accepted_rows),
            "accepted_genes": len(
                {
                    gene
                    for row in accepted_rows
                    for gene in (row["gene_id_a"], row["gene_id_b"])
                }
            ),
            "decision_reason_counts": dict(sorted(reason_counts.items())),
            "skipped_pair_rows": dict(sorted(skipped_counts.items())),
        },
        "outputs": {
            "pairs": {
                "file_name": output_pair.name,
                "sha256": _file_sha256(output_pair),
                "rows": len(accepted_rows),
            },
            "decisions": {
                "file_name": decisions.name,
                "sha256": _file_sha256(decisions),
                "rows": len(decision_rows),
            },
        },
        "claim_boundary": (
            "syntenic duplicate pair tied to a candidate named WGD event; "
            "not gene-tree-confirmed homeology"
        ),
    }
    with manifest_path.open("x", encoding="utf-8", newline="") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return manifest
