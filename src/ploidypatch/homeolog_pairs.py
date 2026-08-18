from __future__ import annotations

import csv
import hashlib
import json
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable

from . import __version__
from .self_wgd_pairs import _map_wgdi_to_feature_ids
from .wgdi_summary import SOURCE_LABEL, _file_sha256, parse_wgdi_collinearity


HOMEOLOG_PAIR_SCHEMA_VERSION = "ploidypatch.homeolog_pairs.v1"
HOMEOLOG_PAIR_MANIFEST_SCHEMA_VERSION = "ploidypatch.homeolog_pair_manifest.v1"
OUTGROUP_DUPLICATE_PAIR_SCHEMA_VERSION = (
    "ploidypatch.outgroup_duplicated_pairs.v2"
)
OUTGROUP_DUPLICATE_MANIFEST_SCHEMA_VERSION = (
    "ploidypatch.outgroup_duplicated_pair_manifest.v2"
)
WGD_EVENT_LABEL = re.compile(r"^[A-Za-z][A-Za-z0-9_.-]*$")


def _parse_collinearity_inputs(
    values: Iterable[str],
) -> tuple[tuple[str, Path], ...]:
    parsed: list[tuple[str, Path]] = []
    seen: set[str] = set()
    for value in values:
        source, separator, raw_path = value.partition("=")
        if not separator or not source or not raw_path:
            raise ValueError("Each --collinearity value must be SOURCE=PATH")
        if not SOURCE_LABEL.fullmatch(source):
            raise ValueError(f"Unsafe WGDI source label: {source}")
        if source in seen:
            raise ValueError(f"Duplicate WGDI source label: {source}")
        seen.add(source)
        parsed.append((source, Path(raw_path)))
    if not parsed:
        raise ValueError("At least one WGDI collinearity input is required")
    return tuple(parsed)


def _read_source_groups(
    path: str | Path | None, sources: tuple[str, ...]
) -> tuple[dict[str, str], str, dict[str, Any] | None]:
    if path is None:
        return {source: source for source in sources}, "source", None
    mapping: dict[str, str] = {}
    group_path = Path(path)
    with group_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or not {
            "source",
            "support_group",
        } <= set(reader.fieldnames):
            raise ValueError(
                "Source-group map requires source and support_group columns"
            )
        for line_number, row in enumerate(reader, start=2):
            source = row["source"]
            group = row["support_group"]
            if not source or not group or source in mapping:
                raise ValueError(
                    f"Empty or duplicate source-group mapping at line {line_number}"
                )
            mapping[source] = group
    missing = set(sources) - set(mapping)
    extra = set(mapping) - set(sources)
    if missing or extra:
        raise ValueError(
            "Source-group map must match collinearity sources exactly; "
            f"missing={sorted(missing)}, extra={sorted(extra)}"
        )
    return mapping, "explicit_group", {
        "file_name": group_path.name,
        "sha256": _file_sha256(group_path),
    }


def _read_gene_subgenomes(
    path: str | Path, subgenomes: tuple[str, str]
) -> dict[str, str]:
    genes: dict[str, str] = {}
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or not {"gene_id", "subgenome"} <= set(
            reader.fieldnames
        ):
            raise ValueError(
                "Gene evidence requires gene_id and subgenome columns"
            )
        for line_number, row in enumerate(reader, start=2):
            gene_id = row["gene_id"]
            subgenome = row["subgenome"]
            if not gene_id or gene_id in genes:
                raise ValueError(
                    f"Empty or duplicate gene_id at line {line_number}"
                )
            if subgenome in subgenomes:
                genes[gene_id] = subgenome
    if not genes:
        raise ValueError("Gene evidence contains no genes in requested subgenomes")
    return genes


def _pair_id(wgd_event: str, gene_a: str, gene_b: str) -> str:
    digest = hashlib.sha256(
        f"{wgd_event}\0{gene_a}\0{gene_b}".encode("utf-8")
    ).hexdigest()[:20]
    return f"PPW-{digest}"


def infer_wgdi_homeolog_pairs(
    *,
    gene_evidence_tsv_path: str | Path,
    collinearity_inputs: Iterable[str],
    output_pair_tsv_path: str | Path,
    decisions_tsv_path: str | Path,
    wgd_event: str,
    subgenomes: tuple[str, str],
    source_group_map_path: str | Path | None = None,
    min_support_group_count: int = 2,
    require_reciprocal_unique: bool = True,
) -> dict[str, Any]:
    """Infer conservative named-WGD pairs from independent WGDI anchors.

    A source supports a pair only when one counterpart gene has exactly one
    query gene from each requested subgenome. Support is then collapsed by an
    optional source-group map. The default reciprocal-unique gate rejects a
    gene if the surviving evidence connects it to more than one partner.
    """

    if not WGD_EVENT_LABEL.fullmatch(wgd_event):
        raise ValueError(f"Unsafe WGD event label: {wgd_event}")
    if len(subgenomes) != 2 or any(not value for value in subgenomes):
        raise ValueError("Exactly two non-empty subgenome labels are required")
    if subgenomes[0] == subgenomes[1]:
        raise ValueError("Subgenome labels must be distinct")
    if min_support_group_count < 1:
        raise ValueError("min_support_group_count must be positive")

    output_pair = Path(output_pair_tsv_path)
    decisions = Path(decisions_tsv_path)
    manifest_path = Path(str(output_pair) + ".manifest.json")
    existing = [path for path in (output_pair, decisions, manifest_path) if path.exists()]
    if existing:
        raise FileExistsError(
            "Refusing to overwrite homeolog-pair artifact(s): "
            + ", ".join(str(path) for path in existing)
        )

    sources_and_paths = _parse_collinearity_inputs(collinearity_inputs)
    source_names = tuple(source for source, _ in sources_and_paths)
    source_groups, support_unit, group_manifest = _read_source_groups(
        source_group_map_path, source_names
    )
    genes = _read_gene_subgenomes(gene_evidence_tsv_path, subgenomes)
    first_subgenome, second_subgenome = subgenomes

    pair_sources: dict[tuple[str, str], set[str]] = defaultdict(set)
    pair_counterparts: dict[tuple[str, str], dict[str, set[str]]] = defaultdict(
        lambda: defaultdict(set)
    )
    input_collinearity: dict[str, dict[str, Any]] = {}
    promiscuous_counterpart_anchors = 0
    unique_counterpart_anchors = 0
    for source, path in sources_and_paths:
        blocks = parse_wgdi_collinearity(path, source)
        counterpart_queries: dict[str, set[str]] = defaultdict(set)
        unexpected: set[str] = set()
        for block in blocks:
            for query_gene, counterpart_gene in block.pairs:
                if query_gene not in genes:
                    unexpected.add(query_gene)
                    continue
                counterpart_queries[counterpart_gene].add(query_gene)
        if unexpected:
            raise ValueError(
                f"WGDI source {source} contains query genes absent from requested "
                f"gene evidence: {', '.join(sorted(unexpected)[:10])}"
            )
        for counterpart, query_genes in counterpart_queries.items():
            first = sorted(
                gene for gene in query_genes if genes[gene] == first_subgenome
            )
            second = sorted(
                gene for gene in query_genes if genes[gene] == second_subgenome
            )
            if not first or not second:
                continue
            if len(first) != 1 or len(second) != 1:
                promiscuous_counterpart_anchors += 1
                continue
            unique_counterpart_anchors += 1
            pair = (first[0], second[0])
            pair_sources[pair].add(source)
            pair_counterparts[pair][source].add(counterpart)
        input_collinearity[source] = {
            "file_name": path.name,
            "sha256": _file_sha256(path),
            "blocks": len(blocks),
            "pair_rows": sum(len(block.pairs) for block in blocks),
            "counterpart_genes": len(counterpart_queries),
        }

    support_groups_by_pair = {
        pair: {source_groups[source] for source in sources}
        for pair, sources in pair_sources.items()
    }
    support_passing = {
        pair
        for pair, groups in support_groups_by_pair.items()
        if len(groups) >= min_support_group_count
    }
    partners: dict[str, set[str]] = defaultdict(set)
    for gene_a, gene_b in support_passing:
        partners[gene_a].add(gene_b)
        partners[gene_b].add(gene_a)

    fields = (
        "pair_id",
        "gene_id_a",
        "gene_id_b",
        "wgd_event",
        "subgenome_a",
        "subgenome_b",
        "support_source_count",
        "support_sources",
        "support_group_count",
        "support_groups",
        "counterpart_anchors",
    )
    decision_fields = (*fields, "status", "reason")
    accepted_rows: list[dict[str, str]] = []
    decision_rows: list[dict[str, str]] = []
    for pair in sorted(pair_sources):
        gene_a, gene_b = pair
        sources = sorted(pair_sources[pair])
        groups = sorted(support_groups_by_pair[pair])
        counterpart_text = ";".join(
            f"{source}={','.join(sorted(pair_counterparts[pair][source]))}"
            for source in sources
        )
        row = {
            "pair_id": _pair_id(wgd_event, gene_a, gene_b),
            "gene_id_a": gene_a,
            "gene_id_b": gene_b,
            "wgd_event": wgd_event,
            "subgenome_a": first_subgenome,
            "subgenome_b": second_subgenome,
            "support_source_count": str(len(sources)),
            "support_sources": ",".join(sources),
            "support_group_count": str(len(groups)),
            "support_groups": ",".join(groups),
            "counterpart_anchors": counterpart_text,
        }
        if pair not in support_passing:
            status, reason = "rejected", "support_below_threshold"
        elif require_reciprocal_unique and (
            len(partners[gene_a]) != 1 or len(partners[gene_b]) != 1
        ):
            status, reason = "rejected", "nonreciprocal_multiple_partners"
        else:
            status, reason = "accepted", "independent_unique_wgdi_support"
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
        "schema_version": HOMEOLOG_PAIR_MANIFEST_SCHEMA_VERSION,
        "pair_schema_version": HOMEOLOG_PAIR_SCHEMA_VERSION,
        "generator": {"name": "PloidyPatch", "version": __version__},
        "inputs": {
            "gene_evidence": {
                "file_name": Path(gene_evidence_tsv_path).name,
                "sha256": _file_sha256(gene_evidence_tsv_path),
                "requested_subgenome_genes": len(genes),
            },
            "collinearity": input_collinearity,
            "source_group_map": group_manifest,
        },
        "parameters": {
            "wgd_event": wgd_event,
            "subgenomes": list(subgenomes),
            "min_support_group_count": min_support_group_count,
            "support_unit": support_unit,
            "within_source_counterpart_policy": "exactly_one_gene_per_subgenome",
            "require_reciprocal_unique": require_reciprocal_unique,
        },
        "counts": {
            "unique_counterpart_anchors": unique_counterpart_anchors,
            "promiscuous_counterpart_anchors_rejected": (
                promiscuous_counterpart_anchors
            ),
            "candidate_pairs": len(pair_sources),
            "support_passing_pairs": len(support_passing),
            "accepted_pairs": len(accepted_rows),
            "accepted_genes": len(
                {
                    gene
                    for row in accepted_rows
                    for gene in (row["gene_id_a"], row["gene_id_b"])
                }
            ),
            "decision_reason_counts": dict(sorted(reason_counts.items())),
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
        "benchmark_contract": (
            "accepted pair TSV is eligible as evaluator-only input for "
            "annotation_copy_collapse"
        ),
    }
    with manifest_path.open("x", encoding="utf-8", newline="") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return manifest


def _read_wgdi_query_genes(path: str | Path) -> dict[str, tuple[str, int]]:
    genes: dict[str, tuple[str, int]] = {}
    with Path(path).open("r", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            fields = raw_line.rstrip("\r\n").split("\t")
            if len(fields) != 6:
                raise ValueError(f"Malformed WGDI GFF line {line_number}")
            seqid, gene_id, _, _, _, order_text = fields
            if not seqid or not gene_id or gene_id in genes:
                raise ValueError(
                    f"Empty or duplicate query gene on WGDI GFF line {line_number}"
                )
            try:
                order = int(order_text)
            except ValueError as exc:
                raise ValueError(
                    f"Non-integer gene order on WGDI GFF line {line_number}"
                ) from exc
            genes[gene_id] = (seqid, order)
    if not genes:
        raise ValueError("WGDI query GFF contains no genes")
    return genes


def infer_outgroup_duplicated_pairs(
    *,
    query_wgdi_gff_path: str | Path,
    source_gff_path: str | Path | None = None,
    collinearity_inputs: Iterable[str],
    output_pair_tsv_path: str | Path,
    decisions_tsv_path: str | Path,
    wgd_event: str,
    source_group_map_path: str | Path | None = None,
    min_support_group_count: int = 2,
    min_block_pairs: int = 20,
    require_cross_seqid: bool = True,
    require_reciprocal_unique: bool = True,
) -> dict[str, Any]:
    """Infer duplicated query descendants supported one-to-two by outgroups.

    Within each outgroup source, a counterpart supports a pair only when it is
    connected in qualifying collinear blocks to exactly two query genes.  Pair
    support is collapsed by declared source group, then the same pair must pass
    the support threshold and, by default, reciprocal uniqueness.  The
    cross-seqid gate excludes tandem/local duplicates from WGD truth sets.
    """

    if not WGD_EVENT_LABEL.fullmatch(wgd_event):
        raise ValueError(f"Unsafe WGD event label: {wgd_event}")
    if min_support_group_count < 1 or min_block_pairs < 2:
        raise ValueError("Support count must be positive and block size at least two")
    output_pair = Path(output_pair_tsv_path)
    decisions = Path(decisions_tsv_path)
    manifest_path = Path(str(output_pair) + ".manifest.json")
    existing = [path for path in (output_pair, decisions, manifest_path) if path.exists()]
    if existing:
        raise FileExistsError(
            "Refusing to overwrite outgroup-duplicate artifact(s): "
            + ", ".join(str(path) for path in existing)
        )

    genes = _read_wgdi_query_genes(query_wgdi_gff_path)
    feature_id_by_wgdi_id, source_gff_manifest = _map_wgdi_to_feature_ids(
        source_gff_path, set(genes)
    )
    sources_and_paths = _parse_collinearity_inputs(collinearity_inputs)
    source_names = tuple(source for source, _ in sources_and_paths)
    source_groups, support_unit, group_manifest = _read_source_groups(
        source_group_map_path, source_names
    )
    pair_sources: dict[tuple[str, str], set[str]] = defaultdict(set)
    pair_counterparts: dict[tuple[str, str], dict[str, set[str]]] = defaultdict(
        lambda: defaultdict(set)
    )
    input_collinearity: dict[str, dict[str, Any]] = {}
    counts: dict[str, int] = defaultdict(int)
    for source, path in sources_and_paths:
        blocks = parse_wgdi_collinearity(path, source)
        qualifying = [block for block in blocks if len(block.pairs) >= min_block_pairs]
        counterpart_queries: dict[str, set[str]] = defaultdict(set)
        unexpected: set[str] = set()
        for block in qualifying:
            for query_gene, counterpart_gene in block.pairs:
                if query_gene not in genes:
                    unexpected.add(query_gene)
                    continue
                counterpart_queries[counterpart_gene].add(query_gene)
        if unexpected:
            raise ValueError(
                f"WGDI source {source} contains query genes absent from query GFF: "
                + ", ".join(sorted(unexpected)[:10])
            )
        for counterpart, query_genes in counterpart_queries.items():
            if len(query_genes) < 2:
                counts["singleton_counterpart_anchors"] += 1
                continue
            if len(query_genes) > 2:
                counts["promiscuous_counterpart_anchors_rejected"] += 1
                continue
            pair = tuple(sorted(query_genes))
            if require_cross_seqid and genes[pair[0]][0] == genes[pair[1]][0]:
                counts["same_seqid_counterpart_anchors_rejected"] += 1
                continue
            counts["unique_pair_counterpart_anchors"] += 1
            pair_sources[pair].add(source)
            pair_counterparts[pair][source].add(counterpart)
        input_collinearity[source] = {
            "file_name": path.name,
            "sha256": _file_sha256(path),
            "blocks": len(blocks),
            "qualifying_blocks": len(qualifying),
            "pair_rows": sum(len(block.pairs) for block in blocks),
            "qualifying_pair_rows": sum(len(block.pairs) for block in qualifying),
            "counterpart_genes_in_qualifying_blocks": len(counterpart_queries),
        }

    support_groups_by_pair = {
        pair: {source_groups[source] for source in sources}
        for pair, sources in pair_sources.items()
    }
    support_passing = {
        pair
        for pair, groups in support_groups_by_pair.items()
        if len(groups) >= min_support_group_count
    }
    partners: dict[str, set[str]] = defaultdict(set)
    for gene_a, gene_b in support_passing:
        partners[gene_a].add(gene_b)
        partners[gene_b].add(gene_a)

    fields = (
        "pair_id",
        "gene_id_a",
        "gene_id_b",
        "wgdi_gene_id_a",
        "wgdi_gene_id_b",
        "seqid_a",
        "seqid_b",
        "gene_order_a",
        "gene_order_b",
        "wgd_event",
        "support_source_count",
        "support_sources",
        "support_group_count",
        "support_groups",
        "counterpart_anchors",
    )
    decision_fields = (*fields, "status", "reason")
    accepted_rows: list[dict[str, str]] = []
    decision_rows: list[dict[str, str]] = []
    for pair in sorted(pair_sources):
        wgdi_gene_a, wgdi_gene_b = pair
        gene_a = feature_id_by_wgdi_id[wgdi_gene_a]
        gene_b = feature_id_by_wgdi_id[wgdi_gene_b]
        sources = sorted(pair_sources[pair])
        groups = sorted(support_groups_by_pair[pair])
        row = {
            "pair_id": _pair_id(wgd_event, gene_a, gene_b),
            "gene_id_a": gene_a,
            "gene_id_b": gene_b,
            "wgdi_gene_id_a": wgdi_gene_a,
            "wgdi_gene_id_b": wgdi_gene_b,
            "seqid_a": genes[wgdi_gene_a][0],
            "seqid_b": genes[wgdi_gene_b][0],
            "gene_order_a": str(genes[wgdi_gene_a][1]),
            "gene_order_b": str(genes[wgdi_gene_b][1]),
            "wgd_event": wgd_event,
            "support_source_count": str(len(sources)),
            "support_sources": ",".join(sources),
            "support_group_count": str(len(groups)),
            "support_groups": ",".join(groups),
            "counterpart_anchors": ";".join(
                f"{source}={','.join(sorted(pair_counterparts[pair][source]))}"
                for source in sources
            ),
        }
        if pair not in support_passing:
            status, reason = "rejected", "support_below_threshold"
        elif require_reciprocal_unique and (
            len(partners[wgdi_gene_a]) != 1
            or len(partners[wgdi_gene_b]) != 1
        ):
            status, reason = "rejected", "nonreciprocal_multiple_partners"
        else:
            status, reason = "accepted", "independent_one_to_two_outgroup_support"
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
            handle, fieldnames=decision_fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(decision_rows)

    reason_counts = Counter(row["reason"] for row in decision_rows)
    manifest: dict[str, Any] = {
        "schema_version": OUTGROUP_DUPLICATE_MANIFEST_SCHEMA_VERSION,
        "pair_schema_version": OUTGROUP_DUPLICATE_PAIR_SCHEMA_VERSION,
        "generator": {"name": "PloidyPatch", "version": __version__},
        "access": "evaluator_only",
        "inputs": {
            "query_wgdi_gff": {
                "file_name": Path(query_wgdi_gff_path).name,
                "sha256": _file_sha256(query_wgdi_gff_path),
                "genes": len(genes),
            },
            "source_gff_mapping": source_gff_manifest,
            "collinearity": input_collinearity,
            "source_group_map": group_manifest,
        },
        "parameters": {
            "wgd_event": wgd_event,
            "min_support_group_count": min_support_group_count,
            "support_unit": support_unit,
            "min_block_pairs": min_block_pairs,
            "within_source_counterpart_policy": "exactly_two_query_genes",
            "require_cross_seqid": require_cross_seqid,
            "require_reciprocal_unique": require_reciprocal_unique,
        },
        "counts": {
            **dict(sorted(counts.items())),
            "candidate_pairs": len(pair_sources),
            "support_passing_pairs": len(support_passing),
            "accepted_pairs": len(accepted_rows),
            "accepted_genes": len(
                {
                    gene
                    for row in accepted_rows
                    for gene in (row["gene_id_a"], row["gene_id_b"])
                }
            ),
            "decision_reason_counts": dict(sorted(reason_counts.items())),
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
        "benchmark_contract": (
            "accepted pair TSV is evaluator-only input for annotation_copy_collapse"
        ),
    }
    with manifest_path.open("x", encoding="utf-8", newline="") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return manifest
