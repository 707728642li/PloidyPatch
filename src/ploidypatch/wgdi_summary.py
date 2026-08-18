from __future__ import annotations

import csv
import hashlib
import json
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from . import __version__


WGDI_SUMMARY_SCHEMA = "ploidypatch.wgdi_gene_evidence.v1"
WGDI_SUMMARY_MANIFEST_SCHEMA = "ploidypatch.wgdi_summary_manifest.v1"
SOURCE_LABEL = re.compile(r"^[A-Za-z][A-Za-z0-9_]*$")
BLOCK_HEADER = re.compile(
    r"^# Alignment (?P<block_id>[^:]+): score=(?P<score>\S+) "
    r"pvalue=(?P<pvalue>\S+) N=(?P<pair_count>\d+) "
    r"(?P<query_seqid>[^&\s]+)&(?P<target_seqid>\S+) "
    r"(?P<orientation>plus|minus)$"
)


@dataclass(frozen=True)
class WgdiBlock:
    source: str
    block_id: str
    score: float
    pvalue: float
    declared_pairs: int
    query_seqid: str
    target_seqid: str
    orientation: str
    pairs: tuple[tuple[str, str], ...]


@dataclass
class GeneEvidence:
    block_ids: set[str]
    counterpart_genes: set[str]
    best_block_score: float | None = None
    best_block_pvalue: float | None = None
    longest_block_pairs: int = 0

    @classmethod
    def empty(cls) -> GeneEvidence:
        return cls(block_ids=set(), counterpart_genes=set())

    def add(self, block: WgdiBlock, counterpart_gene: str) -> None:
        self.block_ids.add(f"{block.source}:{block.block_id}")
        self.counterpart_genes.add(counterpart_gene)
        if self.best_block_score is None or block.score > self.best_block_score:
            self.best_block_score = block.score
        if self.best_block_pvalue is None or block.pvalue < self.best_block_pvalue:
            self.best_block_pvalue = block.pvalue
        self.longest_block_pairs = max(self.longest_block_pairs, len(block.pairs))


def _file_sha256(path: str | Path, block_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        while block := handle.read(block_size):
            digest.update(block)
    return digest.hexdigest()


def parse_wgdi_collinearity(path: str | Path, source: str) -> list[WgdiBlock]:
    """Parse the text output produced by ``wgdi -icl``."""

    if not SOURCE_LABEL.fullmatch(source):
        raise ValueError(f"Unsafe WGDI source label: {source}")
    blocks: list[WgdiBlock] = []
    header: dict[str, str] | None = None
    pairs: list[tuple[str, str]] = []

    def finish_block() -> None:
        nonlocal header, pairs
        if header is None:
            return
        block = WgdiBlock(
            source=source,
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
                f"WGDI block {source}:{block.block_id} declares "
                f"{block.declared_pairs} pairs but contains {len(block.pairs)}"
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
                finish_block()
                match = BLOCK_HEADER.fullmatch(line)
                if match is None:
                    raise ValueError(
                        f"Malformed WGDI block header on line {line_number}: {line}"
                    )
                header = match.groupdict()
                continue
            if header is None:
                raise ValueError(
                    f"WGDI pair encountered before a block header on line {line_number}"
                )
            fields = line.split()
            if len(fields) != 5:
                raise ValueError(
                    f"Malformed WGDI pair on line {line_number}: expected 5 fields"
                )
            pairs.append((fields[0], fields[2]))
    finish_block()
    return blocks


def _read_query_gff(path: str | Path) -> dict[str, tuple[str, int]]:
    genes: dict[str, tuple[str, int]] = {}
    with Path(path).open("r", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            fields = raw_line.rstrip("\r\n").split("\t")
            if len(fields) != 6:
                raise ValueError(f"Malformed WGDI GFF line {line_number}")
            seqid, gene_id, _, _, _, order_text = fields
            if gene_id in genes:
                raise ValueError(f"Duplicate gene in WGDI GFF: {gene_id}")
            genes[gene_id] = (seqid, int(order_text))
    return genes


def _source_columns(source: str) -> tuple[str, ...]:
    return (
        f"{source}_collinear",
        f"{source}_block_count",
        f"{source}_counterpart_gene_count",
        f"{source}_best_block_score",
        f"{source}_best_block_pvalue",
        f"{source}_longest_block_pairs",
    )


def summarize_wgdi_gene_evidence(
    query_gff_path: str | Path,
    query_input_manifest_path: str | Path,
    collinearity_inputs: list[str],
    expected_source_inputs: list[str],
    output_gene_tsv_path: str | Path,
    output_block_tsv_path: str | Path,
) -> dict[str, Any]:
    """Summarize external WGDI blocks into evaluator-ready gene strata."""

    sources_and_paths: list[tuple[str, Path]] = []
    seen_sources: set[str] = set()
    for item in collinearity_inputs:
        source, separator, raw_path = item.partition("=")
        if not separator or not source or not raw_path:
            raise ValueError("Each --collinearity value must be SOURCE=PATH")
        if not SOURCE_LABEL.fullmatch(source):
            raise ValueError(f"Unsafe WGDI source label: {source}")
        if source in seen_sources:
            raise ValueError(f"Duplicate WGDI source label: {source}")
        seen_sources.add(source)
        sources_and_paths.append((source, Path(raw_path)))
    if not sources_and_paths:
        raise ValueError("At least one WGDI collinearity input is required")
    expected_source_by_subgenome: dict[str, str] = {}
    for item in expected_source_inputs:
        subgenome, separator, source = item.partition("=")
        if not separator or not subgenome or not source:
            raise ValueError("Each --expected-source value must be SUBGENOME=SOURCE")
        if subgenome in expected_source_by_subgenome:
            raise ValueError(f"Duplicate expected-source subgenome: {subgenome}")
        if source not in seen_sources:
            raise ValueError(
                f"Expected source {source} is absent from --collinearity inputs"
            )
        expected_source_by_subgenome[subgenome] = source

    output_gene_path = Path(output_gene_tsv_path)
    output_block_path = Path(output_block_tsv_path)
    manifest_path = Path(str(output_gene_path) + ".manifest.json")
    collisions = [
        path
        for path in (output_gene_path, output_block_path, manifest_path)
        if path.exists()
    ]
    if collisions:
        raise FileExistsError(
            "Refusing to overwrite WGDI summary artifact(s): "
            + ", ".join(str(path) for path in collisions)
        )

    query_genes = _read_query_gff(query_gff_path)
    query_manifest = json.loads(
        Path(query_input_manifest_path).read_text(encoding="utf-8")
    )
    chromosome_labels = query_manifest.get("chromosome_labels")
    if chromosome_labels is None:
        chromosome_labels = query_manifest.get("selection", {}).get(
            "chromosome_labels", {}
        )
    if not isinstance(chromosome_labels, dict):
        raise ValueError("Query input manifest chromosome_labels must be a mapping")
    declared_subgenomes = sorted(
        expected_source_by_subgenome, key=lambda value: (-len(value), value)
    )
    all_blocks: list[WgdiBlock] = []
    evidence: dict[str, dict[str, GeneEvidence]] = {
        source: defaultdict(GeneEvidence.empty) for source, _ in sources_and_paths
    }
    input_manifest: dict[str, Any] = {}
    for source, path in sources_and_paths:
        blocks = parse_wgdi_collinearity(path, source)
        all_blocks.extend(blocks)
        unexpected_genes: set[str] = set()
        for block in blocks:
            for query_gene, counterpart_gene in block.pairs:
                if query_gene not in query_genes:
                    unexpected_genes.add(query_gene)
                    continue
                evidence[source][query_gene].add(block, counterpart_gene)
        if unexpected_genes:
            examples = ", ".join(sorted(unexpected_genes)[:10])
            raise ValueError(
                f"WGDI source {source} contains query genes absent from query GFF: "
                f"{examples}"
            )
        input_manifest[source] = {
            "file_name": path.name,
            "sha256": _file_sha256(path),
            "blocks": len(blocks),
            "pair_rows": sum(len(block.pairs) for block in blocks),
        }

    source_names = [source for source, _ in sources_and_paths]
    fixed_columns = [
        "gene_id",
        "seqid",
        "gene_order",
        "chromosome_label",
        "subgenome",
    ]
    evidence_columns = [
        column for source in source_names for column in _source_columns(source)
    ]
    final_columns = [
        "progenitor_support_class",
        "expected_progenitor",
        "expected_progenitor_supported",
        "cross_progenitor_supported",
        "synteny_stratum",
    ]
    gene_rows: list[dict[str, str]] = []
    for gene_id, (seqid, order) in sorted(
        query_genes.items(), key=lambda item: (item[1][0], item[1][1], item[0])
    ):
        chromosome_label = chromosome_labels.get(seqid, "")
        subgenome = next(
            (
                value
                for value in declared_subgenomes
                if chromosome_label.startswith(value)
            ),
            "",
        )
        row = {
            "gene_id": gene_id,
            "seqid": seqid,
            "gene_order": str(order),
            "chromosome_label": chromosome_label,
            "subgenome": subgenome,
        }
        supported_sources: list[str] = []
        for source in source_names:
            item = evidence[source].get(gene_id, GeneEvidence.empty())
            is_supported = bool(item.block_ids)
            if is_supported:
                supported_sources.append(source)
            row.update(
                {
                    f"{source}_collinear": "true" if is_supported else "false",
                    f"{source}_block_count": str(len(item.block_ids)),
                    f"{source}_counterpart_gene_count": str(
                        len(item.counterpart_genes)
                    ),
                    f"{source}_best_block_score": ""
                    if item.best_block_score is None
                    else str(item.best_block_score),
                    f"{source}_best_block_pvalue": ""
                    if item.best_block_pvalue is None
                    else str(item.best_block_pvalue),
                    f"{source}_longest_block_pairs": str(item.longest_block_pairs),
                }
            )

        support_class = (
            "none" if not supported_sources else "+".join(sorted(supported_sources))
        )
        expected = expected_source_by_subgenome.get(subgenome, "")
        expected_supported = bool(expected and expected in supported_sources)
        cross_supported = bool(
            expected and any(source != expected for source in supported_sources)
        )
        if expected_supported and cross_supported:
            stratum = "expected_and_cross"
        elif expected_supported:
            stratum = "expected_only"
        elif cross_supported:
            stratum = "cross_only"
        else:
            stratum = "no_block"
        row.update(
            {
                "progenitor_support_class": support_class,
                "expected_progenitor": expected,
                "expected_progenitor_supported": "true"
                if expected_supported
                else "false",
                "cross_progenitor_supported": "true" if cross_supported else "false",
                "synteny_stratum": stratum,
            }
        )
        gene_rows.append(row)

    output_gene_path.parent.mkdir(parents=True, exist_ok=True)
    with output_gene_path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=(*fixed_columns, *evidence_columns, *final_columns),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(gene_rows)

    block_columns = (
        "source",
        "block_id",
        "score",
        "pvalue",
        "pair_count",
        "query_seqid",
        "target_seqid",
        "orientation",
    )
    with output_block_path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=block_columns,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for block in all_blocks:
            writer.writerow(
                {
                    "source": block.source,
                    "block_id": block.block_id,
                    "score": block.score,
                    "pvalue": block.pvalue,
                    "pair_count": len(block.pairs),
                    "query_seqid": block.query_seqid,
                    "target_seqid": block.target_seqid,
                    "orientation": block.orientation,
                }
            )

    stratum_counts: dict[str, int] = defaultdict(int)
    for row in gene_rows:
        stratum_counts[row["synteny_stratum"]] += 1
    manifest: dict[str, Any] = {
        "schema_version": WGDI_SUMMARY_MANIFEST_SCHEMA,
        "gene_evidence_schema_version": WGDI_SUMMARY_SCHEMA,
        "generator": {"name": "PloidyPatch", "version": __version__},
        "inputs": {
            "query_gff": {
                "file_name": Path(query_gff_path).name,
                "sha256": _file_sha256(query_gff_path),
            },
            "query_input_manifest": {
                "file_name": Path(query_input_manifest_path).name,
                "sha256": _file_sha256(query_input_manifest_path),
            },
            "collinearity": input_manifest,
        },
        "outputs": {
            "gene_evidence": {
                "file_name": output_gene_path.name,
                "sha256": _file_sha256(output_gene_path),
                "rows": len(gene_rows),
            },
            "blocks": {
                "file_name": output_block_path.name,
                "sha256": _file_sha256(output_block_path),
                "rows": len(all_blocks),
            },
        },
        "synteny_stratum_counts": dict(sorted(stratum_counts.items())),
        "expected_source_by_subgenome": dict(
            sorted(expected_source_by_subgenome.items())
        ),
        "subgenome_assignment": {
            "method": "longest_declared_prefix_of_chromosome_label",
            "declared_subgenomes": declared_subgenomes,
        },
        "access": "evaluator_only",
    }
    with manifest_path.open("x", encoding="utf-8", newline="") as handle:
        handle.write(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return manifest
