from __future__ import annotations

import csv
import hashlib
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

from . import __version__
from .perturb import (
    BOUNDARY_SHIFT_EVENT,
    COPY_COLLAPSE_EVENT,
    FUSED_GENE_EVENT,
    MANIFEST_SCHEMA_VERSION,
    MISSING_INTERNAL_EXON_EVENT,
    SPLIT_GENE_EVENT,
    TRUTH_SCHEMA_VERSION,
    GffRecord,
    MissingGeneCandidate,
    _file_sha256,
    _json_text,
    _text_sha256,
    _write_text_exclusive,
    find_missing_gene_candidates,
    read_gff_document,
)


STRUCTURE_EVENT_TYPES = frozenset(
    {
        MISSING_INTERNAL_EXON_EVENT,
        BOUNDARY_SHIFT_EVENT,
        SPLIT_GENE_EVENT,
        FUSED_GENE_EVENT,
        COPY_COLLAPSE_EVENT,
    }
)


@dataclass(frozen=True)
class SimpleGene:
    hierarchy: MissingGeneCandidate
    transcript: GffRecord
    exons: tuple[GffRecord, ...]
    cds: tuple[GffRecord, ...]
    direct_children: tuple[GffRecord, ...]

    @property
    def gene(self) -> GffRecord:
        return self.hierarchy.gene

    @property
    def gene_id(self) -> str:
        assert self.gene.feature_id is not None
        return self.gene.feature_id

    @property
    def transcript_id(self) -> str:
        assert self.transcript.feature_id is not None
        return self.transcript.feature_id


@dataclass(frozen=True)
class EventOption:
    key: str
    genes: tuple[str, ...]
    payload: Any


def _simple_genes(document: Any) -> dict[str, SimpleGene]:
    simple: dict[str, SimpleGene] = {}
    for hierarchy in find_missing_gene_candidates(document, require_exon=False):
        transcript_records = tuple(
            record
            for record in hierarchy.removed_records
            if record.feature_type in {"mRNA", "transcript"}
        )
        if len(transcript_records) != 1:
            continue
        transcript = transcript_records[0]
        transcript_id = transcript.feature_id
        if transcript_id is None or transcript.parents != (hierarchy.gene.feature_id,):
            continue
        direct_children = tuple(
            record
            for record in hierarchy.removed_records
            if record.parents == (transcript_id,)
        )
        exons = tuple(
            sorted(
                (record for record in direct_children if record.feature_type == "exon"),
                key=lambda record: (record.start, record.end, record.line_number),
            )
        )
        cds = tuple(
            sorted(
                (record for record in direct_children if record.feature_type == "CDS"),
                key=lambda record: (record.start, record.end, record.line_number),
            )
        )
        if len(cds) < 1:
            continue
        simple[hierarchy.gene.feature_id] = SimpleGene(
            hierarchy=hierarchy,
            transcript=transcript,
            exons=exons,
            cds=cds,
            direct_children=direct_children,
        )
    return simple


def _rank(seed: int, event_type: str, key: str) -> str:
    return hashlib.sha256(f"{seed}\0{event_type}\0{key}".encode()).hexdigest()


def copy_collapse_partner_ids(
    gene_a: str, gene_b: str, *, seed: int
) -> tuple[str, str]:
    """Return the deterministically collapsed and retained copy IDs."""

    if not gene_a or not gene_b or gene_a == gene_b:
        raise ValueError("Copy-collapse partners must be distinct non-empty IDs")
    first, second = sorted((gene_a, gene_b))
    key = f"{first}|{second}"
    return min(
        ((first, second), (second, first)),
        key=lambda pair: _rank(
            seed, COPY_COLLAPSE_EVENT, f"{key}|{pair[0]}"
        ),
    )


def _token(source_sha: str, seed: int, event_type: str, key: str) -> str:
    return hashlib.sha256(
        f"{source_sha}\0{seed}\0{event_type}\0{key}".encode()
    ).hexdigest()[:20]


def _event_id(source_sha: str, seed: int, event_type: str, key: str) -> str:
    prefixes = {
        MISSING_INTERNAL_EXON_EVENT: "ME",
        BOUNDARY_SHIFT_EVENT: "BS",
        SPLIT_GENE_EVENT: "SG",
        FUSED_GENE_EVENT: "FG",
        COPY_COLLAPSE_EVENT: "CC",
    }
    return f"PPB-{prefixes[event_type]}-{_token(source_sha, seed, event_type, key)}"


def _line_ending(raw_line: str) -> str:
    if raw_line.endswith("\r\n"):
        return "\r\n"
    if raw_line.endswith("\n"):
        return "\n"
    return ""


def _replace_attribute(text: str, key: str, value: str) -> str:
    fields = text.split(";") if text else []
    replaced = False
    output: list[str] = []
    for field in fields:
        if field.startswith(f"{key}="):
            output.append(f"{key}={value}")
            replaced = True
        elif field:
            output.append(field)
    if not replaced:
        output.append(f"{key}={value}")
    return ";".join(output)


def _rewrite_record(
    record: GffRecord,
    *,
    start: int | None = None,
    end: int | None = None,
    feature_id: str | None = None,
    parent: str | None = None,
) -> str:
    ending = _line_ending(record.raw_line)
    fields = record.raw_line.rstrip("\r\n").split("\t")
    if start is not None:
        fields[3] = str(start)
    if end is not None:
        fields[4] = str(end)
    attributes = fields[8]
    if feature_id is not None:
        attributes = _replace_attribute(attributes, "ID", feature_id)
    if parent is not None:
        attributes = _replace_attribute(attributes, "Parent", parent)
    fields[8] = attributes
    return "\t".join(fields) + ending


def _line_edit(record: GffRecord, replacements: Iterable[str]) -> dict[str, Any]:
    perturbed = []
    for raw_line in replacements:
        perturbed.append(
            {
                "line_sha256": hashlib.sha256(raw_line.encode()).hexdigest(),
                "raw_line": raw_line,
            }
        )
    return {
        "source_line_number": record.line_number,
        "source_line_sha256": hashlib.sha256(record.raw_line.encode()).hexdigest(),
        "source_raw_line": record.raw_line,
        "perturbed_lines": perturbed,
    }


def _target(
    source_genes: tuple[SimpleGene, ...],
    perturbed_gene_ids: tuple[str, ...],
) -> dict[str, Any]:
    starts = [gene.gene.start for gene in source_genes]
    ends = [gene.gene.end for gene in source_genes]
    return {
        "gene_id": source_genes[0].gene_id,
        "gene_ids": [gene.gene_id for gene in source_genes],
        "transcript_ids": [gene.transcript_id for gene in source_genes],
        "perturbed_gene_ids": list(perturbed_gene_ids),
        "seqid": source_genes[0].gene.seqid,
        "start": min(starts),
        "end": max(ends),
        "strand": source_genes[0].gene.strand,
    }


def _missing_exon_options(genes: dict[str, SimpleGene]) -> list[EventOption]:
    options: list[EventOption] = []
    for gene in genes.values():
        if len(gene.exons) < 3 or len(gene.cds) < 3:
            continue
        choices: list[tuple[GffRecord, GffRecord]] = []
        for exon in gene.exons[1:-1]:
            exact_cds = [
                cds
                for cds in gene.cds
                if (cds.start, cds.end) == (exon.start, exon.end)
            ]
            overlaps = [
                child
                for child in gene.direct_children
                if child.line_number != exon.line_number
                and child.start <= exon.end
                and child.end >= exon.start
            ]
            if len(exact_cds) == 1 and overlaps == exact_cds:
                choices.append((exon, exact_cds[0]))
        if choices:
            options.append(EventOption(gene.gene_id, (gene.gene_id,), (gene, choices)))
    return options


def _boundary_options(genes: dict[str, SimpleGene]) -> list[EventOption]:
    options: list[EventOption] = []
    for gene in genes.values():
        if not gene.exons:
            continue
        choices: list[tuple[str, GffRecord, GffRecord]] = []
        for side, exon in (("left", gene.exons[0]), ("right", gene.exons[-1])):
            exact_cds = [
                cds
                for cds in gene.cds
                if (cds.start, cds.end) == (exon.start, exon.end)
            ]
            if len(exact_cds) != 1 or exon.end - exon.start + 1 <= 18:
                continue
            shifted_boundary = (
                exon.start + 15 if side == "left" else exon.end - 15
            )
            conflicting_children = [
                child
                for child in gene.direct_children
                if child.line_number
                not in {exon.line_number, exact_cds[0].line_number}
                and (
                    (side == "left" and child.start < shifted_boundary)
                    or (side == "right" and child.end > shifted_boundary)
                )
            ]
            if conflicting_children:
                continue
            boundary = exon.start if side == "left" else exon.end
            gene_boundary = gene.gene.start if side == "left" else gene.gene.end
            transcript_boundary = (
                gene.transcript.start if side == "left" else gene.transcript.end
            )
            if boundary == gene_boundary == transcript_boundary:
                choices.append((side, exon, exact_cds[0]))
        if choices:
            options.append(EventOption(gene.gene_id, (gene.gene_id,), (gene, choices)))
    return options


def _split_options(genes: dict[str, SimpleGene]) -> list[EventOption]:
    options: list[EventOption] = []
    for gene in genes.values():
        if len(gene.exons) < 4 or len(gene.cds) < 4:
            continue
        split_index = len(gene.exons) // 2
        left_exons = gene.exons[:split_index]
        right_exons = gene.exons[split_index:]
        cut_left = left_exons[-1].end
        cut_right = right_exons[0].start
        if cut_left >= cut_right:
            continue
        if any(
            child.start <= cut_left and child.end >= cut_right
            for child in gene.direct_children
        ):
            continue
        left_children = tuple(
            child for child in gene.direct_children if child.end <= cut_left
        )
        right_children = tuple(
            child for child in gene.direct_children if child.start >= cut_right
        )
        if len(left_children) + len(right_children) != len(gene.direct_children):
            continue
        if not any(record.feature_type == "CDS" for record in left_children):
            continue
        if not any(record.feature_type == "CDS" for record in right_children):
            continue
        options.append(
            EventOption(
                gene.gene_id,
                (gene.gene_id,),
                (gene, left_children, right_children),
            )
        )
    return options


def _fusion_options(genes: dict[str, SimpleGene], document: Any) -> list[EventOption]:
    ordered = sorted(
        (record for record in document.records if record.feature_type == "gene"),
        key=lambda gene: (gene.seqid, gene.start, gene.end, gene.feature_id or ""),
    )
    options: list[EventOption] = []
    for left_record, right_record in zip(ordered, ordered[1:]):
        if left_record.feature_id not in genes or right_record.feature_id not in genes:
            continue
        left = genes[left_record.feature_id]
        right = genes[right_record.feature_id]
        if left.gene.seqid != right.gene.seqid:
            continue
        if left.gene.strand != right.gene.strand:
            continue
        gap = right.gene.start - left.gene.end - 1
        if gap < 1 or gap > 20_000:
            continue
        key = f"{left.gene_id}|{right.gene_id}"
        options.append(
            EventOption(key, (left.gene_id, right.gene_id), (left, right, gap))
        )
    return options


def _copy_options(
    genes: dict[str, SimpleGene], pair_tsv_path: str | Path | None
) -> list[EventOption]:
    if pair_tsv_path is None:
        raise ValueError("annotation_copy_collapse requires --pair-tsv")
    options: list[EventOption] = []
    with Path(pair_tsv_path).open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"gene_id_a", "gene_id_b"}
        if reader.fieldnames is None or not required <= set(reader.fieldnames):
            raise ValueError("Copy-pair TSV requires gene_id_a and gene_id_b columns")
        seen: set[tuple[str, str]] = set()
        for line_number, row in enumerate(reader, start=2):
            gene_a = row["gene_id_a"]
            gene_b = row["gene_id_b"]
            if not gene_a or not gene_b or gene_a == gene_b:
                raise ValueError(f"Invalid copy pair at line {line_number}")
            pair = tuple(sorted((gene_a, gene_b)))
            if pair in seen:
                raise ValueError(f"Duplicate copy pair at line {line_number}")
            seen.add(pair)
            if gene_a not in genes or gene_b not in genes:
                continue
            options.append(
                EventOption(
                    f"{pair[0]}|{pair[1]}",
                    pair,
                    (genes[pair[0]], genes[pair[1]]),
                )
            )
    return options


def _select_options(
    options: list[EventOption], count: int, seed: int, event_type: str
) -> list[EventOption]:
    selected: list[EventOption] = []
    occupied: set[str] = set()
    for option in sorted(options, key=lambda item: _rank(seed, event_type, item.key)):
        if occupied.intersection(option.genes):
            continue
        selected.append(option)
        occupied.update(option.genes)
        if len(selected) == count:
            break
    if len(selected) != count:
        raise ValueError(
            f"Requested {count} non-overlapping {event_type} events but only "
            f"{len(selected)} are eligible"
        )
    return selected


def _event_from_option(
    option: EventOption,
    *,
    event_type: str,
    source_sha: str,
    seed: int,
) -> dict[str, Any]:
    event_id = _event_id(source_sha, seed, event_type, option.key)
    token = event_id.rsplit("-", 1)[-1]
    edits: list[dict[str, Any]] = []
    details: dict[str, Any] = {}

    if event_type == MISSING_INTERNAL_EXON_EVENT:
        gene, choices = option.payload
        exon, cds = min(
            choices,
            key=lambda pair: _rank(seed, event_type, f"{option.key}|{pair[0].start}"),
        )
        edits.extend((_line_edit(exon, ()), _line_edit(cds, ())))
        source_genes = (gene,)
        perturbed_gene_ids = (gene.gene_id,)
        details = {
            "removed_exon": [exon.start, exon.end],
            "removed_cds": [cds.start, cds.end, cds.phase],
            "internal_exon_rank": gene.exons.index(exon),
        }
    elif event_type == BOUNDARY_SHIFT_EVENT:
        gene, choices = option.payload
        side, exon, cds = min(
            choices,
            key=lambda item: _rank(seed, event_type, f"{option.key}|{item[0]}"),
        )
        shift = 15
        if side == "left":
            replacements = {
                gene.gene.line_number: _rewrite_record(gene.gene, start=gene.gene.start + shift),
                gene.transcript.line_number: _rewrite_record(
                    gene.transcript, start=gene.transcript.start + shift
                ),
                exon.line_number: _rewrite_record(exon, start=exon.start + shift),
                cds.line_number: _rewrite_record(cds, start=cds.start + shift),
            }
        else:
            replacements = {
                gene.gene.line_number: _rewrite_record(gene.gene, end=gene.gene.end - shift),
                gene.transcript.line_number: _rewrite_record(
                    gene.transcript, end=gene.transcript.end - shift
                ),
                exon.line_number: _rewrite_record(exon, end=exon.end - shift),
                cds.line_number: _rewrite_record(cds, end=cds.end - shift),
            }
        records = {record.line_number: record for record in (gene.gene, gene.transcript, exon, cds)}
        edits.extend(
            _line_edit(records[line_number], (raw_line,))
            for line_number, raw_line in sorted(replacements.items())
        )
        source_genes = (gene,)
        perturbed_gene_ids = (gene.gene_id,)
        details = {"side_in_genomic_coordinates": side, "shift_bp": shift}
    elif event_type == SPLIT_GENE_EVENT:
        gene, left_children, right_children = option.payload
        left_gene_id = f"PPERR-SG-{token}-A"
        right_gene_id = f"PPERR-SG-{token}-B"
        left_tx_id = f"PPERR-ST-{token}-A"
        right_tx_id = f"PPERR-ST-{token}-B"
        left_start = min(record.start for record in left_children)
        left_end = max(record.end for record in left_children)
        right_start = min(record.start for record in right_children)
        right_end = max(record.end for record in right_children)
        edits.append(
            _line_edit(
                gene.gene,
                (
                    _rewrite_record(gene.gene, start=left_start, end=left_end, feature_id=left_gene_id),
                    _rewrite_record(gene.gene, start=right_start, end=right_end, feature_id=right_gene_id),
                ),
            )
        )
        edits.append(
            _line_edit(
                gene.transcript,
                (
                    _rewrite_record(gene.transcript, start=left_start, end=left_end, feature_id=left_tx_id, parent=left_gene_id),
                    _rewrite_record(gene.transcript, start=right_start, end=right_end, feature_id=right_tx_id, parent=right_gene_id),
                ),
            )
        )
        for child in left_children:
            edits.append(_line_edit(child, (_rewrite_record(child, parent=left_tx_id),)))
        for child in right_children:
            edits.append(_line_edit(child, (_rewrite_record(child, parent=right_tx_id),)))
        source_genes = (gene,)
        perturbed_gene_ids = (left_gene_id, right_gene_id)
        details = {"split_between": [left_end, right_start]}
    elif event_type == FUSED_GENE_EVENT:
        left, right, gap = option.payload
        fused_gene_id = f"PPERR-FG-{token}"
        fused_tx_id = f"PPERR-FT-{token}"
        fused_start = min(left.gene.start, right.gene.start)
        fused_end = max(left.gene.end, right.gene.end)
        edits.append(
            _line_edit(
                left.gene,
                (_rewrite_record(left.gene, start=fused_start, end=fused_end, feature_id=fused_gene_id),),
            )
        )
        edits.append(_line_edit(right.gene, ()))
        edits.append(
            _line_edit(
                left.transcript,
                (_rewrite_record(left.transcript, start=fused_start, end=fused_end, feature_id=fused_tx_id, parent=fused_gene_id),),
            )
        )
        edits.append(_line_edit(right.transcript, ()))
        for child in (*left.direct_children, *right.direct_children):
            edits.append(_line_edit(child, (_rewrite_record(child, parent=fused_tx_id),)))
        source_genes = (left, right)
        perturbed_gene_ids = (fused_gene_id,)
        details = {"intergenic_gap_bp": gap}
    elif event_type == COPY_COLLAPSE_EVENT:
        first, second = option.payload
        collapsed_id, retained_id = copy_collapse_partner_ids(
            first.gene_id, second.gene_id, seed=seed
        )
        by_id = {first.gene_id: first, second.gene_id: second}
        collapsed, retained = by_id[collapsed_id], by_id[retained_id]
        edits.extend(_line_edit(record, ()) for record in collapsed.hierarchy.removed_records)
        source_genes = (collapsed,)
        perturbed_gene_ids = ()
        details = {
            "retained_partner_gene_id": retained.gene_id,
            "pair_gene_ids": sorted((collapsed.gene_id, retained.gene_id)),
        }
    else:
        raise ValueError(f"Unsupported structure event type: {event_type}")

    edits.sort(key=lambda edit: int(edit["source_line_number"]))
    return {
        "event_id": event_id,
        "event_type": event_type,
        "target": _target(source_genes, perturbed_gene_ids),
        "details": details,
        "line_edits": edits,
    }


def generate_structure_benchmark(
    gff_path: str | Path,
    output_dir: str | Path,
    *,
    event_type: str,
    count: int,
    seed: int,
    truth_output_dir: str | Path | None = None,
    pair_tsv_path: str | Path | None = None,
) -> dict[str, Any]:
    """Create an exactly reversible non-missing-gene annotation benchmark."""

    if event_type not in STRUCTURE_EVENT_TYPES:
        raise ValueError(f"Unsupported structure event type: {event_type}")
    if count < 1:
        raise ValueError("count must be at least 1")
    source_path = Path(gff_path)
    output_dir = Path(output_dir)
    truth_dir = Path(truth_output_dir) if truth_output_dir is not None else output_dir
    output_gff = output_dir / "perturbed.gff3"
    truth_path = truth_dir / "hidden_truth.json"
    manifest_path = output_dir / "manifest.json"
    collisions = [path for path in (output_gff, truth_path, manifest_path) if path.exists()]
    if collisions:
        raise FileExistsError(
            "Refusing to overwrite benchmark output(s): "
            + ", ".join(str(path) for path in collisions)
        )

    document = read_gff_document(source_path)
    source_sha = document.text_sha256
    genes = _simple_genes(document)
    if event_type == MISSING_INTERNAL_EXON_EVENT:
        options = _missing_exon_options(genes)
    elif event_type == BOUNDARY_SHIFT_EVENT:
        options = _boundary_options(genes)
    elif event_type == SPLIT_GENE_EVENT:
        options = _split_options(genes)
    elif event_type == FUSED_GENE_EVENT:
        options = _fusion_options(genes, document)
    else:
        options = _copy_options(genes, pair_tsv_path)
    selected = _select_options(options, count, seed, event_type)
    events = [
        _event_from_option(
            option, event_type=event_type, source_sha=source_sha, seed=seed
        )
        for option in selected
    ]

    replacements_by_line: dict[int, tuple[str, ...]] = {}
    for event in events:
        for edit in event["line_edits"]:
            line_number = int(edit["source_line_number"])
            if line_number in replacements_by_line:
                raise AssertionError(f"Selected events overlap at source line {line_number}")
            replacements_by_line[line_number] = tuple(
                replacement["raw_line"] for replacement in edit["perturbed_lines"]
            )
    perturbed_lines: list[str] = []
    for line_number, raw_line in enumerate(document.lines, start=1):
        perturbed_lines.extend(replacements_by_line.get(line_number, (raw_line,)))

    pair_source = None
    if pair_tsv_path is not None:
        pair_source = {
            "file_name": Path(pair_tsv_path).name,
            "sha256": _file_sha256(pair_tsv_path),
        }
    source = {
        "file_name": source_path.name,
        "file_sha256": _file_sha256(source_path),
        "text_sha256": source_sha,
        "line_count": len(document.lines),
    }
    truth: dict[str, Any] = {
        "schema_version": TRUTH_SCHEMA_VERSION,
        "generator": {"name": "PloidyPatch", "version": __version__},
        "source": source,
        "perturbation": {
            "event_type": event_type,
            "seed": seed,
            "selection_algorithm": "sha256_rank_nonoverlap_v1",
            "pair_source": pair_source,
            "perturbed_text_sha256": _text_sha256(perturbed_lines),
        },
        "events": events,
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    truth_dir.mkdir(parents=True, exist_ok=True)
    _write_text_exclusive(output_gff, "".join(perturbed_lines))
    _write_text_exclusive(truth_path, _json_text(truth))
    manifest: dict[str, Any] = {
        "schema_version": MANIFEST_SCHEMA_VERSION,
        "generator": {"name": "PloidyPatch", "version": __version__},
        "source": source,
        "perturbation": {
            "event_type": event_type,
            "seed": seed,
            "selection_algorithm": "sha256_rank_nonoverlap_v1",
            "pair_source": pair_source,
            "requested_events": count,
            "eligible_options": len(options),
            "generated_events": len(events),
        },
        "outputs": {
            "perturbed_gff": {
                "file_name": output_gff.name,
                "sha256": _file_sha256(output_gff),
            },
            "hidden_truth": {
                "file_name": truth_path.name,
                "sha256": _file_sha256(truth_path),
                "storage": (
                    "separate_evaluator_directory"
                    if truth_dir.resolve() != output_dir.resolve()
                    else "benchmark_output_directory"
                ),
            },
        },
    }
    _write_text_exclusive(manifest_path, _json_text(manifest))
    return manifest
