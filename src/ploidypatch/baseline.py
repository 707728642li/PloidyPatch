from __future__ import annotations

import csv
import bisect
import hashlib
import json
import re
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any
from urllib.parse import quote

from . import __version__
from .audit import PROTEIN_ALPHABET
from .gff import TRANSCRIPT_TYPES, parse_attributes
from .io import iter_fasta, open_text
from .perturb import GffRecord, read_gff_document


REFERENCE_PROTEIN_MANIFEST_SCHEMA_VERSION = (
    "ploidypatch.reference_protein_manifest.v1"
)
MINIPROT_ADAPTER_MANIFEST_SCHEMA_VERSION = (
    "ploidypatch.miniprot_adapter_manifest.v2"
)
PROJECTION_SUPPORT_SCHEMA_VERSION = "ploidypatch.projection_support.v1"
GFF_ADAPTER_MANIFEST_SCHEMA_VERSION = "ploidypatch.gff_adapter_manifest.v1"
SAFE_SOURCE = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]*")


def _file_sha256(path: str | Path, block_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        while block := handle.read(block_size):
            digest.update(block)
    return digest.hexdigest()


def _parse_named_paths(values: list[str]) -> list[tuple[str, Path]]:
    parsed: list[tuple[str, Path]] = []
    seen: set[str] = set()
    for value in values:
        source, separator, path_text = value.partition("=")
        if not separator or not SAFE_SOURCE.fullmatch(source) or not path_text:
            raise ValueError(
                f"Invalid SOURCE=PATH protein input: {value!r}"
            )
        if source in seen:
            raise ValueError(f"Duplicate protein source: {source}")
        seen.add(source)
        parsed.append((source, Path(path_text)))
    if not parsed:
        raise ValueError("At least one SOURCE=PATH protein input is required")
    return parsed


def prepare_reference_proteins(
    *,
    protein_inputs: list[str],
    output_fasta_path: str | Path,
    output_map_path: str | Path,
) -> dict[str, Any]:
    """Prefix and combine external proteins with complete source provenance."""

    inputs = _parse_named_paths(protein_inputs)
    output_fasta = Path(output_fasta_path)
    output_map = Path(output_map_path)
    manifest_path = Path(str(output_fasta) + ".manifest.json")
    collisions = [
        path for path in (output_fasta, output_map, manifest_path) if path.exists()
    ]
    if collisions:
        raise FileExistsError(
            "Refusing to overwrite reference protein artifact(s): "
            + ", ".join(str(path) for path in collisions)
        )
    if output_fasta.resolve() == output_map.resolve():
        raise ValueError("Output FASTA and map paths must differ")

    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    output_map.parent.mkdir(parents=True, exist_ok=True)
    source_reports: dict[str, dict[str, Any]] = {}
    query_ids: set[str] = set()
    total_records = 0
    with output_fasta.open("x", encoding="utf-8", newline="") as fasta_handle, (
        output_map.open("x", encoding="utf-8", newline="")
    ) as map_handle:
        writer = csv.writer(map_handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            ("query_id", "source", "source_record_id", "length_aa", "source_header")
        )
        for source, path in inputs:
            source_count = 0
            for record_id, header, sequence in iter_fasta(path):
                query_id = f"{source}__{record_id}"
                if query_id in query_ids:
                    raise ValueError(f"Duplicate prefixed protein ID: {query_id}")
                query_ids.add(query_id)
                protein = sequence.upper()
                if not protein:
                    raise ValueError(f"Empty protein record: {source}:{record_id}")
                invalid = set(protein) - PROTEIN_ALPHABET
                if invalid:
                    raise ValueError(
                        f"Invalid protein character(s) in {source}:{record_id}: "
                        + ", ".join(sorted(invalid))
                    )
                fasta_handle.write(f">{query_id}\n")
                for index in range(0, len(protein), 70):
                    fasta_handle.write(protein[index : index + 70] + "\n")
                writer.writerow((query_id, source, record_id, len(protein), header))
                source_count += 1
                total_records += 1
            source_reports[source] = {
                "file_name": path.name,
                "bytes": path.stat().st_size,
                "sha256": _file_sha256(path),
                "records": source_count,
            }

    manifest: dict[str, Any] = {
        "schema_version": REFERENCE_PROTEIN_MANIFEST_SCHEMA_VERSION,
        "generator": {"name": "PloidyPatch", "version": __version__},
        "inputs": source_reports,
        "outputs": {
            "fasta": {
                "file_name": output_fasta.name,
                "bytes": output_fasta.stat().st_size,
                "sha256": _file_sha256(output_fasta),
                "records": total_records,
            },
            "map": {
                "file_name": output_map.name,
                "bytes": output_map.stat().st_size,
                "sha256": _file_sha256(output_map),
                "rows": total_records,
            },
        },
    }
    with manifest_path.open("x", encoding="utf-8", newline="") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return manifest


@dataclass(frozen=True)
class ProjectedModel:
    model_id: str
    query_id: str
    source: str
    seqid: str
    strand: str
    start: int
    end: int
    score: float
    rank: int
    identity: float
    positive: float
    query_coverage: float
    frameshifts: int
    stop_codons: int
    cds: tuple[tuple[int, int, str], ...]

    @property
    def cds_bp(self) -> int:
        return sum(end - start + 1 for start, end, _ in self.cds)


@dataclass(frozen=True)
class MergedIntervalIndex:
    intervals: tuple[tuple[int, int], ...]
    ends: tuple[int, ...]


@dataclass(frozen=True)
class AnnotationCandidateModel:
    model_id: str
    seqid: str
    strand: str
    start: int
    end: int
    score: float | None
    score_source: str
    exons: tuple[tuple[int, int], ...]
    cds: tuple[tuple[int, int, str], ...]
    phase_normalization: str
    upstream_attributes: str

    @property
    def cds_bp(self) -> int:
        return sum(end - start + 1 for start, end in _model_cds_intervals(self))


def _read_protein_map(path: str | Path) -> dict[str, tuple[str, int]]:
    rows: dict[str, tuple[str, int]] = {}
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"query_id", "source", "length_aa"}
        if reader.fieldnames is None or not required <= set(reader.fieldnames):
            raise ValueError("Protein map lacks query_id, source, or length_aa")
        for line_number, row in enumerate(reader, start=2):
            query_id = row["query_id"]
            if not query_id or query_id in rows:
                raise ValueError(
                    f"Empty or duplicate query_id in protein map line {line_number}"
                )
            try:
                length = int(row["length_aa"])
            except ValueError as exc:
                raise ValueError(
                    f"Invalid protein length in map line {line_number}"
                ) from exc
            if length < 1:
                raise ValueError(f"Non-positive protein length in map line {line_number}")
            rows[query_id] = (row["source"], length)
    return rows


def _parse_target(value: str, line_number: int) -> tuple[str, int, int]:
    fields = value.split()
    if len(fields) != 3:
        raise ValueError(f"Malformed miniprot Target at line {line_number}: {value}")
    try:
        start = int(fields[1])
        end = int(fields[2])
    except ValueError as exc:
        raise ValueError(f"Malformed miniprot Target at line {line_number}: {value}") from exc
    if start < 1 or end < start:
        raise ValueError(f"Invalid miniprot Target interval at line {line_number}")
    return fields[0], start, end


def _parse_miniprot_models(
    path: str | Path,
    protein_map: dict[str, tuple[str, int]],
) -> list[ProjectedModel]:
    models: list[ProjectedModel] = []
    model_ids: set[str] = set()
    current: dict[str, Any] | None = None

    def finish() -> None:
        nonlocal current
        if current is None:
            return
        if current["model_id"] in model_ids:
            raise ValueError(f"Duplicate miniprot model ID: {current['model_id']}")
        model_ids.add(current["model_id"])
        models.append(
            ProjectedModel(
                model_id=current["model_id"],
                query_id=current["query_id"],
                source=current["source"],
                seqid=current["seqid"],
                strand=current["strand"],
                start=current["start"],
                end=current["end"],
                score=current["score"],
                rank=current["rank"],
                identity=current["identity"],
                positive=current["positive"],
                query_coverage=current["query_coverage"],
                frameshifts=current["frameshifts"],
                stop_codons=current["stop_codons"],
                cds=tuple(sorted(set(current["cds"]))),
            )
        )
        current = None

    with open_text(path) as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\r\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) != 9:
                raise ValueError(f"Malformed miniprot GFF3 line {line_number}")
            attributes, malformed = parse_attributes(fields[8])
            if malformed:
                raise ValueError(
                    f"Malformed miniprot attributes at line {line_number}"
                )
            feature_type = fields[2]
            if feature_type == "mRNA":
                finish()
                model_id = attributes.get("ID")
                target = attributes.get("Target")
                if not model_id or not target:
                    raise ValueError(
                        f"Miniprot mRNA lacks ID or Target at line {line_number}"
                    )
                query_id, query_start, query_end = _parse_target(
                    target, line_number
                )
                query_metadata = protein_map.get(query_id)
                if query_metadata is None:
                    raise ValueError(
                        f"Miniprot query is absent from protein map: {query_id}"
                    )
                source, query_length = query_metadata
                current = {
                    "model_id": model_id,
                    "query_id": query_id,
                    "source": source,
                    "seqid": fields[0],
                    "strand": fields[6],
                    "start": int(fields[3]),
                    "end": int(fields[4]),
                    "score": float(fields[5]) if fields[5] != "." else 0.0,
                    "rank": int(attributes.get("Rank", "1")),
                    "identity": float(attributes.get("Identity", "0")),
                    "positive": float(attributes.get("Positive", "0")),
                    "query_coverage": (query_end - query_start + 1) / query_length,
                    "frameshifts": int(attributes.get("Frameshift", "0")),
                    "stop_codons": int(attributes.get("StopCodon", "0")),
                    "cds": [],
                }
                continue
            if feature_type not in {"CDS", "stop_codon"}:
                raise ValueError(
                    f"Unsupported miniprot feature type at line {line_number}: "
                    f"{feature_type}"
                )
            if current is None:
                raise ValueError(
                    f"Miniprot child feature precedes mRNA at line {line_number}"
                )
            if attributes.get("Parent") != current["model_id"]:
                raise ValueError(
                    f"Miniprot child Parent mismatch at line {line_number}"
                )
            if feature_type == "CDS":
                current["cds"].append((int(fields[3]), int(fields[4]), fields[7]))
    finish()
    return models


def _merge_intervals(
    intervals: list[tuple[int, int]],
) -> list[tuple[int, int]]:
    merged: list[list[int]] = []
    for start, end in sorted(intervals):
        if merged and start <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return [(start, end) for start, end in merged]


def _existing_cds_intervals(
    gff_path: str | Path,
) -> dict[tuple[str, str], MergedIntervalIndex]:
    raw: dict[tuple[str, str], list[tuple[int, int]]] = defaultdict(list)
    with open_text(gff_path) as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\r\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) != 9:
                raise ValueError(f"Malformed blind GFF3 line {line_number}")
            if fields[2] != "CDS":
                continue
            raw[(fields[0], fields[6])].append((int(fields[3]), int(fields[4])))
    indexes: dict[tuple[str, str], MergedIntervalIndex] = {}
    for key, values in raw.items():
        intervals = tuple(_merge_intervals(values))
        indexes[key] = MergedIntervalIndex(
            intervals=intervals,
            ends=tuple(end for _, end in intervals),
        )
    return indexes


def _overlap_with_index(
    query: list[tuple[int, int]],
    reference: MergedIntervalIndex | None,
) -> int:
    if reference is None:
        return 0
    total = 0
    for query_start, query_end in query:
        index = bisect.bisect_left(reference.ends, query_start)
        while index < len(reference.intervals):
            reference_start, reference_end = reference.intervals[index]
            if reference_start > query_end:
                break
            start = max(query_start, reference_start)
            end = min(query_end, reference_end)
            if start <= end:
                total += end - start + 1
            index += 1
    return total


def _overlap_bp(
    left: list[tuple[int, int]],
    right: list[tuple[int, int]],
) -> int:
    total = 0
    left_index = 0
    right_index = 0
    while left_index < len(left) and right_index < len(right):
        left_start, left_end = left[left_index]
        right_start, right_end = right[right_index]
        start = max(left_start, right_start)
        end = min(left_end, right_end)
        if start <= end:
            total += end - start + 1
        if left_end <= right_end:
            left_index += 1
        else:
            right_index += 1
    return total


def _model_intervals(model: ProjectedModel) -> list[tuple[int, int]]:
    return _merge_intervals([(start, end) for start, end, _ in model.cds])


def _model_cds_intervals(
    model: AnnotationCandidateModel,
) -> list[tuple[int, int]]:
    return _merge_intervals([(start, end) for start, end, _ in model.cds])


def _annotation_model_score(
    record: GffRecord,
    score_attribute: str | None,
) -> tuple[float | None, str]:
    if score_attribute is not None and score_attribute in record.attributes:
        raw_score = record.attributes[score_attribute]
        try:
            return float(raw_score), f"attribute:{score_attribute}"
        except ValueError as exc:
            raise ValueError(
                f"Non-numeric {score_attribute!r} score for candidate "
                f"{record.feature_id!r}"
            ) from exc
    fields = record.raw_line.rstrip("\r\n").split("\t")
    if fields[5] == ".":
        return None, "none"
    try:
        return float(fields[5]), "gff_score_column"
    except ValueError as exc:
        raise ValueError(
            f"Non-numeric GFF score for candidate {record.feature_id!r}"
        ) from exc


def _parse_annotation_candidate_models(
    path: str | Path,
    score_attribute: str | None,
    infer_missing_cds_phase: bool = False,
) -> list[AnnotationCandidateModel]:
    document = read_gff_document(path)
    transcript_records: list[GffRecord] = []
    seen_transcripts: set[str] = set()
    children: dict[str, list[GffRecord]] = defaultdict(list)
    for record in document.records:
        if record.feature_type in TRANSCRIPT_TYPES:
            if not record.feature_id:
                raise ValueError(
                    f"Candidate transcript lacks ID at line {record.line_number}"
                )
            if record.feature_id in seen_transcripts:
                raise ValueError(
                    f"Duplicate candidate transcript ID: {record.feature_id}"
                )
            seen_transcripts.add(record.feature_id)
            transcript_records.append(record)
        if record.feature_type in {"exon", "CDS"}:
            for parent in record.parents:
                children[parent].append(record)

    models: list[AnnotationCandidateModel] = []
    for record in transcript_records:
        assert record.feature_id is not None
        model_children = children.get(record.feature_id, [])
        for child in model_children:
            if child.seqid != record.seqid or child.strand != record.strand:
                raise ValueError(
                    "Candidate child disagrees with transcript locus: "
                    f"{record.feature_id} line {child.line_number}"
                )
        exons = tuple(
            _merge_intervals(
                [
                    (child.start, child.end)
                    for child in model_children
                    if child.feature_type == "exon"
                ]
            )
        )
        raw_cds = {
            (child.start, child.end, child.phase)
            for child in model_children
            if child.feature_type == "CDS"
        }
        phase_normalization = "upstream"
        if (
            infer_missing_cds_phase
            and raw_cds
            and all(phase == "." for _, _, phase in raw_cds)
            and record.strand in {"+", "-"}
        ):
            ordered_cds = sorted(
                raw_cds,
                key=lambda value: (value[0], value[1]),
                reverse=record.strand == "-",
            )
            cumulative_length = 0
            inferred_cds = []
            for start, end, _ in ordered_cds:
                phase = str((3 - cumulative_length % 3) % 3)
                inferred_cds.append((start, end, phase))
                cumulative_length += end - start + 1
            raw_cds = set(inferred_cds)
            phase_normalization = "inferred_full_cds_first_phase_zero"
        cds = tuple(sorted(raw_cds))
        score, score_source = _annotation_model_score(record, score_attribute)
        models.append(
            AnnotationCandidateModel(
                model_id=record.feature_id,
                seqid=record.seqid,
                strand=record.strand,
                start=record.start,
                end=record.end,
                score=score,
                score_source=score_source,
                exons=exons,
                cds=cds,
                phase_normalization=phase_normalization,
                upstream_attributes=record.raw_line.rstrip("\r\n").split("\t")[8],
            )
        )
    return models


def _write_appended_annotation_model(
    handle: Any,
    model: AnnotationCandidateModel,
    number: int,
    source: str,
) -> None:
    gene_id = f"PPGFF_gene_{number:06d}"
    transcript_id = f"PPGFF_tx_{number:06d}"
    score = "." if model.score is None else f"{model.score:g}"
    attributes = (
        f"baseline_source={quote(source, safe='._:-')};"
        f"upstream_model={quote(model.model_id, safe='._:-')};"
        f"upstream_score_source={quote(model.score_source, safe='._:-')}"
    )
    prefix = f"{model.seqid}\tPloidyPatchBaseline"
    handle.write(
        f"{prefix}\tgene\t{model.start}\t{model.end}\t{score}\t"
        f"{model.strand}\t.\tID={gene_id};{attributes}\n"
    )
    handle.write(
        f"{prefix}\tmRNA\t{model.start}\t{model.end}\t{score}\t"
        f"{model.strand}\t.\tID={transcript_id};Parent={gene_id};{attributes}\n"
    )
    exon_intervals = list(model.exons) or _model_cds_intervals(model)
    for exon_number, (start, end) in enumerate(exon_intervals, start=1):
        handle.write(
            f"{prefix}\texon\t{start}\t{end}\t.\t{model.strand}\t.\t"
            f"ID=PPGFF_exon_{number:06d}_{exon_number};Parent={transcript_id}\n"
        )
    for start, end, phase in model.cds:
        handle.write(
            f"{prefix}\tCDS\t{start}\t{end}\t.\t{model.strand}\t{phase}\t"
            f"Parent={transcript_id}\n"
        )


def adapt_annotation_gff_baseline(
    *,
    perturbed_gff_path: str | Path,
    candidate_gff_path: str | Path,
    source: str,
    output_gff_path: str | Path,
    decisions_tsv_path: str | Path,
    score_attribute: str | None = None,
    max_existing_cds_overlap: float = 0.2,
    max_redundancy_overlap: float = 0.5,
    infer_missing_cds_phase: bool = False,
) -> dict[str, Any]:
    """Append coding models from a general GFF3 as gap-only hypotheses."""

    if not SAFE_SOURCE.fullmatch(source):
        raise ValueError(
            "source must start with an alphanumeric character and contain only "
            "letters, digits, dots, underscores, or hyphens"
        )
    for name, value in (
        ("max_existing_cds_overlap", max_existing_cds_overlap),
        ("max_redundancy_overlap", max_redundancy_overlap),
    ):
        if not 0 <= value <= 1:
            raise ValueError(f"{name} must be between zero and one")
    output_gff = Path(output_gff_path)
    decisions_path = Path(decisions_tsv_path)
    manifest_path = Path(str(output_gff) + ".manifest.json")
    collisions = [
        path for path in (output_gff, decisions_path, manifest_path) if path.exists()
    ]
    if collisions:
        raise FileExistsError(
            "Refusing to overwrite GFF adapter artifact(s): "
            + ", ".join(str(path) for path in collisions)
        )

    models = _parse_annotation_candidate_models(
        candidate_gff_path,
        score_attribute,
        infer_missing_cds_phase=infer_missing_cds_phase,
    )
    blind_document = read_gff_document(perturbed_gff_path)
    blind_seqids = {record.seqid for record in blind_document.records}
    blind_ids = {
        record.feature_id for record in blind_document.records if record.feature_id
    }
    existing = _existing_cds_intervals(perturbed_gff_path)
    decisions: dict[str, tuple[str, str, float, str]] = {}
    eligible: list[AnnotationCandidateModel] = []
    eligible_existing_overlap: dict[str, float] = {}
    for model in models:
        reason = ""
        existing_overlap = 0.0
        raw_cds_bp = sum(end - start + 1 for start, end, _ in model.cds)
        child_intervals = list(model.exons) + [
            (start, end) for start, end, _ in model.cds
        ]
        if not model.cds:
            reason = "no_cds"
        elif model.seqid not in blind_seqids:
            reason = "unknown_target_seqid"
        elif model.strand not in {"+", "-"}:
            reason = "invalid_strand"
        elif any(phase not in {"0", "1", "2"} for _, _, phase in model.cds):
            reason = "invalid_cds_phase"
        elif any(
            start < model.start or end > model.end for start, end in child_intervals
        ):
            reason = "child_outside_transcript"
        elif raw_cds_bp != model.cds_bp:
            reason = "overlapping_cds_segments"
        else:
            intervals = _model_cds_intervals(model)
            existing_overlap = _overlap_with_index(
                intervals, existing.get((model.seqid, model.strand))
            ) / model.cds_bp
            if existing_overlap > max_existing_cds_overlap:
                reason = "overlaps_blind_cds"
        if reason:
            decisions[model.model_id] = (
                "rejected",
                reason,
                existing_overlap,
                "",
            )
        else:
            eligible.append(model)
            eligible_existing_overlap[model.model_id] = existing_overlap

    accepted: list[AnnotationCandidateModel] = []
    redundancy_bin_size = 100_000
    accepted_bins: dict[tuple[str, str, int], list[int]] = defaultdict(list)
    ranked = sorted(
        eligible,
        key=lambda model: (
            model.score is None,
            -(model.score or 0.0),
            -model.cds_bp,
            model.seqid,
            model.start,
            model.model_id,
        ),
    )
    for model in ranked:
        intervals = _model_cds_intervals(model)
        nearby_indices: set[int] = set()
        for bin_number in range(
            model.start // redundancy_bin_size,
            model.end // redundancy_bin_size + 1,
        ):
            nearby_indices.update(
                accepted_bins.get((model.seqid, model.strand, bin_number), [])
            )
        redundant_matches: list[tuple[float, float, str]] = []
        for accepted_index in nearby_indices:
            other = accepted[accepted_index]
            shared = _overlap_bp(intervals, _model_cds_intervals(other))
            denominator = min(model.cds_bp, other.cds_bp)
            if denominator and shared / denominator > max_redundancy_overlap:
                redundant_matches.append(
                    (shared / denominator, other.score or 0.0, other.model_id)
                )
        if redundant_matches:
            _, _, redundant_with = max(redundant_matches)
            decisions[model.model_id] = (
                "rejected",
                "redundant_candidate",
                eligible_existing_overlap[model.model_id],
                redundant_with,
            )
        else:
            accepted_index = len(accepted)
            accepted.append(model)
            for bin_number in range(
                model.start // redundancy_bin_size,
                model.end // redundancy_bin_size + 1,
            ):
                accepted_bins[(model.seqid, model.strand, bin_number)].append(
                    accepted_index
                )
            decisions[model.model_id] = (
                "accepted",
                "accepted",
                eligible_existing_overlap[model.model_id],
                "",
            )

    generated_ids: set[str] = set()
    for number, model in enumerate(accepted, start=1):
        generated_ids.update(
            {
                f"PPGFF_gene_{number:06d}",
                f"PPGFF_tx_{number:06d}",
            }
        )
        exon_count = len(model.exons) or len(_model_cds_intervals(model))
        generated_ids.update(
            f"PPGFF_exon_{number:06d}_{exon_number}"
            for exon_number in range(1, exon_count + 1)
        )
    collisions_with_blind = generated_ids & blind_ids
    if collisions_with_blind:
        raise ValueError(
            "Generated adapter IDs collide with blind annotation: "
            + ", ".join(sorted(collisions_with_blind)[:10])
        )

    output_gff.parent.mkdir(parents=True, exist_ok=True)
    with open_text(perturbed_gff_path) as source_handle, output_gff.open(
        "x", encoding="utf-8", newline=""
    ) as output_handle:
        last_line = ""
        for line in source_handle:
            output_handle.write(line)
            last_line = line
        if last_line and not last_line.endswith(("\n", "\r")):
            output_handle.write("\n")
        output_handle.write("###\n")
        for number, model in enumerate(accepted, start=1):
            _write_appended_annotation_model(output_handle, model, number, source)

    decisions_path.parent.mkdir(parents=True, exist_ok=True)
    reason_counts: Counter[str] = Counter()
    with decisions_path.open("x", encoding="utf-8", newline="") as handle:
        fieldnames = (
            "model_id",
            "source",
            "seqid",
            "start",
            "end",
            "strand",
            "score",
            "score_source",
            "exon_segments",
            "cds_segments",
            "cds_bp",
            "existing_cds_overlap_fraction",
            "status",
            "reason",
            "redundant_with_model_id",
            "upstream_attributes",
            "phase_normalization",
        )
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for model in models:
            status, reason, existing_overlap, redundant_with = decisions[
                model.model_id
            ]
            reason_counts[reason] += 1
            writer.writerow(
                {
                    "model_id": model.model_id,
                    "source": source,
                    "seqid": model.seqid,
                    "start": model.start,
                    "end": model.end,
                    "strand": model.strand,
                    "score": "" if model.score is None else f"{model.score:g}",
                    "score_source": model.score_source,
                    "exon_segments": len(model.exons),
                    "cds_segments": len(model.cds),
                    "cds_bp": model.cds_bp,
                    "existing_cds_overlap_fraction": f"{existing_overlap:.6f}",
                    "status": status,
                    "reason": reason,
                    "redundant_with_model_id": redundant_with,
                    "upstream_attributes": model.upstream_attributes,
                    "phase_normalization": model.phase_normalization,
                }
            )

    manifest: dict[str, Any] = {
        "schema_version": GFF_ADAPTER_MANIFEST_SCHEMA_VERSION,
        "generator": {"name": "PloidyPatch", "version": __version__},
        "inputs": {
            "perturbed_gff": {
                "file_name": Path(perturbed_gff_path).name,
                "sha256": _file_sha256(perturbed_gff_path),
            },
            "candidate_gff": {
                "file_name": Path(candidate_gff_path).name,
                "sha256": _file_sha256(candidate_gff_path),
            },
        },
        "parameters": {
            "source": source,
            "score_attribute": score_attribute,
            "max_existing_cds_overlap": max_existing_cds_overlap,
            "max_redundancy_overlap": max_redundancy_overlap,
            "infer_missing_cds_phase": infer_missing_cds_phase,
        },
        "normalization": {
            "phase_normalization_counts": dict(
                sorted(Counter(model.phase_normalization for model in models).items())
            ),
            "inference_contract": (
                "all_phases_missing_full_cds_first_phase_zero"
                if infer_missing_cds_phase
                else "disabled"
            ),
        },
        "models": {
            "input": len(models),
            "accepted": len(accepted),
            "rejected": len(models) - len(accepted),
            "decision_counts": dict(sorted(reason_counts.items())),
        },
        "outputs": {
            "candidate_gff": {
                "file_name": output_gff.name,
                "sha256": _file_sha256(output_gff),
            },
            "decisions": {
                "file_name": decisions_path.name,
                "sha256": _file_sha256(decisions_path),
                "rows": len(models),
            },
        },
    }
    with manifest_path.open("x", encoding="utf-8", newline="") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return manifest


def _write_appended_model(
    handle: Any,
    model: ProjectedModel,
    number: int,
) -> None:
    gene_id = f"PPBASE_gene_{number:06d}"
    transcript_id = f"PPBASE_tx_{number:06d}"
    attributes = (
        f"Derives_from={model.query_id};miniprot_model={model.model_id};"
        f"reference_source={model.source};Identity={model.identity:.4f};"
        f"QueryCoverage={model.query_coverage:.4f}"
    )
    prefix = f"{model.seqid}\tPloidyPatchBaseline"
    handle.write(
        f"{prefix}\tgene\t{model.start}\t{model.end}\t{model.score:g}\t"
        f"{model.strand}\t.\tID={gene_id};{attributes}\n"
    )
    handle.write(
        f"{prefix}\tmRNA\t{model.start}\t{model.end}\t{model.score:g}\t"
        f"{model.strand}\t.\tID={transcript_id};Parent={gene_id};{attributes}\n"
    )
    for exon_number, (start, end) in enumerate(_model_intervals(model), start=1):
        handle.write(
            f"{prefix}\texon\t{start}\t{end}\t.\t{model.strand}\t.\t"
            f"ID=PPBASE_exon_{number:06d}_{exon_number};Parent={transcript_id}\n"
        )
    for start, end, phase in model.cds:
        handle.write(
            f"{prefix}\tCDS\t{start}\t{end}\t.\t{model.strand}\t{phase}\t"
            f"Parent={transcript_id}\n"
        )


def adapt_miniprot_baseline(
    *,
    perturbed_gff_path: str | Path,
    miniprot_gff_path: str | Path,
    protein_map_path: str | Path,
    output_gff_path: str | Path,
    decisions_tsv_path: str | Path,
    min_identity: float = 0.5,
    min_query_coverage: float = 0.5,
    max_existing_cds_overlap: float = 0.2,
    max_redundancy_overlap: float = 0.5,
    require_intact: bool = True,
) -> dict[str, Any]:
    """Append conservative, annotation-gap miniprot models to a blind GFF3."""

    for name, value in (
        ("min_identity", min_identity),
        ("min_query_coverage", min_query_coverage),
        ("max_existing_cds_overlap", max_existing_cds_overlap),
        ("max_redundancy_overlap", max_redundancy_overlap),
    ):
        if not 0 <= value <= 1:
            raise ValueError(f"{name} must be between zero and one")
    output_gff = Path(output_gff_path)
    decisions_path = Path(decisions_tsv_path)
    manifest_path = Path(str(output_gff) + ".manifest.json")
    collisions = [
        path for path in (output_gff, decisions_path, manifest_path) if path.exists()
    ]
    if collisions:
        raise FileExistsError(
            "Refusing to overwrite miniprot adapter artifact(s): "
            + ", ".join(str(path) for path in collisions)
        )

    protein_map = _read_protein_map(protein_map_path)
    models = _parse_miniprot_models(miniprot_gff_path, protein_map)
    existing = _existing_cds_intervals(perturbed_gff_path)
    decisions: dict[str, tuple[str, str, float, str]] = {}
    eligible: list[ProjectedModel] = []
    eligible_existing_overlap: dict[str, float] = {}
    for model in models:
        reason = ""
        existing_overlap = 0.0
        if not model.cds:
            reason = "no_cds"
        elif model.identity < min_identity:
            reason = "identity_below_threshold"
        elif model.query_coverage < min_query_coverage:
            reason = "query_coverage_below_threshold"
        elif require_intact and (model.frameshifts > 0 or model.stop_codons > 0):
            reason = "non_intact_alignment"
        else:
            intervals = _model_intervals(model)
            existing_overlap = _overlap_with_index(
                intervals, existing.get((model.seqid, model.strand))
            ) / model.cds_bp
            if existing_overlap > max_existing_cds_overlap:
                reason = "overlaps_blind_cds"
        if reason:
            decisions[model.model_id] = (
                "rejected",
                reason,
                existing_overlap,
                "",
            )
        else:
            eligible.append(model)
            eligible_existing_overlap[model.model_id] = existing_overlap

    accepted: list[ProjectedModel] = []
    redundancy_bin_size = 100_000
    accepted_bins: dict[tuple[str, str, int], list[int]] = defaultdict(list)
    ranked = sorted(
        eligible,
        key=lambda model: (
            -model.score,
            -model.identity,
            -model.query_coverage,
            model.rank,
            model.model_id,
        ),
    )
    for model in ranked:
        intervals = _model_intervals(model)
        redundant_matches: list[tuple[float, float, str]] = []
        nearby_indices: set[int] = set()
        for bin_number in range(
            model.start // redundancy_bin_size,
            model.end // redundancy_bin_size + 1,
        ):
            nearby_indices.update(
                accepted_bins.get((model.seqid, model.strand, bin_number), [])
            )
        for accepted_index in nearby_indices:
            other = accepted[accepted_index]
            shared = _overlap_bp(intervals, _model_intervals(other))
            denominator = min(model.cds_bp, other.cds_bp)
            if denominator and shared / denominator > max_redundancy_overlap:
                redundant_matches.append(
                    (shared / denominator, other.score, other.model_id)
                )
        if redundant_matches:
            _, _, redundant_with = max(redundant_matches)
            decisions[model.model_id] = (
                "rejected",
                "redundant_projection",
                eligible_existing_overlap[model.model_id],
                redundant_with,
            )
        else:
            accepted_index = len(accepted)
            accepted.append(model)
            for bin_number in range(
                model.start // redundancy_bin_size,
                model.end // redundancy_bin_size + 1,
            ):
                accepted_bins[(model.seqid, model.strand, bin_number)].append(
                    accepted_index
                )
            decisions[model.model_id] = (
                "accepted",
                "accepted",
                eligible_existing_overlap[model.model_id],
                "",
            )

    output_gff.parent.mkdir(parents=True, exist_ok=True)
    with open_text(perturbed_gff_path) as source_handle, output_gff.open(
        "x", encoding="utf-8", newline=""
    ) as output_handle:
        last_line = ""
        for line in source_handle:
            output_handle.write(line)
            last_line = line
        if last_line and not last_line.endswith(("\n", "\r")):
            output_handle.write("\n")
        output_handle.write("###\n")
        for number, model in enumerate(accepted, start=1):
            _write_appended_model(output_handle, model, number)

    decisions_path.parent.mkdir(parents=True, exist_ok=True)
    reason_counts: Counter[str] = Counter()
    with decisions_path.open("x", encoding="utf-8", newline="") as handle:
        fieldnames = (
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
            "positive",
            "query_coverage",
            "frameshifts",
            "stop_codons",
            "cds_segments",
            "cds_bp",
            "existing_cds_overlap_fraction",
            "status",
            "reason",
            "redundant_with_model_id",
        )
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for model in models:
            status, reason, existing_overlap, redundant_with = decisions[
                model.model_id
            ]
            reason_counts[reason] += 1
            writer.writerow(
                {
                    "model_id": model.model_id,
                    "query_id": model.query_id,
                    "source": model.source,
                    "seqid": model.seqid,
                    "start": model.start,
                    "end": model.end,
                    "strand": model.strand,
                    "score": f"{model.score:g}",
                    "rank": model.rank,
                    "identity": f"{model.identity:.6f}",
                    "positive": f"{model.positive:.6f}",
                    "query_coverage": f"{model.query_coverage:.6f}",
                    "frameshifts": model.frameshifts,
                    "stop_codons": model.stop_codons,
                    "cds_segments": len(model.cds),
                    "cds_bp": model.cds_bp,
                    "existing_cds_overlap_fraction": f"{existing_overlap:.6f}",
                    "status": status,
                    "reason": reason,
                    "redundant_with_model_id": redundant_with,
                }
            )

    manifest: dict[str, Any] = {
        "schema_version": MINIPROT_ADAPTER_MANIFEST_SCHEMA_VERSION,
        "generator": {"name": "PloidyPatch", "version": __version__},
        "inputs": {
            "perturbed_gff": {
                "file_name": Path(perturbed_gff_path).name,
                "sha256": _file_sha256(perturbed_gff_path),
            },
            "miniprot_gff": {
                "file_name": Path(miniprot_gff_path).name,
                "sha256": _file_sha256(miniprot_gff_path),
            },
            "protein_map": {
                "file_name": Path(protein_map_path).name,
                "sha256": _file_sha256(protein_map_path),
                "rows": len(protein_map),
            },
        },
        "parameters": {
            "min_identity": min_identity,
            "min_query_coverage": min_query_coverage,
            "max_existing_cds_overlap": max_existing_cds_overlap,
            "max_redundancy_overlap": max_redundancy_overlap,
            "require_intact": require_intact,
        },
        "models": {
            "input": len(models),
            "accepted": len(accepted),
            "rejected": len(models) - len(accepted),
            "decision_counts": dict(sorted(reason_counts.items())),
        },
        "outputs": {
            "candidate_gff": {
                "file_name": output_gff.name,
                "sha256": _file_sha256(output_gff),
            },
            "decisions": {
                "file_name": decisions_path.name,
                "sha256": _file_sha256(decisions_path),
                "rows": len(models),
            },
        },
    }
    with manifest_path.open("x", encoding="utf-8", newline="") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return manifest


def summarize_projection_support(
    *,
    decisions_tsv_path: str | Path,
    output_tsv_path: str | Path,
) -> dict[str, Any]:
    """Recover multi-query and multi-progenitor support hidden by deduplication."""

    output = Path(output_tsv_path)
    manifest_path = Path(str(output) + ".manifest.json")
    collisions = [path for path in (output, manifest_path) if path.exists()]
    if collisions:
        raise FileExistsError(
            "Refusing to overwrite projection-support artifact(s): "
            + ", ".join(str(path) for path in collisions)
        )
    required = {
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
        "existing_cds_overlap_fraction",
        "status",
        "reason",
        "redundant_with_model_id",
    }
    rows: list[dict[str, str]] = []
    by_model: dict[str, dict[str, str]] = {}
    with Path(decisions_tsv_path).open(
        "r", encoding="utf-8", newline=""
    ) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = required - set(reader.fieldnames or [])
        if reader.fieldnames is None or missing:
            raise ValueError(
                "Decision TSV is missing column(s): " + ", ".join(sorted(missing))
            )
        for line_number, row in enumerate(reader, start=2):
            model_id = row["model_id"]
            if not model_id or model_id in by_model:
                raise ValueError(
                    f"Empty or duplicate decision model at line {line_number}"
                )
            by_model[model_id] = row
            rows.append(row)

    accepted = {
        model_id: [row]
        for model_id, row in by_model.items()
        if row["status"] == "accepted"
    }
    if not accepted:
        raise ValueError("Decision TSV contains no accepted models")
    redundant_rows = 0
    for row in rows:
        if row["reason"] != "redundant_projection":
            continue
        redundant_rows += 1
        parent_id = row["redundant_with_model_id"]
        if not parent_id or parent_id not in accepted:
            raise ValueError(
                f"Redundant model has no accepted parent: {row['model_id']}"
            )
        parent = by_model[parent_id]
        if row["seqid"] != parent["seqid"] or row["strand"] != parent["strand"]:
            raise ValueError(
                f"Redundant model disagrees with parent locus: {row['model_id']}"
            )
        accepted[parent_id].append(row)

    output_fields = (
        "model_id",
        "query_id",
        "source",
        "model_seqid",
        "model_start",
        "model_end",
        "model_span_bp",
        "model_strand",
        "model_score",
        "model_rank",
        "model_identity",
        "model_query_coverage",
        "baseline_existing_cds_overlap_fraction",
        "support_model_count",
        "support_query_count",
        "support_source_count",
        "support_sources",
        "support_rank1_model_count",
        "support_best_identity",
        "support_min_identity",
        "support_best_query_coverage",
        "support_min_query_coverage",
    )
    output_rows = []
    source_count_distribution: Counter[int] = Counter()
    for model_id, support in accepted.items():
        anchor = by_model[model_id]
        sources = sorted({row["source"] for row in support})
        identities = [float(row["identity"]) for row in support]
        coverages = [float(row["query_coverage"]) for row in support]
        source_count_distribution[len(sources)] += 1
        output_rows.append(
            {
                "model_id": model_id,
                "query_id": anchor["query_id"],
                "source": anchor["source"],
                "model_seqid": anchor["seqid"],
                "model_start": anchor["start"],
                "model_end": anchor["end"],
                "model_span_bp": int(anchor["end"]) - int(anchor["start"]) + 1,
                "model_strand": anchor["strand"],
                "model_score": anchor["score"],
                "model_rank": anchor["rank"],
                "model_identity": anchor["identity"],
                "model_query_coverage": anchor["query_coverage"],
                "baseline_existing_cds_overlap_fraction": anchor[
                    "existing_cds_overlap_fraction"
                ],
                "support_model_count": len(support),
                "support_query_count": len({row["query_id"] for row in support}),
                "support_source_count": len(sources),
                "support_sources": ",".join(sources),
                "support_rank1_model_count": sum(
                    int(int(row["rank"]) == 1) for row in support
                ),
                "support_best_identity": f"{max(identities):.6f}",
                "support_min_identity": f"{min(identities):.6f}",
                "support_best_query_coverage": f"{max(coverages):.6f}",
                "support_min_query_coverage": f"{min(coverages):.6f}",
            }
        )
    output_rows.sort(
        key=lambda row: (
            row["model_seqid"],
            int(row["model_start"]),
            row["model_id"],
        )
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=output_fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(output_rows)
    manifest = {
        "schema_version": PROJECTION_SUPPORT_SCHEMA_VERSION,
        "generator": {"name": "PloidyPatch", "version": __version__},
        "input": {
            "file_name": Path(decisions_tsv_path).name,
            "sha256": _file_sha256(decisions_tsv_path),
            "rows": len(rows),
        },
        "counts": {
            "accepted_models": len(accepted),
            "redundant_support_models": redundant_rows,
            "accepted_models_with_redundant_support": sum(
                int(len(support) > 1) for support in accepted.values()
            ),
            "support_source_count_distribution": {
                str(key): value
                for key, value in sorted(source_count_distribution.items())
            },
        },
        "output": {
            "file_name": output.name,
            "rows": len(output_rows),
            "sha256": _file_sha256(output),
        },
    }
    with manifest_path.open("x", encoding="utf-8", newline="") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return manifest
