from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter, defaultdict, deque
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from . import __version__
from .gff import parse_attributes
from .io import open_text


TRUTH_SCHEMA_VERSION = "ploidypatch.hidden_truth.v1"
MANIFEST_SCHEMA_VERSION = "ploidypatch.perturbation_manifest.v1"
MISSING_GENE_EVENT = "annotation_missing_gene"
MISSING_INTERNAL_EXON_EVENT = "annotation_missing_internal_exon"
BOUNDARY_SHIFT_EVENT = "annotation_boundary_shift"
SPLIT_GENE_EVENT = "annotation_split_gene"
FUSED_GENE_EVENT = "annotation_fused_gene"
COPY_COLLAPSE_EVENT = "annotation_copy_collapse"
SUPPORTED_ANNOTATION_EVENTS = frozenset(
    {
        MISSING_GENE_EVENT,
        MISSING_INTERNAL_EXON_EVENT,
        BOUNDARY_SHIFT_EVENT,
        SPLIT_GENE_EVENT,
        FUSED_GENE_EVENT,
        COPY_COLLAPSE_EVENT,
    }
)
SAMPLE_MANIFEST_SCHEMA_VERSION = "ploidypatch.candidate_sample_manifest.v1"


@dataclass(frozen=True)
class GffRecord:
    """A parsed feature line whose original text and position are retained."""

    line_number: int
    raw_line: str
    seqid: str
    feature_type: str
    start: int
    end: int
    strand: str
    phase: str
    attributes: dict[str, str]

    @property
    def feature_id(self) -> str | None:
        return self.attributes.get("ID")

    @property
    def parents(self) -> tuple[str, ...]:
        return tuple(
            value for value in self.attributes.get("Parent", "").split(",") if value
        )


@dataclass(frozen=True)
class GffDocument:
    """A lossless, decompressed GFF3 document plus parsed feature records."""

    lines: tuple[str, ...]
    records: tuple[GffRecord, ...]

    @property
    def text_sha256(self) -> str:
        return _text_sha256(self.lines)


@dataclass(frozen=True)
class MissingGeneCandidate:
    gene: GffRecord
    removed_records: tuple[GffRecord, ...]
    transcript_ids: tuple[str, ...]
    feature_type_counts: dict[str, int]


def _text_sha256(lines: tuple[str, ...] | list[str]) -> str:
    digest = hashlib.sha256()
    for line in lines:
        digest.update(line.encode("utf-8"))
    return digest.hexdigest()


def _file_sha256(path: str | Path, block_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        while block := handle.read(block_size):
            digest.update(block)
    return digest.hexdigest()


def read_gff_document(
    path: str | Path,
    include_seqids: frozenset[str] | None = None,
) -> GffDocument:
    """Parse a GFF3 while preserving every decompressed input line.

    When ``include_seqids`` is supplied, feature records on other sequences are
    not parsed. The original text is still retained so its checksum continues
    to identify the complete source file.
    """

    lines: list[str] = []
    records: list[GffRecord] = []
    with open_text(path) as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            lines.append(raw_line)
            stripped = raw_line.rstrip("\r\n")
            if not stripped or stripped.startswith("#"):
                continue
            fields = stripped.split("\t")
            if include_seqids is not None and fields[0] not in include_seqids:
                continue
            if len(fields) != 9:
                raise ValueError(
                    f"Malformed GFF3 feature line {line_number}: expected 9 columns"
                )
            try:
                start = int(fields[3])
                end = int(fields[4])
            except ValueError as exc:
                raise ValueError(
                    f"Malformed GFF3 coordinates on line {line_number}"
                ) from exc
            if start < 1 or end < start:
                raise ValueError(f"Invalid GFF3 interval on line {line_number}")
            attributes, malformed = parse_attributes(fields[8])
            if malformed:
                raise ValueError(
                    f"Malformed GFF3 attributes on line {line_number}: "
                    f"{malformed} invalid field(s)"
                )
            records.append(
                GffRecord(
                    line_number=line_number,
                    raw_line=raw_line,
                    seqid=fields[0],
                    feature_type=fields[2],
                    start=start,
                    end=end,
                    strand=fields[6],
                    phase=fields[7],
                    attributes=attributes,
                )
            )
    return GffDocument(lines=tuple(lines), records=tuple(records))


def find_missing_gene_candidates(
    document: GffDocument, *, require_exon: bool = True
) -> list[MissingGeneCandidate]:
    """Return coding gene hierarchies removable without touching shared features.

    The first benchmark primitive is deliberately conservative: a gene must
    have a unique ID, at least one transcript and CDS, and every feature in its
    descendant closure must have only one parent. The default also requires an
    explicit exon; callers that only remove complete coding hierarchies may
    disable that requirement. This prevents a
    synthetic missing-gene event from silently damaging a neighboring or
    alternatively shared model.
    """

    records_by_id: dict[str, list[GffRecord]] = defaultdict(list)
    children_by_parent: dict[str, list[GffRecord]] = defaultdict(list)
    for record in document.records:
        if record.feature_id:
            records_by_id[record.feature_id].append(record)
        for parent in record.parents:
            children_by_parent[parent].append(record)

    candidates: list[MissingGeneCandidate] = []
    genes = [record for record in document.records if record.feature_type == "gene"]
    for gene in genes:
        gene_id = gene.feature_id
        if not gene_id or len(records_by_id[gene_id]) != 1:
            continue

        queue: deque[str] = deque([gene_id])
        visited_ids = {gene_id}
        descendant_by_line: dict[int, GffRecord] = {}
        has_shared_feature = False
        while queue:
            parent_id = queue.popleft()
            for child in children_by_parent.get(parent_id, []):
                if len(child.parents) != 1:
                    has_shared_feature = True
                    continue
                descendant_by_line[child.line_number] = child
                child_id = child.feature_id
                if child_id and child_id not in visited_ids:
                    visited_ids.add(child_id)
                    queue.append(child_id)

        if has_shared_feature:
            continue
        descendants = tuple(
            descendant_by_line[key] for key in sorted(descendant_by_line)
        )
        counts = Counter(record.feature_type for record in descendants)
        transcript_ids = tuple(
            sorted(
                record.feature_id
                for record in descendants
                if record.feature_type in {"mRNA", "transcript"}
                and record.feature_id is not None
            )
        )
        if (
            not transcript_ids
            or counts["CDS"] == 0
            or (require_exon and counts["exon"] == 0)
        ):
            continue
        removed_records = tuple(sorted((gene, *descendants), key=lambda x: x.line_number))
        feature_type_counts = dict(
            sorted(Counter(record.feature_type for record in removed_records).items())
        )
        candidates.append(
            MissingGeneCandidate(
                gene=gene,
                removed_records=removed_records,
                transcript_ids=transcript_ids,
                feature_type_counts=feature_type_counts,
            )
        )

    return sorted(
        candidates,
        key=lambda item: (
            item.gene.seqid,
            item.gene.start,
            item.gene.end,
            item.gene.feature_id or "",
        ),
    )


def _candidate_rank(seed: int, candidate: MissingGeneCandidate) -> str:
    gene = candidate.gene
    identity = (
        f"{seed}\0{gene.seqid}\0{gene.start}\0{gene.end}\0"
        f"{gene.strand}\0{gene.feature_id}"
    )
    return hashlib.sha256(identity.encode("utf-8")).hexdigest()


def candidate_id_for_gene(source_text_sha256: str, gene_id: str) -> str:
    payload = f"{source_text_sha256}\0{MISSING_GENE_EVENT}\0{gene_id}"
    token = hashlib.sha256(payload.encode("utf-8")).hexdigest()[:20]
    return f"PPC-MG-{token}"


def _read_evaluator_selection(
    selection_path: str | Path,
    source_text_sha256: str,
    candidates: list[MissingGeneCandidate],
    seed: int,
) -> tuple[list[MissingGeneCandidate], dict[str, Any]]:
    selection_path = Path(selection_path)
    manifest_path = Path(str(selection_path) + ".manifest.json")
    if not manifest_path.is_file():
        raise FileNotFoundError(
            f"Evaluator selection manifest is required: {manifest_path}"
        )
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("schema_version") != SAMPLE_MANIFEST_SCHEMA_VERSION:
        raise ValueError("Unsupported evaluator selection manifest schema")
    selection_sha256 = _file_sha256(selection_path)
    if manifest.get("output", {}).get("sha256") != selection_sha256:
        raise ValueError("Evaluator selection checksum does not match its manifest")
    if manifest.get("source", {}).get("text_sha256") != source_text_sha256:
        raise ValueError("Evaluator selection belongs to a different source GFF3")
    if manifest.get("selection", {}).get("seed") != seed:
        raise ValueError("Evaluator selection was generated with a different seed")

    candidates_by_gene = {
        candidate.gene.feature_id: candidate
        for candidate in candidates
        if candidate.gene.feature_id is not None
    }
    selected: list[MissingGeneCandidate] = []
    seen_gene_ids: set[str] = set()
    with selection_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or not {"candidate_id", "gene_id"} <= set(
            reader.fieldnames
        ):
            raise ValueError("Evaluator selection lacks candidate_id or gene_id")
        for line_number, row in enumerate(reader, start=2):
            gene_id = row["gene_id"]
            if gene_id in seen_gene_ids:
                raise ValueError(f"Duplicate selected gene_id: {gene_id}")
            seen_gene_ids.add(gene_id)
            candidate = candidates_by_gene.get(gene_id)
            if candidate is None:
                raise ValueError(
                    f"Selected gene is not eligible in the source at line "
                    f"{line_number}: {gene_id}"
                )
            expected_candidate_id = candidate_id_for_gene(
                source_text_sha256, gene_id
            )
            if row["candidate_id"] != expected_candidate_id:
                raise ValueError(
                    f"Candidate ID mismatch at evaluator selection line {line_number}"
                )
            selected.append(candidate)
    if not selected:
        raise ValueError("Evaluator selection contains no candidates")
    if manifest.get("output", {}).get("rows") != len(selected):
        raise ValueError("Evaluator selection row count does not match its manifest")
    return selected, {
        "file_name": selection_path.name,
        "sha256": selection_sha256,
        "manifest_file_name": manifest_path.name,
        "manifest_sha256": _file_sha256(manifest_path),
    }


def _event_id(
    source_text_sha256: str,
    seed: int,
    candidate: MissingGeneCandidate,
) -> str:
    identity = (
        f"{source_text_sha256}\0{seed}\0{MISSING_GENE_EVENT}\0"
        f"{candidate.gene.feature_id}"
    )
    token = hashlib.sha256(identity.encode("utf-8")).hexdigest()[:20]
    return f"PPB-MG-{token}"


def _write_text_exclusive(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8", newline="") as handle:
        handle.write(text)


def _json_text(value: Any) -> str:
    return json.dumps(value, indent=2, sort_keys=True, ensure_ascii=False) + "\n"


def generate_missing_gene_benchmark(
    gff_path: str | Path,
    output_dir: str | Path,
    count: int | None,
    seed: int,
    selection_tsv_path: str | Path | None = None,
    truth_output_dir: str | Path | None = None,
) -> dict[str, Any]:
    """Create a deterministic missing-annotation benchmark and hidden truth."""

    if (count is None) == (selection_tsv_path is None):
        raise ValueError("Specify exactly one of count or selection_tsv_path")
    if count is not None and count < 1:
        raise ValueError("count must be at least 1")
    source_path = Path(gff_path)
    output_dir = Path(output_dir)
    truth_dir = Path(truth_output_dir) if truth_output_dir is not None else output_dir
    output_gff = output_dir / "perturbed.gff3"
    truth_path = truth_dir / "hidden_truth.json"
    manifest_path = output_dir / "manifest.json"
    collisions = [path for path in (output_gff, truth_path, manifest_path) if path.exists()]
    if collisions:
        names = ", ".join(str(path) for path in collisions)
        raise FileExistsError(f"Refusing to overwrite benchmark output(s): {names}")

    document = read_gff_document(source_path)
    source_text_sha256 = document.text_sha256
    candidates = find_missing_gene_candidates(document)
    selection_source: dict[str, Any] | None = None
    if selection_tsv_path is not None:
        selected, selection_source = _read_evaluator_selection(
            selection_tsv_path, source_text_sha256, candidates, seed
        )
        selection_algorithm = "evaluator_selection_tsv_v1"
        requested_count = len(selected)
    else:
        assert count is not None
        if len(candidates) < count:
            raise ValueError(
                f"Requested {count} events but only {len(candidates)} genes satisfy "
                "the conservative missing-gene eligibility rules"
            )
        selected = sorted(candidates, key=lambda item: _candidate_rank(seed, item))[
            :count
        ]
        selection_algorithm = "sha256_rank_v1"
        requested_count = count
    selected = sorted(selected, key=lambda item: item.gene.line_number)

    removed_line_numbers: set[int] = set()
    events: list[dict[str, Any]] = []
    for candidate in selected:
        overlap = removed_line_numbers & {
            record.line_number for record in candidate.removed_records
        }
        if overlap:
            raise AssertionError(f"Selected gene hierarchies overlap at lines {overlap}")
        removed_line_numbers.update(
            record.line_number for record in candidate.removed_records
        )
        gene = candidate.gene
        events.append(
            {
                "event_id": _event_id(source_text_sha256, seed, candidate),
                "event_type": MISSING_GENE_EVENT,
                "target": {
                    "gene_id": gene.feature_id,
                    "seqid": gene.seqid,
                    "start": gene.start,
                    "end": gene.end,
                    "strand": gene.strand,
                    "transcript_ids": list(candidate.transcript_ids),
                },
                "feature_type_counts": candidate.feature_type_counts,
                "removed_records": [
                    {
                        "line_number": record.line_number,
                        "line_sha256": hashlib.sha256(
                            record.raw_line.encode("utf-8")
                        ).hexdigest(),
                        "raw_line": record.raw_line,
                    }
                    for record in candidate.removed_records
                ],
            }
        )

    perturbed_lines = [
        line
        for line_number, line in enumerate(document.lines, start=1)
        if line_number not in removed_line_numbers
    ]
    truth: dict[str, Any] = {
        "schema_version": TRUTH_SCHEMA_VERSION,
        "generator": {"name": "PloidyPatch", "version": __version__},
        "source": {
            "file_name": source_path.name,
            "file_sha256": _file_sha256(source_path),
            "text_sha256": source_text_sha256,
            "line_count": len(document.lines),
        },
        "perturbation": {
            "event_type": MISSING_GENE_EVENT,
            "seed": seed,
            "selection_algorithm": selection_algorithm,
            "selection_source": selection_source,
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
        "source": truth["source"],
        "perturbation": {
            "event_type": MISSING_GENE_EVENT,
            "seed": seed,
            "selection_algorithm": selection_algorithm,
            "selection_source": selection_source,
            "requested_events": requested_count,
            "eligible_genes": len(candidates),
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


def restore_gff_from_truth(
    perturbed_gff_path: str | Path,
    truth_path: str | Path,
    output_gff_path: str | Path,
) -> dict[str, Any]:
    """Invert a filtering perturbation and verify the original text checksum."""

    output_path = Path(output_gff_path)
    if output_path.exists():
        raise FileExistsError(f"Refusing to overwrite restored GFF3: {output_path}")
    truth = json.loads(Path(truth_path).read_text(encoding="utf-8"))
    if truth.get("schema_version") != TRUTH_SCHEMA_VERSION:
        raise ValueError("Unsupported or missing hidden-truth schema version")

    edits_by_line: dict[int, tuple[str, tuple[str, ...]]] = {}
    for event in truth.get("events", []):
        line_edits = event.get("line_edits")
        if line_edits is None:
            line_edits = [
                {
                    "source_line_number": record["line_number"],
                    "source_line_sha256": record["line_sha256"],
                    "source_raw_line": record["raw_line"],
                    "perturbed_lines": [],
                }
                for record in event.get("removed_records", [])
            ]
        for edit in line_edits:
            line_number = int(edit["source_line_number"])
            source_raw_line = edit["source_raw_line"]
            expected_source_sha = edit["source_line_sha256"]
            observed_source_sha = hashlib.sha256(
                source_raw_line.encode("utf-8")
            ).hexdigest()
            if observed_source_sha != expected_source_sha:
                raise ValueError(
                    f"Hidden-truth source-line checksum failed at {line_number}"
                )
            perturbed_lines: list[str] = []
            for replacement in edit.get("perturbed_lines", []):
                raw_line = replacement["raw_line"]
                expected_sha = replacement["line_sha256"]
                observed_sha = hashlib.sha256(raw_line.encode("utf-8")).hexdigest()
                if observed_sha != expected_sha:
                    raise ValueError(
                        "Hidden-truth perturbed-line checksum failed at "
                        f"source line {line_number}"
                    )
                perturbed_lines.append(raw_line)
            if line_number in edits_by_line:
                raise ValueError(f"Duplicate hidden-truth line number: {line_number}")
            edits_by_line[line_number] = (
                source_raw_line,
                tuple(perturbed_lines),
            )

    with open_text(perturbed_gff_path) as handle:
        perturbed_lines = list(handle)
    original_line_count = int(truth["source"]["line_count"])
    expected_perturbed_count = original_line_count + sum(
        len(replacements) - 1 for _, replacements in edits_by_line.values()
    )
    if len(perturbed_lines) != expected_perturbed_count:
        raise ValueError(
            "Perturbed GFF3 line count is inconsistent with the hidden truth: "
            f"expected {expected_perturbed_count}, observed {len(perturbed_lines)}"
        )

    restored_lines: list[str] = []
    perturbed_index = 0
    for line_number in range(1, original_line_count + 1):
        if line_number in edits_by_line:
            source_raw_line, expected_replacements = edits_by_line[line_number]
            observed_replacements = tuple(
                perturbed_lines[
                    perturbed_index : perturbed_index + len(expected_replacements)
                ]
            )
            if observed_replacements != expected_replacements:
                raise ValueError(
                    "Perturbed GFF3 does not match hidden-truth replacement at "
                    f"source line {line_number}"
                )
            perturbed_index += len(expected_replacements)
            restored_lines.append(source_raw_line)
        else:
            if perturbed_index >= len(perturbed_lines):
                raise ValueError("Perturbed GFF3 ended before restoration completed")
            restored_lines.append(perturbed_lines[perturbed_index])
            perturbed_index += 1
    if perturbed_index != len(perturbed_lines):
        raise ValueError("Perturbed GFF3 contains unexpected trailing records")

    observed_sha = _text_sha256(restored_lines)
    expected_sha = truth["source"]["text_sha256"]
    if observed_sha != expected_sha:
        raise ValueError(
            "Restored GFF3 does not match the source text checksum; the perturbed "
            "input or hidden truth was modified"
        )
    _write_text_exclusive(output_path, "".join(restored_lines))
    return {
        "restored_gff": str(output_path),
        "text_sha256": observed_sha,
        "restored_records": len(edits_by_line),
        "events": len(truth.get("events", [])),
    }
