from __future__ import annotations

import csv
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable

from . import __version__
from .gff import parse_attributes
from .perturb import MISSING_GENE_EVENT, TRUTH_SCHEMA_VERSION, _file_sha256
from .perturb import read_gff_document
from .score import CdsChainSignature, build_annotation_index


LOCALIZATION_SCORE_SCHEMA_VERSION = "ploidypatch.synteny_localization_score.v1"
MODEL_LABEL_SCHEMA_VERSION = "ploidypatch.synteny_model_labels.v1"


@dataclass(frozen=True)
class Span:
    item_id: str
    seqid: str
    start: int
    end: int
    strand: str

    @property
    def length(self) -> int:
        return self.end - self.start + 1


def _rate(numerator: int, denominator: int) -> float | None:
    return numerator / denominator if denominator else None


def _f1(precision: float | None, recall: float | None) -> float | None:
    if precision is None or recall is None or precision + recall == 0:
        return None
    return 2 * precision * recall / (precision + recall)


def _overlap_bp(left: Span, right: Span) -> int:
    if left.seqid != right.seqid:
        return 0
    return max(0, min(left.end, right.end) - max(left.start, right.start) + 1)


def _evaluate_relationship(
    candidates: list[Span],
    truth: list[Span],
    predicate: Callable[[Span, Span], bool],
) -> tuple[dict[str, int | float | None], set[str]]:
    truth_by_seqid: dict[str, list[Span]] = {}
    for event in truth:
        truth_by_seqid.setdefault(event.seqid, []).append(event)
    candidate_hits = 0
    truth_hits: set[str] = set()
    relationship_pairs = 0
    for candidate in candidates:
        matched = False
        for event in truth_by_seqid.get(candidate.seqid, []):
            if predicate(candidate, event):
                matched = True
                truth_hits.add(event.item_id)
                relationship_pairs += 1
        candidate_hits += int(matched)
    precision = _rate(candidate_hits, len(candidates))
    recall = _rate(len(truth_hits), len(truth))
    return (
        {
            "candidate_items": len(candidates),
            "truth_events": len(truth),
            "relationship_pairs": relationship_pairs,
            "candidate_hits": candidate_hits,
            "event_hits": len(truth_hits),
            "candidate_precision": precision,
            "event_recall": recall,
            "f1": _f1(precision, recall),
        },
        truth_hits,
    )


def _read_tsv(path: str | Path, required: set[str]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = required - set(reader.fieldnames or [])
        if reader.fieldnames is None or missing:
            raise ValueError(
                f"TSV {path} is missing column(s): " + ", ".join(sorted(missing))
            )
        rows.extend(reader)
    return rows


def _validated_span(
    *, item_id: str, seqid: str, start: str, end: str, strand: str
) -> Span:
    if not item_id or not seqid:
        raise ValueError("Span identifiers and sequence IDs must be non-empty")
    try:
        start_int = int(start)
        end_int = int(end)
    except ValueError as exc:
        raise ValueError(f"Non-integer span coordinate for {item_id}") from exc
    if start_int < 1 or end_int < start_int or strand not in {"+", "-", "."}:
        raise ValueError(f"Invalid span for {item_id}")
    return Span(item_id, seqid, start_int, end_int, strand)


def _read_baseline_model_records(
    candidate_gff_path: str | Path,
) -> dict[str, dict[str, Any]]:
    model_records: dict[str, dict[str, Any]] = {}
    current_model_id: str | None = None
    with Path(candidate_gff_path).open("r", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\r\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) != 9:
                raise ValueError(f"Malformed candidate GFF line {line_number}")
            if fields[1] != "PloidyPatchBaseline":
                continue
            feature_type = fields[2]
            if feature_type == "mRNA":
                attributes, malformed = parse_attributes(fields[8])
                model_id = attributes.get("miniprot_model")
                if malformed or not model_id or model_id in model_records:
                    raise ValueError(
                        f"Invalid or duplicate baseline mRNA at line {line_number}"
                    )
                current_model_id = model_id
                model_records[model_id] = {
                    "seqid": fields[0],
                    "start": int(fields[3]),
                    "end": int(fields[4]),
                    "strand": fields[6],
                    "cds": [],
                }
            elif feature_type == "CDS":
                if current_model_id is None:
                    raise ValueError(
                        f"Baseline CDS precedes mRNA at line {line_number}"
                    )
                model_records[current_model_id]["cds"].append(
                    (int(fields[3]), int(fields[4]), fields[7])
                )
    return model_records


def score_synteny_localization(
    *,
    gap_tsv_paths: list[str | Path],
    selection_tsv_path: str | Path,
    truth_path: str | Path,
    include_event_details: bool = False,
) -> dict[str, Any]:
    """Score blind gap localization separately from exact gene-model recovery."""

    if not gap_tsv_paths:
        raise ValueError("At least one upstream synteny-gap TSV is required")
    truth_document = json.loads(Path(truth_path).read_text(encoding="utf-8"))
    if truth_document.get("schema_version") != TRUTH_SCHEMA_VERSION:
        raise ValueError("Unsupported or missing hidden-truth schema version")
    unsupported = {
        event.get("event_type")
        for event in truth_document.get("events", [])
        if event.get("event_type") != MISSING_GENE_EVENT
    }
    if unsupported:
        raise ValueError(f"Unsupported event type(s): {unsupported}")

    truth_spans: list[Span] = []
    truth_events: dict[str, dict[str, Any]] = {}
    for event in truth_document.get("events", []):
        event_id = event["event_id"]
        if event_id in truth_events:
            raise ValueError(f"Duplicate hidden event ID: {event_id}")
        target = event["target"]
        truth_events[event_id] = event
        truth_spans.append(
            _validated_span(
                item_id=event_id,
                seqid=target["seqid"],
                start=str(target["start"]),
                end=str(target["end"]),
                strand=target["strand"],
            )
        )
    if not truth_spans:
        raise ValueError("Hidden truth contains no events")

    gap_required = {
        "candidate_id",
        "gap_id",
        "source_label",
        "query_seqid",
        "locus_start",
        "locus_end",
        "target_gene_id",
    }
    gap_rows: dict[str, dict[str, str]] = {}
    gaps_by_id: dict[str, Span] = {}
    gap_input_manifest = []
    for raw_path in gap_tsv_paths:
        rows = _read_tsv(raw_path, gap_required)
        for row in rows:
            candidate_id = row["candidate_id"]
            if not candidate_id or candidate_id in gap_rows:
                raise ValueError(f"Empty or duplicate gap candidate ID: {candidate_id}")
            gap_rows[candidate_id] = row
            span = _validated_span(
                item_id=row["gap_id"],
                seqid=row["query_seqid"],
                start=row["locus_start"],
                end=row["locus_end"],
                strand=".",
            )
            previous = gaps_by_id.get(span.item_id)
            if previous is not None and previous != span:
                raise ValueError(f"Inconsistent coordinates for gap {span.item_id}")
            gaps_by_id[span.item_id] = span
        gap_input_manifest.append(
            {
                "file_name": Path(raw_path).name,
                "sha256": _file_sha256(raw_path),
                "gene_hypotheses": len(rows),
            }
        )

    selection_required = {
        *gap_required,
        "model_id",
        "model_seqid",
        "model_start",
        "model_end",
        "model_strand",
    }
    selection_rows = _read_tsv(selection_tsv_path, selection_required)
    selected_gap_ids: set[str] = set()
    model_ids: set[str] = set()
    model_spans: list[Span] = []
    consistency_fields = (
        "gap_id",
        "source_label",
        "query_seqid",
        "locus_start",
        "locus_end",
        "target_gene_id",
    )
    for row in selection_rows:
        candidate_id = row["candidate_id"]
        upstream = gap_rows.get(candidate_id)
        if upstream is None:
            raise ValueError(f"Selection references unknown candidate: {candidate_id}")
        if any(row[field] != upstream[field] for field in consistency_fields):
            raise ValueError(f"Selection/upstream gap mismatch for {candidate_id}")
        model_id = row["model_id"]
        if not model_id or model_id in model_ids:
            raise ValueError(f"Empty or duplicate selected model ID: {model_id}")
        model_ids.add(model_id)
        model = _validated_span(
            item_id=model_id,
            seqid=row["model_seqid"],
            start=row["model_start"],
            end=row["model_end"],
            strand=row["model_strand"],
        )
        gap = gaps_by_id[row["gap_id"]]
        if (
            model.seqid != gap.seqid
            or model.start < gap.start
            or model.end > gap.end
        ):
            raise ValueError(f"Selected model is outside its gap: {model_id}")
        model_spans.append(model)
        selected_gap_ids.add(row["gap_id"])

    all_gaps = sorted(gaps_by_id.values(), key=lambda span: span.item_id)
    selected_gaps = sorted(
        (gaps_by_id[gap_id] for gap_id in selected_gap_ids),
        key=lambda span: span.item_id,
    )
    model_spans.sort(key=lambda span: span.item_id)

    def overlaps(candidate: Span, event: Span) -> bool:
        return _overlap_bp(candidate, event) > 0

    def contains(candidate: Span, event: Span) -> bool:
        return candidate.start <= event.start and candidate.end >= event.end

    def same_strand_overlap(candidate: Span, event: Span) -> bool:
        return candidate.strand == event.strand and overlaps(candidate, event)

    def model_center_inside_truth(candidate: Span, event: Span) -> bool:
        center = (candidate.start + candidate.end) / 2
        return (
            candidate.strand == event.strand
            and event.start <= center <= event.end
        )

    def truth_center_inside_model(candidate: Span, event: Span) -> bool:
        center = (event.start + event.end) / 2
        return (
            candidate.strand == event.strand
            and candidate.start <= center <= candidate.end
        )

    def reciprocal_overlap_50(candidate: Span, event: Span) -> bool:
        overlap = _overlap_bp(candidate, event)
        return (
            candidate.strand == event.strand
            and overlap / candidate.length >= 0.5
            and overlap / event.length >= 0.5
        )

    all_gap_overlap, _ = _evaluate_relationship(all_gaps, truth_spans, overlaps)
    all_gap_containment, _ = _evaluate_relationship(all_gaps, truth_spans, contains)
    selected_gap_overlap, _ = _evaluate_relationship(
        selected_gaps, truth_spans, overlaps
    )
    selected_gap_containment, _ = _evaluate_relationship(
        selected_gaps, truth_spans, contains
    )
    model_metrics: dict[str, Any] = {}
    event_hits_by_metric: dict[str, set[str]] = {}
    for label, predicate in (
        ("any_span_overlap", overlaps),
        ("same_strand_overlap", same_strand_overlap),
        ("same_strand_model_center_inside_truth", model_center_inside_truth),
        ("same_strand_truth_center_inside_model", truth_center_inside_model),
        ("same_strand_reciprocal_overlap_at_least_0_5", reciprocal_overlap_50),
    ):
        metrics, event_hits = _evaluate_relationship(
            model_spans, truth_spans, predicate
        )
        model_metrics[label] = metrics
        event_hits_by_metric[label] = event_hits

    report: dict[str, Any] = {
        "schema_version": LOCALIZATION_SCORE_SCHEMA_VERSION,
        "evaluator": {"name": "PloidyPatch", "version": __version__},
        "inputs": {
            "synteny_gaps": gap_input_manifest,
            "selection": {
                "file_name": Path(selection_tsv_path).name,
                "sha256": _file_sha256(selection_tsv_path),
                "models": len(model_spans),
            },
            "hidden_truth": {
                "file_name": Path(truth_path).name,
                "sha256": _file_sha256(truth_path),
                "events": len(truth_spans),
            },
        },
        "all_blind_gap_loci": {
            "any_span_overlap": all_gap_overlap,
            "full_gene_span_containment": all_gap_containment,
        },
        "selected_gap_loci": {
            "any_span_overlap": selected_gap_overlap,
            "full_gene_span_containment": selected_gap_containment,
        },
        "selected_model_spans": model_metrics,
        "quality_gate": {
            "selection_rows_trace_to_upstream_gaps": True,
            "selected_models_fully_contained_in_gaps": True,
            "selected_model_ids_unique": True,
            "grade": "pass",
        },
    }
    if include_event_details:
        report["event_details"] = [
            {
                "event_id": event.item_id,
                "target_gene_id": truth_events[event.item_id]["target"]["gene_id"],
                **{
                    label: event.item_id in event_hits
                    for label, event_hits in event_hits_by_metric.items()
                },
            }
            for event in truth_spans
        ]
    return report


def write_synteny_model_labels(
    *,
    source_gff_path: str | Path,
    candidate_gff_path: str | Path,
    selection_tsv_path: str | Path,
    baseline_decisions_tsv_path: str | Path,
    truth_path: str | Path,
    output_tsv_path: str | Path,
    control_candidate_gff_path: str | Path | None = None,
) -> dict[str, Any]:
    """Write evaluator-only per-model labels for calibration and diagnostics."""

    output = Path(output_tsv_path)
    manifest_path = Path(str(output) + ".manifest.json")
    collisions = [path for path in (output, manifest_path) if path.exists()]
    if collisions:
        raise FileExistsError(
            "Refusing to overwrite model-label artifact(s): "
            + ", ".join(str(path) for path in collisions)
        )

    truth_document = json.loads(Path(truth_path).read_text(encoding="utf-8"))
    if truth_document.get("schema_version") != TRUTH_SCHEMA_VERSION:
        raise ValueError("Unsupported or missing hidden-truth schema version")
    unsupported = {
        event.get("event_type")
        for event in truth_document.get("events", [])
        if event.get("event_type") != MISSING_GENE_EVENT
    }
    if unsupported:
        raise ValueError(f"Unsupported event type(s): {unsupported}")

    source_document = read_gff_document(source_gff_path)
    if source_document.text_sha256 != truth_document["source"]["text_sha256"]:
        raise ValueError("Source GFF3 text checksum does not match hidden truth")
    source_index = build_annotation_index(source_document)
    truth_spans: list[Span] = []
    hidden_cds_chains: set[CdsChainSignature] = set()
    for event in truth_document.get("events", []):
        target = event["target"]
        truth_spans.append(
            _validated_span(
                item_id=event["event_id"],
                seqid=target["seqid"],
                start=str(target["start"]),
                end=str(target["end"]),
                strand=target["strand"],
            )
        )
        for transcript_id in target["transcript_ids"]:
            transcript = source_index.transcripts.get(transcript_id)
            if transcript is None:
                raise ValueError(
                    f"Hidden transcript is absent from source: {transcript_id}"
                )
            signature = transcript.signature
            if signature.cds:
                hidden_cds_chains.add(
                    CdsChainSignature(
                        signature.seqid, signature.strand, signature.cds
                    )
                )

    decision_required = {
        "model_id",
        "status",
        "reason",
        "existing_cds_overlap_fraction",
    }
    decisions = _read_tsv(baseline_decisions_tsv_path, decision_required)
    decision_by_model: dict[str, dict[str, str]] = {}
    for row in decisions:
        model_id = row["model_id"]
        if not model_id or model_id in decision_by_model:
            raise ValueError(f"Empty or duplicate baseline decision: {model_id}")
        decision_by_model[model_id] = row

    model_records = _read_baseline_model_records(candidate_gff_path)
    control_cds_chains: set[CdsChainSignature] = set()
    if control_candidate_gff_path is not None:
        for record in _read_baseline_model_records(
            control_candidate_gff_path
        ).values():
            if record["cds"]:
                control_cds_chains.add(
                    CdsChainSignature(
                        record["seqid"],
                        record["strand"],
                        tuple(sorted(set(record["cds"]))),
                    )
                )

    with Path(selection_tsv_path).open(
        "r", encoding="utf-8", newline=""
    ) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        selection_fields = tuple(reader.fieldnames or ())
        required = {
            "model_id",
            "model_seqid",
            "model_start",
            "model_end",
            "model_strand",
        }
        missing = required - set(selection_fields)
        if not selection_fields or missing:
            raise ValueError(
                "Selection TSV is missing column(s): " + ", ".join(sorted(missing))
            )
        selection_rows = list(reader)

    overlap_field = "baseline_existing_cds_overlap_fraction"
    label_fields = (
        "label_paired_differential_cds_chain",
        "label_any_span_overlap",
        "label_same_strand_overlap",
        "label_same_strand_model_center_inside_truth",
        "label_same_strand_reciprocal_overlap_50",
        "label_exact_cds_chain",
    )
    if set(label_fields) & set(selection_fields):
        raise ValueError("Selection TSV collides with evaluator label columns")
    output_extra_fields = (
        label_fields
        if overlap_field in selection_fields
        else (overlap_field, *label_fields)
    )
    output_rows = []
    label_counts = {field: 0 for field in label_fields if field.startswith("label_")}
    seen_models: set[str] = set()
    truth_by_seqid: dict[str, list[Span]] = {}
    for event in truth_spans:
        truth_by_seqid.setdefault(event.seqid, []).append(event)
    for row in selection_rows:
        model_id = row["model_id"]
        if not model_id or model_id in seen_models:
            raise ValueError(f"Empty or duplicate selected model: {model_id}")
        seen_models.add(model_id)
        decision = decision_by_model.get(model_id)
        if decision is None or decision["status"] != "accepted":
            raise ValueError(f"Selected model lacks an accepted decision: {model_id}")
        if overlap_field in row and float(row[overlap_field]) != float(
            decision["existing_cds_overlap_fraction"]
        ):
            raise ValueError(f"Selection/decision overlap mismatch for {model_id}")
        model_record = model_records.get(model_id)
        if model_record is None or not model_record["cds"]:
            raise ValueError(f"Selected model lacks candidate CDS records: {model_id}")
        model_span = _validated_span(
            item_id=model_id,
            seqid=row["model_seqid"],
            start=row["model_start"],
            end=row["model_end"],
            strand=row["model_strand"],
        )
        if (
            model_record["seqid"] != model_span.seqid
            or model_record["start"] != model_span.start
            or model_record["end"] != model_span.end
            or model_record["strand"] != model_span.strand
        ):
            raise ValueError(f"Selection/candidate GFF mismatch for {model_id}")
        events = truth_by_seqid.get(model_span.seqid, [])
        any_overlap = any(_overlap_bp(model_span, event) > 0 for event in events)
        same_strand = [
            event
            for event in events
            if event.strand == model_span.strand
            and _overlap_bp(model_span, event) > 0
        ]
        model_center = (model_span.start + model_span.end) / 2
        center_inside = any(
            event.start <= model_center <= event.end for event in same_strand
        )
        reciprocal_50 = any(
            _overlap_bp(model_span, event) / model_span.length >= 0.5
            and _overlap_bp(model_span, event) / event.length >= 0.5
            for event in same_strand
        )
        cds_chain = CdsChainSignature(
            model_record["seqid"],
            model_record["strand"],
            tuple(sorted(set(model_record["cds"]))),
        )
        labels = {
            "label_paired_differential_cds_chain": int(
                not control_cds_chains or cds_chain not in control_cds_chains
            ),
            "label_any_span_overlap": int(any_overlap),
            "label_same_strand_overlap": int(bool(same_strand)),
            "label_same_strand_model_center_inside_truth": int(center_inside),
            "label_same_strand_reciprocal_overlap_50": int(reciprocal_50),
            "label_exact_cds_chain": int(cds_chain in hidden_cds_chains),
        }
        for field, value in labels.items():
            label_counts[field] += value
        output_rows.append(
            {
                **row,
                overlap_field: decision[
                    "existing_cds_overlap_fraction"
                ],
                **labels,
            }
        )

    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=(*selection_fields, *output_extra_fields),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(output_rows)
    manifest = {
        "schema_version": MODEL_LABEL_SCHEMA_VERSION,
        "generator": {"name": "PloidyPatch", "version": __version__},
        "access": "evaluator_only",
        "inputs": {
            "source_gff_sha256": _file_sha256(source_gff_path),
            "candidate_gff_sha256": _file_sha256(candidate_gff_path),
            "selection_tsv_sha256": _file_sha256(selection_tsv_path),
            "baseline_decisions_sha256": _file_sha256(
                baseline_decisions_tsv_path
            ),
            "hidden_truth_sha256": _file_sha256(truth_path),
        },
        "counts": {
            "selected_models": len(output_rows),
            "hidden_events": len(truth_spans),
            "hidden_cds_chains": len(hidden_cds_chains),
            **label_counts,
        },
        "output": {
            "file_name": output.name,
            "rows": len(output_rows),
            "sha256": _file_sha256(output),
        },
    }
    if control_candidate_gff_path is not None:
        manifest["inputs"]["control_candidate_gff_sha256"] = _file_sha256(
            control_candidate_gff_path
        )
        manifest["counts"]["control_novel_cds_chains"] = len(control_cds_chains)
    with manifest_path.open("x", encoding="utf-8", newline="") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return manifest
