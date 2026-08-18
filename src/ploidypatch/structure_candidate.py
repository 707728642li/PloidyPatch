from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

from . import __version__
from .baseline import (
    ProjectedModel,
    _existing_cds_intervals,
    _file_sha256,
    _merge_intervals,
    _model_intervals,
    _overlap_bp,
    _overlap_with_index,
    _parse_miniprot_models,
    _read_protein_map,
)
from .io import open_text
from .perturb import read_gff_document
from .score import build_annotation_index


STRUCTURE_CANDIDATE_SCHEMA_VERSION = (
    "ploidypatch.miniprot_structure_candidate_manifest.v1"
)
GENE_BIN_SIZE = 100_000


def _chain_key(model: ProjectedModel) -> tuple[Any, ...]:
    return (model.seqid, model.strand, model.cds)


def _chain_digest(key: tuple[Any, ...]) -> str:
    payload = json.dumps(key, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()[:20]


def _annotation_context(
    annotation_gff_path: str | Path,
) -> tuple[
    set[tuple[Any, ...]],
    dict[tuple[str, str, int], tuple[str, ...]],
    dict[str, tuple[tuple[int, int], ...]],
]:
    document = read_gff_document(annotation_gff_path)
    index = build_annotation_index(document)
    exact_chains = {
        (signature.seqid, signature.strand, signature.cds)
        for signature in index.signatures
        if signature.cds
    }
    gene_intervals: dict[str, tuple[tuple[int, int], ...]] = {}
    mutable_bins: dict[tuple[str, str, int], list[str]] = defaultdict(list)
    for gene_id, signatures in index.gene_signatures.items():
        cds = [
            (start, end)
            for signature in signatures
            for start, end, _ in signature.cds
        ]
        if not cds:
            continue
        seqids = {signature.seqid for signature in signatures if signature.cds}
        strands = {signature.strand for signature in signatures if signature.cds}
        if len(seqids) != 1 or len(strands) != 1:
            raise ValueError(f"Gene spans incompatible loci: {gene_id}")
        intervals = tuple(_merge_intervals(cds))
        gene_intervals[gene_id] = intervals
        seqid = next(iter(seqids))
        strand = next(iter(strands))
        for bin_number in range(
            intervals[0][0] // GENE_BIN_SIZE,
            intervals[-1][1] // GENE_BIN_SIZE + 1,
        ):
            mutable_bins[(seqid, strand, bin_number)].append(gene_id)
    gene_bins = {
        key: tuple(sorted(set(values))) for key, values in mutable_bins.items()
    }
    return exact_chains, gene_bins, gene_intervals


def _overlapping_genes(
    model: ProjectedModel,
    *,
    gene_bins: dict[tuple[str, str, int], tuple[str, ...]],
    gene_intervals: dict[str, tuple[tuple[int, int], ...]],
    min_fraction: float,
) -> tuple[str, ...]:
    candidate = _model_intervals(model)
    nearby: set[str] = set()
    for bin_number in range(
        model.start // GENE_BIN_SIZE,
        model.end // GENE_BIN_SIZE + 1,
    ):
        nearby.update(gene_bins.get((model.seqid, model.strand, bin_number), ()))
    denominator = model.cds_bp
    supported = [
        gene_id
        for gene_id in nearby
        if denominator
        and _overlap_bp(candidate, list(gene_intervals[gene_id])) / denominator
        >= min_fraction
    ]
    return tuple(sorted(supported))


def _write_structure_candidate(
    handle: Any,
    model: ProjectedModel,
    report: dict[str, Any],
) -> None:
    """Write one traceable structural candidate with stable digest-based IDs."""

    digest = report["group_id"].removeprefix("PPSC-")
    gene_id = f"PPSC_gene_{digest}"
    transcript_id = f"PPSC_tx_{digest}"
    attributes = (
        f"candidate_group={report['group_id']};"
        f"candidate_category={report['category']};"
        f"support_source_count={report['support_source_count']};"
        f"support_sources={report['support_sources']};"
        f"support_group_count={report['support_group_count']};"
        f"support_groups={report['support_groups']};"
        f"Derives_from={model.query_id};miniprot_model={model.model_id};"
        f"Identity={model.identity:.4f};QueryCoverage={model.query_coverage:.4f}"
    )
    prefix = f"{model.seqid}\tPloidyPatchCandidate"
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
            f"ID=PPSC_exon_{digest}_{exon_number};Parent={transcript_id}\n"
        )
    for start, end, phase in model.cds:
        handle.write(
            f"{prefix}\tCDS\t{start}\t{end}\t.\t{model.strand}\t{phase}\t"
            f"Parent={transcript_id}\n"
        )


def _read_source_groups(
    path: str | Path | None,
    protein_map: dict[str, tuple[str, int]],
) -> dict[str, str]:
    sources = {source for source, _ in protein_map.values()}
    if path is None:
        return {source: source for source in sources}
    mapping: dict[str, str] = {}
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"source", "support_group"}
        if reader.fieldnames is None or not required <= set(reader.fieldnames):
            raise ValueError(
                "Source-group map requires source and support_group columns"
            )
        for line_number, row in enumerate(reader, start=2):
            source = row["source"].strip()
            support_group = row["support_group"].strip()
            if not source or not support_group:
                raise ValueError(
                    f"Empty source-group value at line {line_number}"
                )
            if source in mapping:
                raise ValueError(
                    f"Duplicate source in source-group map: {source}"
                )
            mapping[source] = support_group
    missing = sorted(sources - set(mapping))
    extra = sorted(set(mapping) - sources)
    if missing or extra:
        raise ValueError(
            "Source-group map must exactly cover protein-map sources; "
            f"missing={missing}, extra={extra}"
        )
    return mapping


def adapt_miniprot_structure_candidates(
    *,
    annotation_gff_path: str | Path,
    miniprot_gff_path: str | Path,
    protein_map_path: str | Path,
    output_gff_path: str | Path,
    decisions_tsv_path: str | Path,
    min_identity: float = 0.5,
    min_query_coverage: float = 0.5,
    min_source_support: int = 1,
    min_gene_overlap_fraction: float = 0.1,
    require_intact: bool = True,
    source_group_map_path: str | Path | None = None,
) -> dict[str, Any]:
    """Append exact-chain projection disagreements as structural candidates."""

    for name, value in (
        ("min_identity", min_identity),
        ("min_query_coverage", min_query_coverage),
        ("min_gene_overlap_fraction", min_gene_overlap_fraction),
    ):
        if not 0 <= value <= 1:
            raise ValueError(f"{name} must be between zero and one")
    if min_source_support < 1:
        raise ValueError("min_source_support must be positive")
    output = Path(output_gff_path)
    decisions = Path(decisions_tsv_path)
    manifest_path = Path(str(output) + ".manifest.json")
    collisions = [path for path in (output, decisions, manifest_path) if path.exists()]
    if collisions:
        raise FileExistsError(
            "Refusing to overwrite structure-candidate artifact(s): "
            + ", ".join(str(path) for path in collisions)
        )

    protein_map = _read_protein_map(protein_map_path)
    source_groups = _read_source_groups(source_group_map_path, protein_map)
    models = _parse_miniprot_models(miniprot_gff_path, protein_map)
    exact_chains, gene_bins, gene_intervals = _annotation_context(
        annotation_gff_path
    )
    existing_cds = _existing_cds_intervals(annotation_gff_path)

    preliminary_reasons: dict[str, str] = {}
    eligible_groups: dict[tuple[Any, ...], list[ProjectedModel]] = defaultdict(list)
    for model in models:
        reason = ""
        if not model.cds:
            reason = "no_cds"
        elif model.strand not in {"+", "-"}:
            reason = "invalid_strand"
        elif any(phase not in {"0", "1", "2"} for _, _, phase in model.cds):
            reason = "invalid_cds_phase"
        elif model.identity < min_identity:
            reason = "identity_below_threshold"
        elif model.query_coverage < min_query_coverage:
            reason = "query_coverage_below_threshold"
        elif require_intact and (model.frameshifts or model.stop_codons):
            reason = "non_intact_alignment"
        if reason:
            preliminary_reasons[model.model_id] = reason
        else:
            eligible_groups[_chain_key(model)].append(model)

    group_reports: dict[tuple[Any, ...], dict[str, Any]] = {}
    accepted: list[ProjectedModel] = []
    for key, support_models in eligible_groups.items():
        representative = min(
            support_models,
            key=lambda model: (
                -model.score,
                -model.identity,
                -model.query_coverage,
                model.rank,
                model.model_id,
            ),
        )
        sources = tuple(sorted({model.source for model in support_models}))
        groups = tuple(sorted({source_groups[source] for source in sources}))
        queries = tuple(sorted({model.query_id for model in support_models}))
        candidate_intervals = _model_intervals(representative)
        existing_overlap_bp = _overlap_with_index(
            candidate_intervals,
            existing_cds.get((representative.seqid, representative.strand)),
        )
        existing_overlap_fraction = existing_overlap_bp / representative.cds_bp
        genes = _overlapping_genes(
            representative,
            gene_bins=gene_bins,
            gene_intervals=gene_intervals,
            min_fraction=min_gene_overlap_fraction,
        )
        if key in exact_chains:
            status = "rejected"
            reason = "exact_cds_chain_present"
            category = "concordant"
        elif len(groups) < min_source_support:
            status = "rejected"
            reason = "source_support_below_threshold"
            category = "insufficient_source_support"
        else:
            status = "accepted"
            reason = "annotation_projection_disagreement"
            if not genes:
                category = "missing_annotation_candidate"
            elif len(genes) >= 2:
                category = "spans_multiple_annotated_genes"
            elif existing_overlap_fraction < 1.0:
                category = "cds_extension_or_missing_segment"
            else:
                category = "alternative_cds_chain_within_gene"
            accepted.append(representative)
        group_reports[key] = {
            "group_id": f"PPSC-{_chain_digest(key)}",
            "representative_model_id": representative.model_id,
            "support_model_count": len(support_models),
            "support_query_count": len(queries),
            "support_source_count": len(sources),
            "support_sources": ",".join(sources),
            "support_group_count": len(groups),
            "support_groups": ",".join(groups),
            "support_query_ids": ",".join(queries),
            "existing_cds_overlap_fraction": existing_overlap_fraction,
            "overlapping_gene_count": len(genes),
            "overlapping_gene_ids": ",".join(genes),
            "category": category,
            "status": status,
            "reason": reason,
        }

    accepted.sort(
        key=lambda model: (
            model.seqid,
            model.start,
            model.end,
            model.strand,
            model.model_id,
        )
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    with open_text(annotation_gff_path) as source_handle, output.open(
        "x", encoding="utf-8", newline=""
    ) as output_handle:
        last_line = ""
        for line in source_handle:
            output_handle.write(line)
            last_line = line
        if last_line and not last_line.endswith(("\n", "\r")):
            output_handle.write("\n")
        output_handle.write("###\n")
        for model in accepted:
            _write_structure_candidate(
                output_handle,
                model,
                group_reports[_chain_key(model)],
            )

    decision_fields = (
        "model_id",
        "query_id",
        "source",
        "seqid",
        "start",
        "end",
        "strand",
        "identity",
        "query_coverage",
        "group_id",
        "representative_model_id",
        "support_model_count",
        "support_query_count",
        "support_source_count",
        "support_sources",
        "support_group_count",
        "support_groups",
        "support_query_ids",
        "existing_cds_overlap_fraction",
        "overlapping_gene_count",
        "overlapping_gene_ids",
        "category",
        "status",
        "reason",
    )
    status_counts: Counter[str] = Counter()
    reason_counts: Counter[str] = Counter()
    category_counts: Counter[str] = Counter()
    decisions.parent.mkdir(parents=True, exist_ok=True)
    with decisions.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=decision_fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for model in models:
            preliminary_reason = preliminary_reasons.get(model.model_id)
            if preliminary_reason:
                report = {
                    "group_id": "",
                    "representative_model_id": "",
                    "support_model_count": 0,
                    "support_query_count": 0,
                    "support_source_count": 0,
                    "support_sources": "",
                    "support_group_count": 0,
                    "support_groups": "",
                    "support_query_ids": "",
                    "existing_cds_overlap_fraction": 0.0,
                    "overlapping_gene_count": 0,
                    "overlapping_gene_ids": "",
                    "category": "ineligible_projection",
                    "status": "rejected",
                    "reason": preliminary_reason,
                }
            else:
                report = dict(group_reports[_chain_key(model)])
                if (
                    report["status"] == "accepted"
                    and model.model_id != report["representative_model_id"]
                ):
                    report["status"] = "supporting"
                    report["reason"] = "exact_chain_support_member"
            status_counts[report["status"]] += 1
            reason_counts[report["reason"]] += 1
            category_counts[report["category"]] += 1
            writer.writerow(
                {
                    "model_id": model.model_id,
                    "query_id": model.query_id,
                    "source": model.source,
                    "seqid": model.seqid,
                    "start": model.start,
                    "end": model.end,
                    "strand": model.strand,
                    "identity": f"{model.identity:.6f}",
                    "query_coverage": f"{model.query_coverage:.6f}",
                    **{
                        field: (
                            f"{report[field]:.6f}"
                            if field == "existing_cds_overlap_fraction"
                            else report[field]
                        )
                        for field in decision_fields[9:]
                    },
                }
            )

    accepted_categories = Counter(
        group_reports[_chain_key(model)]["category"] for model in accepted
    )
    manifest = {
        "schema_version": STRUCTURE_CANDIDATE_SCHEMA_VERSION,
        "generator": {"name": "PloidyPatch", "version": __version__},
        "inputs": {
            "annotation_gff": {
                "file_name": Path(annotation_gff_path).name,
                "sha256": _file_sha256(annotation_gff_path),
            },
            "miniprot_gff": {
                "file_name": Path(miniprot_gff_path).name,
                "sha256": _file_sha256(miniprot_gff_path),
            },
            "protein_map": {
                "file_name": Path(protein_map_path).name,
                "sha256": _file_sha256(protein_map_path),
            },
            "source_group_map": (
                {
                    "file_name": Path(source_group_map_path).name,
                    "sha256": _file_sha256(source_group_map_path),
                }
                if source_group_map_path is not None
                else None
            ),
        },
        "parameters": {
            "min_identity": min_identity,
            "min_query_coverage": min_query_coverage,
            "min_source_support": min_source_support,
            "support_unit": "explicit_group" if source_group_map_path else "source",
            "min_gene_overlap_fraction": min_gene_overlap_fraction,
            "require_intact": require_intact,
            "candidate_rule": "eligible_exact_cds_chain_absent_from_annotation",
        },
        "counts": {
            "input_models": len(models),
            "eligible_chain_groups": len(eligible_groups),
            "accepted_chain_groups": len(accepted),
            "status_counts": dict(sorted(status_counts.items())),
            "reason_counts": dict(sorted(reason_counts.items())),
            "category_counts": dict(sorted(category_counts.items())),
            "accepted_category_counts": dict(sorted(accepted_categories.items())),
        },
        "outputs": {
            "candidate_gff": {
                "file_name": output.name,
                "sha256": _file_sha256(output),
            },
            "decisions": {
                "file_name": decisions.name,
                "rows": len(models),
                "sha256": _file_sha256(decisions),
            },
        },
    }
    with manifest_path.open("x", encoding="utf-8", newline="") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return manifest
