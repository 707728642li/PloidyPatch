from __future__ import annotations

import csv
import json
from collections import Counter
from pathlib import Path
from typing import Any, Iterable

from . import __version__
from .baseline import _file_sha256
from .consensus import SAFE_METHOD, _chain_digest
from .gff import TRANSCRIPT_TYPES, parse_attributes


COPY_FEATURE_SCHEMA_VERSION = "ploidypatch.copy_candidate_features.v1"
COPY_LABEL_SCHEMA_VERSION = "ploidypatch.copy_candidate_labels.v1"
SUPPORTED_METHODS = ("miniprot", "gemoma", "lifton")
FEATURE_COLUMNS = (
    "candidate_digest",
    "seqid",
    "start",
    "end",
    "strand",
    "span_bp",
    "cds_segments",
    "cds_bp",
    "support_method_count",
    "support_methods",
    "has_miniprot",
    "has_gemoma",
    "has_lifton",
    "miniprot_identity",
    "miniprot_query_coverage",
    "miniprot_rank",
    "miniprot_score",
    "miniprot_positive",
    "miniprot_frameshifts",
    "miniprot_stop_codons",
    "gemoma_pAA",
    "gemoma_iAA",
    "gemoma_score_ratio",
    "gemoma_complete_exon_count",
    "gemoma_reference_exon_count",
    "gemoma_start_complete",
    "gemoma_stop_complete",
    "gemoma_nps",
    "lifton_dna_identity",
    "lifton_protein_identity",
    "lifton_frameshift",
    "lifton_stop_missing",
    "lifton_start_missing",
    "lifton_inframe_indel",
    "wgd_existing_partner",
    "wgd_support_block_count",
    "wgd_longest_block_pairs",
)


def _parse_labeled_paths(
    values: Iterable[str],
) -> tuple[tuple[str, Path], ...]:
    parsed: list[tuple[str, Path]] = []
    observed: set[str] = set()
    for value in values:
        label, separator, raw_path = value.partition("=")
        if not separator or not label or not raw_path:
            raise ValueError("Each method decision input must be METHOD=PATH")
        if not SAFE_METHOD.fullmatch(label):
            raise ValueError(f"Unsafe method label: {label}")
        if label not in SUPPORTED_METHODS:
            raise ValueError(
                f"Unsupported copy-feature method {label!r}; expected one of "
                + ", ".join(SUPPORTED_METHODS)
            )
        if label in observed:
            raise ValueError(f"Duplicate method decision label: {label}")
        observed.add(label)
        parsed.append((label, Path(raw_path)))
    if not parsed:
        raise ValueError("At least one method decision input is required")
    return tuple(parsed)


def _read_tsv(path: Path, required: set[str]) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or not required <= set(reader.fieldnames):
            raise ValueError(
                f"{path} lacks required columns: {sorted(required)}"
            )
        return list(reader)


def _read_method_rows(path: Path) -> dict[str, dict[str, str]]:
    rows = _read_tsv(path, {"model_id", "status"})
    indexed: dict[str, dict[str, str]] = {}
    for line_number, row in enumerate(rows, start=2):
        model_id = row["model_id"]
        if not model_id or model_id in indexed:
            raise ValueError(
                f"Empty or duplicate method model_id in {path} line {line_number}"
            )
        indexed[model_id] = row
    return indexed


def _attribute_dict(row: dict[str, str]) -> dict[str, str]:
    raw = row.get("upstream_attributes", "")
    if not raw:
        return {}
    attributes, malformed = parse_attributes(raw)
    if malformed:
        raise ValueError("Malformed upstream_attributes in method decision row")
    return attributes


def _number(value: str | None, default: str = "") -> str:
    if value is None or value == "":
        return default
    try:
        return format(float(value), ".12g")
    except ValueError as exc:
        raise ValueError(f"Expected numeric feature value, observed {value!r}") from exc


def _ratio(numerator: str | None, denominator: str | None) -> str:
    if not numerator or not denominator:
        return ""
    try:
        top = float(numerator)
        bottom = float(denominator)
    except ValueError as exc:
        raise ValueError("Expected numeric ratio inputs") from exc
    return "" if bottom == 0 else format(top / bottom, ".12g")


def _method_features(method: str, row: dict[str, str]) -> dict[str, str]:
    if row.get("status") != "accepted":
        raise ValueError(
            f"Consensus links to non-accepted {method} model {row.get('model_id')!r}"
        )
    if method == "miniprot":
        return {
            "miniprot_identity": _number(row.get("identity")),
            "miniprot_query_coverage": _number(row.get("query_coverage")),
            "miniprot_rank": _number(row.get("rank")),
            "miniprot_score": _number(row.get("score")),
            "miniprot_positive": _number(row.get("positive")),
            "miniprot_frameshifts": _number(row.get("frameshifts"), "0"),
            "miniprot_stop_codons": _number(row.get("stop_codons"), "0"),
        }
    attributes = _attribute_dict(row)
    if method == "gemoma":
        return {
            "gemoma_pAA": _number(attributes.get("pAA")),
            "gemoma_iAA": _number(attributes.get("iAA")),
            "gemoma_score_ratio": _ratio(
                attributes.get("score"), attributes.get("bestScore")
            ),
            "gemoma_complete_exon_count": _number(attributes.get("ce")),
            "gemoma_reference_exon_count": _number(attributes.get("rce")),
            "gemoma_start_complete": str(int(attributes.get("start") == "M")),
            "gemoma_stop_complete": str(int(attributes.get("stop") == "*")),
            "gemoma_nps": _number(attributes.get("nps"), "0"),
        }
    mutations = {
        value.strip()
        for value in attributes.get("mutation", "").split(",")
        if value.strip()
    }
    return {
        "lifton_dna_identity": _number(attributes.get("dna_identity")),
        "lifton_protein_identity": _number(attributes.get("protein_identity")),
        "lifton_frameshift": str(int("frameshift" in mutations)),
        "lifton_stop_missing": str(int("stop_missing" in mutations)),
        "lifton_start_missing": str(int("start_missing" in mutations)),
        "lifton_inframe_indel": str(
            int(bool({"inframe_insertion", "inframe_deletion"} & mutations))
        ),
    }


def _read_wgd_rows(path: Path | None) -> dict[str, dict[str, str]]:
    if path is None:
        return {}
    rows = _read_tsv(
        path,
        {
            "consensus_digest",
            "support_block_count",
            "longest_block_pairs",
            "status",
        },
    )
    indexed: dict[str, dict[str, str]] = {}
    for line_number, row in enumerate(rows, start=2):
        digest = row["consensus_digest"]
        if not digest or digest in indexed:
            raise ValueError(
                f"Empty or duplicate WGD consensus_digest at line {line_number}"
            )
        indexed[digest] = row
    return indexed


def build_copy_candidate_features(
    *,
    consensus_decisions_tsv_path: str | Path,
    method_decision_inputs: Iterable[str],
    output_tsv_path: str | Path,
    wgd_selection_tsv_path: str | Path | None = None,
) -> dict[str, Any]:
    """Build truth-blind copy-candidate features from declared method families."""

    consensus_path = Path(consensus_decisions_tsv_path)
    output_path = Path(output_tsv_path)
    manifest_path = Path(str(output_path) + ".manifest.json")
    collisions = [path for path in (output_path, manifest_path) if path.exists()]
    if collisions:
        raise FileExistsError(
            "Refusing to overwrite copy-feature artifact(s): "
            + ", ".join(str(path) for path in collisions)
        )
    method_paths = _parse_labeled_paths(method_decision_inputs)
    method_rows = {
        method: _read_method_rows(path) for method, path in method_paths
    }
    wgd_path = Path(wgd_selection_tsv_path) if wgd_selection_tsv_path else None
    wgd_rows = _read_wgd_rows(wgd_path)
    consensus_rows = _read_tsv(
        consensus_path,
        {
            "consensus_digest",
            "seqid",
            "start",
            "end",
            "strand",
            "cds_segments",
            "cds_bp",
            "support_method_count",
            "support_methods",
            "upstream_models",
            "status",
        },
    )

    feature_rows: list[dict[str, str]] = []
    pattern_counts: Counter[str] = Counter()
    for row in consensus_rows:
        if row["status"] != "accepted":
            continue
        digest = row["consensus_digest"]
        methods = tuple(value for value in row["support_methods"].split(",") if value)
        if len(methods) != int(row["support_method_count"]):
            raise ValueError(f"Method support count mismatch for {digest}")
        if set(methods) - set(method_rows):
            raise ValueError(f"Missing method decision input for candidate {digest}")
        links: dict[str, str] = {}
        for value in row["upstream_models"].split(","):
            method, separator, model_id = value.partition(":")
            if not separator or not method or not model_id or method in links:
                raise ValueError(f"Malformed or duplicate upstream model link: {value}")
            links[method] = model_id
        if set(links) != set(methods):
            raise ValueError(f"Upstream model links disagree with methods for {digest}")
        start = int(row["start"])
        end = int(row["end"])
        output = {column: "" for column in FEATURE_COLUMNS}
        output.update(
            {
                "candidate_digest": digest,
                "seqid": row["seqid"],
                "start": str(start),
                "end": str(end),
                "strand": row["strand"],
                "span_bp": str(end - start + 1),
                "cds_segments": row["cds_segments"],
                "cds_bp": row["cds_bp"],
                "support_method_count": row["support_method_count"],
                "support_methods": row["support_methods"],
                **{
                    f"has_{method}": str(int(method in methods))
                    for method in SUPPORTED_METHODS
                },
            }
        )
        for method, model_id in links.items():
            decision = method_rows[method].get(model_id)
            if decision is None:
                raise ValueError(
                    f"Unknown {method} upstream model {model_id!r} for {digest}"
                )
            output.update(_method_features(method, decision))
        wgd = wgd_rows.get(digest)
        output["wgd_existing_partner"] = str(
            int(wgd is not None and wgd["status"] == "accepted")
        )
        output["wgd_support_block_count"] = _number(
            None if wgd is None else wgd["support_block_count"], "0"
        )
        output["wgd_longest_block_pairs"] = _number(
            None if wgd is None else wgd["longest_block_pairs"], "0"
        )
        feature_rows.append(output)
        pattern_counts[row["support_methods"]] += 1

    if wgd_rows:
        missing_wgd = {row["candidate_digest"] for row in feature_rows} - set(wgd_rows)
        if missing_wgd:
            raise ValueError(
                "WGD selection lacks candidate digest(s): "
                + ", ".join(sorted(missing_wgd)[:10])
            )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=FEATURE_COLUMNS, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(feature_rows)

    manifest: dict[str, Any] = {
        "schema_version": COPY_FEATURE_SCHEMA_VERSION,
        "generator": {"name": "PloidyPatch", "version": __version__},
        "truth_access": False,
        "inputs": {
            "consensus_decisions": {
                "file_name": consensus_path.name,
                "sha256": _file_sha256(consensus_path),
            },
            "method_decisions": {
                method: {"file_name": path.name, "sha256": _file_sha256(path)}
                for method, path in method_paths
            },
            "wgd_selection": (
                None
                if wgd_path is None
                else {"file_name": wgd_path.name, "sha256": _file_sha256(wgd_path)}
            ),
        },
        "feature_columns": list(FEATURE_COLUMNS),
        "counts": {
            "accepted_candidates": len(feature_rows),
            "support_patterns": dict(sorted(pattern_counts.items())),
            "wgd_partner_candidates": sum(
                int(row["wgd_existing_partner"]) for row in feature_rows
            ),
        },
        "outputs": {
            "features": {
                "file_name": output_path.name,
                "sha256": _file_sha256(output_path),
                "rows": len(feature_rows),
            }
        },
    }
    with manifest_path.open("x", encoding="utf-8", newline="") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return manifest


def _truth_chain_index(
    truth: dict[str, Any],
) -> tuple[dict[str, tuple[str, str]], list[dict[str, Any]]]:
    digest_to_event: dict[str, tuple[str, str]] = {}
    loci: list[dict[str, Any]] = []
    for event in truth.get("events", []):
        event_id = event.get("event_id", "")
        target = event.get("target", {})
        gene_id = target.get("gene_id", "")
        loci.append(
            {
                "event_id": event_id,
                "gene_id": gene_id,
                "seqid": target.get("seqid", ""),
                "start": int(target.get("start", 0)),
                "end": int(target.get("end", 0)),
                "strand": target.get("strand", ""),
            }
        )
        transcripts: dict[str, tuple[str, str]] = {}
        cds_by_parent: dict[str, list[tuple[int, int, str]]] = {}
        for edit in event.get("line_edits", []):
            raw_line = edit.get("source_raw_line", "").rstrip("\r\n")
            fields = raw_line.split("\t")
            if len(fields) != 9:
                raise ValueError(f"Malformed hidden source GFF line in {event_id}")
            attributes, malformed = parse_attributes(fields[8])
            if malformed:
                raise ValueError(f"Malformed hidden source attributes in {event_id}")
            if fields[2] in TRANSCRIPT_TYPES:
                transcript_id = attributes.get("ID", "")
                if transcript_id:
                    transcripts[transcript_id] = (fields[0], fields[6])
            elif fields[2] == "CDS":
                parents = [
                    value
                    for value in attributes.get("Parent", "").split(",")
                    if value
                ]
                if len(parents) != 1:
                    raise ValueError(f"Hidden CDS requires one parent in {event_id}")
                cds_by_parent.setdefault(parents[0], []).append(
                    (int(fields[3]), int(fields[4]), fields[7])
                )
        for transcript_id, cds in cds_by_parent.items():
            if transcript_id not in transcripts:
                raise ValueError(f"Hidden CDS parent is absent in {event_id}")
            seqid, strand = transcripts[transcript_id]
            key = (seqid, strand, tuple(sorted(cds)))
            digest = _chain_digest(key)
            if digest in digest_to_event:
                raise ValueError(f"Hidden CDS digest occurs in multiple events: {digest}")
            digest_to_event[digest] = (event_id, gene_id)
    return digest_to_event, loci


def label_copy_candidate_features(
    *,
    feature_tsv_path: str | Path,
    hidden_truth_json_path: str | Path,
    output_tsv_path: str | Path,
) -> dict[str, Any]:
    """Attach evaluator-only exact-CDS labels to a truth-blind feature table."""

    feature_path = Path(feature_tsv_path)
    truth_path = Path(hidden_truth_json_path)
    output_path = Path(output_tsv_path)
    manifest_path = Path(str(output_path) + ".manifest.json")
    collisions = [path for path in (output_path, manifest_path) if path.exists()]
    if collisions:
        raise FileExistsError(
            "Refusing to overwrite copy-label artifact(s): "
            + ", ".join(str(path) for path in collisions)
        )
    feature_rows = _read_tsv(feature_path, set(FEATURE_COLUMNS))
    with truth_path.open("r", encoding="utf-8") as handle:
        truth = json.load(handle)
    digest_to_event, loci = _truth_chain_index(truth)
    label_columns = (
        *FEATURE_COLUMNS,
        "label_exact_cds",
        "truth_event_id",
        "truth_gene_id",
        "truth_locus_overlap_bp",
        "truth_locus_overlap_fraction_shorter",
    )
    output_rows: list[dict[str, str]] = []
    matched_truth: set[str] = set()
    near_miss_candidates = 0
    for row in feature_rows:
        digest = row["candidate_digest"]
        exact = digest_to_event.get(digest)
        best_overlap = 0
        best_denominator = 0
        start = int(row["start"])
        end = int(row["end"])
        for locus in loci:
            if row["seqid"] != locus["seqid"] or row["strand"] != locus["strand"]:
                continue
            overlap = max(
                0, min(end, locus["end"]) - max(start, locus["start"]) + 1
            )
            if overlap > best_overlap:
                best_overlap = overlap
                best_denominator = min(
                    end - start + 1, locus["end"] - locus["start"] + 1
                )
        if exact is not None:
            matched_truth.add(digest)
        elif best_overlap:
            near_miss_candidates += 1
        output_rows.append(
            {
                **row,
                "label_exact_cds": str(int(exact is not None)),
                "truth_event_id": "" if exact is None else exact[0],
                "truth_gene_id": "" if exact is None else exact[1],
                "truth_locus_overlap_bp": str(best_overlap),
                "truth_locus_overlap_fraction_shorter": (
                    "0"
                    if not best_denominator
                    else format(best_overlap / best_denominator, ".12g")
                ),
            }
        )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=label_columns, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(output_rows)
    manifest: dict[str, Any] = {
        "schema_version": COPY_LABEL_SCHEMA_VERSION,
        "generator": {"name": "PloidyPatch", "version": __version__},
        "evaluator_only": True,
        "inputs": {
            "features": {
                "file_name": feature_path.name,
                "sha256": _file_sha256(feature_path),
            },
            "hidden_truth": {
                "file_name": truth_path.name,
                "sha256": _file_sha256(truth_path),
            },
        },
        "counts": {
            "candidates": len(output_rows),
            "positive_exact_cds": len(matched_truth),
            "negative_candidates": len(output_rows) - len(matched_truth),
            "nonexact_candidates_overlapping_hidden_loci": near_miss_candidates,
            "hidden_cds_chains": len(digest_to_event),
            "hidden_cds_chains_matched": len(matched_truth),
        },
        "outputs": {
            "labels": {
                "file_name": output_path.name,
                "sha256": _file_sha256(output_path),
                "rows": len(output_rows),
            }
        },
    }
    with manifest_path.open("x", encoding="utf-8", newline="") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return manifest
