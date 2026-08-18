from __future__ import annotations

import csv
import datetime as dt
import hashlib
import json
import os
import shutil
import tempfile
from collections import Counter
from pathlib import Path
from typing import Any

from . import __version__
from .baseline import _file_sha256
from .gff import parse_attributes
from .io import open_text
from .wgd_candidate_select import _candidate_blocks


COPY_ADDITION_EDIT_SCHEMA = "ploidypatch.copy_addition_patch_edits.v1"
REVIEW_DECISION_SCHEMA = "ploidypatch.copy_review_decisions.v1"
POOL_SCHEMA = "ploidypatch.method_candidate_pool.v2"
REVIEW_FIELDS = (
    "candidate_digest",
    "review_decision",
    "reviewer",
    "reviewed_at_utc",
    "rationale",
)
REVIEW_DECISIONS = frozenset({"accept_addition", "reject", "defer"})


def _write_json_atomic_no_overwrite(
    payload: dict[str, Any], output: Path, *, artifact_name: str
) -> None:
    """Write JSON through a visible working file and refuse stale outputs."""

    working = Path(f"{output}.working")
    if output.exists() or working.exists():
        raise FileExistsError(f"Refusing to overwrite {artifact_name}: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    try:
        with working.open("x", encoding="utf-8", newline="") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        if output.exists():
            raise FileExistsError(f"Refusing to overwrite {artifact_name}: {output}")
        working.rename(output)
    except Exception:
        if working.exists():
            working.unlink()
        raise


def _feature_ids(path: Path) -> set[str]:
    feature_ids: set[str] = set()
    with open_text(path) as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            stripped = raw_line.rstrip("\r\n")
            if not stripped or stripped.startswith("#"):
                continue
            fields = stripped.split("\t")
            if len(fields) != 9:
                raise ValueError(f"Malformed annotation line {line_number}")
            attributes, malformed = parse_attributes(fields[8])
            if malformed:
                raise ValueError(f"Malformed annotation attributes at line {line_number}")
            feature_id = attributes.get("ID", "")
            if feature_id:
                if feature_id in feature_ids and fields[2] in {
                    "gene",
                    "mRNA",
                    "transcript",
                    "exon",
                }:
                    raise ValueError(f"Duplicate patch-relevant feature ID: {feature_id}")
                feature_ids.add(feature_id)
    return feature_ids


def compile_copy_addition_patch_edits(
    *,
    annotation_gff_path: str | Path,
    candidate_gff_path: str | Path,
    output_edits_json_path: str | Path,
) -> dict[str, Any]:
    """Compile reviewed coding-copy additions into one reversible EOF edit."""

    annotation_path = Path(annotation_gff_path)
    candidate_path = Path(candidate_gff_path)
    output_path = Path(output_edits_json_path)
    if output_path.exists() or Path(f"{output_path}.working").exists():
        raise FileExistsError(f"Refusing to overwrite copy-addition edits: {output_path}")
    for required in (annotation_path, candidate_path):
        if not required.is_file() or required.stat().st_size == 0:
            raise ValueError(f"Missing or empty copy-addition input: {required}")
        if required.suffix == ".gz":
            raise ValueError("Copy-addition compilation requires uncompressed GFF3")

    blocks, digests = _candidate_blocks(candidate_path, annotation_path)
    source_ids = _feature_ids(annotation_path)
    appended_ids: set[str] = set()
    for gene_id, lines in blocks.items():
        for raw_line in lines:
            fields = raw_line.rstrip("\r\n").split("\t")
            attributes, malformed = parse_attributes(fields[8])
            if malformed:
                raise ValueError(f"Malformed candidate block for {gene_id}")
            feature_id = attributes.get("ID", "")
            if feature_id:
                if feature_id in appended_ids:
                    raise ValueError(f"Duplicate appended feature ID: {feature_id}")
                appended_ids.add(feature_id)
    collisions = source_ids & appended_ids
    if collisions:
        raise ValueError(
            "Copy-addition feature IDs collide with the annotation: "
            + ", ".join(sorted(collisions)[:10])
        )

    with open_text(annotation_path) as handle:
        source_lines = list(handle)
    if not source_lines:
        raise ValueError("Annotation GFF3 is empty")
    last_line_number = len(source_lines)
    last_source_line = source_lines[-1]
    retained_last_line = last_source_line
    if not retained_last_line.endswith(("\n", "\r")):
        retained_last_line += "\n"
    ordered_genes = sorted(blocks, key=lambda gene_id: digests[gene_id])
    replacement_lines = [retained_last_line, "###\n"]
    event_ids: list[str] = []
    events: list[dict[str, Any]] = []
    for gene_id in ordered_genes:
        digest = digests[gene_id]
        event_id = f"PPCOPY-{digest[:20]}"
        event_ids.append(event_id)
        replacement_lines.extend(blocks[gene_id])
        events.append(
            {
                "event_id": event_id,
                "event_type": "annotation_copy_collapse",
                "candidate_gene_id": gene_id,
                "consensus_digest": digest,
                "feature_lines": len(blocks[gene_id]),
            }
        )

    payload: dict[str, Any] = {
        "schema_version": COPY_ADDITION_EDIT_SCHEMA,
        "generator": {"name": "PloidyPatch", "version": __version__},
        "inputs": {
            "annotation_gff": {
                "file_name": annotation_path.name,
                "sha256": _file_sha256(annotation_path),
            },
            "candidate_gff": {
                "file_name": candidate_path.name,
                "sha256": _file_sha256(candidate_path),
            },
        },
        "policy": {
            "operation": "append_reviewed_complete_gene_hierarchies",
            "existing_line_edits": 0,
            "existing_feature_deletions": 0,
            "output_scope": "coding_model_only",
            "automatic_approval": False,
        },
        "event_ids": event_ids,
        "events": events,
        "operations": [
            {
                "source_line_number": last_line_number,
                "source_raw_line": last_source_line,
                "replacement_lines": replacement_lines,
            }
        ],
        "counts": {
            "copy_addition_events": len(events),
            "appended_feature_lines": sum(
                len(blocks[gene_id]) for gene_id in ordered_genes
            ),
            "source_feature_ids": len(source_ids),
            "appended_feature_ids": len(appended_ids),
        },
    }
    _write_json_atomic_no_overwrite(
        payload, output_path, artifact_name="copy-addition edits"
    )
    return payload


def _load_json_object(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"JSON root is not an object: {path}")
    return value


def _read_pool_decisions(path: Path) -> dict[str, dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "consensus_digest",
            "status",
            "conflict_set_digest",
            "conflict_member_count",
        }
        if reader.fieldnames is None or not required <= set(reader.fieldnames):
            raise ValueError("Pool decisions lack review provenance fields")
        forbidden = [
            field
            for field in reader.fieldnames
            if "truth" in field.lower()
            or "label" in field.lower()
            or field.lower().startswith("v04_")
        ]
        if forbidden:
            raise ValueError("Pool decisions contain evaluator or reserved fields")
        rows = list(reader)
    indexed = {row["consensus_digest"]: row for row in rows}
    if len(indexed) != len(rows) or "" in indexed:
        raise ValueError("Pool decisions contain empty or duplicate digests")
    return indexed


def _read_review_decisions(path: Path) -> dict[str, dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != REVIEW_FIELDS:
            raise ValueError(
                "Review decisions must use the exact v1 fields: "
                + ",".join(REVIEW_FIELDS)
            )
        rows = list(reader)
    indexed: dict[str, dict[str, str]] = {}
    for row_number, row in enumerate(rows, start=2):
        digest = row["candidate_digest"]
        if not digest or digest in indexed:
            raise ValueError(
                f"Review decisions contain empty/duplicate digest at line {row_number}"
            )
        if row["review_decision"] not in REVIEW_DECISIONS:
            raise ValueError(f"Invalid review decision at line {row_number}")
        if not row["reviewer"].strip() or not row["rationale"].strip():
            raise ValueError(f"Review provenance is incomplete at line {row_number}")
        timestamp = row["reviewed_at_utc"]
        try:
            parsed = dt.datetime.fromisoformat(timestamp.replace("Z", "+00:00"))
        except ValueError as exc:
            raise ValueError(
                f"Invalid reviewed_at_utc at line {row_number}"
            ) from exc
        if parsed.tzinfo is None or parsed.utcoffset() is None:
            raise ValueError(f"Review timestamp lacks timezone at line {row_number}")
        if parsed.utcoffset() != dt.timedelta(0):
            raise ValueError(f"Review timestamp is not UTC at line {row_number}")
        indexed[digest] = row
    if not indexed:
        raise ValueError("Review decisions are empty")
    return indexed


def compile_reviewed_copy_addition_patch_edits(
    *,
    annotation_gff_path: str | Path,
    candidate_gff_path: str | Path,
    pool_decisions_tsv_path: str | Path,
    pool_manifest_json_path: str | Path,
    review_decisions_tsv_path: str | Path,
    output_edits_json_path: str | Path,
) -> dict[str, Any]:
    """Compile only explicit human-accepted copy additions into patch edits."""

    annotation = Path(annotation_gff_path)
    candidate = Path(candidate_gff_path)
    pool_decisions = Path(pool_decisions_tsv_path)
    pool_manifest = Path(pool_manifest_json_path)
    review_decisions = Path(review_decisions_tsv_path)
    output = Path(output_edits_json_path)
    for required in (
        annotation,
        candidate,
        pool_decisions,
        pool_manifest,
        review_decisions,
    ):
        if not required.is_file() or required.stat().st_size == 0:
            raise ValueError(f"Missing or empty reviewed-copy input: {required}")
    if output.exists() or Path(f"{output}.working").exists():
        raise FileExistsError(f"Refusing to overwrite reviewed copy edits: {output}")
    if annotation.suffix == ".gz" or candidate.suffix == ".gz":
        raise ValueError("Reviewed copy compilation requires uncompressed GFF3")

    manifest = _load_json_object(pool_manifest)
    decisions = _read_pool_decisions(pool_decisions)
    if (
        manifest.get("schema_version") != POOL_SCHEMA
        or manifest.get("outputs", {}).get("candidate_gff", {}).get("sha256")
        != _file_sha256(candidate)
        or manifest.get("outputs", {}).get("decisions", {}).get("sha256")
        != _file_sha256(pool_decisions)
    ):
        raise ValueError("Candidate pool manifest fails schema or checksum gate")
    parameters = manifest.get("parameters", {})
    if (
        parameters.get("redundancy_policy") != "retain_distinct_chains"
        or parameters.get("conflict_action")
        != "retain_all_for_ranking_and_review"
        or parameters.get("min_method_support") != 1
        or parameters.get("max_redundancy_overlap") != 0.5
    ):
        raise ValueError("Candidate pool is not the review-only chain-preserving policy")
    accepted_pool = {
        digest: row for digest, row in decisions.items() if row["status"] == "accepted"
    }
    if int(manifest.get("counts", {}).get("accepted_models", -1)) != len(
        accepted_pool
    ):
        raise ValueError("Candidate pool accepted count differs")
    if int(
        manifest.get("outputs", {}).get("decisions", {}).get("rows", -1)
    ) != len(decisions):
        raise ValueError("Candidate pool decision row count differs")
    members_by_conflict: dict[str, list[str]] = {}
    for digest, row in accepted_pool.items():
        conflict = row["conflict_set_digest"]
        declared = int(row["conflict_member_count"])
        if conflict:
            members_by_conflict.setdefault(conflict, []).append(digest)
        elif declared != 1:
            raise ValueError("Non-conflict candidate has non-singleton member count")
    for conflict, members in members_by_conflict.items():
        if any(
            int(accepted_pool[digest]["conflict_member_count"]) != len(members)
            for digest in members
        ):
            raise ValueError(f"Conflict member count differs: {conflict}")
    counts = manifest.get("counts", {})
    if (
        int(counts.get("conflict_sets", -1)) != len(members_by_conflict)
        or int(counts.get("conflicted_chains", -1))
        != sum(len(members) for members in members_by_conflict.values())
    ):
        raise ValueError("Candidate pool conflict counts differ")

    blocks, digest_by_gene = _candidate_blocks(candidate, annotation)
    gene_by_digest = {digest: gene for gene, digest in digest_by_gene.items()}
    if len(gene_by_digest) != len(digest_by_gene) or set(gene_by_digest) != set(
        accepted_pool
    ):
        raise ValueError("Candidate GFF and accepted pool universe differ")
    reviews = _read_review_decisions(review_decisions)
    unknown = set(reviews) - set(accepted_pool)
    if unknown:
        raise ValueError(
            "Review decisions contain unknown candidates: "
            + ", ".join(sorted(unknown)[:10])
        )
    accepted_reviews = {
        digest: row
        for digest, row in reviews.items()
        if row["review_decision"] == "accept_addition"
    }
    if not accepted_reviews:
        raise ValueError("No candidate has explicit accept_addition review")
    accepted_by_conflict: dict[str, list[str]] = {}
    for digest in accepted_reviews:
        conflict = accepted_pool[digest]["conflict_set_digest"]
        if conflict:
            accepted_by_conflict.setdefault(conflict, []).append(digest)
    collisions = {
        conflict: members
        for conflict, members in accepted_by_conflict.items()
        if len(members) > 1
    }
    if collisions:
        raise ValueError(
            "Review accepts mutually exclusive conflict alternatives: "
            + ", ".join(sorted(collisions)[:10])
        )

    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="ploidypatch-reviewed-copy-") as temp:
        temporary = Path(temp)
        selected_candidate = temporary / "selected_reviewed_candidates.gff3"
        with annotation.open("rb") as source, selected_candidate.open("xb") as target:
            shutil.copyfileobj(source, target, length=8 * 1024 * 1024)
            source.seek(0, 2)
            if source.tell():
                source.seek(-1, 2)
                if source.read(1) not in {b"\n", b"\r"}:
                    target.write(b"\n")
            target.write(b"###\n")
            for digest in sorted(accepted_reviews):
                gene_id = gene_by_digest[digest]
                target.write("".join(blocks[gene_id]).encode("utf-8"))
        temporary_edits = temporary / "edits.json"
        payload = compile_copy_addition_patch_edits(
            annotation_gff_path=annotation,
            candidate_gff_path=selected_candidate,
            output_edits_json_path=temporary_edits,
        )
        selected_candidate_sha256 = _file_sha256(selected_candidate)

    payload["inputs"]["candidate_gff"] = {
        "file_name": candidate.name,
        "sha256": _file_sha256(candidate),
    }
    payload["inputs"].update(
        {
            "pool_decisions": {
                "file_name": pool_decisions.name,
                "sha256": _file_sha256(pool_decisions),
            },
            "pool_manifest": {
                "file_name": pool_manifest.name,
                "sha256": _file_sha256(pool_manifest),
            },
            "review_decisions": {
                "file_name": review_decisions.name,
                "sha256": _file_sha256(review_decisions),
                "schema_version": REVIEW_DECISION_SCHEMA,
            },
            "selected_candidate_gff_sha256": selected_candidate_sha256,
        }
    )
    payload["policy"].update(
        {
            "candidate_selection": "explicit_human_accept_addition_only",
            "unreviewed_policy": "exclude",
            "rejected_or_deferred_policy": "exclude",
            "conflict_policy": "at_most_one_accepted_per_conflict_set",
            "review_decision_schema": REVIEW_DECISION_SCHEMA,
        }
    )
    for event in payload["events"]:
        review = accepted_reviews[event["consensus_digest"]]
        event["review"] = {
            "decision": review["review_decision"],
            "reviewer": review["reviewer"],
            "reviewed_at_utc": review["reviewed_at_utc"],
            "rationale": review["rationale"],
        }
    decision_counts = Counter(row["review_decision"] for row in reviews.values())
    payload["counts"].update(
        {
            "pool_candidates": len(accepted_pool),
            "review_rows": len(reviews),
            "review_decision_counts": dict(sorted(decision_counts.items())),
            "unreviewed_candidates": len(accepted_pool) - len(reviews),
        }
    )
    _write_json_atomic_no_overwrite(
        payload, output, artifact_name="reviewed copy edits"
    )
    return payload
