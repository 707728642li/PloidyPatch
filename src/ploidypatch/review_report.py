from __future__ import annotations

import csv
import datetime as dt
import html
from importlib import resources
import json
import math
import os
import re
from collections import Counter, defaultdict
from pathlib import Path
import shutil
from typing import Any

from . import __version__
from .artifact_manifest import sha256_file, verify_sha256sums, write_sha256sums
from .gff import parse_attributes
from .io import open_text


REPORT_SCHEMA_VERSION = "ploidypatch.review_report.v1"
POOL_SCHEMA_VERSION = "ploidypatch.method_candidate_pool.v2"
REVIEW_FIELDS = (
    "candidate_digest",
    "review_decision",
    "reviewer",
    "reviewed_at_utc",
    "rationale",
)
REVIEW_DECISIONS = frozenset({"accept_addition", "reject", "defer"})
METHODS = ("miniprot", "gemoma", "lifton")
SCORE_FIELDS = (
    "v04_primary_rank_score",
    "homeolog_topology_rank_score",
    "v03_primary_rank_score",
    "model_calibrated_probability",
)
SCORE_DIRECTION = "higher_is_higher_review_priority"
PERCENTILE_FIELDS = (
    "v04_primary_rank_percentile",
    "homeolog_topology_rank_percentile",
)


def _regular(path: Path, label: str) -> None:
    if not path.is_file() or path.is_symlink() or path.stat().st_size == 0:
        raise ValueError(f"{label} must be a non-empty regular file: {path}")


def _read_json(path: Path, label: str) -> dict[str, Any]:
    _regular(path, label)
    with path.open(encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"{label} JSON root is not an object")
    return value


def _read_tsv(
    path: Path,
    *,
    label: str,
    required: set[str],
    exact: tuple[str, ...] | None = None,
    reject_truth: bool = False,
) -> tuple[tuple[str, ...], list[dict[str, str]]]:
    _regular(path, label)
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = tuple(reader.fieldnames or ())
        if exact is not None and fields != exact:
            raise ValueError(f"{label} fields differ from the exact schema")
        if not required <= set(fields):
            raise ValueError(f"{label} lacks required fields: {sorted(required - set(fields))}")
        if reject_truth:
            exact_forbidden = {
                "answer",
                "exact_positive",
                "gold",
                "gold_label",
                "ground_truth",
                "is_positive",
                "label",
                "oracle_label",
                "truth",
                "y_true",
            }
            forbidden = [
                field
                for field in fields
                if field.lower() in exact_forbidden
                or field.lower().startswith(("label_", "truth_", "gold_", "oracle_"))
                or "ground_truth" in field.lower()
                or field.lower().endswith("_truth")
            ]
            if forbidden:
                raise ValueError(
                    f"{label} contains evaluator-only fields: {', '.join(forbidden)}"
                )
        rows = list(reader)
    return fields, rows


def _candidate_models(path: Path) -> dict[str, dict[str, Any]]:
    _regular(path, "candidate GFF3")
    rows: list[tuple[str, str, int, int, str, dict[str, str]]] = []
    candidates: dict[str, dict[str, Any]] = {}
    gene_by_id: dict[str, str] = {}
    transcript_to_gene: dict[str, str] = {}
    with open_text(path) as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\r\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) != 9:
                raise ValueError(f"Malformed candidate GFF3 line {line_number}")
            seqid, _, feature_type, start_text, end_text, _, strand, phase, raw_attrs = fields
            attrs, malformed = parse_attributes(raw_attrs)
            if malformed:
                raise ValueError(f"Malformed candidate attributes at line {line_number}")
            try:
                start = int(start_text)
                end = int(end_text)
            except ValueError as exc:
                raise ValueError(f"Invalid candidate coordinates at line {line_number}") from exc
            rows.append((feature_type, seqid, start, end, strand, attrs | {"_phase": phase}))
            if feature_type == "gene" and attrs.get("consensus_digest"):
                digest = attrs["consensus_digest"]
                gene_id = attrs.get("ID", "")
                if not gene_id or not digest or digest in candidates or gene_id in gene_by_id:
                    raise ValueError("Candidate GFF3 contains empty or duplicate candidate identifiers")
                methods = tuple(
                    method
                    for method in attrs.get("support_methods", "").split(",")
                    if method
                )
                if any(method not in METHODS for method in methods):
                    raise ValueError(f"Unsupported candidate method label for {digest}")
                candidates[digest] = {
                    "candidate_digest": digest,
                    "gene_id": gene_id,
                    "seqid": seqid,
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "span_bp": end - start + 1,
                    "support_methods": list(methods),
                    "transcripts": [],
                }
                gene_by_id[gene_id] = digest

    for feature_type, seqid, start, end, strand, attrs in rows:
        if feature_type in {"mRNA", "transcript"}:
            parent = attrs.get("Parent", "")
            transcript_id = attrs.get("ID", "")
            if parent in gene_by_id and transcript_id:
                digest = gene_by_id[parent]
                transcript_to_gene[transcript_id] = digest
                candidates[digest]["transcripts"].append(
                    {
                        "transcript_id": transcript_id,
                        "start": start,
                        "end": end,
                        "cds": [],
                    }
                )

    transcript_index: dict[tuple[str, str], dict[str, Any]] = {}
    for digest, candidate in candidates.items():
        for transcript in candidate["transcripts"]:
            transcript_index[(digest, transcript["transcript_id"])] = transcript
    for feature_type, seqid, start, end, strand, attrs in rows:
        if feature_type != "CDS":
            continue
        for parent in attrs.get("Parent", "").split(","):
            digest = transcript_to_gene.get(parent)
            if digest is not None:
                transcript_index[(digest, parent)]["cds"].append(
                    [start, end, attrs.get("_phase", ".")]
                )
    for candidate in candidates.values():
        for transcript in candidate["transcripts"]:
            transcript["cds"].sort(key=lambda item: item[0])
        candidate["transcripts"].sort(
            key=lambda item: (-sum(end - start + 1 for start, end, _ in item["cds"]), item["transcript_id"])
        )
    if not candidates:
        raise ValueError("Candidate GFF3 contains no consensus_digest gene models")
    return candidates


def _index_rows(
    rows: list[dict[str, str]], key: str, label: str
) -> dict[str, dict[str, str]]:
    indexed: dict[str, dict[str, str]] = {}
    for row_number, row in enumerate(rows, start=2):
        value = row.get(key, "")
        if not value or value in indexed:
            raise ValueError(f"{label} contains empty or duplicate {key} at row {row_number}")
        indexed[value] = row
    return indexed


def _validate_reviews(reviews: dict[str, dict[str, str]]) -> None:
    for row_number, row in enumerate(reviews.values(), start=2):
        if row["review_decision"] not in REVIEW_DECISIONS:
            raise ValueError(f"Review decisions contain an unsupported decision at row {row_number}")
        if not row["reviewer"].strip() or not row["rationale"].strip():
            raise ValueError(f"Review provenance is incomplete at row {row_number}")
        timestamp = row["reviewed_at_utc"]
        try:
            parsed = dt.datetime.fromisoformat(timestamp.replace("Z", "+00:00"))
        except ValueError as exc:
            raise ValueError(f"Invalid reviewed_at_utc at row {row_number}") from exc
        if parsed.tzinfo is None or parsed.utcoffset() != dt.timedelta(0):
            raise ValueError(f"Review timestamp is not UTC at row {row_number}")


def _validate_conflicts(
    *,
    accepted_pool: dict[str, dict[str, str]],
    manifest: dict[str, Any],
) -> dict[str, list[str]]:
    members_by_conflict: dict[str, list[str]] = defaultdict(list)
    for digest, row in accepted_pool.items():
        try:
            declared = int(row["conflict_member_count"])
        except ValueError as exc:
            raise ValueError(f"Invalid conflict member count for {digest}") from exc
        conflict = row["conflict_set_digest"]
        if conflict:
            members_by_conflict[conflict].append(digest)
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
        raise ValueError("Candidate-pool manifest conflict counts differ")
    return members_by_conflict


def _optional_candidate_table(
    path: Path | None, *, label: str
) -> tuple[tuple[str, ...], dict[str, dict[str, str]]]:
    if path is None:
        return (), {}
    fields, rows = _read_tsv(
        path,
        label=label,
        required={"candidate_digest"},
        reject_truth=True,
    )
    return fields, _index_rows(rows, "candidate_digest", label)


def _float(row: dict[str, str], names: tuple[str, ...]) -> tuple[str, float | None]:
    for name in names:
        value = row.get(name, "")
        if value != "":
            try:
                parsed = float(value)
                if not math.isfinite(parsed):
                    raise ValueError(f"Non-finite {name}: {value!r}")
                return name, parsed
            except ValueError as exc:
                raise ValueError(f"Non-numeric {name}: {value!r}") from exc
    return "", None


def _input_record(path: Path, role: str) -> dict[str, Any]:
    return {
        "role": role,
        "file_name": path.name,
        "bytes": path.stat().st_size,
        "sha256": sha256_file(path),
    }


def _score_histogram(values: list[float], bins: int = 10) -> dict[str, Any] | None:
    if not values:
        return None
    low = min(values)
    high = max(values)
    if low == high:
        return {
            "minimum": low,
            "maximum": high,
            "bins": [{"lower": low, "upper": high, "count": len(values)}],
        }
    width = (high - low) / bins
    counts = [0] * bins
    for value in values:
        index = min(int((value - low) / width), bins - 1)
        counts[index] += 1
    return {
        "minimum": low,
        "maximum": high,
        "bins": [
            {
                "lower": low + index * width,
                "upper": low + (index + 1) * width,
                "count": count,
            }
            for index, count in enumerate(counts)
        ],
    }


def _plural(count: int, singular: str, plural: str | None = None) -> str:
    return f"{count:,} {singular if count == 1 else (plural or singular + 's')}"


def _generated_at_utc() -> str:
    epoch = os.environ.get("SOURCE_DATE_EPOCH")
    if epoch is None:
        moment = dt.datetime.now(dt.timezone.utc)
    else:
        try:
            moment = dt.datetime.fromtimestamp(int(epoch), tz=dt.timezone.utc)
        except (ValueError, OverflowError, OSError) as exc:
            raise ValueError(f"Invalid SOURCE_DATE_EPOCH: {epoch!r}") from exc
    return moment.replace(microsecond=0).isoformat().replace("+00:00", "Z")


def _primary_chain(candidate: dict[str, Any]) -> tuple[tuple[int, int, str], ...]:
    if not candidate["transcripts"]:
        return ()
    return tuple(
        (int(start), int(end), str(phase))
        for start, end, phase in candidate["transcripts"][0]["cds"]
    )


def _chain_difference(
    left: dict[str, Any], right: dict[str, Any]
) -> dict[str, Any]:
    left_chain = _primary_chain(left)
    right_chain = _primary_chain(right)
    same_count = len(left_chain) == len(right_chain)
    displacement = None
    phase_mismatches = None
    if same_count:
        displacement = sum(
            abs(left_part[0] - right_part[0]) + abs(left_part[1] - right_part[1])
            for left_part, right_part in zip(left_chain, right_chain, strict=True)
        )
        phase_mismatches = sum(
            left_part[2] != right_part[2]
            for left_part, right_part in zip(left_chain, right_chain, strict=True)
        )
    left_span = int(left["end"]) - int(left["start"]) + 1
    right_span = int(right["end"]) - int(right["start"]) + 1
    overlap = max(
        0,
        min(int(left["end"]), int(right["end"]))
        - max(int(left["start"]), int(right["start"]))
        + 1,
    )
    return {
        "left_candidate_digest": left["candidate_digest"],
        "right_candidate_digest": right["candidate_digest"],
        "same_cds_count": same_count,
        "cds_counts": [len(left_chain), len(right_chain)],
        "summed_boundary_displacement_bp": displacement,
        "phase_mismatch_segments": phase_mismatches,
        "span_overlap_fraction_of_shorter": round(
            overlap / min(left_span, right_span), 6
        ),
        "same_strand": left["strand"] == right["strand"],
    }


def _conflict_summaries(
    conflict_members: dict[str, list[str]],
    candidates_by_digest: dict[str, dict[str, Any]],
) -> list[dict[str, Any]]:
    summaries: list[dict[str, Any]] = []
    for conflict, member_digests in sorted(conflict_members.items()):
        members = [candidates_by_digest[digest] for digest in sorted(member_digests)]
        differences = [
            _chain_difference(members[left], members[right])
            for left in range(len(members))
            for right in range(left + 1, len(members))
        ]
        summaries.append(
            {
                "conflict_set_digest": conflict,
                "member_digests": [member["candidate_digest"] for member in members],
                "member_count": len(members),
                "seqids": sorted({str(member["seqid"]) for member in members}),
                "shared_axis": {
                    "start": min(int(member["start"]) for member in members),
                    "end": max(int(member["end"]) for member in members),
                },
                "pairwise_differences": differences,
            }
        )
    return summaries


def _triage(
    report_candidates: list[dict[str, Any]],
    conflict_members: dict[str, list[str]],
    *,
    limit: int = 20,
) -> dict[str, Any]:
    by_digest = {row["candidate_digest"]: row for row in report_candidates}
    unreviewed = [row for row in report_candidates if row["review"]["decision"] == "unreviewed"]
    deferred = [row for row in report_candidates if row["review"]["decision"] == "defer"]
    unresolved: list[dict[str, Any]] = []
    for conflict, members in sorted(conflict_members.items()):
        rows = [by_digest[digest] for digest in members]
        decisions = {row["review"]["decision"] for row in rows}
        if "accept_addition" in decisions or not decisions & {"unreviewed", "defer"}:
            continue
        lead = min(
            rows,
            key=lambda row: (
                row["score"] is None,
                -(row["score"] or 0.0),
                row["candidate_digest"],
            ),
        )
        unresolved.append(
            {
                "conflict_set_digest": conflict,
                "member_count": len(rows),
                "lead_candidate_digest": lead["candidate_digest"],
                "lead_gene_id": lead["gene_id"],
                "lead_rank_position": lead["rank_position"],
            }
        )
    unresolved.sort(
        key=lambda row: (
            row["lead_rank_position"] is None,
            row["lead_rank_position"] or 0,
            row["conflict_set_digest"],
        )
    )

    def slim(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
        return [
            {
                key: row[key]
                for key in (
                    "candidate_digest",
                    "gene_id",
                    "seqid",
                    "start",
                    "end",
                    "rank_position",
                    "rank_percentile",
                    "rank_band",
                )
            }
            for row in rows[:limit]
        ]

    return {
        "policy": "inspection_order_only_never_acceptance_advice",
        "counts": {
            "unreviewed": len(unreviewed),
            "deferred": len(deferred),
            "unresolved_conflict_sets": len(unresolved),
        },
        "top_unreviewed": slim(unreviewed),
        "deferred": slim(deferred),
        "unresolved_conflict_sets": unresolved[:limit],
    }


def _status_and_interpretation(
    *,
    total: int,
    reviews: dict[str, dict[str, str]],
    accepted_reviews: set[str],
    conflict_members: dict[str, list[str]],
    patch_digests: set[str] | None,
    byte_identical_reversion: bool | None,
    automatic_approval: int,
    automatic_approval_checked: bool,
    automatic_approval_fields: int,
    topology_available: int | None,
) -> tuple[str, str, list[dict[str, str]]]:
    reviewed = len(reviews)
    collisions = sum(
        1
        for members in conflict_members.values()
        if sum(digest in accepted_reviews for digest in members) > 1
    )
    patch_consistent = patch_digests is not None and patch_digests == accepted_reviews
    patch_mismatch = patch_digests is not None and not patch_consistent
    if automatic_approval or collisions or patch_mismatch:
        state = "attention_required"
        headline = "Safety checks require attention"
    elif reviewed < total:
        state = "review_in_progress"
        headline = f"Human review is {reviewed / total:.1%} complete"
    elif patch_consistent and byte_identical_reversion is True:
        state = "validated_reversible_run"
        headline = "Review complete and byte-identical reversion verified"
    elif patch_consistent:
        state = "ready_to_apply"
        headline = "Reviewed patch is ready for controlled application"
    else:
        state = "review_complete_patch_pending"
        headline = "Review is complete; compile the immutable patch next"
    notes: list[dict[str, str]] = []
    notes.append(
        {
            "level": "positive" if reviewed == total else "action",
            "title": "Review coverage",
            "text": (
                f"All {total:,} candidates have explicit decisions."
                if reviewed == total
                else f"{total - reviewed:,} candidates remain unreviewed; missing rows are never treated as acceptance."
            ),
        }
    )
    notes.append(
        {
            "level": "positive" if byte_identical_reversion is True else "action",
            "title": "Patch verification",
            "text": (
                "The applied annotation was reverted to the exact source bytes."
                if byte_identical_reversion is True
                else (
                    "A compiled patch matches the accepted candidates; apply it to a new GFF3 and verify byte-identical reversion."
                    if patch_consistent
                    else (
                        "The supplied patch does not match the accepted review set and must not be applied."
                        if patch_mismatch
                        else "Compile a patch from the final ledger before application and reversion testing."
                    )
                )
            ),
        }
    )
    notes.append(
        {
            "level": "positive" if collisions == 0 else "critical",
            "title": "Conflict safety",
            "text": (
                f"All {_plural(len(conflict_members), 'mutually exclusive conflict set')} obey the at-most-one-acceptance rule."
                if collisions == 0
                else f"{collisions:,} conflict sets contain multiple accepted alternatives and must be resolved."
            ),
        }
    )
    if topology_available is None:
        notes.append(
            {
                "level": "neutral",
                "title": "Evolutionary context",
                "text": "Topology evidence was not supplied. The core review and patch workflow remains valid without it.",
            }
        )
    else:
        notes.append(
            {
                "level": "neutral",
                "title": "Evolutionary context",
                "text": f"Topology is available for {topology_available:,}/{total:,} candidates ({topology_available / total:.1%}) and is interpreted only as applicability-gated supporting evidence.",
            }
        )
    if not automatic_approval_checked:
        notes.append(
            {
                "level": "neutral",
                "title": "Decision boundary",
                "text": (
                    "No score table containing an automatic-approval field was supplied, so this report cannot verify that run-level state. "
                    "PloidyPatch policy still forbids automatic approval; supply the checksum-bound score table to audit its execution."
                ),
            }
        )
    else:
        notes.append(
            {
                "level": "positive" if automatic_approval == 0 else "critical",
                "title": "Decision boundary",
                "text": (
                    f"Checked {total:,} candidates across {_plural(automatic_approval_fields, 'automatic-approval field')}: none was automatically approved. Ranking only prioritizes human review."
                    if automatic_approval == 0
                    else f"{automatic_approval:,} automatic-approval flags were observed and contradict the review-only policy."
                ),
            }
        )
    return state, headline, notes


def _render_html(
    report: dict[str, Any], *, embedded_digests: set[str]
) -> str:
    triage = {
        **report["triage"],
        "top_unreviewed": [
            row
            for row in report["triage"]["top_unreviewed"]
            if row["candidate_digest"] in embedded_digests
        ],
        "deferred": [
            row
            for row in report["triage"]["deferred"]
            if row["candidate_digest"] in embedded_digests
        ],
        "unresolved_conflict_sets": [
            row
            for row in report["triage"]["unresolved_conflict_sets"]
            if row["lead_candidate_digest"] in embedded_digests
        ],
    }
    html_report = {
        **report,
        "triage": triage,
        "candidates": [
            row
            for row in report["candidates"]
            if row["candidate_digest"] in embedded_digests
        ],
        "conflict_sets": [
            conflict
            for conflict in report["conflict_sets"]
            if set(conflict["member_digests"]) <= embedded_digests
        ],
    }
    payload = json.dumps(html_report, ensure_ascii=False, separators=(",", ":"))
    payload = payload.replace("<", "\\u003c").replace(">", "\\u003e").replace("&", "\\u0026")
    title = html.escape(str(report["run"]["title"]))
    asset_root = resources.files("ploidypatch").joinpath("report_assets")
    template = asset_root.joinpath("report.html").read_text(encoding="utf-8")
    report_css = asset_root.joinpath("report.css").read_text(encoding="utf-8")
    report_js = asset_root.joinpath("report.js").read_text(encoding="utf-8")
    replacements = {
        "__REPORT_TITLE__": title,
        "__REPORT_CSS__": report_css,
        "__REPORT_JS__": report_js,
        "__PLOIDYPATCH_VERSION__": html.escape(__version__),
        "__REPORT_PAYLOAD__": payload,
    }
    token = re.compile("|".join(re.escape(value) for value in replacements))
    return token.sub(lambda match: replacements[match.group(0)], template)


def build_review_report(
    *,
    candidate_gff_path: str | Path,
    pool_decisions_tsv_path: str | Path,
    pool_manifest_json_path: str | Path,
    output_dir: str | Path,
    review_decisions_tsv_path: str | Path | None = None,
    scores_tsv_path: str | Path | None = None,
    copy_features_tsv_path: str | Path | None = None,
    topology_features_tsv_path: str | Path | None = None,
    patch_edits_json_path: str | Path | None = None,
    run_summary_json_path: str | Path | None = None,
    title: str = "Annotation repair review",
    max_embedded_candidates: int = 5000,
) -> dict[str, Any]:
    title = title.strip()
    if not title or len(title) > 160 or any(character in title for character in "\r\n"):
        raise ValueError("Report title must contain 1-160 characters on one line")
    if max_embedded_candidates < 1:
        raise ValueError("max_embedded_candidates must be positive")
    candidate_path = Path(candidate_gff_path)
    pool_path = Path(pool_decisions_tsv_path)
    manifest_path = Path(pool_manifest_json_path)
    output = Path(output_dir)
    working = Path(f"{output}.working")
    if output.exists() or output.is_symlink() or working.exists() or working.is_symlink():
        raise FileExistsError(f"Refusing to overwrite report output: {output}")
    candidates = _candidate_models(candidate_path)
    manifest = _read_json(manifest_path, "candidate-pool manifest")
    pool_fields, pool_rows = _read_tsv(
        pool_path,
        label="candidate-pool decisions",
        required={"consensus_digest", "status", "conflict_set_digest", "conflict_member_count"},
        reject_truth=True,
    )
    pool = _index_rows(pool_rows, "consensus_digest", "candidate-pool decisions")
    accepted_pool = {digest: row for digest, row in pool.items() if row["status"] == "accepted"}
    if (
        manifest.get("schema_version") != POOL_SCHEMA_VERSION
        or manifest.get("outputs", {}).get("candidate_gff", {}).get("sha256") != sha256_file(candidate_path)
        or manifest.get("outputs", {}).get("decisions", {}).get("sha256") != sha256_file(pool_path)
        or set(candidates) != set(accepted_pool)
    ):
        raise ValueError("Candidate pool, manifest, decisions, and GFF3 fail exact-universe validation")
    if int(manifest.get("counts", {}).get("accepted_models", -1)) != len(candidates):
        raise ValueError("Candidate-pool manifest accepted count differs")
    if int(manifest.get("outputs", {}).get("decisions", {}).get("rows", -1)) != len(pool):
        raise ValueError("Candidate-pool manifest decision row count differs")
    parameters = manifest.get("parameters", {})
    if (
        parameters.get("redundancy_policy") != "retain_distinct_chains"
        or parameters.get("conflict_action") != "retain_all_for_ranking_and_review"
        or parameters.get("min_method_support") != 1
        or parameters.get("max_redundancy_overlap") != 0.5
    ):
        raise ValueError("Candidate pool is not the review-only chain-preserving policy")
    conflict_members = _validate_conflicts(accepted_pool=accepted_pool, manifest=manifest)

    reviews: dict[str, dict[str, str]] = {}
    review_path = Path(review_decisions_tsv_path) if review_decisions_tsv_path else None
    if review_path is not None:
        _, review_rows = _read_tsv(
            review_path,
            label="review decisions",
            required=set(REVIEW_FIELDS),
            exact=REVIEW_FIELDS,
        )
        reviews = _index_rows(review_rows, "candidate_digest", "review decisions")
        if not set(reviews) <= set(candidates):
            raise ValueError("Review decisions contain candidates outside the pool")
        _validate_reviews(reviews)

    score_path = Path(scores_tsv_path) if scores_tsv_path else None
    score_fields, scores = _optional_candidate_table(score_path, label="candidate scores")
    feature_path = Path(copy_features_tsv_path) if copy_features_tsv_path else None
    feature_fields, features = _optional_candidate_table(feature_path, label="copy features")
    topology_path = Path(topology_features_tsv_path) if topology_features_tsv_path else None
    topology_fields, topologies = _optional_candidate_table(topology_path, label="topology features")
    for label, table in (("candidate scores", scores), ("copy features", features), ("topology features", topologies)):
        if table and set(table) != set(candidates):
            raise ValueError(f"{label} candidate universe differs from the pool")

    patch_path = Path(patch_edits_json_path) if patch_edits_json_path else None
    patch_digests: set[str] | None = None
    if patch_path is not None:
        patch = _read_json(patch_path, "patch edits")
        events = patch.get("events")
        if not isinstance(events, list):
            raise ValueError("Patch edits lack an events list")
        raw_patch_digests = [str(event.get("consensus_digest", "")) for event in events]
        patch_digests = set(raw_patch_digests)
        if len(patch_digests) != len(raw_patch_digests):
            raise ValueError("Patch edits contain duplicate candidate digests")
        if "" in patch_digests or not patch_digests <= set(candidates):
            raise ValueError("Patch edits contain unknown or empty candidate digests")

    summary_path = Path(run_summary_json_path) if run_summary_json_path else None
    run_summary = _read_json(summary_path, "run summary") if summary_path else None
    if run_summary is not None and run_summary.get("byte_identical_reversion") is not True:
        raise ValueError("Run summary does not prove byte-identical reversion")
    accepted_reviews = {
        digest for digest, row in reviews.items() if row["review_decision"] == "accept_addition"
    }
    if run_summary is not None:
        if run_summary.get("automatic_approval") is not False:
            raise ValueError("Run summary does not prove automatic approval was disabled")
        if run_summary.get("source_sha256") != run_summary.get("reverted_sha256"):
            raise ValueError("Run summary source and reverted SHA256 differ")
        if int(run_summary.get("accepted_additions", -1)) != len(accepted_reviews):
            raise ValueError("Run summary accepted count differs from the review ledger")
        selected = run_summary.get("selected_candidate_digests")
        if not isinstance(selected, list) or set(map(str, selected)) != accepted_reviews:
            raise ValueError("Run summary selected candidates differ from the review ledger")
        if len(selected) != len(set(map(str, selected))):
            raise ValueError("Run summary selected candidates contain duplicates")
        if patch_digests is None or patch_digests != accepted_reviews:
            raise ValueError("A validated run summary requires a matching patch artifact")
    automatic_fields = [field for field in score_fields if field.endswith("automatic_approval")]
    automatic_approval_checked = bool(scores) and bool(automatic_fields)
    automatic_approval = 0
    if automatic_approval_checked:
        for digest, row in scores.items():
            values: list[float] = []
            for field in automatic_fields:
                try:
                    value = float(row.get(field, "0") or "0")
                except ValueError as exc:
                    raise ValueError(
                        f"Candidate automatic-approval flag is not numeric: {digest}"
                    ) from exc
                if not math.isfinite(value) or value not in {0.0, 1.0}:
                    raise ValueError(
                        f"Candidate automatic-approval flag must be zero or one: {digest}"
                    )
                values.append(value)
            automatic_approval += int(any(value == 1.0 for value in values))
    topology_available: int | None = None
    if topologies:
        topology_available = sum(
            int(float(row.get("topology_available", "0") or "0"))
            for row in topologies.values()
        )

    support_combinations: Counter[str] = Counter()
    report_candidates: list[dict[str, Any]] = []
    for digest in sorted(candidates):
        base = candidates[digest]
        pool_row = accepted_pool[digest]
        feature = features.get(digest, {})
        methods = base["support_methods"] or [
            method for method in METHODS if feature.get(f"has_{method}") == "1"
        ]
        if methods:
            support_combinations["+".join(sorted(methods))] += 1
        score_field, score = _float(scores.get(digest, {}), SCORE_FIELDS)
        percentile_field, percentile = _float(scores.get(digest, {}), PERCENTILE_FIELDS)
        review = reviews.get(digest)
        topology = topologies.get(digest)
        report_candidates.append(
            {
                **base,
                "support_methods": sorted(methods),
                "conflict_set_digest": pool_row["conflict_set_digest"],
                "conflict_member_count": int(pool_row["conflict_member_count"]),
                "score": score,
                "score_field": score_field,
                "rank_percentile": percentile,
                "rank_percentile_field": percentile_field,
                "review": {
                    "decision": review["review_decision"] if review else "unreviewed",
                    "reviewer": review["reviewer"] if review else "",
                    "reviewed_at_utc": review["reviewed_at_utc"] if review else "",
                    "rationale": review["rationale"] if review else "",
                },
                "topology": {
                    "available": (
                        None
                        if topology is None
                        else bool(int(float(topology.get("topology_available", "0") or "0")))
                    ),
                    "reason": "" if topology is None else topology.get("topology_reason", ""),
                    "partner_gene_id": (
                        ""
                        if topology is None
                        or topology.get("partner_gene_id", "") in {"", ".", "NA"}
                        else topology["partner_gene_id"]
                    ),
                    "coherence_score": (
                        None
                        if topology is None
                        or topology.get("topology_coherence_score", "") in {"", ".", "NA"}
                        else float(topology["topology_coherence_score"])
                    ),
                },
                "patch_included": patch_digests is not None and digest in patch_digests,
            }
        )
    report_candidates.sort(
        key=lambda row: (
            row["score"] is None,
            -(row["score"] or 0.0),
            row["candidate_digest"],
        )
    )
    scored_candidates = [row for row in report_candidates if row["score"] is not None]
    scored_total = len(scored_candidates)
    for position, candidate in enumerate(scored_candidates, start=1):
        candidate["rank_position"] = position
        if candidate["rank_percentile"] is None:
            candidate["rank_percentile"] = (
                1.0 if scored_total == 1 else 1.0 - (position - 1) / (scored_total - 1)
            )
            candidate["rank_percentile_field"] = "derived_from_descending_review_rank"
        percentile = float(candidate["rank_percentile"])
        if not 0.0 <= percentile <= 1.0:
            raise ValueError("Candidate rank percentile must be between zero and one")
        candidate["rank_band"] = (
            "top_1_percent"
            if percentile >= 0.99
            else "top_5_percent"
            if percentile >= 0.95
            else "top_20_percent"
            if percentile >= 0.80
            else "remainder"
        )
    for row in report_candidates:
        if row["score"] is None:
            row["rank_position"] = None
            row["rank_band"] = None
    selected_score_fields = {row["score_field"] for row in report_candidates if row["score_field"]}
    if len(selected_score_fields) > 1:
        raise ValueError("Candidate score rows do not use one consistent report score field")
    selected_score_field = next(iter(selected_score_fields), "")
    score_interpretation = (
        "calibrated_model_output_not_automatic_decision"
        if selected_score_field == "model_calibrated_probability"
        else "review_rank_not_probability"
    )
    score_values = [
        float(row["score"])
        for row in report_candidates
        if row["score"] is not None
    ]
    by_seqid = Counter(str(row["seqid"]) for row in report_candidates)

    state, headline, interpretation = _status_and_interpretation(
        total=len(candidates),
        reviews=reviews,
        accepted_reviews=accepted_reviews,
        conflict_members=conflict_members,
        patch_digests=patch_digests,
        byte_identical_reversion=(
            None
            if run_summary is None
            else bool(run_summary["byte_identical_reversion"])
        ),
        automatic_approval=automatic_approval,
        automatic_approval_checked=automatic_approval_checked,
        automatic_approval_fields=len(automatic_fields),
        topology_available=topology_available,
    )
    interpretation.insert(
        1,
        {
            "level": "neutral",
            "title": "Review priority",
            "text": (
                f"{selected_score_field} spans {min(score_values):.4g} to {max(score_values):.4g}. It controls inspection order, not acceptance."
                if score_values
                else "No score table was supplied; candidates use deterministic digest order and still require explicit review."
            ),
        },
    )
    decision_counts = Counter(
        row["review"]["decision"] for row in report_candidates
    )
    for decision in (*sorted(REVIEW_DECISIONS), "unreviewed"):
        decision_counts.setdefault(decision, 0)
    inputs = [
        _input_record(candidate_path, "candidate_gff"),
        _input_record(pool_path, "pool_decisions"),
        _input_record(manifest_path, "pool_manifest"),
    ]
    for path, role in (
        (review_path, "review_decisions"),
        (score_path, "scores"),
        (feature_path, "copy_features"),
        (topology_path, "topology_features"),
        (patch_path, "patch_edits"),
        (summary_path, "run_summary"),
    ):
        if path is not None:
            inputs.append(_input_record(path, role))
    generated = _generated_at_utc()
    state_detail = {
        "review_in_progress": f"Review {len(candidates) - len(reviews):,} remaining candidates; unreviewed rows are excluded.",
        "review_complete_patch_pending": "Compile the immutable patch from the completed human ledger.",
        "ready_to_apply": "Apply the matching patch to a new output and run the reversion check.",
        "validated_reversible_run": "Archive this checksum-bound report with the patch and run summary.",
        "attention_required": "Resolve every safety finding before applying or distributing a patch.",
    }[state]
    report: dict[str, Any] = {
        "schema_version": REPORT_SCHEMA_VERSION,
        "generator": {"name": "PloidyPatch", "version": __version__},
        "truth_access": False,
        "run": {"title": title, "generated_at_utc": generated},
        "summary": {
            "state": state,
            "headline": headline,
            "detail": state_detail,
        },
        "counts": {
            "candidates": len(candidates),
            "conflict_sets": len(conflict_members),
            "conflicted_candidates": sum(len(value) for value in conflict_members.values()),
            "nonconflicted": len(candidates) - sum(len(value) for value in conflict_members.values()),
            "reviewed": len(reviews),
            "accepted_review": len(accepted_reviews),
            "patch_events": 0 if patch_digests is None else len(patch_digests),
            "byte_identical_reversion": (
                None
                if run_summary is None
                else bool(run_summary["byte_identical_reversion"])
            ),
            "automatic_approval": automatic_approval,
            "automatic_approval_checked": automatic_approval_checked,
            "topology_available": topology_available,
            "review_decisions": dict(sorted(decision_counts.items())),
        },
        "method_overlap": {
            "available": bool(support_combinations),
            "by_method": {
                method: sum(
                    count
                    for combination, count in support_combinations.items()
                    if method in combination.split("+")
                )
                for method in METHODS
            },
            "all_three": support_combinations.get("gemoma+lifton+miniprot", 0),
            "combinations": dict(sorted(support_combinations.items())),
        },
        "distributions": {
            "by_seqid": [
                {"seqid": seqid, "count": count}
                for seqid, count in sorted(
                    by_seqid.items(), key=lambda item: (-item[1], item[0])
                )
            ],
            "score_histogram": _score_histogram(score_values),
        },
        "interpretation": interpretation,
        "triage": _triage(report_candidates, conflict_members),
        "conflict_sets": _conflict_summaries(
            conflict_members,
            {row["candidate_digest"]: row for row in report_candidates},
        ),
        "policy": {
            "score_interpretation": score_interpretation,
            "displayed_score_field": selected_score_field or None,
            "score_direction": SCORE_DIRECTION,
            "automatic_approval": False,
            "automatic_approval_verified": automatic_approval_checked,
            "unreviewed_policy": "exclude",
            "topology_role": "optional_applicability_gated_support",
            "conflict_rule": "at_most_one_accepted_alternative_per_conflict_set",
        },
        "inputs": inputs,
        "candidates": report_candidates,
    }

    initially_embedded = report_candidates[:max_embedded_candidates]
    embedded_digests = {row["candidate_digest"] for row in initially_embedded}
    conflict_closure_extra = 0
    for members in conflict_members.values():
        if embedded_digests.intersection(members):
            for digest in members:
                if digest not in embedded_digests:
                    embedded_digests.add(digest)
                    conflict_closure_extra += 1
    report["embedding"] = {
        "maximum_requested": max_embedded_candidates,
        "embedded": len(embedded_digests),
        "total": len(report_candidates),
        "truncated": len(embedded_digests) < len(report_candidates),
        "conflict_closure_extra": conflict_closure_extra,
        "order": "descending_review_priority_with_conflict_set_closure",
        "full_universe_files": ["report.json", "candidates.tsv", "review_ledger_template.tsv"],
    }

    output.parent.mkdir(parents=True, exist_ok=True)
    working.mkdir()
    try:
        (working / "report.json").write_text(
            json.dumps(report, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
            newline="\n",
        )
        (working / "index.html").write_text(
            _render_html(report, embedded_digests=embedded_digests),
            encoding="utf-8",
            newline="\n",
        )
        with (working / "candidates.tsv").open("x", encoding="utf-8", newline="") as handle:
            fields = (
                "candidate_digest",
                "gene_id",
                "seqid",
                "start",
                "end",
                "strand",
                "support_methods",
                "conflict_set_digest",
                "conflict_member_count",
                "score",
                "score_field",
                "rank_percentile",
                "rank_position",
                "rank_band",
                "review_decision",
                "reviewer",
                "rationale",
                "topology_available",
                "topology_reason",
                "patch_included",
            )
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for candidate in report_candidates:
                writer.writerow(
                    {
                        "candidate_digest": candidate["candidate_digest"],
                        "gene_id": candidate["gene_id"],
                        "seqid": candidate["seqid"],
                        "start": candidate["start"],
                        "end": candidate["end"],
                        "strand": candidate["strand"],
                        "support_methods": ",".join(candidate["support_methods"]),
                        "conflict_set_digest": candidate["conflict_set_digest"],
                        "conflict_member_count": candidate["conflict_member_count"],
                        "score": "" if candidate["score"] is None else repr(candidate["score"]),
                        "score_field": candidate["score_field"],
                        "rank_percentile": "" if candidate["rank_percentile"] is None else repr(candidate["rank_percentile"]),
                        "rank_position": "" if candidate["rank_position"] is None else candidate["rank_position"],
                        "rank_band": candidate["rank_band"] or "",
                        "review_decision": candidate["review"]["decision"],
                        "reviewer": candidate["review"]["reviewer"],
                        "rationale": candidate["review"]["rationale"],
                        "topology_available": "" if candidate["topology"]["available"] is None else int(candidate["topology"]["available"]),
                        "topology_reason": candidate["topology"]["reason"],
                        "patch_included": int(candidate["patch_included"]),
                    }
                )
        with (working / "review_ledger_template.tsv").open(
            "x", encoding="utf-8", newline=""
        ) as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=REVIEW_FIELDS,
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            for row in report_candidates:
                existing = reviews.get(row["candidate_digest"])
                writer.writerow(
                    {
                        "candidate_digest": row["candidate_digest"],
                        "review_decision": "" if existing is None else existing["review_decision"],
                        "reviewer": "" if existing is None else existing["reviewer"],
                        "reviewed_at_utc": "" if existing is None else existing["reviewed_at_utc"],
                        "rationale": "" if existing is None else existing["rationale"],
                    }
                )
        write_sha256sums(working)
        verify_sha256sums(working, ignore_checksum_file=True)
        os.replace(working, output)
    except BaseException:
        shutil.rmtree(working, ignore_errors=True)
        raise
    return report
