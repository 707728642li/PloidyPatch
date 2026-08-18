from __future__ import annotations

import csv
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path
from typing import Any, Sequence

from . import __version__
from .baseline import _file_sha256
from .homeolog_ranker import _load_json, _percentile_ranks
from .support_ranker import SUPPORT_CONDITIONED_SCORE_SCHEMA_VERSION


CONFLICT_GUARD_SCORE_SCHEMA_VERSION = (
    "ploidypatch.conflict_winner_guard_scores.v1"
)
BASELINE_FIELD = "v03_baseline_logit"
PRIMARY_FIELD = "v03_primary_rank_score"
POOL_SCHEMA_VERSION = "ploidypatch.method_candidate_pool.v2"
REQUIRED_POOL_PARAMETERS: dict[str, Any] = {
    "conflict_action": "retain_all_for_ranking_and_review",
    "max_redundancy_overlap": 0.5,
    "min_method_support": 1,
    "redundancy_policy": "retain_distinct_chains",
}


def _read_indexed_tsv(
    path: Path, key_fields: Sequence[str]
) -> tuple[list[str], list[dict[str, str]], dict[str, dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Missing TSV header: {path}")
        fields = list(reader.fieldnames)
        key = next((field for field in key_fields if field in fields), None)
        if key is None:
            raise ValueError(f"Missing candidate key {key_fields}: {path}")
        rows = list(reader)
    indexed = {row[key]: row for row in rows}
    if len(indexed) != len(rows) or "" in indexed:
        raise ValueError(f"Empty or duplicate candidate key: {path}")
    return fields, rows, indexed


def _winner(
    digests: Sequence[str], scores: dict[str, float]
) -> str:
    return min(digests, key=lambda digest: (-scores[digest], digest))


def _validate_blind_fields(fields: Sequence[str], path: Path) -> None:
    forbidden = [
        field
        for field in fields
        if "label" in field.lower()
        or "truth" in field.lower()
        or field.lower().startswith("v04_")
    ]
    if forbidden:
        raise ValueError(
            f"Truth/label or reserved v0.4 field in blind input {path}: {forbidden}"
        )


def _winner_mapping_sha256(mapping: dict[str, str]) -> str:
    payload = "".join(
        f"{conflict}\t{mapping[conflict]}\n" for conflict in sorted(mapping)
    )
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def compute_conflict_winner_guard(
    *,
    digests: Sequence[str],
    conflict_by_digest: dict[str, str],
    baseline_scores: dict[str, float],
    primary_scores: dict[str, float],
) -> dict[str, Any]:
    """Compute the deterministic label-free guard using one production primitive."""

    ordered = list(digests)
    universe = set(ordered)
    if not ordered or len(universe) != len(ordered) or "" in universe:
        raise ValueError("Guard digests must be nonempty and unique")
    for name, mapping in (
        ("conflict", conflict_by_digest),
        ("baseline", baseline_scores),
        ("primary", primary_scores),
    ):
        if set(mapping) != universe:
            raise ValueError(f"Guard {name} universe differs from candidate digests")
    if any(
        not math.isfinite(float(scores[digest]))
        for scores in (baseline_scores, primary_scores)
        for digest in ordered
    ):
        raise ValueError("Guard scores must be finite")

    conflicts: dict[str, list[str]] = defaultdict(list)
    for digest in ordered:
        conflict = conflict_by_digest[digest]
        if conflict:
            conflicts[conflict].append(digest)
    if any(len(members) < 2 for members in conflicts.values()):
        raise ValueError("Conflict annotations contain singleton set(s)")

    baseline_winners = {
        conflict: _winner(members, baseline_scores)
        for conflict, members in conflicts.items()
    }
    primary_winners = {
        conflict: _winner(members, primary_scores)
        for conflict, members in conflicts.items()
    }
    guarded_sets = {
        conflict
        for conflict in conflicts
        if baseline_winners[conflict] != primary_winners[conflict]
    }
    guarded_digests = {
        digest for conflict in guarded_sets for digest in conflicts[conflict]
    }
    guard_scores = {
        digest: (
            float(baseline_scores[digest])
            if digest in guarded_digests
            else float(primary_scores[digest])
        )
        for digest in ordered
    }
    guard_winners = {
        conflict: _winner(members, guard_scores)
        for conflict, members in conflicts.items()
    }
    mismatches = [
        conflict
        for conflict in conflicts
        if baseline_winners[conflict] != guard_winners[conflict]
    ]
    if mismatches:
        raise AssertionError("Conflict guard failed to preserve baseline winner")
    mapping_text = "".join(
        f"{conflict}\t{baseline_winners[conflict]}\t"
        f"{primary_winners[conflict]}\t{guard_winners[conflict]}\n"
        for conflict in sorted(conflicts)
    )
    return {
        "scores": guard_scores,
        "conflicts": conflicts,
        "guarded_sets": guarded_sets,
        "guarded_digests": guarded_digests,
        "baseline_winners": baseline_winners,
        "primary_winners": primary_winners,
        "guard_winners": guard_winners,
        "baseline_winner_mapping_sha256": _winner_mapping_sha256(
            baseline_winners
        ),
        "primary_winner_mapping_sha256": _winner_mapping_sha256(primary_winners),
        "guard_winner_mapping_sha256": _winner_mapping_sha256(guard_winners),
        "winner_mapping_sha256": hashlib.sha256(
            mapping_text.encode("utf-8")
        ).hexdigest(),
        "winner_mismatch_count": 0,
    }


def apply_conflict_winner_guard(
    *,
    v03_score_tsv_path: str | Path,
    pool_decisions_tsv_path: str | Path,
    pool_manifest_json_path: str | Path,
    output_tsv_path: str | Path,
) -> dict[str, Any]:
    """Apply a label-free fallback when v0.3 changes a conflict-set winner.

    The fallback is deliberately local: candidates outside winner-disagreement
    conflict sets retain their v0.3 primary score. All members of a guarded
    set receive their frozen baseline logits, which makes its deterministic
    top-1 candidate exactly identical to the baseline winner.
    """

    score_path = Path(v03_score_tsv_path)
    decisions_path = Path(pool_decisions_tsv_path)
    pool_manifest_path = Path(pool_manifest_json_path)
    output_path = Path(output_tsv_path)
    score_manifest_path = Path(str(score_path) + ".manifest.json")
    output_manifest_path = Path(str(output_path) + ".manifest.json")
    for required in (
        score_path,
        decisions_path,
        score_manifest_path,
        pool_manifest_path,
    ):
        if not required.is_file() or required.stat().st_size == 0:
            raise ValueError(f"Missing or empty conflict-guard input: {required}")
    collisions = [
        path for path in (output_path, output_manifest_path) if path.exists()
    ]
    if collisions:
        raise FileExistsError(
            "Refusing to overwrite conflict-guard artifact(s): "
            + ", ".join(str(path) for path in collisions)
        )

    score_manifest = _load_json(score_manifest_path)
    if (
        score_manifest.get("schema_version")
        != SUPPORT_CONDITIONED_SCORE_SCHEMA_VERSION
        or score_manifest.get("truth_access") is not False
        or score_manifest.get("outputs", {}).get("scores", {}).get("sha256")
        != _file_sha256(score_path)
    ):
        raise ValueError("v0.3 score manifest fails schema, truth, or checksum gate")

    score_fields, score_rows, scores_by_digest = _read_indexed_tsv(
        score_path, ("candidate_digest",)
    )
    _validate_blind_fields(score_fields, score_path)
    required_score_fields = {"candidate_digest", BASELINE_FIELD, PRIMARY_FIELD}
    if not required_score_fields <= set(score_fields):
        raise ValueError("v0.3 scores lack conflict-guard fields")
    decision_fields, decision_rows, decisions_by_digest = _read_indexed_tsv(
        decisions_path, ("candidate_digest", "consensus_digest")
    )
    _validate_blind_fields(decision_fields, decisions_path)
    required_decision_fields = {
        "conflict_set_digest",
        "conflict_member_count",
        "status",
    }
    if not required_decision_fields <= set(decision_fields):
        raise ValueError("Pool decisions lack conflict provenance fields")
    pool_manifest = _load_json(pool_manifest_path)
    if (
        pool_manifest.get("schema_version") != POOL_SCHEMA_VERSION
        or pool_manifest.get("outputs", {}).get("decisions", {}).get("sha256")
        != _file_sha256(decisions_path)
    ):
        raise ValueError("Pool manifest fails schema or decisions checksum gate")
    parameters = pool_manifest.get("parameters", {})
    if any(parameters.get(key) != value for key, value in REQUIRED_POOL_PARAMETERS.items()):
        raise ValueError("Pool manifest is not the frozen chain-preserving policy")
    declared_rows = pool_manifest.get("outputs", {}).get("decisions", {}).get("rows")
    if int(declared_rows) != len(decision_rows):
        raise ValueError("Pool manifest decisions row count differs")
    accepted = {
        digest: row
        for digest, row in decisions_by_digest.items()
        if row.get("status") == "accepted"
    }
    if set(scores_by_digest) != set(accepted):
        raise ValueError("v0.3 scores and accepted pool decisions differ")
    declared_accepted = pool_manifest.get("counts", {}).get("accepted_models")
    if int(declared_accepted) != len(accepted):
        raise ValueError("Pool manifest accepted-model count differs")

    baseline: dict[str, float] = {}
    primary: dict[str, float] = {}
    for digest, row in scores_by_digest.items():
        for field, target in ((BASELINE_FIELD, baseline), (PRIMARY_FIELD, primary)):
            value = float(row[field])
            if not math.isfinite(value):
                raise ValueError(f"Non-finite {field} for {digest}")
            target[digest] = value

    conflict_by_digest = {
        digest: row.get("conflict_set_digest", "") for digest, row in accepted.items()
    }
    computation = compute_conflict_winner_guard(
        digests=[row["candidate_digest"] for row in score_rows],
        conflict_by_digest=conflict_by_digest,
        baseline_scores=baseline,
        primary_scores=primary,
    )
    conflicts: dict[str, list[str]] = computation["conflicts"]
    guarded_sets: set[str] = computation["guarded_sets"]
    guarded_digests: set[str] = computation["guarded_digests"]
    actual_conflicted = sum(len(members) for members in conflicts.values())
    if int(pool_manifest.get("counts", {}).get("conflict_sets")) != len(conflicts):
        raise ValueError("Pool manifest conflict-set count differs")
    if int(pool_manifest.get("counts", {}).get("conflicted_chains")) != actual_conflicted:
        raise ValueError("Pool manifest conflicted-chain count differs")
    for digest, row in accepted.items():
        conflict = conflict_by_digest[digest]
        declared_size = int(row["conflict_member_count"])
        expected_size = len(conflicts[conflict]) if conflict else 1
        if declared_size != expected_size:
            raise ValueError(f"Conflict member count differs for {digest}")

    output_rows: list[dict[str, Any]] = []
    v04_scores: list[float] = []
    for row in score_rows:
        digest = row["candidate_digest"]
        conflict = accepted[digest].get("conflict_set_digest", "")
        guarded = digest in guarded_digests
        score = computation["scores"][digest]
        if guarded:
            source = "v03_baseline_conflict_winner_fallback"
        elif conflict:
            source = "v03_primary_same_conflict_winner"
        else:
            source = "v03_primary_no_conflict"
        v04_scores.append(score)
        output_rows.append(
            {
                **row,
                "v04_primary_rank_score": format(score, ".17g"),
                "v04_conflict_guard_applied": int(guarded),
                "v04_topology_abstained": int(guarded),
                "v04_score_source": source,
                "v04_automatic_approval": 0,
            }
        )

    percentiles = _percentile_ranks(v04_scores)
    output_fields = [
        *score_fields,
        "v04_primary_rank_score",
        "v04_primary_rank_percentile",
        "v04_conflict_guard_applied",
        "v04_topology_abstained",
        "v04_score_source",
        "v04_automatic_approval",
    ]
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=output_fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for row, percentile in zip(output_rows, percentiles, strict=True):
            row["v04_primary_rank_percentile"] = format(percentile, ".17g")
            writer.writerow(row)

    manifest: dict[str, Any] = {
        "schema_version": CONFLICT_GUARD_SCORE_SCHEMA_VERSION,
        "generator": {"name": "PloidyPatch", "version": __version__},
        "truth_access": False,
        "inputs": {
            "v03_scores": _file_sha256(score_path),
            "v03_score_manifest": _file_sha256(score_manifest_path),
            "pool_decisions": _file_sha256(decisions_path),
            "pool_manifest": _file_sha256(pool_manifest_path),
        },
        "policy": {
            "name": "baseline_fallback_on_conflict_winner_disagreement",
            "winner_tie_break": "descending_score_then_candidate_digest",
            "outside_guarded_sets": "retain_v03_primary_score",
            "inside_guarded_sets": "use_v03_baseline_logit_and_abstain_topology",
            "guarantee": "every_conflict_top1_identical_to_v03_baseline",
            "interpretation": "uncalibrated_review_rank_not_probability",
            "automatic_approval": False,
        },
        "counts": {
            "candidates": len(score_rows),
            "conflict_sets": len(conflicts),
            "winner_disagreement_sets": len(guarded_sets),
            "guarded_candidates": len(guarded_digests),
            "topology_abstained_candidates": len(guarded_digests),
            "automatic_approved": 0,
            "winner_mismatch_count": computation["winner_mismatch_count"],
        },
        "winner_audit": {
            "mapping_columns": [
                "conflict_set_digest",
                "baseline_winner",
                "v03_primary_winner",
                "v04_guard_winner",
            ],
            "mapping_sha256": computation["winner_mapping_sha256"],
            "baseline_mapping_sha256": computation[
                "baseline_winner_mapping_sha256"
            ],
            "v03_primary_mapping_sha256": computation[
                "primary_winner_mapping_sha256"
            ],
            "v04_guard_mapping_sha256": computation[
                "guard_winner_mapping_sha256"
            ],
            "mismatch_count": computation["winner_mismatch_count"],
        },
        "outputs": {
            "scores": {
                "file_name": output_path.name,
                "sha256": _file_sha256(output_path),
                "rows": len(output_rows),
            }
        },
    }
    with output_manifest_path.open("x", encoding="utf-8", newline="") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return manifest
