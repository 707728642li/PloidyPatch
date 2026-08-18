"""Species-neutral evaluation for a no-ranker, paired-event core H1."""
from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any, Mapping

from .artifact_manifest import sha256_file
from .bootstrap import paired_event_bootstrap


SCOPES = ("combined", "bua_only", "mauritiana_only")
ARMS = ("retain_distinct", "suppress_overlap")
EXPECTED_SCORE_KEYS = frozenset(f"{scope}_{arm}" for scope in SCOPES for arm in ARMS)


def load_json_object(path: str | Path) -> dict[str, Any]:
    source = Path(path)
    if not source.is_file() or source.is_symlink() or source.stat().st_size == 0:
        raise ValueError(f"Missing, empty or symlinked JSON: {source}")
    value = json.loads(source.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"Expected JSON object: {source}")
    return value


def score_collateral_gate(score: Mapping[str, Any]) -> bool:
    return (
        score.get("quality_gate", {}).get("grade") == "pass"
        and score.get("collateral_changes", {}).get(
            "baseline_transcript_structures_missing_from_candidate"
        )
        == 0
    )


def _arm_summary(score: Mapping[str, Any]) -> dict[str, Any]:
    recovery = score.get("event_recovery", {})
    collateral = score.get("collateral_changes", {})
    return {
        "events": recovery.get("events"),
        "exact_phased_cds_chain_recovered": recovery.get(
            "complete_cds_chain_recovery"
        ),
        "exact_phased_cds_chain_recall": recovery.get(
            "complete_cds_chain_recall"
        ),
        "baseline_transcript_structure_loss": collateral.get(
            "baseline_transcript_structures_missing_from_candidate"
        ),
    }


def evaluate_core_h1_scores(
    *,
    score_paths: Mapping[str, str | Path],
    holdout_id: str,
    policy_id: str,
    schema_version: str,
    primary_scope: str,
    bootstrap_seed: int,
    bootstrap_replicates: int,
    output: str | Path,
) -> dict[str, Any]:
    if set(score_paths) != EXPECTED_SCORE_KEYS:
        raise ValueError("Core-H1 evaluation requires exact three-scope/two-arm scores")
    if primary_scope != "combined" or bootstrap_replicates != 20_000:
        raise ValueError("Coffea formal H1 requires combined scope and 20,000 replicates")
    destination = Path(output)
    bootstrap_path = destination.with_name("paired_event_bootstrap.json")
    if any(path.exists() or path.is_symlink() for path in (destination, bootstrap_path)):
        raise FileExistsError("Refusing to overwrite core-H1 evaluation")
    scores = {name: load_json_object(path) for name, path in score_paths.items()}
    for name, score in scores.items():
        if not score_collateral_gate(score):
            raise ValueError(f"Core-H1 collateral/quality gate failed: {name}")
    event_counts = {
        score.get("event_recovery", {}).get("events") for score in scores.values()
    }
    if len(event_counts) != 1 or next(iter(event_counts), 0) in {None, 0}:
        raise ValueError("Core-H1 score event universes differ or are empty")
    retain_key = f"{primary_scope}_retain_distinct"
    suppress_key = f"{primary_scope}_suppress_overlap"
    bootstrap = paired_event_bootstrap(
        score_inputs=(
            ("retain_distinct", score_paths[retain_key]),
            ("suppress_overlap", score_paths[suppress_key]),
        ),
        output_json_path=bootstrap_path,
        metric="complete_cds_chain_recovery",
        replicates=bootstrap_replicates,
        seed=bootstrap_seed,
    )
    difference = bootstrap["paired_differences"][0]
    if (
        difference.get("left") != "retain_distinct"
        or difference.get("right") != "suppress_overlap"
    ):
        raise ValueError("Core-H1 bootstrap contrast direction differs")
    success = difference["observed_delta"] > 0 and difference["ci_lower"] > 0
    arms = {
        scope: {
            arm: _arm_summary(scores[f"{scope}_{arm}"])
            for arm in ARMS
        }
        for scope in SCOPES
    }
    result = {
        "schema_version": schema_version,
        "holdout_id": holdout_id,
        "policy_id": policy_id,
        "status": "ready",
        "reason_codes": [],
        "formal_outcome": (
            "formal_positive_external_result"
            if success
            else "formal_negative_external_result"
        ),
        "hypothesis": "H1_retain_distinct_vs_suppress_overlap_only",
        "primary_reference_scope": primary_scope,
        "descriptive_reference_scopes": ["bua_only", "mauritiana_only"],
        "metric": "event_exact_phased_CDS_recall_retain_distinct_minus_suppress_overlap",
        "arms": arms,
        "paired_event_bootstrap": difference,
        "bootstrap_parameters": bootstrap["parameters"],
        "bootstrap_replicates": bootstrap_replicates,
        "bootstrap_unit": "paired_event",
        "success": success,
        "gates": {
            "H1_delta_positive": difference["observed_delta"] > 0,
            "H1_CI_lower_positive": difference["ci_lower"] > 0,
            "all_arm_collateral_loss_zero": True,
        },
        "sentinels": {
            "all_arm_collateral_loss_zero": True,
            "automatic_approval": False,
        },
        "protocol_profile": "core_H1_known_subgenome_no_ranker",
        "ranker_or_model_executed": False,
        "h2_or_topology_ranking_executed": False,
        "all_arm_collateral_loss": 0,
        "inputs": {
            name: {"sha256": sha256_file(path)}
            for name, path in sorted(score_paths.items())
        },
    }
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(f".{destination.name}.working.{os.getpid()}")
    if temporary.exists() or temporary.is_symlink():
        raise FileExistsError(f"Stale core-H1 evaluation output: {temporary}")
    temporary.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    os.replace(temporary, destination)
    return result
