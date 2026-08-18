"""Exact engineering contract for the Coffea core-H1 external validation."""
from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

from .artifact_manifest import sha256_file, verify_sha256sums
from .holdout_contract import KNOWN_SUBGENOME_CORE_H1_MODEL_VERSION, load_holdout_contract


HOLDOUT_ID = "coffea_et39_v1.0"
POLICY_ID = "ploidypatch_coffea_external_core_h1_v1.0"
PROTOCOL_SCHEMA = "ploidypatch.coffea_core_h1_protocol_freeze.v1.0"
EXECUTION_SCHEMA = "ploidypatch.coffea_core_h1_execution_freeze.v1.0"
ROLE_SCHEMA = "ploidypatch.coffea_core_h1_blind_role_manifest.v1.0"
EVALUATOR_ROLE_SCHEMA = "ploidypatch.coffea_core_h1_evaluator_role_manifest.v1.0"
CUSTODY_SCHEMA = "ploidypatch.coffea_core_h1_blind_custody.v1.0"
MOUNT_SCHEMA = "ploidypatch.coffea_core_h1_mount_manifest.v1.0"
NAMESPACE_SCHEMA = "ploidypatch.coffea_core_h1_namespace_validation.v1.0"
REVEAL_STATUS_SCHEMA = "ploidypatch.coffea_h1_reveal_status.v1.0"
REVEAL_INPUT_SCHEMA = "ploidypatch.coffea_h1_reveal_inputs.v1.0"
EVALUATION_STATUS_SCHEMA = "ploidypatch.coffea_h1_evaluation_status.v1.0"
EVALUATION_SCHEMA = "ploidypatch.coffea_external_core_h1_evaluation.v1.0"
RAW_SCHEMA = "ploidypatch.core_h1_raw_predictions.v1"
STATUS_VALUES = frozenset({"ready", "not_evaluable", "invalid"})
PATCH_STAGE = "post_evaluator_truth_failed_blind_pre_candidate_pre_label_execution_patch"
PATCH2_STAGE = (
    "post_evaluator_truth_failed_blind_partial_candidate_pre_label_execution_patch_2"
)
PATCH3_STAGE = (
    "post_evaluator_truth_two_complete_blind_runs_pre_label_"
    "reproducibility_patch_3"
)
PATCH4_STAGE = (
    "post_blind_custody_reveal_authorized_pre_evaluator_environment_patch_4"
)
PATCH5_STAGE = (
    "post_blind_custody_pre_truth_authorization_custody_lineage_patch_5"
)

PIPELINE_ENTRIES = {
    "blind_pipeline": "scripts/run_coffea_blind_pipeline_v1.0.sh",
    "reveal_input_builder": "scripts/build_coffea_complete_control_reveal_inputs_v1.0.py",
    "evaluator": "scripts/evaluate_coffea_external_h1_v1.0.py",
}


def _blind_outputs() -> dict[str, str]:
    root = "results/copy_collapse/external/coffea_v1.0_h1"
    result = {"raw_predictions_manifest": f"{root}/raw_predictions.manifest.json"}
    for scope in ("combined", "bua_only", "mauritiana_only"):
        for arm in ("retain_distinct", "suppress_overlap"):
            prefix = f"{scope}_{arm}"
            base = f"{root}/{scope}/{arm}/blind"
            result[f"{prefix}_pool"] = f"{base}/candidate.gff3"
            result[f"{prefix}_decisions"] = f"{base}/decisions.tsv"
            result[f"{prefix}_manifest"] = f"{base}/candidate.gff3.manifest.json"
    result["command_log"] = "pipeline_commands.tsv"
    return result


BLIND_OUTPUTS = _blind_outputs()
REQUIRED_ENVIRONMENTS = frozenset(
    {
        "ploidypatch-dev",
        "ploidypatch-baseline",
        "ploidypatch-synteny",
        "ploidypatch-syngap",
        "ploidypatch-gemoma",
        "ploidypatch-lifton",
    }
)
RAW_PREDICTION_TREE_KEYS = frozenset(
    f"{method}__{bundle}"
    for method in ("miniprot", "gemoma", "lifton")
    for bundle in ("candidate_bua", "candidate_mauritiana")
)


def load_json(path: str | Path) -> dict[str, Any]:
    source = Path(path)
    if not source.is_file() or source.is_symlink() or source.stat().st_size == 0:
        raise ValueError(f"Missing, empty or symlinked JSON: {source}")

    def reject_duplicates(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        value: dict[str, Any] = {}
        for key, item in pairs:
            if key in value:
                raise ValueError(f"Duplicate JSON key in {source}: {key}")
            value[key] = item
        return value

    value = json.loads(
        source.read_text(encoding="utf-8"), object_pairs_hook=reject_duplicates
    )
    if not isinstance(value, dict):
        raise ValueError(f"Expected JSON object: {source}")
    return value


def read_tsv(path: str | Path) -> tuple[list[str], list[dict[str, str]]]:
    source = Path(path)
    if not source.is_file() or source.is_symlink() or source.stat().st_size == 0:
        raise ValueError(f"Missing, empty or symlinked TSV: {source}")
    with source.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or ())
        if not fields or len(fields) != len(set(fields)):
            raise ValueError(f"Malformed TSV header: {source}")
        return fields, list(reader)


def verify_protocol(protocol: str | Path) -> tuple[dict[str, Any], Any]:
    root = Path(protocol)
    verify_sha256sums(root, ignore_checksum_file=True)
    manifest = load_json(root / "protocol_manifest.json")
    contract = load_holdout_contract(root / "contract.json")
    if (
        manifest.get("schema_version") != PROTOCOL_SCHEMA
        or manifest.get("holdout_id") != HOLDOUT_ID
        or manifest.get("policy_id") != POLICY_ID
        or contract.holdout_id != HOLDOUT_ID
        or contract.policy_id != POLICY_ID
        or contract.model_version != KNOWN_SUBGENOME_CORE_H1_MODEL_VERSION
        or manifest.get("model_version") != KNOWN_SUBGENOME_CORE_H1_MODEL_VERSION
        or manifest.get("contract_sha256") != sha256_file(root / "contract.json")
        or manifest.get("protocol_profile")
        != "core_H1_known_subgenome_no_ranker"
        or manifest.get("ranker_enabled") is not False
        or manifest.get("h2_or_topology_ranking_enabled") is not False
        or manifest.get("truth_access") is not False
        or manifest.get("wgd_pairs_enumerated") is not False
        or manifest.get("candidate_counts_computed") is not False
        or manifest.get("truth_labels_accessed") is not False
    ):
        raise ValueError("Protocol is not the exact Coffea core-H1 freeze")
    return manifest, contract


def verify_execution(
    execution: str | Path, protocol: str | Path
) -> tuple[dict[str, Any], dict[str, Any], Any]:
    root = Path(execution)
    protocol_root = Path(protocol)
    verify_sha256sums(root, ignore_checksum_file=True)
    protocol_manifest, contract = verify_protocol(protocol_root)
    manifest = load_json(root / "execution_manifest.json")
    patch = manifest.get("execution_patch")
    stage = manifest.get("freeze_stage")
    original_stage = "post_metadata_pre_pair_pre_candidate_pre_label"
    stage_valid = False
    if stage == original_stage and patch is None:
        stage_valid = (
            manifest.get("truth_pairs_enumerated_before_execution_freeze") is False
            and manifest.get("candidate_counts_computed_before_execution_freeze") is False
            and manifest.get("truth_labels_accessed_before_execution_freeze") is False
        )
    elif stage in {
        PATCH_STAGE, PATCH2_STAGE, PATCH3_STAGE, PATCH4_STAGE, PATCH5_STAGE
    } and isinstance(patch, dict):
        patch_sequence = patch.get("patch_sequence", 1)
        common_patch_valid = (
            patch.get("schema_version")
            == "ploidypatch.coffea_core_h1_execution_patch.v1.0"
            and patch.get("freeze_stage") == stage
            and patch.get("patch_code_commit") == manifest.get("code_commit")
            and isinstance(patch.get("base_code_commit"), str)
            and len(patch["base_code_commit"]) == 40
            and patch.get("evaluator_truth_construction_completed_before_patch") is True
            and patch.get("formal_scores_generated_before_patch") is False
            and patch.get("truth_labels_accessed_before_patch") is False
            and manifest.get("truth_pairs_enumerated_before_execution_freeze") is True
            and manifest.get("truth_labels_accessed_before_execution_freeze") is False
        )
        if stage == PATCH_STAGE:
            stage_valid = (
                common_patch_valid
                and patch_sequence == 1
                and patch.get("failed_attempt_exit_status") not in (None, 0)
                and patch.get("candidate_generation_completed_before_patch") is False
                and patch.get("scientific_rules_or_thresholds_changed") is False
                and manifest.get("candidate_counts_computed_before_execution_freeze")
                is False
            )
        elif stage == PATCH2_STAGE:
            replay = patch.get("primary_combined_replay_outputs")
            stage_valid = (
                common_patch_valid
                and patch_sequence == 2
                and patch.get("failed_attempt_exit_status") not in (None, 0)
                and patch.get("candidate_generation_completed_before_patch") is False
                and patch.get("scientific_rules_or_thresholds_changed") is False
                and patch.get("candidate_generation_completed_before_patch") is False
                and patch.get("partial_candidate_generation_before_patch") is True
                and patch.get("formal_blind_outputs_frozen_before_patch") is False
                and patch.get("blind_custody_completed_before_patch") is False
                and manifest.get("candidate_counts_computed_before_execution_freeze")
                is True
                and isinstance(replay, dict)
                and set(replay)
                == {
                    "combined_retain_distinct_pool",
                    "combined_retain_distinct_decisions",
                    "combined_suppress_overlap_pool",
                    "combined_suppress_overlap_decisions",
                }
                and all(
                    isinstance(value, str) and len(value) == 64
                    for value in replay.values()
                )
            )
        elif stage == PATCH3_STAGE:
            run_manifest_keys = (
                "reproducibility_run_a_manifest_sha256",
                "reproducibility_run_b_manifest_sha256",
            )
            stage_valid = (
                common_patch_valid
                and patch_sequence == 3
                and patch.get("failed_attempt_exit_status") == 0
                and patch.get("candidate_generation_completed_before_patch") is True
                and patch.get("two_blind_candidate_executions_completed_before_patch")
                is True
                and patch.get("blind_custody_completed_before_patch") is False
                and patch.get("scientific_rules_or_thresholds_changed") is True
                and patch.get("projection_reproducibility_abstention_rule_added")
                is True
                and patch.get("biological_rules_or_performance_thresholds_changed")
                is False
                and patch.get("label_informed_selection") is False
                and patch.get("unstable_projection_policy") == "abstain"
                and all(
                    isinstance(patch.get(key), str) and len(patch[key]) == 64
                    for key in run_manifest_keys
                )
                and all(
                    isinstance(patch.get(key), str) and patch[key]
                    for key in (
                        "reproducibility_run_a_manifest_relative_path",
                        "reproducibility_run_b_manifest_relative_path",
                    )
                )
                and manifest.get("candidate_counts_computed_before_execution_freeze")
                is True
            )
        elif stage == PATCH4_STAGE:
            stage_valid = (
                common_patch_valid
                and patch_sequence == 4
                and patch.get("failed_attempt_exit_status") == 1
                and patch.get("candidate_generation_completed_before_patch") is True
                and patch.get("blind_custody_completed_before_patch") is True
                and patch.get("formal_blind_outputs_frozen_before_patch") is True
                and patch.get("truth_reveal_authorized_before_patch") is True
                and patch.get("evaluator_truth_bytes_hashed_before_patch") is True
                and patch.get("evaluator_invoked_before_patch") is False
                and patch.get("performance_metrics_computed_before_patch") is False
                and patch.get("scientific_rules_or_thresholds_changed") is False
                and patch.get("environment_interpreter_symlink_validation_fixed")
                is True
                and patch.get("canonical_nested_manifest_writer_fixed") is True
                and isinstance(patch.get("blind_custody_manifest_sha256"), str)
                and len(patch["blind_custody_manifest_sha256"]) == 64
                and isinstance(patch.get("blind_root_SHA256SUMS_sha256"), str)
                and len(patch["blind_root_SHA256SUMS_sha256"]) == 64
                and isinstance(patch.get("failed_reveal_error_log_sha256"), str)
                and len(patch["failed_reveal_error_log_sha256"]) == 64
                and manifest.get("candidate_counts_computed_before_execution_freeze")
                is True
            )
        else:
            stage_valid = (
                common_patch_valid
                and patch_sequence == 5
                and patch.get("failed_attempt_exit_status") == 1
                and patch.get("candidate_generation_completed_before_patch") is True
                and patch.get("blind_custody_completed_before_patch") is True
                and patch.get("formal_blind_outputs_frozen_before_patch") is True
                and patch.get("truth_reveal_authorized_before_patch") is False
                and patch.get("evaluator_truth_bytes_hashed_before_patch") is False
                and patch.get("evaluator_invoked_before_patch") is False
                and patch.get("performance_metrics_computed_before_patch") is False
                and patch.get("scientific_rules_or_thresholds_changed") is False
                and patch.get("custody_execution_chain_validation_fixed") is True
                and all(
                    isinstance(patch.get(key), str) and len(patch[key]) == 64
                    for key in (
                        "blind_custody_manifest_sha256",
                        "blind_custody_execution_SHA256SUMS_sha256",
                        "blind_root_SHA256SUMS_sha256",
                        "failed_reveal_error_log_sha256",
                    )
                )
                and manifest.get("candidate_counts_computed_before_execution_freeze")
                is True
            )
    if (
        manifest.get("schema_version") != EXECUTION_SCHEMA
        or manifest.get("holdout_id") != HOLDOUT_ID
        or manifest.get("policy_id") != POLICY_ID
        or manifest.get("protocol_SHA256SUMS_sha256")
        != sha256_file(protocol_root / "SHA256SUMS")
        or manifest.get("contract_sha256")
        != sha256_file(protocol_root / "contract.json")
        or manifest.get("pipeline_entries") != PIPELINE_ENTRIES
        or manifest.get("blind_outputs") != BLIND_OUTPUTS
        or not stage_valid
        or manifest.get("ranker_or_model_execution") is not False
        or manifest.get("h2_or_topology_ranking_enabled") is not False
        or any(
            manifest.get(name) is not False
            for name in (
                "network_access_in_blind_runner",
                "nas_data_mount_in_blind_runner",
                "complete_target_annotation_mount_in_blind_runner",
                "evaluator_only_mount_in_blind_runner",
                "truth_or_label_mount_in_blind_runner",
            )
        )
    ):
        raise ValueError("Execution is not the exact Coffea core-H1 freeze")
    return manifest, protocol_manifest, contract


def validate_status(value: dict[str, Any], *, expected_schema: str) -> str:
    status = value.get("status")
    reasons = value.get("reason_codes")
    if (
        value.get("schema_version") != expected_schema
        or status not in STATUS_VALUES
        or not isinstance(reasons, list)
        or not all(isinstance(reason, str) and reason for reason in reasons)
        or (status == "ready" and reasons)
        or (status != "ready" and not reasons)
    ):
        raise ValueError("Malformed Coffea tri-state status")
    return str(status)
