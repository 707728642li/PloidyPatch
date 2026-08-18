"""Fail-closed engineering contract for the Walnut core-H1 execution chain."""
from __future__ import annotations

import csv
import json
import os
from pathlib import Path
from typing import Any, Iterable

from .artifact_manifest import sha256_file, verify_sha256sums
from .holdout_contract import CORE_H1_MODEL_VERSION, load_holdout_contract, safe_relative_path


HOLDOUT_ID = "walnut_walnut2_v0.8"
POLICY_ID = "ploidypatch_walnut_external_core_h1_v0.8"
PROTOCOL_SCHEMA = "ploidypatch.walnut_core_h1_protocol_freeze.v0.8"
EXECUTION_SCHEMA = "ploidypatch.walnut_core_h1_execution_freeze.v0.8"
ROLE_SCHEMA = "ploidypatch.walnut_core_h1_blind_role_manifest.v0.8"
CUSTODY_SCHEMA = "ploidypatch.walnut_core_h1_blind_custody.v0.8"
MOUNT_SCHEMA = "ploidypatch.walnut_core_h1_mount_manifest.v0.8"
NAMESPACE_SCHEMA = "ploidypatch.walnut_core_h1_namespace_validation.v0.8"
REVEAL_STATUS_SCHEMA = "ploidypatch.walnut_h1_reveal_status.v0.8"
EVALUATION_STATUS_SCHEMA = "ploidypatch.walnut_h1_evaluation_status.v0.8"
EVALUATION_SCHEMA = "ploidypatch.walnut_external_core_h1_evaluation.v0.8"
STATUS_VALUES = frozenset({"ready", "not_evaluable", "invalid"})

PIPELINE_ENTRIES = {
    "blind_pipeline": "scripts/run_walnut_blind_pipeline_v0.8.sh",
    "reveal_input_builder": "scripts/build_walnut_complete_control_reveal_inputs_v0.8.py",
    "evaluator": "scripts/evaluate_walnut_external_h1_v0.8.py",
}
BLIND_OUTPUTS = {
    "raw_predictions_manifest": (
        "results/copy_collapse/external/walnut_v0.8_h1/"
        "raw_predictions.manifest.json"
    ),
    "retain_pool": (
        "results/copy_collapse/external/walnut_v0.8_h1/retain_distinct/"
        "blind/candidate.gff3"
    ),
    "retain_decisions": (
        "results/copy_collapse/external/walnut_v0.8_h1/retain_distinct/"
        "blind/decisions.tsv"
    ),
    "retain_manifest": (
        "results/copy_collapse/external/walnut_v0.8_h1/retain_distinct/"
        "blind/candidate.gff3.manifest.json"
    ),
    "suppress_pool": (
        "results/copy_collapse/external/walnut_v0.8_h1/suppress_overlap/"
        "blind/candidate.gff3"
    ),
    "suppress_decisions": (
        "results/copy_collapse/external/walnut_v0.8_h1/suppress_overlap/"
        "blind/decisions.tsv"
    ),
    "suppress_manifest": (
        "results/copy_collapse/external/walnut_v0.8_h1/suppress_overlap/"
        "blind/candidate.gff3.manifest.json"
    ),
    "command_log": "pipeline_commands.tsv",
}
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
    for bundle in ("candidate_mandshurica", "candidate_carya")
)
FORBIDDEN_PIPELINE_PATH_TOKENS = (
    "ranker", "stable_ranker", "topology_features", "score_candidates", "h2",
)
FORBIDDEN_BLIND_TEXT_TOKENS = (
    "/nas_data", "evaluator_only", "target_complete", "/truth", "/labels",
    "composite_model", "ranker", "topology_features", "candidate_labels",
)


def _reject_duplicate_keys(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"Duplicate JSON key: {key}")
        result[key] = value
    return result


def load_json(path: Path) -> dict[str, Any]:
    if not path.is_file() or path.is_symlink() or path.stat().st_size == 0:
        raise ValueError(f"Missing, empty or symlinked JSON: {path}")
    value = json.loads(
        path.read_text(encoding="utf-8"), object_pairs_hook=_reject_duplicate_keys
    )
    if not isinstance(value, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return value


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    if not path.is_file() or path.is_symlink() or path.stat().st_size == 0:
        raise ValueError(f"Missing, empty or symlinked TSV: {path}")
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or ())
        if not fields or len(fields) != len(set(fields)):
            raise ValueError(f"Malformed TSV header: {path}")
        return fields, list(reader)


def safe_join(root: Path, raw: str, context: str) -> Path:
    relative = safe_relative_path(raw, context)
    path = root.joinpath(*relative.parts)
    try:
        path.resolve().relative_to(root.resolve())
    except ValueError:
        raise ValueError(f"{context} escapes root") from None
    return path


def verify_protocol(protocol: Path) -> tuple[dict[str, Any], Any]:
    verify_sha256sums(protocol, ignore_checksum_file=True)
    manifest = load_json(protocol / "protocol_manifest.json")
    contract = load_holdout_contract(protocol / "contract.json")
    if (
        manifest.get("schema_version") != PROTOCOL_SCHEMA
        or manifest.get("holdout_id") != HOLDOUT_ID
        or manifest.get("policy_id") != POLICY_ID
        or contract.holdout_id != HOLDOUT_ID
        or contract.policy_id != POLICY_ID
        or contract.model_version != CORE_H1_MODEL_VERSION
        or manifest.get("model_version") != CORE_H1_MODEL_VERSION
        or manifest.get("contract_sha256") != sha256_file(protocol / "contract.json")
        or manifest.get("ranker_enabled") is not False
        or manifest.get("h2_or_topology_ranking_enabled") is not False
        or manifest.get("truth_access") is not False
        or manifest.get("wgd_pairs_enumerated") is not False
        or manifest.get("candidate_counts_computed") is not False
        or manifest.get("truth_labels_accessed") is not False
        or manifest.get("all_arm_collateral_loss_maximum") != 0
    ):
        raise ValueError("Protocol is not the exact Walnut core-H1 no-ranker freeze")
    return manifest, contract


def verify_execution(execution: Path, protocol: Path) -> tuple[dict[str, Any], dict[str, Any], Any]:
    verify_sha256sums(execution, ignore_checksum_file=True)
    protocol_manifest, contract = verify_protocol(protocol)
    manifest = load_json(execution / "execution_manifest.json")
    if (
        manifest.get("schema_version") != EXECUTION_SCHEMA
        or manifest.get("holdout_id") != HOLDOUT_ID
        or manifest.get("policy_id") != POLICY_ID
        or manifest.get("protocol_SHA256SUMS_sha256")
        != sha256_file(protocol / "SHA256SUMS")
        or manifest.get("contract_sha256") != sha256_file(protocol / "contract.json")
        or manifest.get("pipeline_entries") != PIPELINE_ENTRIES
        or manifest.get("blind_outputs") != BLIND_OUTPUTS
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
        raise ValueError("Execution is not the exact Walnut core-H1 freeze")
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
        raise ValueError("Malformed tri-state status")
    return str(status)


def atomic_write_json(path: Path, value: dict[str, Any]) -> None:
    if path.exists() or path.is_symlink():
        raise FileExistsError(f"Refusing to overwrite JSON: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.partial.{os.getpid()}")
    if temporary.exists() or temporary.is_symlink():
        raise FileExistsError(f"Refusing to overwrite partial JSON: {temporary}")
    try:
        temporary.write_text(
            json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        os.replace(temporary, path)
    finally:
        if temporary.exists() and not temporary.is_symlink():
            temporary.unlink()


def reject_forbidden_text(path: Path, extra: Iterable[str] = ()) -> None:
    text = path.read_text(encoding="utf-8", errors="strict").casefold()
    if not text.strip():
        raise ValueError(f"Empty command/audit text: {path}")
    forbidden = tuple(token.casefold() for token in (*FORBIDDEN_BLIND_TEXT_TOKENS, *extra))
    hits = [token for token in forbidden if token in text]
    if hits:
        raise ValueError(f"Forbidden blind text in {path}: {hits}")
