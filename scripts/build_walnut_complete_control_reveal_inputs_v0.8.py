#!/usr/bin/env python3
"""Build Walnut complete controls only after blind custody and reveal authorization."""
from __future__ import annotations

import importlib.util
import json
import os
from pathlib import Path
import shutil
import sys
import tempfile
from typing import Any

from ploidypatch.artifact_manifest import sha256_file, verify_sha256sums, write_sha256sums
from ploidypatch.holdout_contract import load_holdout_contract
from ploidypatch.score import score_annotation_repair
from ploidypatch.walnut_h1 import (
    HOLDOUT_ID,
    POLICY_ID,
    load_json_object,
    score_collateral_gate,
)
from ploidypatch.walnut_h1_framework import BLIND_OUTPUTS, CUSTODY_SCHEMA


STATUS_SCHEMA = "ploidypatch.walnut_h1_reveal_status.v0.8"
REVEAL_SCHEMA = "ploidypatch.walnut_h1_reveal_inputs.v0.8"
NOT_EVALUABLE_REASONS = frozenset(
    {
        "formal_event_count_below_500",
        "target_primary_chromosome_count_below_12",
        "complexity_bin_one_below_20",
        "complexity_bin_two_to_three_below_20",
        "complexity_bin_four_to_six_below_20",
        "complexity_bin_seven_plus_below_20",
    }
)


def required(name: str) -> Path:
    value = os.environ.get(name)
    if not value:
        raise ValueError(f"{name} is required")
    return Path(value).resolve()


def write_json(path: Path, value: dict[str, Any]) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def load_pool_builder(execution: Path) -> Any:
    path = execution / "source/scripts/build_walnut_h1_candidate_pools_v0.8.py"
    spec = importlib.util.spec_from_file_location("walnut_h1_candidate_pools_v0_8", path)
    if spec is None or spec.loader is None:
        raise ValueError("Frozen Walnut pool builder cannot be imported")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def validate_blind_custody(blind_run: Path, custody_path: Path) -> dict[str, Path]:
    verify_sha256sums(blind_run, ignore_checksum_file=True)
    custody = load_json_object(custody_path)
    if (
        custody.get("schema_version") != CUSTODY_SCHEMA
        or custody.get("holdout_id") != HOLDOUT_ID
        or custody.get("policy_id") != POLICY_ID
        or custody.get("ranker_or_model_executed") is not False
        or custody.get("h2_or_topology_ranking_executed") is not False
        or any(
            custody.get(field) is not False
            for field in (
                "truth_mounted",
                "complete_target_annotation_mounted",
                "evaluator_references_mounted",
                "nas_data_mounted",
                "network_access",
            )
        )
    ):
        raise ValueError("Walnut blind custody firewall differs")
    project = blind_run / "project"
    values = custody.get("blind_outputs")
    if not isinstance(values, dict) or set(values) != set(BLIND_OUTPUTS):
        raise ValueError("Walnut custody lacks exact eight blind outputs")
    resolved: dict[str, Path] = {}
    for name, expected_relative in BLIND_OUTPUTS.items():
        item = values.get(name)
        if not isinstance(item, dict) or set(item) != {
            "relative_path", "sha256", "bytes"
        }:
            raise ValueError(f"Malformed Walnut custody binding: {name}")
        if item["relative_path"] != expected_relative:
            raise ValueError(f"Walnut custody path differs: {name}")
        path = project / Path(expected_relative)
        if (
            not path.is_file()
            or path.is_symlink()
            or path.stat().st_size != item["bytes"]
            or sha256_file(path) != item["sha256"]
        ):
            raise ValueError(f"Walnut custody hash differs: {name}")
        resolved[name] = path
    return resolved


def publish_status(
    working: Path,
    *,
    status: str,
    reasons: list[str],
    custody_sha: str,
    inputs: dict[str, dict[str, str]] | None = None,
) -> None:
    write_json(
        working / "status.json",
        {
            "schema_version": STATUS_SCHEMA,
            "holdout_id": HOLDOUT_ID,
            "policy_id": POLICY_ID,
            "status": status,
            "reason_codes": reasons,
        },
    )
    write_json(
        working / "reveal_input_manifest.json",
        {
            "schema_version": REVEAL_SCHEMA,
            "holdout_id": HOLDOUT_ID,
            "policy_id": POLICY_ID,
            "formal_status": status,
            "reason_codes": reasons,
            "custody_manifest_sha256": custody_sha,
            "ranker_or_model_access": False,
            "h2_or_topology_ranking_access": False,
            "raw_predictions_reused_without_rerun": True,
            "evaluation_inputs": inputs or {},
        },
    )


def main(argv: list[str] | None = None) -> int:
    values = sys.argv[1:] if argv is None else argv
    if len(values) != 1:
        raise SystemExit("usage: build_walnut_complete_control_reveal_inputs_v0.8.py PROJECT_ROOT")
    project_root = Path(values[0]).resolve()
    contract_path = required("PLOIDYPATCH_HOLDOUT_CONTRACT")
    protocol = required("PLOIDYPATCH_PROTOCOL_FREEZE")
    execution = required("PLOIDYPATCH_EXECUTION_FREEZE")
    blind_run = required("PLOIDYPATCH_BLIND_RUN_ROOT")
    custody_path = required("PLOIDYPATCH_CUSTODY_MANIFEST")
    authorization_path = required("PLOIDYPATCH_REVEAL_AUTHORIZATION")
    evaluator = required("PLOIDYPATCH_EVALUATOR_ONLY_ROOT")
    output = required("PLOIDYPATCH_REVEAL_INPUTS_OUTPUT")
    if output.exists() or output.is_symlink():
        raise FileExistsError(f"Refusing to overwrite Walnut reveal inputs: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    working = Path(tempfile.mkdtemp(prefix=f".{output.name}.working.", dir=output.parent))
    exit_status = 0
    custody_sha = sha256_file(custody_path)
    try:
        for root in (protocol, execution, evaluator):
            verify_sha256sums(root, ignore_checksum_file=True)
        contract = load_holdout_contract(contract_path)
        if contract.holdout_id != HOLDOUT_ID or contract.policy_id != POLICY_ID:
            raise ValueError("Walnut reveal builder received a different contract")
        authorization = load_json_object(authorization_path)
        if (
            authorization.get("schema_version")
            != "ploidypatch.walnut_h1_reveal_authorization.v0.8"
            or authorization.get("holdout_id") != HOLDOUT_ID
            or authorization.get("truth_reveal_authorized") is not True
            or authorization.get("custody_manifest_sha256") != custody_sha
            or authorization.get("ranker_or_model_authorized") is not False
            or authorization.get("h2_or_topology_ranking_authorized") is not False
        ):
            raise ValueError("Walnut reveal authorization differs")
        blind = validate_blind_custody(blind_run, custody_path)
        raw_manifest = load_json_object(blind["raw_predictions_manifest"])
        raw_hashes = raw_manifest.get("input_hashes")
        if not isinstance(raw_hashes, dict):
            raise ValueError("Walnut raw-prediction input hashes are missing")
        source_genome = evaluator / "normalized/target_complete/primary_chromosomes.genome.fa"
        evaluator_perturbed = evaluator / "benchmark/inputs/perturbed.gff3"
        if (
            not source_genome.is_file()
            or source_genome.is_symlink()
            or raw_hashes.get("target_genome_sha256") != sha256_file(source_genome)
            or not evaluator_perturbed.is_file()
            or evaluator_perturbed.is_symlink()
            or raw_hashes.get("perturbed_gff3_sha256")
            != sha256_file(evaluator_perturbed)
        ):
            raise ValueError("Walnut blind/evaluator genome or perturbation sentinel differs")
        evaluability_path = evaluator / "benchmark/pair_selection/evaluability.json"
        evaluability = load_json_object(evaluability_path)
        status = evaluability.get("status")
        reasons = evaluability.get("reason_codes")
        if status not in {"ready", "not_evaluable", "invalid"} or not isinstance(reasons, list):
            raise ValueError("Walnut evaluability status is malformed")
        if status == "not_evaluable" and (
            not reasons or not set(reasons) <= NOT_EVALUABLE_REASONS
        ):
            raise ValueError("Walnut not-evaluable reasons are not fixed data gates")
        if status == "ready" and reasons:
            raise ValueError("Ready Walnut evaluability has reasons")
        sentinels = evaluability.get("sentinels")
        event_count = evaluability.get("events")
        if (
            not isinstance(sentinels, dict)
            or sentinels.get("blind_noop_exact_recovery") != 0
            or sentinels.get("complete_oracle_exact_recovery") != event_count
            or sentinels.get("restoration_byte_identical") is not True
            or sentinels.get("blind_complete_genome_sha256_identical") is not True
            or sentinels.get("noop_quality_grade_pass") is not True
            or sentinels.get("oracle_quality_grade_pass") is not True
        ):
            raise ValueError("Walnut no-op/oracle/restoration/genome sentinel differs")
        shutil.copyfile(evaluability_path, working / "evaluability.json")
        if status != "ready":
            publish_status(
                working,
                status=status,
                reasons=reasons or ["evaluator_holdout_invalid"],
                custody_sha=custody_sha,
            )
            exit_status = 1 if status == "invalid" else 0
        else:
            source_gff = evaluator / "normalized/target_complete/primary_chromosomes.gff3"
            # Evaluator-owned duplicate of the exact perturbed bytes is sealed
            # before blind execution; the host reveal interface need not remount
            # the blind-role root.
            blind_gff = evaluator_perturbed
            truth = evaluator / "benchmark/truth/hidden_truth.json"
            for path in (source_gff, blind_gff, truth):
                if not path.is_file() or path.is_symlink() or path.stat().st_size == 0:
                    raise ValueError(f"Walnut reveal science input is missing: {path}")
            builder = load_pool_builder(execution)
            complete = working / "complete_control"
            builder.build_pools(
                project_root=blind_run / "project",
                base_gff=source_gff,
                output=complete,
                include_raw_manifest=False,
                raw_manifest_path=blind["raw_predictions_manifest"],
                seal_manifests=False,
            )
            scores = working / "scores"
            scores.mkdir()
            bindings: dict[str, dict[str, str]] = {}
            for arm, key in (("retain_distinct", "retain_score"), ("suppress_overlap", "suppress_score")):
                candidate_key = "retain_pool" if arm == "retain_distinct" else "suppress_pool"
                control = complete / arm / "complete_control/candidate.gff3"
                score = score_annotation_repair(
                    source_gff_path=source_gff,
                    perturbed_gff_path=blind_gff,
                    candidate_gff_path=blind[candidate_key],
                    truth_path=truth,
                    include_event_details=True,
                    control_candidate_gff_path=control,
                )
                if not score_collateral_gate(score):
                    raise ValueError(f"Walnut {arm} collateral loss is nonzero")
                score_path = scores / f"{arm}.json"
                write_json(score_path, score)
                bindings[key] = {
                    "relative_path": score_path.relative_to(working).as_posix(),
                    "sha256": sha256_file(score_path),
                }
            bindings["evaluability"] = {
                "relative_path": "evaluability.json",
                "sha256": sha256_file(working / "evaluability.json"),
            }
            publish_status(
                working,
                status="ready",
                reasons=[],
                custody_sha=custody_sha,
                inputs=bindings,
            )
        write_sha256sums(working)
        verify_sha256sums(working, ignore_checksum_file=True)
        os.replace(working, output)
        return exit_status
    except BaseException as error:
        exit_status = 1
        for child in list(working.iterdir()):
            if child.is_dir():
                shutil.rmtree(child)
            else:
                child.unlink()
        publish_status(
            working,
            status="invalid",
            reasons=[f"{type(error).__name__}:{error}"],
            custody_sha=custody_sha,
        )
        write_sha256sums(working)
        verify_sha256sums(working, ignore_checksum_file=True)
        os.replace(working, output)
        return exit_status


if __name__ == "__main__":
    raise SystemExit(main())
