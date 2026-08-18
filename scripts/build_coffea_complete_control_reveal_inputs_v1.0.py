#!/usr/bin/env python3
"""Build Coffea complete controls only after frozen blind custody."""
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
from ploidypatch.coffea_h1_framework import (
    BLIND_OUTPUTS,
    CUSTODY_SCHEMA,
    EVALUATOR_ROLE_SCHEMA,
    HOLDOUT_ID,
    POLICY_ID,
    REVEAL_INPUT_SCHEMA,
    REVEAL_STATUS_SCHEMA,
    load_json,
    verify_execution,
)
from ploidypatch.core_h1_evaluation import score_collateral_gate
from ploidypatch.holdout_contract import load_holdout_contract
from ploidypatch.score import score_annotation_repair


NOT_EVALUABLE_REASONS = frozenset(
    {
        "formal_event_count_below_500",
        "target_primary_chromosome_count_below_17",
        "homoeolog_group_count_below_9",
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
    path = execution / "source/scripts/build_coffea_h1_candidate_pools_v1.0.py"
    spec = importlib.util.spec_from_file_location("coffea_h1_candidate_pools_v1_0", path)
    if spec is None or spec.loader is None:
        raise ValueError("Frozen Coffea pool builder cannot be imported")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def validate_blind_custody(blind_run: Path, custody_path: Path) -> dict[str, Path]:
    verify_sha256sums(blind_run, ignore_checksum_file=True)
    custody = load_json(custody_path)
    if (
        custody.get("schema_version") != CUSTODY_SCHEMA
        or custody.get("holdout_id") != HOLDOUT_ID
        or custody.get("policy_id") != POLICY_ID
        or custody.get("ranker_or_model_executed") is not False
        or custody.get("h2_or_topology_ranking_executed") is not False
        or any(custody.get(field) is not False for field in (
            "truth_mounted", "complete_target_annotation_mounted",
            "evaluator_references_mounted", "nas_data_mounted", "network_access",
        ))
    ):
        raise ValueError("Coffea blind custody firewall differs")
    project = blind_run / "project"
    values = custody.get("blind_outputs")
    if not isinstance(values, dict) or set(values) != set(BLIND_OUTPUTS):
        raise ValueError("Coffea custody lacks exact blind output universe")
    resolved: dict[str, Path] = {}
    for name, expected_relative in BLIND_OUTPUTS.items():
        item = values.get(name)
        if not isinstance(item, dict) or set(item) != {"relative_path", "sha256", "bytes"}:
            raise ValueError(f"Malformed Coffea custody binding: {name}")
        path = project / Path(expected_relative)
        if (
            item["relative_path"] != expected_relative
            or not path.is_file()
            or path.is_symlink()
            or path.stat().st_size != item["bytes"]
            or sha256_file(path) != item["sha256"]
        ):
            raise ValueError(f"Coffea custody output differs: {name}")
        resolved[name] = path
    return resolved


def publish_status(
    working: Path, *, status: str, reasons: list[str], custody_sha: str,
    inputs: dict[str, dict[str, str]] | None = None
) -> None:
    write_json(
        working / "status.json",
        {
            "schema_version": REVEAL_STATUS_SCHEMA,
            "holdout_id": HOLDOUT_ID,
            "policy_id": POLICY_ID,
            "status": status,
            "reason_codes": reasons,
        },
    )
    write_json(
        working / "reveal_input_manifest.json",
        {
            "schema_version": REVEAL_INPUT_SCHEMA,
            "holdout_id": HOLDOUT_ID,
            "policy_id": POLICY_ID,
            "formal_status": status,
            "reason_codes": reasons,
            "custody_manifest_sha256": custody_sha,
            "ranker_or_model_access": False,
            "h2_or_topology_ranking_access": False,
            "raw_predictions_reused_without_rerun": True,
            "primary_reference_scope": "combined",
            "evaluation_inputs": inputs or {},
        },
    )


def main(argv: list[str] | None = None) -> int:
    values = sys.argv[1:] if argv is None else argv
    if len(values) != 1:
        raise SystemExit("usage: build_coffea_complete_control_reveal_inputs_v1.0.py PROJECT_ROOT")
    project_root = Path(values[0]).resolve()
    contract_path = required("PLOIDYPATCH_HOLDOUT_CONTRACT")
    protocol = required("PLOIDYPATCH_PROTOCOL_FREEZE")
    execution = required("PLOIDYPATCH_EXECUTION_FREEZE")
    blind_run = required("PLOIDYPATCH_BLIND_RUN_ROOT")
    custody_path = required("PLOIDYPATCH_CUSTODY_MANIFEST")
    authorization_path = required("PLOIDYPATCH_REVEAL_AUTHORIZATION")
    evaluator_input = required("PLOIDYPATCH_EVALUATOR_INPUT_ROOT")
    evaluator = required("PLOIDYPATCH_EVALUATOR_ONLY_ROOT")
    output = required("PLOIDYPATCH_REVEAL_INPUTS_OUTPUT")
    if output.exists() or output.is_symlink():
        raise FileExistsError(f"Refusing to overwrite Coffea reveal inputs: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    working = Path(tempfile.mkdtemp(prefix=f".{output.name}.working.", dir=output.parent))
    custody_sha = sha256_file(custody_path)
    exit_status = 0
    try:
        verify_execution(execution, protocol)
        verify_sha256sums(evaluator_input, ignore_checksum_file=True)
        verify_sha256sums(evaluator, ignore_checksum_file=True)
        contract = load_holdout_contract(contract_path)
        evaluator_role = load_json(evaluator_input / "role_manifest.json")
        if (
            contract.holdout_id != HOLDOUT_ID
            or contract.policy_id != POLICY_ID
            or evaluator_role.get("schema_version") != EVALUATOR_ROLE_SCHEMA
            or evaluator_role.get("candidate_reference_access") is not False
            or evaluator_role.get("blind_candidate_outputs_access") is not False
        ):
            raise ValueError("Coffea reveal evaluator role differs")
        authorization = load_json(authorization_path)
        if (
            authorization.get("schema_version")
            != "ploidypatch.coffea_h1_reveal_authorization.v1.0"
            or authorization.get("holdout_id") != HOLDOUT_ID
            or authorization.get("truth_reveal_authorized") is not True
            or authorization.get("custody_manifest_sha256") != custody_sha
            or authorization.get("ranker_or_model_authorized") is not False
            or authorization.get("h2_or_topology_ranking_authorized") is not False
        ):
            raise ValueError("Coffea reveal authorization differs")
        blind = validate_blind_custody(blind_run, custody_path)
        raw = load_json(blind["raw_predictions_manifest"])
        raw_hashes = raw.get("input_hashes")
        source_genome = evaluator_input / "normalized/target_complete/primary_chromosomes.genome.fa"
        evaluator_perturbed = evaluator / "benchmark/inputs/perturbed.gff3"
        if (
            not isinstance(raw_hashes, dict)
            or raw_hashes.get("target_genome_sha256") != sha256_file(source_genome)
            or raw_hashes.get("perturbed_gff3_sha256") != sha256_file(evaluator_perturbed)
        ):
            raise ValueError("Coffea blind/evaluator genome or perturbation sentinel differs")
        evaluability_path = evaluator / "benchmark/pair_selection/evaluability.json"
        evaluability = load_json(evaluability_path)
        status = evaluability.get("status")
        reasons = evaluability.get("reason_codes")
        if status not in {"ready", "not_evaluable", "invalid"} or not isinstance(reasons, list):
            raise ValueError("Coffea evaluability status is malformed")
        if status == "not_evaluable" and (not reasons or not set(reasons) <= NOT_EVALUABLE_REASONS):
            raise ValueError("Coffea not-evaluable reasons differ from frozen gates")
        if status == "ready" and reasons:
            raise ValueError("Ready Coffea evaluability contains reason codes")
        sentinels = evaluability.get("sentinels", {})
        event_count = evaluability.get("events")
        if (
            sentinels.get("blind_noop_exact_recovery") != 0
            or sentinels.get("complete_oracle_exact_recovery") != event_count
            or sentinels.get("restoration_byte_identical") is not True
            or sentinels.get("blind_complete_genome_sha256_identical") is not True
            or sentinels.get("noop_quality_grade_pass") is not True
            or sentinels.get("oracle_quality_grade_pass") is not True
        ):
            raise ValueError("Coffea no-op/oracle/restoration/genome sentinel differs")
        shutil.copyfile(evaluability_path, working / "evaluability.json")
        if status != "ready":
            publish_status(
                working, status=status,
                reasons=reasons or ["evaluator_holdout_invalid"], custody_sha=custody_sha,
            )
            exit_status = 1 if status == "invalid" else 0
        else:
            source_gff = evaluator_input / "normalized/target_complete/primary_chromosomes.gff3"
            truth = evaluator / "benchmark/truth/hidden_truth.json"
            for path in (source_gff, evaluator_perturbed, truth):
                if not path.is_file() or path.is_symlink() or path.stat().st_size == 0:
                    raise ValueError(f"Coffea reveal science input is missing: {path}")
            builder = load_pool_builder(execution)
            builder.build_pools(
                project_root=blind_run / "project",
                raw_baselines=blind_run / "project/results/baselines/coffea/v1.0",
                base_gff=source_gff,
                output=working / "complete_control",
                include_raw_manifest=False,
                raw_manifest_path=blind["raw_predictions_manifest"],
                seal_manifests=False,
            )
            scores = working / "scores"; scores.mkdir()
            bindings: dict[str, dict[str, str]] = {}
            for scope in ("combined", "bua_only", "mauritiana_only"):
                for arm in ("retain_distinct", "suppress_overlap"):
                    key = f"{scope}_{arm}"
                    blind_pool = blind[f"{key}_pool"]
                    control = (
                        working / "complete_control" / scope / arm
                        / "complete_control/candidate.gff3"
                    )
                    score = score_annotation_repair(
                        source_gff_path=source_gff,
                        perturbed_gff_path=evaluator_perturbed,
                        candidate_gff_path=blind_pool,
                        truth_path=truth,
                        include_event_details=True,
                        control_candidate_gff_path=control,
                    )
                    if not score_collateral_gate(score):
                        raise ValueError(f"Coffea {key} collateral loss is nonzero")
                    score_path = scores / f"{key}.json"
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
                working, status="ready", reasons=[], custody_sha=custody_sha,
                inputs=bindings,
            )
        write_sha256sums(working)
        verify_sha256sums(working, ignore_checksum_file=True)
        os.replace(working, output)
        return exit_status
    except BaseException as error:
        for child in list(working.iterdir()):
            shutil.rmtree(child) if child.is_dir() else child.unlink()
        publish_status(
            working, status="invalid", reasons=[f"{type(error).__name__}:{error}"],
            custody_sha=custody_sha,
        )
        write_sha256sums(working)
        verify_sha256sums(working, ignore_checksum_file=True)
        os.replace(working, output)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
