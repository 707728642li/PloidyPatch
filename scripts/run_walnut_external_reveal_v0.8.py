#!/usr/bin/env python3
"""Authorize and orchestrate Walnut core-H1 reveal with a sealed tri-state result."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
from typing import Any, Iterable

from ploidypatch.artifact_manifest import sha256_file, verify_sha256sums, write_sha256sums
from ploidypatch.walnut_h1_framework import (
    CUSTODY_SCHEMA,
    EVALUATION_SCHEMA,
    EVALUATION_STATUS_SCHEMA,
    PIPELINE_ENTRIES,
    REVEAL_STATUS_SCHEMA,
    atomic_write_json,
    load_json,
    validate_status,
    verify_execution,
)


ORCHESTRATION_SCHEMA = "ploidypatch.walnut_h1_reveal_orchestration.v0.8"
REVEAL_INPUT_SCHEMA = "ploidypatch.walnut_h1_reveal_inputs.v0.8"
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


def validate_reason_class(status: str, reasons: list[str]) -> None:
    if status == "not_evaluable" and not set(reasons) <= NOT_EVALUABLE_REASONS:
        raise ValueError("not_evaluable contains a non-data-gate reason")


def validate_custody(
    blind_run: Path, execution: Path, protocol: Path
) -> dict[str, Any]:
    verify_sha256sums(blind_run, ignore_checksum_file=True)
    custody = load_json(blind_run / "custody_manifest.json")
    if (
        custody.get("schema_version") != CUSTODY_SCHEMA
        or custody.get("protocol_SHA256SUMS_sha256")
        != sha256_file(protocol / "SHA256SUMS")
        or custody.get("execution_SHA256SUMS_sha256")
        != sha256_file(execution / "SHA256SUMS")
        or custody.get("ranker_or_model_executed") is not False
        or custody.get("h2_or_topology_ranking_executed") is not False
        or any(custody.get(key) is not False for key in (
            "truth_mounted", "complete_target_annotation_mounted",
            "evaluator_references_mounted", "nas_data_mounted", "network_access",
        ))
    ):
        raise ValueError("Blind custody does not authorize Walnut reveal")
    return custody


def run_frozen_entry(
    *, python: Path, entry: Path, project_root: Path, environment: dict[str, str],
    stdout: Path, stderr: Path
) -> int:
    if not entry.is_file() or entry.is_symlink() or not python.is_file():
        raise ValueError("Frozen reveal entry or Python is missing")
    with stdout.open("x", encoding="utf-8") as out, stderr.open(
        "x", encoding="utf-8"
    ) as err:
        return subprocess.run(
            [str(python), str(entry), str(project_root)], env=environment,
            stdout=out, stderr=err, check=False,
        ).returncode


def sealed_status_or_invalid(root: Path, expected_schema: str, reason: str) -> dict[str, Any]:
    try:
        verify_sha256sums(root, ignore_checksum_file=True)
        value = load_json(root / "status.json")
        validate_status(value, expected_schema=expected_schema)
        return value
    except (OSError, ValueError, json.JSONDecodeError):
        return {"schema_version": expected_schema, "status": "invalid",
                "reason_codes": [reason]}


def forbidden_evaluation_keys(value: Any) -> list[str]:
    forbidden = ("h2", "average_precision", "topology", "rank", "score")
    allowed_firewall_keys = {
        "ranker_or_model_executed", "h2_or_topology_ranking_executed"
    }
    hits: list[str] = []
    if isinstance(value, dict):
        for key, item in value.items():
            if key not in allowed_firewall_keys and any(
                token in key.casefold() for token in forbidden
            ):
                hits.append(key)
            hits.extend(forbidden_evaluation_keys(item))
    elif isinstance(value, list):
        for item in value:
            hits.extend(forbidden_evaluation_keys(item))
    return hits


def validate_evaluation(root: Path) -> dict[str, Any]:
    verify_sha256sums(root, ignore_checksum_file=True)
    status_value = load_json(root / "status.json")
    status = validate_status(status_value, expected_schema=EVALUATION_STATUS_SCHEMA)
    validate_reason_class(status, status_value["reason_codes"])
    evaluation = load_json(root / "evaluation.json")
    if (
        evaluation.get("schema_version") != EVALUATION_SCHEMA
        or evaluation.get("status") != status
        or evaluation.get("reason_codes") != status_value["reason_codes"]
        or evaluation.get("ranker_or_model_executed") is not False
        or evaluation.get("h2_or_topology_ranking_executed") is not False
        or forbidden_evaluation_keys(evaluation)
    ):
        raise ValueError("Evaluation contains forbidden ranker/H2 semantics")
    if status == "ready" and (
        evaluation.get("bootstrap_replicates") != 20_000
        or evaluation.get("bootstrap_unit") != "paired_event"
        or evaluation.get("metric")
        != "event_exact_phased_CDS_recall_retain_distinct_minus_suppress_overlap"
        or evaluation.get("all_arm_collateral_loss") != 0
    ):
        raise ValueError("Ready H1 evaluation differs from fixed estimand/safety gate")
    return status_value


def orchestrate_reveal(
    *, project_root: Path, protocol: Path, execution: Path, blind_run: Path,
    evaluator_only: Path, output: Path
) -> Path:
    working = Path(str(output) + ".working")
    if output.exists() or output.is_symlink() or working.exists() or working.is_symlink():
        raise FileExistsError("Refusing to overwrite reveal output or working directory")
    execution_manifest, protocol_manifest, contract = verify_execution(execution, protocol)
    custody = validate_custody(blind_run, execution, protocol)
    working.mkdir(parents=True)
    try:
        authorization = working / "reveal_authorization.json"
        atomic_write_json(
            authorization,
            {"schema_version": "ploidypatch.walnut_h1_reveal_authorization.v0.8",
             "holdout_id": contract.holdout_id,
             "protocol_SHA256SUMS_sha256": sha256_file(protocol / "SHA256SUMS"),
             "execution_SHA256SUMS_sha256": sha256_file(execution / "SHA256SUMS"),
             "blind_run_SHA256SUMS_sha256": sha256_file(blind_run / "SHA256SUMS"),
             "custody_manifest_sha256": sha256_file(blind_run / "custody_manifest.json"),
             "truth_reveal_authorized": True,
             "ranker_or_model_authorized": False,
             "h2_or_topology_ranking_authorized": False},
        )
        # Evaluator-owned bytes are resolved only after immutable authorization.
        if not evaluator_only.is_dir() or evaluator_only.is_symlink():
            raise ValueError("Evaluator-only root is missing or symlinked")
        verify_sha256sums(evaluator_only, ignore_checksum_file=True)
        dev = next(
            Path(row["host_prefix"]) for row in execution_manifest["environments"]
            if row["name"] == "ploidypatch-dev"
        )
        python = dev / ("python.exe" if os.name == "nt" else "bin/python")
        environment = os.environ.copy()
        for key in list(environment):
            if "MODEL" in key or "RANK" in key:
                del environment[key]
        environment.update(
            {"PLOIDYPATCH_NETWORK_ACCESS": "none",
             "PLOIDYPATCH_HOLDOUT_CONTRACT": str(protocol / "contract.json"),
             "PLOIDYPATCH_PROTOCOL_FREEZE": str(protocol),
             "PLOIDYPATCH_EXECUTION_FREEZE": str(execution),
             "PLOIDYPATCH_BLIND_RUN_ROOT": str(blind_run),
             "PLOIDYPATCH_CUSTODY_MANIFEST": str(blind_run / "custody_manifest.json"),
             "PLOIDYPATCH_REVEAL_AUTHORIZATION": str(authorization),
             "PLOIDYPATCH_EVALUATOR_ONLY_ROOT": str(evaluator_only),
             "PYTHONPATH": str(execution / "source/src")},
        )
        reveal_inputs = working / "reveal_inputs"
        environment["PLOIDYPATCH_REVEAL_INPUTS_OUTPUT"] = str(reveal_inputs)
        builder = execution / "source" / PIPELINE_ENTRIES["reveal_input_builder"]
        builder_status = run_frozen_entry(
            python=python, entry=builder, project_root=project_root,
            environment=environment, stdout=working / "builder.stdout.log",
            stderr=working / "builder.stderr.log",
        )
        reveal_status = sealed_status_or_invalid(
            reveal_inputs, REVEAL_STATUS_SCHEMA, "reveal_builder_unsealed_failure"
        )
        status = validate_status(reveal_status, expected_schema=REVEAL_STATUS_SCHEMA)
        validate_reason_class(status, reveal_status["reason_codes"])
        if builder_status != 0 and status != "invalid":
            status = "invalid"; reveal_status = {
                "schema_version": REVEAL_STATUS_SCHEMA, "status": status,
                "reason_codes": ["reveal_builder_nonzero_with_noninvalid_status"]}
        if status == "ready":
            manifest = load_json(reveal_inputs / "reveal_input_manifest.json")
            if (
                manifest.get("schema_version") != REVEAL_INPUT_SCHEMA
                or manifest.get("formal_status") != "ready"
                or manifest.get("custody_manifest_sha256")
                != sha256_file(blind_run / "custody_manifest.json")
                or manifest.get("ranker_or_model_access") is not False
                or manifest.get("h2_or_topology_ranking_access") is not False
            ):
                status = "invalid"; reveal_status = {
                    "schema_version": REVEAL_STATUS_SCHEMA, "status": status,
                    "reason_codes": ["reveal_input_manifest_lineage_violation"]}
        evaluation_status: dict[str, Any] | None = None
        if status == "ready":
            evaluation_root = working / "evaluation"
            environment["PLOIDYPATCH_REVEAL_INPUTS_ROOT"] = str(reveal_inputs)
            environment["PLOIDYPATCH_EVALUATION_OUTPUT"] = str(evaluation_root)
            evaluator = execution / "source" / PIPELINE_ENTRIES["evaluator"]
            evaluator_rc = run_frozen_entry(
                python=python, entry=evaluator, project_root=project_root,
                environment=environment, stdout=working / "evaluator.stdout.log",
                stderr=working / "evaluator.stderr.log",
            )
            try:
                evaluation_status = validate_evaluation(evaluation_root)
                status = evaluation_status["status"]
                if evaluator_rc != 0 and status != "invalid":
                    raise ValueError("Evaluator returned nonzero with noninvalid status")
            except (OSError, ValueError, json.JSONDecodeError):
                status = "invalid"; evaluation_status = {
                    "schema_version": EVALUATION_STATUS_SCHEMA, "status": "invalid",
                    "reason_codes": ["evaluation_unsealed_or_invalid"]}
        reasons = (
            evaluation_status["reason_codes"] if evaluation_status is not None
            else reveal_status["reason_codes"]
        )
        final_status = {"schema_version": ORCHESTRATION_SCHEMA, "status": status,
                        "reason_codes": reasons,
                        "evaluator_invoked": evaluation_status is not None}
        atomic_write_json(working / "status.json", final_status)
        atomic_write_json(
            working / "orchestration_manifest.json",
            {**final_status, "holdout_id": contract.holdout_id,
             "protocol_SHA256SUMS_sha256": sha256_file(protocol / "SHA256SUMS"),
             "execution_SHA256SUMS_sha256": sha256_file(execution / "SHA256SUMS"),
             "blind_run_SHA256SUMS_sha256": sha256_file(blind_run / "SHA256SUMS"),
             "evaluator_only_SHA256SUMS_sha256": sha256_file(evaluator_only / "SHA256SUMS"),
             "authorization_sha256": sha256_file(authorization)},
        )
        write_sha256sums(working)
        verify_sha256sums(working, ignore_checksum_file=True)
        os.replace(working, output)
    except BaseException:
        if working.exists() and not working.is_symlink():
            shutil.rmtree(working)
        raise
    return output


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", required=True)
    parser.add_argument("--protocol-freeze", required=True)
    parser.add_argument("--execution-freeze", required=True)
    parser.add_argument("--blind-run-root", required=True)
    parser.add_argument("--evaluator-only-root", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args(argv)
    orchestrate_reveal(
        project_root=Path(args.project_root).resolve(),
        protocol=Path(args.protocol_freeze).resolve(),
        execution=Path(args.execution_freeze).resolve(),
        blind_run=Path(args.blind_run_root).resolve(),
        evaluator_only=Path(args.evaluator_only_root),
        output=Path(args.output_dir).resolve(),
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
