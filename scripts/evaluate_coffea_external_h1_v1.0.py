#!/usr/bin/env python3
"""Evaluate the frozen Coffea combined-reference core H1 after reveal."""
from __future__ import annotations

import json
import os
from pathlib import Path, PurePosixPath
import shutil
import sys
import tempfile
from typing import Any

from ploidypatch.artifact_manifest import sha256_file, verify_sha256sums, write_sha256sums
from ploidypatch.coffea_h1_framework import (
    CUSTODY_SCHEMA,
    EVALUATION_SCHEMA,
    EVALUATION_STATUS_SCHEMA,
    HOLDOUT_ID,
    POLICY_ID,
    REVEAL_INPUT_SCHEMA,
    load_json,
)
from ploidypatch.core_h1_evaluation import EXPECTED_SCORE_KEYS, evaluate_core_h1_scores
from ploidypatch.holdout_contract import load_holdout_contract


def required(name: str) -> Path:
    value = os.environ.get(name)
    if not value:
        raise ValueError(f"{name} is required")
    return Path(value).resolve()


def bound_input(root: Path, value: Any, field: str) -> Path:
    if not isinstance(value, dict) or set(value) != {"relative_path", "sha256"}:
        raise ValueError(f"Malformed Coffea reveal binding: {field}")
    relative = PurePosixPath(value["relative_path"])
    if relative.is_absolute() or any(part in {"", ".", ".."} for part in relative.parts):
        raise ValueError(f"Unsafe Coffea reveal path: {field}")
    path = root.joinpath(*relative.parts)
    if not path.is_file() or path.is_symlink() or sha256_file(path) != value["sha256"]:
        raise ValueError(f"Coffea reveal input differs: {field}")
    return path


def write_json(path: Path, value: dict[str, Any]) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def main(argv: list[str] | None = None) -> int:
    values = sys.argv[1:] if argv is None else argv
    if len(values) != 1:
        raise SystemExit("usage: evaluate_coffea_external_h1_v1.0.py PROJECT_ROOT")
    project_root = Path(values[0]).resolve()
    reveal_root = required("PLOIDYPATCH_REVEAL_INPUTS_ROOT")
    output = required("PLOIDYPATCH_EVALUATION_OUTPUT")
    contract_path = required("PLOIDYPATCH_HOLDOUT_CONTRACT")
    custody_path = required("PLOIDYPATCH_CUSTODY_MANIFEST")
    authorization_path = required("PLOIDYPATCH_REVEAL_AUTHORIZATION")
    for name in (
        "PLOIDYPATCH_PROTOCOL_FREEZE",
        "PLOIDYPATCH_EXECUTION_FREEZE",
        "PLOIDYPATCH_BLIND_RUN_ROOT",
        "PLOIDYPATCH_EVALUATOR_INPUT_ROOT",
        "PLOIDYPATCH_EVALUATOR_ONLY_ROOT",
    ):
        required(name)
    if output.exists() or output.is_symlink():
        raise FileExistsError(f"Refusing to overwrite Coffea evaluation: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    working = Path(tempfile.mkdtemp(prefix=f".{output.name}.working.", dir=output.parent))
    exit_status = 0
    try:
        verify_sha256sums(reveal_root, ignore_checksum_file=True)
        contract = load_holdout_contract(contract_path)
        if contract.holdout_id != HOLDOUT_ID or contract.policy_id != POLICY_ID:
            raise ValueError("Coffea evaluator received a different contract")
        custody = load_json(custody_path)
        authorization = load_json(authorization_path)
        if (
            custody.get("schema_version") != CUSTODY_SCHEMA
            or custody.get("holdout_id") != HOLDOUT_ID
            or authorization.get("truth_reveal_authorized") is not True
            or authorization.get("custody_manifest_sha256") != sha256_file(custody_path)
            or authorization.get("ranker_or_model_authorized") is not False
        ):
            raise ValueError("Coffea authorization or custody binding differs")
        reveal = load_json(reveal_root / "reveal_input_manifest.json")
        status_record = load_json(reveal_root / "status.json")
        status = reveal.get("formal_status")
        reasons = status_record.get("reason_codes")
        if (
            reveal.get("schema_version") != REVEAL_INPUT_SCHEMA
            or reveal.get("holdout_id") != HOLDOUT_ID
            or reveal.get("policy_id") != POLICY_ID
            or status not in {"ready", "not_evaluable", "invalid"}
            or status_record.get("status") != status
            or not isinstance(reasons, list)
            or (status == "ready" and reasons)
            or (status != "ready" and not reasons)
            or reveal.get("custody_manifest_sha256") != sha256_file(custody_path)
        ):
            raise ValueError("Coffea reveal input manifest differs")
        if status == "ready":
            bindings = reveal.get("evaluation_inputs")
            expected = {*EXPECTED_SCORE_KEYS, "evaluability"}
            if not isinstance(bindings, dict) or set(bindings) != expected:
                raise ValueError("Coffea ready reveal lacks exact six scores")
            score_paths = {
                name: bound_input(reveal_root, bindings[name], name)
                for name in EXPECTED_SCORE_KEYS
            }
            bound_input(reveal_root, bindings["evaluability"], "evaluability")
            evaluation = evaluate_core_h1_scores(
                score_paths=score_paths,
                holdout_id=HOLDOUT_ID,
                policy_id=POLICY_ID,
                schema_version=EVALUATION_SCHEMA,
                primary_scope="combined",
                bootstrap_seed=20260912,
                bootstrap_replicates=20_000,
                output=working / "evaluation.json",
            )
            if evaluation.get("schema_version") != EVALUATION_SCHEMA:
                raise ValueError("Coffea H1 evaluation schema differs")
        else:
            evaluation = {
                "schema_version": EVALUATION_SCHEMA,
                "holdout_id": HOLDOUT_ID,
                "policy_id": POLICY_ID,
                "status": status,
                "formal_outcome": status,
                "hypothesis": "H1_retain_distinct_vs_suppress_overlap_only",
                "reason_codes": reasons,
                "formal_test_executed": False,
                "automatic_approval": False,
                "ranker_or_model_executed": False,
                "h2_or_topology_ranking_executed": False,
            }
            write_json(working / "evaluation.json", evaluation)
        write_json(
            working / "status.json",
            {
                "schema_version": EVALUATION_STATUS_SCHEMA,
                "holdout_id": HOLDOUT_ID,
                "policy_id": POLICY_ID,
                "status": status,
                "reason_codes": reasons,
            },
        )
    except BaseException as error:
        exit_status = 1
        for child in list(working.iterdir()):
            shutil.rmtree(child) if child.is_dir() else child.unlink()
        reason = f"{type(error).__name__}:{error}"
        write_json(
            working / "evaluation.json",
            {
                "schema_version": EVALUATION_SCHEMA,
                "holdout_id": HOLDOUT_ID,
                "policy_id": POLICY_ID,
                "status": "invalid",
                "formal_outcome": "invalid",
                "reason_codes": [reason],
                "formal_test_executed": False,
                "automatic_approval": False,
                "ranker_or_model_executed": False,
                "h2_or_topology_ranking_executed": False,
            },
        )
        write_json(
            working / "status.json",
            {
                "schema_version": EVALUATION_STATUS_SCHEMA,
                "holdout_id": HOLDOUT_ID,
                "policy_id": POLICY_ID,
                "status": "invalid",
                "reason_codes": [reason],
            },
        )
    write_sha256sums(working)
    verify_sha256sums(working, ignore_checksum_file=True)
    os.replace(working, output)
    return exit_status


if __name__ == "__main__":
    raise SystemExit(main())
