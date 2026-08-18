#!/usr/bin/env python3
from __future__ import annotations

import json
import os
from pathlib import Path, PurePosixPath
import shutil
import sys
import tempfile
from typing import Any

from ploidypatch.artifact_manifest import (
    sha256_file,
    verify_sha256sums,
    write_sha256sums,
)
from ploidypatch.holdout_contract import load_holdout_contract
from ploidypatch.walnut_h1 import (
    EVALUATION_SCHEMA,
    FORMAL_STATUSES,
    HOLDOUT_ID,
    POLICY_ID,
    evaluate_h1_scores,
    load_json_object,
)


STATUS_SCHEMA = "ploidypatch.walnut_h1_evaluation_status.v0.8"
REVEAL_SCHEMA = "ploidypatch.walnut_h1_reveal_inputs.v0.8"


def required_environment(name: str) -> Path:
    value = os.environ.get(name)
    if not value:
        raise ValueError(f"{name} is required")
    return Path(value).resolve()


def safe_relative(value: str, field: str) -> PurePosixPath:
    path = PurePosixPath(value)
    if path.is_absolute() or not path.parts or any(part in {"", ".", ".."} for part in path.parts):
        raise ValueError(f"Unsafe relative path for {field}: {value!r}")
    return path


def resolve_bound_input(root: Path, value: Any, field: str) -> Path:
    if not isinstance(value, dict) or set(value) != {"relative_path", "sha256"}:
        raise ValueError(f"Malformed reveal input binding: {field}")
    relative = safe_relative(value["relative_path"], field)
    path = root.joinpath(*relative.parts)
    if (
        not path.is_file()
        or path.is_symlink()
        or sha256_file(path) != value["sha256"]
    ):
        raise ValueError(f"Reveal input differs from manifest: {field}")
    return path


def write_record(path: Path, value: dict[str, Any]) -> None:
    path.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def main(argv: list[str] | None = None) -> int:
    values = sys.argv[1:] if argv is None else argv
    if len(values) != 1:
        raise SystemExit("usage: evaluate_walnut_external_h1_v0.8.py PROJECT_ROOT")
    project_root = Path(values[0]).resolve()
    if not project_root.is_dir():
        raise SystemExit("PROJECT_ROOT must be an existing directory")

    reveal_root = required_environment("PLOIDYPATCH_REVEAL_INPUTS_ROOT")
    output = required_environment("PLOIDYPATCH_EVALUATION_OUTPUT")
    contract_path = required_environment("PLOIDYPATCH_HOLDOUT_CONTRACT")
    custody_path = required_environment("PLOIDYPATCH_CUSTODY_MANIFEST")
    authorization_path = required_environment("PLOIDYPATCH_REVEAL_AUTHORIZATION")
    for name in (
        "PLOIDYPATCH_PROTOCOL_FREEZE",
        "PLOIDYPATCH_EXECUTION_FREEZE",
        "PLOIDYPATCH_BLIND_RUN_ROOT",
        "PLOIDYPATCH_EVALUATOR_ONLY_ROOT",
    ):
        required_environment(name)

    if output.exists() or output.is_symlink():
        raise FileExistsError(f"Refusing to overwrite Walnut evaluation: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    working = Path(tempfile.mkdtemp(prefix=f".{output.name}.working.", dir=output.parent))
    exit_status = 0
    try:
        verify_sha256sums(reveal_root, ignore_checksum_file=True)
        contract = load_holdout_contract(contract_path)
        if contract.holdout_id != HOLDOUT_ID or contract.policy_id != POLICY_ID:
            raise ValueError("Walnut evaluator received a different holdout contract")
        custody = load_json_object(custody_path)
        authorization = load_json_object(authorization_path)
        if (
            authorization.get("truth_reveal_authorized") is not True
            or authorization.get("custody_manifest_sha256") != sha256_file(custody_path)
            or custody.get("holdout_id") != HOLDOUT_ID
        ):
            raise ValueError("Walnut reveal authorization or custody binding differs")
        reveal = load_json_object(reveal_root / "reveal_input_manifest.json")
        status_record = load_json_object(reveal_root / "status.json")
        formal_status = reveal.get("formal_status")
        reason_codes = status_record.get("reason_codes")
        if (
            reveal.get("schema_version") != REVEAL_SCHEMA
            or reveal.get("holdout_id") != HOLDOUT_ID
            or reveal.get("policy_id") != POLICY_ID
            or formal_status not in FORMAL_STATUSES
            or status_record.get("status") != formal_status
            or not isinstance(reason_codes, list)
            or (formal_status == "ready" and reason_codes)
            or (formal_status != "ready" and not reason_codes)
            or reveal.get("custody_manifest_sha256") != sha256_file(custody_path)
        ):
            raise ValueError("Walnut reveal-input manifest differs")

        if formal_status == "ready":
            bindings = reveal.get("evaluation_inputs")
            if not isinstance(bindings, dict) or set(bindings) != {
                "retain_score", "suppress_score", "evaluability"
            }:
                raise ValueError("Walnut ready reveal lacks the exact H1 inputs")
            retain = resolve_bound_input(reveal_root, bindings["retain_score"], "retain_score")
            suppress = resolve_bound_input(reveal_root, bindings["suppress_score"], "suppress_score")
            resolve_bound_input(reveal_root, bindings["evaluability"], "evaluability")
            evaluation = evaluate_h1_scores(
                retain_score_path=retain,
                suppress_score_path=suppress,
                output=working / "evaluation.json",
            )
            if evaluation.get("schema_version") != EVALUATION_SCHEMA:
                raise ValueError("Walnut H1 evaluation schema differs")
        else:
            evaluation = {
                "schema_version": EVALUATION_SCHEMA,
                "holdout_id": HOLDOUT_ID,
                "policy_id": POLICY_ID,
                "status": formal_status,
                "hypothesis": "H1_retain_distinct_vs_suppress_overlap_only",
                "reason_codes": reason_codes,
                "formal_test_executed": False,
                "automatic_approval": False,
                "ranker_or_model_executed": False,
                "h2_or_topology_ranking_executed": False,
            }
            write_record(working / "evaluation.json", evaluation)
        write_record(
            working / "status.json",
            {
                "schema_version": STATUS_SCHEMA,
                "holdout_id": HOLDOUT_ID,
                "policy_id": POLICY_ID,
                "status": formal_status,
                "reason_codes": reason_codes,
            },
        )
    except BaseException as error:
        exit_status = 1
        for child in list(working.iterdir()):
            if child.is_dir():
                shutil.rmtree(child)
            else:
                child.unlink()
        formal_status = "invalid"
        reason = f"{type(error).__name__}:{error}"
        write_record(
            working / "evaluation.json",
            {
                "schema_version": EVALUATION_SCHEMA,
                "holdout_id": HOLDOUT_ID,
                "policy_id": POLICY_ID,
                "status": formal_status,
                "hypothesis": "H1_retain_distinct_vs_suppress_overlap_only",
                "reason_codes": [reason],
                "formal_test_executed": False,
                "automatic_approval": False,
                "ranker_or_model_executed": False,
                "h2_or_topology_ranking_executed": False,
            },
        )
        write_record(
            working / "status.json",
            {
                "schema_version": STATUS_SCHEMA,
                "holdout_id": HOLDOUT_ID,
                "policy_id": POLICY_ID,
                "status": formal_status,
                "reason_codes": [reason],
            },
        )
    write_sha256sums(working)
    verify_sha256sums(working, ignore_checksum_file=True)
    os.replace(working, output)
    return exit_status


if __name__ == "__main__":
    raise SystemExit(main())
