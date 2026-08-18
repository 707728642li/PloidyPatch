#!/usr/bin/env python3
"""Authorize and execute the Coffea H1 reveal after blind custody."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
from typing import Any, Iterable

from ploidypatch.artifact_manifest import sha256_file, verify_sha256sums, write_sha256sums
from ploidypatch.coffea_h1_framework import (
    BLIND_OUTPUTS,
    CUSTODY_SCHEMA,
    EVALUATION_STATUS_SCHEMA,
    HOLDOUT_ID,
    POLICY_ID,
    REVEAL_STATUS_SCHEMA,
    load_json,
    validate_status,
    verify_execution,
)


AUTHORIZATION_SCHEMA = "ploidypatch.coffea_h1_reveal_authorization.v1.0"
RESULT_SCHEMA = "ploidypatch.coffea_h1_reveal_result.v1.0"


def write_json(path: Path, value: dict[str, Any]) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def verify_custody_before_truth(
    *, blind_run: Path, execution: Path, protocol: Path
) -> Path:
    execution_manifest, _, _ = verify_execution(execution, protocol)
    verify_sha256sums(blind_run, ignore_checksum_file=True)
    custody_path = blind_run / "custody_manifest.json"
    custody = load_json(custody_path)
    expected_execution_sha256 = sha256_file(execution / "SHA256SUMS")
    patch = execution_manifest.get("execution_patch")
    if isinstance(patch, dict) and patch.get("patch_sequence") == 5:
        expected_execution_sha256 = patch.get(
            "blind_custody_execution_SHA256SUMS_sha256"
        )
        if (
            not isinstance(expected_execution_sha256, str)
            or len(expected_execution_sha256) != 64
            or patch.get("blind_custody_manifest_sha256")
            != sha256_file(custody_path)
        ):
            raise ValueError("Coffea custody patch binding differs")
    if (
        custody.get("schema_version") != CUSTODY_SCHEMA
        or custody.get("holdout_id") != HOLDOUT_ID
        or custody.get("policy_id") != POLICY_ID
        or custody.get("protocol_SHA256SUMS_sha256")
        != sha256_file(protocol / "SHA256SUMS")
        or custody.get("execution_SHA256SUMS_sha256")
        != expected_execution_sha256
        or custody.get("ranker_or_model_executed") is not False
        or custody.get("h2_or_topology_ranking_executed") is not False
        or any(custody.get(field) is not False for field in (
            "truth_mounted", "complete_target_annotation_mounted",
            "evaluator_references_mounted", "nas_data_mounted", "network_access",
        ))
    ):
        raise ValueError("Coffea custody failed pre-truth authorization")
    bindings = custody.get("blind_outputs")
    if not isinstance(bindings, dict) or set(bindings) != set(BLIND_OUTPUTS):
        raise ValueError("Coffea custody lacks exact blind output universe")
    project = blind_run / "project"
    for name, relative in BLIND_OUTPUTS.items():
        item = bindings[name]
        path = project / Path(relative)
        if (
            not isinstance(item, dict)
            or item.get("relative_path") != relative
            or item.get("sha256") != sha256_file(path)
            or item.get("bytes") != path.stat().st_size
        ):
            raise ValueError(f"Coffea blind output differs before reveal: {name}")
    return custody_path


def dev_environment(execution: Path) -> Path:
    manifest = load_json(execution / "execution_manifest.json")
    rows = manifest.get("environments")
    if not isinstance(rows, list):
        raise ValueError("Coffea execution lacks environment bindings")
    matches = [row for row in rows if row.get("name") == "ploidypatch-dev"]
    if len(matches) != 1:
        raise ValueError("Coffea execution lacks one dev environment")
    row = matches[0]
    prefix = Path(row["host_prefix"]).resolve()
    python = prefix / "bin/python"
    try:
        resolved_python = python.resolve(strict=True)
        resolved_python.relative_to(prefix)
    except (FileNotFoundError, ValueError):
        raise ValueError("Coffea frozen dev Python escapes its environment") from None
    if not resolved_python.is_file() or not os.access(resolved_python, os.X_OK):
        raise ValueError("Coffea frozen dev Python is missing")
    explicit = subprocess.run(
        ["conda", "list", "--explicit", "-p", str(prefix)],
        check=True, capture_output=True,
    ).stdout
    if sha256_file(execution / row["explicit_relative_path"]) != row["explicit_sha256"]:
        raise ValueError("Coffea frozen explicit lock differs")
    import hashlib
    if hashlib.sha256(explicit).hexdigest() != row["explicit_sha256"]:
        raise ValueError("Coffea current dev environment differs from freeze")
    return resolved_python


def run_entry(
    *, python: Path, source: Path, entry: str, project_root: Path,
    environment: dict[str, str], stdout: Path, stderr: Path
) -> int:
    env = {"PATH": os.environ.get("PATH", ""), "LANG": "C.UTF-8", "LC_ALL": "C.UTF-8"}
    env.update(environment)
    env["PYTHONPATH"] = str(source / "src")
    with stdout.open("x", encoding="utf-8") as out, stderr.open("x", encoding="utf-8") as err:
        return subprocess.run(
            [str(python), str(source / entry), str(project_root)],
            env=env, stdout=out, stderr=err, check=False,
        ).returncode


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", required=True)
    parser.add_argument("--execution-freeze", required=True)
    parser.add_argument("--protocol-freeze", required=True)
    parser.add_argument("--blind-run", required=True)
    # Deliberately retained as an opaque string until the custody barrier passes.
    parser.add_argument("--evaluator-input-root", required=True)
    parser.add_argument("--evaluator-run-root", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args(argv)
    project_root = Path(args.project_root).resolve()
    execution = Path(args.execution_freeze).resolve()
    protocol = Path(args.protocol_freeze).resolve()
    blind_run = Path(args.blind_run).resolve()
    output = Path(args.output_dir).resolve()
    allowed = project_root / "results/evaluator"
    try:
        output.relative_to(allowed.resolve())
    except ValueError:
        raise ValueError("Coffea reveal output must be under results/evaluator") from None
    if output.exists() or output.is_symlink():
        raise FileExistsError(f"Refusing to overwrite Coffea reveal: {output}")
    # Hard barrier: no stat/resolve/read of evaluator-only argument above this line.
    custody_path = verify_custody_before_truth(
        blind_run=blind_run, execution=execution, protocol=protocol
    )
    custody_sha = sha256_file(custody_path)
    evaluator_input = Path(args.evaluator_input_root).resolve()
    evaluator = Path(args.evaluator_run_root).resolve()
    verify_sha256sums(evaluator_input, ignore_checksum_file=True)
    verify_sha256sums(evaluator, ignore_checksum_file=True)
    output.parent.mkdir(parents=True, exist_ok=True)
    working = Path(tempfile.mkdtemp(prefix=f".{output.name}.working.", dir=output.parent))
    exit_status = 0
    try:
        authorization = working / "reveal_authorization.json"
        write_json(
            authorization,
            {
                "schema_version": AUTHORIZATION_SCHEMA,
                "holdout_id": HOLDOUT_ID,
                "policy_id": POLICY_ID,
                "custody_manifest_sha256": custody_sha,
                "truth_reveal_authorized": True,
                "ranker_or_model_authorized": False,
                "h2_or_topology_ranking_authorized": False,
            },
        )
        manifest, _, _ = verify_execution(execution, protocol)
        python = dev_environment(execution)
        source = execution / "source"
        common = {
            "PLOIDYPATCH_HOLDOUT_CONTRACT": str(protocol / "contract.json"),
            "PLOIDYPATCH_PROTOCOL_FREEZE": str(protocol),
            "PLOIDYPATCH_EXECUTION_FREEZE": str(execution),
            "PLOIDYPATCH_BLIND_RUN_ROOT": str(blind_run),
            "PLOIDYPATCH_CUSTODY_MANIFEST": str(custody_path),
            "PLOIDYPATCH_REVEAL_AUTHORIZATION": str(authorization),
            "PLOIDYPATCH_EVALUATOR_INPUT_ROOT": str(evaluator_input),
            "PLOIDYPATCH_EVALUATOR_ONLY_ROOT": str(evaluator),
        }
        reveal_inputs = working / "reveal_inputs"
        builder_env = {
            **common,
            "PLOIDYPATCH_REVEAL_INPUTS_OUTPUT": str(reveal_inputs),
        }
        builder_status = run_entry(
            python=python, source=source,
            entry=manifest["pipeline_entries"]["reveal_input_builder"],
            project_root=project_root, environment=builder_env,
            stdout=working / "builder.stdout.log", stderr=working / "builder.stderr.log",
        )
        if not reveal_inputs.is_dir():
            raise RuntimeError(f"Coffea reveal builder produced no result ({builder_status})")
        evaluation = working / "evaluation"
        evaluator_env = {
            **common,
            "PLOIDYPATCH_REVEAL_INPUTS_ROOT": str(reveal_inputs),
            "PLOIDYPATCH_EVALUATION_OUTPUT": str(evaluation),
        }
        evaluator_status = run_entry(
            python=python, source=source,
            entry=manifest["pipeline_entries"]["evaluator"],
            project_root=project_root, environment=evaluator_env,
            stdout=working / "evaluator.stdout.log", stderr=working / "evaluator.stderr.log",
        )
        if not evaluation.is_dir():
            raise RuntimeError(f"Coffea evaluator produced no result ({evaluator_status})")
        reveal_status = load_json(reveal_inputs / "status.json")
        evaluation_status = load_json(evaluation / "status.json")
        status = validate_status(reveal_status, expected_schema=REVEAL_STATUS_SCHEMA)
        evaluated = validate_status(
            evaluation_status, expected_schema=EVALUATION_STATUS_SCHEMA
        )
        if status != evaluated:
            raise ValueError("Coffea reveal/evaluation tri-state differs")
        exit_status = 1 if status == "invalid" else 0
        write_json(
            working / "result_manifest.json",
            {
                "schema_version": RESULT_SCHEMA,
                "holdout_id": HOLDOUT_ID,
                "policy_id": POLICY_ID,
                "status": status,
                "reason_codes": evaluation_status["reason_codes"],
                "custody_manifest_sha256": custody_sha,
                "execution_SHA256SUMS_sha256": sha256_file(execution / "SHA256SUMS"),
                "protocol_SHA256SUMS_sha256": sha256_file(protocol / "SHA256SUMS"),
                "evaluator_only_SHA256SUMS_sha256": sha256_file(
                    evaluator / "SHA256SUMS"
                ),
                "evaluator_input_SHA256SUMS_sha256": sha256_file(
                    evaluator_input / "SHA256SUMS"
                ),
                "builder_exit_status": builder_status,
                "evaluator_exit_status": evaluator_status,
                "ranker_or_model_executed": False,
                "h2_or_topology_ranking_executed": False,
            },
        )
        write_sha256sums(working)
        verify_sha256sums(working, ignore_checksum_file=True)
        if os.name != "nt":
            for path in sorted(working.rglob("*"), reverse=True):
                path.chmod(0o550 if path.is_dir() else 0o440)
            working.chmod(0o550)
        os.replace(working, output)
        return exit_status
    except BaseException:
        # Retain the non-overwritable working tree and logs for audit.
        failed = output.with_name(output.name + ".invalid_run")
        if failed.exists() or failed.is_symlink():
            raise RuntimeError(f"Coffea invalid reveal output exists: {failed}")
        os.replace(working, failed)
        raise


if __name__ == "__main__":
    raise SystemExit(main())
