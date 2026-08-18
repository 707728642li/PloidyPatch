#!/usr/bin/env python3
"""Resolve frozen Actinidia reveal inputs and invoke the generic v0.5 evaluator."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess
import sys
from typing import Any, Iterable

from ploidypatch.artifact_manifest import sha256_file, verify_sha256sums
from ploidypatch.holdout_contract import safe_relative_path


PROTOCOL_SCHEMA = "ploidypatch.external_holdout_protocol_freeze.v0.5"
EXECUTION_SCHEMA = "ploidypatch.external_holdout_execution_freeze.v0.5"
CUSTODY_SCHEMA = "ploidypatch.external_holdout_blind_custody.v0.5"
REVEAL_SCHEMA = "ploidypatch.actinidia_reveal_inputs.v0.5"
POLICY_ID = "ploidypatch_actinidia_external_validation_v0.5"
HOLDOUT_ID = "actinidia_red5_v0.5"


def _load_json(path: Path) -> dict[str, Any]:
    if not path.is_file() or path.is_symlink() or path.stat().st_size == 0:
        raise ValueError(f"Missing, empty or symlinked JSON input: {path}")
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return value


def _safe_join(root: Path, value: str, context: str) -> Path:
    relative = safe_relative_path(value, context)
    path = root.joinpath(*relative.parts)
    try:
        path.resolve().relative_to(root.resolve())
    except ValueError:
        raise ValueError(f"{context} escapes its frozen root") from None
    if not path.is_file() or path.is_symlink() or path.stat().st_size == 0:
        raise ValueError(f"Missing, empty or symlinked {context}: {path}")
    return path


def _manifest_path(root: Path, item: Any, context: str) -> Path:
    if not isinstance(item, dict) or set(item) != {"relative_path", "sha256"}:
        raise ValueError(f"Malformed {context} binding")
    path = _safe_join(root, item["relative_path"], context)
    if sha256_file(path) != item["sha256"]:
        raise ValueError(f"SHA-256 differs for {context}")
    return path


def evaluate(
    *,
    execution: Path,
    protocol: Path,
    model: Path,
    blind_run: Path,
    custody_path: Path,
    reveal_inputs: Path,
    output: Path,
) -> None:
    for root in (execution, protocol, model, reveal_inputs):
        verify_sha256sums(root, ignore_checksum_file=True)
    protocol_manifest = _load_json(protocol / "protocol_manifest.json")
    execution_manifest = _load_json(execution / "execution_manifest.json")
    custody = _load_json(custody_path)
    reveal = _load_json(reveal_inputs / "reveal_input_manifest.json")
    if (
        protocol_manifest.get("schema_version") != PROTOCOL_SCHEMA
        or execution_manifest.get("schema_version") != EXECUTION_SCHEMA
        or custody.get("schema_version") != CUSTODY_SCHEMA
        or reveal.get("schema_version") != REVEAL_SCHEMA
        or any(
            value.get("holdout_id") != HOLDOUT_ID
            for value in (protocol_manifest, execution_manifest, custody, reveal)
        )
        or any(
            value.get("policy_id") != POLICY_ID
            for value in (protocol_manifest, custody, reveal)
        )
        or execution_manifest.get("protocol_SHA256SUMS_sha256")
        != sha256_file(protocol / "SHA256SUMS")
        or custody.get("frozen_inputs", {}).get("execution_SHA256SUMS_sha256")
        != sha256_file(execution / "SHA256SUMS")
        or custody.get("frozen_inputs", {}).get("protocol_SHA256SUMS_sha256")
        != sha256_file(protocol / "SHA256SUMS")
        or custody.get("frozen_inputs", {}).get(
            "composite_model_SHA256SUMS_sha256"
        )
        != sha256_file(model / "SHA256SUMS")
        or reveal.get("formal_status") != "ready_for_evaluation"
        or reveal.get("custody_manifest_sha256") != sha256_file(custody_path)
    ):
        raise ValueError("Actinidia execution, custody and reveal bindings differ")

    blind_outputs = custody.get("blind_outputs")
    if not isinstance(blind_outputs, dict):
        raise ValueError("Custody lacks blind outputs")
    blind_paths: dict[str, Path] = {}
    for name in ("scores", "score_manifest", "pool_decisions", "pool_manifest"):
        item = blind_outputs.get(name)
        if not isinstance(item, dict) or set(item) != {"relative_path", "sha256"}:
            raise ValueError(f"Malformed custody blind output: {name}")
        path = _safe_join(blind_run, item["relative_path"], f"blind output {name}")
        if sha256_file(path) != item["sha256"]:
            raise ValueError(f"Custody digest differs for blind output: {name}")
        blind_paths[name] = path
    expected_manifest = Path(str(blind_paths["scores"]) + ".manifest.json")
    if expected_manifest.resolve() != blind_paths["score_manifest"].resolve():
        raise ValueError("Custody score manifest is not adjacent to the blind scores")

    inputs = reveal.get("evaluation_inputs")
    if not isinstance(inputs, dict):
        raise ValueError("Reveal manifest lacks evaluation inputs")
    required = {"labels", "primary_pool_score", "legacy_pool_score"}
    if not required <= set(inputs):
        raise ValueError("Reveal manifest lacks primary evaluation inputs")
    resolved = {
        name: _manifest_path(reveal_inputs, item, f"reveal input {name}")
        for name, item in inputs.items()
    }
    policy = _safe_join(
        protocol,
        "protocol_artifacts/config/actinidia_external_validation_policy_v0.5.tsv",
        "frozen Actinidia policy",
    )
    evaluator = Path(__file__).with_name("evaluate_external_v0.5.py")
    if not evaluator.is_file() or evaluator.is_symlink():
        raise ValueError("Frozen generic v0.5 evaluator is missing")

    command = [
        sys.executable,
        str(evaluator),
        "--scores",
        str(blind_paths["scores"]),
        "--labels",
        str(resolved["labels"]),
        "--pool-decisions",
        str(blind_paths["pool_decisions"]),
        "--pool-manifest",
        str(blind_paths["pool_manifest"]),
        "--primary-pool-score",
        str(resolved["primary_pool_score"]),
        "--legacy-pool-score",
        str(resolved["legacy_pool_score"]),
        "--evaluability",
        str(_manifest_path(reveal_inputs, reveal["evaluability"], "evaluability report")),
        "--custody-manifest",
        str(custody_path),
        "--protocol-freeze",
        str(protocol),
        "--composite-model-freeze",
        str(model),
        "--policy",
        str(policy),
        "--input-root",
        f"blind_run={blind_run}",
        "--input-root",
        f"reveal_inputs={reveal_inputs}",
        "--input-root",
        f"protocol_freeze={protocol}",
        "--output-dir",
        str(output),
    ]
    for name in sorted(inputs):
        if not name.startswith("secondary:"):
            continue
        comparator = name.split(":", 1)[1]
        if not comparator:
            raise ValueError("Empty secondary comparator name")
        command.extend(("--secondary-score", f"{comparator}={resolved[name]}"))
    subprocess.run(command, check=True)


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--execution-freeze", required=True)
    parser.add_argument("--protocol-freeze", required=True)
    parser.add_argument("--composite-model-freeze", required=True)
    parser.add_argument("--blind-run", required=True)
    parser.add_argument("--custody-manifest", required=True)
    parser.add_argument("--reveal-inputs", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args(argv)
    evaluate(
        execution=Path(args.execution_freeze),
        protocol=Path(args.protocol_freeze),
        model=Path(args.composite_model_freeze),
        blind_run=Path(args.blind_run),
        custody_path=Path(args.custody_manifest),
        reveal_inputs=Path(args.reveal_inputs),
        output=Path(args.output_dir),
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
