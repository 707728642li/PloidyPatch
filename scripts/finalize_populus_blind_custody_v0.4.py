#!/usr/bin/env python3
"""Validate and seal a Populus blind run before any evaluator reveal."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import shlex
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable


SCHEMA_VERSION = "ploidypatch.blind_run_custody.v1"
SCORE_SCHEMA = "ploidypatch.conflict_winner_guard_scores.v1"
POOL_SCHEMA = "ploidypatch.method_candidate_pool.v2"
EXECUTION_SCHEMA = "ploidypatch.populus_execution_freeze.v0.4"

RANKING_RELATIVE = Path(
    "results/copy_collapse/external/populus_v0.4_blind_rankings"
)
POOL_RELATIVE = Path(
    "results/copy_collapse/external/populus_v0.4_method_trio/consensus/primary_union/blind"
)
SCORES_RELATIVE = RANKING_RELATIVE / "scores/v04.tsv"
SCORE_MANIFEST_RELATIVE = RANKING_RELATIVE / "scores/v04.tsv.manifest.json"
POOL_DECISIONS_RELATIVE = POOL_RELATIVE / "decisions.tsv"
POOL_MANIFEST_RELATIVE = POOL_RELATIVE / "candidate.gff3.manifest.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return value


def verify_sha256sums(root: Path) -> None:
    checksum = root / "SHA256SUMS"
    if not checksum.is_file() or checksum.stat().st_size == 0:
        raise ValueError(f"Missing SHA256SUMS: {root}")
    for number, line in enumerate(checksum.read_text(encoding="utf-8").splitlines(), 1):
        digest, separator, relative = line.partition("  ")
        path = root / relative
        if (
            not separator
            or not re.fullmatch(r"[0-9a-f]{64}", digest)
            or not path.is_file()
            or sha256(path) != digest
        ):
            raise ValueError(f"Checksum failure at line {number}: {root}")


def accepted_decision_digests(path: Path) -> set[str]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or not {
            "consensus_digest",
            "status",
        } <= set(reader.fieldnames):
            raise ValueError("Pool decisions lack consensus_digest/status")
        rows = list(reader)
    digests = {
        row["consensus_digest"] for row in rows if row["status"] == "accepted"
    }
    if not digests or len(digests) != sum(row["status"] == "accepted" for row in rows):
        raise ValueError("Accepted pool decisions are empty or duplicate")
    return digests


def score_digests(path: Path) -> set[str]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        forbidden = [field for field in fields if "truth" in field.lower() or "label" in field.lower()]
        if forbidden or "candidate_digest" not in fields:
            raise ValueError("Blind score columns contain truth/label data or lack digest")
        rows = list(reader)
    digests = {row["candidate_digest"] for row in rows}
    if not digests or len(digests) != len(rows):
        raise ValueError("Blind scores are empty or contain duplicate candidate digests")
    return digests


def tree_digest(root: Path, excluded: set[Path]) -> tuple[str, int, int]:
    digest = hashlib.sha256()
    files = 0
    total_bytes = 0
    for path in sorted(item for item in root.rglob("*") if item.is_file()):
        if path in excluded:
            continue
        relative = path.relative_to(root).as_posix()
        item_hash = sha256(path)
        digest.update(f"{item_hash}  {relative}\n".encode("utf-8"))
        files += 1
        total_bytes += path.stat().st_size
    if not files:
        raise ValueError("Blind output tree is empty")
    return digest.hexdigest(), files, total_bytes


def reject_forbidden_outputs(root: Path) -> None:
    forbidden_parts = {"evaluator_only", "truth", "labels", "target_complete"}
    for path in root.rglob("*"):
        lowered = {part.lower() for part in path.relative_to(root).parts}
        if lowered & forbidden_parts:
            raise ValueError(f"Forbidden evaluator artifact in blind output: {path}")


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--blind-project-root", required=True)
    parser.add_argument("--execution-freeze", required=True)
    parser.add_argument("--protocol-freeze", required=True)
    parser.add_argument("--composite-model-freeze", required=True)
    parser.add_argument("--shared-target", required=True)
    parser.add_argument("--candidate-only", required=True)
    parser.add_argument("--blind-annotation-root", required=True)
    parser.add_argument("--blind-role-checksums", required=True)
    parser.add_argument("--runner-command", required=True)
    parser.add_argument("--bwrap-version", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args(argv)

    project = Path(args.blind_project_root).resolve()
    execution = Path(args.execution_freeze).resolve()
    protocol = Path(args.protocol_freeze).resolve()
    model = Path(args.composite_model_freeze).resolve()
    output = Path(args.output).resolve()
    blind_role_checksums = Path(args.blind_role_checksums).resolve()
    if not blind_role_checksums.is_file():
        raise FileNotFoundError("Missing frozen blind-role checksum list")
    if output.exists():
        raise FileExistsError("Refusing to overwrite blind custody manifest")
    for root in (execution, protocol, model):
        verify_sha256sums(root)
    execution_manifest = load_json(execution / "execution_manifest.json")
    if execution_manifest.get("schema_version") != EXECUTION_SCHEMA:
        raise ValueError("Wrong Populus execution freeze schema")
    runner_command = Path(args.runner_command)
    command_text = runner_command.read_text(encoding="utf-8")
    command_tokens = shlex.split(command_text)
    if (
        "--unshare-all" not in command_tokens
        or "--unshare-net" not in command_tokens
        or "--share-net" in command_tokens
        or any(
            marker in token.lower()
            for token in command_tokens
            for marker in ("/nas_data", "/evaluator_only", "/target_complete", "/truth_references")
        )
    ):
        raise ValueError("Recorded bubblewrap command violates blind isolation")

    safe_mounts = {
        "shared_target": Path(args.shared_target).resolve(),
        "candidate_only": Path(args.candidate_only).resolve(),
        "blind_perturbed_annotation": Path(args.blind_annotation_root).resolve(),
        "frozen_code": (execution / "source").resolve(),
        "frozen_model": model,
        "frozen_execution_contract": execution,
        "frozen_protocol": protocol,
    }
    if (
        safe_mounts["shared_target"].name != "shared_target"
        or safe_mounts["candidate_only"].name != "candidate_only"
        or safe_mounts["blind_perturbed_annotation"].name != "blind"
    ):
        raise ValueError("Blind data mounts do not use the frozen role directories")
    environment_rows = execution_manifest.get("environments")
    if not isinstance(environment_rows, list) or not environment_rows:
        raise ValueError("Execution freeze lacks environment bindings")
    for row in environment_rows:
        if not isinstance(row, dict) or not row.get("name") or not row.get("host_prefix"):
            raise ValueError("Malformed frozen environment binding")
        safe_mounts[f"frozen_environment:{row['name']}"] = Path(
            str(row["host_prefix"])
        ).resolve()
    for role, path in safe_mounts.items():
        lowered_parts = {part.lower() for part in path.parts}
        if (
            not path.exists()
            or str(path) == "/nas_data"
            or str(path).startswith("/nas_data/")
            or lowered_parts
            & {"evaluator_only", "target_complete", "truth", "truth_references", "labels"}
        ):
            raise ValueError(f"Unsafe or missing blind mount {role}: {path}")

    reject_forbidden_outputs(project)
    scores = project / SCORES_RELATIVE
    score_manifest_path = project / SCORE_MANIFEST_RELATIVE
    decisions = project / POOL_DECISIONS_RELATIVE
    pool_manifest_path = project / POOL_MANIFEST_RELATIVE
    for path in (scores, score_manifest_path, decisions, pool_manifest_path):
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(f"Missing blind output: {path}")
    pool_manifest = load_json(pool_manifest_path)
    score_manifest = load_json(score_manifest_path)
    if (
        pool_manifest.get("schema_version") != POOL_SCHEMA
        or pool_manifest.get("outputs", {}).get("decisions", {}).get("sha256")
        != sha256(decisions)
    ):
        raise ValueError("Blind pool manifest does not bind decisions")
    if (
        score_manifest.get("schema_version") != SCORE_SCHEMA
        or score_manifest.get("truth_access") is not False
        or score_manifest.get("inputs", {}).get("pool_decisions") != sha256(decisions)
        or score_manifest.get("inputs", {}).get("pool_manifest") != sha256(pool_manifest_path)
        or score_manifest.get("outputs", {}).get("scores", {}).get("sha256")
        != sha256(scores)
    ):
        raise ValueError("Blind score manifest does not bind pool/scores")
    if score_digests(scores) != accepted_decision_digests(decisions):
        raise ValueError("Blind score and accepted-pool candidate universes differ")

    tree_hash, file_count, total_bytes = tree_digest(project, {output})
    custody = {
        "schema_version": SCHEMA_VERSION,
        "runner_identity": "bubblewrap_populus_v0.4_blind_runner",
        "frozen_before_truth_reveal_at": datetime.now(timezone.utc)
        .replace(microsecond=0)
        .isoformat()
        .replace("+00:00", "Z"),
        "bubblewrap": {
            "version": args.bwrap_version,
            "required_flags": ["--unshare-all", "--unshare-net"],
            "command_sha256": sha256(runner_command),
        },
        "truth_mounted": False,
        "complete_target_annotation_mounted": False,
        "evaluator_references_mounted": False,
        "nas_data_mounted": False,
        "network_access": False,
        "mounts": {
            role: {"host_path": str(path), "read_only": True}
            for role, path in sorted(safe_mounts.items())
        },
        "frozen_inputs": {
            "execution_SHA256SUMS_sha256": sha256(execution / "SHA256SUMS"),
            "protocol_SHA256SUMS_sha256": sha256(protocol / "SHA256SUMS"),
            "composite_model_SHA256SUMS_sha256": sha256(model / "SHA256SUMS"),
            "blind_role_SHA256SUMS_sha256": sha256(blind_role_checksums),
        },
        "blind_outputs": {
            "scores_relative_path": SCORES_RELATIVE.as_posix(),
            "scores_sha256": sha256(scores),
            "score_manifest_relative_path": SCORE_MANIFEST_RELATIVE.as_posix(),
            "score_manifest_sha256": sha256(score_manifest_path),
            "pool_decisions_relative_path": POOL_DECISIONS_RELATIVE.as_posix(),
            "pool_decisions_sha256": sha256(decisions),
            "pool_manifest_relative_path": POOL_MANIFEST_RELATIVE.as_posix(),
            "pool_manifest_sha256": sha256(pool_manifest_path),
            "output_tree_sha256": tree_hash,
            "file_count": file_count,
            "bytes": total_bytes,
        },
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL
    descriptor = os.open(output, flags, 0o440)
    with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as handle:
        json.dump(custody, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
