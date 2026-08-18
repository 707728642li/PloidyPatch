#!/usr/bin/env python3
"""Seal candidate-safe staged data plus one truth-free blind perturbation."""
from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path, PurePosixPath
import shutil
from typing import Any, Iterable

from ploidypatch.artifact_manifest import (
    sha256_file,
    verify_sha256sums,
    write_sha256sums,
)
from ploidypatch.holdout_contract import load_holdout_contract, safe_relative_path


SCHEMA_VERSION = "ploidypatch.blind_role_manifest.v0.5"
BENCHMARK_SCHEMA = "ploidypatch.blind_benchmark_input.v0.5"
PROTOCOL_SCHEMA = "ploidypatch.external_holdout_protocol_freeze.v0.5"
EXPECTED_BENCHMARK_FILES = frozenset(
    {"perturbed.gff3", "blind_manifest.json", "SHA256SUMS"}
)


def load_json(path: Path) -> dict[str, Any]:
    if not path.is_file() or path.is_symlink() or path.stat().st_size == 0:
        raise ValueError(f"Missing, empty or symlinked JSON input: {path}")
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return value


def target_genome_path(stage: Path, contract: Any) -> Path:
    with (stage / "role_manifest.tsv").open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("Staged role manifest has no header")
        rows = [
            row
            for row in reader
            if row.get("role") == "target"
            and row.get("species_id") == contract.target.species_id
            and row.get("artifact") == "genome"
        ]
    if len(rows) != 1:
        raise ValueError("Staged role manifest lacks one unique target genome")
    relative = safe_relative_path(
        rows[0]["staged_relative_path"], "staged target genome"
    )
    if not relative.parts or relative.parts[0] != "shared_target":
        raise ValueError("Target genome is not in the shared_target role")
    path = stage.joinpath(*relative.parts)
    if not path.is_file() or path.is_symlink():
        raise ValueError("Staged target genome is missing or symlinked")
    if sha256_file(path) != rows[0]["staged_sha256"]:
        raise ValueError("Staged target genome differs from role manifest")
    return path


def validate_benchmark(
    benchmark: Path, stage: Path, contract: Any
) -> dict[str, Any]:
    verify_sha256sums(benchmark, ignore_checksum_file=True)
    files = {
        path.relative_to(benchmark).as_posix()
        for path in benchmark.rglob("*")
        if path.is_file()
    }
    if files != EXPECTED_BENCHMARK_FILES:
        raise ValueError(
            "Blind benchmark must contain exactly perturbed.gff3, "
            "blind_manifest.json and SHA256SUMS"
        )
    manifest = load_json(benchmark / "blind_manifest.json")
    perturbed = benchmark / "perturbed.gff3"
    target_genome_path(stage, contract)
    expected_target_sha = manifest.get("target_genome", {}).get("sha256")
    if (
        manifest.get("schema_version") != BENCHMARK_SCHEMA
        or manifest.get("truth_access") is not False
        or manifest.get("complete_target_annotation_access") is not False
        or manifest.get("perturbed_annotation")
        != {"file_name": "perturbed.gff3", "sha256": sha256_file(perturbed)}
        or manifest.get("target_genome", {}).get("mount_role")
        != "shared_target_genome"
        or not isinstance(expected_target_sha, str)
        or len(expected_target_sha) != 64
        or any(character not in "0123456789abcdef" for character in expected_target_sha)
    ):
        raise ValueError("Blind benchmark manifest or target-genome binding differs")
    return manifest


def reject_forbidden_names(root: Path) -> None:
    forbidden = {"evaluator", "evaluator_only", "truth", "labels", "target_complete"}
    for path in root.rglob("*"):
        if path.is_symlink():
            raise ValueError(f"Symlink is forbidden in blind role input: {path}")
        folded = {part.casefold() for part in path.relative_to(root).parts}
        if folded & forbidden:
            raise ValueError(f"Forbidden evaluator/truth path in blind role: {path}")


def chmod_read_only(root: Path) -> None:
    for path in sorted(root.rglob("*"), reverse=True):
        path.chmod(0o550 if path.is_dir() else 0o440)
    root.chmod(0o550)


def build_role_root(
    *, staged_inputs: Path, blind_benchmark: Path, protocol: Path, output: Path
) -> Path:
    partial = Path(str(output) + ".partial")
    if output.exists() or output.is_symlink() or partial.exists() or partial.is_symlink():
        raise FileExistsError("Refusing to overwrite blind role root or partial output")
    for root in (staged_inputs, protocol):
        verify_sha256sums(root, ignore_checksum_file=True)
    contract = load_holdout_contract(protocol / "contract.json")
    protocol_manifest = load_json(protocol / "protocol_manifest.json")
    if (
        protocol_manifest.get("schema_version") != PROTOCOL_SCHEMA
        or protocol_manifest.get("holdout_id") != contract.holdout_id
        or protocol_manifest.get("contract_sha256")
        != sha256_file(protocol / "contract.json")
        or protocol_manifest.get("staged_input_SHA256SUMS_sha256")
        != sha256_file(staged_inputs / "SHA256SUMS")
        or (staged_inputs / "role_manifest.tsv").read_bytes()
        != (protocol / "role_manifest.tsv").read_bytes()
        or (staged_inputs / "role_contract.json").read_bytes()
        != (protocol / "role_contract.json").read_bytes()
    ):
        raise ValueError("Staged input, protocol and contract bindings differ")
    benchmark_manifest = validate_benchmark(blind_benchmark, staged_inputs, contract)
    for role in ("shared_target", "candidate_only"):
        source = staged_inputs / role
        if not source.is_dir() or source.is_symlink():
            raise ValueError(f"Staged inputs lack safe role directory: {role}")
        reject_forbidden_names(source)

    partial.mkdir(parents=True)
    try:
        for role in ("shared_target", "candidate_only"):
            shutil.copytree(
                staged_inputs / role,
                partial / role,
                copy_function=shutil.copyfile,
            )
        shutil.copytree(
            blind_benchmark,
            partial / "blind_benchmark",
            copy_function=shutil.copyfile,
        )
        role_manifest = {
            "schema_version": SCHEMA_VERSION,
            "holdout_id": contract.holdout_id,
            "policy_id": contract.policy_id,
            "model_version": contract.model_version,
            "code_commit": protocol_manifest["code_commit"],
            "contract_sha256": sha256_file(protocol / "contract.json"),
            "protocol_SHA256SUMS_sha256": sha256_file(protocol / "SHA256SUMS"),
            "staged_input_SHA256SUMS_sha256": sha256_file(
                staged_inputs / "SHA256SUMS"
            ),
            "blind_benchmark_SHA256SUMS_sha256": sha256_file(
                blind_benchmark / "SHA256SUMS"
            ),
            "blind_benchmark_manifest_sha256": sha256_file(
                blind_benchmark / "blind_manifest.json"
            ),
            "target_genome_sha256": benchmark_manifest["target_genome"]["sha256"],
            "roles": ["shared_target", "candidate_only", "blind_benchmark"],
            "truth_access": False,
            "complete_target_annotation_present": False,
            "evaluator_references_present": False,
            "network_access": False,
        }
        (partial / "role_manifest.json").write_text(
            json.dumps(role_manifest, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        reject_forbidden_names(partial)
        write_sha256sums(partial)
        verify_sha256sums(partial, ignore_checksum_file=True)
        chmod_read_only(partial)
        os.replace(partial, output)
    except BaseException:
        if partial.exists() and not partial.is_symlink():
            shutil.rmtree(partial)
        raise
    return output


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--staged-inputs", required=True)
    parser.add_argument("--blind-benchmark-root", required=True)
    parser.add_argument("--protocol-freeze", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args(argv)
    build_role_root(
        staged_inputs=Path(args.staged_inputs).resolve(),
        blind_benchmark=Path(args.blind_benchmark_root).resolve(),
        protocol=Path(args.protocol_freeze).resolve(),
        output=Path(args.output_dir).resolve(),
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
