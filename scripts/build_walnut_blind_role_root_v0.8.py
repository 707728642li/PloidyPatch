#!/usr/bin/env python3
"""Seal exactly three candidate-safe roles for Walnut core-H1 execution."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import shutil
from typing import Any, Iterable

from ploidypatch.artifact_manifest import sha256_file, verify_sha256sums, write_sha256sums
from ploidypatch.walnut_h1_framework import ROLE_SCHEMA, load_json, verify_protocol


BENCHMARK_SCHEMA = "ploidypatch.walnut_h1_blind_benchmark_input.v0.8"
EXPECTED_BENCHMARK_FILES = frozenset(
    {"perturbed.gff3", "blind_manifest.json", "SHA256SUMS"}
)


def reject_role_tree(root: Path) -> None:
    forbidden_names = {"evaluator_only", "target_complete", "truth", "labels"}
    for path in root.rglob("*"):
        if path.is_symlink():
            raise ValueError(f"Symlink is forbidden in blind role: {path}")
        if set(part.casefold() for part in path.relative_to(root).parts) & forbidden_names:
            raise ValueError(f"Forbidden evaluator/truth path in blind role: {path}")


def validate_benchmark(root: Path) -> dict[str, Any]:
    verify_sha256sums(root, ignore_checksum_file=True)
    files = {
        path.relative_to(root).as_posix()
        for path in root.rglob("*") if path.is_file()
    }
    if files != EXPECTED_BENCHMARK_FILES:
        raise ValueError("Blind benchmark file universe differs")
    manifest = load_json(root / "blind_manifest.json")
    perturbed = root / "perturbed.gff3"
    target_sha = manifest.get("target_genome", {}).get("sha256")
    if (
        manifest.get("schema_version") != BENCHMARK_SCHEMA
        or manifest.get("truth_access") is not False
        or manifest.get("complete_target_annotation_access") is not False
        or manifest.get("ranker_access") is not False
        or manifest.get("h2_or_topology_ranking_access") is not False
        or manifest.get("perturbed_annotation")
        != {"file_name": "perturbed.gff3", "sha256": sha256_file(perturbed)}
        or manifest.get("target_genome", {}).get("mount_role")
        != "shared_target_genome"
        or not isinstance(target_sha, str)
        or len(target_sha) != 64
        or any(character not in "0123456789abcdef" for character in target_sha)
    ):
        raise ValueError("Blind benchmark is not exact truth-free Walnut input")
    return manifest


def chmod_read_only(root: Path) -> None:
    if os.name == "nt":
        return
    for path in sorted(root.rglob("*"), reverse=True):
        path.chmod(0o550 if path.is_dir() else 0o440)
    root.chmod(0o550)


def build_role_root(
    *, staged_inputs: Path, blind_benchmark: Path, protocol: Path, output: Path
) -> Path:
    partial = Path(str(output) + ".partial")
    if output.exists() or output.is_symlink() or partial.exists() or partial.is_symlink():
        raise FileExistsError("Refusing to overwrite blind role root or partial")
    protocol_manifest, contract = verify_protocol(protocol)
    verify_sha256sums(staged_inputs, ignore_checksum_file=True)
    if (
        protocol_manifest.get("staged_input_SHA256SUMS_sha256")
        != sha256_file(staged_inputs / "SHA256SUMS")
        or (protocol / "role_manifest.tsv").read_bytes()
        != (staged_inputs / "role_manifest.tsv").read_bytes()
        or (protocol / "role_contract.json").read_bytes()
        != (staged_inputs / "role_contract.json").read_bytes()
    ):
        raise ValueError("Staged inputs differ from protocol lineage")
    benchmark_manifest = validate_benchmark(blind_benchmark)
    for role in ("shared_target", "candidate_only"):
        role_path = staged_inputs / role
        if not role_path.is_dir() or role_path.is_symlink():
            raise ValueError(f"Staged inputs lack safe role: {role}")
        reject_role_tree(role_path)

    partial.mkdir(parents=True)
    try:
        for role in ("shared_target", "candidate_only"):
            shutil.copytree(staged_inputs / role, partial / role,
                            copy_function=shutil.copyfile)
        shutil.copytree(blind_benchmark, partial / "blind_benchmark",
                        copy_function=shutil.copyfile)
        manifest = {
            "schema_version": ROLE_SCHEMA,
            "holdout_id": contract.holdout_id,
            "policy_id": contract.policy_id,
            "protocol_profile": "core_H1_only_no_ranker",
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
            "nas_data_present": False,
            "network_access": False,
            "ranker_or_model_present": False,
            "h2_or_topology_ranking_present": False,
        }
        (partial / "role_manifest.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        reject_role_tree(partial)
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
