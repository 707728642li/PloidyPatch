#!/usr/bin/env python3
"""Seal candidate-safe Coffea data and exact provider-protein subsets."""
from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path
import shutil
import tempfile
from typing import Any, Iterable

from ploidypatch.artifact_manifest import sha256_file, verify_sha256sums, write_sha256sums
from ploidypatch.coffea_h1_framework import ROLE_SCHEMA, verify_protocol


BENCHMARK_SCHEMA = "ploidypatch.coffea_h1_blind_benchmark_input.v1.0"
SPECIES = ("Coffea_eugenioides_BuA", "Coffea_mauritiana")
EXPECTED_BENCHMARK_FILES = frozenset(
    {"perturbed.gff3", "blind_manifest.json", "SHA256SUMS"}
)


def reject_tree(root: Path) -> None:
    forbidden = {"evaluator_only", "target_complete", "truth", "labels"}
    for path in root.rglob("*"):
        if path.is_symlink():
            raise ValueError(f"Symlink is forbidden in Coffea blind role: {path}")
        if set(part.casefold() for part in path.relative_to(root).parts) & forbidden:
            raise ValueError(f"Forbidden Coffea evaluator/truth path: {path}")


def validate_benchmark(root: Path) -> dict[str, Any]:
    verify_sha256sums(root, ignore_checksum_file=True)
    files = {
        path.relative_to(root).as_posix()
        for path in root.rglob("*")
        if path.is_file()
    }
    if files != EXPECTED_BENCHMARK_FILES:
        raise ValueError("Coffea blind benchmark file universe differs")
    manifest = json.loads((root / "blind_manifest.json").read_text(encoding="utf-8"))
    target_sha = manifest.get("target_genome", {}).get("sha256")
    if (
        manifest.get("schema_version") != BENCHMARK_SCHEMA
        or manifest.get("truth_access") is not False
        or manifest.get("complete_target_annotation_access") is not False
        or manifest.get("ranker_access") is not False
        or manifest.get("h2_or_topology_ranking_access") is not False
        or manifest.get("perturbed_annotation")
        != {
            "file_name": "perturbed.gff3",
            "sha256": sha256_file(root / "perturbed.gff3"),
        }
        or manifest.get("target_genome", {}).get("mount_role")
        != "shared_target_genome"
        or not isinstance(target_sha, str)
        or len(target_sha) != 64
    ):
        raise ValueError("Coffea blind benchmark manifest differs")
    return manifest


def filtered_role_manifest(source: Path, output: Path) -> None:
    with source.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or ())
        rows = [
            row
            for row in reader
            if row.get("role") == "candidate_reference"
            or (row.get("role") == "target" and row.get("artifact") == "genome")
        ]
    expected = 7
    if len(rows) != expected or len(fields) != len(set(fields)):
        raise ValueError("Coffea candidate-safe role manifest is not exact 1+2x3")
    with output.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def chmod_read_only(root: Path) -> None:
    if os.name == "nt":
        return
    for path in sorted(root.rglob("*"), reverse=True):
        path.chmod(0o550 if path.is_dir() else 0o440)
    root.chmod(0o550)


def build_role_root(
    *, staged_inputs: Path, protein_universes: Path, blind_benchmark: Path,
    protocol: Path, output: Path
) -> Path:
    if output.exists() or output.is_symlink():
        raise FileExistsError(f"Refusing to overwrite Coffea blind role root: {output}")
    verify_sha256sums(staged_inputs, ignore_checksum_file=True)
    verify_sha256sums(protein_universes, ignore_checksum_file=True)
    protocol_manifest, contract = verify_protocol(protocol)
    if (
        protocol_manifest.get("staged_input_SHA256SUMS_sha256")
        != sha256_file(staged_inputs / "SHA256SUMS")
        or protocol_manifest.get("protein_universes_SHA256SUMS_sha256")
        != sha256_file(protein_universes / "SHA256SUMS")
        or (protocol / "role_manifest.tsv").read_bytes()
        != (staged_inputs / "role_manifest.tsv").read_bytes()
    ):
        raise ValueError("Coffea staged/protein inputs differ from protocol lineage")
    benchmark_manifest = validate_benchmark(blind_benchmark)
    for role in ("shared_target", "candidate_only"):
        if not (staged_inputs / role).is_dir():
            raise ValueError(f"Missing Coffea candidate-safe role: {role}")
        reject_tree(staged_inputs / role)
    output.parent.mkdir(parents=True, exist_ok=True)
    working = Path(tempfile.mkdtemp(prefix=f".{output.name}.working.", dir=output.parent))
    try:
        shutil.copytree(staged_inputs / "shared_target", working / "shared_target")
        shutil.copytree(staged_inputs / "candidate_only", working / "candidate_only")
        shutil.copytree(blind_benchmark, working / "blind_benchmark")
        filtered_role_manifest(staged_inputs / "role_manifest.tsv", working / "role_manifest.tsv")
        protein_root = working / "candidate_only/protein_universes"
        protein_root.mkdir()
        for species in SPECIES:
            universe = protein_universes / species
            verify_sha256sums(universe, ignore_checksum_file=True)
            shutil.copytree(universe, protein_root / species)
        write_sha256sums(protein_root)
        manifest = {
            "schema_version": ROLE_SCHEMA,
            "holdout_id": contract.holdout_id,
            "policy_id": contract.policy_id,
            "protocol_profile": "core_H1_known_subgenome_no_ranker",
            "contract_sha256": sha256_file(protocol / "contract.json"),
            "protocol_SHA256SUMS_sha256": sha256_file(protocol / "SHA256SUMS"),
            "staged_input_SHA256SUMS_sha256": sha256_file(staged_inputs / "SHA256SUMS"),
            "protein_universes_SHA256SUMS_sha256": sha256_file(
                protein_universes / "SHA256SUMS"
            ),
            "candidate_protein_universes_SHA256SUMS_sha256": sha256_file(
                protein_root / "SHA256SUMS"
            ),
            "blind_benchmark_SHA256SUMS_sha256": sha256_file(
                blind_benchmark / "SHA256SUMS"
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
        (working / "role_manifest.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        reject_tree(working)
        write_sha256sums(working)
        verify_sha256sums(working, ignore_checksum_file=True)
        chmod_read_only(working)
        os.replace(working, output)
    except BaseException:
        shutil.rmtree(working, ignore_errors=True)
        raise
    return output


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--staged-inputs", required=True)
    parser.add_argument("--protein-universes", required=True)
    parser.add_argument("--blind-benchmark-root", required=True)
    parser.add_argument("--protocol-freeze", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args(argv)
    build_role_root(
        staged_inputs=Path(args.staged_inputs).resolve(),
        protein_universes=Path(args.protein_universes).resolve(),
        blind_benchmark=Path(args.blind_benchmark_root).resolve(),
        protocol=Path(args.protocol_freeze).resolve(),
        output=Path(args.output_dir).resolve(),
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
