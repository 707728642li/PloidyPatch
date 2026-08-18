#!/usr/bin/env python3
"""Seal complete target and evaluator-only Coffea data with no candidates."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import shutil
import tempfile
from typing import Iterable

from ploidypatch.artifact_manifest import sha256_file, verify_sha256sums, write_sha256sums
from ploidypatch.coffea_h1_framework import EVALUATOR_ROLE_SCHEMA, verify_protocol


SPECIES = (
    "Coffea_arabica_ET39",
    "Gardenia_jasminoides",
    "Ophiorrhiza_pumila",
)


def reject_candidate_content(root: Path) -> None:
    forbidden = {"candidate_only", "candidate_bua", "candidate_mauritiana"}
    for path in root.rglob("*"):
        if path.is_symlink():
            raise ValueError(f"Symlink is forbidden in Coffea evaluator role: {path}")
        if set(part.casefold() for part in path.relative_to(root).parts) & forbidden:
            raise ValueError(f"Candidate-only path leaked into evaluator role: {path}")


def chmod_read_only(root: Path) -> None:
    if os.name == "nt":
        return
    for path in sorted(root.rglob("*"), reverse=True):
        path.chmod(0o550 if path.is_dir() else 0o440)
    root.chmod(0o550)


def build_evaluator_root(
    *, normalized_root: Path, protein_universes: Path, protocol: Path, output: Path
) -> Path:
    if output.exists() or output.is_symlink():
        raise FileExistsError(f"Refusing to overwrite Coffea evaluator role: {output}")
    verify_sha256sums(normalized_root, ignore_checksum_file=True)
    verify_sha256sums(protein_universes, ignore_checksum_file=True)
    protocol_manifest, contract = verify_protocol(protocol)
    if (
        protocol_manifest.get("normalized_input_SHA256SUMS_sha256")
        != sha256_file(normalized_root / "SHA256SUMS")
        or protocol_manifest.get("protein_universes_SHA256SUMS_sha256")
        != sha256_file(protein_universes / "SHA256SUMS")
    ):
        raise ValueError("Coffea normalized/protein evaluator lineage differs")
    source = normalized_root / "evaluator_only"
    if not source.is_dir() or source.is_symlink():
        raise ValueError("Coffea normalized evaluator-only tree is missing")
    reject_candidate_content(source)
    output.parent.mkdir(parents=True, exist_ok=True)
    working = Path(tempfile.mkdtemp(prefix=f".{output.name}.working.", dir=output.parent))
    try:
        for child in source.iterdir():
            if child.is_dir():
                shutil.copytree(child, working / child.name)
            elif child.is_file() and not child.is_symlink():
                shutil.copyfile(child, working / child.name)
            else:
                raise ValueError(f"Unsafe Coffea evaluator input: {child}")
        proteins = working / "protein_universes"
        proteins.mkdir()
        for species in SPECIES:
            universe = protein_universes / species
            verify_sha256sums(universe, ignore_checksum_file=True)
            shutil.copytree(universe, proteins / species)
        write_sha256sums(proteins)
        manifest = {
            "schema_version": EVALUATOR_ROLE_SCHEMA,
            "holdout_id": contract.holdout_id,
            "policy_id": contract.policy_id,
            "protocol_SHA256SUMS_sha256": sha256_file(protocol / "SHA256SUMS"),
            "normalized_input_SHA256SUMS_sha256": sha256_file(
                normalized_root / "SHA256SUMS"
            ),
            "protein_universes_SHA256SUMS_sha256": sha256_file(
                protein_universes / "SHA256SUMS"
            ),
            "evaluator_protein_universes_SHA256SUMS_sha256": sha256_file(
                proteins / "SHA256SUMS"
            ),
            "roles": ["target_complete", "evaluator_only"],
            "candidate_reference_access": False,
            "blind_candidate_outputs_access": False,
            "ranker_or_model_access": False,
            "h2_or_topology_ranking_access": False,
        }
        (working / "role_manifest.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        reject_candidate_content(working)
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
    parser.add_argument("--normalized-root", required=True)
    parser.add_argument("--protein-universes", required=True)
    parser.add_argument("--protocol-freeze", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args(argv)
    build_evaluator_root(
        normalized_root=Path(args.normalized_root).resolve(),
        protein_universes=Path(args.protein_universes).resolve(),
        protocol=Path(args.protocol_freeze).resolve(),
        output=Path(args.output_dir).resolve(),
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
