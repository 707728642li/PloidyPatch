#!/usr/bin/env python3
"""Stage a frozen external-holdout contract into disjoint role directories."""
from __future__ import annotations

import argparse
import csv
from concurrent.futures import ThreadPoolExecutor
import json
import os
from pathlib import Path, PurePosixPath
import re
import shutil
from typing import Any, Iterable

from ploidypatch.artifact_manifest import (
    sha256_file,
    verify_sha256sums,
    write_sha256sums,
)
from ploidypatch.holdout_contract import (
    ArtifactSource,
    HoldoutContract,
    ReferenceContract,
    load_holdout_contract,
    staged_relative_path,
)
from ploidypatch.safe_tar_fasta import extract_single_member_tar_fasta


SCHEMA_VERSION = "ploidypatch.external_holdout_input_stage.v0.5"


def _local_path(root: Path, relative: PurePosixPath) -> Path:
    return root.joinpath(*relative.parts)


def _source_without_symlinks(root: Path, relative: PurePosixPath) -> Path:
    if root.is_symlink():
        raise ValueError(f"Source root may not be a symlink: {root}")
    current = root
    for part in relative.parts:
        current = current / part
        if current.is_symlink():
            raise ValueError(f"Source path contains a symlink: {current}")
    if not current.is_file() or current.stat().st_size == 0:
        raise FileNotFoundError(f"Missing or empty holdout source artifact: {current}")
    try:
        current.resolve().relative_to(root.resolve())
    except ValueError:
        raise ValueError(f"Holdout source escapes source root: {current}") from None
    return current


def _validate_source(
    source_root: Path,
    reference: ReferenceContract,
    artifact_name: str,
    artifact: ArtifactSource,
) -> tuple[ReferenceContract, str, ArtifactSource, Path, PurePosixPath]:
    source = _source_without_symlinks(source_root, artifact.source_relative_path)
    if source.stat().st_size != artifact.bytes:
        raise ValueError(
            f"Frozen source byte count differs for {reference.species_id}/{artifact_name}"
        )
    observed = sha256_file(source)
    if observed != artifact.sha256:
        raise ValueError(
            f"Frozen source SHA-256 differs for {reference.species_id}/{artifact_name}"
        )
    return (
        reference,
        artifact_name,
        artifact,
        source,
        staged_relative_path(reference, artifact_name),
    )


def _copy_one(
    partial: Path,
    item: tuple[ReferenceContract, str, ArtifactSource, Path, PurePosixPath],
) -> dict[str, Any]:
    reference, artifact_name, artifact, source, relative = item
    destination = _local_path(partial, relative)
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() or destination.is_symlink():
        raise FileExistsError(f"Refusing to overwrite staged artifact: {destination}")
    if artifact.container is None:
        shutil.copyfile(source, destination, follow_symlinks=False)
    else:
        extract_single_member_tar_fasta(
            source,
            destination,
            expected_outer_bytes=artifact.bytes,
            expected_outer_sha256=artifact.sha256,
            expected_member_name=artifact.container.member_name.as_posix(),
            expected_member_bytes=artifact.container.member_bytes,
            expected_member_sha256=artifact.container.member_sha256,
        )
    if destination.is_symlink() or not destination.is_file():
        raise ValueError(f"Staged artifact is not a regular file: {destination}")
    if destination.stat().st_size != artifact.staged_bytes:
        raise IOError(f"Staged byte count differs: {destination}")
    observed = sha256_file(destination)
    if observed != artifact.staged_sha256:
        raise IOError(f"Staged SHA-256 differs: {destination}")
    return {
        "role": reference.role,
        "species_id": reference.species_id,
        "release": reference.release,
        "bundle_id": reference.bundle_id,
        "wgdi_prefix": reference.wgdi_prefix,
        "artifact": artifact_name,
        "bytes": artifact.staged_bytes,
        "sha256": artifact.staged_sha256,
        "source_relative_path": artifact.source_relative_path.as_posix(),
        "staged_relative_path": relative.as_posix(),
        "staged_sha256": observed,
    }


def _write_role_manifest(path: Path, rows: list[dict[str, Any]]) -> None:
    fields = (
        "role",
        "species_id",
        "release",
        "bundle_id",
        "wgdi_prefix",
        "artifact",
        "bytes",
        "sha256",
        "source_relative_path",
        "staged_relative_path",
        "staged_sha256",
    )
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(
            sorted(rows, key=lambda row: (row["role"], row["species_id"], row["artifact"]))
        )


def _write_stage_contract(
    path: Path,
    *,
    contract: HoldoutContract,
    contract_path: Path,
    source_root: Path,
    code_commit: str,
    rows: list[dict[str, Any]],
) -> None:
    payload = {
        "schema_version": SCHEMA_VERSION,
        "holdout_contract_schema": contract.schema_version,
        "holdout_id": contract.holdout_id,
        "policy_id": contract.policy_id,
        "test_role": contract.test_role,
        "model_version": contract.model_version,
        "code_commit": code_commit,
        "contract": {
            "path": str(contract_path),
            "bytes": contract_path.stat().st_size,
            "sha256": sha256_file(contract_path),
        },
        "source_root": str(source_root),
        "truth_blind": dict(contract.truth_blind),
        "target_resolved_parameters": {
            "primary_chromosome_count": (
                contract.target_resolved_parameters.primary_chromosome_count
            ),
            "minimum_target_chromosomes_fraction": (
                contract.target_resolved_parameters.minimum_target_chromosomes_fraction
            ),
            "minimum_target_chromosomes": (
                contract.target_resolved_parameters.minimum_target_chromosomes
            ),
        },
        "role_boundaries": {
            "shared_target": "target_genome_only",
            "candidate_only": "candidate_generation_only",
            "evaluator_only": (
                "complete_target_annotation_and_evaluator_truth_references"
            ),
            "candidate_evaluator_species_overlap": False,
        },
        "counts": {
            "references": len(contract.references),
            "artifacts": len(rows),
            "bytes": sum(int(row["bytes"]) for row in rows),
            "target_references": 1,
            "candidate_references": 2,
            "evaluator_references": 2,
        },
    }
    with path.open("x", encoding="utf-8", newline="") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")


def _assert_staged_tree_has_no_symlinks(root: Path) -> None:
    for path in root.rglob("*"):
        if path.is_symlink():
            raise ValueError(f"Symlink is forbidden in staged holdout inputs: {path}")


def stage_inputs(
    *,
    contract_path: str | Path,
    source_root: str | Path,
    output_dir: str | Path,
    code_commit: str,
    copy_workers: int = 4,
) -> Path:
    """Validate, copy, hash and atomically publish one holdout input stage."""

    if re.fullmatch(r"[0-9a-f]{40}", code_commit) is None:
        raise ValueError("code_commit must be a full lowercase Git SHA")
    if isinstance(copy_workers, bool) or not isinstance(copy_workers, int) or copy_workers < 1:
        raise ValueError("copy_workers must be a positive integer")
    contract_file = Path(contract_path)
    contract = load_holdout_contract(contract_file)
    raw_source_root = Path(source_root)
    if not raw_source_root.is_dir() or raw_source_root.is_symlink():
        raise ValueError(f"Source root must be a real directory: {raw_source_root}")
    source = raw_source_root.resolve()
    output = Path(output_dir)
    partial = Path(str(output) + ".partial")
    if output.exists() or output.is_symlink() or partial.exists() or partial.is_symlink():
        raise FileExistsError("Refusing to overwrite staged holdout inputs or partial stage")

    validated: list[
        tuple[ReferenceContract, str, ArtifactSource, Path, PurePosixPath]
    ] = []
    for reference in contract.references:
        for artifact_name, artifact in reference.artifact_items():
            validated.append(
                _validate_source(source, reference, artifact_name, artifact)
            )
    destinations = [item[-1].as_posix().casefold() for item in validated]
    if len(destinations) != len(set(destinations)):
        raise ValueError("Two source artifacts map to one staged destination")

    partial.mkdir(parents=True)
    try:
        workers = min(copy_workers, len(validated))
        with ThreadPoolExecutor(max_workers=workers) as executor:
            rows = list(executor.map(lambda item: _copy_one(partial, item), validated))
        _write_role_manifest(partial / "role_manifest.tsv", rows)
        _write_stage_contract(
            partial / "role_contract.json",
            contract=contract,
            contract_path=contract_file,
            source_root=source,
            code_commit=code_commit,
            rows=rows,
        )
        _assert_staged_tree_has_no_symlinks(partial)
        checksum = write_sha256sums(partial)
        verify_sha256sums(
            partial, checksum, ignore_checksum_file=True
        )
        os.replace(partial, output)
    except BaseException:
        if partial.exists() and not partial.is_symlink():
            shutil.rmtree(partial)
        raise
    return output


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", required=True)
    parser.add_argument("--source-root", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--code-commit", required=True)
    parser.add_argument("--copy-workers", type=int, default=4)
    args = parser.parse_args(argv)
    stage_inputs(
        contract_path=args.contract,
        source_root=args.source_root,
        output_dir=args.output_dir,
        code_commit=args.code_commit,
        copy_workers=args.copy_workers,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
