#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import shutil
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any


SCHEMA_VERSION = "ploidypatch.populus_external_input_stage.v0.4"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def verify_sha256sums(root: Path) -> None:
    checksum_path = root / "SHA256SUMS"
    if not checksum_path.is_file() or checksum_path.stat().st_size == 0:
        raise ValueError(f"Missing protocol SHA256SUMS: {root}")
    for line_number, line in enumerate(
        checksum_path.read_text(encoding="utf-8").splitlines(), start=1
    ):
        digest, separator, name = line.partition("  ")
        path = root / name
        if (
            not separator
            or not re.fullmatch(r"[0-9a-f]{64}", digest)
            or not path.is_file()
            or sha256(path) != digest
        ):
            raise ValueError(f"Protocol checksum failure at line {line_number}")


def destination_for(row: dict[str, str]) -> Path:
    species = row["species_id"]
    artifact = row["artifact"]
    source_name = Path(row["source_path"]).name
    if row["role"] == "target" and artifact == "genome":
        return Path("shared_target") / species / source_name
    if row["role"] == "target":
        return Path("evaluator_only/target_complete") / species / source_name
    if row["role"] == "candidate_reference":
        return Path("candidate_only") / species / source_name
    if row["role"] == "evaluator_reference":
        return Path("evaluator_only/truth_references") / species / source_name
    raise ValueError(f"Unknown staged role: {row['role']}")


def read_manifest(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "role",
            "species_id",
            "release",
            "artifact",
            "bytes",
            "sha256",
            "source_path",
        }
        if reader.fieldnames is None or not required <= set(reader.fieldnames):
            raise ValueError("Preflight input manifest lacks staged-role fields")
        rows = list(reader)
    if not rows:
        raise ValueError("Preflight input manifest is empty")
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Copy frozen Populus inputs into disjoint blind/evaluator roles"
    )
    parser.add_argument("--protocol-freeze", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--code-commit", required=True)
    args = parser.parse_args()
    if not re.fullmatch(r"[0-9a-f]{40}", args.code_commit):
        raise ValueError("--code-commit must be a full lowercase Git SHA")
    protocol_root = Path(args.protocol_freeze)
    output = Path(args.output_dir)
    partial = Path(str(output) + ".partial")
    if output.exists() or partial.exists():
        raise FileExistsError("Refusing to overwrite staged Populus inputs")
    verify_sha256sums(protocol_root)
    manifest_path = protocol_root / "preflight_input_manifest.tsv"
    rows = read_manifest(manifest_path)
    destinations = [destination_for(row) for row in rows]
    if len(destinations) != len(set(destinations)):
        raise ValueError("Two source artifacts map to one staged destination")
    partial.mkdir(parents=True)

    def copy_one(item: tuple[dict[str, str], Path]) -> dict[str, Any]:
        row, relative = item
        source = Path(row["source_path"])
        destination = partial / relative
        if (
            not source.is_file()
            or source.stat().st_size != int(row["bytes"])
            or sha256(source) != row["sha256"]
        ):
            raise ValueError(f"Public source differs from frozen preflight: {source}")
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, destination)
        if destination.stat().st_size != int(row["bytes"]) or sha256(destination) != row["sha256"]:
            raise IOError(f"Staged copy verification failed: {destination}")
        return {
            **row,
            "staged_relative_path": relative.as_posix(),
            "staged_sha256": sha256(destination),
        }

    with ThreadPoolExecutor(max_workers=min(8, len(rows))) as executor:
        staged_rows = list(executor.map(copy_one, zip(rows, destinations, strict=True)))
    role_manifest = partial / "role_manifest.tsv"
    fields = [
        "role",
        "species_id",
        "release",
        "artifact",
        "bytes",
        "sha256",
        "source_path",
        "staged_relative_path",
        "staged_sha256",
    ]
    with role_manifest.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(
            sorted(
                staged_rows,
                key=lambda row: (row["role"], row["species_id"], row["artifact"]),
            )
        )
    contract = {
        "schema_version": SCHEMA_VERSION,
        "code_commit": args.code_commit,
        "protocol_SHA256SUMS_sha256": sha256(protocol_root / "SHA256SUMS"),
        "wgd_pairs_enumerated": False,
        "candidate_counts_computed": False,
        "truth_labels_accessed": False,
        "role_boundaries": {
            "shared_target": "target_genome_only",
            "candidate_only": "Salix_candidate_generation_only",
            "evaluator_only": "complete_target_annotation_and_evaluator_references",
            "candidate_evaluator_species_overlap": False,
        },
        "counts": {
            "artifacts": len(staged_rows),
            "bytes": sum(int(row["bytes"]) for row in staged_rows),
        },
    }
    with (partial / "role_contract.json").open(
        "x", encoding="utf-8", newline=""
    ) as handle:
        json.dump(contract, handle, indent=2, sort_keys=True)
        handle.write("\n")
    checksum_path = partial / "SHA256SUMS"
    with checksum_path.open("x", encoding="utf-8", newline="") as handle:
        for path in sorted(item for item in partial.rglob("*") if item.is_file()):
            if path != checksum_path:
                handle.write(f"{sha256(path)}  {path.relative_to(partial).as_posix()}\n")
    os.replace(partial, output)
    for path in output.rglob("*"):
        path.chmod(0o550 if path.is_dir() else 0o440)
    output.chmod(0o550)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
