#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil
import subprocess
from pathlib import Path
from typing import Any


SCHEMA_VERSION = "ploidypatch.species_applicability_input_stage.v0.6"
FIELDS = ("role", "artifact", "release", "bytes", "sha256", "source_path")
STAGED_NAMES = {
    "genome_fasta_gz": "genome.fa",
    "annotation_gff3_gz": "annotation.gff3",
    "protein_fasta_gz": "protein.fa",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(8 * 1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def reject_symlink_tree(path: Path, label: str) -> None:
    absolute = path.absolute()
    for component in (absolute, *absolute.parents):
        if component.exists() and component.is_symlink():
            raise ValueError(f"{label} traverses a symbolic link: {component}")


def read_sources(path: Path) -> list[dict[str, str]]:
    if not path.is_file() or path.is_symlink() or path.stat().st_size == 0:
        raise ValueError(f"Missing, empty or symlinked source table: {path}")
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != FIELDS:
            raise ValueError("Source table columns differ")
        rows = list(reader)
    if {row["artifact"] for row in rows} != set(STAGED_NAMES) or len(rows) != 3:
        raise ValueError("Source table must contain exactly genome, annotation and protein")
    if any(row["role"] != "target_quality_only" for row in rows):
        raise ValueError("Only the target_quality_only role is permitted")
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Stage byte-locked target-only inputs for a species-applicability audit"
    )
    parser.add_argument("--source-table", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--pigz", default="/usr/bin/pigz", type=Path)
    args = parser.parse_args()

    rows = read_sources(args.source_table)
    output = args.output_dir
    working = Path(f"{output}.working")
    reject_symlink_tree(output, "output directory")
    reject_symlink_tree(working, "working directory")
    if output.exists() or working.exists():
        raise FileExistsError(f"Refusing to overwrite applicability inputs: {output}")
    if not args.pigz.is_file() or args.pigz.is_symlink():
        raise ValueError(f"Missing or symlinked pigz: {args.pigz}")

    try:
        output.parent.mkdir(parents=True, exist_ok=True)
        working.mkdir()
        staged: dict[str, Any] = {}
        for row in sorted(rows, key=lambda item: item["artifact"]):
            source = Path(row["source_path"])
            reject_symlink_tree(source, f"source {row['artifact']}")
            expected_bytes = int(row["bytes"])
            if (
                not source.is_file()
                or source.stat().st_size != expected_bytes
                or sha256(source) != row["sha256"]
            ):
                raise ValueError(f"Source bytes or SHA-256 differ: {source}")
            compressed = working / (STAGED_NAMES[row["artifact"]] + ".gz")
            decompressed = working / STAGED_NAMES[row["artifact"]]
            shutil.copyfile(source, compressed)
            if compressed.stat().st_size != expected_bytes or sha256(compressed) != row["sha256"]:
                raise ValueError(f"Staged compressed bytes differ: {compressed}")
            subprocess.run(
                [str(args.pigz), "-t", str(compressed)], check=True,
                stdout=subprocess.DEVNULL,
            )
            with decompressed.open("xb") as handle:
                subprocess.run(
                    [str(args.pigz), "-dc", str(compressed)],
                    check=True,
                    stdout=handle,
                )
                handle.flush()
                os.fsync(handle.fileno())
            if decompressed.stat().st_size == 0:
                raise ValueError(f"Decompressed staged input is empty: {decompressed}")
            staged[row["artifact"]] = {
                "release": row["release"],
                "source_path": str(source),
                "compressed_file": compressed.name,
                "compressed_bytes": compressed.stat().st_size,
                "compressed_sha256": sha256(compressed),
                "file": decompressed.name,
                "bytes": decompressed.stat().st_size,
                "sha256": sha256(decompressed),
            }

        manifest = {
            "schema_version": SCHEMA_VERSION,
            "role": "target_quality_only",
            "candidate_access": False,
            "truth_or_label_access": False,
            "inputs": {
                "source_table_sha256": sha256(args.source_table),
                "pigz_sha256": sha256(args.pigz),
            },
            "artifacts": staged,
        }
        manifest_path = working / "manifest.json"
        manifest_path.write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
            newline="",
        )
        checksum_path = working / "SHA256SUMS"
        with checksum_path.open("x", encoding="utf-8", newline="") as handle:
            for path in sorted(working.iterdir(), key=lambda item: item.name):
                if path == checksum_path or not path.is_file() or path.is_symlink():
                    continue
                handle.write(f"{sha256(path)}  {path.name}\n")
        os.replace(working, output)
    except Exception:
        if working.exists():
            shutil.rmtree(working)
        raise
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
