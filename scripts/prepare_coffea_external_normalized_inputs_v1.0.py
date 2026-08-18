#!/usr/bin/env python3
"""Build role-separated, primary-only normalized Coffea validation inputs."""
from __future__ import annotations

import csv
import json
import os
from pathlib import Path
import shutil
import sys
import tempfile
from typing import Any

from ploidypatch.artifact_manifest import (
    sha256_file,
    verify_sha256sums,
    write_sha256sums,
)
from ploidypatch.holdout_contract import load_holdout_contract
from ploidypatch.io import iter_fasta, open_text
from ploidypatch.known_subgenome_h1 import HOLDOUT_ID, POLICY_ID
from ploidypatch.normalize import prepare_primary_annotation_bundle, read_primary_seqid_table
from ploidypatch.seqid_alias import rewrite_gff3_seqids_exact


ROLE_LAYOUT = {
    ("target", "Coffea_arabica_ET39"): "evaluator_only/normalized/target_complete",
    (
        "candidate_reference",
        "Coffea_eugenioides_BuA",
    ): "blind/candidate_only/candidate_bua",
    (
        "candidate_reference",
        "Coffea_mauritiana",
    ): "blind/candidate_only/candidate_mauritiana",
    (
        "evaluator_reference",
        "Gardenia_jasminoides",
    ): "evaluator_only/normalized/truth_references/evaluator_gardenia",
    (
        "evaluator_reference",
        "Ophiorrhiza_pumila",
    ): "evaluator_only/normalized/truth_references/evaluator_ophiorrhiza",
}


def write_primary_genome(source: Path, table: Path, output: Path) -> None:
    seqids, _ = read_primary_seqid_table(table)
    seen: set[str] = set()
    fai: list[str] = []
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("xb") as handle:
        for seqid, _header, sequence in iter_fasta(source):
            if seqid not in seqids:
                continue
            if seqid in seen or not sequence:
                raise ValueError(f"Duplicate or empty primary target seqid: {seqid}")
            seen.add(seqid)
            handle.write(f">{seqid}\n".encode())
            offset = handle.tell()
            encoded = sequence.encode("ascii")
            for index in range(0, len(encoded), 60):
                handle.write(encoded[index : index + 60] + b"\n")
            line_bases = min(60, len(encoded))
            fai.append(
                f"{seqid}\t{len(encoded)}\t{offset}\t{line_bases}\t{line_bases + 1}\n"
            )
    if seen != set(seqids):
        raise ValueError(
            "Target primary genome seqids differ: "
            + ",".join(sorted(set(seqids) - seen))
        )
    output.with_suffix(output.suffix + ".fai").write_text(
        "".join(fai), encoding="utf-8"
    )


def read_role_manifest(root: Path) -> dict[tuple[str, str, str], dict[str, str]]:
    with (root / "role_manifest.tsv").open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "role",
            "species_id",
            "artifact",
            "staged_relative_path",
            "bytes",
            "sha256",
        }
        if not required <= set(reader.fieldnames or ()):
            raise ValueError("Coffea staged role manifest fields differ")
        rows = list(reader)
    indexed = {
        (row["role"], row["species_id"], row["artifact"]): row for row in rows
    }
    if len(rows) != 15 or len(indexed) != 15:
        raise ValueError("Coffea staged role manifest must contain 15 artifacts")
    return indexed


def main(argv: list[str] | None = None) -> int:
    values = sys.argv[1:] if argv is None else argv
    if len(values) != 1:
        raise SystemExit(
            "usage: prepare_coffea_external_normalized_inputs_v1.0.py PROJECT_ROOT"
        )
    project_root = Path(values[0]).resolve()
    staged = Path(os.environ["PLOIDYPATCH_STAGED_INPUT_ROOT"]).resolve()
    execution = Path(os.environ["PLOIDYPATCH_EXECUTION_FREEZE"]).resolve()
    contract_path = Path(os.environ["PLOIDYPATCH_HOLDOUT_CONTRACT"]).resolve()
    output = project_root / "data/derived/external_inputs/coffea/v1.0"
    if output.exists() or output.is_symlink():
        raise FileExistsError(f"Refusing to overwrite Coffea normalized inputs: {output}")
    for root in (staged, execution):
        verify_sha256sums(root, ignore_checksum_file=True)
    contract = load_holdout_contract(contract_path)
    if contract.holdout_id != HOLDOUT_ID or contract.policy_id != POLICY_ID:
        raise ValueError("Coffea normalization received a different contract")
    rows = read_role_manifest(staged)
    output.parent.mkdir(parents=True, exist_ok=True)
    working = Path(tempfile.mkdtemp(prefix=".v1.0.working.", dir=output.parent))
    try:
        input_records: dict[str, Any] = {}
        for reference in contract.references:
            key = (reference.role, reference.species_id)
            if key not in ROLE_LAYOUT:
                raise ValueError(f"Unexpected Coffea role/species: {key}")
            bundle = working / ROLE_LAYOUT[key]
            artifacts: dict[str, Path] = {}
            for artifact in ("genome", "gff3", "protein"):
                row = rows[(reference.role, reference.species_id, artifact)]
                relative = Path(row["staged_relative_path"])
                if relative.is_absolute() or ".." in relative.parts:
                    raise ValueError(f"Unsafe staged Coffea path: {relative}")
                path = staged / relative
                if (
                    not path.is_file()
                    or path.is_symlink()
                    or path.stat().st_size != int(row["bytes"])
                    or sha256_file(path) != row["sha256"]
                ):
                    raise ValueError(f"Staged Coffea artifact differs: {path}")
                artifacts[artifact] = path
                input_records[f"{reference.bundle_id}:{artifact}"] = {
                    "bytes": path.stat().st_size,
                    "sha256": sha256_file(path),
                }
            primary = execution / "source" / reference.primary_seqid_table
            if not primary.is_file() or primary.is_symlink():
                raise ValueError(f"Frozen primary seqid table is missing: {primary}")

            gff_input = artifacts["gff3"]
            if reference.species_id == "Gardenia_jasminoides":
                compatibility = working / "evaluator_only/compatibility/gardenia_seqids"
                compatibility.mkdir(parents=True)
                gff_input = compatibility / "primary_alias_normalized.gff3"
                rewrite_gff3_seqids_exact(
                    input_gff_path=artifacts["gff3"],
                    alias_table_path=(
                        execution
                        / "source/config/coffea_gardenia_seqid_aliases_v1.0.tsv"
                    ),
                    primary_seqid_table_path=primary,
                    output_gff_path=gff_input,
                    manifest_path=compatibility / "manifest.json",
                )

            prepare_primary_annotation_bundle(
                gff_path=gff_input,
                genome_path=artifacts["genome"],
                primary_seqid_table_path=primary,
                output_dir=bundle,
                canonical_fasta_headers=True,
            )
            with open_text(artifacts["protein"]) as source_handle, (
                bundle / "provider.protein.fa"
            ).open("x", encoding="utf-8", newline="") as destination:
                shutil.copyfileobj(source_handle, destination)
            if not (bundle / "provider.protein.fa").stat().st_size:
                raise ValueError(f"Empty provider protein output for {reference.species_id}")

            if reference.role == "target":
                shared = working / "blind/shared_target/target_et39"
                write_primary_genome(
                    artifacts["genome"],
                    primary,
                    shared / "primary_chromosomes.genome.fa",
                )
                (shared / "manifest.json").write_text(
                    json.dumps(
                        {
                            "schema_version": (
                                "ploidypatch.coffea_blind_target_genome.v1.0"
                            ),
                            "holdout_id": HOLDOUT_ID,
                            "truth_access": False,
                            "complete_target_annotation_access": False,
                            "primary_seqid_table_sha256": sha256_file(primary),
                            "genome_sha256": sha256_file(
                                shared / "primary_chromosomes.genome.fa"
                            ),
                        },
                        indent=2,
                        sort_keys=True,
                    )
                    + "\n",
                    encoding="utf-8",
                )

        blind = working / "blind"
        forbidden = ("target_complete", "evaluator_gardenia", "evaluator_ophiorrhiza")
        if any(
            any(token in path.as_posix() for token in forbidden)
            for path in blind.rglob("*")
        ):
            raise ValueError("Complete annotation or evaluator reference leaked blind")
        manifest = {
            "schema_version": "ploidypatch.coffea_external_normalized_inputs.v1.0",
            "holdout_id": HOLDOUT_ID,
            "policy_id": POLICY_ID,
            "role_layout": {
                "|".join(key): value for key, value in sorted(ROLE_LAYOUT.items())
            },
            "inputs": input_records,
            "firewall": {
                "blind_complete_target_annotation": False,
                "blind_evaluator_references": False,
                "blind_truth": False,
                "blind_labels": False,
                "blind_homoeolog_group_table": False,
            },
        }
        (working / "manifest.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        write_sha256sums(working)
        verify_sha256sums(working, ignore_checksum_file=True)
        os.replace(working, output)
        return 0
    except BaseException:
        shutil.rmtree(working, ignore_errors=True)
        raise


if __name__ == "__main__":
    raise SystemExit(main())
