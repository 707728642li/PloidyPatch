#!/usr/bin/env python3
"""Prepare only candidate-safe Walnut inputs inside the blind namespace."""
from __future__ import annotations

import csv
import json
import os
from pathlib import Path
import shutil
import sys
import tempfile

from ploidypatch.artifact_manifest import sha256_file, verify_sha256sums, write_sha256sums
from ploidypatch.gff_compat import synthesize_missing_transcript_exons
from ploidypatch.holdout_contract import load_holdout_contract
from ploidypatch.io import iter_fasta
from ploidypatch.normalize import prepare_primary_annotation_bundle, read_primary_seqid_table
from ploidypatch.walnut_h1 import (
    HOLDOUT_ID,
    POLICY_ID,
    write_primary_candidate_proteins,
)


CANDIDATES = {
    "Juglans_mandshurica": ("candidate_mandshurica", "config/primary_seqids/juglans_mandshurica_gwhbeun_v1.tsv"),
    "Carya_illinoinensis": ("candidate_carya", "config/primary_seqids/carya_illinoinensis_pawnee_v1.tsv"),
}


def read_roles(root: Path) -> dict[tuple[str, str, str], dict[str, str]]:
    with (root / "role_manifest.tsv").open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"role", "species_id", "artifact", "staged_relative_path", "bytes", "sha256"}
        if not required <= set(reader.fieldnames or ()):
            raise ValueError("Walnut blind role manifest fields differ")
        rows = list(reader)
    indexed = {(row["role"], row["species_id"], row["artifact"]): row for row in rows}
    if len(rows) != 15 or len(indexed) != 15:
        raise ValueError("Walnut role manifest must contain 15 unique artifacts")
    return indexed


def checked_artifact(root: Path, row: dict[str, str]) -> Path:
    relative = Path(row["staged_relative_path"])
    if relative.is_absolute() or ".." in relative.parts:
        raise ValueError(f"Unsafe staged Walnut path: {relative}")
    path = root / relative
    if (
        not path.is_file()
        or path.is_symlink()
        or path.stat().st_size != int(row["bytes"])
        or sha256_file(path) != row["sha256"]
    ):
        raise ValueError(f"Walnut blind artifact differs from role manifest: {path}")
    return path


def write_primary_genome(source: Path, table: Path, output: Path) -> list[str]:
    seqids, labels = read_primary_seqid_table(table)
    allowed = set(seqids)
    seen: set[str] = set()
    fai_rows: list[str] = []
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("xb") as handle:
        for seqid, header, sequence in iter_fasta(source):
            if seqid not in allowed:
                continue
            if seqid in seen or not sequence:
                raise ValueError(f"Duplicate or empty target primary seqid: {seqid}")
            seen.add(seqid)
            handle.write(f">{header}\n".encode("utf-8"))
            sequence_offset = handle.tell()
            encoded = sequence.encode("ascii")
            for offset in range(0, len(encoded), 60):
                handle.write(encoded[offset:offset + 60] + b"\n")
            line_bases = min(60, len(encoded))
            fai_rows.append(
                f"{seqid}\t{len(encoded)}\t{sequence_offset}\t{line_bases}\t{line_bases + 1}\n"
            )
    if seen != allowed:
        raise ValueError(f"Walnut target primary seqids differ: {sorted(allowed - seen)}")
    output.with_suffix(output.suffix + ".fai").write_text(
        "".join(fai_rows), encoding="utf-8"
    )
    return list(labels)


def main(argv: list[str] | None = None) -> int:
    values = sys.argv[1:] if argv is None else argv
    if len(values) != 1:
        raise SystemExit("usage: prepare_walnut_blind_candidate_inputs_v0.8.py PROJECT_ROOT")
    project_root = Path(values[0]).resolve()
    staged = Path(os.environ["PLOIDYPATCH_STAGED_INPUT_ROOT"]).resolve()
    protocol = Path(os.environ["PLOIDYPATCH_PROTOCOL_FREEZE"]).resolve()
    execution = Path(os.environ["PLOIDYPATCH_EXECUTION_FREEZE"]).resolve()
    contract_path = Path(os.environ["PLOIDYPATCH_HOLDOUT_CONTRACT"]).resolve()
    benchmark = Path(os.environ["PLOIDYPATCH_BLIND_BENCHMARK_ROOT"]).resolve()
    output = project_root / "results/baselines/walnut/v0.8/normalized"
    if output.exists() or output.is_symlink():
        raise FileExistsError(f"Refusing to overwrite Walnut blind normalized inputs: {output}")
    for root in (protocol, execution, benchmark):
        verify_sha256sums(root, ignore_checksum_file=True)
    contract = load_holdout_contract(contract_path)
    if contract.holdout_id != HOLDOUT_ID or contract.policy_id != POLICY_ID:
        raise ValueError("Walnut blind normalization received a different contract")
    if (staged / "evaluator_only").exists():
        raise ValueError("Evaluator-only role is visible to Walnut blind normalization")
    roles = read_roles(staged)
    output.parent.mkdir(parents=True, exist_ok=True)
    working = Path(tempfile.mkdtemp(prefix=".normalized.working.", dir=output.parent))
    try:
        target_genome = checked_artifact(
            staged, roles[("target", "Juglans_regia", "genome")]
        )
        target_table = execution / "source/config/primary_seqids/juglans_regia_walnut_2.0.tsv"
        target_output = working / "target_walnut2/primary_chromosomes.genome.fa"
        seqids = write_primary_genome(target_genome, target_table, target_output)
        blind_manifest = json.loads(
            (benchmark / "blind_manifest.json").read_text(encoding="utf-8")
        )
        expected_genome = blind_manifest.get("target_genome", {}).get("sha256")
        if expected_genome != sha256_file(target_output):
            raise ValueError("Walnut blind target genome differs from benchmark sentinel")

        inputs: dict[str, dict[str, str | int]] = {
            "target_genome": {
                "bytes": target_genome.stat().st_size,
                "sha256": sha256_file(target_genome),
            }
        }
        for species, (bundle_name, table_relative) in CANDIDATES.items():
            artifacts = {
                name: checked_artifact(
                    staged, roles[("candidate_reference", species, name)]
                )
                for name in ("genome", "gff3", "protein")
            }
            bundle = working / bundle_name
            prepare_primary_annotation_bundle(
                gff_path=artifacts["gff3"],
                genome_path=artifacts["genome"],
                primary_seqid_table_path=execution / "source" / table_relative,
                output_dir=bundle,
                canonical_fasta_headers=True,
            )
            lifton_gff = bundle / "primary_chromosomes.lifton.gff3"
            report = synthesize_missing_transcript_exons(
                bundle / "primary_chromosomes.gff3",
                lifton_gff,
                repair_parent_bounds=True,
            )
            counts = report.get("counts", {})
            if (
                report.get("child_coordinate_or_cds_changes") is not False
                or counts.get("unresolved_transcripts") != 0
                or counts.get("input_cds_records") != counts.get("output_cds_records")
            ):
                raise ValueError(f"Walnut LiftOn compatibility invariant failed: {species}")
            write_primary_candidate_proteins(
                primary_gff_path=bundle / "primary_chromosomes.gff3",
                provider_protein_path=artifacts["protein"],
                output_fasta_path=bundle / "provider.protein.fa",
                whitelist_tsv_path=bundle / "primary_protein_whitelist.tsv",
                manifest_path=bundle / "primary_protein_whitelist.manifest.json",
            )
            for artifact, path in artifacts.items():
                inputs[f"{bundle_name}:{artifact}"] = {
                    "bytes": path.stat().st_size,
                    "sha256": sha256_file(path),
                }
        manifest = {
            "schema_version": "ploidypatch.walnut_h1_blind_normalized_inputs.v0.8",
            "holdout_id": HOLDOUT_ID,
            "policy_id": POLICY_ID,
            "target_primary_seqids": seqids,
            "roles": ["shared_target_genome", "candidate_reference"],
            "inputs": inputs,
            "truth_access": False,
            "complete_target_annotation_access": False,
            "evaluator_reference_access": False,
            "ranker_access": False,
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
