#!/usr/bin/env python3
"""Prepare only candidate-safe Coffea inputs inside the blind namespace."""
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
from ploidypatch.known_subgenome_h1 import HOLDOUT_ID, POLICY_ID
from ploidypatch.normalize import prepare_primary_annotation_bundle, read_primary_seqid_table


CANDIDATES = {
    "Coffea_eugenioides_BuA": (
        "candidate_bua",
        "config/primary_seqids/coffea_eugenioides_bua_v1.tsv",
    ),
    "Coffea_mauritiana": (
        "candidate_mauritiana",
        "config/primary_seqids/coffea_mauritiana_v1.tsv",
    ),
}


def read_roles(root: Path) -> dict[tuple[str, str, str], dict[str, str]]:
    with (root / "role_manifest.tsv").open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "role", "species_id", "artifact", "staged_relative_path", "bytes", "sha256"
        }
        if not required <= set(reader.fieldnames or ()):
            raise ValueError("Coffea blind role manifest fields differ")
        rows = list(reader)
    indexed = {(row["role"], row["species_id"], row["artifact"]): row for row in rows}
    expected = {
        ("target", "Coffea_arabica_ET39", "genome"),
        *{
            ("candidate_reference", species, artifact)
            for species in CANDIDATES
            for artifact in ("genome", "gff3", "protein")
        },
    }
    if set(indexed) != expected or len(rows) != len(expected):
        raise ValueError("Coffea blind role manifest must contain exact target+candidate artifacts")
    return indexed


def checked_artifact(root: Path, row: dict[str, str]) -> Path:
    relative = Path(row["staged_relative_path"])
    if relative.is_absolute() or ".." in relative.parts:
        raise ValueError(f"Unsafe staged Coffea path: {relative}")
    path = root / relative
    if (
        not path.is_file()
        or path.is_symlink()
        or path.stat().st_size != int(row["bytes"])
        or sha256_file(path) != row["sha256"]
    ):
        raise ValueError(f"Coffea blind artifact differs from role manifest: {path}")
    return path


def write_primary_genome(source: Path, table: Path, output: Path) -> list[str]:
    seqids, labels = read_primary_seqid_table(table)
    allowed = set(seqids)
    seen: set[str] = set()
    fai_rows: list[str] = []
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("xb") as handle:
        for seqid, _header, sequence in iter_fasta(source):
            if seqid not in allowed:
                continue
            if seqid in seen or not sequence:
                raise ValueError(f"Duplicate or empty Coffea target primary seqid: {seqid}")
            seen.add(seqid)
            handle.write(f">{seqid}\n".encode("utf-8"))
            sequence_offset = handle.tell()
            encoded = sequence.encode("ascii")
            for offset in range(0, len(encoded), 60):
                handle.write(encoded[offset : offset + 60] + b"\n")
            line_bases = min(60, len(encoded))
            fai_rows.append(
                f"{seqid}\t{len(encoded)}\t{sequence_offset}\t{line_bases}\t{line_bases + 1}\n"
            )
    if seen != allowed:
        raise ValueError(f"Coffea target primary seqids differ: {sorted(allowed - seen)}")
    output.with_suffix(output.suffix + ".fai").write_text(
        "".join(fai_rows), encoding="utf-8"
    )
    return list(labels)


def main(argv: list[str] | None = None) -> int:
    values = sys.argv[1:] if argv is None else argv
    if len(values) != 1:
        raise SystemExit("usage: prepare_coffea_blind_candidate_inputs_v1.0.py PROJECT_ROOT")
    project_root = Path(values[0]).resolve()
    staged = Path(os.environ["PLOIDYPATCH_STAGED_INPUT_ROOT"]).resolve()
    protocol = Path(os.environ["PLOIDYPATCH_PROTOCOL_FREEZE"]).resolve()
    execution = Path(os.environ["PLOIDYPATCH_EXECUTION_FREEZE"]).resolve()
    contract_path = Path(os.environ["PLOIDYPATCH_HOLDOUT_CONTRACT"]).resolve()
    benchmark = Path(os.environ["PLOIDYPATCH_BLIND_BENCHMARK_ROOT"]).resolve()
    output = project_root / "results/baselines/coffea/v1.0/normalized"
    if output.exists() or output.is_symlink():
        raise FileExistsError(f"Refusing to overwrite Coffea blind normalized inputs: {output}")
    for root in (staged, protocol, execution, benchmark):
        verify_sha256sums(root, ignore_checksum_file=True)
    contract = load_holdout_contract(contract_path)
    if contract.holdout_id != HOLDOUT_ID or contract.policy_id != POLICY_ID:
        raise ValueError("Coffea blind normalization received a different contract")
    forbidden_names = {"evaluator_only", "target_complete", "truth", "labels"}
    if any(
        set(part.casefold() for part in path.relative_to(staged).parts) & forbidden_names
        for path in staged.rglob("*")
    ):
        raise ValueError("Evaluator/truth role is visible to Coffea blind normalization")
    roles = read_roles(staged)
    output.parent.mkdir(parents=True, exist_ok=True)
    working = Path(tempfile.mkdtemp(prefix=".normalized.working.", dir=output.parent))
    try:
        target_genome = checked_artifact(
            staged, roles[("target", "Coffea_arabica_ET39", "genome")]
        )
        target_table = execution / "source/config/primary_seqids/coffea_arabica_et39_hifi.tsv"
        target_output = working / "target_et39/primary_chromosomes.genome.fa"
        seqids = write_primary_genome(target_genome, target_table, target_output)
        blind_manifest = json.loads(
            (benchmark / "blind_manifest.json").read_text(encoding="utf-8")
        )
        if blind_manifest.get("target_genome", {}).get("sha256") != sha256_file(
            target_output
        ):
            raise ValueError("Coffea blind target genome differs from benchmark sentinel")
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
            universe = staged / "candidate_only/protein_universes" / species
            verify_sha256sums(universe, ignore_checksum_file=True)
            representative = universe / "representative.protein.fa"
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
                raise ValueError(f"Coffea LiftOn compatibility invariant failed: {species}")
            shutil.copyfile(representative, bundle / "provider.protein.fa")
            if sha256_file(bundle / "provider.protein.fa") != sha256_file(representative):
                raise ValueError(f"Coffea candidate protein copy differs: {species}")
            for artifact, path in artifacts.items():
                inputs[f"{bundle_name}:{artifact}"] = {
                    "bytes": path.stat().st_size,
                    "sha256": sha256_file(path),
                }
        manifest = {
            "schema_version": "ploidypatch.coffea_h1_blind_normalized_inputs.v1.0",
            "holdout_id": HOLDOUT_ID,
            "policy_id": POLICY_ID,
            "target_primary_seqids": seqids,
            "roles": ["shared_target_genome", "candidate_reference"],
            "candidate_protein_policy": "exact_provider_supported_subset",
            "truth_access": False,
            "complete_target_annotation_access": False,
            "evaluator_reference_access": False,
            "ranker_access": False,
            "inputs": inputs,
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
