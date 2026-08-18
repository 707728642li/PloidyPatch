#!/usr/bin/env python3
"""Prepare role-isolated Coffea target and evaluator-only WGDI inputs."""
from __future__ import annotations

import csv
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

from ploidypatch.artifact_manifest import sha256_file, verify_sha256sums, write_sha256sums
from ploidypatch.io import iter_fasta, normalize_feature_id
from ploidypatch.known_subgenome_h1 import HOLDOUT_ID, POLICY_ID
from ploidypatch.perturb import read_gff_document
from ploidypatch.synteny_io import prepare_wgdi_inputs


BUNDLES = {
    "target": ("normalized/target_complete", "Coffea_arabica_ET39", "car"),
    "gardenia": (
        "normalized/truth_references/evaluator_gardenia",
        "Gardenia_jasminoides",
        "gja",
    ),
    "ophiorrhiza": (
        "normalized/truth_references/evaluator_ophiorrhiza",
        "Ophiorrhiza_pumila",
        "opu",
    ),
}


def build_target_gene_cds(
    *, supported_gff: Path, representatives: Path, transcript_cds: Path, output: Path
) -> None:
    document = read_gff_document(supported_gff)
    protein_transcripts: dict[str, set[str]] = {}
    for record in document.records:
        if record.feature_type != "CDS":
            continue
        protein_id = record.attributes.get("protein_id")
        if not protein_id:
            continue
        parents = {normalize_feature_id(parent) for parent in record.parents}
        if len(parents) != 1:
            raise ValueError("Target CDS protein relation is not single-transcript")
        protein_transcripts.setdefault(normalize_feature_id(protein_id), set()).update(
            parents
        )
    cds = {identifier: sequence for identifier, _, sequence in iter_fasta(transcript_cds)}
    rows: list[tuple[str, str]] = []
    with representatives.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not {"gene_id", "protein_id"} <= set(reader.fieldnames or ()):
            raise ValueError("Target representatives lack gene_id/protein_id")
        for row in reader:
            transcripts = protein_transcripts.get(normalize_feature_id(row["protein_id"]), set())
            if len(transcripts) != 1:
                raise ValueError(
                    f"Target representative protein lacks one transcript: {row['gene_id']}"
                )
            transcript = next(iter(transcripts))
            if transcript not in cds:
                raise ValueError(f"Target representative lacks exact CDS: {row['gene_id']}")
            rows.append((row["gene_id"], cds[transcript]))
    if not rows or len(rows) != len({gene for gene, _ in rows}):
        raise ValueError("Target representative CDS universe is empty or duplicated")
    with output.open("x", encoding="utf-8", newline="") as handle:
        for gene, sequence in rows:
            handle.write(f">{gene}\n")
            for offset in range(0, len(sequence), 60):
                handle.write(sequence[offset : offset + 60] + "\n")


def main(argv: list[str] | None = None) -> int:
    values = sys.argv[1:] if argv is None else argv
    if len(values) != 1:
        raise SystemExit("usage: prepare_coffea_evaluator_wgdi_inputs_v1.0.py PROJECT_ROOT")
    project_root = Path(values[0]).resolve()
    input_root = Path(os.environ["PLOIDYPATCH_EVALUATOR_INPUT_ROOT"]).resolve()
    evaluator_root = Path(os.environ["PLOIDYPATCH_EVALUATOR_ONLY_ROOT"]).resolve()
    execution = Path(os.environ["PLOIDYPATCH_EXECUTION_FREEZE"]).resolve()
    output = evaluator_root / "wgdi/input"
    gffread = project_root / "envs/ploidypatch-syngap/bin/gffread"
    alias_builder = execution / "source/scripts/build_wgdi_source_alias_gff.py"
    for path in (project_root, input_root, evaluator_root, execution):
        if not path.is_dir() or path.is_symlink():
            raise ValueError(f"Coffea evaluator WGDI root is missing or symlinked: {path}")
    verify_sha256sums(execution, ignore_checksum_file=True)
    verify_sha256sums(input_root, ignore_checksum_file=True)
    for forbidden in ("candidate_bua", "candidate_mauritiana", "candidate_only"):
        if any(forbidden in path.as_posix() for path in input_root.rglob("*")):
            raise ValueError(f"Candidate-only role leaked to Coffea evaluator: {forbidden}")
    for executable in (gffread, alias_builder):
        if not executable.is_file() or executable.is_symlink():
            raise ValueError(f"Coffea WGDI executable/helper is missing: {executable}")
    if output.exists() or output.is_symlink():
        raise FileExistsError(f"Refusing to overwrite Coffea evaluator WGDI inputs: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    working = Path(tempfile.mkdtemp(prefix=".input.working.", dir=output.parent))
    try:
        for label, (relative, species, prefix) in BUNDLES.items():
            bundle = input_root / relative
            universe = input_root / "protein_universes" / species
            genome = bundle / "primary_chromosomes.genome.fa"
            fai = bundle / "primary_chromosomes.genome.fa.fai"
            gff = universe / "protein_supported.gff3"
            protein = universe / "representative.protein.fa"
            for path in (genome, fai, gff, protein, universe / "manifest.json"):
                if not path.is_file() or path.is_symlink() or path.stat().st_size == 0:
                    raise ValueError(f"Coffea evaluator WGDI artifact is missing: {path}")
            verify_sha256sums(universe, ignore_checksum_file=True)
            destination = working / label
            destination.mkdir()
            prepare_wgdi_inputs(
                gff_path=gff,
                protein_path=protein,
                fai_path=fai,
                output_dir=destination,
                prefix=prefix,
                min_genes_per_seqid=100,
            )
            if label == "target":
                subprocess.run(
                    [
                        sys.executable,
                        str(alias_builder),
                        "--source-gff",
                        str(gff),
                        "--representatives",
                        str(destination / f"{prefix}.representatives.tsv"),
                        "--output-gff",
                        str(destination / f"{prefix}.source_alias.gff3"),
                    ],
                    check=True,
                )
                transcript_cds = destination / "target.transcript.cds.fa"
                subprocess.run(
                    [str(gffread), str(gff), "-g", str(genome), "-x", str(transcript_cds)],
                    check=True,
                )
                build_target_gene_cds(
                    supported_gff=gff,
                    representatives=destination / f"{prefix}.representatives.tsv",
                    transcript_cds=transcript_cds,
                    output=destination / f"{prefix}.wgdi.cds.fa",
                )
        manifest = {
            "schema_version": "ploidypatch.coffea_evaluator_wgdi_inputs.v1.0",
            "holdout_id": HOLDOUT_ID,
            "policy_id": POLICY_ID,
            "role": "evaluator_only",
            "candidate_references_used": False,
            "protein_universe_policy": "exact_provider_supported_subset",
            "bundles": {
                label: {
                    "species_id": species,
                    "prefix": prefix,
                    "wgdi_manifest_sha256": sha256_file(
                        working / label / f"{prefix}.wgdi_inputs.manifest.json"
                    ),
                }
                for label, (_relative, species, prefix) in BUNDLES.items()
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
