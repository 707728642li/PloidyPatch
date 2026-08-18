#!/usr/bin/env python3
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
from ploidypatch.io import fasta_relation_id, iter_fasta
from ploidypatch.synteny_io import prepare_wgdi_inputs
from ploidypatch.walnut_h1 import (
    HOLDOUT_ID,
    POLICY_ID,
    write_primary_candidate_proteins,
)


BUNDLES = {
    "target": ("normalized/target_complete", "jre"),
    "corylus": ("normalized/truth_references/evaluator_corylus", "cav"),
    "castanea": ("normalized/truth_references/evaluator_castanea", "cmo"),
}


def build_target_gene_cds(
    *, representatives: Path, provider_protein: Path, transcript_cds: Path, output: Path
) -> None:
    protein_to_transcript: dict[str, str] = {}
    for protein_id, header, _ in iter_fasta(provider_protein):
        transcript_id, _ = fasta_relation_id(protein_id, header, relation="transcript")
        protein_to_transcript[protein_id] = transcript_id
    cds = {identifier: sequence for identifier, _, sequence in iter_fasta(transcript_cds)}
    rows: list[tuple[str, str]] = []
    with representatives.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            transcript = protein_to_transcript.get(row["protein_id"])
            if not transcript or transcript not in cds:
                raise ValueError(
                    f"Target representative lacks exact transcript CDS: {row['gene_id']}"
                )
            rows.append((row["gene_id"], cds[transcript]))
    if not rows:
        raise ValueError("Target representative CDS set is empty")
    with output.open("x", encoding="utf-8", newline="") as handle:
        for gene, sequence in rows:
            handle.write(f">{gene}\n")
            for index in range(0, len(sequence), 60):
                handle.write(sequence[index:index + 60] + "\n")


def main(argv: list[str] | None = None) -> int:
    values = sys.argv[1:] if argv is None else argv
    if len(values) != 1:
        raise SystemExit("usage: prepare_walnut_evaluator_wgdi_inputs_v0.8.py PROJECT_ROOT")
    project_root = Path(values[0]).resolve()
    evaluator_root = Path(os.environ["PLOIDYPATCH_EVALUATOR_ONLY_ROOT"]).resolve()
    execution = Path(os.environ["PLOIDYPATCH_EXECUTION_FREEZE"]).resolve()
    output = evaluator_root / "wgdi/input"
    gffread = project_root / "envs/ploidypatch-syngap/bin/gffread"
    alias_builder = execution / "source/scripts/build_wgdi_source_alias_gff.py"
    for path in (project_root, evaluator_root, execution):
        if not path.is_dir() or path.is_symlink():
            raise ValueError(f"Walnut WGDI root is missing or symlinked: {path}")
    verify_sha256sums(execution, ignore_checksum_file=True)
    for executable in (gffread, alias_builder):
        if not executable.is_file() or executable.is_symlink():
            raise ValueError(f"Walnut WGDI executable/helper is missing: {executable}")
    if output.exists() or output.is_symlink():
        raise FileExistsError(f"Refusing to overwrite Walnut WGDI inputs: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    working = Path(tempfile.mkdtemp(prefix=".input.working.", dir=output.parent))
    try:
        for label, (relative, prefix) in BUNDLES.items():
            bundle = evaluator_root / relative
            gff = bundle / "primary_chromosomes.gff3"
            genome = bundle / "primary_chromosomes.genome.fa"
            fai = bundle / "primary_chromosomes.genome.fa.fai"
            protein = bundle / "provider.protein.fa"
            for path in (gff, genome, fai, protein):
                if not path.is_file() or path.is_symlink() or path.stat().st_size == 0:
                    raise ValueError(f"Walnut normalized evaluator artifact is missing: {path}")
            destination = working / label
            destination.mkdir()
            wgdi_protein = protein
            if label != "target":
                wgdi_protein = destination / f"{prefix}.primary.provider.protein.fa"
                write_primary_candidate_proteins(
                    primary_gff_path=gff,
                    provider_protein_path=protein,
                    output_fasta_path=wgdi_protein,
                    whitelist_tsv_path=destination / f"{prefix}.primary.protein_whitelist.tsv",
                    manifest_path=destination / f"{prefix}.primary.protein_manifest.json",
                )
            prepare_wgdi_inputs(
                gff_path=gff,
                protein_path=wgdi_protein,
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
                        "--source-gff", str(gff),
                        "--representatives", str(destination / f"{prefix}.representatives.tsv"),
                        "--output-gff", str(destination / f"{prefix}.source_alias.gff3"),
                    ],
                    check=True,
                )
                transcript_cds = destination / "target.transcript.cds.fa"
                subprocess.run(
                    [str(gffread), str(gff), "-g", str(genome), "-x", str(transcript_cds)],
                    check=True,
                )
                build_target_gene_cds(
                    representatives=destination / f"{prefix}.representatives.tsv",
                    provider_protein=protein,
                    transcript_cds=transcript_cds,
                    output=destination / f"{prefix}.wgdi.cds.fa",
                )
        manifest = {
            "schema_version": "ploidypatch.walnut_evaluator_wgdi_inputs.v0.8",
            "holdout_id": HOLDOUT_ID,
            "policy_id": POLICY_ID,
            "role": "evaluator_only",
            "candidate_references_used": False,
            "bundles": {
                label: {
                    "prefix": prefix,
                    "wgdi_manifest_sha256": sha256_file(
                        working / label / f"{prefix}.wgdi_inputs.manifest.json"
                    ),
                    "primary_protein_manifest_sha256": (
                        sha256_file(
                            working / label / f"{prefix}.primary.protein_manifest.json"
                        )
                        if label != "target"
                        else None
                    ),
                }
                for label, (_, prefix) in BUNDLES.items()
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
