#!/usr/bin/env python3
"""Build the five metadata-only Coffea protein-supported WGDI universes."""
from __future__ import annotations

import json
import os
from pathlib import Path
import sys
import tempfile

from ploidypatch.artifact_manifest import verify_sha256sums, write_sha256sums
from ploidypatch.known_subgenome_h1 import HOLDOUT_ID, POLICY_ID
from ploidypatch.protein_universe import build_protein_supported_gene_universe


BUNDLES = {
    "Coffea_arabica_ET39": "evaluator_only/normalized/target_complete",
    "Coffea_eugenioides_BuA": "blind/candidate_only/candidate_bua",
    "Coffea_mauritiana": "blind/candidate_only/candidate_mauritiana",
    "Gardenia_jasminoides": (
        "evaluator_only/normalized/truth_references/evaluator_gardenia"
    ),
    "Ophiorrhiza_pumila": (
        "evaluator_only/normalized/truth_references/evaluator_ophiorrhiza"
    ),
}


def main(argv: list[str] | None = None) -> int:
    values = sys.argv[1:] if argv is None else argv
    if len(values) != 2:
        raise SystemExit(
            "usage: build_coffea_protein_supported_universes_v1.0.py "
            "NORMALIZED_ROOT OUTPUT_ROOT"
        )
    normalized = Path(values[0]).resolve()
    output = Path(values[1]).resolve()
    verify_sha256sums(normalized, ignore_checksum_file=True)
    if output.exists() or output.is_symlink():
        raise FileExistsError(f"Refusing to overwrite Coffea protein universes: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    working = Path(tempfile.mkdtemp(prefix=f".{output.name}.working.", dir=output.parent))
    try:
        manifests = {}
        for species_id, relative in sorted(BUNDLES.items()):
            bundle = normalized / relative
            manifests[species_id] = build_protein_supported_gene_universe(
                primary_gff_path=bundle / "primary_chromosomes.gff3",
                provider_protein_path=bundle / "provider.protein.fa",
                output_dir=working / species_id,
                species_id=species_id,
                holdout_id=HOLDOUT_ID,
                policy_id=POLICY_ID,
            )
        (working / "manifest.json").write_text(
            json.dumps(
                {
                    "schema_version": "ploidypatch.coffea_protein_universes.v1.0",
                    "holdout_id": HOLDOUT_ID,
                    "policy_id": POLICY_ID,
                    "labels_used": False,
                    "truth_pairs_used": False,
                    "candidate_predictions_used": False,
                    "species": {
                        species_id: {
                            "coding_genes": manifest["counts"]["coding_genes"],
                            "genes_with_exact_provider_protein": manifest["counts"][
                                "genes_with_exact_provider_protein"
                            ],
                            "genes_excluded_without_exact_provider_protein": manifest[
                                "counts"
                            ]["genes_excluded_without_exact_provider_protein"],
                        }
                        for species_id, manifest in sorted(manifests.items())
                    },
                },
                indent=2,
                sort_keys=True,
            )
            + "\n",
            encoding="utf-8",
        )
        write_sha256sums(working)
        verify_sha256sums(working, ignore_checksum_file=True)
        os.replace(working, output)
        return 0
    except BaseException:
        # A working tree is intentionally retained for a diagnosable failed run.
        raise


if __name__ == "__main__":
    raise SystemExit(main())

