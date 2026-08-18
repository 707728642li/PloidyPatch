#!/usr/bin/env python3
"""Seal the complete Walnut evaluator-only scientific input tree."""
from __future__ import annotations

import json
import os
from pathlib import Path
import sys

from ploidypatch.artifact_manifest import sha256_file, verify_sha256sums, write_sha256sums
from ploidypatch.walnut_h1 import HOLDOUT_ID, POLICY_ID, load_json_object


def main(argv: list[str] | None = None) -> int:
    values = sys.argv[1:] if argv is None else argv
    if len(values) != 1:
        raise SystemExit("usage: finalize_walnut_evaluator_only_inputs_v0.8.py PROJECT_ROOT")
    project_root = Path(values[0]).resolve()
    root = Path(os.environ["PLOIDYPATCH_EVALUATOR_ONLY_ROOT"]).resolve()
    if not project_root.is_dir() or not root.is_dir() or root.is_symlink():
        raise ValueError("Walnut evaluator-only root is missing or symlinked")
    if (root / "SHA256SUMS").exists() or (root / "manifest.json").exists():
        raise FileExistsError("Refusing to overwrite finalized Walnut evaluator-only root")
    required = {
        "target_gff": root / "normalized/target_complete/primary_chromosomes.gff3",
        "target_genome": root / "normalized/target_complete/primary_chromosomes.genome.fa",
        "corylus_gff": root / "normalized/truth_references/evaluator_corylus/primary_chromosomes.gff3",
        "castanea_gff": root / "normalized/truth_references/evaluator_castanea/primary_chromosomes.gff3",
        "wgdi_inputs": root / "wgdi/input/SHA256SUMS",
        "wgdi_evidence": root / "wgdi/evidence/SHA256SUMS",
        "truth_pairs": root / "truth_pairs/SHA256SUMS",
        "benchmark": root / "benchmark/SHA256SUMS",
        "evaluability": root / "benchmark/pair_selection/evaluability.json",
        "truth": root / "benchmark/truth/hidden_truth.json",
        "perturbed_copy": root / "benchmark/inputs/perturbed.gff3",
    }
    for label, path in required.items():
        if not path.is_file() or path.is_symlink() or path.stat().st_size == 0:
            raise ValueError(f"Walnut evaluator artifact is missing: {label}")
    for subroot in (root / "wgdi/input", root / "wgdi/evidence", root / "truth_pairs", root / "benchmark"):
        verify_sha256sums(subroot, ignore_checksum_file=True)
    evaluability = load_json_object(required["evaluability"])
    if (
        evaluability.get("holdout_id") != HOLDOUT_ID
        or evaluability.get("policy_id") != POLICY_ID
        or evaluability.get("status") not in {"ready", "not_evaluable", "invalid"}
    ):
        raise ValueError("Walnut evaluator evaluability record differs")
    manifest = {
        "schema_version": "ploidypatch.walnut_h1_evaluator_only_inputs.v0.8",
        "holdout_id": HOLDOUT_ID,
        "policy_id": POLICY_ID,
        "formal_status": evaluability["status"],
        "role": "evaluator_only_never_mounted_to_blind_runner",
        "truth_reveal_requires_sealed_blind_custody": True,
        "candidate_references_used_for_truth": False,
        "closed_yn00_ks": [0.10, 0.75],
        "evaluator_groups": ["Corylus_avellana", "Castanea_mollissima"],
        "exact_target_counterpart_multiplicity": 2,
        "minimum_block_pairs": 20,
        "inputs": {
            label: {"relative_path": path.relative_to(root).as_posix(), "sha256": sha256_file(path)}
            for label, path in required.items()
        },
        "ranker_or_model_present": False,
        "h2_or_topology_ranking_present": False,
        "automatic_approval": False,
    }
    (root / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    write_sha256sums(root)
    verify_sha256sums(root, ignore_checksum_file=True)
    return 1 if evaluability["status"] == "invalid" else 0


if __name__ == "__main__":
    raise SystemExit(main())
