#!/usr/bin/env python3
from __future__ import annotations

import json
import os
from pathlib import Path
import shutil
import sys
import tempfile

from ploidypatch.artifact_manifest import sha256_file, verify_sha256sums, write_sha256sums
from ploidypatch.homeolog_pairs import infer_outgroup_duplicated_pairs
from ploidypatch.self_wgd_pairs import infer_self_wgdi_pairs
from ploidypatch.walnut_h1 import (
    HOLDOUT_ID,
    POLICY_ID,
    TARGET_EVENT,
    filter_self_pairs_by_closed_yn00_ks,
    infer_exact_two_outgroup_pair_consistent_truth,
)


def main(argv: list[str] | None = None) -> int:
    values = sys.argv[1:] if argv is None else argv
    if len(values) != 1:
        raise SystemExit("usage: infer_walnut_external_pairs_v0.8.py PROJECT_ROOT")
    project_root = Path(values[0]).resolve()
    evaluator_root = Path(os.environ["PLOIDYPATCH_EVALUATOR_ONLY_ROOT"]).resolve()
    wgdi = evaluator_root / "wgdi"
    evidence = wgdi / "evidence"
    output = evaluator_root / "truth_pairs"
    query_gff = wgdi / "input/target/jre.wgdi.gff"
    source_alias = wgdi / "input/target/jre.source_alias.gff3"
    self_collinearity = evidence / "collinearity/jre_self.tsv"
    corylus_collinearity = evidence / "collinearity/jre_vs_cav.tsv"
    castanea_collinearity = evidence / "collinearity/jre_vs_cmo.tsv"
    ks = evidence / "ks/target/ks_merged.tsv"
    for path in (
        project_root, evaluator_root, query_gff, source_alias, self_collinearity,
        corylus_collinearity, castanea_collinearity, ks,
    ):
        if not path.exists() or path.is_symlink():
            raise ValueError(f"Walnut pair prerequisite is missing or symlinked: {path}")
    if output.exists() or output.is_symlink():
        raise FileExistsError(f"Refusing to overwrite Walnut truth pairs: {output}")
    working = Path(tempfile.mkdtemp(prefix=".truth_pairs.working.", dir=evaluator_root))
    try:
        for name in ("self", "ks", "corylus", "castanea"):
            (working / name).mkdir()
        infer_self_wgdi_pairs(
            query_wgdi_gff_path=query_gff,
            collinearity_path=self_collinearity,
            source_gff_path=source_alias,
            output_pair_tsv_path=working / "self/pairs.tsv",
            decisions_tsv_path=working / "self/decisions.tsv",
            wgd_event="Juglandaceae_ancestral_juglandoid_WGD_self",
            min_block_pairs=20,
            require_different_seqids=True,
            require_reciprocal_unique=True,
        )
        filter_self_pairs_by_closed_yn00_ks(
            self_pairs_path=working / "self/pairs.tsv",
            ks_path=ks,
            output_pairs_path=working / "ks/pairs.tsv",
            decisions_path=working / "ks/decisions.tsv",
        )
        for label, collinearity in (
            ("corylus", corylus_collinearity),
            ("castanea", castanea_collinearity),
        ):
            infer_outgroup_duplicated_pairs(
                query_wgdi_gff_path=query_gff,
                source_gff_path=source_alias,
                collinearity_inputs=(f"{label}={collinearity}",),
                output_pair_tsv_path=working / label / "pairs.tsv",
                decisions_tsv_path=working / label / "decisions.tsv",
                wgd_event=f"Juglandaceae_ancestral_juglandoid_WGD_{label}",
                min_support_group_count=1,
                min_block_pairs=20,
                require_cross_seqid=True,
                require_reciprocal_unique=True,
            )
        final_manifest = infer_exact_two_outgroup_pair_consistent_truth(
            ks_filtered_self_pairs_path=working / "ks/pairs.tsv",
            evaluator_group_decisions={
                "corylus": working / "corylus/decisions.tsv",
                "castanea": working / "castanea/decisions.tsv",
            },
            output_dir=working / "final",
        )
        shutil.copyfile(working / "final/pairs.tsv", working / "pairs.tsv")
        shutil.copyfile(working / "final/decisions.tsv", working / "decisions.tsv")
        manifest = {
            "schema_version": "ploidypatch.walnut_truth_pair_bundle.v0.8",
            "holdout_id": HOLDOUT_ID,
            "policy_id": POLICY_ID,
            "target_event": TARGET_EVENT,
            "formal_rule": final_manifest["formal_rule"],
            "candidate_references_used_for_truth": False,
            "evaluator_groups": ["Corylus_avellana", "Castanea_mollissima"],
            "parameters": {
                "minimum_block_pairs": 20,
                "cross_primary_chromosome": True,
                "reciprocal_unique": True,
                "closed_yn00_ks": [0.10, 0.75],
                "exact_counterpart_target_multiplicity": 2,
                "post_enumeration_relaxation": False,
            },
            "counts": final_manifest["counts"],
            "outputs": {
                "pairs_sha256": sha256_file(working / "pairs.tsv"),
                "decisions_sha256": sha256_file(working / "decisions.tsv"),
            },
            "truth_labels_accessed": False,
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
