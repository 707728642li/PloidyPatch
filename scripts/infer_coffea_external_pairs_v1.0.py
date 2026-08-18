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
from ploidypatch.known_subgenome_h1 import (
    HOLDOUT_ID,
    POLICY_ID,
    TARGET_EVENT,
    attach_descriptive_yn00_ks,
    filter_self_pairs_by_known_subgenome,
    infer_exact_two_outgroup_pair_consistent_truth,
)
from ploidypatch.self_wgd_pairs import infer_self_wgdi_pairs


def main(argv: list[str] | None = None) -> int:
    values = sys.argv[1:] if argv is None else argv
    if len(values) != 1:
        raise SystemExit("usage: infer_coffea_external_pairs_v1.0.py PROJECT_ROOT")
    project_root = Path(values[0]).resolve()
    evaluator_root = Path(os.environ["PLOIDYPATCH_EVALUATOR_ONLY_ROOT"]).resolve()
    execution = Path(os.environ["PLOIDYPATCH_EXECUTION_FREEZE"]).resolve()
    wgdi = evaluator_root / "wgdi"
    evidence = wgdi / "evidence"
    output = evaluator_root / "truth_pairs"
    query_gff = wgdi / "input/target/car.wgdi.gff"
    source_alias = wgdi / "input/target/car.source_alias.gff3"
    self_collinearity = evidence / "collinearity/car_self.tsv"
    gardenia_collinearity = evidence / "collinearity/car_vs_gja.tsv"
    ophiorrhiza_collinearity = evidence / "collinearity/car_vs_opu.tsv"
    ks = evidence / "ks/target/ks_merged.tsv"
    groups = execution / "source/config/coffea_et39_homoeolog_groups_v1.0.tsv"
    for path in (
        project_root,
        evaluator_root,
        execution / "SHA256SUMS",
        query_gff,
        source_alias,
        self_collinearity,
        gardenia_collinearity,
        ophiorrhiza_collinearity,
        ks,
        groups,
    ):
        if not path.exists() or path.is_symlink():
            raise ValueError(f"Coffea pair prerequisite is missing or symlinked: {path}")
    verify_sha256sums(execution, ignore_checksum_file=True)
    if output.exists() or output.is_symlink():
        raise FileExistsError(f"Refusing to overwrite Coffea truth pairs: {output}")
    working = Path(tempfile.mkdtemp(prefix=".truth_pairs.working.", dir=evaluator_root))
    try:
        for name in ("self", "known_subgenome", "gardenia", "ophiorrhiza"):
            (working / name).mkdir()
        infer_self_wgdi_pairs(
            query_wgdi_gff_path=query_gff,
            collinearity_path=self_collinearity,
            source_gff_path=source_alias,
            output_pair_tsv_path=working / "self/pairs.tsv",
            decisions_tsv_path=working / "self/decisions.tsv",
            wgd_event=f"{TARGET_EVENT}_target_self",
            min_block_pairs=20,
            require_different_seqids=True,
            require_reciprocal_unique=True,
        )
        attach_descriptive_yn00_ks(
            self_pairs_path=working / "self/pairs.tsv",
            ks_path=ks,
            output_path=working / "self/pairs_with_descriptive_ks.tsv",
        )
        filter_self_pairs_by_known_subgenome(
            self_pairs_path=working / "self/pairs_with_descriptive_ks.tsv",
            homoeolog_groups_path=groups,
            output_dir=working / "known_subgenome/filter",
        )
        for label, collinearity in (
            ("gardenia", gardenia_collinearity),
            ("ophiorrhiza", ophiorrhiza_collinearity),
        ):
            infer_outgroup_duplicated_pairs(
                query_wgdi_gff_path=query_gff,
                source_gff_path=source_alias,
                collinearity_inputs=(f"{label}={collinearity}",),
                output_pair_tsv_path=working / label / "pairs.tsv",
                decisions_tsv_path=working / label / "decisions.tsv",
                wgd_event=f"{TARGET_EVENT}_{label}",
                min_support_group_count=1,
                min_block_pairs=20,
                require_cross_seqid=True,
                require_reciprocal_unique=True,
            )
        final_manifest = infer_exact_two_outgroup_pair_consistent_truth(
            known_subgenome_self_pairs_path=working / "known_subgenome/filter/pairs.tsv",
            evaluator_group_decisions={
                "gardenia": working / "gardenia/decisions.tsv",
                "ophiorrhiza": working / "ophiorrhiza/decisions.tsv",
            },
            output_dir=working / "final",
        )
        shutil.copyfile(working / "final/pairs.tsv", working / "pairs.tsv")
        shutil.copyfile(working / "final/decisions.tsv", working / "decisions.tsv")
        manifest = {
            "schema_version": "ploidypatch.coffea_truth_pair_bundle.v1.0",
            "holdout_id": HOLDOUT_ID,
            "policy_id": POLICY_ID,
            "target_event": TARGET_EVENT,
            "formal_rule": final_manifest["formal_rule"],
            "candidate_references_used_for_truth": False,
            "evaluator_groups": ["Gardenia_jasminoides", "Ophiorrhiza_pumila"],
            "parameters": {
                "minimum_block_pairs": 20,
                "cross_primary_chromosome": True,
                "reciprocal_unique": True,
                "predeclared_same_group_exact_C_E": True,
                "yn00_ks_selection_use": False,
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

