from __future__ import annotations

import json
from pathlib import Path

from ploidypatch.core_h1_pools import (
    build_raw_predictions_manifest,
    collapse_duplicate_exact_chains,
)


def _model(gene: str, transcript: str, start: int) -> str:
    return (
        f"chr1\tx\tgene\t{start}\t{start+20}\t.\t+\t.\tID={gene}\n"
        f"chr1\tx\tmRNA\t{start}\t{start+20}\t.\t+\t.\tID={transcript};Parent={gene}\n"
        f"chr1\tx\tCDS\t{start}\t{start+20}\t.\t+\t0\tID={transcript}.cds;Parent={transcript}\n"
    )


def test_exact_chain_collapse_is_method_vote_not_reference_vote(tmp_path: Path) -> None:
    source = tmp_path / "candidate.gff3"
    source.write_text(
        "##gff-version 3\n" + _model("g1", "t1", 1) + _model("g2", "t2", 1),
        encoding="utf-8",
    )
    output = tmp_path / "deduplicated.gff3"
    audit = collapse_duplicate_exact_chains(candidate_gff=source, output=output)
    assert audit["input_transcripts"] == 2
    assert audit["retained_transcripts"] == 1
    assert audit["collapsed_duplicate_transcripts"] == 1
    assert "ID=t1" in output.read_text(encoding="utf-8")
    assert "ID=t2" not in output.read_text(encoding="utf-8")
    manifest = json.loads(
        Path(str(output) + ".manifest.json").read_text(encoding="utf-8")
    )
    assert manifest["truth_access"] is False
    assert manifest["ranker_access"] is False


def test_raw_prediction_manifest_requires_exact_three_by_reference_tree(
    tmp_path: Path,
) -> None:
    trees = {}
    for method in ("miniprot", "gemoma", "lifton"):
        for reference in ("bua", "mauritiana"):
            tree = tmp_path / method / reference
            tree.mkdir(parents=True)
            (tree / "artifact.txt").write_text(f"{method}:{reference}\n", encoding="utf-8")
            trees[f"{method}__{reference}"] = tree
    output = tmp_path / "raw.json"
    manifest = build_raw_predictions_manifest(
        holdout_id="coffea",
        policy_id="policy",
        project_root=tmp_path,
        raw_prediction_trees=trees,
        candidate_references=("bua", "mauritiana"),
        input_hashes={"input": "a" * 64},
        output=output,
    )
    assert len(manifest["raw_prediction_trees"]) == 6
    assert manifest["within_method_reference_vote_count"] == 1
    assert manifest["truth_access"] is False
    assert manifest["ranker_access"] is False

