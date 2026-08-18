from __future__ import annotations

import csv
from pathlib import Path

import pytest

from ploidypatch.pair_consensus import intersect_copy_pair_evidence


def _write_pairs(path: Path, rows: list[tuple[str, str, str]]) -> Path:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("pair_id", "gene_id_a", "gene_id_b"))
        writer.writerows(rows)
    return path


def test_intersect_copy_pair_evidence_is_exact_unordered_and_audited(
    tmp_path: Path,
) -> None:
    self_pairs = _write_pairs(
        tmp_path / "self.tsv",
        [("s1", "A", "B"), ("s2", "C", "D"), ("s3", "E", "F")],
    )
    outgroup_pairs = _write_pairs(
        tmp_path / "outgroup.tsv",
        [("o1", "B", "A"), ("o2", "D", "C"), ("o3", "E", "G")],
    )
    output = tmp_path / "accepted.tsv"
    decisions = tmp_path / "decisions.tsv"

    manifest = intersect_copy_pair_evidence(
        pair_inputs=[f"self={self_pairs}", f"outgroups={outgroup_pairs}"],
        output_pair_tsv_path=output,
        decisions_tsv_path=decisions,
        pair_set_label="apple_wgd",
    )

    with output.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [(row["gene_id_a"], row["gene_id_b"]) for row in rows] == [
        ("A", "B"),
        ("C", "D"),
    ]
    assert manifest["counts"] == {
        "union_pairs": 4,
        "exact_intersection_pairs_before_reciprocal_gate": 2,
        "accepted_pairs": 2,
        "decision_reason_counts": {
            "all_required_pair_evidence_exact_match": 2,
            "missing_required_evidence": 2,
        },
    }
    assert manifest["access"] == "evaluator_only"

    with pytest.raises(FileExistsError, match="Refusing to overwrite"):
        intersect_copy_pair_evidence(
            pair_inputs=[f"self={self_pairs}", f"outgroups={outgroup_pairs}"],
            output_pair_tsv_path=output,
            decisions_tsv_path=decisions,
            pair_set_label="apple_wgd",
        )


def test_intersection_rejects_nonreciprocal_partner_sets(tmp_path: Path) -> None:
    left = _write_pairs(
        tmp_path / "left.tsv", [("l1", "A", "B"), ("l2", "A", "C")]
    )
    right = _write_pairs(
        tmp_path / "right.tsv", [("r1", "A", "B"), ("r2", "A", "C")]
    )
    manifest = intersect_copy_pair_evidence(
        pair_inputs=[f"left={left}", f"right={right}"],
        output_pair_tsv_path=tmp_path / "accepted.tsv",
        decisions_tsv_path=tmp_path / "decisions.tsv",
        pair_set_label="test",
    )
    assert manifest["counts"]["accepted_pairs"] == 0
    assert manifest["counts"]["decision_reason_counts"] == {
        "nonreciprocal_intersection_partner": 2
    }
