from __future__ import annotations

import csv
from pathlib import Path

import pytest

from ploidypatch.consensus import select_method_consensus
from ploidypatch.perturb import read_gff_document
from ploidypatch.score import build_annotation_index


def write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8", newline="")
    return path


def candidate(base_text: str, models: list[tuple[str, int, int, str]]) -> str:
    rows = [base_text, "###\n"]
    for number, (model, start, end, phase) in enumerate(models, start=1):
        gene = f"PPGFF_gene_{number:06d}"
        transcript = f"PPGFF_tx_{number:06d}"
        attrs = f"baseline_source=tool;upstream_model={model}"
        rows.extend(
            [
                f"chr1\tPloidyPatchBaseline\tgene\t{start}\t{end}\t.\t+\t.\tID={gene};{attrs}\n",
                f"chr1\tPloidyPatchBaseline\tmRNA\t{start}\t{end}\t.\t+\t.\tID={transcript};Parent={gene};{attrs}\n",
                f"chr1\tPloidyPatchBaseline\texon\t{start}\t{end}\t.\t+\t.\tParent={transcript}\n",
                f"chr1\tPloidyPatchBaseline\tCDS\t{start}\t{end}\t.\t+\t{phase}\tParent={transcript}\n",
            ]
        )
    return "".join(rows)


def test_method_consensus_selects_exact_chain_support(tmp_path: Path) -> None:
    base_text = (
        "##gff-version 3\n"
        "chr1\ttest\tgene\t1\t30\t.\t+\t.\tID=G1\n"
        "chr1\ttest\tmRNA\t1\t30\t.\t+\t.\tID=T1;Parent=G1\n"
        "chr1\ttest\tCDS\t1\t30\t.\t+\t0\tParent=T1\n"
    )
    base = write(tmp_path / "base.gff3", base_text)
    left = write(
        tmp_path / "left.gff3",
        candidate(base_text, [("L_SHARED", 101, 130, "0"), ("L_ONLY", 201, 230, "0")]),
    )
    right = write(
        tmp_path / "right.gff3",
        candidate(base_text, [("R_SHARED", 101, 130, "0"), ("R_ONLY", 301, 330, "0")]),
    )
    output = tmp_path / "consensus.gff3"
    decisions = tmp_path / "decisions.tsv"

    manifest = select_method_consensus(
        base_gff_path=base,
        candidate_inputs=(("left", left), ("right", right)),
        output_gff_path=output,
        decisions_tsv_path=decisions,
        min_method_support=2,
    )

    assert output.read_bytes().startswith(base.read_bytes() + b"###\n")
    index = build_annotation_index(read_gff_document(output))
    added = [
        tx.signature
        for tx_id, tx in index.transcripts.items()
        if tx_id.startswith("PPCONS_tx_")
    ]
    assert len(added) == 1
    assert added[0].cds == ((101, 130, "0"),)
    assert manifest["counts"]["accepted_models"] == 1
    assert manifest["counts"]["decision_counts"] == {
        "method_support_below_threshold": 2,
        "method_support_pass": 1,
    }
    with decisions.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    accepted = [row for row in rows if row["status"] == "accepted"]
    assert accepted[0]["support_methods"] == "left,right"
    assert accepted[0]["upstream_models"] == "left:L_SHARED,right:R_SHARED"


def test_method_consensus_rejects_wrong_base_prefix(tmp_path: Path) -> None:
    base = write(tmp_path / "base.gff3", "##gff-version 3\n")
    candidate_path = write(
        tmp_path / "candidate.gff3", "##gff-version 2\n###\n"
    )
    with pytest.raises(ValueError, match="exact base GFF prefix"):
        select_method_consensus(
            base_gff_path=base,
            candidate_inputs=(("a", candidate_path), ("b", candidate_path)),
            output_gff_path=tmp_path / "out.gff3",
            decisions_tsv_path=tmp_path / "decisions.tsv",
        )


def test_method_consensus_requires_unique_method_labels(tmp_path: Path) -> None:
    base_text = "##gff-version 3\n"
    base = write(tmp_path / "base.gff3", base_text)
    candidate_path = write(tmp_path / "candidate.gff3", base_text + "###\n")
    with pytest.raises(ValueError, match="unique"):
        select_method_consensus(
            base_gff_path=base,
            candidate_inputs=(("same", candidate_path), ("same", candidate_path)),
            output_gff_path=tmp_path / "out.gff3",
            decisions_tsv_path=tmp_path / "decisions.tsv",
        )


def test_method_consensus_interval_index_preserves_cross_bin_redundancy(
    tmp_path: Path,
) -> None:
    base_text = "##gff-version 3\n"
    base = write(tmp_path / "base.gff3", base_text)
    models = [
        ("LONG", 999_900, 1_000_100, "0"),
        ("OVERLAP", 1_000_050, 1_000_150, "0"),
        ("DISTANT", 2_000_100, 2_000_200, "0"),
    ]
    left = write(tmp_path / "left.gff3", candidate(base_text, models))
    right = write(tmp_path / "right.gff3", candidate(base_text, models))
    decisions = tmp_path / "decisions.tsv"
    manifest = select_method_consensus(
        base_gff_path=base,
        candidate_inputs=(("left", left), ("right", right)),
        output_gff_path=tmp_path / "consensus.gff3",
        decisions_tsv_path=decisions,
        min_method_support=2,
        max_redundancy_overlap=0.5,
    )
    assert manifest["counts"]["accepted_models"] == 2
    assert manifest["parameters"]["redundancy_interval_bin_bp"] == 1_000_000
    with decisions.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    reason_by_start = {int(row["start"]): row["reason"] for row in rows}
    assert reason_by_start == {
        999_900: "method_support_pass",
        1_000_050: "redundant_consensus_candidate",
        2_000_100: "method_support_pass",
    }


def test_chain_preserving_pool_retains_and_labels_alternative_chains(
    tmp_path: Path,
) -> None:
    base_text = "##gff-version 3\n"
    base = write(tmp_path / "base.gff3", base_text)
    models = [
        ("LONG", 999_900, 1_000_100, "0"),
        ("OVERLAP", 1_000_050, 1_000_150, "0"),
        ("DISTANT", 2_000_100, 2_000_200, "0"),
    ]
    left = write(tmp_path / "left.gff3", candidate(base_text, models))
    right = write(tmp_path / "right.gff3", candidate(base_text, models))
    decisions = tmp_path / "decisions.tsv"

    manifest = select_method_consensus(
        base_gff_path=base,
        candidate_inputs=(("left", left), ("right", right)),
        output_gff_path=tmp_path / "pool.gff3",
        decisions_tsv_path=decisions,
        min_method_support=2,
        max_redundancy_overlap=0.5,
        redundancy_policy="retain_distinct_chains",
    )

    assert manifest["schema_version"] == "ploidypatch.method_candidate_pool.v2"
    assert manifest["counts"]["accepted_models"] == 3
    assert manifest["counts"]["conflict_sets"] == 1
    assert manifest["counts"]["conflicted_chains"] == 2
    with decisions.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    by_start = {int(row["start"]): row for row in rows}
    assert {row["status"] for row in rows} == {"accepted"}
    assert by_start[999_900]["conflict_set_digest"]
    assert (
        by_start[999_900]["conflict_set_digest"]
        == by_start[1_000_050]["conflict_set_digest"]
    )
    assert by_start[999_900]["conflict_member_count"] == "2"
    assert by_start[2_000_100]["conflict_set_digest"] == ""
    assert by_start[2_000_100]["conflict_member_count"] == "1"


def test_chain_preserving_pool_still_deduplicates_exact_chains(
    tmp_path: Path,
) -> None:
    base_text = "##gff-version 3\n"
    base = write(tmp_path / "base.gff3", base_text)
    left = write(
        tmp_path / "left.gff3",
        candidate(base_text, [("LEFT", 101, 130, "0")]),
    )
    right = write(
        tmp_path / "right.gff3",
        candidate(base_text, [("RIGHT", 101, 130, "0")]),
    )
    decisions = tmp_path / "decisions.tsv"

    manifest = select_method_consensus(
        base_gff_path=base,
        candidate_inputs=(("left", left), ("right", right)),
        output_gff_path=tmp_path / "pool.gff3",
        decisions_tsv_path=decisions,
        min_method_support=1,
        redundancy_policy="retain_distinct_chains",
    )

    assert manifest["counts"]["input_models"] == 2
    assert manifest["counts"]["unique_cds_chains"] == 1
    assert manifest["counts"]["accepted_models"] == 1
    with decisions.open("r", encoding="utf-8", newline="") as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    assert row["support_method_count"] == "2"
    assert row["support_methods"] == "left,right"
    assert row["conflict_member_count"] == "1"
