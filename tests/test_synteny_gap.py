from __future__ import annotations

import csv
from pathlib import Path

from ploidypatch.synteny_gap import (
    infer_wgdi_synteny_gaps,
    select_synteny_gap_models,
)


def test_infers_and_deduplicates_strict_adjacent_anchor_gaps(tmp_path: Path) -> None:
    query_gff = tmp_path / "query.gff"
    query_gff.write_text(
        "qchr\tQ1\t100\t200\t+\t1\n"
        "qchr\tQ2\t500\t600\t+\t2\n"
        "qchr\tQ3\t1000\t1100\t-\t3\n",
        encoding="utf-8",
    )
    target_gff = tmp_path / "target.gff"
    target_gff.write_text(
        "tchr\tT1\t100\t200\t+\t1\n"
        "tchr\tTM\t300\t400\t-\t2\n"
        "tchr\tT2\t500\t600\t+\t3\n"
        "tchr\tT3\t700\t800\t+\t4\n"
        "t2\tX1\t100\t200\t+\t1\n"
        "t2\tXM\t300\t400\t+\t2\n"
        "t2\tX2\t500\t600\t+\t3\n",
        encoding="utf-8",
    )
    collinearity = tmp_path / "blocks.tsv"
    collinearity.write_text(
        "# Alignment 1: score=100 pvalue=0.01 N=3 qchr&tchr plus\n"
        "Q1 1 T1 1 1\n"
        "Q2 2 T2 3 1\n"
        "Q3 3 T3 4 1\n"
        "# Alignment 2: score=90 pvalue=0.02 N=2 qchr&tchr plus\n"
        "Q1 1 T1 1 1\n"
        "Q2 2 T2 3 1\n"
        "# Alignment 3: score=80 pvalue=0.03 N=2 qchr&tchr plus\n"
        "Q1 1 T1 1 1\n"
        "Q3 3 T3 4 1\n"
        "# Alignment 4: score=70 pvalue=0.04 N=2 qchr&t2 plus\n"
        "Q1 1 X1 1 1\n"
        "Q2 2 X2 3 1\n",
        encoding="utf-8",
    )
    expected_pairs = tmp_path / "expected_pairs.tsv"
    expected_pairs.write_text(
        "query_seqid\ttarget_seqid\tsource_label\nqchr\ttchr\tbra_a\n",
        encoding="utf-8",
    )
    output = tmp_path / "gaps.tsv"
    manifest = infer_wgdi_synteny_gaps(
        query_wgdi_gff_path=query_gff,
        target_wgdi_gff_path=target_gff,
        collinearity_path=collinearity,
        source_label="bra_a",
        output_tsv_path=output,
        expected_chromosome_pair_tsv_path=expected_pairs,
        max_query_intervening_genes=0,
        min_target_excess_genes=1,
        max_target_gap_genes=3,
        max_query_locus_bp=1000,
    )
    with output.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 1
    row = rows[0]
    assert row["target_gene_id"] == "TM"
    assert row["query_left_anchor_gene"] == "Q1"
    assert row["query_right_anchor_gene"] == "Q2"
    assert row["locus_start"] == "201"
    assert row["locus_end"] == "499"
    assert row["target_gap_genes"] == "1"
    assert row["target_excess_genes"] == "1"
    assert row["supporting_block_count"] == "2"
    assert row["best_block_id"] == "bra_a:1"
    assert manifest["counts"]["raw_gene_hypotheses"] == 2
    assert manifest["counts"]["unique_gene_hypotheses"] == 1
    assert manifest["counts"]["anchor_pairs_query_gap_too_large"] == 1
    assert manifest["counts"]["blocks_expected_chromosome_pair"] == 3
    assert manifest["counts"]["blocks_outside_expected_chromosome_pair"] == 1
    assert manifest["inputs"]["expected_chromosome_pairs"]["pairs"] == 1


def test_selects_only_local_accepted_models_and_retains_blind_gff(
    tmp_path: Path,
) -> None:
    gap = tmp_path / "gaps.tsv"
    gap.write_text(
        "candidate_id\tgap_id\tsource_label\tquery_seqid\tlocus_start\t"
        "locus_end\tlocus_span_bp\tquery_left_anchor_gene\t"
        "query_right_anchor_gene\ttarget_seqid\ttarget_gene_id\t"
        "target_gene_order\ttarget_gap_genes\ttarget_excess_genes\t"
        "best_block_id\tbest_block_score\tbest_block_pvalue\t"
        "best_block_pairs\tsupporting_block_count\tsupporting_block_ids\n"
        "C1\tG1\tbra_a\tchrA01\t201\t499\t299\tQ1\tQ2\tA01\tT1\t2\t"
        "1\t1\tbra_a:1\t100\t0.01\t3\t2\tbra_a:1,bra_a:2\n",
        encoding="utf-8",
    )
    decisions = tmp_path / "decisions.tsv"
    decisions.write_text(
        "model_id\tquery_id\tsource\tseqid\tstart\tend\tstrand\tscore\trank\t"
        "identity\tquery_coverage\tframeshifts\tstop_codons\tstatus\treason\n"
        "M1\tbra_a__T1\tbra_a\tchrA01\t250\t450\t+\t90\t1\t0.90\t0.95\t"
        "0\t0\taccepted\taccepted\n"
        "M2\tbra_a__T1\tbra_a\tchrA01\t150\t450\t+\t99\t1\t0.99\t0.99\t"
        "0\t0\taccepted\taccepted\n"
        "M3\tbra_a__T1\tbra_a\tchrA01\t260\t440\t+\t100\t1\t0.99\t0.99\t"
        "0\t0\trejected\toverlaps_blind_cds\n",
        encoding="utf-8",
    )
    candidate = tmp_path / "adapted.gff3"
    candidate.write_text(
        "##gff-version 3\n"
        "chrA01\tblind\tgene\t1\t100\t.\t+\t.\tID=blind_gene\n"
        "###\n"
        "chrA01\tPloidyPatchBaseline\tgene\t250\t450\t90\t+\t.\t"
        "ID=PPBASE_gene_000001;miniprot_model=M1\n"
        "chrA01\tPloidyPatchBaseline\tmRNA\t250\t450\t90\t+\t.\t"
        "ID=PPBASE_tx_000001;Parent=PPBASE_gene_000001;miniprot_model=M1\n"
        "chrA01\tPloidyPatchBaseline\tCDS\t250\t450\t.\t+\t0\t"
        "Parent=PPBASE_tx_000001\n"
        "chrA01\tPloidyPatchBaseline\tgene\t150\t450\t99\t+\t.\t"
        "ID=PPBASE_gene_000002;miniprot_model=M2\n"
        "chrA01\tPloidyPatchBaseline\tmRNA\t150\t450\t99\t+\t.\t"
        "ID=PPBASE_tx_000002;Parent=PPBASE_gene_000002;miniprot_model=M2\n",
        encoding="utf-8",
    )
    selection = tmp_path / "selection.tsv"
    output_gff = tmp_path / "selected.gff3"
    manifest = select_synteny_gap_models(
        gap_tsv_paths=[gap],
        baseline_decisions_tsv_path=decisions,
        adapted_candidate_gff_path=candidate,
        output_selection_tsv_path=selection,
        output_candidate_gff_path=output_gff,
    )

    with selection.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    output_text = output_gff.read_text(encoding="utf-8")
    assert len(rows) == 1
    assert rows[0]["model_id"] == "M1"
    assert "ID=blind_gene" in output_text
    assert "miniprot_model=M1" in output_text
    assert "miniprot_model=M2" not in output_text
    assert manifest["counts"]["selected_models"] == 1
    assert manifest["counts"]["adapted_baseline_features_selected"] == 3
