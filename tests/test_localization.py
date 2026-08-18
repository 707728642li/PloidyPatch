from __future__ import annotations

import json
from pathlib import Path

from ploidypatch.localization import (
    score_synteny_localization,
    write_synteny_model_labels,
)
from ploidypatch.perturb import generate_missing_gene_benchmark


SOURCE_GFF = (
    "##gff-version 3\n"
    "chr1\ttest\tgene\t101\t190\t.\t+\t.\tID=gene:G1\n"
    "chr1\ttest\tmRNA\t101\t190\t.\t+\t.\tID=transcript:T1;Parent=gene:G1\n"
    "chr1\ttest\texon\t101\t190\t.\t+\t.\tParent=transcript:T1\n"
    "chr1\ttest\tCDS\t110\t180\t.\t+\t0\tParent=transcript:T1\n"
)


def test_scores_gap_and_model_localization_independently(tmp_path: Path) -> None:
    source = tmp_path / "source.gff3"
    source.write_text(SOURCE_GFF, encoding="utf-8")
    benchmark = tmp_path / "benchmark"
    generate_missing_gene_benchmark(source, benchmark, count=1, seed=1)

    gaps = tmp_path / "gaps.tsv"
    gaps.write_text(
        "candidate_id\tgap_id\tsource_label\tquery_seqid\tlocus_start\t"
        "locus_end\ttarget_gene_id\n"
        "C1\tGAP1\tbra_a\tchr1\t90\t200\tT1\n"
        "C2\tGAP2\tbra_a\tchr1\t300\t400\tT2\n",
        encoding="utf-8",
    )
    selection = tmp_path / "selection.tsv"
    selection.write_text(
        "candidate_id\tgap_id\tsource_label\tquery_seqid\tlocus_start\t"
        "locus_end\ttarget_gene_id\tmodel_id\tmodel_seqid\tmodel_start\t"
        "model_end\tmodel_strand\n"
        "C1\tGAP1\tbra_a\tchr1\t90\t200\tT1\tM1\tchr1\t110\t180\t+\n"
        "C2\tGAP2\tbra_a\tchr1\t300\t400\tT2\tM2\tchr1\t320\t380\t+\n",
        encoding="utf-8",
    )
    report = score_synteny_localization(
        gap_tsv_paths=[gaps],
        selection_tsv_path=selection,
        truth_path=benchmark / "hidden_truth.json",
        include_event_details=True,
    )

    contained = report["all_blind_gap_loci"]["full_gene_span_containment"]
    model_overlap = report["selected_model_spans"]["same_strand_overlap"]
    reciprocal = report["selected_model_spans"][
        "same_strand_reciprocal_overlap_at_least_0_5"
    ]
    assert contained["candidate_precision"] == 0.5
    assert contained["event_recall"] == 1.0
    assert model_overlap["candidate_precision"] == 0.5
    assert model_overlap["event_recall"] == 1.0
    assert reciprocal["candidate_precision"] == 0.5
    assert reciprocal["event_recall"] == 1.0
    assert report["quality_gate"]["grade"] == "pass"
    assert report["event_details"][0]["same_strand_overlap"] is True
    assert "gene:G1" in json.dumps(report["event_details"])


def test_writes_exact_cds_labels_for_selected_models(tmp_path: Path) -> None:
    source = tmp_path / "source.gff3"
    source.write_text(SOURCE_GFF, encoding="utf-8")
    benchmark = tmp_path / "benchmark"
    generate_missing_gene_benchmark(source, benchmark, count=1, seed=1)
    selection = tmp_path / "selection.tsv"
    selection.write_text(
        "candidate_id\tgap_id\tsource_label\tquery_seqid\tlocus_start\t"
        "locus_end\ttarget_gene_id\tmodel_id\tmodel_seqid\tmodel_start\t"
        "model_end\tmodel_strand\n"
        "C1\tGAP1\tbra_a\tchr1\t90\t200\tT1\tM1\tchr1\t110\t180\t+\n",
        encoding="utf-8",
    )
    decisions = tmp_path / "decisions.tsv"
    decisions.write_text(
        "model_id\tstatus\treason\texisting_cds_overlap_fraction\n"
        "M1\taccepted\taccepted\t0\n",
        encoding="utf-8",
    )
    candidate = tmp_path / "candidate.gff3"
    candidate.write_text(
        (benchmark / "perturbed.gff3").read_text(encoding="utf-8")
        + "chr1\tPloidyPatchBaseline\tgene\t110\t180\t.\t+\t.\t"
        "ID=PG1;miniprot_model=M1\n"
        + "chr1\tPloidyPatchBaseline\tmRNA\t110\t180\t.\t+\t.\t"
        "ID=PT1;Parent=PG1;miniprot_model=M1\n"
        + "chr1\tPloidyPatchBaseline\tCDS\t110\t180\t.\t+\t0\tParent=PT1\n",
        encoding="utf-8",
    )
    labels = tmp_path / "labels.tsv"
    manifest = write_synteny_model_labels(
        source_gff_path=source,
        candidate_gff_path=candidate,
        selection_tsv_path=selection,
        baseline_decisions_tsv_path=decisions,
        truth_path=benchmark / "hidden_truth.json",
        output_tsv_path=labels,
        control_candidate_gff_path=source,
    )

    header, values = labels.read_text(encoding="utf-8").splitlines()
    row = dict(zip(header.split("\t"), values.split("\t"), strict=True))
    assert row["label_same_strand_overlap"] == "1"
    assert row["label_exact_cds_chain"] == "1"
    assert row["label_paired_differential_cds_chain"] == "1"
    assert row["baseline_existing_cds_overlap_fraction"] == "0"
    assert manifest["counts"]["label_exact_cds_chain"] == 1
