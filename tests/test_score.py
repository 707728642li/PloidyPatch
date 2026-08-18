from __future__ import annotations

import json
from pathlib import Path

import pytest

from ploidypatch.perturb import (
    generate_missing_gene_benchmark,
    restore_gff_from_truth,
)
from ploidypatch.score import score_annotation_repair


SOURCE_GFF = (
    "##gff-version 3\n"
    "chr1\ttest\tgene\t1\t90\t.\t+\t.\tID=gene:G1\n"
    "chr1\ttest\tmRNA\t1\t90\t.\t+\t.\tID=transcript:T1;Parent=gene:G1\n"
    "chr1\ttest\texon\t1\t30\t.\t+\t.\tID=exon:E1;Parent=transcript:T1\n"
    "chr1\ttest\tCDS\t1\t30\t.\t+\t0\tID=CDS:C1;Parent=transcript:T1\n"
    "chr1\ttest\tgene\t101\t190\t.\t-\t.\tID=gene:G2\n"
    "chr1\ttest\tmRNA\t101\t190\t.\t-\t.\tID=transcript:T2;Parent=gene:G2\n"
    "chr1\ttest\texon\t161\t190\t.\t-\t.\tID=exon:E2;Parent=transcript:T2\n"
    "chr1\ttest\tCDS\t161\t190\t.\t-\t0\tID=CDS:C2;Parent=transcript:T2\n"
)


def write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8", newline="")
    return path


def make_benchmark(tmp_path: Path) -> tuple[Path, Path]:
    source = write(tmp_path / "source.gff3", SOURCE_GFF)
    run = tmp_path / "run"
    generate_missing_gene_benchmark(source, run, count=1, seed=5)
    return source, run


def test_perfect_repair_scores_one_without_identifier_dependence(tmp_path: Path) -> None:
    source, run = make_benchmark(tmp_path)
    restored = run / "restored.gff3"
    restore_gff_from_truth(run / "perturbed.gff3", run / "hidden_truth.json", restored)

    report = score_annotation_repair(
        source,
        run / "perturbed.gff3",
        restored,
        run / "hidden_truth.json",
        include_event_details=True,
    )

    assert report["strict_transcript_structure"]["precision"] == 1.0
    assert report["strict_transcript_structure"]["recall"] == 1.0
    assert report["strict_transcript_structure"]["f1"] == 1.0
    assert report["strict_cds_chain"]["f1"] == 1.0
    assert report["partial_structure"]["exon_segments"]["f1"] == 1.0
    assert report["partial_structure"]["cds_segments_with_phase"]["f1"] == 1.0
    assert report["partial_structure"]["cds_nucleotide_coverage"]["f1"] == 1.0
    assert report["event_recovery"]["complete_transcript_recall"] == 1.0
    assert report["event_recovery"]["exact_gene_recall"] == 1.0
    assert report["event_recovery"]["complete_cds_chain_recall"] == 1.0
    assert report["event_recovery"]["exact_cds_gene_recall"] == 1.0
    assert report["collateral_changes"][
        "baseline_transcript_structures_missing_from_candidate"
    ] == 0
    assert len(report["event_details"]) == 1


def test_no_repair_has_zero_recall_and_no_defined_precision(tmp_path: Path) -> None:
    source, run = make_benchmark(tmp_path)
    report = score_annotation_repair(
        source,
        run / "perturbed.gff3",
        run / "perturbed.gff3",
        run / "hidden_truth.json",
    )

    metrics = report["strict_transcript_structure"]
    assert metrics["true_positive"] == 0
    assert metrics["false_positive"] == 0
    assert metrics["false_negative"] == 1
    assert metrics["precision"] is None
    assert metrics["recall"] == 0.0
    assert report["strict_cds_chain"]["recall"] == 0.0
    assert report["event_recovery"]["complete_transcript_recall"] == 0.0
    assert report["partial_structure"]["exon_segments"]["recall"] == 0.0
    assert report["partial_structure"]["cds_nucleotide_coverage"]["recall"] == 0.0


def test_spurious_novel_model_reduces_precision(tmp_path: Path) -> None:
    source, run = make_benchmark(tmp_path)
    candidate = run / "candidate_with_false_positive.gff3"
    candidate.write_text(
        SOURCE_GFF
        + "chr1\ttest\tgene\t300\t330\t.\t+\t.\tID=gene:FP\n"
        + "chr1\ttest\tmRNA\t300\t330\t.\t+\t.\tID=transcript:FP;Parent=gene:FP\n"
        + "chr1\ttest\texon\t300\t330\t.\t+\t.\tParent=transcript:FP\n"
        + "chr1\ttest\tCDS\t300\t329\t.\t+\t0\tParent=transcript:FP\n",
        encoding="utf-8",
        newline="",
    )
    report = score_annotation_repair(
        source,
        run / "perturbed.gff3",
        candidate,
        run / "hidden_truth.json",
    )

    metrics = report["strict_transcript_structure"]
    assert metrics["true_positive"] == 1
    assert metrics["false_positive"] == 1
    assert metrics["precision"] == 0.5
    assert metrics["recall"] == 1.0


def test_paired_control_subtracts_reproducible_background_models(
    tmp_path: Path,
) -> None:
    source, run = make_benchmark(tmp_path)
    restored = run / "restored.gff3"
    restore_gff_from_truth(run / "perturbed.gff3", run / "hidden_truth.json", restored)
    false_model = (
        "chr1\ttest\tgene\t300\t330\t.\t+\t.\tID=gene:FP\n"
        "chr1\ttest\tmRNA\t300\t330\t.\t+\t.\tID=transcript:FP;Parent=gene:FP\n"
        "chr1\ttest\texon\t300\t330\t.\t+\t.\tParent=transcript:FP\n"
        "chr1\ttest\tCDS\t300\t329\t.\t+\t0\tParent=transcript:FP\n"
    )
    candidate = run / "paired_candidate.gff3"
    candidate.write_text(
        restored.read_text(encoding="utf-8") + false_model,
        encoding="utf-8",
        newline="",
    )
    control = run / "complete_control_candidate.gff3"
    control.write_text(
        source.read_text(encoding="utf-8") + false_model,
        encoding="utf-8",
        newline="",
    )

    report = score_annotation_repair(
        source,
        run / "perturbed.gff3",
        candidate,
        run / "hidden_truth.json",
        control_candidate_gff_path=control,
    )

    assert report["evaluation_mode"] == "paired_complete_annotation_difference"
    assert report["strict_transcript_structure"]["precision"] == 1.0
    assert report["strict_cds_chain"]["precision"] == 1.0
    assert report["background_subtraction"][
        "candidate_novel_transcript_structures_before_subtraction"
    ] == 2
    assert report["background_subtraction"][
        "differential_candidate_transcript_structures"
    ] == 1
    assert report["quality_gate"]["control_retains_all_source_structures"] is True


def test_partial_structure_scores_inexact_model_components(tmp_path: Path) -> None:
    source, run = make_benchmark(tmp_path)
    truth_path = run / "hidden_truth.json"
    truth = json.loads(truth_path.read_text(encoding="utf-8"))
    removed = [
        record["raw_line"].rstrip("\n").split("\t")
        for record in truth["events"][0]["removed_records"]
    ]
    target = truth["events"][0]["target"]
    exon = next(fields for fields in removed if fields[2] == "exon")
    cds = next(fields for fields in removed if fields[2] == "CDS")
    changed_exon_start = int(exon[3]) + 1
    candidate = run / "partial_candidate.gff3"
    with (run / "perturbed.gff3").open("r", encoding="utf-8") as handle:
        perturbed_text = handle.read()
    candidate.write_text(
        perturbed_text
        + f"{target['seqid']}\ttest\tgene\t{target['start']}\t{target['end']}\t.\t"
        f"{target['strand']}\t.\tID=gene:partial\n"
        + f"{target['seqid']}\ttest\tmRNA\t{target['start']}\t{target['end']}\t.\t"
        f"{target['strand']}\t.\tID=transcript:partial;Parent=gene:partial\n"
        + f"{exon[0]}\ttest\texon\t{changed_exon_start}\t{exon[4]}\t.\t"
        f"{exon[6]}\t.\tParent=transcript:partial\n"
        + f"{cds[0]}\ttest\tCDS\t{cds[3]}\t{cds[4]}\t.\t{cds[6]}\t"
        f"{cds[7]}\tParent=transcript:partial\n",
        encoding="utf-8",
        newline="",
    )

    report = score_annotation_repair(
        source,
        run / "perturbed.gff3",
        candidate,
        truth_path,
    )

    assert report["strict_transcript_structure"]["true_positive"] == 0
    assert report["strict_transcript_structure"]["false_positive"] == 1
    assert report["partial_structure"]["exon_segments"]["true_positive"] == 0
    assert report["partial_structure"]["cds_segments_with_phase"]["recall"] == 1.0
    assert report["partial_structure"]["cds_nucleotide_coverage"]["recall"] == 1.0
    assert report["strict_cds_chain"]["recall"] == 1.0


def test_score_reports_anonymous_stratified_recall(tmp_path: Path) -> None:
    source, run = make_benchmark(tmp_path)
    truth = json.loads((run / "hidden_truth.json").read_text(encoding="utf-8"))
    target_gene = truth["events"][0]["target"]["gene_id"]
    strata = write(
        tmp_path / "strata.tsv",
        "gene_id\tsubgenome\tsynteny\n"
        f"{target_gene}\tA\texpected_only\n",
    )

    report = score_annotation_repair(
        source,
        run / "perturbed.gff3",
        source,
        run / "hidden_truth.json",
        event_strata_path=strata,
        stratum_columns=("subgenome", "synteny"),
    )

    subgenome = report["stratified_recall"]["marginal"]["subgenome"][0]
    assert subgenome["subgenome"] == "A"
    assert subgenome["events"] == 1
    assert subgenome["strict_transcript_structures"]["recall"] == 1.0
    assert subgenome["strict_cds_chains"]["recall"] == 1.0
    assert report["stratified_recall"]["joint"][0]["synteny"] == "expected_only"
    assert "gene_id" not in json.dumps(report["stratified_recall"])


def test_score_fails_when_event_strata_are_incomplete(tmp_path: Path) -> None:
    source, run = make_benchmark(tmp_path)
    strata = write(tmp_path / "strata.tsv", "gene_id\tsubgenome\nother\tA\n")
    with pytest.raises(ValueError, match="absent from event strata"):
        score_annotation_repair(
            source,
            run / "perturbed.gff3",
            source,
            run / "hidden_truth.json",
            event_strata_path=strata,
            stratum_columns=("subgenome",),
        )
