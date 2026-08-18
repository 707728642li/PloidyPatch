from __future__ import annotations

import csv
import json
from pathlib import Path

from ploidypatch.isoseq_validation import (
    bootstrap_isoseq_review_yield,
    filter_candidate_query_paf,
    join_isoseq_review_rankings,
    prepare_b73_isoseq_transcripts,
    validate_isoseq_candidate_chains,
)


def write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8", newline="")
    return path


def test_prepare_b73_isoseq_uses_only_declared_pure_b73_columns(tmp_path: Path) -> None:
    counts = write(
        tmp_path / "counts.csv",
        "id,EM1,R2,END1,R1\n"
        "PB.1.1,1,100,0,1\n"
        "PB.1.2,0,100,0,1\n",
    )
    fasta = write(tmp_path / "all.fa", ">PB.1.1\nACGT\n>PB.1.2\nTTTT\n")
    output_fasta = tmp_path / "b73.fa"
    output_counts = tmp_path / "b73.tsv"

    manifest = prepare_b73_isoseq_transcripts(
        fasta_path=fasta,
        count_csv_path=counts,
        output_fasta_path=output_fasta,
        output_count_tsv_path=output_counts,
        minimum_b73_full_length_reads=2,
    )

    assert output_fasta.read_text(encoding="ascii") == ">PB.1.1\nACGT\n"
    assert manifest["counts"]["selected_transcripts"] == 1
    assert manifest["parameters"]["other_genotypes_excluded"] is True
    with output_counts.open("r", encoding="utf-8", newline="") as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    assert row["b73_full_length_reads"] == "2"


def test_validate_isoseq_requires_full_cds_coverage_and_exact_junctions(
    tmp_path: Path,
) -> None:
    sequence = list("A" * 200)
    sequence[40:42] = "GT"
    sequence[57:59] = "AG"
    genome = write(tmp_path / "genome.fa", ">chr1\n" + "".join(sequence) + "\n")
    write(tmp_path / "genome.fa.fai", "chr1\t200\t6\t200\t201\n")
    digest = "a" * 64
    candidate = write(
        tmp_path / "candidate.gff3",
        "##gff-version 3\n"
        f"chr1\tPloidyPatchConsensus\tgene\t21\t80\t.\t+\t.\tID=g;consensus_digest={digest}\n"
        f"chr1\tPloidyPatchConsensus\tmRNA\t21\t80\t.\t+\t.\tID=t;Parent=g;consensus_digest={digest}\n"
        "chr1\tPloidyPatchConsensus\tCDS\t21\t40\t.\t+\t0\tParent=t\n"
        "chr1\tPloidyPatchConsensus\tCDS\t60\t80\t.\t+\t0\tParent=t\n",
    )
    selected_counts = write(
        tmp_path / "selected.tsv",
        "transcript_id\tEM1\tR1\tEND1\tb73_full_length_reads\n"
        "PB.1.1\t2\t0\t0\t2\n",
    )
    paf = write(
        tmp_path / "alignments.paf",
        "PB.1.1\t60\t0\t60\t+\tchr1\t200\t10\t89\t60\t60\t60"
        "\ttp:A:P\tcg:Z:30M19N30M\tAS:i:120\n",
    )
    output = tmp_path / "evidence.tsv"

    manifest = validate_isoseq_candidate_chains(
        candidate_gff_path=candidate,
        paf_path=paf,
        selected_count_tsv_path=selected_counts,
        genome_fasta_path=genome,
        output_evidence_tsv_path=output,
        flank_bp=5,
    )

    with output.open("r", encoding="utf-8", newline="") as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    assert row["evidence_state"] == "full_chain_supported"
    assert row["matched_candidate_junctions"] == "1"
    assert float(row["best_candidate_cds_coverage"]) == 1.0
    assert row["supporting_full_length_reads"] == "2"
    assert row["supporting_query_groups"] == ""
    assert manifest["counts"]["alignment_decisions"]["usable"] == 1
    assert manifest["counts"]["evidence_states"] == {"full_chain_supported": 1}
    assert manifest["parameters"]["selected_count_field"] == "b73_full_length_reads"

    generic_counts = write(
        tmp_path / "raw_read_counts.tsv",
        "transcript_id\tfull_length_reads\nCRR1429911|read1\t1\n",
    )
    generic_paf = write(
        tmp_path / "raw_read_alignments.paf",
        "CRR1429911|read1\t60\t0\t60\t+\tchr1\t200\t10\t89\t60\t60\t60"
        "\ttp:A:P\tcg:Z:30M19N30M\tAS:i:120\n",
    )
    generic_output = tmp_path / "raw_read_evidence.tsv"
    generic_manifest = validate_isoseq_candidate_chains(
        candidate_gff_path=candidate,
        paf_path=generic_paf,
        selected_count_tsv_path=generic_counts,
        genome_fasta_path=genome,
        output_evidence_tsv_path=generic_output,
        flank_bp=5,
    )
    with generic_output.open("r", encoding="utf-8", newline="") as handle:
        generic_row = next(csv.DictReader(handle, delimiter="\t"))
    assert generic_row["supporting_full_length_reads"] == "1"
    assert generic_row["supporting_query_groups"] == "CRR1429911"
    assert generic_manifest["parameters"]["selected_count_field"] == "full_length_reads"


def test_minimap2_ts_mode_combines_paf_orientation_and_ts_relation(
    tmp_path: Path,
) -> None:
    sequence = list("A" * 200)
    sequence[40:42] = "GT"
    sequence[57:59] = "AG"
    genome = write(tmp_path / "genome.fa", ">chr1\n" + "".join(sequence) + "\n")
    write(tmp_path / "genome.fa.fai", "chr1\t200\t6\t200\t201\n")
    digest = "c" * 64
    candidate = write(
        tmp_path / "candidate.gff3",
        "##gff-version 3\n"
        f"chr1\tPloidyPatchConsensus\tmRNA\t21\t80\t.\t+\t.\t"
        f"ID=t;consensus_digest={digest}\n"
        "chr1\tPloidyPatchConsensus\tCDS\t21\t40\t.\t+\t0\tParent=t\n"
        "chr1\tPloidyPatchConsensus\tCDS\t60\t80\t.\t+\t0\tParent=t\n",
    )
    counts = write(
        tmp_path / "counts.tsv",
        "transcript_id\tfull_length_reads\nread_with_ts\t1\nread_without_ts\t1\n",
    )
    paf = write(
        tmp_path / "alignments.paf",
        # PAF '-' times ts '-' gives a '+' transcript relative to reference.
        "read_with_ts\t60\t0\t60\t-\tchr1\t200\t10\t89\t60\t60\t60"
        "\ttp:A:P\tts:A:-\tcg:Z:30M19N30M\tAS:i:120\n"
        "read_without_ts\t60\t0\t60\t+\tchr1\t200\t10\t89\t60\t60\t60"
        "\ttp:A:P\tcg:Z:30M19N30M\tAS:i:120\n",
    )
    output = tmp_path / "evidence.tsv"

    manifest = validate_isoseq_candidate_chains(
        candidate_gff_path=candidate,
        paf_path=paf,
        selected_count_tsv_path=counts,
        genome_fasta_path=genome,
        output_evidence_tsv_path=output,
        alignment_strand_source="minimap2_ts",
        flank_bp=5,
    )

    with output.open("r", encoding="utf-8", newline="") as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    assert row["evidence_state"] == "full_chain_supported"
    assert row["supporting_transcripts"] == "read_with_ts"
    assert manifest["parameters"]["alignment_strand_source"] == "minimap2_ts"
    assert manifest["counts"]["paf_strand_availability"] == {
        "available": 1,
        "unavailable": 1,
    }
    assert manifest["counts"]["alignment_decisions"] == {
        "transcript_strand_unavailable": 1,
        "usable": 1,
    }


def test_candidate_query_filter_is_ts_aware_and_retains_all_query_alignments(
    tmp_path: Path,
) -> None:
    digest = "d" * 64
    candidate = write(
        tmp_path / "candidate.gff3",
        "##gff-version 3\n"
        f"chr1\tPloidyPatchConsensus\tmRNA\t21\t80\t.\t+\t.\t"
        f"ID=t;consensus_digest={digest}\n"
        "chr1\tPloidyPatchConsensus\tCDS\t21\t40\t.\t+\t0\tParent=t\n"
        "chr1\tPloidyPatchConsensus\tCDS\t60\t80\t.\t+\t0\tParent=t\n",
    )
    paf = write(
        tmp_path / "run.paf",
        "q1\t60\t0\t60\t-\tchr1\t200\t10\t89\t60\t60\t60"
        "\ttp:A:P\tts:A:-\tcg:Z:30M19N30M\tAS:i:120\n"
        # q1 is retained losslessly even though this secondary is opposite.
        "q1\t60\t0\t60\t+\tchr1\t200\t10\t89\t59\t60\t10"
        "\ttp:A:S\tts:A:-\tcg:Z:30M19N30M\tAS:i:80\n"
        "q2\t60\t0\t60\t+\tchr1\t200\t10\t89\t60\t60\t60"
        "\ttp:A:P\tts:A:-\tcg:Z:30M19N30M\tAS:i:120\n"
        "q3\t60\t0\t60\t+\tchr1\t200\t10\t89\t60\t60\t60"
        "\ttp:A:P\tcg:Z:30M19N30M\tAS:i:120\n",
    )

    manifest = filter_candidate_query_paf(
        candidate_gff_path=candidate,
        paf_inputs=(("RUN1", paf),),
        output_paf_path=tmp_path / "filtered.paf",
        output_count_tsv_path=tmp_path / "counts.tsv",
        output_summary_tsv_path=tmp_path / "summary.tsv",
        output_manifest_json_path=tmp_path / "manifest.json",
        alignment_strand_source="minimap2_ts",
    )

    filtered = (tmp_path / "filtered.paf").read_text(encoding="utf-8")
    assert filtered.count("\n") == 2
    assert all(line.startswith("RUN1|q1\t") for line in filtered.splitlines())
    assert (tmp_path / "counts.tsv").read_text(encoding="utf-8") == (
        "transcript_id\tfull_length_reads\nRUN1|q1\t1\n"
    )
    assert manifest["counts"]["retained_queries"] == 1
    assert manifest["counts"]["retained_alignments"] == 2
    assert manifest["counts"]["strand_unavailable_alignments"] == 1


def test_join_isoseq_review_rankings_reports_predeclared_topology_delta(
    tmp_path: Path,
) -> None:
    digest_a = "a" * 64
    digest_b = "b" * 64
    oversized_support_field = "read," * 30_000
    evidence = write(
        tmp_path / "evidence.tsv",
        "candidate_digest\tevidence_state\tsupporting_transcripts\n"
        f"{digest_a}\tfull_chain_supported\t{oversized_support_field}\n"
        f"{digest_b}\tno_qualifying_observation\t\n",
    )
    rankings = write(
        tmp_path / "rankings.tsv",
        "estimator\treview_rank\tcandidate_digest\n"
        f"baseline\t1\t{digest_b}\n"
        f"baseline\t2\t{digest_a}\n"
        f"topology\t1\t{digest_a}\n"
        f"topology\t2\t{digest_b}\n",
    )

    summary = join_isoseq_review_rankings(
        evidence_tsv_path=evidence,
        review_rankings_tsv_path=rankings,
        output_tsv_path=tmp_path / "joined.tsv",
        output_summary_json_path=tmp_path / "summary.json",
        review_budgets=(1, 2, 100),
    )

    assert summary["primary"]["full_chain_supported"] == {
        "baseline": 1,
        "topology": 1,
    }
    assert summary["primary"]["topology_minus_baseline"] == 0
    disk_summary = json.loads((tmp_path / "summary.json").read_text(encoding="utf-8"))
    topology_at_one = next(
        row
        for row in disk_summary["yield"]
        if row["estimator"] == "topology" and row["review_budget"] == 1
    )
    assert topology_at_one["full_chain_supported"] == 1

    generic_rankings = write(
        tmp_path / "generic_rankings.tsv",
        rankings.read_text(encoding="utf-8").replace("topology", "v04_guard"),
    )
    generic_summary = join_isoseq_review_rankings(
        evidence_tsv_path=evidence,
        review_rankings_tsv_path=generic_rankings,
        output_tsv_path=tmp_path / "generic_joined.tsv",
        output_summary_json_path=tmp_path / "generic_summary.json",
        review_budgets=(1, 2, 100),
        comparator_estimator="baseline",
        primary_estimator="v04_guard",
    )
    assert generic_summary["primary"]["estimator_delta"] == {
        "primary_estimator": "v04_guard",
        "comparator_estimator": "baseline",
        "delta_full_chain_supported": 0,
    }
    assert "topology_minus_baseline" not in generic_summary["primary"]


def test_isoseq_review_bootstrap_preserves_top_k_chromosome_allocation(
    tmp_path: Path,
) -> None:
    evidence = tmp_path / "evidence.tsv"
    rankings = tmp_path / "rankings.tsv"
    # The production function requires the preregistered K=100, so the toy
    # universe is large enough to exercise that contract without duplicates.
    evidence_rows = ["candidate_digest\tseqid\tevidence_state\n"]
    ranking_rows = ["estimator\treview_rank\tcandidate_digest\n"]
    large_digests = [f"{index:064x}" for index in range(120)]
    for index, digest in enumerate(large_digests):
        evidence_rows.append(
            f"{digest}\tchr{index % 2 + 1}\t"
            + ("full_chain_supported\n" if index < 10 else "no_qualifying_observation\n")
        )
    for estimator, order in (
        ("baseline", list(range(120))),
        ("topology", list(range(9, -1, -1)) + list(range(10, 120))),
    ):
        for rank, index in enumerate(order, start=1):
            ranking_rows.append(f"{estimator}\t{rank}\t{large_digests[index]}\n")
    evidence.write_text("".join(evidence_rows), encoding="utf-8", newline="")
    rankings.write_text("".join(ranking_rows), encoding="utf-8", newline="")

    report = bootstrap_isoseq_review_yield(
        evidence_tsv_path=evidence,
        review_rankings_tsv_path=rankings,
        output_json_path=tmp_path / "bootstrap.json",
        review_budgets=(25, 50, 100),
        replicates=200,
        seed=7,
    )

    assert report["counts"]["target_chromosomes"] == 2
    baseline_25 = next(
        row
        for row in report["random_ranking_enrichment"]
        if row["estimator"] == "baseline" and row["review_budget"] == 25
    )
    assert sum(baseline_25["chromosome_allocation"].values()) == 25

    generic_rankings = tmp_path / "generic_bootstrap_rankings.tsv"
    generic_rankings.write_text(
        rankings.read_text(encoding="utf-8").replace("topology", "v04_guard"),
        encoding="utf-8",
        newline="",
    )
    generic_report = bootstrap_isoseq_review_yield(
        evidence_tsv_path=evidence,
        review_rankings_tsv_path=generic_rankings,
        output_json_path=tmp_path / "generic_bootstrap.json",
        review_budgets=(25, 50, 100),
        replicates=200,
        seed=7,
        comparator_estimator="baseline",
        primary_estimator="v04_guard",
    )
    assert generic_report["parameters"]["primary_estimator"] == "v04_guard"
    assert generic_report["primary_estimator_delta"]["primary_estimator"] == "v04_guard"
    assert "primary_topology_minus_baseline" not in generic_report
