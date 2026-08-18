from __future__ import annotations

import csv
import json
from pathlib import Path

from ploidypatch.natural_audit import (
    audit_natural_candidates,
    export_natural_candidate_cds,
)


def write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8", newline="")
    return path


def _toy_inputs(tmp_path: Path) -> tuple[Path, Path, Path, Path, Path, str]:
    digest = "a" * 64
    cds = "ATGGCCGTTACGCGACCTGACGTACGCTAA"
    assert len(cds) % 3 == 0
    sequence = (
        "A" * 10
        + cds[:15]
        + "CCCCC"
        + cds[15:]
        + "C" * (300 - 10 - len(cds) - 5)
    )
    genome = write(tmp_path / "genome.fa", f">chr1\n{sequence}\n")
    write(tmp_path / "genome.fa.fai", "chr1\t300\t6\t300\t301\n")
    candidate = write(
        tmp_path / "candidate.gff3",
        "##gff-version 3\n"
        "chr1\tPloidyPatchConsensus\tgene\t11\t45\t.\t+\t."
        f"\tID=g;consensus_digest={digest}\n"
        "chr1\tPloidyPatchConsensus\tmRNA\t11\t45\t.\t+\t."
        f"\tID=t;Parent=g;consensus_digest={digest}\n"
        "chr1\tPloidyPatchConsensus\tCDS\t11\t25\t.\t+\t0\tParent=t\n"
        "chr1\tPloidyPatchConsensus\tCDS\t31\t45\t.\t+\t0\tParent=t\n",
    )
    base = write(
        tmp_path / "base.gff3",
        "##gff-version 3\n"
        "chr1\tbase\tgene\t100\t120\t.\t+\t.\tID=existing\n"
        "chr1\tbase\tmRNA\t100\t120\t.\t+\t.\tID=existing_t;Parent=existing\n"
        "chr1\tbase\tCDS\t100\t120\t.\t+\t0\tParent=existing_t\n",
    )
    rankings = write(
        tmp_path / "rankings.tsv",
        "estimator\treview_rank\tcandidate_digest\n"
        f"baseline\t1\t{digest}\n"
        f"topology\t1\t{digest}\n",
    )
    isoseq = write(
        tmp_path / "isoseq.tsv",
        "candidate_digest\tevidence_state\tsupporting_transcripts"
        "\tsupporting_b73_full_length_reads\n"
        f"{digest}\tfull_chain_supported\tPB.1.1\t3\n",
    )
    paf = write(
        tmp_path / "self.paf",
        f"{digest}\t{len(cds)}\t0\t{len(cds)}\t+\tchr1\t300\t10"
        f"\t45\t{len(cds)}\t{len(cds)}\t60\tAS:i:{2 * len(cds)}\n",
    )
    return candidate, base, genome, rankings, isoseq, digest


def test_natural_audit_integrates_orf_collision_self_map_and_rna(
    tmp_path: Path,
) -> None:
    candidate, base, genome, rankings, isoseq, digest = _toy_inputs(tmp_path)
    isoseq.write_text(
        isoseq.read_text(encoding="utf-8").replace("PB.1.1", "read," * 30_000),
        encoding="utf-8",
        newline="",
    )
    paf = tmp_path / "self.paf"
    output = tmp_path / "audit.tsv"
    summary_path = tmp_path / "summary.json"

    summary = audit_natural_candidates(
        candidate_gff_path=candidate,
        base_gff_path=base,
        genome_fasta_path=genome,
        review_rankings_tsv_path=rankings,
        isoseq_evidence_tsv_path=isoseq,
        self_map_paf_path=paf,
        output_tsv_path=output,
        output_summary_json_path=summary_path,
        review_budgets=(1,),
    )

    with output.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    baseline = next(row for row in rows if row["estimator"] == "baseline")
    assert baseline["candidate_digest"] == digest
    assert baseline["complete_orf"] == "1"
    assert baseline["strict_collision_free"] == "1"
    assert baseline["mappability_class"] == "unique_locus"
    assert baseline["transcript_chain_supported"] == "1"
    assert baseline["case_study_ready"] == "1"
    assert summary["counts"]["case_study_ready"] == 1
    assert summary["review_yield"][0]["case_study_ready"] == 1
    manifest = json.loads(
        Path(str(output) + ".manifest.json").read_text(encoding="utf-8")
    )
    assert manifest["truth_access"] is True
    assert manifest["automatic_approval"] is False


def test_export_candidate_cds_is_digest_keyed_and_checksummed(tmp_path: Path) -> None:
    candidate, _, genome, _, _, digest = _toy_inputs(tmp_path)
    output = tmp_path / "candidate_cds.fa"

    manifest = export_natural_candidate_cds(
        candidate_gff_path=candidate,
        genome_fasta_path=genome,
        output_fasta_path=output,
    )

    text = output.read_text(encoding="ascii")
    assert text.startswith(f">{digest} chr1:11-45(+)\nATG")
    assert manifest["counts"] == {"candidates": 1, "cds_bases": 30}
    assert manifest["truth_access"] is False


def test_official_repeat_cds_overlap_blocks_case_study_ready(tmp_path: Path) -> None:
    candidate, base, genome, rankings, isoseq, _ = _toy_inputs(tmp_path)
    repeat = write(
        tmp_path / "repeat.gff3",
        "##gff-version 3\n"
        "source_chr1\tEDTA\tLTR_retrotransposon\t31\t40\t.\t+\t.\tID=te1\n"
        "off_target_scaffold\tEDTA\trepeat_region\t1\t20\t.\t+\t.\tID=te2\n",
    )
    seqid_map = write(
        tmp_path / "seqid_map.tsv",
        "source_seqid\ttarget_seqid\nsource_chr1\tchr1\n",
    )
    output = tmp_path / "repeat_audit.tsv"

    summary = audit_natural_candidates(
        candidate_gff_path=candidate,
        base_gff_path=base,
        genome_fasta_path=genome,
        review_rankings_tsv_path=rankings,
        isoseq_evidence_tsv_path=isoseq,
        self_map_paf_path=tmp_path / "self.paf",
        repeat_gff_path=repeat,
        repeat_seqid_map_tsv_path=seqid_map,
        output_tsv_path=output,
        output_summary_json_path=tmp_path / "repeat_summary.json",
        review_budgets=(1,),
    )

    with output.open(encoding="utf-8", newline="") as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    assert row["repeat_context_class"] == "cds_repeat_overlap"
    assert int(row["repeat_cds_overlap_bp"]) == 10
    assert row["repeat_risk_flag"] == "1"
    assert row["case_study_ready"] == "0"
    assert summary["counts"]["repeat_features"] == 1
    assert summary["counts"]["repeat_features_excluded_by_seqid_map"] == 1
    assert summary["counts"]["repeat_source_seqids_excluded_by_seqid_map"] == [
        "off_target_scaffold"
    ]


def test_minimum_raw_full_length_read_support_is_enforced(tmp_path: Path) -> None:
    candidate, base, genome, rankings, isoseq, _ = _toy_inputs(tmp_path)
    text = isoseq.read_text(encoding="utf-8").replace(
        "full_chain_supported\tPB.1.1\t3",
        "full_chain_supported\tPB.1.1\t1",
    )
    isoseq.write_text(text, encoding="utf-8", newline="")
    output = tmp_path / "minimum_reads.tsv"

    summary = audit_natural_candidates(
        candidate_gff_path=candidate,
        base_gff_path=base,
        genome_fasta_path=genome,
        review_rankings_tsv_path=rankings,
        isoseq_evidence_tsv_path=isoseq,
        self_map_paf_path=tmp_path / "self.paf",
        minimum_full_length_read_support=2,
        output_tsv_path=output,
        output_summary_json_path=tmp_path / "minimum_reads.json",
        review_budgets=(1,),
    )

    with output.open(encoding="utf-8", newline="") as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    assert row["transcript_chain_supported"] == "0"
    assert row["case_study_ready"] == "0"
    assert summary["parameters"]["minimum_full_length_read_support"] == 2
