from __future__ import annotations

import csv
from pathlib import Path

from ploidypatch.perturb import read_gff_document
from ploidypatch.score import build_annotation_index
from ploidypatch.wgd_candidate_select import (
    propagate_wgd_selection_to_conflict_pool,
    select_wgd_supported_candidates,
)


def write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8", newline="")
    return path


def test_wgd_candidate_selection_requires_existing_partner(tmp_path: Path) -> None:
    base_text = (
        "##gff-version 3\n"
        "chr1\ttest\tgene\t1\t30\t.\t+\t.\tID=G1\n"
        "chr1\ttest\tmRNA\t1\t30\t.\t+\t.\tID=T1;Parent=G1\n"
        "chr1\ttest\tCDS\t1\t30\t.\t+\t0\tParent=T1\n"
    )
    base = write(tmp_path / "base.gff3", base_text)
    digest_a = "a" * 64
    digest_b = "b" * 64
    candidate = write(
        tmp_path / "candidate.gff3",
        base_text
        + "###\n"
        + f"chr2\tPloidyPatchConsensus\tgene\t101\t130\t.\t+\t.\tID=PPCONS_gene_a;consensus_digest={digest_a}\n"
        + "chr2\tPloidyPatchConsensus\tmRNA\t101\t130\t.\t+\t.\tID=PPCONS_tx_a;Parent=PPCONS_gene_a\n"
        + "chr2\tPloidyPatchConsensus\texon\t101\t130\t.\t+\t.\tParent=PPCONS_tx_a\n"
        + "chr2\tPloidyPatchConsensus\tCDS\t101\t130\t.\t+\t0\tParent=PPCONS_tx_a\n"
        + f"chr3\tPloidyPatchConsensus\tgene\t201\t230\t.\t+\t.\tID=PPCONS_gene_b;consensus_digest={digest_b}\n"
        + "chr3\tPloidyPatchConsensus\tmRNA\t201\t230\t.\t+\t.\tID=PPCONS_tx_b;Parent=PPCONS_gene_b\n"
        + "chr3\tPloidyPatchConsensus\texon\t201\t230\t.\t+\t.\tParent=PPCONS_tx_b\n"
        + "chr3\tPloidyPatchConsensus\tCDS\t201\t230\t.\t+\t0\tParent=PPCONS_tx_b\n",
    )
    pairs = write(
        tmp_path / "pairs.tsv",
        "pair_id\tgene_id_a\tgene_id_b\tsupport_block_count\tlongest_block_pairs\n"
        "P1\tG1\tPPCONS_gene_a\t1\t30\n"
        "P2\tPPCONS_gene_b\tPPCONS_gene_a\t1\t25\n",
    )
    output = tmp_path / "selected.gff3"
    selection = tmp_path / "selection.tsv"

    manifest = select_wgd_supported_candidates(
        base_gff_path=base,
        candidate_gff_path=candidate,
        pair_tsv_path=pairs,
        output_gff_path=output,
        selection_tsv_path=selection,
    )

    # Candidate A has an existing-gene partner but also a candidate partner;
    # circular support forces abstention for both candidates.
    assert manifest["counts"]["selected_models"] == 0
    index = build_annotation_index(read_gff_document(output))
    assert set(index.gene_signatures) == {"G1"}
    with selection.open("r", encoding="utf-8", newline="") as handle:
        rows = {row["gene_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
    assert rows["PPCONS_gene_a"]["reason"] == "candidate_partner_only"
    assert rows["PPCONS_gene_b"]["reason"] == "candidate_partner_only"


def test_wgd_candidate_selection_retains_unique_base_partner(tmp_path: Path) -> None:
    base_text = "##gff-version 3\nchr1\ttest\tgene\t1\t30\t.\t+\t.\tID=G1\n"
    base = write(tmp_path / "base.gff3", base_text)
    digest = "c" * 64
    candidate = write(
        tmp_path / "candidate.gff3",
        base_text
        + "###\n"
        + f"chr2\tPloidyPatchConsensus\tgene\t101\t130\t.\t+\t.\tID=PPCONS_gene_c;consensus_digest={digest}\n"
        + "chr2\tPloidyPatchConsensus\tmRNA\t101\t130\t.\t+\t.\tID=PPCONS_tx_c;Parent=PPCONS_gene_c\n"
        + "chr2\tPloidyPatchConsensus\tCDS\t101\t130\t.\t+\t0\tParent=PPCONS_tx_c\n",
    )
    pairs = write(
        tmp_path / "pairs.tsv",
        "pair_id\tgene_id_a\tgene_id_b\tsupport_block_count\tlongest_block_pairs\n"
        "P1\tG1\tPPCONS_gene_c\t2\t40\n",
    )
    output = tmp_path / "selected.gff3"
    manifest = select_wgd_supported_candidates(
        base_gff_path=base,
        candidate_gff_path=candidate,
        pair_tsv_path=pairs,
        output_gff_path=output,
        selection_tsv_path=tmp_path / "selection.tsv",
    )

    assert manifest["counts"]["selected_models"] == 1
    index = build_annotation_index(read_gff_document(output))
    assert "PPCONS_gene_c" in index.gene_signatures


def test_wgd_partner_propagates_only_within_unambiguous_conflict_set(
    tmp_path: Path,
) -> None:
    base_text = (
        "##gff-version 3\n"
        "chr1\ttest\tgene\t1\t30\t.\t+\t.\tID=base_gene\n"
        "chr1\ttest\tmRNA\t1\t30\t.\t+\t.\tID=base_tx;Parent=base_gene\n"
        "chr1\ttest\tCDS\t1\t30\t.\t+\t0\tParent=base_tx\n"
    )
    digest_a = "a" * 64
    digest_b = "b" * 64
    digest_c = "c" * 64
    conflict = "d" * 64
    candidate_text = base_text + "###\n"
    for index, digest in enumerate((digest_a, digest_b, digest_c), start=1):
        start = 100 + index * 20
        gene = f"PPCONS_gene_{digest[:20]}"
        transcript = f"PPCONS_tx_{digest[:20]}"
        candidate_text += (
            f"chr1\tPloidyPatchConsensus\tgene\t{start}\t{start + 9}\t.\t+\t.\tID={gene};consensus_digest={digest}\n"
            f"chr1\tPloidyPatchConsensus\tmRNA\t{start}\t{start + 9}\t.\t+\t.\tID={transcript};Parent={gene};consensus_digest={digest}\n"
            f"chr1\tPloidyPatchConsensus\tCDS\t{start}\t{start + 9}\t.\t+\t0\tParent={transcript}\n"
        )
    base = write(tmp_path / "base.gff3", base_text)
    candidate_gff = write(tmp_path / "pool.gff3", candidate_text)
    decisions = write(
        tmp_path / "pool.tsv",
        "consensus_digest\tstatus\tconflict_set_digest\tconflict_member_count\n"
        f"{digest_a}\taccepted\t{conflict}\t2\n"
        f"{digest_b}\taccepted\t{conflict}\t2\n"
        f"{digest_c}\taccepted\t\t1\n",
    )
    prior = write(
        tmp_path / "prior.tsv",
        "gene_id\tconsensus_digest\tpair_id\tpartner_gene_id\t"
        "support_block_count\tlongest_block_pairs\tstatus\treason\n"
        f"PPCONS_gene_{digest_a[:20]}\t{digest_a}\tpair1\tbase_gene\t3\t40\taccepted\tprior_pass\n",
    )
    output = tmp_path / "propagated.tsv"

    manifest = propagate_wgd_selection_to_conflict_pool(
        base_gff_path=base,
        candidate_gff_path=candidate_gff,
        pool_decisions_tsv_path=decisions,
        prior_wgd_selection_tsv_path=prior,
        output_selection_tsv_path=output,
    )

    with output.open("r", encoding="utf-8", newline="") as handle:
        rows = {
            row["consensus_digest"]: row
            for row in csv.DictReader(handle, delimiter="\t")
        }
    assert rows[digest_a]["evidence_origin"] == "prior_exact_candidate"
    assert rows[digest_b]["status"] == "accepted"
    assert rows[digest_b]["partner_gene_id"] == "base_gene"
    assert rows[digest_b]["source_candidate_digest"] == digest_a
    assert rows[digest_c]["status"] == "rejected"
    assert manifest["counts"]["inherited_accepted"] == 1
