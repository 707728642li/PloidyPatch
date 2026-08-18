from __future__ import annotations

import csv
from pathlib import Path

import pytest

from ploidypatch.cli import main
from ploidypatch.homeolog_topology import build_homeolog_topology_features


def _write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8", newline="")
    return path


def test_build_homeolog_topology_features_is_truth_blind_and_isoform_aware(
    tmp_path: Path,
) -> None:
    features = _write(
        tmp_path / "features.tsv",
        "candidate_digest\tsupport_method_count\nD1\t2\nD2\t1\n",
    )
    selection = _write(
        tmp_path / "selection.tsv",
        "gene_id\tconsensus_digest\tpair_id\tpartner_gene_id\t"
        "support_block_count\tlongest_block_pairs\tstatus\treason\n"
        "Cgene\tD1\tPAIR1\tPgene\t2\t25\taccepted\t"
        "reciprocal_wgd_partner_is_existing_gene\n"
        "Rgene\tD2\t\t\t\t\trejected\tno_accepted_wgd_partner\n",
    )
    candidate = _write(
        tmp_path / "candidate.gff3",
        "##gff-version 3\n"
        "chr1\ttest\tgene\t1\t210\t.\t+\t.\tID=Cgene\n"
        "chr1\ttest\tmRNA\t1\t210\t.\t+\t.\tID=Ctx;Parent=Cgene\n"
        "chr1\ttest\tCDS\t1\t90\t.\t+\t0\tID=Ccds1;Parent=Ctx\n"
        "chr1\ttest\tCDS\t151\t210\t.\t+\t1\tID=Ccds2;Parent=Ctx\n",
    )
    base = _write(
        tmp_path / "base.gff3",
        "##gff-version 3\n"
        "chr2\ttest\tgene\t1001\t1799\t.\t+\t.\tID=Pgene\n"
        "chr2\ttest\tmRNA\t1001\t1360\t.\t+\t.\tID=Ptx1;Parent=Pgene\n"
        "chr2\ttest\tCDS\t1001\t1180\t.\t+\t0\tID=P1;Parent=Ptx1\n"
        "chr2\ttest\tCDS\t1241\t1360\t.\t+\t1\tID=P2;Parent=Ptx1\n"
        "chr2\ttest\tmRNA\t1500\t1799\t.\t+\t.\tID=Ptx2;Parent=Pgene\n"
        "chr2\ttest\tCDS\t1500\t1799\t.\t+\t0\tID=P3;Parent=Ptx2\n",
    )
    output = tmp_path / "topology.tsv"

    manifest = build_homeolog_topology_features(
        copy_feature_tsv_path=features,
        wgd_selection_tsv_path=selection,
        candidate_gff_path=candidate,
        base_gff_path=base,
        output_tsv_path=output,
    )

    assert manifest["truth_access"] is False
    assert manifest["counts"] == {
        "candidates": 2,
        "topology_available": 1,
        "topology_unavailable": 1,
    }
    with output.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    accepted, rejected = rows
    assert accepted["candidate_transcript_id"] == "Ctx"
    assert accepted["partner_transcript_id"] == "Ptx1"
    assert float(accepted["cds_bp_ratio"]) == pytest.approx(0.5)
    assert float(accepted["cds_segment_count_ratio"]) == pytest.approx(1.0)
    assert float(accepted["phase_lcs_similarity"]) == pytest.approx(1.0)
    assert float(accepted["junction_fraction_similarity"]) == pytest.approx(1.0)
    assert float(accepted["coding_span_ratio"]) == pytest.approx(210 / 360)
    assert float(accepted["topology_coherence_score"]) == pytest.approx(
        (0.5 + 1 + 1 + 1 + 210 / 360) / 5
    )
    assert rejected["topology_available"] == "0"
    assert rejected["topology_coherence_score"] == ""

    with pytest.raises(FileExistsError, match="Refusing to overwrite"):
        build_homeolog_topology_features(
            copy_feature_tsv_path=features,
            wgd_selection_tsv_path=selection,
            candidate_gff_path=candidate,
            base_gff_path=base,
            output_tsv_path=output,
        )


def test_homeolog_topology_cli_routes(tmp_path: Path) -> None:
    features = _write(tmp_path / "features.tsv", "candidate_digest\nD1\n")
    selection = _write(
        tmp_path / "selection.tsv",
        "gene_id\tconsensus_digest\tpair_id\tpartner_gene_id\t"
        "support_block_count\tlongest_block_pairs\tstatus\treason\n"
        "C1\tD1\t\t\t\t\trejected\tno_accepted_wgd_partner\n",
    )
    empty_gff = _write(tmp_path / "empty.gff3", "##gff-version 3\n")
    output = tmp_path / "output.tsv"

    assert main(
        [
            "evidence",
            "build-homeolog-topology-features",
            "--copy-features",
            str(features),
            "--wgd-selection",
            str(selection),
            "--candidate-gff",
            str(empty_gff),
            "--base-gff",
            str(empty_gff),
            "--output-tsv",
            str(output),
        ]
    ) == 0
    assert output.exists()
