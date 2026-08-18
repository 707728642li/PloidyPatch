from __future__ import annotations

import csv
from pathlib import Path

from ploidypatch.cli import main
from ploidypatch.self_wgd_pairs import infer_self_wgdi_pairs


def write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8", newline="")
    return path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_self_wgdi_pairs_require_cross_seqid_reciprocal_uniqueness(
    tmp_path: Path,
) -> None:
    gff = write(
        tmp_path / "self.gff",
        "chr1\tA1\t1\t10\t+\t1\n"
        "chr1\tA2\t20\t30\t+\t2\n"
        "chr2\tB1\t1\t10\t+\t1\n"
        "chr2\tB2\t20\t30\t+\t2\n"
        "chr3\tC2\t20\t30\t+\t1\n",
    )
    collinearity = write(
        tmp_path / "self.tsv",
        "# Alignment 1: score=100 pvalue=0.01 N=4 chr1&chr2 plus\n"
        "A1 1 B1 1 1\n"
        "A2 2 B2 2 1\n"
        "A2 2 C2 1 1\n"
        "A1 1 A1 1 1\n",
    )
    pairs = tmp_path / "pairs.tsv"
    decisions = tmp_path / "decisions.tsv"
    manifest = infer_self_wgdi_pairs(
        query_wgdi_gff_path=gff,
        collinearity_path=collinearity,
        output_pair_tsv_path=pairs,
        decisions_tsv_path=decisions,
        wgd_event="test_paleopolyploid_wgd",
        min_block_pairs=4,
    )

    assert [(row["gene_id_a"], row["gene_id_b"]) for row in read_tsv(pairs)] == [
        ("A1", "B1")
    ]
    rejected = {
        (row["gene_id_a"], row["gene_id_b"]): row["reason"]
        for row in read_tsv(decisions)
        if row["status"] == "rejected"
    }
    assert rejected == {
        ("A2", "B2"): "nonreciprocal_multiple_partners",
        ("A2", "C2"): "nonreciprocal_multiple_partners",
    }
    assert manifest["counts"]["skipped_pair_rows"]["self_pairs"] == 1


def test_self_wgdi_pair_cli_routes(tmp_path: Path) -> None:
    gff = write(
        tmp_path / "self.gff",
        "chr1\tA\t1\t10\t+\t1\n"
        "chr2\tB\t1\t10\t+\t1\n",
    )
    collinearity = write(
        tmp_path / "self.tsv",
        "# Alignment 1: score=100 pvalue=0.01 N=2 chr1&chr2 plus\n"
        "A 1 B 1 1\n"
        "B 1 A 1 1\n",
    )
    source_gff = write(
        tmp_path / "source.gff3",
        "##gff-version 3\n"
        "chr1\ttest\tgene\t1\t10\t.\t+\t.\tID=gene:A;gene_id=A\n"
        "chr2\ttest\tgene\t1\t10\t.\t+\t.\tID=gene:B;gene_id=B\n",
    )
    pairs = tmp_path / "pairs.tsv"
    assert (
        main(
            [
                "evidence",
                "infer-self-wgd-pairs",
                "--query-wgdi-gff",
                str(gff),
                "--collinearity",
                str(collinearity),
                "--source-gff",
                str(source_gff),
                "--wgd-event",
                "test_wgd",
                "--min-block-pairs",
                "2",
                "--output-pairs",
                str(pairs),
                "--decisions-tsv",
                str(tmp_path / "decisions.tsv"),
            ]
        )
        == 0
    )
    rows = read_tsv(pairs)
    assert len(rows) == 1
    assert (rows[0]["gene_id_a"], rows[0]["gene_id_b"]) == (
        "gene:A",
        "gene:B",
    )
    assert (rows[0]["wgdi_gene_id_a"], rows[0]["wgdi_gene_id_b"]) == (
        "A",
        "B",
    )
    assert rows[0]["support_block_count"] == "1"
