from __future__ import annotations

import csv
from pathlib import Path

import pytest

from ploidypatch.cli import main
from ploidypatch.homeolog_pairs import (
    infer_outgroup_duplicated_pairs,
    infer_wgdi_homeolog_pairs,
)


def write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8", newline="")
    return path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def collinearity(path: Path, rows: list[tuple[str, str]]) -> Path:
    text = (
        f"# Alignment 1: score=100 pvalue=0.01 N={len(rows)} chrA&ref plus\n"
    )
    for number, (query, target) in enumerate(rows, start=1):
        text += f"{query} {number} {target} {number} 1\n"
    return write(path, text)


def make_inputs(tmp_path: Path) -> tuple[Path, Path, Path]:
    genes = write(
        tmp_path / "genes.tsv",
        "gene_id\tsubgenome\n"
        "A1\tA\n"
        "C1\tC\n"
        "A2\tA\n"
        "C2\tC\n"
        "C3\tC\n",
    )
    source_one = collinearity(
        tmp_path / "one.tsv",
        [
            ("A1", "X1"),
            ("C1", "X1"),
            ("A2", "X2"),
            ("C2", "X2"),
            ("A2", "X3"),
            ("C3", "X3"),
        ],
    )
    source_two = collinearity(
        tmp_path / "two.tsv",
        [
            ("A1", "Y1"),
            ("C1", "Y1"),
            ("A2", "Y2"),
            ("C2", "Y2"),
            ("A2", "Y3"),
            ("C3", "Y3"),
        ],
    )
    return genes, source_one, source_two


def test_infer_homeolog_pairs_requires_independent_reciprocal_support(
    tmp_path: Path,
) -> None:
    genes, source_one, source_two = make_inputs(tmp_path)
    pairs = tmp_path / "pairs.tsv"
    decisions = tmp_path / "decisions.tsv"
    manifest = infer_wgdi_homeolog_pairs(
        gene_evidence_tsv_path=genes,
        collinearity_inputs=(f"one={source_one}", f"two={source_two}"),
        output_pair_tsv_path=pairs,
        decisions_tsv_path=decisions,
        wgd_event="test_allopolyploid_AC",
        subgenomes=("A", "C"),
    )

    accepted = read_tsv(pairs)
    assert [(row["gene_id_a"], row["gene_id_b"]) for row in accepted] == [
        ("A1", "C1")
    ]
    assert accepted[0]["support_group_count"] == "2"
    rejected = {
        (row["gene_id_a"], row["gene_id_b"]): row["reason"]
        for row in read_tsv(decisions)
        if row["status"] == "rejected"
    }
    assert rejected == {
        ("A2", "C2"): "nonreciprocal_multiple_partners",
        ("A2", "C3"): "nonreciprocal_multiple_partners",
    }
    assert manifest["counts"]["accepted_pairs"] == 1
    assert manifest["parameters"]["support_unit"] == "source"


def test_homeolog_pairs_collapse_correlated_sources_and_cli_routes(
    tmp_path: Path,
) -> None:
    genes, source_one, source_two = make_inputs(tmp_path)
    groups = write(
        tmp_path / "groups.tsv",
        "source\tsupport_group\n"
        "one\tcorrelated\n"
        "two\tcorrelated\n",
    )
    pairs = tmp_path / "pairs.tsv"
    assert (
        main(
            [
                "evidence",
                "infer-homeolog-pairs",
                "--gene-evidence",
                str(genes),
                "--collinearity",
                f"one={source_one}",
                "--collinearity",
                f"two={source_two}",
                "--source-group-map",
                str(groups),
                "--wgd-event",
                "test_allopolyploid_AC",
                "--subgenome",
                "A",
                "--subgenome",
                "C",
                "--min-support-group-count",
                "2",
                "--output-pairs",
                str(pairs),
                "--decisions-tsv",
                str(tmp_path / "decisions.tsv"),
            ]
        )
        == 0
    )
    assert read_tsv(pairs) == []
    assert all(
        row["reason"] == "support_below_threshold"
        for row in read_tsv(tmp_path / "decisions.tsv")
    )

    with pytest.raises(FileExistsError, match="Refusing to overwrite"):
        infer_wgdi_homeolog_pairs(
            gene_evidence_tsv_path=genes,
            collinearity_inputs=(f"one={source_one}", f"two={source_two}"),
            output_pair_tsv_path=pairs,
            decisions_tsv_path=tmp_path / "unused.tsv",
            wgd_event="test_allopolyploid_AC",
            subgenomes=("A", "C"),
        )


def test_outgroup_duplicated_pairs_require_matching_one_to_two_support(
    tmp_path: Path,
) -> None:
    query = write(
        tmp_path / "query.gff",
        "chr1\tM1\t1\t100\t+\t1\n"
        "chr2\tM2\t1\t100\t+\t1\n"
        "chr3\tN1\t1\t100\t+\t1\n"
        "chr4\tN2\t1\t100\t+\t1\n"
        "chr5\tN3\t1\t100\t+\t1\n"
        "chr6\tP1\t1\t100\t+\t1\n"
        "chr6\tP2\t201\t300\t+\t2\n"
        "chr7\tQ1\t1\t100\t+\t1\n"
        "chr8\tQ2\t1\t100\t+\t1\n"
        "chr9\tQ3\t1\t100\t+\t1\n",
    )
    rows_one = [
        ("M1", "S1"),
        ("M2", "S1"),
        ("N1", "S2"),
        ("N2", "S2"),
        ("N1", "S3"),
        ("N3", "S3"),
        ("P1", "S4"),
        ("P2", "S4"),
        ("Q1", "S5"),
        ("Q2", "S5"),
        ("Q3", "S5"),
    ]
    rows_two = [
        (query_gene, target.replace("S", "T"))
        for query_gene, target in rows_one
    ]
    source_one = collinearity(tmp_path / "sorghum.tsv", rows_one)
    source_two = collinearity(tmp_path / "setaria.tsv", rows_two)
    pairs = tmp_path / "pairs.tsv"
    decisions = tmp_path / "decisions.tsv"
    manifest = infer_outgroup_duplicated_pairs(
        query_wgdi_gff_path=query,
        collinearity_inputs=(
            f"sorghum={source_one}",
            f"setaria={source_two}",
        ),
        output_pair_tsv_path=pairs,
        decisions_tsv_path=decisions,
        wgd_event="maize1",
        min_block_pairs=2,
    )

    accepted = read_tsv(pairs)
    assert [(row["gene_id_a"], row["gene_id_b"]) for row in accepted] == [
        ("M1", "M2")
    ]
    assert accepted[0]["seqid_a"] == "chr1"
    assert accepted[0]["seqid_b"] == "chr2"
    rejected = {
        (row["gene_id_a"], row["gene_id_b"]): row["reason"]
        for row in read_tsv(decisions)
        if row["status"] == "rejected"
    }
    assert rejected == {
        ("N1", "N2"): "nonreciprocal_multiple_partners",
        ("N1", "N3"): "nonreciprocal_multiple_partners",
    }
    assert manifest["counts"]["same_seqid_counterpart_anchors_rejected"] == 2
    assert manifest["counts"]["promiscuous_counterpart_anchors_rejected"] == 2
    assert manifest["counts"]["accepted_pairs"] == 1


def test_outgroup_duplicated_pairs_cli_routes(tmp_path: Path) -> None:
    query = write(
        tmp_path / "query.gff",
        "1\tM1\t1\t100\t+\t1\n2\tM2\t1\t100\t+\t1\n",
    )
    source_one = collinearity(tmp_path / "one.tsv", [("M1", "X"), ("M2", "X")])
    source_two = collinearity(tmp_path / "two.tsv", [("M1", "Y"), ("M2", "Y")])
    source_gff = write(
        tmp_path / "source.gff3",
        "1\ttest\tgene\t1\t100\t.\t+\t.\tID=gene:M1;gene_id=M1\n"
        "2\ttest\tgene\t1\t100\t.\t+\t.\tID=gene:M2;gene_id=M2\n",
    )
    output = tmp_path / "pairs.tsv"
    assert main(
        [
            "evidence",
            "infer-outgroup-duplicated-pairs",
            "--query-wgdi-gff",
            str(query),
            "--source-gff",
            str(source_gff),
            "--collinearity",
            f"one={source_one}",
            "--collinearity",
            f"two={source_two}",
            "--wgd-event",
            "maize1",
            "--min-block-pairs",
            "2",
            "--output-pairs",
            str(output),
            "--decisions-tsv",
            str(tmp_path / "decisions.tsv"),
        ]
    ) == 0
    rows = read_tsv(output)
    assert len(rows) == 1
    assert rows[0]["gene_id_a"] == "gene:M1"
    assert rows[0]["gene_id_b"] == "gene:M2"
    assert rows[0]["wgdi_gene_id_a"] == "M1"
