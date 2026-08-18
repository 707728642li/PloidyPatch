from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from ploidypatch.wgdi_summary import (
    parse_wgdi_collinearity,
    summarize_wgdi_gene_evidence,
)


def write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8", newline="")
    return path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def make_inputs(tmp_path: Path) -> tuple[Path, Path, Path, Path]:
    query_gff = write(
        tmp_path / "query.gff",
        "chrA\tG1\t1\t100\t+\t1\n"
        "chrA\tG2\t200\t300\t+\t2\n"
        "chrC\tG3\t1\t100\t-\t1\n"
        "chrC\tG4\t200\t300\t-\t2\n",
    )
    query_manifest = write(
        tmp_path / "query.manifest.json",
        json.dumps({"chromosome_labels": {"chrA": "A1", "chrC": "C1"}}),
    )
    bra = write(
        tmp_path / "bra.collinearity.tsv",
        "# Alignment 1: score=200 pvalue=0.01 N=2 chrA&A01 plus\n"
        "G1 1 Bra1 10 1\n"
        "G2 2 Bra2 11 1\n"
        "# Alignment 2: score=100 pvalue=0.1 N=1 chrC&A02 minus\n"
        "G3 1 Bra3 20 2\n",
    )
    bol = write(
        tmp_path / "bol.collinearity.tsv",
        "# Alignment 1: score=180 pvalue=0.02 N=2 chrC&C01 plus\n"
        "G3 1 Bo1 10 1\n"
        "G4 2 Bo2 11 1\n"
        "# Alignment 2: score=90 pvalue=0.2 N=1 chrA&C02 minus\n"
        "G2 2 Bo3 20 2\n",
    )
    return query_gff, query_manifest, bra, bol


def test_parse_wgdi_collinearity_validates_declared_pair_count(
    tmp_path: Path,
) -> None:
    bad = write(
        tmp_path / "bad.tsv",
        "# Alignment 1: score=10 pvalue=0.1 N=2 chr1&chr2 plus\n"
        "G1 1 H1 1 1\n",
    )
    with pytest.raises(ValueError, match="declares 2 pairs"):
        parse_wgdi_collinearity(bad, "source")


def test_wgdi_summary_builds_expected_and_cross_support_strata(
    tmp_path: Path,
) -> None:
    query_gff, query_manifest, bra, bol = make_inputs(tmp_path)
    genes = tmp_path / "gene_evidence.tsv"
    blocks = tmp_path / "blocks.tsv"
    manifest = summarize_wgdi_gene_evidence(
        query_gff,
        query_manifest,
        [f"bra_a={bra}", f"bol_c={bol}"],
        ["A=bra_a", "C=bol_c"],
        genes,
        blocks,
    )
    rows = {row["gene_id"]: row for row in read_tsv(genes)}

    assert rows["G1"]["synteny_stratum"] == "expected_only"
    assert rows["G2"]["synteny_stratum"] == "expected_and_cross"
    assert rows["G3"]["synteny_stratum"] == "expected_and_cross"
    assert rows["G4"]["synteny_stratum"] == "expected_only"
    assert rows["G1"]["bra_a_best_block_score"] == "200.0"
    assert rows["G4"]["bol_c_counterpart_gene_count"] == "1"
    assert len(read_tsv(blocks)) == 4
    assert manifest["synteny_stratum_counts"] == {
        "expected_and_cross": 2,
        "expected_only": 2,
    }


def test_wgdi_summary_refuses_unknown_query_gene_and_overwrite(
    tmp_path: Path,
) -> None:
    query_gff, query_manifest, bra, _ = make_inputs(tmp_path)
    unknown = write(
        tmp_path / "unknown.tsv",
        "# Alignment 1: score=10 pvalue=0.1 N=1 chrA&A01 plus\n"
        "missing 1 Bra1 1 1\n",
    )
    with pytest.raises(ValueError, match="absent from query GFF"):
        summarize_wgdi_gene_evidence(
            query_gff,
            query_manifest,
            [f"bra_a={unknown}"],
            ["A=bra_a"],
            tmp_path / "unknown_genes.tsv",
            tmp_path / "unknown_blocks.tsv",
        )

    genes = tmp_path / "genes.tsv"
    blocks = tmp_path / "blocks.tsv"
    summarize_wgdi_gene_evidence(
        query_gff,
        query_manifest,
        [f"bra_a={bra}"],
        ["A=bra_a"],
        genes,
        blocks,
    )
    with pytest.raises(FileExistsError, match="Refusing to overwrite"):
        summarize_wgdi_gene_evidence(
            query_gff,
            query_manifest,
            [f"bra_a={bra}"],
            ["A=bra_a"],
            genes,
            blocks,
        )


def test_wgdi_summary_supports_nested_manifest_and_ad_subgenomes(
    tmp_path: Path,
) -> None:
    query_gff = write(
        tmp_path / "cotton.gff",
        "A01\tGA\t1\t100\t+\t1\nD01\tGD\t1\t100\t+\t1\n",
    )
    query_manifest = write(
        tmp_path / "cotton.manifest.json",
        json.dumps(
            {"selection": {"chromosome_labels": {"A01": "A01", "D01": "D01"}}}
        ),
    )
    arb = write(
        tmp_path / "arb.collinearity.tsv",
        "# Alignment 1: score=10 pvalue=0.1 N=1 A01&A01 plus\n"
        "GA 1 Arb1 1 1\n",
    )
    rai = write(
        tmp_path / "rai.collinearity.tsv",
        "# Alignment 1: score=10 pvalue=0.1 N=1 D01&D01 plus\n"
        "GD 1 Rai1 1 1\n",
    )
    genes = tmp_path / "cotton_genes.tsv"
    manifest = summarize_wgdi_gene_evidence(
        query_gff,
        query_manifest,
        [f"arb_a={arb}", f"rai_d={rai}"],
        ["A=arb_a", "D=rai_d"],
        genes,
        tmp_path / "cotton_blocks.tsv",
    )

    rows = {row["gene_id"]: row for row in read_tsv(genes)}
    assert rows["GA"]["subgenome"] == "A"
    assert rows["GD"]["subgenome"] == "D"
    assert rows["GA"]["synteny_stratum"] == "expected_only"
    assert rows["GD"]["synteny_stratum"] == "expected_only"
    assert manifest["subgenome_assignment"] == {
        "method": "longest_declared_prefix_of_chromosome_label",
        "declared_subgenomes": ["A", "D"],
    }
