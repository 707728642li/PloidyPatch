from __future__ import annotations

import csv
from pathlib import Path

import pytest

from ploidypatch.catalog import write_missing_gene_candidate_catalog


GFF = (
    "##gff-version 3\n"
    "chr1\ttest\tgene\t1\t900\t.\t+\t.\tID=gene:G1;biotype=protein_coding\n"
    "chr1\ttest\tmRNA\t1\t900\t.\t+\t.\tID=transcript:T1;Parent=gene:G1\n"
    "chr1\ttest\texon\t1\t90\t.\t+\t.\tID=exon:E1;Parent=transcript:T1\n"
    "chr1\ttest\tCDS\t1\t90\t.\t+\t0\tID=CDS:C1;Parent=transcript:T1\n"
    "chr1\ttest\tgene\t1001\t7000\t.\t-\t.\tID=gene:G2;biotype=protein_coding\n"
    "chr1\ttest\tmRNA\t1001\t7000\t.\t-\t.\tID=transcript:T2;Parent=gene:G2\n"
    "chr1\ttest\texon\t6901\t7000\t.\t-\t.\tID=exon:E2;Parent=transcript:T2\n"
    "chr1\ttest\tCDS\t6901\t7000\t.\t-\t0\tID=CDS:C2;Parent=transcript:T2\n"
)


def write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8", newline="")
    return path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_candidate_catalog_has_structural_bins_and_external_strata(
    tmp_path: Path,
) -> None:
    source = write(tmp_path / "source.gff3", GFF)
    strata = write(
        tmp_path / "strata.tsv",
        "gene_id\tduplication_class\tsubgenome\n"
        "gene:G1\tWGD\tA\n"
        "gene:not_eligible\ttandem\t\n",
    )
    output = tmp_path / "catalog.tsv"

    manifest = write_missing_gene_candidate_catalog(source, output, strata)
    rows = read_tsv(output)

    assert len(rows) == 2
    assert rows[0]["gene_id"] == "gene:G1"
    assert rows[0]["gene_span_bin"] == "lt_1kb"
    assert rows[0]["duplication_class"] == "WGD"
    assert rows[0]["subgenome"] == "A"
    assert rows[0]["structurally_unique"] == "true"
    assert rows[1]["gene_span_bin"] == "5_to_lt_20kb"
    assert rows[1]["duplication_class"] == ""
    assert manifest["catalog"]["eligible_candidates"] == 2
    assert manifest["external_strata"]["matched_candidates"] == 1
    assert manifest["external_strata"]["external_rows_not_eligible"] == 1


def test_catalog_flags_identical_structure_assigned_to_another_gene(
    tmp_path: Path,
) -> None:
    duplicated = GFF + (
        "chr1\ttest\tgene\t1\t900\t.\t+\t.\tID=gene:G3\n"
        "chr1\ttest\tmRNA\t1\t900\t.\t+\t.\tID=transcript:T3;Parent=gene:G3\n"
        "chr1\ttest\texon\t1\t90\t.\t+\t.\tParent=transcript:T3\n"
        "chr1\ttest\tCDS\t1\t90\t.\t+\t0\tParent=transcript:T3\n"
    )
    source = write(tmp_path / "source.gff3", duplicated)
    output = tmp_path / "catalog.tsv"
    manifest = write_missing_gene_candidate_catalog(source, output)
    rows = {row["gene_id"]: row for row in read_tsv(output)}

    assert rows["gene:G1"]["structurally_unique"] == "false"
    assert rows["gene:G3"]["structurally_unique"] == "false"
    assert manifest["catalog"]["structurally_ambiguous_candidates"] == 2


def test_catalog_prefixes_external_columns_to_avoid_core_collision(
    tmp_path: Path,
) -> None:
    source = write(tmp_path / "source.gff3", GFF)
    strata = write(
        tmp_path / "strata.tsv",
        "gene_id\tseqid\tsubgenome\n"
        "gene:G1\tchr1\tA\n",
    )
    output = tmp_path / "catalog.tsv"
    manifest = write_missing_gene_candidate_catalog(
        source,
        output,
        strata,
        external_strata_prefix="wgdi_",
    )
    rows = read_tsv(output)

    assert rows[0]["seqid"] == "chr1"
    assert rows[0]["wgdi_seqid"] == "chr1"
    assert rows[0]["wgdi_subgenome"] == "A"
    assert manifest["external_strata"]["column_prefix"] == "wgdi_"


def test_catalog_refuses_overwrite_and_column_collision(tmp_path: Path) -> None:
    source = write(tmp_path / "source.gff3", GFF)
    output = tmp_path / "catalog.tsv"
    write_missing_gene_candidate_catalog(source, output)
    with pytest.raises(FileExistsError, match="Refusing to overwrite"):
        write_missing_gene_candidate_catalog(source, output)

    collision = write(
        tmp_path / "collision.tsv",
        "gene_id\tgene_span_bp\n"
        "gene:G1\t900\n",
    )
    with pytest.raises(ValueError, match="collide"):
        write_missing_gene_candidate_catalog(
            source, tmp_path / "collision_catalog.tsv", collision
        )
