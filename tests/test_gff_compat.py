from __future__ import annotations

import csv
from pathlib import Path

import pytest

from ploidypatch.gff_compat import (
    synthesize_missing_transcript_exons,
    synthesize_root_transcript_genes,
)


def test_synthesizes_stable_gene_parents_for_root_transcripts(tmp_path: Path) -> None:
    source = tmp_path / "source.gff3"
    source.write_text(
        "chr1\ttest\tmRNA\t10\t90\t.\t+\t.\tID=T1;Name=one;\n"
        "chr1\ttest\tCDS\t10\t90\t.\t+\t0\tID=C1;Parent=T1;\n",
        encoding="utf-8",
        newline="",
    )
    output = tmp_path / "adapted.gff3"
    mapping = tmp_path / "mapping.tsv"

    report = synthesize_root_transcript_genes(source, output, mapping)

    lines = output.read_text(encoding="utf-8").splitlines()
    assert lines[0].split("\t")[2] == "gene"
    gene_id = lines[0].split("ID=", 1)[1].split(";", 1)[0]
    assert lines[1].endswith(f"Parent={gene_id};")
    assert lines[2].endswith("ID=C1;Parent=T1;")
    with mapping.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert rows == [{"transcript_id": "T1", "synthesized_gene_id": gene_id}]
    assert report["counts"]["synthesized_gene_records"] == 1


def test_refuses_annotations_that_already_have_gene_records(tmp_path: Path) -> None:
    source = tmp_path / "source.gff3"
    source.write_text(
        "chr1\ttest\tgene\t1\t10\t.\t+\t.\tID=G1\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="zero gene records"):
        synthesize_root_transcript_genes(
            source, tmp_path / "adapted.gff3", tmp_path / "mapping.tsv"
        )


def test_groups_gffread_transcripts_by_gene_id(tmp_path: Path) -> None:
    source = tmp_path / "gffread.gff3"
    source.write_text(
        "chr1\ttest\tmRNA\t10\t40\t.\t+\t.\tID=T1;geneID=G1\n"
        "chr1\ttest\texon\t10\t40\t.\t+\t.\tParent=T1\n"
        "chr1\ttest\tmRNA\t20\t90\t.\t+\t.\tID=T2;geneID=G1\n"
        "chr1\ttest\texon\t20\t90\t.\t+\t.\tParent=T2\n",
        encoding="utf-8",
        newline="",
    )
    output = tmp_path / "adapted.gff3"
    report = synthesize_root_transcript_genes(
        source, output, tmp_path / "mapping.tsv"
    )

    lines = output.read_text(encoding="utf-8").splitlines()
    genes = [line for line in lines if line.split("\t")[2] == "gene"]
    transcripts = [line for line in lines if line.split("\t")[2] == "mRNA"]
    assert len(genes) == 1
    assert genes[0].split("\t")[3:5] == ["10", "90"]
    assert "ID=G1;" in genes[0]
    assert all(line.endswith(";Parent=G1;") for line in transcripts)
    assert report["counts"] == {
        "input_gene_records": 0,
        "root_transcripts": 2,
        "synthesized_gene_records": 1,
    }


def test_synthesizes_missing_exons_without_changing_cds(tmp_path: Path) -> None:
    source = tmp_path / "source.gff3"
    source.write_text(
        "##gff-version 3\n"
        "chr1\ttest\tgene\t10\t100\t.\t+\t.\tID=G1\n"
        "chr1\ttest\tmRNA\t10\t100\t.\t+\t.\tID=T1;Parent=G1\n"
        "chr1\ttest\tfive_prime_UTR\t10\t19\t.\t+\t.\tParent=T1\n"
        "chr1\ttest\tCDS\t20\t40\t.\t+\t0\tID=C1;Parent=T1\n"
        "chr1\ttest\tCDS\t70\t90\t.\t+\t2\tID=C2;Parent=T1\n"
        "chr1\ttest\tthree_prime_UTR\t91\t100\t.\t+\t.\tParent=T1\n",
        encoding="utf-8",
        newline="",
    )
    output = tmp_path / "adapted.gff3"

    report = synthesize_missing_transcript_exons(source, output)

    text = output.read_text(encoding="utf-8")
    assert "chr1\ttest\tCDS\t20\t40\t.\t+\t0\tID=C1;Parent=T1\n" in text
    exons = [line.split("\t") for line in text.splitlines() if "\texon\t" in line]
    assert [(row[3], row[4]) for row in exons] == [("10", "40"), ("70", "100")]
    assert all("Parent=T1;" in row[8] for row in exons)
    exon_ids = [row[8].split("ID=", 1)[1].split(";", 1)[0] for row in exons]
    assert all(exon_id.startswith("PPX") and "-" not in exon_id for exon_id in exon_ids)
    assert len(set(exon_ids)) == len(exon_ids)
    assert report["coordinate_or_cds_changes"] is False
    assert report["cds_rows_sha256"]["input"] == report["cds_rows_sha256"]["output"]
    assert report["counts"]["transcripts_with_synthesized_exons"] == 1
    assert report["counts"]["unresolved_transcripts"] == 0


def test_preserves_transcripts_that_already_have_exons(tmp_path: Path) -> None:
    source = tmp_path / "source.gff3"
    source.write_text(
        "chr1\ttest\tmRNA\t1\t30\t.\t-\t.\tID=T1\n"
        "chr1\ttest\texon\t1\t30\t.\t-\t.\tID=E1;Parent=T1\n"
        "chr1\ttest\tCDS\t4\t30\t.\t-\t0\tID=C1;Parent=T1\n",
        encoding="utf-8",
        newline="",
    )
    output = tmp_path / "adapted.gff3"

    report = synthesize_missing_transcript_exons(source, output)

    assert output.read_bytes() == source.read_bytes()
    assert report["counts"]["transcripts_with_existing_exons"] == 1
    assert report["counts"]["synthesized_exon_records"] == 0


def test_reports_unresolved_transcript_without_guessing_span(tmp_path: Path) -> None:
    source = tmp_path / "source.gff3"
    source.write_text(
        "chr1\ttest\tmRNA\t1\t30\t.\t+\t.\tID=T1\n",
        encoding="utf-8",
    )

    report = synthesize_missing_transcript_exons(source, tmp_path / "out.gff3")

    assert report["counts"]["unresolved_transcripts"] == 1
    assert report["unresolved_transcript_ids"] == ["T1"]


def test_out_of_bounds_child_requires_explicit_parent_repair(tmp_path: Path) -> None:
    source = tmp_path / "provider.gff3"
    source.write_text(
        "chr5\ttest\tgene\t100\t300\t.\t+\t.\tID=G1\n"
        "chr5\ttest\tmRNA\t100\t300\t.\t+\t.\tID=T1;Parent=G1\n"
        "chr5\ttest\tCDS\t80\t150\t.\t+\t0\tID=C1;Parent=T1\n",
        encoding="utf-8",
        newline="",
    )
    with pytest.raises(ValueError, match="outside parent bounds"):
        synthesize_missing_transcript_exons(source, tmp_path / "strict.gff3")


def test_explicit_parent_repair_expands_only_gene_and_transcript_bounds(
    tmp_path: Path,
) -> None:
    source = tmp_path / "provider.gff3"
    unchanged_children = (
        "chr5\ttest\tfive_prime_UTR\t75\t79\t.\t+\t.\tID=U1;Parent=T1\n"
        "chr5\ttest\tCDS\t80\t150\t7.5\t+\t2\tID=C1;Parent=T1\n"
        "chr5\ttest\tCDS\t280\t325\t.\t+\t0\tID=C2;Parent=T1\n"
    )
    source.write_text(
        "##gff-version 3\n"
        "chr5\ttest\tgene\t100\t300\t.\t+\t.\tID=G1\n"
        "chr5\ttest\tmRNA\t100\t300\t.\t+\t.\tID=T1;Parent=G1\n"
        + unchanged_children,
        encoding="utf-8",
        newline="",
    )
    output = tmp_path / "repaired.gff3"
    report = synthesize_missing_transcript_exons(
        source, output, repair_parent_bounds=True
    )
    text = output.read_text(encoding="utf-8")
    assert "chr5\ttest\tgene\t75\t325\t.\t+\t.\tID=G1\n" in text
    assert "chr5\ttest\tmRNA\t75\t325\t.\t+\t.\tID=T1;Parent=G1\n" in text
    for child in unchanged_children.splitlines(keepends=True):
        assert child in text
    assert report["schema_version"] == "ploidypatch.missing_transcript_exon_compat.v2"
    assert report["repair_parent_bounds"] is True
    assert report["coordinate_or_cds_changes"] is True
    assert report["child_coordinate_or_cds_changes"] is False
    assert report["counts"]["transcript_parent_bounds_repaired"] == 1
    assert report["counts"]["gene_parent_bounds_repaired"] == 1
    assert report["counts"]["parent_bounds_repaired"] == 2
    assert report["counts"]["protected_input_child_records_verified_unchanged"] == 3
    assert report["cds_rows_sha256"]["input"] == report["cds_rows_sha256"]["output"]
    assert {item["feature_type"] for item in report["parent_bound_repairs"]} == {
        "gene",
        "transcript",
    }


@pytest.mark.parametrize(
    ("child_seqid", "child_strand", "message"),
    [
        ("chr6", "+", "changes seqid"),
        ("chr5", "-", "changes strand"),
    ],
)
def test_parent_repair_rejects_cross_seqid_or_strand_children(
    tmp_path: Path, child_seqid: str, child_strand: str, message: str
) -> None:
    source = tmp_path / "provider.gff3"
    source.write_text(
        "chr5\ttest\tgene\t100\t300\t.\t+\t.\tID=G1\n"
        "chr5\ttest\tmRNA\t100\t300\t.\t+\t.\tID=T1;Parent=G1\n"
        f"{child_seqid}\ttest\tCDS\t80\t150\t.\t{child_strand}\t0\t"
        "ID=C1;Parent=T1\n",
        encoding="utf-8",
        newline="",
    )
    output = tmp_path / "repaired.gff3"
    with pytest.raises(ValueError, match=message):
        synthesize_missing_transcript_exons(
            source, output, repair_parent_bounds=True
        )
    assert not output.exists()


def test_parent_repair_rejects_unresolved_structural_child_parent(tmp_path: Path) -> None:
    source = tmp_path / "provider.gff3"
    source.write_text(
        "chr5\ttest\tgene\t100\t300\t.\t+\t.\tID=G1\n"
        "chr5\ttest\tmRNA\t100\t300\t.\t+\t.\tID=T1;Parent=G1\n"
        "chr5\ttest\tCDS\t100\t150\t.\t+\t0\tID=C1;Parent=missing\n",
        encoding="utf-8",
        newline="",
    )
    with pytest.raises(ValueError, match="not a transcript"):
        synthesize_missing_transcript_exons(
            source, tmp_path / "repaired.gff3", repair_parent_bounds=True
        )


def test_parent_repair_preserves_unique_exon_only_feature_containers(
    tmp_path: Path,
) -> None:
    source = tmp_path / "provider.gff3"
    body = (
        "chr5\ttest\tgene\t100\t300\t.\t+\t.\tID=G1\n"
        "chr5\ttest\tlnc_RNA\t100\t300\t.\t+\t.\tID=R1;Parent=G1\n"
        "chr5\ttest\texon\t110\t290\t.\t+\t.\tID=E1;Parent=R1\n"
        "chr5\ttest\tmRNA\t120\t280\t.\t+\t.\tID=T1;Parent=G1\n"
        "chr5\ttest\texon\t120\t280\t.\t+\t.\tID=E2;Parent=T1\n"
        "chr5\ttest\tpseudogene\t400\t500\t.\t-\t.\tID=P1\n"
        "chr5\ttest\texon\t420\t480\t.\t-\t.\tID=E3;Parent=P1\n"
        "chr5\ttest\tpseudogene\t600\t800\t.\t+\t.\tID=P2\n"
        "chr5\ttest\ttranscript\t620\t780\t.\t+\t.\tID=T2;Parent=P2\n"
        "chr5\ttest\texon\t630\t770\t.\t+\t.\tID=E4;Parent=T2\n"
    )
    source.write_text(
        body,
        encoding="utf-8",
        newline="",
    )
    output = tmp_path / "repaired.gff3"
    report = synthesize_missing_transcript_exons(
        source, output, repair_parent_bounds=True
    )
    observed = output.read_text(encoding="utf-8")
    assert observed == body
    assert report["promoted_transcript_like_feature_types"] == {}
    assert report["promoted_gene_like_root_feature_types"] == {"pseudogene": 1}
    assert report["preserved_exon_only_feature_container_types"] == {
        "lnc_RNA": 1,
        "pseudogene": 1,
    }
    assert report["counts"]["promoted_transcript_like_features"] == 0
    assert report["counts"]["promoted_gene_like_root_features"] == 1
    assert report["counts"]["preserved_exon_only_feature_containers"] == 2
    assert report["counts"]["transcripts_with_existing_exons"] == 2
    assert report["counts"]["synthesized_exon_records"] == 0
    assert report["counts"]["parent_bounds_repaired"] == 0
    assert report["child_coordinate_or_cds_changes"] is False


def test_parent_repair_preserves_existing_exon_row_byte_for_byte(tmp_path: Path) -> None:
    source = tmp_path / "provider.gff3"
    exon = "chr5\ttest\texon\t90\t120\t9.0\t-\t.\tID=E1;Parent=T1;Note=keep_me\n"
    source.write_text(
        "chr5\ttest\tgene\t100\t300\t.\t-\t.\tID=G1\n"
        "chr5\ttest\tmRNA\t100\t300\t.\t-\t.\tID=T1;Parent=G1\n"
        + exon
        + "chr5\ttest\tCDS\t90\t120\t.\t-\t0\tID=C1;Parent=T1\n",
        encoding="utf-8",
        newline="",
    )
    output = tmp_path / "repaired.gff3"
    report = synthesize_missing_transcript_exons(
        source, output, repair_parent_bounds=True
    )
    assert output.read_text(encoding="utf-8").count(exon) == 1
    assert report["counts"]["transcripts_with_existing_exons"] == 1
    assert report["counts"]["synthesized_exon_records"] == 0
    assert report["cds_rows_sha256"]["input"] == report["cds_rows_sha256"]["output"]


def test_eriantha_dtz79_parent_bound_regression(tmp_path: Path) -> None:
    source = tmp_path / "eriantha.gff3"
    cds = (
        "DTZ79_05\tmaker\tCDS\t2960679\t2960853\t.\t+\t0\t"
        "ID=DTZ79_05g01990:2960679-2960853;Parent=DTZ79_05g01990\n"
    )
    source.write_text(
        "DTZ79_05\tmaker\tgene\t2960873\t2966950\t.\t+\t.\t"
        "ID=gene:DTZ79_05g01990\n"
        "DTZ79_05\tmaker\tmRNA\t2960873\t2966950\t.\t+\t.\t"
        "ID=DTZ79_05g01990;Parent=gene:DTZ79_05g01990\n"
        + cds,
        encoding="utf-8",
        newline="",
    )
    output = tmp_path / "repaired.gff3"
    report = synthesize_missing_transcript_exons(
        source, output, repair_parent_bounds=True
    )
    text = output.read_text(encoding="utf-8")
    assert "\tgene\t2960679\t2966950\t" in text
    assert "\tmRNA\t2960679\t2966950\t" in text
    assert cds in text
    assert report["counts"]["parent_bounds_repaired"] == 2
