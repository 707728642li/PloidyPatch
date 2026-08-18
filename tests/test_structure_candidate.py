from __future__ import annotations

import csv
import json
from pathlib import Path

from ploidypatch.cli import main
from ploidypatch.perturb import read_gff_document
from ploidypatch.score import build_annotation_index
from ploidypatch.structure_candidate import adapt_miniprot_structure_candidates


def write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8", newline="")
    return path


def test_structure_candidate_requires_independent_source_support(
    tmp_path: Path,
) -> None:
    annotation = write(
        tmp_path / "annotation.gff3",
        "##gff-version 3\n"
        "chr1\ttest\tgene\t1\t30\t.\t+\t.\tID=G1\n"
        "chr1\ttest\tmRNA\t1\t30\t.\t+\t.\tID=T1;Parent=G1\n"
        "chr1\ttest\tCDS\t1\t30\t.\t+\t0\tParent=T1\n",
    )
    protein_map = write(
        tmp_path / "protein_map.tsv",
        "query_id\tsource\tsource_record_id\tlength_aa\tsource_header\n"
        "a__p1\ta\tp1\t10\tp1\n"
        "a__p2\ta\tp2\t7\tp2\n"
        "b__p3\tb\tp3\t7\tp3\n"
        "a__p4\ta\tp4\t7\tp4\n",
    )
    miniprot = write(
        tmp_path / "miniprot.gff3",
        "##gff-version 3\n"
        "chr1\tminiprot\tmRNA\t1\t30\t100\t+\t.\tID=M1;Rank=1;Identity=0.9;Positive=0.9;Target=a__p1 1 10\n"
        "chr1\tminiprot\tCDS\t1\t30\t100\t+\t0\tParent=M1\n"
        "chr1\tminiprot\tmRNA\t40\t60\t90\t+\t.\tID=M2;Rank=1;Identity=0.9;Positive=0.9;Target=a__p2 1 7\n"
        "chr1\tminiprot\tCDS\t40\t60\t90\t+\t0\tParent=M2\n"
        "chr1\tminiprot\tmRNA\t40\t60\t80\t+\t.\tID=M3;Rank=1;Identity=0.9;Positive=0.9;Target=b__p3 1 7\n"
        "chr1\tminiprot\tCDS\t40\t60\t80\t+\t0\tParent=M3\n"
        "chr1\tminiprot\tmRNA\t70\t90\t70\t+\t.\tID=M4;Rank=1;Identity=0.9;Positive=0.9;Target=a__p4 1 7\n"
        "chr1\tminiprot\tCDS\t70\t90\t70\t+\t0\tParent=M4\n",
    )
    output = tmp_path / "candidates.gff3"
    decisions = tmp_path / "decisions.tsv"
    manifest = adapt_miniprot_structure_candidates(
        annotation_gff_path=annotation,
        miniprot_gff_path=miniprot,
        protein_map_path=protein_map,
        output_gff_path=output,
        decisions_tsv_path=decisions,
        min_source_support=2,
    )

    index = build_annotation_index(read_gff_document(output))
    assert len(index.transcripts) == 2
    candidate_ids = [
        transcript_id
        for transcript_id in index.transcripts
        if transcript_id.startswith("PPSC_tx_")
    ]
    assert len(candidate_ids) == 1
    assert index.transcripts[candidate_ids[0]].signature.cds == (
        (40, 60, "0"),
    )
    output_text = output.read_text(encoding="utf-8")
    assert "candidate_category=missing_annotation_candidate" in output_text
    assert "support_source_count=2" in output_text
    assert "support_sources=a,b" in output_text
    assert "support_group_count=2" in output_text
    assert "support_groups=a,b" in output_text
    assert manifest["counts"]["accepted_chain_groups"] == 1
    assert manifest["counts"]["accepted_category_counts"] == {
        "missing_annotation_candidate": 1
    }
    with decisions.open("r", encoding="utf-8", newline="") as handle:
        rows = {row["model_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
    assert rows["M1"]["reason"] == "exact_cds_chain_present"
    assert rows["M2"]["status"] == "accepted"
    assert rows["M2"]["support_source_count"] == "2"
    assert rows["M3"]["status"] == "supporting"
    assert rows["M4"]["reason"] == "source_support_below_threshold"


def test_structure_candidate_collapses_correlated_sources(
    tmp_path: Path,
) -> None:
    annotation = write(
        tmp_path / "annotation.gff3",
        "##gff-version 3\n",
    )
    protein_map = write(
        tmp_path / "protein_map.tsv",
        "query_id\tsource\tsource_record_id\tlength_aa\tsource_header\n"
        "a__p1\ta\tp1\t7\tp1\n"
        "b__p2\tb\tp2\t7\tp2\n",
    )
    source_groups = write(
        tmp_path / "source_groups.tsv",
        "source\tsupport_group\n"
        "a\tsame_species\n"
        "b\tsame_species\n",
    )
    miniprot = write(
        tmp_path / "miniprot.gff3",
        "##gff-version 3\n"
        "chr1\tminiprot\tmRNA\t40\t60\t90\t+\t.\tID=M1;Rank=1;Identity=0.9;Positive=0.9;Target=a__p1 1 7\n"
        "chr1\tminiprot\tCDS\t40\t60\t90\t+\t0\tParent=M1\n"
        "chr1\tminiprot\tmRNA\t40\t60\t80\t+\t.\tID=M2;Rank=1;Identity=0.9;Positive=0.9;Target=b__p2 1 7\n"
        "chr1\tminiprot\tCDS\t40\t60\t80\t+\t0\tParent=M2\n",
    )
    manifest = adapt_miniprot_structure_candidates(
        annotation_gff_path=annotation,
        miniprot_gff_path=miniprot,
        protein_map_path=protein_map,
        source_group_map_path=source_groups,
        output_gff_path=tmp_path / "candidate.gff3",
        decisions_tsv_path=tmp_path / "decisions.tsv",
        min_source_support=2,
    )

    assert manifest["counts"]["accepted_chain_groups"] == 0
    assert manifest["parameters"]["support_unit"] == "explicit_group"
    with (tmp_path / "decisions.tsv").open(
        "r", encoding="utf-8", newline=""
    ) as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert {row["support_source_count"] for row in rows} == {"2"}
    assert {row["support_group_count"] for row in rows} == {"1"}
    assert {row["reason"] for row in rows} == {
        "source_support_below_threshold"
    }

    cli_output = tmp_path / "cli_candidate.gff3"
    assert (
        main(
            [
                "baseline",
                "adapt-miniprot-structure",
                "--annotation-gff",
                str(annotation),
                "--miniprot-gff",
                str(miniprot),
                "--protein-map",
                str(protein_map),
                "--source-group-map",
                str(source_groups),
                "--min-source-support",
                "2",
                "--output-gff",
                str(cli_output),
                "--decisions-tsv",
                str(tmp_path / "cli_decisions.tsv"),
            ]
        )
        == 0
    )
    cli_manifest = json.loads(
        Path(str(cli_output) + ".manifest.json").read_text(encoding="utf-8")
    )
    assert cli_manifest["parameters"]["support_unit"] == "explicit_group"
