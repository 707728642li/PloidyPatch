from __future__ import annotations

import csv
from pathlib import Path

from ploidypatch.structure_hypothesis import infer_structure_hypotheses


def write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8", newline="")
    return path


def candidate_model(
    group: str,
    number: int,
    category: str,
    cds: tuple[tuple[int, int, str], ...],
) -> str:
    start = min(item[0] for item in cds)
    end = max(item[1] for item in cds)
    attributes = (
        f"candidate_group={group};candidate_category={category};"
        "support_source_count=2;support_sources=a,b"
    )
    lines = [
        f"chr1\tPloidyPatchCandidate\tgene\t{start}\t{end}\t.\t+\t.\t"
        f"ID=CG{number};{attributes}\n",
        f"chr1\tPloidyPatchCandidate\tmRNA\t{start}\t{end}\t.\t+\t.\t"
        f"ID=CT{number};Parent=CG{number};{attributes}\n",
    ]
    for segment_start, segment_end, phase in cds:
        lines.append(
            f"chr1\tPloidyPatchCandidate\tCDS\t{segment_start}\t"
            f"{segment_end}\t.\t+\t{phase}\tParent=CT{number}\n"
        )
    return "".join(lines)


def test_exact_topology_infers_four_structure_events(tmp_path: Path) -> None:
    annotation_text = (
        "##gff-version 3\n"
        "chr1\ttest\tgene\t1\t80\t.\t+\t.\tID=GB\n"
        "chr1\ttest\tmRNA\t1\t80\t.\t+\t.\tID=TB;Parent=GB\n"
        "chr1\ttest\tCDS\t1\t30\t.\t+\t0\tParent=TB\n"
        "chr1\ttest\tCDS\t50\t80\t.\t+\t0\tParent=TB\n"
        "chr1\ttest\tgene\t200\t280\t.\t+\t.\tID=GM\n"
        "chr1\ttest\tmRNA\t200\t280\t.\t+\t.\tID=TM;Parent=GM\n"
        "chr1\ttest\tCDS\t200\t220\t.\t+\t0\tParent=TM\n"
        "chr1\ttest\tCDS\t260\t280\t.\t+\t0\tParent=TM\n"
        "chr1\ttest\tgene\t400\t420\t.\t+\t.\tID=GS1\n"
        "chr1\ttest\tmRNA\t400\t420\t.\t+\t.\tID=TS1;Parent=GS1\n"
        "chr1\ttest\tCDS\t400\t420\t.\t+\t0\tParent=TS1\n"
        "chr1\ttest\tgene\t450\t470\t.\t+\t.\tID=GS2\n"
        "chr1\ttest\tmRNA\t450\t470\t.\t+\t.\tID=TS2;Parent=GS2\n"
        "chr1\ttest\tCDS\t450\t470\t.\t+\t0\tParent=TS2\n"
        "chr1\ttest\tgene\t600\t670\t.\t+\t.\tID=GF\n"
        "chr1\ttest\tmRNA\t600\t670\t.\t+\t.\tID=TF;Parent=GF\n"
        "chr1\ttest\tCDS\t600\t620\t.\t+\t0\tParent=TF\n"
        "chr1\ttest\tCDS\t650\t670\t.\t+\t0\tParent=TF\n"
    )
    annotation = write(tmp_path / "annotation.gff3", annotation_text)
    candidates = write(
        tmp_path / "candidates.gff3",
        annotation_text
        + "###\n"
        + candidate_model(
            "PPSC-boundary",
            1,
            "cds_extension_or_missing_segment",
            ((1, 30, "0"), (50, 95, "0")),
        )
        + candidate_model(
            "PPSC-missing",
            2,
            "cds_extension_or_missing_segment",
            ((200, 220, "0"), (230, 240, "0"), (260, 280, "0")),
        )
        + candidate_model(
            "PPSC-split",
            3,
            "spans_multiple_annotated_genes",
            ((400, 420, "0"), (450, 470, "0")),
        )
        + candidate_model(
            "PPSC-fusion-a",
            4,
            "alternative_cds_chain_within_gene",
            ((600, 620, "0"),),
        )
        + candidate_model(
            "PPSC-fusion-b",
            5,
            "alternative_cds_chain_within_gene",
            ((650, 670, "0"),),
        )
        + candidate_model(
            "PPSC-background",
            6,
            "alternative_cds_chain_within_gene",
            ((10, 30, "0"), (50, 70, "0")),
        ),
    )
    output = tmp_path / "hypotheses.tsv"
    topology = tmp_path / "topology.tsv"
    manifest = infer_structure_hypotheses(
        annotation_gff_path=annotation,
        candidate_gff_path=candidates,
        output_tsv_path=output,
        candidate_topology_tsv_path=topology,
    )

    with output.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert {row["event_type"] for row in rows} == {
        "annotation_boundary_shift",
        "annotation_fused_gene",
        "annotation_missing_internal_exon",
        "annotation_split_gene",
    }
    assert len(rows) == 4
    fusion = next(row for row in rows if row["event_type"] == "annotation_fused_gene")
    assert fusion["candidate_group_ids"] == "PPSC-fusion-a,PPSC-fusion-b"
    assert fusion["annotation_gene_ids"] == "GF"
    assert manifest["counts"]["event_type_counts"] == {
        "annotation_boundary_shift": 1,
        "annotation_fused_gene": 1,
        "annotation_missing_internal_exon": 1,
        "annotation_split_gene": 1,
    }
    with topology.open("r", encoding="utf-8", newline="") as handle:
        topology_rows = {
            row["candidate_group_id"]: row
            for row in csv.DictReader(handle, delimiter="\t")
        }
    assert topology_rows["PPSC-background"]["status"] == "no_exact_topology"
