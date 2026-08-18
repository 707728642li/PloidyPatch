from __future__ import annotations

import csv
from pathlib import Path

from ploidypatch.patch import (
    apply_annotation_patch,
    create_annotation_patch,
    revert_annotation_patch,
)
from ploidypatch.perturb import read_gff_document
from ploidypatch.score import build_annotation_index
from ploidypatch.structure_patch import compile_structure_patch_edits


def write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8", newline="")
    return path


def test_compile_split_and_fusion_hypotheses_to_reversible_patch(
    tmp_path: Path,
) -> None:
    annotation_text = (
        "##gff-version 3\n"
        "chr1\ttest\tgene\t1\t20\t.\t+\t.\tID=G1\n"
        "chr1\ttest\tmRNA\t1\t20\t.\t+\t.\tID=T1;Parent=G1\n"
        "chr1\ttest\tCDS\t1\t10\t.\t+\t0\tID=CDS:T1;Parent=T1\n"
        "chr1\ttest\tCDS\t11\t20\t.\t+\t0\tID=CDS:T1;Parent=T1\n"
        "chr1\ttest\tgene\t40\t60\t.\t+\t.\tID=G2\n"
        "chr1\ttest\tmRNA\t40\t60\t.\t+\t.\tID=T2;Parent=G2\n"
        "chr1\ttest\tCDS\t40\t60\t.\t+\t0\tParent=T2\n"
        "chr1\ttest\tgene\t100\t160\t.\t+\t.\tID=GF\n"
        "chr1\ttest\tmRNA\t100\t160\t.\t+\t.\tID=TF;Parent=GF\n"
        "chr1\ttest\tCDS\t100\t120\t.\t+\t0\tID=CDS:TF;Parent=TF\n"
        "chr1\ttest\tCDS\t140\t160\t.\t+\t0\tID=CDS:TF;Parent=TF\n"
    )
    annotation = write(tmp_path / "annotation.gff3", annotation_text)
    candidate_text = annotation_text + "###\n"
    candidate_rows = (
        ("S", "split", ((1, 10, "0"), (11, 20, "0"), (40, 60, "0"))),
        ("FA", "alternative", ((100, 120, "0"),)),
        ("FB", "alternative", ((140, 160, "0"),)),
    )
    for number, (group, category, cds) in enumerate(candidate_rows, start=1):
        start = min(item[0] for item in cds)
        end = max(item[1] for item in cds)
        attrs = (
            f"candidate_group={group};candidate_category={category};"
            "support_source_count=2;support_sources=a,b;"
            "support_group_count=2;support_groups=ga,gb"
        )
        candidate_text += (
            f"chr1\tPloidyPatchCandidate\tgene\t{start}\t{end}\t.\t+\t.\t"
            f"ID=CG{number};{attrs}\n"
            f"chr1\tPloidyPatchCandidate\tmRNA\t{start}\t{end}\t.\t+\t.\t"
            f"ID=CT{number};Parent=CG{number};{attrs}\n"
        )
        for segment_start, segment_end, phase in cds:
            candidate_text += (
                f"chr1\tPloidyPatchCandidate\tCDS\t{segment_start}\t"
                f"{segment_end}\t.\t+\t{phase}\tParent=CT{number}\n"
            )
    candidates = write(tmp_path / "candidates.gff3", candidate_text)
    hypotheses = tmp_path / "hypotheses.tsv"
    fields = (
        "hypothesis_id",
        "event_type",
        "candidate_group_ids",
        "annotation_gene_ids",
        "annotation_transcript_ids",
        "support_group_count_min",
    )
    with hypotheses.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerow(
            {
                "hypothesis_id": "PPH-split",
                "event_type": "annotation_split_gene",
                "candidate_group_ids": "S",
                "annotation_gene_ids": "G1,G2",
                "annotation_transcript_ids": "T1,T2",
                "support_group_count_min": 2,
            }
        )
        writer.writerow(
            {
                "hypothesis_id": "PPH-fusion",
                "event_type": "annotation_fused_gene",
                "candidate_group_ids": "FA,FB",
                "annotation_gene_ids": "GF",
                "annotation_transcript_ids": "TF",
                "support_group_count_min": 2,
            }
        )

    edits = tmp_path / "edits.json"
    report = compile_structure_patch_edits(
        annotation_gff_path=annotation,
        candidate_gff_path=candidates,
        hypotheses_tsv_path=hypotheses,
        output_edits_json_path=edits,
        allowed_event_types=("annotation_split_gene", "annotation_fused_gene"),
        min_support_group_count=2,
    )
    assert report["counts"]["compiled_hypotheses"] == 2
    patch = tmp_path / "patch.json"
    create_annotation_patch(annotation, edits, patch)
    patched = tmp_path / "patched.gff3"
    apply_annotation_patch(annotation, patch, patched)
    index = build_annotation_index(read_gff_document(patched))
    assert set(index.transcripts) == {
        "PPR_tx_split",
        "PPR_tx_fusion_1",
        "PPR_tx_fusion_2",
    }
    assert {model.signature.cds for model in index.transcripts.values()} == {
        ((1, 10, "0"), (11, 20, "0"), (40, 60, "0")),
        ((100, 120, "0"),),
        ((140, 160, "0"),),
    }
    reverted = tmp_path / "reverted.gff3"
    revert_annotation_patch(patched, patch, reverted)
    assert reverted.read_bytes() == annotation.read_bytes()
