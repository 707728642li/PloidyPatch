from __future__ import annotations

import csv
from pathlib import Path

from ploidypatch.cli import main
from ploidypatch.projection_select import select_projection_support_models


def write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8", newline="")
    return path


def candidate_text() -> str:
    return (
        "##gff-version 3\n"
        "chr1\ttest\tgene\t1\t20\t.\t+\t.\tID=G0\n"
        "###\n"
        "chr1\tPloidyPatchBaseline\tgene\t40\t60\t.\t+\t.\tID=PG1;miniprot_model=M1\n"
        "chr1\tPloidyPatchBaseline\tmRNA\t40\t60\t.\t+\t.\tID=PT1;Parent=PG1;miniprot_model=M1\n"
        "chr1\tPloidyPatchBaseline\texon\t40\t60\t.\t+\t.\tID=PE1;Parent=PT1\n"
        "chr1\tPloidyPatchBaseline\tCDS\t40\t60\t.\t+\t0\tParent=PT1\n"
        "chr1\tPloidyPatchBaseline\tgene\t80\t100\t.\t+\t.\tID=PG2;miniprot_model=M2\n"
        "chr1\tPloidyPatchBaseline\tmRNA\t80\t100\t.\t+\t.\tID=PT2;Parent=PG2;miniprot_model=M2\n"
        "chr1\tPloidyPatchBaseline\texon\t80\t100\t.\t+\t.\tID=PE2;Parent=PT2\n"
        "chr1\tPloidyPatchBaseline\tCDS\t80\t100\t.\t+\t0\tParent=PT2\n"
    )


def support_text() -> str:
    return (
        "model_id\tsupport_source_count\tsupport_sources\n"
        "M1\t2\ta,b\n"
        "M2\t1\ta\n"
    )


def test_select_projection_models_preserves_annotation_and_hierarchies(
    tmp_path: Path,
) -> None:
    candidate = write(tmp_path / "candidate.gff3", candidate_text())
    support = write(tmp_path / "support.tsv", support_text())
    output = tmp_path / "selected.gff3"
    selection = tmp_path / "selection.tsv"
    manifest = select_projection_support_models(
        candidate_gff_path=candidate,
        projection_support_tsv_path=support,
        output_gff_path=output,
        selection_tsv_path=selection,
    )

    text = output.read_text(encoding="utf-8")
    assert "ID=G0" in text
    assert "miniprot_model=M1" in text
    assert "Parent=PT1" in text
    assert "miniprot_model=M2" not in text
    assert "Parent=PT2" not in text
    assert manifest["counts"]["selected_models"] == 1
    with selection.open("r", encoding="utf-8", newline="") as handle:
        rows = {row["model_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
    assert rows["M1"]["support_group_count"] == "2"
    assert rows["M2"]["reason"] == "support_below_threshold"


def test_projection_selection_collapses_correlated_sources_through_cli(
    tmp_path: Path,
) -> None:
    candidate = write(tmp_path / "candidate.gff3", candidate_text())
    support = write(tmp_path / "support.tsv", support_text())
    groups = write(
        tmp_path / "groups.tsv",
        "source\tsupport_group\n"
        "a\tsame\n"
        "b\tsame\n",
    )
    output = tmp_path / "selected.gff3"
    assert (
        main(
            [
                "evidence",
                "select-projection-support",
                "--candidate-gff",
                str(candidate),
                "--projection-support",
                str(support),
                "--source-group-map",
                str(groups),
                "--min-support-group-count",
                "2",
                "--output-gff",
                str(output),
                "--selection-tsv",
                str(tmp_path / "selection.tsv"),
            ]
        )
        == 0
    )
    text = output.read_text(encoding="utf-8")
    assert "ID=G0" in text
    assert "PloidyPatchBaseline" not in text
