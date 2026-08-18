from __future__ import annotations

import json
from pathlib import Path

import pytest

from ploidypatch.perturb import (
    BOUNDARY_SHIFT_EVENT,
    COPY_COLLAPSE_EVENT,
    FUSED_GENE_EVENT,
    MISSING_INTERNAL_EXON_EVENT,
    SPLIT_GENE_EVENT,
    restore_gff_from_truth,
)
from ploidypatch.score import build_annotation_index, score_annotation_repair
from ploidypatch.structure_perturb import generate_structure_benchmark
from ploidypatch.perturb import read_gff_document


def _gene(gene: str, transcript: str, start: int) -> str:
    intervals = [(start, start + 29), (start + 50, start + 79),
                 (start + 100, start + 129), (start + 150, start + 179)]
    lines = [
        f"chr1\ttest\tgene\t{start}\t{start + 179}\t.\t+\t.\tID={gene}\n",
        f"chr1\ttest\tmRNA\t{start}\t{start + 179}\t.\t+\t.\t"
        f"ID={transcript};Parent={gene}\n",
    ]
    for index, (left, right) in enumerate(intervals, start=1):
        lines.extend(
            [
                f"chr1\ttest\texon\t{left}\t{right}\t.\t+\t.\t"
                f"ID={gene}.exon{index};Parent={transcript}\n",
                f"chr1\ttest\tCDS\t{left}\t{right}\t.\t+\t0\t"
                f"ID={gene}.cds{index};Parent={transcript}\n",
            ]
        )
    return "".join(lines)


SOURCE_GFF = (
    "##gff-version 3\n"
    + _gene("gene:G1", "transcript:T1", 1)
    + _gene("gene:G2", "transcript:T2", 201)
    + _gene("gene:G3", "transcript:T3", 501)
)


def _write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8", newline="")
    return path


@pytest.mark.parametrize(
    "event_type",
    [
        MISSING_INTERNAL_EXON_EVENT,
        BOUNDARY_SHIFT_EVENT,
        SPLIT_GENE_EVENT,
        FUSED_GENE_EVENT,
    ],
)
def test_structure_perturbations_are_reversible_and_perfectly_scorable(
    tmp_path: Path, event_type: str
) -> None:
    source = _write(tmp_path / "source.gff3", SOURCE_GFF)
    run = tmp_path / event_type
    manifest = generate_structure_benchmark(
        source, run, event_type=event_type, count=1, seed=17
    )

    assert manifest["perturbation"]["generated_events"] == 1
    truth = json.loads((run / "hidden_truth.json").read_text(encoding="utf-8"))
    assert truth["events"][0]["event_type"] == event_type
    assert truth["events"][0]["line_edits"]
    assert (run / "perturbed.gff3").read_bytes() != source.read_bytes()

    restored = run / "restored.gff3"
    restore = restore_gff_from_truth(
        run / "perturbed.gff3", run / "hidden_truth.json", restored
    )
    assert restored.read_bytes() == source.read_bytes()
    assert restore["events"] == 1

    report = score_annotation_repair(
        source,
        run / "perturbed.gff3",
        source,
        run / "hidden_truth.json",
        include_event_details=True,
    )
    assert report["schema_version"] == "ploidypatch.annotation_repair_score.v5"
    assert report["strict_transcript_structure"]["f1"] == 1.0
    assert report["strict_cds_chain"]["f1"] == 1.0
    assert report["event_recovery"]["complete_transcript_recall"] == 1.0
    assert report["event_recovery"]["exact_gene_recall"] == 1.0
    assert report["event_recovery"]["complete_error_removal_recall"] == 1.0
    assert report["event_details"][0]["complete_error_removal"] is True


def test_no_repair_does_not_get_credit_for_retaining_introduced_error(
    tmp_path: Path,
) -> None:
    source = _write(tmp_path / "source.gff3", SOURCE_GFF)
    run = tmp_path / "boundary"
    generate_structure_benchmark(
        source, run, event_type=BOUNDARY_SHIFT_EVENT, count=1, seed=3
    )

    report = score_annotation_repair(
        source,
        run / "perturbed.gff3",
        run / "perturbed.gff3",
        run / "hidden_truth.json",
        include_event_details=True,
    )

    assert report["strict_transcript_structure"]["recall"] == 0.0
    assert report["event_recovery"]["complete_transcript_recall"] == 0.0
    assert report["event_recovery"]["complete_error_removal_recall"] == 0.0
    assert report["event_details"][0]["introduced_error_structures_remaining"] == 1


def test_split_and_fusion_change_gene_grouping_in_opposite_directions(
    tmp_path: Path,
) -> None:
    source = _write(tmp_path / "source.gff3", SOURCE_GFF)
    source_index = build_annotation_index(read_gff_document(source))

    split_run = tmp_path / "split"
    generate_structure_benchmark(
        source, split_run, event_type=SPLIT_GENE_EVENT, count=1, seed=1
    )
    split_index = build_annotation_index(
        read_gff_document(split_run / "perturbed.gff3")
    )
    assert len(split_index.gene_signatures) == len(source_index.gene_signatures) + 1

    fusion_run = tmp_path / "fusion"
    generate_structure_benchmark(
        source, fusion_run, event_type=FUSED_GENE_EVENT, count=1, seed=1
    )
    fusion_index = build_annotation_index(
        read_gff_document(fusion_run / "perturbed.gff3")
    )
    assert len(fusion_index.gene_signatures) == len(source_index.gene_signatures) - 1


def test_copy_collapse_requires_external_pairs_and_is_reversible(
    tmp_path: Path,
) -> None:
    source = _write(tmp_path / "source.gff3", SOURCE_GFF)
    with pytest.raises(ValueError, match="requires --pair-tsv"):
        generate_structure_benchmark(
            source,
            tmp_path / "missing_pairs",
            event_type=COPY_COLLAPSE_EVENT,
            count=1,
            seed=9,
        )

    pairs = _write(
        tmp_path / "pairs.tsv",
        "gene_id_a\tgene_id_b\n"
        "gene:G1\tgene:G3\n",
    )
    run = tmp_path / "copy"
    generate_structure_benchmark(
        source,
        run,
        event_type=COPY_COLLAPSE_EVENT,
        count=1,
        seed=9,
        pair_tsv_path=pairs,
    )
    truth = json.loads((run / "hidden_truth.json").read_text(encoding="utf-8"))
    details = truth["events"][0]["details"]
    assert len(details["pair_gene_ids"]) == 2
    assert details["retained_partner_gene_id"] in details["pair_gene_ids"]

    restored = run / "restored.gff3"
    restore_gff_from_truth(run / "perturbed.gff3", run / "hidden_truth.json", restored)
    assert restored.read_bytes() == source.read_bytes()
    report = score_annotation_repair(
        source, run / "perturbed.gff3", source, run / "hidden_truth.json"
    )
    assert report["event_recovery"]["complete_transcript_recall"] == 1.0
    assert report["event_recovery"]["events_with_introduced_errors"] == 0


def test_copy_collapse_supports_cds_only_transcript_annotations(
    tmp_path: Path,
) -> None:
    source = _write(
        tmp_path / "cds_only.gff3",
        SOURCE_GFF.replace("##gff-version 3\n", "").replace(
            "\ttest\texon\t", "\ttest\tignored_exon\t"
        ),
    )
    pairs = _write(
        tmp_path / "cds_only_pairs.tsv",
        "gene_id_a\tgene_id_b\ngene:G1\tgene:G3\n",
    )
    run = tmp_path / "cds_only_copy"

    manifest = generate_structure_benchmark(
        source,
        run,
        event_type=COPY_COLLAPSE_EVENT,
        count=1,
        seed=19,
        pair_tsv_path=pairs,
    )

    assert manifest["perturbation"]["generated_events"] == 1
    restored = run / "restored.gff3"
    restore_gff_from_truth(run / "perturbed.gff3", run / "hidden_truth.json", restored)
    assert restored.read_bytes() == source.read_bytes()


def test_restore_detects_modified_structural_replacement(tmp_path: Path) -> None:
    source = _write(tmp_path / "source.gff3", SOURCE_GFF)
    run = tmp_path / "shift"
    generate_structure_benchmark(
        source, run, event_type=BOUNDARY_SHIFT_EVENT, count=1, seed=2
    )
    text = (run / "perturbed.gff3").read_text(encoding="utf-8")
    (run / "perturbed.gff3").write_text(
        text.replace("\ttest\tgene\t", "\tchanged\tgene\t", 1),
        encoding="utf-8",
        newline="",
    )
    with pytest.raises(ValueError, match="replacement|checksum"):
        restore_gff_from_truth(
            run / "perturbed.gff3",
            run / "hidden_truth.json",
            run / "restored.gff3",
        )
    with pytest.raises(ValueError, match="Perturbed GFF3 text checksum"):
        score_annotation_repair(
            source,
            run / "perturbed.gff3",
            source,
            run / "hidden_truth.json",
        )
