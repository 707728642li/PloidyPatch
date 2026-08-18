from __future__ import annotations

import csv
from pathlib import Path

import pytest

from ploidypatch.baseline import (
    MergedIntervalIndex,
    _overlap_with_index,
    adapt_annotation_gff_baseline,
    adapt_miniprot_baseline,
    prepare_reference_proteins,
    summarize_projection_support,
)
from ploidypatch.perturb import read_gff_document
from ploidypatch.score import build_annotation_index


def _write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8", newline="")
    return path


def test_prepare_reference_proteins_prefixes_and_maps_sources(tmp_path: Path) -> None:
    bra = _write(tmp_path / "bra.fa", ">Bra1 old description\nMAAA\n")
    bol = _write(tmp_path / "bol.fa", ">Bo1\nMCCC\n")
    output = tmp_path / "combined.fa"
    mapping = tmp_path / "combined.tsv"

    manifest = prepare_reference_proteins(
        protein_inputs=[f"bra_a={bra}", f"bol_c={bol}"],
        output_fasta_path=output,
        output_map_path=mapping,
    )

    assert output.read_text(encoding="utf-8") == (
        ">bra_a__Bra1\nMAAA\n>bol_c__Bo1\nMCCC\n"
    )
    with mapping.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert rows[0] == {
        "query_id": "bra_a__Bra1",
        "source": "bra_a",
        "source_record_id": "Bra1",
        "length_aa": "4",
        "source_header": "Bra1 old description",
    }
    assert manifest["outputs"]["fasta"]["records"] == 2
    assert manifest["inputs"]["bra_a"]["records"] == 1


def test_prepare_reference_proteins_fails_on_invalid_or_existing_output(
    tmp_path: Path,
) -> None:
    invalid = _write(tmp_path / "invalid.fa", ">p1\nMA.\n")
    with pytest.raises(ValueError, match="Invalid protein"):
        prepare_reference_proteins(
            protein_inputs=[f"source={invalid}"],
            output_fasta_path=tmp_path / "combined.fa",
            output_map_path=tmp_path / "combined.tsv",
        )

    valid = _write(tmp_path / "valid.fa", ">p1\nMAA\n")
    output = _write(tmp_path / "existing.fa", "occupied\n")
    with pytest.raises(FileExistsError, match="Refusing to overwrite"):
        prepare_reference_proteins(
            protein_inputs=[f"source={valid}"],
            output_fasta_path=output,
            output_map_path=tmp_path / "map.tsv",
        )


def test_adapt_miniprot_filters_existing_and_redundant_models(tmp_path: Path) -> None:
    blind = _write(
        tmp_path / "blind.gff3",
        "##gff-version 3\n"
        "chr1\ttest\tgene\t1\t30\t.\t+\t.\tID=gene:G1\n"
        "chr1\ttest\tmRNA\t1\t30\t.\t+\t.\tID=transcript:T1;Parent=gene:G1\n"
        "chr1\ttest\texon\t1\t30\t.\t+\t.\tParent=transcript:T1\n"
        "chr1\ttest\tCDS\t1\t30\t.\t+\t0\tParent=transcript:T1\n",
    )
    protein_map = _write(
        tmp_path / "proteins.tsv",
        "query_id\tsource\tsource_record_id\tlength_aa\tsource_header\n"
        "bra__p1\tbra\tp1\t10\tp1\n"
        "bra__p2\tbra\tp2\t10\tp2\n"
        "bol__p3\tbol\tp3\t10\tp3\n"
        "bol__p4\tbol\tp4\t10\tp4\n",
    )
    miniprot = _write(
        tmp_path / "miniprot.gff3",
        "##gff-version 3\n"
        "chr1\tminiprot\tmRNA\t1\t30\t100\t+\t.\tID=MP1;Rank=1;Identity=0.9;Positive=0.9;Target=bra__p1 1 10\n"
        "chr1\tminiprot\tCDS\t1\t30\t100\t+\t0\tParent=MP1\n"
        "chr1\tminiprot\tmRNA\t25\t54\t90\t+\t.\tID=MP2;Rank=1;Identity=0.9;Positive=0.9;Target=bra__p2 1 10\n"
        "chr1\tminiprot\tCDS\t25\t54\t90\t+\t0\tParent=MP2\n"
        "chr1\tminiprot\tmRNA\t26\t53\t80\t+\t.\tID=MP3;Rank=1;Identity=0.9;Positive=0.9;Target=bol__p3 1 10\n"
        "chr1\tminiprot\tCDS\t26\t53\t80\t+\t0\tParent=MP3\n"
        "chr1\tminiprot\tmRNA\t201\t230\t70\t+\t.\tID=MP4;Rank=1;Identity=0.3;Positive=0.4;Target=bol__p4 1 10\n"
        "chr1\tminiprot\tCDS\t201\t230\t70\t+\t0\tParent=MP4\n",
    )
    output = tmp_path / "candidate.gff3"
    decisions = tmp_path / "decisions.tsv"

    manifest = adapt_miniprot_baseline(
        perturbed_gff_path=blind,
        miniprot_gff_path=miniprot,
        protein_map_path=protein_map,
        output_gff_path=output,
        decisions_tsv_path=decisions,
    )

    assert manifest["models"]["input"] == 4
    assert manifest["models"]["accepted"] == 1
    assert manifest["models"]["decision_counts"] == {
        "accepted": 1,
        "identity_below_threshold": 1,
        "overlaps_blind_cds": 1,
        "redundant_projection": 1,
    }
    index = build_annotation_index(read_gff_document(output))
    assert len(index.transcripts) == 2
    added = index.transcripts["PPBASE_tx_000001"].signature
    assert added.exons == ((25, 54),)
    assert added.cds == ((25, 54, "0"),)
    with decisions.open("r", encoding="utf-8", newline="") as handle:
        decision_rows = {
            row["model_id"]: row for row in csv.DictReader(handle, delimiter="\t")
        }
    assert decision_rows["MP2"]["existing_cds_overlap_fraction"] == "0.200000"
    assert decision_rows["MP3"]["existing_cds_overlap_fraction"] == "0.178571"
    assert decision_rows["MP3"]["redundant_with_model_id"] == "MP2"
    assert decision_rows["MP2"]["redundant_with_model_id"] == ""

    support = tmp_path / "projection_support.tsv"
    support_manifest = summarize_projection_support(
        decisions_tsv_path=decisions,
        output_tsv_path=support,
    )
    with support.open("r", encoding="utf-8", newline="") as handle:
        support_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(support_rows) == 1
    assert support_rows[0]["model_id"] == "MP2"
    assert support_rows[0]["support_model_count"] == "2"
    assert support_rows[0]["support_query_count"] == "2"
    assert support_rows[0]["support_source_count"] == "2"
    assert support_rows[0]["support_sources"] == "bol,bra"
    assert support_manifest["counts"]["redundant_support_models"] == 1


def test_adapt_general_gff_preserves_blind_and_filters_candidates(
    tmp_path: Path,
) -> None:
    blind = _write(
        tmp_path / "blind.gff3",
        "##gff-version 3\n"
        "chr1\ttest\tgene\t1\t30\t.\t+\t.\tID=gene:G1\n"
        "chr1\ttest\tmRNA\t1\t30\t.\t+\t.\tID=transcript:T1;Parent=gene:G1\n"
        "chr1\ttest\texon\t1\t30\t.\t+\t.\tParent=transcript:T1\n"
        "chr1\ttest\tCDS\t1\t30\t.\t+\t0\tParent=transcript:T1\n",
    )
    candidates = _write(
        tmp_path / "upstream.gff3",
        "##gff-version 3\n"
        "chr1\ttool\tmRNA\t1\t30\t100\t+\t.\tID=C1\n"
        "chr1\ttool\tCDS\t1\t30\t.\t+\t0\tParent=C1\n"
        "chr1\ttool\tmRNA\t20\t60\t90\t+\t.\tID=C2;Evidence=strong\n"
        "chr1\ttool\texon\t23\t56\t.\t+\t.\tParent=C2\n"
        "chr1\ttool\tCDS\t25\t54\t.\t+\t0\tParent=C2\n"
        "chr1\ttool\tmRNA\t26\t53\t80\t+\t.\tID=C3\n"
        "chr1\ttool\tCDS\t26\t53\t.\t+\t0\tParent=C3\n"
        "chr1\ttool\tmRNA\t101\t120\t70\t+\t.\tID=N1\n"
        "chr2\ttool\tmRNA\t201\t230\t60\t+\t.\tID=U1\n"
        "chr2\ttool\tCDS\t201\t230\t.\t+\t0\tParent=U1\n"
        "chr1\ttool\tmRNA\t301\t330\t50\t+\t.\tID=P1\n"
        "chr1\ttool\tCDS\t301\t330\t.\t+\t.\tParent=P1\n",
    )
    output = tmp_path / "candidate.gff3"
    decisions = tmp_path / "decisions.tsv"

    manifest = adapt_annotation_gff_baseline(
        perturbed_gff_path=blind,
        candidate_gff_path=candidates,
        source="SynGAP-1.2.5",
        output_gff_path=output,
        decisions_tsv_path=decisions,
    )

    assert manifest["models"] == {
        "input": 6,
        "accepted": 1,
        "rejected": 5,
        "decision_counts": {
            "accepted": 1,
            "invalid_cds_phase": 1,
            "no_cds": 1,
            "overlaps_blind_cds": 1,
            "redundant_candidate": 1,
            "unknown_target_seqid": 1,
        },
    }
    index = build_annotation_index(read_gff_document(output))
    assert len(index.transcripts) == 2
    added = index.transcripts["PPGFF_tx_000001"].signature
    assert added.exons == ((23, 56),)
    assert added.cds == ((25, 54, "0"),)
    with decisions.open("r", encoding="utf-8", newline="") as handle:
        rows = {row["model_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
    assert rows["C2"]["status"] == "accepted"
    assert rows["C2"]["existing_cds_overlap_fraction"] == "0.200000"
    assert rows["C2"]["upstream_attributes"] == "ID=C2;Evidence=strong"
    assert rows["C3"]["redundant_with_model_id"] == "C2"
    assert output.read_text(encoding="utf-8").startswith(
        blind.read_text(encoding="utf-8") + "###\n"
    )

    with pytest.raises(FileExistsError, match="Refusing to overwrite"):
        adapt_annotation_gff_baseline(
            perturbed_gff_path=blind,
            candidate_gff_path=candidates,
            source="SynGAP-1.2.5",
            output_gff_path=output,
            decisions_tsv_path=decisions,
        )


def test_general_gff_adapter_can_explicitly_infer_all_missing_phases(
    tmp_path: Path,
) -> None:
    blind = _write(
        tmp_path / "blind.gff3",
        "##gff-version 3\n"
        "chr1\ttest\tgene\t1\t30\t.\t+\t.\tID=gene:G1\n"
        "chr1\ttest\tmRNA\t1\t30\t.\t+\t.\tID=transcript:T1;Parent=gene:G1\n"
        "chr1\ttest\tCDS\t1\t30\t.\t+\t0\tParent=transcript:T1\n",
    )
    candidates = _write(
        tmp_path / "upstream.gff3",
        "##gff-version 3\n"
        "chr1\ttool\tmRNA\t100\t200\t50\t-\t.\tID=MISSING_PHASE\n"
        "chr1\ttool\tCDS\t190\t200\t.\t-\t.\tParent=MISSING_PHASE\n"
        "chr1\ttool\tCDS\t100\t110\t.\t-\t.\tParent=MISSING_PHASE\n",
    )
    output = tmp_path / "candidate.gff3"
    decisions = tmp_path / "decisions.tsv"
    manifest = adapt_annotation_gff_baseline(
        perturbed_gff_path=blind,
        candidate_gff_path=candidates,
        source="phase-normalized",
        output_gff_path=output,
        decisions_tsv_path=decisions,
        infer_missing_cds_phase=True,
    )

    index = build_annotation_index(read_gff_document(output))
    added = index.transcripts["PPGFF_tx_000001"].signature
    assert added.cds == ((100, 110, "1"), (190, 200, "0"))
    assert manifest["parameters"]["infer_missing_cds_phase"] is True
    assert manifest["normalization"]["phase_normalization_counts"] == {
        "inferred_full_cds_first_phase_zero": 1
    }
    with decisions.open("r", encoding="utf-8", newline="") as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    assert row["status"] == "accepted"
    assert row["phase_normalization"] == (
        "inferred_full_cds_first_phase_zero"
    )


def test_indexed_interval_overlap_skips_nonoverlapping_prefix() -> None:
    reference = MergedIntervalIndex(
        intervals=((1, 10), (101, 120), (201, 220)),
        ends=(10, 120, 220),
    )
    assert _overlap_with_index([(105, 110), (115, 205)], reference) == 6 + 6 + 5
    assert _overlap_with_index([(50, 60)], reference) == 0
    assert _overlap_with_index([(1, 220)], reference) == 50
    assert _overlap_with_index([(1, 10)], None) == 0
