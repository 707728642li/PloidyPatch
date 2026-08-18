from __future__ import annotations

import json
from pathlib import Path

import pytest

from ploidypatch.synteny_io import (
    discover_primary_chromosome_seqids,
    prepare_wgdi_inputs,
)


GFF = (
    "##gff-version 3\n"
    "chr1\ttest\tregion\t1\t1000\t.\t+\t.\t"
    "ID=chr1;Name=A1;chromosome=A1;genome=chromosome\n"
    "chr1\ttest\tgene\t10\t90\t.\t+\t.\tID=gene:G1\n"
    "chr1\ttest\tmRNA\t10\t90\t.\t+\t.\tID=transcript:T1a;Parent=gene:G1\n"
    "chr1\ttest\tmRNA\t10\t90\t.\t+\t.\tID=transcript:T1b;Parent=gene:G1\n"
    "chr1\ttest\tgene\t110\t190\t.\t-\t.\tID=gene:G2\n"
    "chr1\ttest\tmRNA\t110\t190\t.\t-\t.\tID=transcript:T2;Parent=gene:G2\n"
    "scaffold1\ttest\tgene\t1\t50\t.\t+\t.\tID=gene:G3\n"
    "scaffold1\ttest\tmRNA\t1\t50\t.\t+\t.\tID=transcript:T3;Parent=gene:G3\n"
    "scaffold1\ttest\tregion\t1\t100\t.\t+\t.\t"
    "ID=scaffold1;Name=unplaced;genome=genomic\n"
)

PROTEINS = (
    ">P1short gene:G1 transcript:T1a\nMAA*\n"
    ">P1long gene:G1 transcript:T1b\nMAAAAA*\n"
    ">P2 transcript:T2\nMCCCC\n"
    ">P3 gene:G3 transcript:T3\nMDDD\n"
    ">P1invalid gene:G1 transcript:T1a\nMA.AAAAA\n"
    ">unmapped\nMXXX\n"
)


def write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8", newline="")
    return path


def test_prepare_wgdi_inputs_selects_longest_and_filters_small_seqids(
    tmp_path: Path,
) -> None:
    gff = write(tmp_path / "source.gff3", GFF)
    proteins = write(tmp_path / "proteins.fa", PROTEINS)
    fai = write(
        tmp_path / "genome.fai",
        "chr1\t1000\t0\t80\t81\nscaffold1\t100\t0\t80\t81\n",
    )
    output = tmp_path / "wgdi"

    manifest = prepare_wgdi_inputs(
        gff, proteins, fai, output, prefix="fixture", min_genes_per_seqid=2
    )

    assert (output / "fixture.wgdi.gff").read_text(encoding="utf-8") == (
        "chr1\tG1\t10\t90\t+\t1\n"
        "chr1\tG2\t110\t190\t-\t2\n"
    )
    assert (output / "fixture.wgdi.lens").read_text(encoding="utf-8") == (
        "chr1\t1000\t2\n"
    )
    peptide = (output / "fixture.wgdi.pep.fa").read_text(encoding="utf-8")
    assert peptide == ">G1\nMAAAAA\n>G2\nMCCCC\n"
    representatives = (output / "fixture.representatives.tsv").read_text(
        encoding="utf-8"
    )
    assert "G1\tP1long\theader:gene\t6\tchr1\t1" in representatives
    assert "G2\tP2\theader:transcript->gff:Parent\t5\tchr1\t2" in representatives
    assert manifest["counts"]["retained_genes"] == 2
    assert manifest["counts"]["genes_excluded_by_seqid_threshold"] == 1
    assert manifest["counts"]["unmapped_protein_records"] == 1
    assert manifest["counts"]["invalid_protein_records"] == 1
    assert manifest["counts"]["invalid_protein_character:."] == 1
    stored = json.loads(
        (output / "fixture.wgdi_inputs.manifest.json").read_text(encoding="utf-8")
    )
    assert stored == manifest


def test_prepare_wgdi_inputs_refuses_missing_fai_seqid_and_overwrite(
    tmp_path: Path,
) -> None:
    gff = write(tmp_path / "source.gff3", GFF)
    proteins = write(tmp_path / "proteins.fa", PROTEINS)
    incomplete_fai = write(tmp_path / "incomplete.fai", "chr1\t1000\t0\t80\t81\n")
    with pytest.raises(ValueError, match="absent from FAI"):
        prepare_wgdi_inputs(
            gff,
            proteins,
            incomplete_fai,
            tmp_path / "missing",
            prefix="fixture",
        )

    complete_fai = write(
        tmp_path / "complete.fai",
        "chr1\t1000\t0\t80\t81\nscaffold1\t100\t0\t80\t81\n",
    )
    output = tmp_path / "complete"
    prepare_wgdi_inputs(gff, proteins, complete_fai, output, prefix="fixture")
    with pytest.raises(FileExistsError, match="Refusing to overwrite"):
        prepare_wgdi_inputs(gff, proteins, complete_fai, output, prefix="fixture")


def test_prepare_wgdi_inputs_rejects_unsafe_prefix(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="prefix"):
        prepare_wgdi_inputs(
            tmp_path / "unused.gff3",
            tmp_path / "unused.fa",
            tmp_path / "unused.fai",
            tmp_path / "out",
            prefix="../escape",
        )


def test_prepare_wgdi_inputs_can_select_ncbi_primary_chromosomes(
    tmp_path: Path,
) -> None:
    gff = write(
        tmp_path / "source.gff3",
        GFF.replace(
            "scaffold1\ttest\tgene\t1\t50",
            "scaffold1\ttest\tgene\t50\t1",
        ),
    )
    proteins = write(tmp_path / "proteins.fa", PROTEINS)
    fai = write(
        tmp_path / "genome.fai",
        "chr1\t1000\t0\t80\t81\nscaffold1\t100\t0\t80\t81\n",
    )
    output = tmp_path / "primary"
    manifest = prepare_wgdi_inputs(
        gff,
        proteins,
        fai,
        output,
        prefix="primary",
        primary_chromosomes_only=True,
    )

    assert manifest["chromosome_labels"] == {"chr1": "A1"}
    assert manifest["counts"]["retained_seqids"] == 1
    assert "scaffold1" not in (output / "primary.wgdi.gff").read_text(
        encoding="utf-8"
    )


def test_prepare_wgdi_inputs_supports_root_transcript_gene_annotations(
    tmp_path: Path,
) -> None:
    gff = write(
        tmp_path / "root_transcripts.gff3",
        "chr1\ttest\tmRNA\t10\t90\t.\t+\t.\tID=R1\n"
        "chr1\ttest\tCDS\t10\t90\t.\t+\t0\tParent=R1\n",
    )
    proteins = write(tmp_path / "root_transcripts.fa", ">R1\nMAAAAA\n")
    fai = write(tmp_path / "root_transcripts.fai", "chr1\t1000\t0\t80\t81\n")
    output = tmp_path / "root_transcripts"

    manifest = prepare_wgdi_inputs(
        gff, proteins, fai, output, prefix="root", min_genes_per_seqid=1
    )

    assert manifest["parameters"]["gene_record_policy"] == "root_transcript_as_gene"
    assert manifest["counts"]["retained_genes"] == 1
    assert (output / "root.wgdi.gff").read_text(encoding="utf-8") == (
        "chr1\tR1\t10\t90\t+\t1\n"
    )
    assert (output / "root.wgdi.pep.fa").read_text(encoding="utf-8") == (
        ">R1\nMAAAAA\n"
    )


def test_discovers_ensembl_chromosome_features(tmp_path: Path) -> None:
    gff = write(
        tmp_path / "ensembl.gff3",
        "##gff-version 3\n"
        "1\tassembly\tchromosome\t1\t1000\t.\t.\t.\t"
        "ID=chromosome:1;Name=1\n"
        "scaffold1\tassembly\tsupercontig\t1\t100\t.\t.\t.\t"
        "ID=supercontig:scaffold1;Name=scaffold1\n",
    )
    seqids, labels = discover_primary_chromosome_seqids(gff)
    assert seqids == frozenset({"1"})
    assert labels == {"1": "1"}


def test_discovers_explicit_chromosome_aliases_without_guessing_scaffolds(
    tmp_path: Path,
) -> None:
    gff = write(
        tmp_path / "ensembl_refseq_import.gff3",
        "##gff-version 3\n"
        "CM014357.1\tassembly\tregion\t1\t1000\t.\t.\t.\t"
        "ID=region:CM014357.1;Alias=NC_041002,1,Chr01,NC_041002.1\n"
        "QZWG01000021.1\tassembly\tregion\t1\t100\t.\t.\t.\t"
        "ID=region:QZWG01000021.1;Alias=NW_021143573,Super-Scaffold_78\n"
        "NC_039768.1\tassembly\tregion\t1\t80\t.\t.\t.\t"
        "ID=region:NC_039768.1;Alias=NC_039768,MT,NC_039768.1\n",
    )
    seqids, labels = discover_primary_chromosome_seqids(gff)
    assert seqids == frozenset({"CM014357.1"})
    assert labels == {"CM014357.1": "1"}
