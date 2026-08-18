from __future__ import annotations

from pathlib import Path

import pytest

from ploidypatch.audit import audit_bundle
from ploidypatch.normalize import (
    OUTPUT_NAMES,
    PRIMARY_ANNOTATION_OUTPUT_NAMES,
    PROVIDER_GFF3_OUTPUT_NAMES,
    normalize_provider_gff3,
    prepare_ncbi_primary_bundle,
    prepare_primary_annotation_bundle,
    read_primary_seqid_table,
)
from ploidypatch.perturb import read_gff_document


def _write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8", newline="")
    return path


def test_prepare_ncbi_primary_bundle_filters_scaffolds_and_bad_translation(
    tmp_path: Path,
) -> None:
    gff = _write(
        tmp_path / "source.gff3",
        "##gff-version 3\n"
        "##sequence-region chr1 1 500\n"
        "chr1\tRefSeq\tregion\t1\t500\t.\t+\t.\tID=chr1;Name=A1;chromosome=A1;genome=chromosome\n"
        "chr1\tRefSeq\tgene\t1\t6\t.\t+\t.\tID=gene:G1\n"
        "chr1\tRefSeq\tmRNA\t1\t6\t.\t+\t.\tID=transcript:T1;Parent=gene:G1\n"
        "chr1\tRefSeq\texon\t1\t6\t.\t+\t.\tParent=transcript:T1\n"
        "chr1\tRefSeq\tCDS\t1\t6\t.\t+\t0\tParent=transcript:T1\n"
        "chr1\tRefSeq\tgene\t101\t106\t.\t+\t.\tID=gene:G2\n"
        "chr1\tRefSeq\tmRNA\t101\t106\t.\t+\t.\tID=transcript:T2;Parent=gene:G2\n"
        "chr1\tRefSeq\texon\t101\t106\t.\t+\t.\tParent=transcript:T2\n"
        "chr1\tRefSeq\tCDS\t101\t106\t.\t+\t0\tParent=transcript:T2\n"
        "##sequence-region scaffold1 1 500\n"
        "scaffold1\tRefSeq\tregion\t1\t500\t.\t+\t.\tID=scaffold1;Name=scaffold1;genome=genomic\n"
        "scaffold1\tRefSeq\tgene\t1\t6\t.\t+\t.\tID=gene:G3\n"
        "scaffold1\tRefSeq\tmRNA\t1\t6\t.\t+\t.\tID=transcript:T3;Parent=gene:G3\n"
        "scaffold1\tRefSeq\texon\t1\t6\t.\t+\t.\tParent=transcript:T3\n"
        "scaffold1\tRefSeq\tCDS\t1\t6\t.\t+\t0\tParent=transcript:T3\n",
    )
    protein = _write(
        tmp_path / "source.protein.fa",
        ">transcript:T1\nMA\n>transcript:T2\nM.\n>transcript:T3\nMA\n",
    )
    cds = _write(
        tmp_path / "source.cds.fa",
        ">transcript:T1\nATGGCT\n>transcript:T2\nATGGCT\n>transcript:T3\nATGGCT\n",
    )
    genome = _write(
        tmp_path / "source.genome.fa",
        ">chr1\n" + "A" * 500 + "\n>scaffold1\n" + "C" * 500 + "\n",
    )
    output = tmp_path / "primary"

    manifest = prepare_ncbi_primary_bundle(
        gff_path=gff,
        protein_path=protein,
        cds_path=cds,
        genome_path=genome,
        output_dir=output,
    )

    assert manifest["selection"]["primary_chromosomes"] == 1
    assert manifest["selection"]["retained_translation_transcripts"] == 1
    assert manifest["selection"]["excluded_translation_transcripts"] == 1
    assert manifest["selection"]["excluded_genes_without_retained_transcript"] == 1
    assert manifest["selection"]["exclusion_reason_counts"] == {
        "invalid_protein_character": 1
    }
    filtered_gff = (output / OUTPUT_NAMES["gff3"]).read_text(encoding="utf-8")
    assert "transcript:T1" in filtered_gff
    assert "transcript:T2" not in filtered_gff
    assert "transcript:T3" not in filtered_gff
    assert "##sequence-region scaffold1" not in filtered_gff
    filtered_genome = (output / OUTPUT_NAMES["genome"]).read_text(encoding="utf-8")
    assert ">chr1\n" in filtered_genome
    assert ">scaffold1\n" not in filtered_genome

    report = audit_bundle(
        gff_path=output / OUTPUT_NAMES["gff3"],
        protein_path=output / OUTPUT_NAMES["protein"],
        cds_path=output / OUTPUT_NAMES["cds"],
        fai_path=output / OUTPUT_NAMES["fai"],
    )
    assert report["quality_gate"]["grade"] == "pass"

    with pytest.raises(FileExistsError, match="Refusing to overwrite"):
        prepare_ncbi_primary_bundle(
            gff_path=gff,
            protein_path=protein,
            cds_path=cds,
            genome_path=genome,
            output_dir=output,
        )


def test_prepare_primary_annotation_bundle_preserves_ensembl_models(
    tmp_path: Path,
) -> None:
    gff = _write(
        tmp_path / "source.gff3",
        "##gff-version 3\n"
        "##sequence-region chr1 1 200\n"
        "chr1\tEnsembl\tchromosome\t1\t200\t.\t.\t.\tID=chromosome:chr1;Name=1\n"
        "chr1\tEnsembl\tgene\t1\t9\t.\t+\t.\tID=gene:G1\n"
        "chr1\tEnsembl\tmRNA\t1\t9\t.\t+\t.\tID=transcript:T1;Parent=gene:G1\n"
        "chr1\tEnsembl\tCDS\t1\t9\t.\t+\t0\tParent=transcript:T1\n"
        "chr1\tEnsembl\tgene\t101\t109\t.\t+\t.\tID=gene:Gbad\n"
        "chr1\tEnsembl\tmRNA\t101\t109\t.\t+\t.\tID=transcript:Tbad;Parent=gene:Gbad\n"
        "chr1\tEnsembl\tCDS\t101\t109\t.\t+\t0\tParent=transcript:Tbad\n"
        "##sequence-region scaffold1 1 100\n"
        "scaffold1\tEnsembl\tgene\t1\t9\t.\t+\t.\tID=gene:G2\n"
        "scaffold1\tEnsembl\tmRNA\t1\t9\t.\t+\t.\tID=transcript:T2;Parent=gene:G2\n"
        "scaffold1\tEnsembl\tCDS\t1\t9\t.\t+\t0\tParent=transcript:T2\n",
    )
    genome = _write(
        tmp_path / "source.genome.fa",
        ">chr1 chromosome 1\n" + "A" * 200 + "\n"
        ">scaffold1\n" + "C" * 100 + "\n",
    )
    output = tmp_path / "primary_annotation"

    manifest = prepare_primary_annotation_bundle(
        gff_path=gff,
        genome_path=genome,
        output_dir=output,
    )

    assert manifest["selection"]["primary_chromosomes"] == 1
    assert manifest["selection"]["translation_filtering"] is False
    assert manifest["selection"]["chromosome_labels"] == {"chr1": "1"}
    assert manifest["selection"]["retained_feature_counts"] == {
        "CDS": 2,
        "chromosome": 1,
        "gene": 2,
        "mRNA": 2,
    }
    filtered_gff = (
        output / PRIMARY_ANNOTATION_OUTPUT_NAMES["gff3"]
    ).read_text(encoding="utf-8")
    assert "transcript:Tbad" in filtered_gff
    assert "transcript:T2" not in filtered_gff
    assert "##sequence-region scaffold1" not in filtered_gff
    filtered_genome = (
        output / PRIMARY_ANNOTATION_OUTPUT_NAMES["genome"]
    ).read_text(encoding="utf-8")
    assert ">chr1 chromosome 1\n" in filtered_genome
    assert ">scaffold1\n" not in filtered_genome

    with pytest.raises(FileExistsError, match="Refusing to overwrite"):
        prepare_primary_annotation_bundle(
            gff_path=gff,
            genome_path=genome,
            output_dir=output,
        )


def test_prepare_primary_annotation_bundle_accepts_explicit_seqid_table(
    tmp_path: Path,
) -> None:
    gff = _write(
        tmp_path / "community.gff3",
        "A01\tPhytozome\tgene\t1\t9\t.\t+\t.\tID=GA1\n"
        "A01\tPhytozome\tmRNA\t1\t9\t.\t+\t.\tID=TA1;Parent=GA1\n"
        "A01\tPhytozome\tCDS\t1\t9\t.\t+\t0\tParent=TA1\n"
        "D01\tPhytozome\tgene\t1\t9\t.\t+\t.\tID=GD1\n"
        "D01\tPhytozome\tmRNA\t1\t9\t.\t+\t.\tID=TD1;Parent=GD1\n"
        "D01\tPhytozome\tCDS\t1\t9\t.\t+\t0\tParent=TD1\n"
        "scaffold_1\tPhytozome\tgene\t1\t9\t.\t+\t.\tID=GS1\n",
    )
    genome = _write(
        tmp_path / "community.fa",
        ">A01\nAAAAAAAAA\n>D01\nCCCCCCCCC\n>scaffold_1\nGGGGGGGGG\n",
    )
    seqids = _write(
        tmp_path / "primary_seqids.tsv",
        "seqid\tchromosome_label\nA01\tA01\nD01\tD01\n",
    )
    output = tmp_path / "primary"

    parsed_seqids, labels = read_primary_seqid_table(seqids)
    assert parsed_seqids == frozenset({"A01", "D01"})
    assert labels == {"A01": "A01", "D01": "D01"}

    manifest = prepare_primary_annotation_bundle(
        gff_path=gff,
        genome_path=genome,
        output_dir=output,
        primary_seqid_table_path=seqids,
    )

    assert manifest["selection"]["method"] == "explicit primary seqid table"
    assert manifest["selection"]["primary_chromosomes"] == 2
    assert manifest["sources"]["primary_seqid_table"]["sha256"]
    filtered_gff = (
        output / PRIMARY_ANNOTATION_OUTPUT_NAMES["gff3"]
    ).read_text(encoding="utf-8")
    assert "GA1" in filtered_gff
    assert "GD1" in filtered_gff
    assert "GS1" not in filtered_gff
    filtered_genome = (
        output / PRIMARY_ANNOTATION_OUTPUT_NAMES["genome"]
    ).read_text(encoding="utf-8")
    assert ">A01\n" in filtered_genome
    assert ">D01\n" in filtered_genome
    assert ">scaffold_1\n" not in filtered_genome


def test_prepare_primary_annotation_bundle_can_canonicalize_fasta_headers(
    tmp_path: Path,
) -> None:
    gff = _write(
        tmp_path / "source.gff3",
        "chr1\tsource\tgene\t1\t9\t.\t+\t.\tID=G1\n"
        "chr1\tsource\tmRNA\t1\t9\t.\t+\t.\tID=T1;Parent=G1\n"
        "chr1\tsource\tCDS\t1\t9\t.\t+\t0\tParent=T1\n",
    )
    genome = _write(tmp_path / "source.fa", ">chr1 chromosome description\nAAAAAAAAA\n")
    table = _write(tmp_path / "primary.tsv", "seqid\tchromosome_label\nchr1\t1\n")
    output = tmp_path / "canonical"

    manifest = prepare_primary_annotation_bundle(
        gff_path=gff,
        genome_path=genome,
        primary_seqid_table_path=table,
        output_dir=output,
        canonical_fasta_headers=True,
    )

    assert (output / PRIMARY_ANNOTATION_OUTPUT_NAMES["genome"]).read_text(
        encoding="utf-8"
    ) == ">chr1\nAAAAAAAAA\n"
    assert manifest["selection"]["genome_fasta_headers"] == "exact_seqid"


def test_read_primary_seqid_table_rejects_duplicate_seqids(tmp_path: Path) -> None:
    table = _write(
        tmp_path / "duplicate.tsv",
        "seqid\tchromosome_label\nA01\tA01\nA01\tA1\n",
    )
    with pytest.raises(ValueError, match="duplicate seqid"):
        read_primary_seqid_table(table)


def test_normalize_provider_gff3_applies_only_enabled_audited_repairs(
    tmp_path: Path,
) -> None:
    source = _write(
        tmp_path / "provider.gff3",
        "##gff-version 3\n"
        "chr1\tLab\tgene\t1\t20\t.\t+\t.\t"
        "ID=G1;Note=DNA binding;sequence-specific DNA binding\n"
        "chr1\tLab\tmRNA\t1\t20\t.\t+\t.\tID=T1;Parent=G1\n"
        "chr1\tLab\texon\t1\t10\t.\t+\t.\tID=E1;Parent=T1\n"
        "chr1\tLab\texon\t10\t20\t.\t+\t.\tID=E2;Parent=T1\n"
        "chr1\tLab\tintron\t11\t9\t.\t+\t.\tID=I1;Parent=T1\n"
        "chr1\tLab\tCDS\t1\t20\t.\t+\t0\tID=C1;Parent=T1\n"
        "##FASTA\n"
        ">chr1\n"
        "AAAAAAAAAAAAAAAAAAAA\n",
    )

    with pytest.raises(ValueError, match="Malformed GFF3 attributes"):
        normalize_provider_gff3(gff_path=source, output_dir=tmp_path / "strict")

    output = tmp_path / "compatible"
    manifest = normalize_provider_gff3(
        gff_path=source,
        output_dir=output,
        repair_unescaped_note_semicolons=True,
        drop_invalid_intron_intervals=True,
        strip_embedded_fasta=True,
    )

    assert manifest["policy"] == {
        "repair_unescaped_note_semicolons": True,
        "drop_invalid_intron_intervals": True,
        "strip_embedded_fasta": True,
    }
    assert manifest["observed"]["repaired_note_records"] == 1
    assert manifest["observed"]["repaired_note_semicolons"] == 1
    assert manifest["observed"]["dropped_invalid_intron_intervals"] == 1
    assert manifest["observed"]["stripped_embedded_fasta_directives"] == 1
    assert manifest["observed"]["stripped_embedded_fasta_lines"] == 2
    compatible = (output / PROVIDER_GFF3_OUTPUT_NAMES["gff3"]).read_text(
        encoding="utf-8"
    )
    assert "Note=DNA binding%3Bsequence-specific DNA binding" in compatible
    assert "\tintron\t" not in compatible
    assert "##FASTA" not in compatible
    assert len(read_gff_document(output / "sanitized.gff3").records) == 5

    with pytest.raises(FileExistsError, match="Refusing to overwrite"):
        normalize_provider_gff3(gff_path=source, output_dir=output)


def test_normalize_provider_gff3_never_drops_invalid_structural_records(
    tmp_path: Path,
) -> None:
    source = _write(
        tmp_path / "invalid_cds.gff3",
        "chr1\tLab\tCDS\t20\t10\t.\t+\t0\tID=C1;Parent=T1\n",
    )
    with pytest.raises(ValueError, match="Invalid GFF3 interval"):
        normalize_provider_gff3(
            gff_path=source,
            output_dir=tmp_path / "out",
            drop_invalid_intron_intervals=True,
        )
