from __future__ import annotations

from pathlib import Path

from ploidypatch.audit import audit_bundle, audit_fasta
from ploidypatch.gff import audit_gff, parse_attributes
from ploidypatch.io import (
    fasta_relation_id,
    normalize_feature_id,
    parse_fasta_header_fields,
    read_fai,
)


def write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


def test_parse_attributes_and_normalize_ids() -> None:
    attrs, malformed = parse_attributes("ID=transcript%3AAT1;Parent=gene%3AAG1")
    assert malformed == 0
    assert attrs == {"ID": "transcript:AT1", "Parent": "gene:AG1"}
    assert normalize_feature_id("transcript:AT1") == "AT1"
    assert normalize_feature_id("custom:AT1") == "custom:AT1"


def test_fasta_header_relation_prefers_explicit_transcript() -> None:
    header = (
        "Bra021453.1-P pep chromosome:Brapa_1.0:A01:1:99:1 "
        "gene:Bra021453 transcript:Bra021453.1 transcript_biotype:protein_coding"
    )
    fields = parse_fasta_header_fields(header)
    assert fields["gene"] == "Bra021453"
    assert fields["transcript"] == "Bra021453.1"
    assert fasta_relation_id("Bra021453.1-P", header) == (
        "Bra021453.1",
        "header:transcript",
    )


def test_fasta_header_relation_has_auditable_fallback() -> None:
    assert fasta_relation_id("transcript:T1", "transcript:T1 cds") == (
        "T1",
        "primary_fallback",
    )


def test_audit_bundle_passes_consistent_fixture(tmp_path: Path) -> None:
    fai = write(tmp_path / "genome.fa.fai", "chr1\t100\t0\t50\t51\n")
    gff = write(
        tmp_path / "genes.gff3",
        "##gff-version 3\n"
        "chr1\ttest\tgene\t1\t12\t.\t+\t.\tID=gene:G1\n"
        "chr1\ttest\tmRNA\t1\t12\t.\t+\t.\tID=transcript:T1;Parent=gene:G1\n"
        "chr1\ttest\texon\t1\t12\t.\t+\t.\tID=exon:E1;Parent=transcript:T1\n"
        "chr1\ttest\tCDS\t1\t12\t.\t+\t0\tID=CDS:C1;Parent=transcript:T1\n",
    )
    protein = write(
        tmp_path / "protein.fa",
        ">P1 pep gene:G1 transcript:T1 transcript_biotype:protein_coding\nMKT*\n",
    )
    cds = write(tmp_path / "cds.fa", ">transcript:T1\nATGAAAACCTAA\n")

    report = audit_bundle(
        gff_path=gff, protein_path=protein, cds_path=cds, fai_path=fai
    )

    assert report["quality_gate"]["grade"] == "pass"
    assert report["gff3"]["gene_ids"] == 1
    assert report["cross_checks"]["transcripts_missing_protein"] == 0
    assert report["protein"]["relation_id_source_counts"] == {
        "header:transcript": 1
    }


def test_gff_audit_detects_hierarchy_and_coordinate_errors(tmp_path: Path) -> None:
    gff = write(
        tmp_path / "bad.gff3",
        "chrX\ttest\tgene\t1\t500\t.\t+\t.\tID=G1\n"
        "chrX\ttest\tmRNA\t5\t4\t.\t+\t.\tID=T1;Parent=missing\n",
    )
    report, identifiers = audit_gff(gff, {"chr1": 100})
    assert report["invalid_coordinates"] == 1
    assert report["unknown_seqid_features"] == 2
    assert report["orphan_parent_ids"] == 1
    assert report["transcripts_without_gene_parent"] == 1
    assert identifiers.transcripts == frozenset({"T1"})
    assert report["examples"]["invalid_coordinates"] == [
        {
            "line_number": 2,
            "seqid": "chrX",
            "start": 5,
            "end": 4,
            "feature_type": "mRNA",
        }
    ]


def test_fasta_and_fai_validation(tmp_path: Path) -> None:
    protein = write(tmp_path / "p.fa", ">P1\nMA*Q\n>P1 duplicate\nMX-\n")
    report, index = audit_fasta(protein, "protein")
    assert index.relation_ids == frozenset({"P1"})
    assert report["duplicate_ids"] == 1
    assert report["internal_stop_records"] == 1
    assert report["invalid_character_records"] == 1
    assert report["invalid_character_counts"] == {"-": 1}
    assert report["examples"]["invalid_character_ids"] == ["P1"]
    assert report["relation_id_source_counts"] == {"primary_fallback": 2}

    fai = write(tmp_path / "x.fai", "chr1\t10\t0\t10\t11\n")
    assert read_fai(fai) == {"chr1": 10}


def test_partial_cds_is_observed_but_translation_consistency_passes(
    tmp_path: Path,
) -> None:
    gff = write(
        tmp_path / "partial.gff3",
        "chr1\ttest\tgene\t1\t5\t.\t+\t.\tID=gene:G1\n"
        "chr1\ttest\tmRNA\t1\t5\t.\t+\t.\tID=transcript:T1;Parent=gene:G1;biotype=protein_coding\n"
        "chr1\ttest\texon\t1\t5\t.\t+\t.\tParent=transcript:T1\n"
        "chr1\ttest\tCDS\t1\t5\t.\t+\t0\tParent=transcript:T1\n",
    )
    protein = write(tmp_path / "partial.pep.fa", ">P1 transcript:T1\nM\n")
    cds = write(tmp_path / "partial.cds.fa", ">T1\nATGAA\n")

    report = audit_bundle(gff_path=gff, protein_path=protein, cds_path=cds)

    assert report["quality_gate"]["grade"] == "pass"
    assert report["quality_gate"]["observations"]["cds_not_multiple_of_three"] == 1
    assert report["cross_checks"]["cds_frame_remainder_translation_consistent"] == 1
