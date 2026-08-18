from __future__ import annotations

import csv
import hashlib
from pathlib import Path

import pytest

from ploidypatch.seqid_alias import SCHEMA_VERSION, rewrite_gff3_seqids_exact


def write_tsv(
    path: Path, fields: tuple[str, ...], rows: list[tuple[str, ...]]
) -> Path:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(fields)
        writer.writerows(rows)
    return path


def test_exact_seqid_rewrite_preserves_coordinates_attributes_and_primary_only(
    tmp_path: Path,
) -> None:
    source = tmp_path / "source.gff3"
    source.write_text(
        "##gff-version 3\n"
        "##sequence-region Gardenia1 1 100\n"
        "##sequence-region contig0 1 50\n"
        "# provider comment\n"
        "Gardenia1\tx\tgene\t2\t90\t.\t+\t.\tID=G1;Name=G1\n"
        "Gardenia2\tx\tgene\t4\t80\t.\t-\t.\tID=G2;Name=G2\n"
        "contig0\tx\tgene\t1\t20\t.\t+\t.\tID=DROP\n",
        encoding="utf-8",
    )
    aliases = write_tsv(
        tmp_path / "aliases.tsv",
        ("gff_seqid", "genome_seqid", "chromosome_label"),
        [("Gardenia1", "CM1.1", "1"), ("Gardenia2", "CM2.1", "2")],
    )
    primary = write_tsv(
        tmp_path / "primary.tsv",
        ("seqid", "chromosome_label"),
        [("CM1.1", "1"), ("CM2.1", "2")],
    )
    output = tmp_path / "normalized.gff3"
    manifest_path = tmp_path / "manifest.json"
    manifest = rewrite_gff3_seqids_exact(
        input_gff_path=source,
        alias_table_path=aliases,
        primary_seqid_table_path=primary,
        output_gff_path=output,
        manifest_path=manifest_path,
    )

    assert manifest["schema_version"] == SCHEMA_VERSION
    assert manifest["counts"]["skipped_nonprimary_features"] == 1
    assert manifest["counts"]["skipped_nonprimary_sequence_regions"] == 1
    assert manifest["coordinate_or_attribute_change"] is False
    text = output.read_text(encoding="utf-8")
    assert "##sequence-region CM1.1 1 100\n" in text
    assert "CM1.1\tx\tgene\t2\t90\t.\t+\t.\tID=G1;Name=G1\n" in text
    assert "CM2.1\tx\tgene\t4\t80\t.\t-\t.\tID=G2;Name=G2\n" in text
    assert "contig0" not in text
    assert "DROP" not in text
    assert manifest["outputs"]["gff3_sha256"] == hashlib.sha256(
        output.read_bytes()
    ).hexdigest()

    with pytest.raises(FileExistsError, match="overwrite"):
        rewrite_gff3_seqids_exact(
            input_gff_path=source,
            alias_table_path=aliases,
            primary_seqid_table_path=primary,
            output_gff_path=output,
            manifest_path=manifest_path,
        )


def test_exact_seqid_rewrite_rejects_nonbijective_or_incomplete_aliases(
    tmp_path: Path,
) -> None:
    source = tmp_path / "source.gff3"
    source.write_text(
        "Gardenia1\tx\tgene\t1\t10\t.\t+\t.\tID=G1\n",
        encoding="utf-8",
    )
    aliases = write_tsv(
        tmp_path / "aliases.tsv",
        ("gff_seqid", "genome_seqid", "chromosome_label"),
        [("Gardenia1", "CM1.1", "1"), ("Gardenia2", "CM1.1", "2")],
    )
    primary = write_tsv(
        tmp_path / "primary.tsv",
        ("seqid", "chromosome_label"),
        [("CM1.1", "1"), ("CM2.1", "2")],
    )
    with pytest.raises(ValueError, match="Target genome seqid alias is not unique"):
        rewrite_gff3_seqids_exact(
            input_gff_path=source,
            alias_table_path=aliases,
            primary_seqid_table_path=primary,
            output_gff_path=tmp_path / "out.gff3",
            manifest_path=tmp_path / "manifest.json",
        )
    assert not (tmp_path / "out.gff3").exists()


def test_exact_seqid_rewrite_rejects_malformed_gff_without_partial_output(
    tmp_path: Path,
) -> None:
    source = tmp_path / "source.gff3"
    source.write_text("Gardenia1\tmalformed\n", encoding="utf-8")
    aliases = write_tsv(
        tmp_path / "aliases.tsv",
        ("gff_seqid", "genome_seqid", "chromosome_label"),
        [("Gardenia1", "CM1.1", "1")],
    )
    primary = write_tsv(
        tmp_path / "primary.tsv",
        ("seqid", "chromosome_label"),
        [("CM1.1", "1")],
    )
    output = tmp_path / "out.gff3"
    with pytest.raises(ValueError, match="expected 9 fields"):
        rewrite_gff3_seqids_exact(
            input_gff_path=source,
            alias_table_path=aliases,
            primary_seqid_table_path=primary,
            output_gff_path=output,
            manifest_path=tmp_path / "manifest.json",
        )
    assert not output.exists()
