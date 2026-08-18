from __future__ import annotations

import json
from pathlib import Path

import pytest

from ploidypatch.artifact_manifest import verify_sha256sums
from ploidypatch.io import iter_fasta
from ploidypatch.perturb import read_gff_document
from ploidypatch.protein_universe import (
    build_protein_supported_gene_universe,
    collapse_exact_structural_duplicates,
)


def _write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


def _gff(*, duplicate: bool = False, drift: bool = False) -> str:
    rows = ["##gff-version 3\n"]
    models = [
        ("g1", "t1", "p1", 1, 100, ((1, 20, "0"), (41, 80, "2"))),
        ("g2", "t2", "p2", 201, 300, ((201, 260, "0"),)),
        ("g3", "t3", "p3", 401, 500, ((401, 450, "0"),)),
    ]
    for gene, transcript, protein, start, end, parts in models:
        gene_rows = [
            f"chr1\tx\tgene\t{start}\t{end}\t.\t+\t.\tID=gene:{gene};coge_fid=2\n"
        ]
        if duplicate:
            gene_rows.append(
                f"chr1\tx\tgene\t{start}\t{end}\t.\t+\t.\tID=gene:{gene};coge_fid=1\n"
            )
        rows.extend(gene_rows)
        transcript_rows = [
            f"chr1\tx\tmRNA\t{start}\t{end}\t.\t+\t.\tID=transcript:{transcript};Parent=gene:{gene};coge_fid=2\n"
        ]
        if duplicate:
            transcript_rows.append(
                f"chr1\tx\tmRNA\t{start}\t{end}\t.\t+\t.\tID=transcript:{transcript};Parent=gene:{gene};coge_fid=1\n"
            )
        rows.extend(transcript_rows)
        for index, (cds_start, cds_end, phase) in enumerate(parts, start=1):
            cds = (
                f"chr1\tx\tCDS\t{cds_start}\t{cds_end}\t.\t+\t{phase}\t"
                f"ID=cds:{protein};Parent=transcript:{transcript};protein_id={protein};coge_fid=2\n"
            )
            rows.append(cds)
            if duplicate:
                rows.append(cds.replace("coge_fid=2", "coge_fid=1"))
    if drift:
        rows.append(
            "chr1\tx\tgene\t700\t800\t.\t+\t.\tID=gene:g1;coge_fid=3\n"
        )
    return "".join(rows)


def test_exact_duplicates_collapse_and_unsupported_genes_are_explicit(tmp_path: Path) -> None:
    gff = _write(tmp_path / "input.gff3", _gff(duplicate=True))
    proteins = _write(
        tmp_path / "proteins.fa",
        ">p1 transcript:t1\n" + "M" * 30 + "\n>p2 transcript:t2\n" + "M" * 10 + "\n",
    )
    output = tmp_path / "universe"
    manifest = build_protein_supported_gene_universe(
        primary_gff_path=gff,
        provider_protein_path=proteins,
        output_dir=output,
        species_id="fixture",
        holdout_id="holdout",
        policy_id="policy",
    )
    assert manifest["counts"]["normalized_genes"] == 3
    assert manifest["counts"]["coding_genes"] == 3
    assert manifest["counts"]["genes_with_exact_provider_protein"] == 2
    assert manifest["counts"]["genes_excluded_without_exact_provider_protein"] == 1
    assert manifest["counts"]["collapsed_duplicate_rows"] == 10
    assert manifest["counts"]["duplicate_groups_with_attribute_variants"] == 10
    retained = read_gff_document(output / "protein_supported.gff3")
    assert sum(record.feature_type == "gene" for record in retained.records) == 2
    # Multipart CDS is preserved, not mistaken for an ID collision.
    assert sum(record.feature_type == "CDS" for record in retained.records) == 3
    assert [item[0] for item in iter_fasta(output / "representative.protein.fa")] == [
        "p1",
        "p2",
    ]
    assert (output / "excluded_genes.tsv").read_text(encoding="utf-8").splitlines() == [
        "gene_id\treason",
        "g3\tno_exact_provider_protein",
    ]
    verify_sha256sums(output, ignore_checksum_file=True)
    with pytest.raises(FileExistsError):
        build_protein_supported_gene_universe(
            primary_gff_path=gff,
            provider_protein_path=proteins,
            output_dir=output,
            species_id="fixture",
            holdout_id="holdout",
            policy_id="policy",
        )


def test_gene_multilocus_collision_is_rejected(tmp_path: Path) -> None:
    document = read_gff_document(_write(tmp_path / "drift.gff3", _gff(drift=True)))
    with pytest.raises(ValueError, match="multiple structures"):
        collapse_exact_structural_duplicates(document.records)


def test_discordant_exact_protein_relations_are_rejected(tmp_path: Path) -> None:
    gff = _write(tmp_path / "input.gff3", _gff())
    proteins = _write(
        tmp_path / "proteins.fa",
        ">p1 gene:g2 transcript:t1\nMMMM\n",
    )
    with pytest.raises(ValueError, match="discordant exact gene mappings"):
        build_protein_supported_gene_universe(
            primary_gff_path=gff,
            provider_protein_path=proteins,
            output_dir=tmp_path / "output",
            species_id="fixture",
            holdout_id="holdout",
            policy_id="policy",
        )


def test_no_prefix_or_fuzzy_protein_rescue(tmp_path: Path) -> None:
    gff = _write(tmp_path / "input.gff3", _gff())
    proteins = _write(tmp_path / "proteins.fa", ">t1_suffix\nMMMM\n")
    with pytest.raises(ValueError, match="No coding gene has an exact provider protein"):
        build_protein_supported_gene_universe(
            primary_gff_path=gff,
            provider_protein_path=proteins,
            output_dir=tmp_path / "output",
            species_id="fixture",
            holdout_id="holdout",
            policy_id="policy",
        )


def test_manifest_is_truth_and_candidate_free(tmp_path: Path) -> None:
    gff = _write(tmp_path / "input.gff3", _gff())
    proteins = _write(tmp_path / "proteins.fa", ">p1 transcript:t1\nMMMM\n")
    output = tmp_path / "output"
    build_protein_supported_gene_universe(
        primary_gff_path=gff,
        provider_protein_path=proteins,
        output_dir=output,
        species_id="fixture",
        holdout_id="holdout",
        policy_id="policy",
    )
    manifest = json.loads((output / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["truth_access"] is False
    assert manifest["candidate_access"] is False
    assert manifest["fuzzy_mapping_used"] is False
    assert manifest["source_annotation_preserved"] is True


def test_parentless_cds_fragment_is_reported_not_guessed(tmp_path: Path) -> None:
    gff = _write(
        tmp_path / "input.gff3",
        _gff()
        + "chr1\tx\tCDS\t700\t750\t.\t+\t0\tID=orphan;Name=orphan\n",
    )
    proteins = _write(tmp_path / "proteins.fa", ">p1 transcript:t1\nMMMM\n")
    manifest = build_protein_supported_gene_universe(
        primary_gff_path=gff,
        provider_protein_path=proteins,
        output_dir=tmp_path / "output",
        species_id="fixture",
        holdout_id="holdout",
        policy_id="policy",
    )
    assert manifest["counts"]["orphan_CDS_rows_without_Parent_excluded"] == 1
    assert "orphan" not in (tmp_path / "output/protein_supported.gff3").read_text(
        encoding="utf-8"
    )


def test_cds_with_unknown_parent_remains_invalid(tmp_path: Path) -> None:
    gff = _write(
        tmp_path / "input.gff3",
        _gff() + "chr1\tx\tCDS\t700\t750\t.\t+\t0\tID=bad;Parent=missing\n",
    )
    proteins = _write(tmp_path / "proteins.fa", ">p1 transcript:t1\nMMMM\n")
    with pytest.raises(ValueError, match="unresolved_CDS_with_parent=1"):
        build_protein_supported_gene_universe(
            primary_gff_path=gff,
            provider_protein_path=proteins,
            output_dir=tmp_path / "output",
            species_id="fixture",
            holdout_id="holdout",
            policy_id="policy",
        )
