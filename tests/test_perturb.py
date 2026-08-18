from __future__ import annotations

import json
from pathlib import Path

import pytest

from ploidypatch.perturb import (
    find_missing_gene_candidates,
    generate_missing_gene_benchmark,
    read_gff_document,
    restore_gff_from_truth,
)


GFF_FIXTURE = (
    "##gff-version 3\n"
    "# retained comment\n"
    "chr1\ttest\tgene\t1\t90\t.\t+\t.\tID=gene:G1\n"
    "chr1\ttest\tmRNA\t1\t90\t.\t+\t.\tID=transcript:T1;Parent=gene:G1\n"
    "chr1\ttest\texon\t1\t30\t.\t+\t.\tID=exon:E1;Parent=transcript:T1\n"
    "chr1\ttest\tCDS\t1\t30\t.\t+\t0\tID=CDS:C1;Parent=transcript:T1\n"
    "chr1\ttest\tgene\t101\t190\t.\t-\t.\tID=gene:G2\n"
    "chr1\ttest\tmRNA\t101\t190\t.\t-\t.\tID=transcript:T2;Parent=gene:G2\n"
    "chr1\ttest\texon\t161\t190\t.\t-\t.\tID=exon:E2;Parent=transcript:T2\n"
    "chr1\ttest\tCDS\t161\t190\t.\t-\t0\tID=CDS:C2;Parent=transcript:T2\n"
)


def write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8", newline="")
    return path


def test_missing_gene_candidates_require_unshared_coding_hierarchy(
    tmp_path: Path,
) -> None:
    gff = write(tmp_path / "source.gff3", GFF_FIXTURE)
    candidates = find_missing_gene_candidates(read_gff_document(gff))

    assert [item.gene.feature_id for item in candidates] == ["gene:G1", "gene:G2"]
    assert candidates[0].feature_type_counts == {
        "CDS": 1,
        "exon": 1,
        "gene": 1,
        "mRNA": 1,
    }


def test_perturbation_is_deterministic_hidden_and_exactly_reversible(
    tmp_path: Path,
) -> None:
    source = write(tmp_path / "source.gff3", GFF_FIXTURE)
    first = tmp_path / "first"
    second = tmp_path / "second"

    manifest1 = generate_missing_gene_benchmark(source, first, count=1, seed=17)
    manifest2 = generate_missing_gene_benchmark(source, second, count=1, seed=17)

    assert (first / "perturbed.gff3").read_bytes() == (
        second / "perturbed.gff3"
    ).read_bytes()
    assert (first / "hidden_truth.json").read_bytes() == (
        second / "hidden_truth.json"
    ).read_bytes()
    assert manifest1["perturbation"] == manifest2["perturbation"]

    truth = json.loads((first / "hidden_truth.json").read_text(encoding="utf-8"))
    assert len(truth["events"]) == 1
    target_gene = truth["events"][0]["target"]["gene_id"]
    perturbed_text = (first / "perturbed.gff3").read_text(encoding="utf-8")
    assert target_gene not in perturbed_text
    assert "# retained comment" in perturbed_text

    restored = first / "restored.gff3"
    report = restore_gff_from_truth(
        first / "perturbed.gff3", first / "hidden_truth.json", restored
    )
    assert restored.read_bytes() == source.read_bytes()
    assert report["events"] == 1
    assert report["restored_records"] == 4


def test_shared_feature_makes_both_gene_hierarchies_ineligible(tmp_path: Path) -> None:
    source = write(
        tmp_path / "shared.gff3",
        GFF_FIXTURE
        + "chr1\ttest\texon\t50\t60\t.\t+\t.\t"
        "ID=exon:shared;Parent=transcript:T1,transcript:T2\n",
    )
    candidates = find_missing_gene_candidates(read_gff_document(source))
    assert candidates == []


def test_generation_refuses_overwrite_and_insufficient_candidates(
    tmp_path: Path,
) -> None:
    source = write(tmp_path / "source.gff3", GFF_FIXTURE)
    with pytest.raises(ValueError, match="only 2 genes"):
        generate_missing_gene_benchmark(source, tmp_path / "too_many", count=3, seed=1)

    output = tmp_path / "run"
    generate_missing_gene_benchmark(source, output, count=1, seed=1)
    with pytest.raises(FileExistsError, match="Refusing to overwrite"):
        generate_missing_gene_benchmark(source, output, count=1, seed=1)


def test_restore_detects_modified_perturbed_input(tmp_path: Path) -> None:
    source = write(tmp_path / "source.gff3", GFF_FIXTURE)
    output = tmp_path / "run"
    generate_missing_gene_benchmark(source, output, count=1, seed=1)
    with (output / "perturbed.gff3").open("a", encoding="utf-8") as handle:
        handle.write("# unexpected\n")

    with pytest.raises(ValueError, match="line count"):
        restore_gff_from_truth(
            output / "perturbed.gff3",
            output / "hidden_truth.json",
            output / "restored.gff3",
        )
