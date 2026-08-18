from __future__ import annotations

import csv
from pathlib import Path

import pytest

from ploidypatch.candidate_merge import merge_candidate_gffs
from ploidypatch.gff import parse_attributes


def _write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


def test_merge_candidate_gffs_namespaces_hierarchies_and_preserves_votes(
    tmp_path: Path,
) -> None:
    first = _write(
        tmp_path / "a.gff3",
        "##gff-version 3\n"
        "chr1\ta\tgene\t1\t30\t.\t+\t.\tID=g1\n"
        "chr1\ta\tmRNA\t1\t30\t.\t+\t.\tID=t1;Parent=g1;pAA=0.9;stop=*\n"
        "chr1\ta\tCDS\t1\t30\t.\t+\t0\tParent=t1\n",
    )
    second = _write(
        tmp_path / "b.gff3",
        "##gff-version 3\n"
        "chr2\tb\tgene\t5\t34\t.\t-\t.\tID=g1\n"
        "chr2\tb\tmRNA\t5\t34\t.\t-\t.\tID=t1;Parent=g1;mutation=frameshift\n"
        "chr2\tb\tCDS\t5\t34\t.\t-\t0\tParent=t1\n"
        "##FASTA\n>chr2\nACGT\n",
    )
    output = tmp_path / "merged.gff3"
    provenance = tmp_path / "provenance.tsv"
    manifest = merge_candidate_gffs(
        candidate_inputs=(f"refA={first}", f"refB={second}"),
        output_gff_path=output,
        provenance_tsv_path=provenance,
    )
    observed = output.read_text(encoding="utf-8")
    assert "ID=refA__g1" in observed and "Parent=refA__g1" in observed
    assert "ID=refB__g1" in observed and "Parent=refB__g1" in observed
    assert "pAA=0.9" in observed and "mutation=frameshift" in observed
    assert ">chr2" not in observed
    assert manifest["method_family_policy"]["vote_count_after_merge"] == 1
    assert manifest["counts"] == {
        "features": 6,
        "ids": 4,
        "id_occurrences": 4,
        "repeated_id_lines": 0,
        "duplicate_feature_lines_skipped": 0,
        "reference_sources": 2,
    }
    with provenance.open(encoding="utf-8", newline="") as handle:
        assert len(list(csv.DictReader(handle, delimiter="\t"))) == 4
    for raw_line in observed.splitlines():
        if raw_line.startswith("#"):
            continue
        attributes, malformed = parse_attributes(raw_line.split("\t")[8])
        assert malformed == 0
        assert attributes


def test_merge_candidate_gffs_namespaces_one_reference_without_extra_vote(
    tmp_path: Path,
) -> None:
    candidate = _write(
        tmp_path / "single.gff3",
        "##gff-version 3\n"
        "chr1\tx\tgene\t1\t30\t.\t+\t.\tID=g1\n"
        "chr1\tx\tmRNA\t1\t30\t.\t+\t.\tID=t1;Parent=g1\n"
        "chr1\tx\tCDS\t1\t30\t.\t+\t0\tParent=t1\n",
    )
    output = tmp_path / "namespaced.gff3"
    provenance = tmp_path / "provenance.tsv"
    manifest = merge_candidate_gffs(
        candidate_inputs=(f"refA={candidate}",),
        output_gff_path=output,
        provenance_tsv_path=provenance,
    )
    observed = output.read_text(encoding="utf-8")
    assert "ID=refA__g1" in observed
    assert "ID=refA__t1;Parent=refA__g1" in observed
    assert "Parent=refA__t1" in observed
    assert manifest["counts"]["reference_sources"] == 1
    assert manifest["method_family_policy"]["vote_count_after_merge"] == 1


def test_merge_candidate_gffs_rejects_empty_input_list(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="At least one reference GFF input"):
        merge_candidate_gffs(
            candidate_inputs=(),
            output_gff_path=tmp_path / "empty.gff3",
            provenance_tsv_path=tmp_path / "empty.tsv",
        )


def test_merge_candidate_gffs_rejects_orphan_parent_without_partial_outputs(
    tmp_path: Path,
) -> None:
    orphan = _write(
        tmp_path / "orphan.gff3",
        "chr1\tx\tCDS\t1\t3\t.\t+\t0\tParent=missing\n",
    )
    valid = _write(
        tmp_path / "valid.gff3",
        "chr1\tx\tmRNA\t1\t3\t.\t+\t.\tID=t\n"
        "chr1\tx\tCDS\t1\t3\t.\t+\t0\tParent=t\n",
    )
    output = tmp_path / "merged.gff3"
    provenance = tmp_path / "provenance.tsv"
    with pytest.raises(ValueError, match="Parent references without IDs"):
        merge_candidate_gffs(
            candidate_inputs=(f"bad={orphan}", f"good={valid}"),
            output_gff_path=output,
            provenance_tsv_path=provenance,
        )
    assert not output.exists()
    assert not provenance.exists()


def test_merge_candidate_gffs_allows_discontinuous_ids_but_not_conflicts(
    tmp_path: Path,
) -> None:
    discontinuous = _write(
        tmp_path / "discontinuous.gff3",
        "chr1\tx\tgene\t1\t30\t.\t+\t.\tID=g1\n"
        "chr1\tx\tmRNA\t1\t30\t.\t+\t.\tID=t1;Parent=g1\n"
        "chr1\tx\tCDS\t1\t9\t.\t+\t0\tID=cds1;Parent=t1\n"
        "chr1\tx\tCDS\t20\t30\t.\t+\t0\tID=cds1;Parent=t1\n"
        "chr1\tx\tCDS\t20\t30\t.\t+\t0\tID=cds1;Parent=t1\n",
    )
    second = _write(
        tmp_path / "second.gff3",
        "chr2\ty\tgene\t1\t9\t.\t+\t.\tID=g1\n"
        "chr2\ty\tmRNA\t1\t9\t.\t+\t.\tID=t1;Parent=g1\n"
        "chr2\ty\tCDS\t1\t9\t.\t+\t0\tID=cds1;Parent=t1\n",
    )
    manifest = merge_candidate_gffs(
        candidate_inputs=(f"a={discontinuous}", f"b={second}"),
        output_gff_path=tmp_path / "merged.gff3",
        provenance_tsv_path=tmp_path / "provenance.tsv",
    )
    assert manifest["counts"]["ids"] == 6
    assert manifest["counts"]["id_occurrences"] == 7
    assert manifest["counts"]["repeated_id_lines"] == 1
    assert manifest["counts"]["duplicate_feature_lines_skipped"] == 1
    assert manifest["inputs"]["a"]["repeated_id_lines"] == 1

    conflict = _write(
        tmp_path / "conflict.gff3",
        "chr3\tz\tmRNA\t1\t9\t.\t+\t.\tID=t1\n"
        "chr3\tz\tmRNA\t20\t30\t.\t+\t.\tID=t2\n"
        "chr3\tz\tCDS\t1\t9\t.\t+\t0\tID=cds1;Parent=t1\n"
        "chr3\tz\tCDS\t20\t30\t.\t+\t0\tID=cds1;Parent=t2\n",
    )
    with pytest.raises(ValueError, match="Conflicting repeated ID"):
        merge_candidate_gffs(
            candidate_inputs=(f"bad={conflict}", f"good={second}"),
            output_gff_path=tmp_path / "conflict_merged.gff3",
            provenance_tsv_path=tmp_path / "conflict_provenance.tsv",
        )
