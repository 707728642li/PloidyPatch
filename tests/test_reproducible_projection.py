from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from ploidypatch.reproducible_projection import (
    compare_decision_tables,
    filter_candidate_gff_by_upstream_models,
    write_comparison_audit,
    verify_tree_manifest,
)
from ploidypatch.artifact_manifest import sha256_file


FIELDS = ("model_id", "source", "seqid", "status", "reason")


def _write_decisions(path: Path, rows: list[tuple[str, ...]], fields=FIELDS) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(fields)
        writer.writerows(rows)


def test_exact_complete_rows_define_stable_models(tmp_path: Path) -> None:
    a = tmp_path / "a.tsv"
    b = tmp_path / "b.tsv"
    _write_decisions(
        a,
        [
            ("stable", "lifton", "chr1", "accepted", "accepted"),
            ("changed", "lifton", "chr1", "accepted", "accepted"),
            ("only_a", "lifton", "chr1", "accepted", "accepted"),
            ("rejected", "lifton", "chr1", "rejected", "overlap"),
        ],
    )
    _write_decisions(
        b,
        [
            ("stable", "lifton", "chr1", "accepted", "accepted"),
            ("changed", "lifton", "chr2", "accepted", "accepted"),
            ("only_b", "lifton", "chr1", "accepted", "accepted"),
            ("rejected", "lifton", "chr1", "rejected", "overlap"),
        ],
    )
    comparison = compare_decision_tables(a, b)
    assert comparison.stable_models == {"stable", "rejected"}
    assert comparison.stable_accepted_models == {"stable"}
    audit = tmp_path / "audit.tsv"
    manifest = tmp_path / "audit.json"
    payload = write_comparison_audit(
        comparison=comparison,
        decisions_a=a,
        decisions_b=b,
        output_tsv=audit,
        output_json=manifest,
        method="lifton",
        reference="ref",
    )
    assert payload["unstable_models"] == 3
    assert json.loads(manifest.read_text())["label_access"] is False
    assert len(list(csv.DictReader(audit.open(), delimiter="\t"))) == 5


def test_truth_bearing_decision_table_is_rejected(tmp_path: Path) -> None:
    a = tmp_path / "a.tsv"
    b = tmp_path / "b.tsv"
    fields = ("model_id", "status", "label")
    _write_decisions(a, [("x", "accepted", "1")], fields)
    _write_decisions(b, [("x", "accepted", "1")], fields)
    with pytest.raises(ValueError, match="Truth-bearing"):
        compare_decision_tables(a, b)


def test_filter_preserves_only_complete_allowed_gene_blocks(tmp_path: Path) -> None:
    source = tmp_path / "candidate_only.gff3"
    source.write_text(
        "##gff-version 3\n"
        "chr1\tP\tgene\t1\t20\t.\t+\t.\tID=g1;upstream_model=keep\n"
        "chr1\tP\tmRNA\t1\t20\t.\t+\t.\tID=t1;Parent=g1;upstream_model=keep\n"
        "chr1\tP\texon\t1\t20\t.\t+\t.\tID=e1;Parent=t1\n"
        "chr1\tP\tCDS\t1\t18\t.\t+\t0\tID=c1;Parent=t1\n"
        "chr1\tP\tgene\t30\t50\t.\t+\t.\tID=g2;upstream_model=drop\n"
        "chr1\tP\tmRNA\t30\t50\t.\t+\t.\tID=t2;Parent=g2;upstream_model=drop\n"
        "chr1\tP\texon\t30\t50\t.\t+\t.\tID=e2;Parent=t2\n"
        "chr1\tP\tCDS\t30\t48\t.\t+\t0\tID=c2;Parent=t2\n",
        encoding="utf-8",
    )
    output = tmp_path / "stable.gff3"
    counts = filter_candidate_gff_by_upstream_models(
        source=source, allowed_models={"keep"}, output=output
    )
    text = output.read_text()
    assert "upstream_model=keep" in text
    assert "upstream_model=drop" not in text
    assert counts == {
        "source_models": 2,
        "allowed_models": 1,
        "kept_models": 1,
        "dropped_models": 1,
    }


def test_filter_rejects_missing_allowed_model(tmp_path: Path) -> None:
    source = tmp_path / "candidate_only.gff3"
    source.write_text(
        "##gff-version 3\n"
        "chr1\tP\tgene\t1\t20\t.\t+\t.\tID=g1;upstream_model=present\n"
        "chr1\tP\tmRNA\t1\t20\t.\t+\t.\tID=t1;Parent=g1;upstream_model=present\n"
        "chr1\tP\tCDS\t1\t18\t.\t+\t0\tID=c1;Parent=t1\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="absent"):
        filter_candidate_gff_by_upstream_models(
            source=source, allowed_models={"missing"}, output=tmp_path / "out.gff3"
        )


def test_filter_accepts_miniprot_lineage_attribute(tmp_path: Path) -> None:
    source = tmp_path / "candidate_only.gff3"
    source.write_text(
        "##gff-version 3\n"
        "chr1\tP\tgene\t1\t20\t.\t+\t.\tID=g1;miniprot_model=MP1\n"
        "chr1\tP\tmRNA\t1\t20\t.\t+\t.\tID=t1;Parent=g1;miniprot_model=MP1\n"
        "chr1\tP\tCDS\t1\t18\t.\t+\t0\tID=c1;Parent=t1\n",
        encoding="utf-8",
    )
    output = tmp_path / "out.gff3"
    counts = filter_candidate_gff_by_upstream_models(
        source=source, allowed_models={"MP1"}, output=output
    )
    assert counts["kept_models"] == 1
    assert "miniprot_model=MP1" in output.read_text()


def test_exact_tree_manifest_rejects_extra_or_changed_files(tmp_path: Path) -> None:
    root = tmp_path / "tree"
    root.mkdir()
    (root / "a.txt").write_bytes(b"a\n")
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text(
        "relative_path\tbytes\tsha256\n"
        f"a.txt\t2\t{sha256_file(root / 'a.txt')}\n",
        encoding="utf-8",
    )
    assert verify_tree_manifest(root=root, manifest=manifest) == {
        "files": 1,
        "bytes": 2,
    }
    (root / "extra.txt").write_bytes(b"x\n")
    with pytest.raises(ValueError, match="Unexpected"):
        verify_tree_manifest(root=root, manifest=manifest)
