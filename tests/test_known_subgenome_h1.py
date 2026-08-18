from __future__ import annotations

import csv
from pathlib import Path

import pytest

from ploidypatch.artifact_manifest import verify_sha256sums
from ploidypatch.known_subgenome_h1 import (
    PAIR_SCHEMA,
    SELF_FILTER_SCHEMA,
    attach_descriptive_yn00_ks,
    filter_self_pairs_by_known_subgenome,
    infer_exact_two_outgroup_pair_consistent_truth,
)


def write_tsv(
    path: Path, fields: tuple[str, ...], rows: list[tuple[object, ...]]
) -> Path:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(fields)
        writer.writerows(rows)
    return path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def homoeolog_table(path: Path) -> Path:
    return write_tsv(
        path,
        (
            "homoeolog_group",
            "c_subgenome_seqid",
            "e_subgenome_seqid",
            "c_label",
            "e_label",
        ),
        [("1", "C1", "E1", "1c", "1e"), ("2", "C2", "E2", "2c", "2e")],
    )


def test_known_subgenome_filter_is_exact_and_ks_independent(tmp_path: Path) -> None:
    groups = homoeolog_table(tmp_path / "groups.tsv")
    self_pairs = write_tsv(
        tmp_path / "self.tsv",
        (
            "pair_id",
            "gene_id_a",
            "gene_id_b",
            "seqid_a",
            "seqid_b",
            "support_block_count",
            "longest_block_pairs",
            "yn00_ks",
        ),
        [
            ("p1", "B", "A", "E1", "C1", 1, 22, "0.12"),
            ("p2", "C", "D", "C1", "E2", 2, 30, "0.20"),
            ("p3", "E", "F", "C1", "C1", 1, 25, "0.30"),
            ("p4", "G", "H", "unknown", "E1", 1, 24, "0.40"),
            ("p5", "I", "J", "C2", "E2", 3, 35, "99.0"),
        ],
    )
    output = tmp_path / "filtered"
    manifest = filter_self_pairs_by_known_subgenome(
        self_pairs_path=self_pairs,
        homoeolog_groups_path=groups,
        output_dir=output,
    )

    assert manifest["schema_version"] == SELF_FILTER_SCHEMA
    assert manifest["parameters"]["yn00_ks_policy"] == (
        "descriptive_only_not_used_for_selection"
    )
    assert manifest["counts"]["accepted_known_subgenome_pairs"] == 2
    accepted = read_tsv(output / "pairs.tsv")
    assert [(row["gene_id_a"], row["gene_id_b"]) for row in accepted] == [
        ("A", "B"),
        ("I", "J"),
    ]
    assert accepted[0]["subgenome_a"] == "c"
    assert accepted[0]["subgenome_b"] == "e"
    assert accepted[1]["yn00_ks"] == "99.0"
    reasons = {row["reason"] for row in read_tsv(output / "decisions.tsv")}
    assert reasons == {
        "predeclared_same_group_exact_C_E_pair",
        "different_predeclared_homoeolog_group",
        "same_subgenome_or_duplicate_chromosome",
        "nonprimary_or_unregistered_seqid",
    }
    verify_sha256sums(output, ignore_checksum_file=True)

    with pytest.raises(FileExistsError, match="overwrite"):
        filter_self_pairs_by_known_subgenome(
            self_pairs_path=self_pairs,
            homoeolog_groups_path=groups,
            output_dir=output,
        )


def test_known_subgenome_table_rejects_duplicate_chromosome_assignments(
    tmp_path: Path,
) -> None:
    groups = write_tsv(
        tmp_path / "bad_groups.tsv",
        (
            "homoeolog_group",
            "c_subgenome_seqid",
            "e_subgenome_seqid",
            "c_label",
            "e_label",
        ),
        [("1", "C1", "E1", "1c", "1e"), ("2", "C1", "E2", "2c", "2e")],
    )
    self_pairs = write_tsv(
        tmp_path / "self.tsv",
        (
            "pair_id",
            "gene_id_a",
            "gene_id_b",
            "seqid_a",
            "seqid_b",
            "support_block_count",
            "longest_block_pairs",
        ),
        [("p1", "A", "B", "C1", "E1", 1, 22)],
    )
    with pytest.raises(ValueError, match="Duplicate homoeolog"):
        filter_self_pairs_by_known_subgenome(
            self_pairs_path=self_pairs,
            homoeolog_groups_path=groups,
            output_dir=tmp_path / "out",
        )


def test_two_event_outgroups_require_pair_consistency_and_veto_alternates(
    tmp_path: Path,
) -> None:
    self_pairs = write_tsv(
        tmp_path / "self.tsv",
        ("pair_id", "gene_id_a", "gene_id_b", "homoeolog_group"),
        [("p1", "A", "B", "1"), ("p2", "I", "J", "2")],
    )
    gardenia = write_tsv(
        tmp_path / "gardenia.tsv",
        ("gene_id_a", "gene_id_b", "status", "reason"),
        [
            ("A", "B", "accepted", "exact_1_to_2"),
            ("B", "A", "accepted", "second_counterpart_same_pair"),
            ("I", "J", "accepted", "exact_1_to_2"),
        ],
    )
    ophiorrhiza = write_tsv(
        tmp_path / "ophiorrhiza.tsv",
        ("gene_id_a", "gene_id_b", "status", "reason"),
        [
            ("A", "B", "accepted", "exact_1_to_2"),
            ("A", "X", "rejected", "nonreciprocal_multiple_partners"),
            ("I", "J", "accepted", "exact_1_to_2"),
        ],
    )
    output = tmp_path / "truth"
    manifest = infer_exact_two_outgroup_pair_consistent_truth(
        known_subgenome_self_pairs_path=self_pairs,
        evaluator_group_decisions={
            "gardenia": gardenia,
            "ophiorrhiza": ophiorrhiza,
        },
        output_dir=output,
    )

    assert manifest["schema_version"] == PAIR_SCHEMA
    assert manifest["counts"]["accepted_truth_pairs"] == 1
    assert [(row["gene_id_a"], row["gene_id_b"]) for row in read_tsv(output / "pairs.tsv")] == [
        ("I", "J")
    ]
    rejected = next(
        row for row in read_tsv(output / "decisions.tsv") if row["gene_id_a"] == "A"
    )
    assert rejected["reason"] == "missing_or_discordant_group:ophiorrhiza"
    ophio_audit = manifest["counts"]["evaluator_groups"][1]
    assert ophio_audit["discordant_target_members"] == 1
    assert manifest["candidate_references_used_for_truth"] is False
    assert manifest["truth_labels_accessed"] is False
    verify_sha256sums(output, ignore_checksum_file=True)


def test_truth_rejects_missing_or_extra_evaluator_group(tmp_path: Path) -> None:
    self_pairs = write_tsv(
        tmp_path / "self.tsv",
        ("pair_id", "gene_id_a", "gene_id_b", "homoeolog_group"),
        [("p1", "A", "B", "1")],
    )
    group = write_tsv(
        tmp_path / "group.tsv",
        ("gene_id_a", "gene_id_b", "status", "reason"),
        [("A", "B", "accepted", "exact_1_to_2")],
    )
    with pytest.raises(ValueError, match="Gardenia and Ophiorrhiza"):
        infer_exact_two_outgroup_pair_consistent_truth(
            known_subgenome_self_pairs_path=self_pairs,
            evaluator_group_decisions={"gardenia": group},
            output_dir=tmp_path / "truth",
        )


def test_descriptive_yn00_attachment_never_filters_or_imputes(tmp_path: Path) -> None:
    pairs = write_tsv(
        tmp_path / "pairs.tsv",
        (
            "pair_id",
            "gene_id_a",
            "gene_id_b",
            "wgdi_gene_id_a",
            "wgdi_gene_id_b",
        ),
        [("p1", "A", "B", "wa", "wb"), ("p2", "C", "D", "wc", "wd")],
    )
    ks = write_tsv(
        tmp_path / "ks.tsv",
        ("id1", "id2", "ka_YN00", "ks_YN00"),
        [("wb", "wa", "0.1", "0.2")],
    )
    output = tmp_path / "with_ks.tsv"
    manifest = attach_descriptive_yn00_ks(
        self_pairs_path=pairs, ks_path=ks, output_path=output
    )
    rows = read_tsv(output)
    assert [row["pair_id"] for row in rows] == ["p1", "p2"]
    assert [row["yn00_ks"] for row in rows] == ["0.2", ""]
    assert manifest["selection_use"] is False
    assert manifest["counts"]["matched_self_pairs"] == 1
    assert manifest["counts"]["missing_self_pairs"] == 1
