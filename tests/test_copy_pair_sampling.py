from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from ploidypatch.cli import main
from ploidypatch.copy_pair_sampling import sample_balanced_copy_pairs
from ploidypatch.perturb import COPY_COLLAPSE_EVENT
from ploidypatch.structure_perturb import (
    copy_collapse_partner_ids,
    generate_structure_benchmark,
)


def _gene(gene: str, seqid: str, start: int, segments: int) -> str:
    transcript = f"{gene}.t1"
    end = start + segments * 60 - 1
    lines = [
        f"{seqid}\ttest\tgene\t{start}\t{end}\t.\t+\t.\tID={gene}\n",
        f"{seqid}\ttest\tmRNA\t{start}\t{end}\t.\t+\t.\t"
        f"ID={transcript};Parent={gene}\n",
    ]
    for index in range(segments):
        left = start + index * 60
        right = left + 29
        lines.append(
            f"{seqid}\ttest\tCDS\t{left}\t{right}\t.\t+\t0\t"
            f"ID={gene}.cds{index + 1};Parent={transcript}\n"
        )
    return "".join(lines)


def _write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8", newline="")
    return path


def _read(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_balanced_pair_sample_matches_perturbation_partner_choice(
    tmp_path: Path,
) -> None:
    source = _write(
        tmp_path / "source.gff3",
        "##gff-version 3\n"
        + _gene("A1", "chr1", 1, 1)
        + _gene("A2", "chr2", 1, 1)
        + _gene("B1", "chr3", 1, 2)
        + _gene("B2", "chr4", 1, 2)
        + _gene("C1", "chr5", 1, 4)
        + _gene("C2", "chr6", 1, 4)
        + _gene("D1", "chr7", 1, 7)
        + _gene("D2", "chr8", 1, 7),
    )
    pairs = _write(
        tmp_path / "pairs.tsv",
        "pair_id\tgene_id_a\tgene_id_b\n"
        "P1\tA1\tA2\n"
        "P2\tB1\tB2\n"
        "P3\tC1\tC2\n"
        "P4\tD1\tD2\n",
    )
    selected_path = tmp_path / "selected.tsv"
    decisions = tmp_path / "decisions.tsv"
    manifest = sample_balanced_copy_pairs(
        source_gff_path=source,
        pair_tsv_path=pairs,
        output_pair_tsv_path=selected_path,
        decisions_tsv_path=decisions,
        count=3,
        seed=29,
    )

    selected = _read(selected_path)
    assert len(selected) == 3
    assert manifest["counts"]["selected_pairs"] == 3
    assert manifest["counts"]["selection_shortfall"] == 0
    assert len(manifest["counts"]["selected_by_complexity"]) == 3
    assert all(row["collapsed_gene_id"] != row["retained_gene_id"] for row in selected)

    benchmark = tmp_path / "benchmark"
    generate_structure_benchmark(
        source,
        benchmark,
        event_type=COPY_COLLAPSE_EVENT,
        count=3,
        seed=29,
        pair_tsv_path=selected_path,
    )
    truth = json.loads((benchmark / "hidden_truth.json").read_text(encoding="utf-8"))
    assert {event["target"]["gene_id"] for event in truth["events"]} == {
        row["collapsed_gene_id"] for row in selected
    }


def test_balanced_pair_sample_cli_reports_shortfall(tmp_path: Path) -> None:
    source = _write(
        tmp_path / "source.gff3",
        "##gff-version 3\n"
        + _gene("A1", "chr1", 1, 1)
        + _gene("A2", "chr2", 1, 1),
    )
    pairs = _write(
        tmp_path / "pairs.tsv",
        "gene_id_a\tgene_id_b\nA1\tA2\n",
    )
    output = tmp_path / "selected.tsv"
    assert main(
        [
            "benchmark",
            "sample-copy-pairs",
            "--source-gff",
            str(source),
            "--pairs",
            str(pairs),
            "--count",
            "2",
            "--seed",
            "29",
            "--output-pairs",
            str(output),
            "--decisions-tsv",
            str(tmp_path / "decisions.tsv"),
        ]
    ) == 0
    manifest = json.loads(
        Path(str(output) + ".manifest.json").read_text(encoding="utf-8")
    )
    assert manifest["counts"]["selected_pairs"] == 1
    assert manifest["counts"]["selection_shortfall"] == 1


def test_balanced_pair_sample_rejects_nonunique_collapsed_phased_cds_chain(
    tmp_path: Path,
) -> None:
    collapsed, _retained = copy_collapse_partner_ids("A1", "A2", seed=29)
    collapsed_seqid = {"A1": "chr1", "A2": "chr2"}[collapsed]
    source = _write(
        tmp_path / "source.gff3",
        "##gff-version 3\n"
        + _gene("A1", "chr1", 1, 1)
        + _gene("A2", "chr2", 1, 1)
        + _gene("DUP", collapsed_seqid, 1, 1),
    )
    pairs = _write(
        tmp_path / "pairs.tsv",
        "pair_id\tgene_id_a\tgene_id_b\nP1\tA1\tA2\n",
    )
    selected_path = tmp_path / "selected.tsv"
    decisions_path = tmp_path / "decisions.tsv"
    manifest = sample_balanced_copy_pairs(
        source_gff_path=source,
        pair_tsv_path=pairs,
        output_pair_tsv_path=selected_path,
        decisions_tsv_path=decisions_path,
        count=1,
        seed=29,
    )

    assert _read(selected_path) == []
    assert _read(decisions_path)[0]["reason"] == (
        "collapsed_phased_cds_chain_not_unique_in_source"
    )
    assert manifest["parameters"][
        "require_unique_collapsed_phased_cds_chain_in_source"
    ] is True
    assert manifest["counts"]["operable_cross_seqid_pairs"] == 0
    assert manifest["counts"]["selected_pairs"] == 0
    assert manifest["counts"]["selection_shortfall"] == 1


def test_balanced_pair_sample_can_round_robin_predeclared_groups(
    tmp_path: Path,
) -> None:
    source = _write(
        tmp_path / "source.gff3",
        "##gff-version 3\n"
        + _gene("A1", "chr1", 1, 1)
        + _gene("A2", "chr2", 1, 1)
        + _gene("B1", "chr3", 1, 2)
        + _gene("B2", "chr4", 1, 2)
        + _gene("C1", "chr5", 1, 4)
        + _gene("C2", "chr6", 1, 4)
        + _gene("D1", "chr7", 1, 7)
        + _gene("D2", "chr8", 1, 7),
    )
    pairs = _write(
        tmp_path / "pairs.tsv",
        "pair_id\tgene_id_a\tgene_id_b\thomoeolog_group\n"
        "P1\tA1\tA2\t1\n"
        "P2\tB1\tB2\t1\n"
        "P3\tC1\tC2\t2\n"
        "P4\tD1\tD2\t2\n",
    )
    selected = tmp_path / "selected.tsv"
    manifest = sample_balanced_copy_pairs(
        source_gff_path=source,
        pair_tsv_path=pairs,
        output_pair_tsv_path=selected,
        decisions_tsv_path=tmp_path / "decisions.tsv",
        count=2,
        seed=29,
        balance_group_field="homoeolog_group",
    )
    assert {row["homoeolog_group"] for row in _read(selected)} == {"1", "2"}
    assert manifest["counts"]["selected_by_balance_group"] == {"1": 1, "2": 1}
    assert manifest["parameters"]["balance_group_policy"] == (
        "declared_group_global_round_robin"
    )


def test_balanced_pair_sample_rejects_missing_declared_group_field(
    tmp_path: Path,
) -> None:
    source = _write(
        tmp_path / "source.gff3",
        "##gff-version 3\n" + _gene("A1", "chr1", 1, 1) + _gene("A2", "chr2", 1, 1),
    )
    pairs = _write(
        tmp_path / "pairs.tsv", "gene_id_a\tgene_id_b\nA1\tA2\n"
    )
    with pytest.raises(ValueError, match="lacks balance group field"):
        sample_balanced_copy_pairs(
            source_gff_path=source,
            pair_tsv_path=pairs,
            output_pair_tsv_path=tmp_path / "selected.tsv",
            decisions_tsv_path=tmp_path / "decisions.tsv",
            count=1,
            seed=29,
            balance_group_field="homoeolog_group",
        )
