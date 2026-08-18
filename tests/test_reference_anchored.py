from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from ploidypatch.cli import build_parser, main
from ploidypatch.reference_anchored import (
    CANDIDATE_SOURCE_FIELDS,
    HOMELOG_BASE_FIELDS,
    WGD_PAIR_FIELDS,
    aggregate_reference_anchored_projection,
)


def write_tsv(
    path: Path, fields: tuple[str, ...], rows: list[dict[str, object]]
) -> Path:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)
    return path


def inputs(tmp_path: Path) -> tuple[Path, Path, Path, Path, Path]:
    candidates = write_tsv(
        tmp_path / "candidate_sources.tsv",
        CANDIDATE_SOURCE_FIELDS,
        [
            {"candidate_digest": "d1", "candidate_gene_id": "c1", "arm_id": "m1", "reference_id": "refA", "source_gene_id": "s1"},
            {"candidate_digest": "d1", "candidate_gene_id": "c1", "arm_id": "m2", "reference_id": "refA", "source_gene_id": "s2"},
            {"candidate_digest": "d2", "candidate_gene_id": "c2", "arm_id": "m5", "reference_id": "refA", "source_gene_id": "s4"},
            {"candidate_digest": "d2", "candidate_gene_id": "c2", "arm_id": "m3", "reference_id": "refB", "source_gene_id": "t1"},
            {"candidate_digest": "d3", "candidate_gene_id": "c3", "arm_id": "m6", "reference_id": "refA", "source_gene_id": "s5"},
            {"candidate_digest": "d4", "candidate_gene_id": "c4", "arm_id": "m4", "reference_id": "refA", "source_gene_id": "s3"},
        ],
    )
    pairs_a = write_tsv(
        tmp_path / "pairs_a.tsv",
        WGD_PAIR_FIELDS,
        [
            {"source_gene_id": "s1", "partner_source_gene_id": "h1", "support_block_count": 2, "longest_block_pairs": 30},
            {"source_gene_id": "h1", "partner_source_gene_id": "s1", "support_block_count": 2, "longest_block_pairs": 30},
            {"source_gene_id": "s3", "partner_source_gene_id": "h3", "support_block_count": 1, "longest_block_pairs": 25},
            {"source_gene_id": "h3", "partner_source_gene_id": "s3", "support_block_count": 1, "longest_block_pairs": 25},
            {"source_gene_id": "s4", "partner_source_gene_id": "h4", "support_block_count": 3, "longest_block_pairs": 35},
            {"source_gene_id": "h4", "partner_source_gene_id": "s4", "support_block_count": 3, "longest_block_pairs": 35},
        ],
    )
    pairs_b = write_tsv(
        tmp_path / "pairs_b.tsv",
        WGD_PAIR_FIELDS,
        [
            {"source_gene_id": "t1", "partner_source_gene_id": "u1", "support_block_count": 4, "longest_block_pairs": 40},
            {"source_gene_id": "u1", "partner_source_gene_id": "t1", "support_block_count": 4, "longest_block_pairs": 40},
        ],
    )
    evidence_a = write_tsv(
        tmp_path / "evidence_a.tsv",
        HOMELOG_BASE_FIELDS,
        [
            {"arm_id": "m1", "source_homeolog_gene_id": "h1", "evidence_status": "evidence", "target_base_gene_id": "base1", "evidence_reason": "unique_projection"},
            {"arm_id": "m4", "source_homeolog_gene_id": "h3", "evidence_status": "conflict", "target_base_gene_id": "", "evidence_reason": "multiple_base_genes"},
            {"arm_id": "m5", "source_homeolog_gene_id": "h4", "evidence_status": "evidence", "target_base_gene_id": "base1", "evidence_reason": "unique_projection"},
        ],
    )
    evidence_b = write_tsv(
        tmp_path / "evidence_b.tsv",
        HOMELOG_BASE_FIELDS,
        [
            {"arm_id": "m3", "source_homeolog_gene_id": "u1", "evidence_status": "evidence", "target_base_gene_id": "base2", "evidence_reason": "unique_projection"},
        ],
    )
    return candidates, pairs_a, pairs_b, evidence_a, evidence_b


def read_index(path: Path, key: str) -> dict[str, dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return {row[key]: row for row in csv.DictReader(handle, delimiter="\t")}


def test_reference_anchored_missing_does_not_veto_and_conflicts_do(
    tmp_path: Path,
) -> None:
    candidates, pairs_a, pairs_b, evidence_a, evidence_b = inputs(tmp_path)
    output = tmp_path / "output"
    report = aggregate_reference_anchored_projection(
        candidate_source_provenance_path=candidates,
        accepted_wgd_pair_inputs=[f"refA={pairs_a}", f"refB={pairs_b}"],
        homeolog_base_evidence_inputs=[
            f"refA={evidence_a}", f"refB={evidence_b}"
        ],
        output_dir_path=output,
    )
    assert report["truth_access"] is False
    assert report["labels_used"] is False
    assert report["candidate_or_target_union_graph_used"] is False
    assert report["counts"] == {
        "candidates": 4,
        "candidate_source_arms": 6,
        "accepted": 1,
        "rejected": 3,
        "reasons": {
            "abstain_conflicting_evidence_arm": 1,
            "abstain_discordant_target_base_partners": 1,
            "abstain_no_evidence_arm": 1,
            "accepted_reference_anchored_unique_consensus": 1,
        },
    }
    selection = read_index(output / "selection.tsv", "consensus_digest")
    assert selection["d1"]["status"] == "accepted"
    assert selection["d1"]["partner_gene_id"] == "base1"
    assert selection["d1"]["anchor_cell_count"] == "1"
    assert selection["d1"]["compatible_partner_count"] == "1"
    assert selection["d2"]["reason"] == "abstain_discordant_target_base_partners"
    assert selection["d2"]["compatible_partner_count"] == "2"
    assert selection["d3"]["reason"] == "abstain_no_evidence_arm"
    assert selection["d3"]["compatible_partner_count"] == "0"
    assert selection["d4"]["reason"] == "abstain_conflicting_evidence_arm"
    assert selection["d4"]["compatible_partner_count"] == "0"
    checksum_rows = (output / "SHA256SUMS").read_text().splitlines()
    assert [row.split("  ", 1)[1] for row in checksum_rows] == [
        "arm_decisions.tsv",
        "manifest.json",
        "selection.tsv",
    ]
    manifest = json.loads((output / "manifest.json").read_text())
    assert manifest["backbone_id"].startswith("PPRA-BB-")
    assert manifest["backbone_id"] == (
        "PPRA-BB-" + manifest["backbone_definition_sha256"][:24]
    )
    assert {row["backbone_id"] for row in selection.values()} == {
        manifest["backbone_id"]
    }
    assert manifest["policy"]["missing_evidence_arms_veto"] is False
    assert manifest["policy"]["conflicting_evidence_arm_veto"] is True


def test_reference_anchored_backbone_and_pair_ids_bind_inputs(tmp_path: Path) -> None:
    candidates, pairs_a, pairs_b, evidence_a, evidence_b = inputs(tmp_path)
    output1 = tmp_path / "output1"
    aggregate_reference_anchored_projection(
        candidate_source_provenance_path=candidates,
        accepted_wgd_pair_inputs=[f"refA={pairs_a}", f"refB={pairs_b}"],
        homeolog_base_evidence_inputs=[f"refA={evidence_a}", f"refB={evidence_b}"],
        output_dir_path=output1,
    )
    manifest1 = json.loads((output1 / "manifest.json").read_text())
    selection1 = read_index(output1 / "selection.tsv", "consensus_digest")

    evidence_rows = list(csv.DictReader(evidence_a.open(newline=""), delimiter="\t"))
    evidence_rows[0]["target_base_gene_id"] = "base9"
    changed_evidence_a = write_tsv(
        tmp_path / "evidence_a_changed.tsv", HOMELOG_BASE_FIELDS, evidence_rows
    )
    output2 = tmp_path / "output2"
    aggregate_reference_anchored_projection(
        candidate_source_provenance_path=candidates,
        accepted_wgd_pair_inputs=[f"refA={pairs_a}", f"refB={pairs_b}"],
        homeolog_base_evidence_inputs=[
            f"refA={changed_evidence_a}",
            f"refB={evidence_b}",
        ],
        output_dir_path=output2,
    )
    manifest2 = json.loads((output2 / "manifest.json").read_text())
    selection2 = read_index(output2 / "selection.tsv", "consensus_digest")

    assert manifest1["backbone_id"] != manifest2["backbone_id"]
    assert manifest1["backbone_definition_sha256"] != manifest2["backbone_definition_sha256"]
    assert selection1["d1"]["pair_id"] != selection2["d1"]["pair_id"]


def test_reference_anchored_rejects_source_arm_reused_across_candidates(
    tmp_path: Path,
) -> None:
    candidates, pairs_a, pairs_b, evidence_a, evidence_b = inputs(tmp_path)
    candidate_rows = list(csv.DictReader(candidates.open(newline=""), delimiter="\t"))
    candidate_rows.append(
        {
            "candidate_digest": "d5",
            "candidate_gene_id": "c5",
            "arm_id": "m1",
            "reference_id": "refA",
            "source_gene_id": "s1",
        }
    )
    write_tsv(candidates, CANDIDATE_SOURCE_FIELDS, candidate_rows)
    with pytest.raises(ValueError, match="reused by different candidates"):
        aggregate_reference_anchored_projection(
            candidate_source_provenance_path=candidates,
            accepted_wgd_pair_inputs=[f"refA={pairs_a}", f"refB={pairs_b}"],
            homeolog_base_evidence_inputs=[
                f"refA={evidence_a}", f"refB={evidence_b}"
            ],
            output_dir_path=tmp_path / "output",
        )


def test_reference_anchored_is_nonoverwriting(tmp_path: Path) -> None:
    candidates, pairs_a, pairs_b, evidence_a, evidence_b = inputs(tmp_path)
    kwargs = dict(
        candidate_source_provenance_path=candidates,
        accepted_wgd_pair_inputs=[f"refA={pairs_a}", f"refB={pairs_b}"],
        homeolog_base_evidence_inputs=[f"refA={evidence_a}", f"refB={evidence_b}"],
        output_dir_path=tmp_path / "output",
    )
    aggregate_reference_anchored_projection(**kwargs)
    with pytest.raises(FileExistsError, match="refusing to overwrite"):
        aggregate_reference_anchored_projection(**kwargs)


def test_reference_anchored_rejects_nonreciprocal_pairs(tmp_path: Path) -> None:
    candidates, pairs_a, pairs_b, evidence_a, evidence_b = inputs(tmp_path)
    rows = list(csv.DictReader(pairs_a.open(newline=""), delimiter="\t"))[:-1]
    write_tsv(pairs_a, WGD_PAIR_FIELDS, rows)
    with pytest.raises(ValueError, match="not exact reciprocal"):
        aggregate_reference_anchored_projection(
            candidate_source_provenance_path=candidates,
            accepted_wgd_pair_inputs=[f"refA={pairs_a}", f"refB={pairs_b}"],
            homeolog_base_evidence_inputs=[f"refA={evidence_a}", f"refB={evidence_b}"],
            output_dir_path=tmp_path / "output",
        )


def test_reference_anchored_cli_contract() -> None:
    args = build_parser().parse_args(
        [
            "evidence",
            "aggregate-reference-anchored",
            "--candidate-source-provenance",
            "candidates.tsv",
            "--accepted-wgd-pairs",
            "refA=pairs.tsv",
            "--homeolog-base-evidence",
            "refA=evidence.tsv",
            "--output-dir",
            "output",
        ]
    )
    assert args.evidence_command == "aggregate-reference-anchored"
    assert args.evidence_type == "reference_anchored_projection"


def test_reference_anchored_cli_executes(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    candidates, pairs_a, pairs_b, evidence_a, evidence_b = inputs(tmp_path)
    output = tmp_path / "cli_output"
    status = main(
        [
            "evidence",
            "aggregate-reference-anchored",
            "--candidate-source-provenance",
            str(candidates),
            "--accepted-wgd-pairs",
            f"refA={pairs_a}",
            "--accepted-wgd-pairs",
            f"refB={pairs_b}",
            "--homeolog-base-evidence",
            f"refA={evidence_a}",
            "--homeolog-base-evidence",
            f"refB={evidence_b}",
            "--output-dir",
            str(output),
        ]
    )
    assert status == 0
    printed = json.loads(capsys.readouterr().out)
    assert printed["counts"]["accepted"] == 1
    assert (output / "manifest.json").is_file()
