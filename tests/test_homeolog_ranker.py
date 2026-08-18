from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

from ploidypatch.cli import main
from ploidypatch.copy_features import COPY_FEATURE_SCHEMA_VERSION, FEATURE_COLUMNS
from ploidypatch.copy_model import fit_copy_feature_contract
from ploidypatch.homeolog_ranker import (
    HOMEOLOG_RANKER_SCHEMA_VERSION,
    HOMEOLOG_REVIEW_RANKING_SCHEMA_VERSION,
    TOPOLOGY_ADDON_FIELDS,
    freeze_homeolog_review_rankings,
    score_homeolog_copy_candidates,
)
from ploidypatch.homeolog_topology import HOMEOLOG_TOPOLOGY_SCHEMA_VERSION


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _write_tsv(
    path: Path, fieldnames: list[str] | tuple[str, ...], rows: list[dict[str, str]]
) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def _fixture(tmp_path: Path) -> tuple[Path, Path, Path]:
    copy_path = tmp_path / "copy.tsv"
    copy_rows: list[dict[str, str]] = []
    for digest, start, wgd in (("D1", "1", "1"), ("D2", "101", "0")):
        row = {field: "" for field in FEATURE_COLUMNS}
        row.update(
            {
                "candidate_digest": digest,
                "seqid": "chr1",
                "start": start,
                "end": str(int(start) + 89),
                "strand": "+",
                "span_bp": "90",
                "cds_segments": "1",
                "cds_bp": "90",
                "support_method_count": "1",
                "support_methods": "miniprot",
                "has_miniprot": "1",
                "has_gemoma": "0",
                "has_lifton": "0",
                "wgd_existing_partner": wgd,
                "wgd_support_block_count": "1" if wgd == "1" else "",
                "wgd_longest_block_pairs": "20" if wgd == "1" else "",
            }
        )
        copy_rows.append(row)
    _write_tsv(copy_path, FEATURE_COLUMNS, copy_rows)
    copy_manifest = {
        "schema_version": COPY_FEATURE_SCHEMA_VERSION,
        "truth_access": False,
        "outputs": {"features": {"sha256": _sha256(copy_path), "rows": 2}},
    }
    Path(str(copy_path) + ".manifest.json").write_text(
        json.dumps(copy_manifest), encoding="utf-8"
    )

    topology_path = tmp_path / "topology.tsv"
    topology_fields = ["candidate_digest", *TOPOLOGY_ADDON_FIELDS]
    topology_rows = [
        {
            "candidate_digest": "D1",
            "topology_available": "1",
            "cds_bp_ratio": "0.9",
            "cds_segment_count_ratio": "1",
            "phase_lcs_similarity": "1",
            "junction_fraction_similarity": "1",
            "coding_span_ratio": "0.8",
        },
        {
            "candidate_digest": "D2",
            "topology_available": "0",
            "cds_bp_ratio": "",
            "cds_segment_count_ratio": "",
            "phase_lcs_similarity": "",
            "junction_fraction_similarity": "",
            "coding_span_ratio": "",
        },
    ]
    _write_tsv(topology_path, topology_fields, topology_rows)
    topology_manifest = {
        "schema_version": HOMEOLOG_TOPOLOGY_SCHEMA_VERSION,
        "truth_access": False,
        "inputs": {"copy_features": _sha256(copy_path)},
        "outputs": {
            "features": {"sha256": _sha256(topology_path), "rows": 2}
        },
    }
    Path(str(topology_path) + ".manifest.json").write_text(
        json.dumps(topology_manifest), encoding="utf-8"
    )

    contract = fit_copy_feature_contract(copy_rows, feature_set="full")
    base_names = contract["expanded_feature_names"]
    topology_coefficients = [0.0] * (len(base_names) + len(TOPOLOGY_ADDON_FIELDS))
    topology_coefficients[-1] = 2.0
    model = {
        "schema_version": HOMEOLOG_RANKER_SCHEMA_VERSION,
        "base_feature_contract": contract,
        "topology_addon_fields": list(TOPOLOGY_ADDON_FIELDS),
        "estimators": {
            "baseline": {
                "intercept": 0.0,
                "coefficients": [0.0] * len(base_names),
                "coefficient_feature_order": base_names,
            },
            "topology": {
                "intercept": 0.0,
                "coefficients": topology_coefficients,
                "coefficient_feature_order": [
                    *base_names,
                    *TOPOLOGY_ADDON_FIELDS,
                ],
            },
        },
        "claim_boundary": {"automatic_approval": False},
    }
    model_path = tmp_path / "model.json"
    model_path.write_text(json.dumps(model), encoding="utf-8")
    return copy_path, topology_path, model_path


def test_score_homeolog_copy_candidates_is_rank_only_and_truth_blind(
    tmp_path: Path,
) -> None:
    copy_path, topology_path, model_path = _fixture(tmp_path)
    output = tmp_path / "scores.tsv"
    manifest = score_homeolog_copy_candidates(
        copy_feature_tsv_path=copy_path,
        topology_tsv_path=topology_path,
        model_json_path=model_path,
        output_tsv_path=output,
    )

    assert manifest["truth_access"] is False
    assert manifest["counts"] == {
        "candidates": 2,
        "topology_available": 1,
        "topology_unavailable": 1,
        "automatic_approved": 0,
    }
    with output.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert float(rows[0]["homeolog_topology_rank_score"]) > float(
        rows[1]["homeolog_topology_rank_score"]
    )
    assert rows[0]["homeolog_topology_rank_percentile"] == "1"
    assert rows[1]["homeolog_topology_rank_percentile"] == "0"
    assert {row["homeolog_automatic_approval"] for row in rows} == {"0"}
    assert rows[1]["homeolog_evidence_tier"].endswith("topology_unavailable")


def test_homeolog_ranker_cli_routes(tmp_path: Path) -> None:
    copy_path, topology_path, model_path = _fixture(tmp_path)
    output = tmp_path / "scores.tsv"
    assert main(
        [
            "evidence",
            "score-homeolog-copy-candidates",
            "--copy-features",
            str(copy_path),
            "--topology-features",
            str(topology_path),
            "--model-json",
            str(model_path),
            "--output-tsv",
            str(output),
        ]
    ) == 0
    assert output.exists()


def test_freeze_homeolog_review_rankings_is_deterministic_and_truth_blind(
    tmp_path: Path,
) -> None:
    copy_path, topology_path, model_path = _fixture(tmp_path)
    scores = tmp_path / "scores.tsv"
    score_homeolog_copy_candidates(
        copy_feature_tsv_path=copy_path,
        topology_tsv_path=topology_path,
        model_json_path=model_path,
        output_tsv_path=scores,
    )
    rankings = tmp_path / "review.tsv"
    manifest = freeze_homeolog_review_rankings(
        score_tsv_path=scores,
        output_tsv_path=rankings,
        review_budgets=(1, 2, 2),
    )

    assert manifest["schema_version"] == HOMEOLOG_REVIEW_RANKING_SCHEMA_VERSION
    assert manifest["truth_access"] is False
    assert manifest["policy"]["review_budgets"] == [1, 2]
    assert manifest["counts"]["candidates"] == 2
    assert manifest["counts"]["ranking_rows"] == 4
    with rankings.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    baseline = [row for row in rows if row["estimator"] == "baseline"]
    topology = [row for row in rows if row["estimator"] == "topology"]
    assert [row["candidate_digest"] for row in baseline] == ["D1", "D2"]
    assert [row["candidate_digest"] for row in topology] == ["D1", "D2"]
    assert baseline[0]["within_top_1"] == "1"
    assert baseline[1]["within_top_1"] == "0"
    assert {row["automatic_approval"] for row in rows} == {"0"}


def test_homeolog_review_ranking_cli_routes(tmp_path: Path) -> None:
    copy_path, topology_path, model_path = _fixture(tmp_path)
    scores = tmp_path / "scores.tsv"
    score_homeolog_copy_candidates(
        copy_feature_tsv_path=copy_path,
        topology_tsv_path=topology_path,
        model_json_path=model_path,
        output_tsv_path=scores,
    )
    rankings = tmp_path / "review.tsv"
    assert main(
        [
            "evidence",
            "freeze-homeolog-review-rankings",
            "--scores",
            str(scores),
            "--review-budget",
            "1",
            "--review-budget",
            "2",
            "--output-tsv",
            str(rankings),
        ]
    ) == 0
    assert rankings.exists()
