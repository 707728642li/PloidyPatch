from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from ploidypatch.copy_features import COPY_FEATURE_SCHEMA_VERSION, FEATURE_COLUMNS
from ploidypatch.copy_model import (
    COPY_MODEL_SCHEMA_VERSION,
    fit_copy_feature_contract,
    predict_copy_candidate,
    score_copy_candidate_features,
    select_scored_copy_candidates,
)


def _row(digest: str, *, wgd: str) -> dict[str, str]:
    row = {field: "" for field in FEATURE_COLUMNS}
    row.update(
        {
            "candidate_digest": digest,
            "seqid": "chr1",
            "start": "10",
            "end": "30",
            "strand": "+",
            "span_bp": "21",
            "cds_segments": "1",
            "cds_bp": "21",
            "support_method_count": "1",
            "support_methods": "miniprot",
            "has_miniprot": "1",
            "has_gemoma": "0",
            "has_lifton": "0",
            "miniprot_identity": "0.9",
            "miniprot_query_coverage": "0.8",
            "miniprot_rank": "1",
            "miniprot_score": "100",
            "miniprot_positive": "1",
            "miniprot_frameshifts": "0",
            "miniprot_stop_codons": "0",
            "wgd_existing_partner": wgd,
            "wgd_support_block_count": "3" if wgd == "1" else "0",
            "wgd_longest_block_pairs": "12" if wgd == "1" else "0",
        }
    )
    return row


def _model(rows: list[dict[str, str]]) -> dict[str, object]:
    contract = fit_copy_feature_contract(rows)
    coefficients = [0.0] * len(contract["expanded_feature_names"])
    coefficients[contract["expanded_feature_names"].index("binary:wgd_existing_partner")] = 4.0
    return {
        "schema_version": COPY_MODEL_SCHEMA_VERSION,
        "feature_contract": contract,
        "estimator": {
            "family": "logistic_regression_l2",
            "intercept": -2.0,
            "coefficients": coefficients,
        },
        "calibration": {"method": "sigmoid_on_raw_logit", "slope": 1.0, "intercept": 0.0},
        "thresholds": {
            "review": {"value": 0.5},
            "high_confidence": {"value": 0.8},
        },
    }


def _write_features(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FEATURE_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    import hashlib

    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    Path(str(path) + ".manifest.json").write_text(
        json.dumps(
            {
                "schema_version": COPY_FEATURE_SCHEMA_VERSION,
                "truth_access": False,
                "outputs": {"features": {"sha256": digest}},
            }
        ),
        encoding="utf-8",
    )


def test_copy_model_prediction_and_truth_blind_scoring(tmp_path: Path) -> None:
    rows = [_row("a" * 64, wgd="0"), _row("b" * 64, wgd="1")]
    model = _model(rows)
    assert predict_copy_candidate(rows[0], model)["calibrated_probability"] < 0.5
    assert predict_copy_candidate(rows[1], model)["calibrated_probability"] > 0.8

    features = tmp_path / "features.tsv"
    model_json = tmp_path / "model.json"
    scores = tmp_path / "scores.tsv"
    _write_features(features, rows)
    model_json.write_text(json.dumps(model), encoding="utf-8")
    manifest = score_copy_candidate_features(
        feature_tsv_path=features,
        model_json_path=model_json,
        output_tsv_path=scores,
    )
    assert manifest["truth_access"] is False
    assert manifest["counts"] == {
        "candidates": 2,
        "review_selected": 1,
        "high_confidence_selected": 1,
    }
    with scores.open(encoding="utf-8", newline="") as handle:
        scored = list(csv.DictReader(handle, delimiter="\t"))
    assert [row["model_review_decision"] for row in scored] == ["0", "1"]


def test_scored_copy_selection_preserves_prefix_and_whole_hierarchy(tmp_path: Path) -> None:
    digest_a = "a" * 64
    digest_b = "b" * 64
    rows = [_row(digest_a, wgd="0"), _row(digest_b, wgd="1")]
    model = _model(rows)
    base = tmp_path / "base.gff3"
    candidate = tmp_path / "candidate.gff3"
    features = tmp_path / "features.tsv"
    model_json = tmp_path / "model.json"
    scores = tmp_path / "scores.tsv"
    output = tmp_path / "selected.gff3"
    decisions = tmp_path / "selection.tsv"
    base_text = "##gff-version 3\nchr1\ttest\tgene\t1\t5\t.\t+\t.\tID=existing\n"
    base.write_text(base_text, encoding="utf-8")
    appended = []
    for index, digest in enumerate((digest_a, digest_b), start=1):
        appended.extend(
            (
                f"chr1\tPloidyPatchConsensus\tgene\t{index * 10}\t{index * 10 + 2}\t.\t+\t.\tID=g{index};consensus_digest={digest}\n",
                f"chr1\tPloidyPatchConsensus\tmRNA\t{index * 10}\t{index * 10 + 2}\t.\t+\t.\tID=t{index};Parent=g{index};consensus_digest={digest}\n",
                f"chr1\tPloidyPatchConsensus\texon\t{index * 10}\t{index * 10 + 2}\t.\t+\t.\tID=e{index};Parent=t{index}\n",
                f"chr1\tPloidyPatchConsensus\tCDS\t{index * 10}\t{index * 10 + 2}\t.\t+\t0\tParent=t{index}\n",
            )
        )
    candidate.write_text(base_text + "###\n" + "".join(appended), encoding="utf-8")
    _write_features(features, rows)
    model_json.write_text(json.dumps(model), encoding="utf-8")
    score_copy_candidate_features(
        feature_tsv_path=features,
        model_json_path=model_json,
        output_tsv_path=scores,
    )
    manifest = select_scored_copy_candidates(
        base_gff_path=base,
        candidate_gff_path=candidate,
        scored_tsv_path=scores,
        model_json_path=model_json,
        output_gff_path=output,
        selection_tsv_path=decisions,
    )
    observed = output.read_text(encoding="utf-8")
    assert observed.startswith(base_text + "###\n")
    assert digest_a not in observed
    assert digest_b in observed
    assert "ID=g2" in observed and "ID=t2" in observed and "ID=e2" in observed
    assert manifest["counts"]["selected_models"] == 1


def test_scored_copy_selection_rejects_digest_universe_mismatch(tmp_path: Path) -> None:
    rows = [_row("a" * 64, wgd="0"), _row("b" * 64, wgd="1")]
    model = _model(rows)
    base = tmp_path / "base.gff3"
    candidate = tmp_path / "candidate.gff3"
    features = tmp_path / "features.tsv"
    model_json = tmp_path / "model.json"
    scores = tmp_path / "scores.tsv"
    base.write_text("##gff-version 3\n", encoding="utf-8")
    candidate.write_text(
        base.read_text(encoding="utf-8")
        + "###\n"
        + f"chr1\tPloidyPatchConsensus\tgene\t1\t3\t.\t+\t.\tID=g;consensus_digest={'a' * 64}\n"
        + f"chr1\tPloidyPatchConsensus\tmRNA\t1\t3\t.\t+\t.\tID=t;Parent=g;consensus_digest={'a' * 64}\n"
        + "chr1\tPloidyPatchConsensus\tCDS\t1\t3\t.\t+\t0\tParent=t\n",
        encoding="utf-8",
    )
    _write_features(features, rows)
    model_json.write_text(json.dumps(model), encoding="utf-8")
    score_copy_candidate_features(
        feature_tsv_path=features,
        model_json_path=model_json,
        output_tsv_path=scores,
    )
    with pytest.raises(ValueError, match="digest universes differ"):
        select_scored_copy_candidates(
            base_gff_path=base,
            candidate_gff_path=candidate,
            scored_tsv_path=scores,
            model_json_path=model_json,
            output_gff_path=tmp_path / "out.gff3",
            selection_tsv_path=tmp_path / "decisions.tsv",
        )


def test_counterfactual_feature_masks_keep_candidate_universe_and_are_audited(
    tmp_path: Path,
) -> None:
    features = tmp_path / "features.tsv"
    rows = [_row("a" * 64, wgd="0"), _row("b" * 64, wgd="1")]
    _write_features(features, rows)
    model_path = tmp_path / "model.json"
    model = _model(rows)
    model_path.write_text(json.dumps(model), encoding="utf-8")

    unmasked = tmp_path / "unmasked.tsv"
    masked = tmp_path / "masked.tsv"
    score_copy_candidate_features(
        feature_tsv_path=features,
        model_json_path=model_path,
        output_tsv_path=unmasked,
    )
    manifest = score_copy_candidate_features(
        feature_tsv_path=features,
        model_json_path=model_path,
        output_tsv_path=masked,
        mask_feature_groups=("wgd_context", "method_quality"),
    )

    with unmasked.open("r", encoding="utf-8", newline="") as handle:
        unmasked_rows = list(csv.DictReader(handle, delimiter="\t"))
    with masked.open("r", encoding="utf-8", newline="") as handle:
        masked_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [row["candidate_digest"] for row in masked_rows] == [
        row["candidate_digest"] for row in unmasked_rows
    ]
    assert manifest["counterfactual_feature_masks"] == [
        "wgd_context",
        "method_quality",
    ]
    assert manifest["truth_access"] is False

    with pytest.raises(ValueError, match="must be unique"):
        score_copy_candidate_features(
            feature_tsv_path=features,
            model_json_path=model_path,
            output_tsv_path=tmp_path / "duplicate.tsv",
            mask_feature_groups=("wgd_context", "wgd_context"),
        )
