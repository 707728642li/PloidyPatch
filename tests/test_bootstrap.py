from __future__ import annotations

import json
from pathlib import Path

import pytest

from ploidypatch.bootstrap import paired_event_bootstrap


def write_score(path: Path, values: list[tuple[str, str, bool]]) -> Path:
    report = {
        "schema_version": "ploidypatch.annotation_repair_score.v5",
        "quality_gate": {"grade": "pass"},
        "event_details": [
            {
                "event_id": event_id,
                "event_type": event_type,
                "complete_cds_chain_recovery": value,
            }
            for event_id, event_type, value in values
        ],
    }
    path.write_text(json.dumps(report), encoding="utf-8", newline="")
    return path


def test_paired_event_bootstrap_is_deterministic_and_stratified(
    tmp_path: Path,
) -> None:
    left = write_score(
        tmp_path / "left.json",
        [
            ("E1", "missing_exon", True),
            ("E2", "missing_exon", True),
            ("E3", "split", True),
            ("E4", "split", False),
        ],
    )
    right = write_score(
        tmp_path / "right.json",
        [
            ("E1", "missing_exon", False),
            ("E2", "missing_exon", True),
            ("E3", "split", False),
            ("E4", "split", False),
        ],
    )
    output = tmp_path / "bootstrap.json"
    report = paired_event_bootstrap(
        score_inputs=[("left", left), ("right", right)],
        output_json_path=output,
        replicates=1000,
        seed=17,
    )
    assert report["counts"]["event_types"] == {"missing_exon": 2, "split": 2}
    assert report["methods"][0]["observed_rate"] == 0.75
    assert report["methods"][1]["observed_rate"] == 0.25
    assert report["paired_differences"][0]["observed_delta"] == 0.5
    assert report["paired_differences"][0]["probability_delta_gt_zero"] > 0.8
    with pytest.raises(FileExistsError, match="Refusing to overwrite"):
        paired_event_bootstrap(
            score_inputs=[("left", left), ("right", right)],
            output_json_path=output,
            replicates=100,
        )


def test_paired_event_bootstrap_rejects_unpaired_events(tmp_path: Path) -> None:
    left = write_score(tmp_path / "left.json", [("E1", "split", True)])
    right = write_score(tmp_path / "right.json", [("E2", "split", True)])
    with pytest.raises(ValueError, match="Paired event IDs differ"):
        paired_event_bootstrap(
            score_inputs=[("left", left), ("right", right)],
            output_json_path=tmp_path / "out.json",
            replicates=100,
        )
