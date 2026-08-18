from __future__ import annotations

import json
from pathlib import Path

from ploidypatch.structure_bootstrap import (
    paired_structure_hypothesis_bootstrap,
)
from ploidypatch.structure_hypothesis_score import (
    STRUCTURE_HYPOTHESIS_SCORE_SCHEMA_VERSION,
)


def score(path: Path, event_type: str, values: tuple[bool, ...]) -> Path:
    report = {
        "schema_version": STRUCTURE_HYPOTHESIS_SCORE_SCHEMA_VERSION,
        "quality_gate": {"grade": "pass"},
        "event_details": [
            {
                "event_id": f"{event_type}-{index}",
                "event_type": event_type,
                "recovered": value,
            }
            for index, value in enumerate(values)
        ],
    }
    path.write_text(json.dumps(report), encoding="utf-8")
    return path


def test_structure_bootstrap_merges_event_files_per_tier(tmp_path: Path) -> None:
    inputs = [
        ("tier_b", score(tmp_path / "b1.json", "boundary", (True, False))),
        ("tier_b", score(tmp_path / "b2.json", "split", (True, True))),
        ("tier_a", score(tmp_path / "a1.json", "boundary", (False, False))),
        ("tier_a", score(tmp_path / "a2.json", "split", (True, False))),
    ]
    report = paired_structure_hypothesis_bootstrap(
        score_inputs=inputs,
        output_json_path=tmp_path / "bootstrap.json",
        replicates=200,
        seed=17,
    )

    methods = {row["label"]: row for row in report["methods"]}
    assert methods["tier_b"]["successes"] == 3
    assert methods["tier_a"]["successes"] == 1
    assert report["counts"]["event_types"] == {"boundary": 2, "split": 2}
    assert report["paired_differences"][0]["observed_delta"] == 0.5
