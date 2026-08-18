from __future__ import annotations

import csv
import importlib.util
import json
import math
from pathlib import Path

import pytest


MODULE_PATH = (
    Path(__file__).parents[1]
    / "scripts"
    / "verify_conflict_guard_production_replay_v0.4.py"
)
SPEC = importlib.util.spec_from_file_location("guard_replay", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def toy_rows() -> tuple[dict[str, str], dict[str, str]]:
    development = {
        "dataset": "toy",
        "candidate_digest": "a",
        "baseline": "1.0",
        "primary": "2.0",
        "v04_guard": "1.0",
        "v04_guard_applied": "1",
        "v04_topology_abstained": "1",
        "v04_automatic_approval": "0",
    }
    production = {
        "candidate_digest": "a",
        "v03_baseline_logit": "1",
        "v03_primary_rank_score": "2",
        "v04_primary_rank_score": "1",
        "v04_conflict_guard_applied": "1",
        "v04_topology_abstained": "1",
        "v04_automatic_approval": "0",
    }
    return development, production


def write_manifest(path: Path) -> None:
    manifest = {
        "schema_version": MODULE.GUARD_MANIFEST_SCHEMA,
        "truth_access": False,
        "outputs": {"scores": {"sha256": MODULE.sha256(path)}},
        "winner_audit": {
            "mismatch_count": 0,
            "baseline_mapping_sha256": "same",
            "v04_guard_mapping_sha256": "same",
        },
    }
    Path(str(path) + ".manifest.json").write_text(
        json.dumps(manifest), encoding="utf-8"
    )


def test_verify_dataset_requires_exact_production_replay(tmp_path: Path) -> None:
    development, production = toy_rows()
    output = tmp_path / "v04.tsv"
    write_tsv(output, [production])
    write_manifest(output)
    report = MODULE.verify_dataset("toy", output, [development])
    assert report["mismatch_count"] == 0
    assert report["candidates"] == 1
    assert report["score_tolerance"]["bitwise_difference_count"] == 0


def test_verify_dataset_rejects_score_drift(tmp_path: Path) -> None:
    development, production = toy_rows()
    production["v04_primary_rank_score"] = "0.999"
    output = tmp_path / "v04.tsv"
    write_tsv(output, [production])
    write_manifest(output)
    with pytest.raises(AssertionError, match="differs"):
        MODULE.verify_dataset("toy", output, [development])


def test_verify_dataset_rejects_truth_column(tmp_path: Path) -> None:
    development, production = toy_rows()
    production["label_exact_cds"] = "1"
    output = tmp_path / "v04.tsv"
    write_tsv(output, [production])
    write_manifest(output)
    with pytest.raises(ValueError, match="truth"):
        MODULE.verify_dataset("toy", output, [development])


def test_verify_dataset_allows_reported_single_ulp_replay_noise(
    tmp_path: Path,
) -> None:
    development, production = toy_rows()
    production["v04_primary_rank_score"] = repr(math.nextafter(1.0, 2.0))
    output = tmp_path / "v04.tsv"
    write_tsv(output, [production])
    write_manifest(output)
    report = MODULE.verify_dataset("toy", output, [development])
    assert report["mismatch_count"] == 0
    assert report["score_tolerance"]["bitwise_difference_count"] == 1
    assert report["score_tolerance"]["maximum_observed_ulp_difference"] == 1.0
