from __future__ import annotations

import importlib.util
import csv
import json
from pathlib import Path
import sys

import numpy as np
import pytest


MODULE_PATH = Path(__file__).parents[1] / "scripts" / "evaluate_external_v0.5.py"
SPEC = importlib.util.spec_from_file_location("external_v05", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def policy() -> dict[str, str]:
    return {
        "minimum_formal_event_count": "500",
        "minimum_target_chromosomes": "22",
        "minimum_events_per_complexity_bin": "20",
    }


def evaluability(events: int = 500) -> dict[str, object]:
    return {
        "events": events,
        "target_chromosomes": 29,
        "complexity_bins": {
            "one": 125,
            "two_to_three": 125,
            "four_to_six": 125,
            "seven_plus": 125,
        },
        "sentinels": {
            "blind_noop_exact_recovery": 0,
            "complete_oracle_exact_recovery": events,
            "restoration_byte_identical": True,
            "blind_complete_genome_sha256_identical": True,
        },
    }


def test_review_metrics_freeze_six_budgets_ties_and_digest() -> None:
    labels = np.asarray([1] + [0] * 299, dtype=np.uint8)
    scores = np.zeros(300, dtype=float)
    digests = [f"d{index:03d}" for index in reversed(range(300))]
    review = MODULE.review_metrics(labels, scores, digests)
    assert list(review) == [
        "top_0.5pct",
        "top_1pct",
        "top_2pct",
        "top_100",
        "top_250",
        "top_500",
    ]
    expected = [f"d{index:03d}" for index in range(250)]
    assert review["top_250"]["reviewed"] == 250
    assert review["top_250"]["selected_digest_sha256"] == (
        MODULE.selected_digest_sha256(expected)
    )


def test_evaluability_distinguishes_invalid_sentinel_and_insufficient_data() -> None:
    formal, gates = MODULE.validate_evaluability(evaluability(), policy())
    assert formal is True
    assert all(gates.values())

    insufficient = evaluability(events=499)
    insufficient["sentinels"]["complete_oracle_exact_recovery"] = 499  # type: ignore[index]
    formal, gates = MODULE.validate_evaluability(insufficient, policy())
    assert formal is False
    assert gates["minimum_events"] is False

    invalid = evaluability()
    invalid["sentinels"]["restoration_byte_identical"] = False  # type: ignore[index]
    with pytest.raises(ValueError, match="invalidates run"):
        MODULE.validate_evaluability(invalid, policy())


def test_custody_rejects_any_truth_or_network_access(tmp_path: Path) -> None:
    scores = tmp_path / "scores.tsv"
    score_manifest = tmp_path / "scores.tsv.manifest.json"
    scores.write_text("candidate_digest\na\n", encoding="utf-8")
    score_manifest.write_text("{}\n", encoding="utf-8")
    custody = {
        "schema_version": MODULE.CUSTODY_SCHEMA,
        "runner_identity": "isolated_blind_runner",
        "truth_mounted": False,
        "complete_target_annotation_mounted": False,
        "evaluator_references_mounted": False,
        "nas_data_mounted": False,
        "network_access": False,
        "frozen_before_truth_reveal_at": "2026-08-08T00:00:00Z",
        "blind_outputs": {
            "scores": {"relative_path": "scores.tsv", "sha256": MODULE.sha256(scores)},
            "score_manifest": {
                "relative_path": "scores.tsv.manifest.json",
                "sha256": MODULE.sha256(score_manifest),
            },
        },
    }
    MODULE.validate_custody(custody, scores, score_manifest)
    custody["network_access"] = True
    with pytest.raises(ValueError, match="custody"):
        MODULE.validate_custody(custody, scores, score_manifest)


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def freeze_dir(root: Path, names: list[str]) -> None:
    with (root / "SHA256SUMS").open("x", encoding="utf-8", newline="") as handle:
        for name in names:
            handle.write(f"{MODULE.sha256(root / name)}  {name}\n")


def test_toy_reveal_writes_audited_not_evaluable_result(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(MODULE, "BOOTSTRAP_REPLICATES", 100)
    composite = tmp_path / "composite"
    composite.mkdir()
    (composite / "composite_manifest.json").write_text(
        json.dumps(
            {
                "schema_version": "ploidypatch.composite_ranker.v0.4",
                "automatic_approval": False,
            }
        ),
        encoding="utf-8",
    )
    freeze_dir(composite, ["composite_manifest.json"])

    protocol = tmp_path / "protocol"
    protocol.mkdir()
    policy_path = protocol / "policy.tsv"
    policy_rows = {
        "policy_id": MODULE.POLICY_ID,
        "test_role": "target_level_predeclared_untouched_secondary_replication",
        "model_version": "PloidyPatch_ranker_v0.4",
        "automatic_copy_addition_approval": "false",
        "H1_bootstrap_seed": str(MODULE.H1_BOOTSTRAP_SEED),
        "H2_bootstrap_seed": str(MODULE.H2_BOOTSTRAP_SEED),
        "guard_v03_bootstrap_seed": str(MODULE.GUARD_V03_BOOTSTRAP_SEED),
        "bootstrap_replicates": str(MODULE.BOOTSTRAP_REPLICATES),
        "minimum_chromosome_bootstrap_valid_replicates": str(
            MODULE.MINIMUM_VALID_BOOTSTRAPS
        ),
        "minimum_formal_event_count": "500",
        "minimum_target_chromosomes": "22",
        "minimum_events_per_complexity_bin": "20",
        "minimum_topology_coverage_among_positive_candidates": "0.70",
        "minimum_v03_AP_gain_retained_fraction": "0.90",
    }
    write_tsv(
        policy_path,
        [{"field": key, "value": value} for key, value in policy_rows.items()],
    )
    (protocol / "protocol_manifest.json").write_text(
        json.dumps(
            {
                "schema_version": "ploidypatch.external_holdout_protocol_freeze.v0.5",
                "policy_id": MODULE.POLICY_ID,
                "model_version": "PloidyPatch_ranker_v0.4",
                "composite_model_SHA256SUMS_sha256": MODULE.sha256(
                    composite / "SHA256SUMS"
                ),
            }
        ),
        encoding="utf-8",
    )
    freeze_dir(protocol, ["policy.tsv", "protocol_manifest.json"])

    pool_decisions = tmp_path / "decisions.tsv"
    decision_rows = [
        {
            "consensus_digest": digest,
            "conflict_set_digest": f"x{index // 2}",
            "status": "accepted",
        }
        for index, digest in enumerate(("a", "b", "c", "d"))
    ]
    write_tsv(pool_decisions, decision_rows)
    pool_manifest = tmp_path / "pool.manifest.json"
    pool_manifest.write_text(
        json.dumps(
            {
                "schema_version": MODULE.POOL_SCHEMA,
                "outputs": {
                    "decisions": {"sha256": MODULE.sha256(pool_decisions)}
                },
            }
        ),
        encoding="utf-8",
    )

    scores = tmp_path / "scores.tsv"
    score_rows = []
    for index, digest in enumerate(("a", "b", "c", "d")):
        score_rows.append(
            {
                "candidate_digest": digest,
                "seqid": str(index // 2 + 1),
                "support_method_count": 2,
                "v03_baseline_logit": (4, 3, 2, 1)[index],
                "v03_primary_rank_score": (4, 3, 2, 1)[index],
                "v03_topology_available": 1,
                "v04_primary_rank_score": (4, 3, 2, 1)[index],
                "v04_conflict_guard_applied": 0,
                "v04_topology_abstained": 0,
                "v04_automatic_approval": 0,
            }
        )
    write_tsv(scores, score_rows)
    score_manifest = Path(str(scores) + ".manifest.json")
    score_manifest.write_text(
        json.dumps(
            {
                "schema_version": MODULE.GUARD_SCORE_SCHEMA,
                "truth_access": False,
                "inputs": {
                    "pool_decisions": MODULE.sha256(pool_decisions),
                    "pool_manifest": MODULE.sha256(pool_manifest),
                },
                "outputs": {"scores": {"sha256": MODULE.sha256(scores)}},
                "winner_audit": {
                    "mismatch_count": 0,
                    "baseline_mapping_sha256": "same",
                    "v04_guard_mapping_sha256": "same",
                },
            }
        ),
        encoding="utf-8",
    )

    labels = tmp_path / "labels.tsv"
    write_tsv(
        labels,
        [
            {"candidate_digest": digest, "label_exact_cds": label}
            for digest, label in zip(("a", "b", "c", "d"), (1, 0, 1, 0))
        ],
    )
    Path(str(labels) + ".manifest.json").write_text(
        json.dumps(
            {
                "schema_version": MODULE.LABEL_SCHEMA,
                "evaluator_only": True,
                "blind_scores_sha256": MODULE.sha256(scores),
                "pool_manifest_sha256": MODULE.sha256(pool_manifest),
                "outputs": {"labels": {"sha256": MODULE.sha256(labels)}},
            }
        ),
        encoding="utf-8",
    )

    events = [
        {"event_id": f"e{index}", "complete_cds_chain_recovery": value}
        for index, value in enumerate((1, 1, 1, 0))
    ]
    pool_score_common = {
        "quality_gate": {"grade": "pass"},
        "collateral_changes": {
            "baseline_transcript_structures_missing_from_candidate": 0
        },
        "event_details": events,
    }
    primary_score = tmp_path / "primary.json"
    primary_score.write_text(json.dumps(pool_score_common), encoding="utf-8")
    legacy_score = tmp_path / "legacy.json"
    legacy_score.write_text(
        json.dumps(
            {
                **pool_score_common,
                "event_details": [
                    {
                        "event_id": f"e{index}",
                        "complete_cds_chain_recovery": value,
                    }
                    for index, value in enumerate((1, 1, 0, 0))
                ],
            }
        ),
        encoding="utf-8",
    )
    evaluability_path = tmp_path / "evaluability.json"
    evaluability_path.write_text(json.dumps(evaluability(events=4)), encoding="utf-8")
    custody_path = tmp_path / "custody.json"
    custody_path.write_text(
        json.dumps(
            {
                "schema_version": MODULE.CUSTODY_SCHEMA,
                "runner_identity": "isolated_blind_runner",
                "truth_mounted": False,
                "complete_target_annotation_mounted": False,
                "evaluator_references_mounted": False,
                "nas_data_mounted": False,
                "network_access": False,
                "frozen_before_truth_reveal_at": "2026-08-08T00:00:00Z",
                "blind_outputs": {
                    "scores": {
                        "relative_path": "scores.tsv",
                        "sha256": MODULE.sha256(scores),
                    },
                    "score_manifest": {
                        "relative_path": "scores.tsv.manifest.json",
                        "sha256": MODULE.sha256(score_manifest),
                    },
                },
            }
        ),
        encoding="utf-8",
    )
    output = tmp_path / "evaluation"
    monkeypatch.setattr(
        sys,
        "argv",
        [
            str(MODULE_PATH),
            "--scores",
            str(scores),
            "--labels",
            str(labels),
            "--pool-decisions",
            str(pool_decisions),
            "--pool-manifest",
            str(pool_manifest),
            "--primary-pool-score",
            str(primary_score),
            "--legacy-pool-score",
            str(legacy_score),
            "--evaluability",
            str(evaluability_path),
            "--custody-manifest",
            str(custody_path),
            "--protocol-freeze",
            str(protocol),
            "--composite-model-freeze",
            str(composite),
            "--policy",
            str(policy_path),
            "--input-root",
            f"fixture={tmp_path}",
            "--output-dir",
            str(output),
        ],
    )
    assert MODULE.main() == 0
    report = json.loads((output / "evaluation.json").read_text(encoding="utf-8"))
    assert report["formal_outcome"] == "not_evaluable_without_rule_relaxation"
    assert report["confirmatory_pass"] is False
    assert report["winner_audit"]["mismatch_count"] == 0
    assert all(item["root_role"] == "fixture" for item in report["inputs"].values())
    assert all(item["path"].startswith("@fixture/") for item in report["inputs"].values())
    assert str(tmp_path) not in json.dumps(report["inputs"], sort_keys=True)
    assert (output / "SHA256SUMS").is_file()
