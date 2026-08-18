from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import pytest

from ploidypatch.cli import main as cli_main
from ploidypatch.conflict_guard import (
    apply_conflict_winner_guard,
    compute_conflict_winner_guard,
)
from ploidypatch.support_ranker import SUPPORT_CONDITIONED_SCORE_SCHEMA_VERSION


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _write_inputs(tmp_path: Path) -> tuple[Path, Path, Path]:
    scores = tmp_path / "v03.tsv"
    rows = [
        ("a", 4.0, 1.0),
        ("b", 3.0, 5.0),
        ("c", 2.0, 6.0),
        ("d", 1.0, 3.0),
        ("e", 0.0, 2.0),
    ]
    with scores.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("candidate_digest", "v03_baseline_logit", "v03_primary_rank_score"))
        writer.writerows(rows)
    Path(str(scores) + ".manifest.json").write_text(
        json.dumps(
            {
                "schema_version": SUPPORT_CONDITIONED_SCORE_SCHEMA_VERSION,
                "truth_access": False,
                "outputs": {"scores": {"sha256": _sha256(scores)}},
            }
        ),
        encoding="utf-8",
    )
    decisions = tmp_path / "decisions.tsv"
    with decisions.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            (
                "consensus_digest",
                "conflict_set_digest",
                "conflict_member_count",
                "status",
            )
        )
        writer.writerows(
            (
                ("a", "x", 2, "accepted"),
                ("b", "x", 2, "accepted"),
                ("c", "y", 2, "accepted"),
                ("d", "y", 2, "accepted"),
                ("e", "", 1, "accepted"),
                ("rejected", "", 1, "rejected"),
            )
        )
    pool_manifest = tmp_path / "candidate.gff3.manifest.json"
    pool_manifest.write_text(
        json.dumps(
            {
                "schema_version": "ploidypatch.method_candidate_pool.v2",
                "counts": {
                    "accepted_models": 5,
                    "conflict_sets": 2,
                    "conflicted_chains": 4,
                },
                "outputs": {
                    "decisions": {
                        "rows": 6,
                        "sha256": _sha256(decisions),
                    }
                },
                "parameters": {
                    "conflict_action": "retain_all_for_ranking_and_review",
                    "max_redundancy_overlap": 0.5,
                    "min_method_support": 1,
                    "redundancy_policy": "retain_distinct_chains",
                },
            }
        ),
        encoding="utf-8",
    )
    return scores, decisions, pool_manifest


def test_guard_falls_back_only_when_conflict_winner_changes(tmp_path: Path) -> None:
    scores, decisions, pool_manifest = _write_inputs(tmp_path)
    output = tmp_path / "v04.tsv"
    manifest = apply_conflict_winner_guard(
        v03_score_tsv_path=scores,
        pool_decisions_tsv_path=decisions,
        pool_manifest_json_path=pool_manifest,
        output_tsv_path=output,
    )
    with output.open(encoding="utf-8", newline="") as handle:
        rows = {row["candidate_digest"]: row for row in csv.DictReader(handle, delimiter="\t")}
    assert float(rows["a"]["v04_primary_rank_score"]) == 4.0
    assert float(rows["b"]["v04_primary_rank_score"]) == 3.0
    assert rows["a"]["v04_conflict_guard_applied"] == "1"
    assert float(rows["c"]["v04_primary_rank_score"]) == 6.0
    assert float(rows["d"]["v04_primary_rank_score"]) == 3.0
    assert rows["c"]["v04_conflict_guard_applied"] == "0"
    assert float(rows["e"]["v04_primary_rank_score"]) == 2.0
    assert all(row["v04_automatic_approval"] == "0" for row in rows.values())
    assert manifest["truth_access"] is False
    assert manifest["counts"]["winner_disagreement_sets"] == 1
    assert manifest["counts"]["guarded_candidates"] == 2
    assert manifest["winner_audit"]["mismatch_count"] == 0


def test_guard_rejects_nonblind_score_manifest(tmp_path: Path) -> None:
    scores, decisions, pool_manifest = _write_inputs(tmp_path)
    manifest_path = Path(str(scores) + ".manifest.json")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["truth_access"] = True
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    with pytest.raises(ValueError, match="schema, truth, or checksum"):
        apply_conflict_winner_guard(
            v03_score_tsv_path=scores,
            pool_decisions_tsv_path=decisions,
            pool_manifest_json_path=pool_manifest,
            output_tsv_path=tmp_path / "v04.tsv",
        )


def test_guard_uses_digest_tie_break_and_reports_exact_winner_mapping() -> None:
    result = compute_conflict_winner_guard(
        digests=("b", "a", "c"),
        conflict_by_digest={"a": "x", "b": "x", "c": ""},
        baseline_scores={"a": 1.0, "b": 1.0, "c": 0.0},
        primary_scores={"a": 0.0, "b": 2.0, "c": 3.0},
    )
    assert result["baseline_winners"] == {"x": "a"}
    assert result["primary_winners"] == {"x": "b"}
    assert result["guard_winners"] == {"x": "a"}
    assert result["scores"]["a"] == result["scores"]["b"] == 1.0
    assert result["winner_mismatch_count"] == 0
    assert (
        result["baseline_winner_mapping_sha256"]
        == result["guard_winner_mapping_sha256"]
    )


def test_guard_rejects_pool_manifest_checksum_mismatch(tmp_path: Path) -> None:
    scores, decisions, pool_manifest = _write_inputs(tmp_path)
    manifest = json.loads(pool_manifest.read_text(encoding="utf-8"))
    manifest["outputs"]["decisions"]["sha256"] = "0" * 64
    pool_manifest.write_text(json.dumps(manifest), encoding="utf-8")
    with pytest.raises(ValueError, match="schema or decisions checksum"):
        apply_conflict_winner_guard(
            v03_score_tsv_path=scores,
            pool_decisions_tsv_path=decisions,
            pool_manifest_json_path=pool_manifest,
            output_tsv_path=tmp_path / "v04.tsv",
        )


def test_guard_rejects_truth_or_reserved_score_columns(tmp_path: Path) -> None:
    scores, decisions, pool_manifest = _write_inputs(tmp_path)
    text = scores.read_text(encoding="utf-8")
    lines = text.splitlines()
    lines[0] += "\tlabel_exact_cds"
    lines[1:] = [line + "\t0" for line in lines[1:]]
    scores.write_text("\n".join(lines) + "\n", encoding="utf-8")
    score_manifest = Path(str(scores) + ".manifest.json")
    manifest = json.loads(score_manifest.read_text(encoding="utf-8"))
    manifest["outputs"]["scores"]["sha256"] = _sha256(scores)
    score_manifest.write_text(json.dumps(manifest), encoding="utf-8")
    with pytest.raises(ValueError, match="Truth/label"):
        apply_conflict_winner_guard(
            v03_score_tsv_path=scores,
            pool_decisions_tsv_path=decisions,
            pool_manifest_json_path=pool_manifest,
            output_tsv_path=tmp_path / "v04.tsv",
        )


def test_guard_rejects_conflict_member_count_mismatch(tmp_path: Path) -> None:
    scores, decisions, pool_manifest = _write_inputs(tmp_path)
    text = decisions.read_text(encoding="utf-8").replace(
        "a\tx\t2\taccepted", "a\tx\t3\taccepted"
    )
    decisions.write_text(text, encoding="utf-8")
    manifest = json.loads(pool_manifest.read_text(encoding="utf-8"))
    manifest["outputs"]["decisions"]["sha256"] = _sha256(decisions)
    pool_manifest.write_text(json.dumps(manifest), encoding="utf-8")
    with pytest.raises(ValueError, match="Conflict member count differs"):
        apply_conflict_winner_guard(
            v03_score_tsv_path=scores,
            pool_decisions_tsv_path=decisions,
            pool_manifest_json_path=pool_manifest,
            output_tsv_path=tmp_path / "v04.tsv",
        )


@pytest.mark.parametrize(
    ("digests", "conflicts", "baseline", "primary", "message"),
    (
        (
            ("a", "a"),
            {"a": ""},
            {"a": 0.0},
            {"a": 0.0},
            "nonempty and unique",
        ),
        (
            ("a", "b"),
            {"a": "x", "b": ""},
            {"a": 0.0, "b": 0.0},
            {"a": 0.0, "b": 0.0},
            "singleton",
        ),
        (
            ("a", "b"),
            {"a": "x", "b": "x"},
            {"a": float("nan"), "b": 0.0},
            {"a": 0.0, "b": 0.0},
            "finite",
        ),
    ),
)
def test_guard_rejects_invalid_primitive_inputs(
    digests: tuple[str, ...],
    conflicts: dict[str, str],
    baseline: dict[str, float],
    primary: dict[str, float],
    message: str,
) -> None:
    with pytest.raises(ValueError, match=message):
        compute_conflict_winner_guard(
            digests=digests,
            conflict_by_digest=conflicts,
            baseline_scores=baseline,
            primary_scores=primary,
        )


def test_guard_rejects_reserved_pool_columns(tmp_path: Path) -> None:
    scores, decisions, pool_manifest = _write_inputs(tmp_path)
    lines = decisions.read_text(encoding="utf-8").splitlines()
    lines[0] += "\tv04_primary_rank_score"
    lines[1:] = [line + "\t0" for line in lines[1:]]
    decisions.write_text("\n".join(lines) + "\n", encoding="utf-8")
    manifest = json.loads(pool_manifest.read_text(encoding="utf-8"))
    manifest["outputs"]["decisions"]["sha256"] = _sha256(decisions)
    pool_manifest.write_text(json.dumps(manifest), encoding="utf-8")
    with pytest.raises(ValueError, match="reserved v0.4"):
        apply_conflict_winner_guard(
            v03_score_tsv_path=scores,
            pool_decisions_tsv_path=decisions,
            pool_manifest_json_path=pool_manifest,
            output_tsv_path=tmp_path / "v04.tsv",
        )


def test_guard_rejects_nonfrozen_pool_policy(tmp_path: Path) -> None:
    scores, decisions, pool_manifest = _write_inputs(tmp_path)
    manifest = json.loads(pool_manifest.read_text(encoding="utf-8"))
    manifest["parameters"]["redundancy_policy"] = "suppress_overlap"
    pool_manifest.write_text(json.dumps(manifest), encoding="utf-8")
    with pytest.raises(ValueError, match="chain-preserving policy"):
        apply_conflict_winner_guard(
            v03_score_tsv_path=scores,
            pool_decisions_tsv_path=decisions,
            pool_manifest_json_path=pool_manifest,
            output_tsv_path=tmp_path / "v04.tsv",
        )


def test_guard_refuses_output_overwrite(tmp_path: Path) -> None:
    scores, decisions, pool_manifest = _write_inputs(tmp_path)
    output = tmp_path / "v04.tsv"
    apply_conflict_winner_guard(
        v03_score_tsv_path=scores,
        pool_decisions_tsv_path=decisions,
        pool_manifest_json_path=pool_manifest,
        output_tsv_path=output,
    )
    with pytest.raises(FileExistsError, match="Refusing to overwrite"):
        apply_conflict_winner_guard(
            v03_score_tsv_path=scores,
            pool_decisions_tsv_path=decisions,
            pool_manifest_json_path=pool_manifest,
            output_tsv_path=output,
        )


def test_conflict_guard_cli_requires_and_uses_pool_manifest(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    scores, decisions, pool_manifest = _write_inputs(tmp_path)
    output = tmp_path / "v04.tsv"
    assert (
        cli_main(
            (
                "evidence",
                "apply-conflict-winner-guard",
                "--v03-scores",
                str(scores),
                "--pool-decisions",
                str(decisions),
                "--pool-manifest",
                str(pool_manifest),
                "--output-tsv",
                str(output),
            )
        )
        == 0
    )
    manifest = json.loads(capsys.readouterr().out)
    assert manifest["inputs"]["pool_manifest"] == _sha256(pool_manifest)
    assert output.is_file()
