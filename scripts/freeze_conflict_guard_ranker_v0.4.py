#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
from pathlib import Path
from typing import Any


SCHEMA_VERSION = "ploidypatch.composite_ranker.v0.4"
EXPECTED_V03_MODEL_SHA256 = (
    "9dcafcc3294bfd85ccfbe824c33debd0e967578b89cf7f08b56076d32efd4390"
)
CODE_FILES = (
    "src/ploidypatch/conflict_guard.py",
    "src/ploidypatch/support_ranker.py",
    "src/ploidypatch/cli.py",
    "scripts/evaluate_conflict_winner_guard_v0.4.py",
    "scripts/verify_conflict_guard_production_replay_v0.4.py",
    "scripts/freeze_conflict_guard_ranker_v0.4.py",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"JSON root is not an object: {path}")
    return value


def verify_sha256sums(root: Path) -> dict[str, str]:
    checksum_path = root / "SHA256SUMS"
    if not checksum_path.is_file() or checksum_path.stat().st_size == 0:
        raise ValueError(f"Missing SHA256SUMS: {root}")
    entries: dict[str, str] = {}
    for line_number, line in enumerate(
        checksum_path.read_text(encoding="utf-8").splitlines(), start=1
    ):
        digest, separator, name = line.partition("  ")
        if not separator or not re.fullmatch(r"[0-9a-f]{64}", digest):
            raise ValueError(f"Malformed SHA256SUMS line {line_number}: {root}")
        path = root / name
        if not path.is_file() or sha256(path) != digest:
            raise ValueError(f"SHA256SUMS verification failed: {path}")
        entries[name] = digest
    return entries


def guard_policy(code_commit: str) -> dict[str, Any]:
    return {
        "schema_version": "ploidypatch.conflict_winner_guard_policy.v1",
        "code_commit": code_commit,
        "base_model": {
            "version": "v0.3",
            "sha256": EXPECTED_V03_MODEL_SHA256,
        },
        "policy": {
            "name": "baseline_fallback_on_conflict_winner_disagreement",
            "winner_tie_break": "descending_score_then_candidate_digest",
            "outside_guarded_sets": "retain_v03_primary_score",
            "inside_guarded_sets": "use_v03_baseline_logit_and_abstain_topology",
            "automatic_approval": False,
            "calibrated_probability": False,
            "decision_threshold": None,
            "interpretation": "uncalibrated_review_rank_not_probability",
        },
        "required_pool": {
            "schema_version": "ploidypatch.method_candidate_pool.v2",
            "conflict_action": "retain_all_for_ranking_and_review",
            "max_redundancy_overlap": 0.5,
            "min_method_support": 1,
            "redundancy_policy": "retain_distinct_chains",
        },
        "production_replay": {
            "score_acceptance": (
                "absolute_difference_at_most_1e-12_or_ulp_difference_at_most_32"
            ),
            "candidate_universe_exact": True,
            "ranking_exact": True,
            "guard_flags_exact": True,
            "conflict_winner_mapping_exact": True,
        },
        "external_safety_gates": {
            "all_conflict_winners_identical_to_baseline": True,
            "top_1pct_true_positives_not_below_baseline": True,
            "minimum_positive_topology_coverage": 0.70,
            "minimum_v03_ap_gain_retained_when_positive": 0.90,
            "automatic_approvals": 0,
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Freeze the composite PloidyPatch v0.4 ranker")
    parser.add_argument("--code-root", required=True)
    parser.add_argument("--v03-model-root", required=True)
    parser.add_argument("--development-evaluation-root", required=True)
    parser.add_argument("--production-replay-json", required=True)
    parser.add_argument("--source-archive", required=True)
    parser.add_argument("--environment-lock", required=True)
    parser.add_argument("--code-commit", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()
    if not re.fullmatch(r"[0-9a-f]{40}", args.code_commit):
        raise ValueError("--code-commit must be a full lowercase Git SHA")

    code_root = Path(args.code_root)
    model_root = Path(args.v03_model_root)
    evaluation_root = Path(args.development_evaluation_root)
    replay_path = Path(args.production_replay_json)
    source_archive = Path(args.source_archive)
    environment_lock = Path(args.environment_lock)
    output = Path(args.output_dir)
    partial = Path(str(output) + ".partial")
    required = [
        replay_path,
        source_archive,
        environment_lock,
        *(code_root / relative for relative in CODE_FILES),
    ]
    for path in required:
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f"Missing composite-freeze input: {path}")
    if output.exists() or partial.exists():
        raise FileExistsError("Refusing to overwrite composite v0.4 freeze")

    model_sums = verify_sha256sums(model_root)
    evaluation_sums = verify_sha256sums(evaluation_root)
    model_path = model_root / "model.json"
    if (
        sha256(model_path) != EXPECTED_V03_MODEL_SHA256
        or model_sums.get("model.json") != EXPECTED_V03_MODEL_SHA256
    ):
        raise ValueError("Frozen v0.3 model checksum differs")
    evaluation = load_json(evaluation_root / "evaluation.json")
    if (
        evaluation.get("schema_version")
        != "ploidypatch.conflict_winner_guard_evaluation.v1"
        or evaluation.get("code_commit") != args.code_commit
        or evaluation.get("all_development_gates_pass") is not True
    ):
        raise ValueError("Development evaluation is not an eligible committed freeze")
    for name, report in evaluation.get("datasets", {}).items():
        identity = report.get("conflict_winner_identity", {})
        reviews = report.get("review_budgets", {})
        if (
            identity.get("mismatch_count") != 0
            or identity.get("identical_to_baseline") is not True
            or any("top_250" not in method for method in reviews.values())
        ):
            raise ValueError(f"Development safety or review-budget gate failed: {name}")

    replay = load_json(replay_path)
    if (
        replay.get("schema_version")
        != "ploidypatch.conflict_winner_guard_production_replay.v1"
        or replay.get("code_commit") != args.code_commit
        or replay.get("labels_used") is not False
        or replay.get("all_replay_gates_pass") is not True
        or replay.get("all_rankings_exact") is not True
    ):
        raise ValueError("Production replay is not an eligible committed freeze")
    for name, report in replay.get("datasets", {}).items():
        if (
            report.get("mismatch_count") != 0
            or report.get("winner_audit", {}).get("mismatch_count") != 0
        ):
            raise ValueError(f"Production replay safety gate failed: {name}")

    partial.mkdir(parents=True)
    shutil.copyfile(model_path, partial / "model_v0.3.json")
    shutil.copyfile(replay_path, partial / "replay_audit.json")
    shutil.copyfile(evaluation_root / "evaluation.json", partial / "development_evaluation.json")
    shutil.copyfile(evaluation_root / "SHA256SUMS", partial / "development_SHA256SUMS")
    shutil.copyfile(source_archive, partial / "source.tar.gz")
    shutil.copyfile(environment_lock, partial / "environment.explicit.txt")
    policy = guard_policy(args.code_commit)
    with (partial / "guard_policy.json").open("x", encoding="utf-8", newline="") as handle:
        json.dump(policy, handle, indent=2, sort_keys=True)
        handle.write("\n")

    manifest = {
        "schema_version": SCHEMA_VERSION,
        "model_id": "PloidyPatch_ranker_v0.4",
        "code_commit": args.code_commit,
        "truth_access": False,
        "automatic_approval": False,
        "interpretation": "uncalibrated_review_rank_not_probability",
        "components": {
            "model_v0.3": sha256(partial / "model_v0.3.json"),
            "guard_policy": sha256(partial / "guard_policy.json"),
            "replay_audit": sha256(partial / "replay_audit.json"),
            "development_evaluation": sha256(
                partial / "development_evaluation.json"
            ),
            "development_predictions": evaluation_sums.get("predictions.tsv"),
            "source_archive": sha256(partial / "source.tar.gz"),
            "environment_explicit": sha256(partial / "environment.explicit.txt"),
        },
        "code_files": {
            relative: sha256(code_root / relative) for relative in CODE_FILES
        },
        "input_freezes": {
            "v03_model_SHA256SUMS": sha256(model_root / "SHA256SUMS"),
            "development_SHA256SUMS": sha256(evaluation_root / "SHA256SUMS"),
        },
    }
    with (partial / "composite_manifest.json").open(
        "x", encoding="utf-8", newline=""
    ) as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
    checksum_path = partial / "SHA256SUMS"
    with checksum_path.open("x", encoding="utf-8", newline="") as handle:
        for path in sorted(partial.iterdir()):
            if path != checksum_path:
                handle.write(f"{sha256(path)}  {path.name}\n")
    os.replace(partial, output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
