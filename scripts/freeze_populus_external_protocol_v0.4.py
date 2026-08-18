#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import shutil
from pathlib import Path
from typing import Any


SCHEMA_VERSION = "ploidypatch.populus_external_protocol_freeze.v0.4"
POLICY_ID = "ploidypatch_populus_external_validation_v0.4"
IMPLEMENTATION_FILES = (
    "scripts/evaluate_external_v0.4.py",
    "scripts/evaluate_apple_external_v0.3.py",
    "scripts/preflight_external_inputs_v0.4.py",
    "scripts/freeze_populus_external_protocol_v0.4.py",
    "src/ploidypatch/conflict_guard.py",
    "src/ploidypatch/support_ranker.py",
    "src/ploidypatch/cli.py",
)
SELECTION_FILES = (
    "config/populus_external_input_sources_v0.4.tsv",
    "config/v0.4_external_species_contamination_registry.tsv",
    "config/v0.4_external_species_eligibility_registry.tsv",
    "config/v0.4_external_species_selection_rule.tsv",
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
        raise ValueError(f"Expected JSON object: {path}")
    return value


def read_policy(path: Path) -> dict[str, str]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ["field", "value"]:
            raise ValueError("Policy must have field and value columns")
        rows = list(reader)
    output = {row["field"]: row["value"] for row in rows}
    if len(output) != len(rows) or "" in output:
        raise ValueError("Policy contains empty or duplicate fields")
    return output


def verify_sha256sums(root: Path) -> None:
    checksum_path = root / "SHA256SUMS"
    if not checksum_path.is_file() or checksum_path.stat().st_size == 0:
        raise ValueError(f"Missing SHA256SUMS: {root}")
    for line_number, line in enumerate(
        checksum_path.read_text(encoding="utf-8").splitlines(), start=1
    ):
        digest, separator, name = line.partition("  ")
        path = root / name
        if (
            not separator
            or not re.fullmatch(r"[0-9a-f]{64}", digest)
            or not path.is_file()
            or sha256(path) != digest
        ):
            raise ValueError(f"Checksum failure at line {line_number}: {root}")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Freeze Populus selection, data roles, evaluator and failure policy"
    )
    parser.add_argument("--project-root", required=True)
    parser.add_argument("--code-root", required=True)
    parser.add_argument("--composite-model-root", required=True)
    parser.add_argument("--preflight-root", required=True)
    parser.add_argument("--source-archive", required=True)
    parser.add_argument("--environment-lock", required=True)
    parser.add_argument("--code-commit", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()
    if not re.fullmatch(r"[0-9a-f]{40}", args.code_commit):
        raise ValueError("--code-commit must be a full lowercase Git SHA")

    project_root = Path(args.project_root).resolve()
    code_root = Path(args.code_root).resolve()
    composite_root = Path(args.composite_model_root).resolve()
    preflight_root = Path(args.preflight_root).resolve()
    source_archive = Path(args.source_archive)
    environment_lock = Path(args.environment_lock)
    output = Path(args.output_dir)
    partial = Path(str(output) + ".partial")
    if output.exists() or partial.exists():
        raise FileExistsError("Refusing to overwrite Populus protocol freeze")
    for forbidden in (
        project_root / "data/derived/holdout_inputs/populus_v0.4",
        project_root / "results/copy_collapse/external/populus_v0.4",
        project_root / "results/holdouts/populus_v0.4",
        project_root / "work/populus_external_v0.4",
    ):
        if forbidden.exists():
            raise ValueError(f"Protocol freeze is too late; target artifact exists: {forbidden}")

    protocol_path = code_root / "docs/POPULUS_EXTERNAL_VALIDATION_PROTOCOL_v0.4.md"
    policy_path = code_root / "config/populus_external_validation_policy_v0.4.tsv"
    required = [
        protocol_path,
        policy_path,
        source_archive,
        environment_lock,
        *(code_root / path for path in IMPLEMENTATION_FILES),
        *(code_root / path for path in SELECTION_FILES),
        preflight_root / "metadata.json",
        preflight_root / "input_manifest.tsv",
        preflight_root / "SHA256SUMS",
        composite_root / "composite_manifest.json",
        composite_root / "SHA256SUMS",
    ]
    for path in required:
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f"Missing protocol-freeze input: {path}")
    verify_sha256sums(preflight_root)
    verify_sha256sums(composite_root)

    policy = read_policy(policy_path)
    expected_hashes = {
        "input_source_table_sha256": sha256(
            code_root / "config/populus_external_input_sources_v0.4.tsv"
        ),
        "input_preflight_sha256sums_sha256": sha256(
            preflight_root / "SHA256SUMS"
        ),
        "contamination_registry_sha256": sha256(
            code_root / "config/v0.4_external_species_contamination_registry.tsv"
        ),
        "eligibility_registry_sha256": sha256(
            code_root / "config/v0.4_external_species_eligibility_registry.tsv"
        ),
        "selection_rule_sha256": sha256(
            code_root / "config/v0.4_external_species_selection_rule.tsv"
        ),
        "composite_model_sha256sums_sha256": sha256(
            composite_root / "SHA256SUMS"
        ),
    }
    if policy.get("policy_id") != POLICY_ID or any(
        policy.get(key) != value for key, value in expected_hashes.items()
    ):
        raise ValueError("Policy identity or frozen input hashes differ")
    if (
        policy.get("test_role") != "untouched_confirmatory_external_species"
        or policy.get("selection_by_pair_yield_candidate_count_or_model_performance")
        != "forbidden"
        or policy.get("automatic_copy_addition_approval") != "false"
    ):
        raise ValueError("Policy claim boundary differs")
    preflight = load_json(preflight_root / "metadata.json")
    if (
        preflight.get("schema_version")
        != "ploidypatch.external_input_preflight.v0.4"
        or preflight.get("wgd_pairs_enumerated") is not False
        or preflight.get("candidate_counts_computed") is not False
        or preflight.get("truth_labels_accessed") is not False
        or preflight.get("selection_by_pair_yield_or_performance") is not False
        or preflight.get("source_table", {}).get("sha256")
        != expected_hashes["input_source_table_sha256"]
    ):
        raise ValueError("Preflight exceeded metadata-only selection scope")

    partial.mkdir(parents=True)
    shutil.copyfile(protocol_path, partial / "protocol.md")
    shutil.copyfile(policy_path, partial / "policy.tsv")
    shutil.copyfile(preflight_root / "metadata.json", partial / "preflight_metadata.json")
    shutil.copyfile(
        preflight_root / "input_manifest.tsv", partial / "preflight_input_manifest.tsv"
    )
    shutil.copyfile(composite_root / "composite_manifest.json", partial / "composite_model_manifest.json")
    shutil.copyfile(composite_root / "SHA256SUMS", partial / "composite_model_SHA256SUMS")
    shutil.copyfile(source_archive, partial / "source.tar.gz")
    shutil.copyfile(environment_lock, partial / "environment.explicit.txt")
    selection_dir = partial / "selection"
    implementation_dir = partial / "implementation"
    selection_dir.mkdir()
    implementation_dir.mkdir()
    for relative in SELECTION_FILES:
        shutil.copyfile(code_root / relative, selection_dir / Path(relative).name)
    for relative in IMPLEMENTATION_FILES:
        destination = implementation_dir / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(code_root / relative, destination)

    with (partial / "run_contract.tsv").open("x", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("field", "value"))
        for key, value in (
            ("schema_version", SCHEMA_VERSION),
            ("code_commit", args.code_commit),
            ("freeze_stage", "post_metadata_pre_pair_pre_candidate_pre_label"),
            ("wgd_pairs_enumerated_before_freeze", "false"),
            ("candidate_counts_computed_before_freeze", "false"),
            ("truth_labels_accessed_before_freeze", "false"),
            ("primary_target", "Populus_trichocarpa"),
            ("primary_failure_replacement", "forbidden_under_v0.4"),
        ):
            writer.writerow((key, value))
    with (partial / "implementation_manifest.tsv").open(
        "x", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("relative_path", "bytes", "sha256"))
        for relative in (*IMPLEMENTATION_FILES, *SELECTION_FILES):
            path = code_root / relative
            writer.writerow((relative, path.stat().st_size, sha256(path)))
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "code_commit": args.code_commit,
        "policy_id": POLICY_ID,
        "truth_access": False,
        "wgd_pairs_enumerated": False,
        "candidate_counts_computed": False,
        "source_archive_sha256": sha256(partial / "source.tar.gz"),
        "environment_lock_sha256": sha256(partial / "environment.explicit.txt"),
        "preflight_SHA256SUMS_sha256": sha256(preflight_root / "SHA256SUMS"),
        "composite_model_SHA256SUMS_sha256": sha256(composite_root / "SHA256SUMS"),
    }
    with (partial / "protocol_manifest.json").open(
        "x", encoding="utf-8", newline=""
    ) as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
    checksum_path = partial / "SHA256SUMS"
    with checksum_path.open("x", encoding="utf-8", newline="") as handle:
        for path in sorted(item for item in partial.rglob("*") if item.is_file()):
            if path != checksum_path:
                handle.write(f"{sha256(path)}  {path.relative_to(partial).as_posix()}\n")
    os.replace(partial, output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
