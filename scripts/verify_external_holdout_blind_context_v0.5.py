#!/usr/bin/env python3
"""Fail-closed validation shared by v0.5 blind candidate-side scripts."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
from typing import Any, Iterable

from ploidypatch.artifact_manifest import verify_sha256sums
from ploidypatch.holdout_contract import (
    CANDIDATE_ROLE,
    EVALUATOR_ROLE,
    TARGET_ROLE,
    load_holdout_contract,
    staged_relative_path,
)


STAGE_SCHEMA = "ploidypatch.external_holdout_input_stage.v0.5"
PROTOCOL_SCHEMA = "ploidypatch.external_holdout_protocol_freeze.v0.5"
EXECUTION_SCHEMA = "ploidypatch.external_holdout_execution_freeze.v0.5"
MODEL_SCHEMA = "ploidypatch.composite_ranker.v0.4"
ROLE_FIELDS = (
    "role",
    "species_id",
    "release",
    "bundle_id",
    "wgdi_prefix",
    "artifact",
    "bytes",
    "sha256",
    "source_relative_path",
    "staged_relative_path",
    "staged_sha256",
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(8 * 1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"JSON root must be an object: {path}")
    return value


def _assert_regular(path: Path, description: str) -> None:
    if not path.is_file() or path.is_symlink() or path.stat().st_size == 0:
        raise ValueError(f"{description} must be a non-empty non-symlink file: {path}")
    parent = path.parent
    while parent != parent.parent:
        if parent.is_symlink():
            raise ValueError(f"{description} has a symlinked parent: {parent}")
        parent = parent.parent


def _read_role_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != ROLE_FIELDS:
            raise ValueError("Role manifest header differs from the strict v0.5 schema")
        rows = list(reader)
    if len(rows) != 15:
        raise ValueError("Role manifest must contain exactly five three-artifact references")
    keys = [(row["species_id"], row["artifact"]) for row in rows]
    if len(keys) != len(set(keys)):
        raise ValueError("Role manifest contains duplicate species/artifact rows")
    return rows


def verify_context(
    *,
    input_root: Path,
    contract_path: Path,
    protocol_root: Path,
    execution_root: Path,
    blind_benchmark_root: Path,
    expected_holdout_id: str,
    expected_primary_chromosomes: int,
    model_root: Path | None = None,
) -> dict[str, Any]:
    for path, description in (
        (contract_path, "holdout contract"),
        (input_root / "role_manifest.tsv", "staged role manifest"),
        (input_root / "role_contract.json", "staged role contract"),
        (protocol_root / "protocol_manifest.json", "protocol manifest"),
        (protocol_root / "SHA256SUMS", "protocol checksum manifest"),
        (execution_root / "execution_manifest.json", "execution manifest"),
        (execution_root / "SHA256SUMS", "execution checksum manifest"),
        (blind_benchmark_root / "perturbed.gff3", "sealed blind annotation"),
        (blind_benchmark_root / "blind_manifest.json", "sealed blind manifest"),
        (blind_benchmark_root / "SHA256SUMS", "sealed blind checksum manifest"),
    ):
        _assert_regular(path, description)

    if os.environ.get("PLOIDYPATCH_BLIND_RUNNER") != "1":
        raise ValueError("Blind context verification requires PLOIDYPATCH_BLIND_RUNNER=1")
    if os.environ.get("PLOIDYPATCH_NETWORK_ACCESS") != "none":
        raise ValueError("Blind context verification requires PLOIDYPATCH_NETWORK_ACCESS=none")
    code_commit = os.environ.get("PLOIDYPATCH_CODE_COMMIT", "")
    if len(code_commit) != 40 or any(character not in "0123456789abcdef" for character in code_commit):
        raise ValueError("Blind context verification requires a frozen 40-hex code commit")

    verify_sha256sums(blind_benchmark_root, ignore_checksum_file=True)
    benchmark_files = {
        path.relative_to(blind_benchmark_root).as_posix()
        for path in blind_benchmark_root.rglob("*")
        if path.is_file()
    }
    if benchmark_files != {"perturbed.gff3", "blind_manifest.json", "SHA256SUMS"}:
        raise ValueError("Blind benchmark role contains an unexpected file universe")
    benchmark = load_json(blind_benchmark_root / "blind_manifest.json")
    perturbed_sha = sha256_file(blind_benchmark_root / "perturbed.gff3")
    target_genome = benchmark.get("target_genome", {})
    if (
        benchmark.get("schema_version") != "ploidypatch.blind_benchmark_input.v0.5"
        or benchmark.get("truth_access") is not False
        or benchmark.get("complete_target_annotation_access") is not False
        or benchmark.get("perturbed_annotation")
        != {"file_name": "perturbed.gff3", "sha256": perturbed_sha}
        or target_genome.get("mount_role") != "shared_target_genome"
        or not isinstance(target_genome.get("sha256"), str)
        or len(target_genome["sha256"]) != 64
    ):
        raise ValueError("Blind benchmark manifest violates its sealed truth-free contract")

    contract = load_holdout_contract(contract_path)
    if contract.holdout_id != expected_holdout_id:
        raise ValueError("Unexpected holdout identity")
    resolved = contract.target_resolved_parameters
    if resolved.primary_chromosome_count != expected_primary_chromosomes:
        raise ValueError("Unexpected target primary chromosome count")
    role_counts = {
        role: sum(reference.role == role for reference in contract.references)
        for role in (TARGET_ROLE, CANDIDATE_ROLE, EVALUATOR_ROLE)
    }
    if role_counts != {TARGET_ROLE: 1, CANDIDATE_ROLE: 2, EVALUATOR_ROLE: 2}:
        raise ValueError("Holdout reference roles differ from 1 target/2 candidate/2 evaluator")

    role_manifest = input_root / "role_manifest.tsv"
    role_contract_path = input_root / "role_contract.json"
    if role_manifest.read_bytes() != (protocol_root / "role_manifest.tsv").read_bytes():
        raise ValueError("Mounted role manifest differs from protocol freeze")
    if role_contract_path.read_bytes() != (protocol_root / "role_contract.json").read_bytes():
        raise ValueError("Mounted role contract differs from protocol freeze")

    rows = _read_role_rows(role_manifest)
    indexed = {(row["species_id"], row["artifact"]): row for row in rows}
    visible_artifacts: list[dict[str, Any]] = []
    for reference in contract.references:
        for artifact_name, artifact in reference.artifact_items():
            row = indexed.get((reference.species_id, artifact_name))
            if row is None:
                raise ValueError(f"Missing role row: {reference.species_id}/{artifact_name}")
            expected_relative = staged_relative_path(reference, artifact_name).as_posix()
            expected = {
                "role": reference.role,
                "release": reference.release,
                "bundle_id": reference.bundle_id,
                "wgdi_prefix": reference.wgdi_prefix,
                "bytes": str(artifact.bytes),
                "sha256": artifact.sha256,
                "source_relative_path": artifact.source_relative_path.as_posix(),
                "staged_relative_path": expected_relative,
                "staged_sha256": artifact.sha256,
            }
            if any(row[key] != value for key, value in expected.items()):
                raise ValueError(f"Role row differs from contract: {reference.species_id}/{artifact_name}")
            is_visible = reference.role == CANDIDATE_ROLE or (
                reference.role == TARGET_ROLE and artifact_name == "genome"
            )
            if not is_visible:
                continue
            path = input_root / Path(expected_relative)
            _assert_regular(path, "candidate-safe staged artifact")
            if path.stat().st_size != artifact.bytes or sha256_file(path) != artifact.sha256:
                raise ValueError(f"Visible staged artifact differs: {expected_relative}")
            visible_artifacts.append(
                {"species_id": reference.species_id, "artifact": artifact_name, "path": expected_relative}
            )
    if len(visible_artifacts) != 7:
        raise ValueError("Blind input root must expose one target genome and six candidate files")

    for forbidden in (
        input_root / "evaluator_only",
        input_root / "target_complete",
        input_root / "truth",
        input_root / "labels",
        Path("/nas_data"),
    ):
        if forbidden.exists():
            raise ValueError(f"Forbidden blind path is visible: {forbidden}")

    role_contract = load_json(role_contract_path)
    contract_sha = sha256_file(contract_path)
    if (
        role_contract.get("schema_version") != STAGE_SCHEMA
        or role_contract.get("holdout_id") != contract.holdout_id
        or role_contract.get("policy_id") != contract.policy_id
        or role_contract.get("test_role") != contract.test_role
        or role_contract.get("model_version") != contract.model_version
        or role_contract.get("contract", {}).get("sha256") != contract_sha
        or role_contract.get("truth_blind") != dict(contract.truth_blind)
        or role_contract.get("target_resolved_parameters")
        != {
            "primary_chromosome_count": resolved.primary_chromosome_count,
            "minimum_target_chromosomes_fraction": resolved.minimum_target_chromosomes_fraction,
            "minimum_target_chromosomes": resolved.minimum_target_chromosomes,
        }
        or role_contract.get("role_boundaries", {}).get("candidate_only")
        != "candidate_generation_only"
        or role_contract.get("role_boundaries", {}).get("evaluator_only")
        != "complete_target_annotation_and_evaluator_truth_references"
    ):
        raise ValueError("Staged role contract differs from holdout contract")

    protocol = load_json(protocol_root / "protocol_manifest.json")
    execution = load_json(execution_root / "execution_manifest.json")
    protocol_sums_sha = sha256_file(protocol_root / "SHA256SUMS")
    patch = execution.get("execution_patch")
    patch_mode = isinstance(patch, dict)
    science_commit = patch.get("base_code_commit") if patch_mode else code_commit
    expected_stage = (
        "post_evaluator_truth_failed_blind_pre_candidate_pre_score_pre_label_execution_patch"
        if patch_mode
        else "post_metadata_pre_pair_pre_candidate_pre_label"
    )
    expected_created_before = {
        "wgd_pair_enumeration": not patch_mode,
        "candidate_generation": True,
        "candidate_labels": True,
        "candidate_scores": True,
    }
    patch_valid = True
    if patch_mode:
        changed = patch.get("changed_files")
        failed = patch.get("failed_attempt")
        patch_valid = (
            patch.get("schema_version")
            == "ploidypatch.external_holdout_execution_patch.v0.5"
            and patch.get("freeze_stage") == expected_stage
            and patch.get("patch_code_commit") == code_commit
            and patch.get("base_protocol_SHA256SUMS_sha256") == protocol_sums_sha
            and isinstance(
                patch.get("superseded_execution_SHA256SUMS_sha256"), str
            )
            and len(patch["superseded_execution_SHA256SUMS_sha256"]) == 64
            and patch.get("contract_sha256") == contract_sha
            and patch.get("staged_input_SHA256SUMS_sha256")
            == protocol.get("staged_input_SHA256SUMS_sha256")
            and patch.get("scientific_protocol_changed") is False
            and patch.get("contract_or_policy_changed") is False
            and patch.get("model_or_threshold_changed") is False
            and patch.get("staged_inputs_changed") is False
            and patch.get("truth_or_benchmark_regenerated") is False
            and patch.get("evaluator_truth_construction_completed_before_patch") is True
            and patch.get("blind_candidate_wgd_completed_before_patch") is False
            and patch.get("candidate_generation_completed_before_patch") is False
            and patch.get("formal_scores_generated_before_patch") is False
            and patch.get("truth_labels_accessed_before_patch") is False
            and patch.get("automatic_approval") is False
            and isinstance(changed, list)
            and bool(changed)
            and all(
                isinstance(row, dict)
                and row.get("status") in {"A", "M"}
                and isinstance(row.get("relative_path"), str)
                and isinstance(row.get("patch_sha256"), str)
                and len(row["patch_sha256"]) == 64
                for row in changed
            )
            and isinstance(failed, dict)
            and isinstance(failed.get("exit_status"), int)
            and failed["exit_status"] != 0
            and isinstance(failed.get("tree_sha256"), str)
            and len(failed["tree_sha256"]) == 64
            and (execution_root / "superseded_failed_attempt_manifest.tsv").is_file()
            and (execution_root / "patch_reason.md").is_file()
        )
    common = (
        protocol.get("holdout_id") == contract.holdout_id
        and execution.get("holdout_id") == contract.holdout_id
        and protocol.get("policy_id") == contract.policy_id
        and execution.get("policy_id") == contract.policy_id
        and protocol.get("test_role") == contract.test_role
        and execution.get("test_role") == contract.test_role
        and protocol.get("model_version") == contract.model_version
        and execution.get("model_version") == contract.model_version
        and protocol.get("code_commit") == science_commit
        and execution.get("code_commit") == code_commit
        and role_contract.get("code_commit") == science_commit
        and protocol.get("contract_sha256") == contract_sha
        and execution.get("contract_sha256") == contract_sha
    )
    if (
        protocol.get("schema_version") != PROTOCOL_SCHEMA
        or execution.get("schema_version") != EXECUTION_SCHEMA
        or not common
        or protocol.get("truth_access") is not False
        or protocol.get("freeze_stage")
        != "post_metadata_pre_pair_pre_candidate_pre_label"
        or protocol.get("wgd_pairs_enumerated") is not False
        or protocol.get("candidate_counts_computed") is not False
        or protocol.get("truth_labels_accessed") is not False
        or execution.get("freeze_stage") != expected_stage
        or execution.get("created_before") != expected_created_before
        or not patch_valid
        or execution.get("network_access_in_blind_runner") is not False
        or execution.get("nas_data_mount_in_blind_runner") is not False
        or execution.get("complete_target_annotation_mount_in_blind_runner") is not False
        or execution.get("evaluator_only_mount_in_blind_runner") is not False
        or execution.get("truth_or_label_mount_in_blind_runner") is not False
        or execution.get("protocol_SHA256SUMS_sha256") != protocol_sums_sha
        or protocol.get("staged_role_contract_sha256") != sha256_file(role_contract_path)
    ):
        raise ValueError("Protocol/execution freeze violates the blind context contract")

    if model_root is not None:
        for name in ("SHA256SUMS", "composite_manifest.json"):
            _assert_regular(model_root / name, "frozen composite model artifact")
        model_manifest = load_json(model_root / "composite_manifest.json")
        model_sha = sha256_file(model_root / "SHA256SUMS")
        if (
            contract.model_version != "PloidyPatch_ranker_v0.4"
            or model_manifest.get("schema_version") != MODEL_SCHEMA
            or model_manifest.get("automatic_approval") is not False
            or protocol.get("composite_model_SHA256SUMS_sha256") != model_sha
            or execution.get("composite_model_SHA256SUMS_sha256") != model_sha
            or (
                patch_mode
                and patch.get("composite_model_SHA256SUMS_sha256") != model_sha
            )
        ):
            raise ValueError("Frozen composite model differs from the exact safe v0.4 model")

    return {
        "schema_version": "ploidypatch.external_holdout_blind_context_check.v0.5",
        "holdout_id": contract.holdout_id,
        "contract_sha256": contract_sha,
        "protocol_SHA256SUMS_sha256": protocol_sums_sha,
        "execution_SHA256SUMS_sha256": sha256_file(execution_root / "SHA256SUMS"),
        "blind_benchmark_SHA256SUMS_sha256": sha256_file(
            blind_benchmark_root / "SHA256SUMS"
        ),
        "visible_candidate_safe_artifacts": visible_artifacts,
        "target_primary_chromosomes": resolved.primary_chromosome_count,
        "candidate_truth_access": False,
        "evaluator_artifact_access": False,
        "complete_target_annotation_access": False,
        "nas_data_access": False,
        "network_access": False,
        "execution_patch": (
            {
                "active": True,
                "freeze_stage": patch["freeze_stage"],
                "base_code_commit": patch["base_code_commit"],
                "patch_code_commit": patch["patch_code_commit"],
                "failed_attempt_tree_sha256": patch["failed_attempt"]["tree_sha256"],
            }
            if patch_mode
            else {"active": False}
        ),
    }


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", required=True, type=Path)
    parser.add_argument("--contract", required=True, type=Path)
    parser.add_argument("--protocol-freeze", required=True, type=Path)
    parser.add_argument("--execution-freeze", required=True, type=Path)
    parser.add_argument("--blind-benchmark-root", required=True, type=Path)
    parser.add_argument("--expected-holdout-id", required=True)
    parser.add_argument("--expected-primary-chromosomes", required=True, type=int)
    parser.add_argument("--model-freeze", type=Path)
    parser.add_argument("--output-json", type=Path)
    args = parser.parse_args(argv)
    report = verify_context(
        input_root=args.input_root,
        contract_path=args.contract,
        protocol_root=args.protocol_freeze,
        execution_root=args.execution_freeze,
        blind_benchmark_root=args.blind_benchmark_root,
        expected_holdout_id=args.expected_holdout_id,
        expected_primary_chromosomes=args.expected_primary_chromosomes,
        model_root=args.model_freeze,
    )
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output_json is None:
        print(rendered, end="")
    else:
        if args.output_json.exists() or args.output_json.is_symlink():
            raise FileExistsError("Refusing to overwrite blind context report")
        args.output_json.parent.mkdir(parents=True, exist_ok=True)
        args.output_json.write_text(rendered, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
