#!/usr/bin/env python3
"""Freeze Walnut core-H1 inputs, protocol bytes and implementation skeleton."""
from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path
import runpy
import shutil
import subprocess
from typing import Any, Iterable

from ploidypatch.artifact_manifest import (
    sha256_file,
    verify_sha256sums,
    write_sha256sums,
)
from ploidypatch.holdout_contract import (
    CORE_H1_MODEL_VERSION,
    CoreH1ScientificParameters,
    load_holdout_contract,
)


SCHEMA_VERSION = "ploidypatch.walnut_core_h1_protocol_freeze.v0.8"
HOLDOUT_ID = "walnut_walnut2_v0.8"
FULL_SHA_LENGTH = 40
PROTOCOL_FILES = (
    "config/holdouts/walnut_walnut2_v0.8/contract.json",
    "config/holdouts/walnut_walnut2_v0.8/no_ranker/composite_manifest.json",
    "config/holdouts/walnut_walnut2_v0.8/no_ranker/SHA256SUMS",
    "config/walnut_external_input_sources_v0.8.tsv",
    "config/walnut_external_validation_policy_v0.8.tsv",
    "config/walnut_external_event_definition_v0.8.tsv",
    "config/primary_seqids/juglans_regia_walnut_2.0.tsv",
    "config/primary_seqids/juglans_mandshurica_gwhbeun_v1.tsv",
    "config/primary_seqids/carya_illinoinensis_pawnee_v1.tsv",
    "config/primary_seqids/corylus_avellana_cavtom2pms_1.0.tsv",
    "config/primary_seqids/castanea_mollissima_nanking_hap2_v1.0.tsv",
    "docs/WALNUT_EXTERNAL_SELECTION_AND_CONTAMINATION_RATIONALE_v0.8.md",
    "docs/WALNUT_EXTERNAL_CORE_H1_PROTOCOL_v0.8.md",
)
IMPLEMENTATION_FILES = (
    "src/ploidypatch/artifact_manifest.py",
    "src/ploidypatch/holdout_contract.py",
    "src/ploidypatch/safe_tar_fasta.py",
    "scripts/stage_external_holdout_inputs_v0.5.py",
    "scripts/preflight_walnut_external_inputs_v0.8.py",
    "scripts/freeze_walnut_external_protocol_v0.8.py",
)
FORBIDDEN_TARGET_ARTIFACTS = (
    "data/derived/external_evaluator/walnut_walnut2_v0.8",
    "data/derived/external_inputs/walnut/v0.8",
    "benchmark/structure/walnut_walnut2_v0.8",
    "results/evaluator/walnut_walnut2_v0.8",
    "results/copy_collapse/external/walnut_walnut2_v0.8",
    "results/holdouts/walnut_walnut2_v0.8",
    "work/walnut_walnut2_v0.8",
)


def verify_git_state(code_root: Path, code_commit: str) -> None:
    if (
        len(code_commit) != FULL_SHA_LENGTH
        or any(character not in "0123456789abcdef" for character in code_commit)
    ):
        raise ValueError("code_commit must be a full lowercase Git SHA")
    head = subprocess.run(
        ["git", "-C", str(code_root), "rev-parse", "HEAD"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    if head != code_commit:
        raise ValueError(f"code_commit differs from code-root HEAD: {head}")
    status = subprocess.run(
        ["git", "-C", str(code_root), "status", "--porcelain"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    if status:
        raise ValueError("Code root must be completely clean before protocol freeze")


def require_regular_files(code_root: Path, relatives: Iterable[str]) -> None:
    for relative in relatives:
        path = code_root / relative
        if not path.is_file() or path.is_symlink() or path.stat().st_size == 0:
            raise ValueError(f"Missing, empty or symlinked frozen file: {relative}")


def reject_too_late(project_root: Path) -> None:
    for relative in FORBIDDEN_TARGET_ARTIFACTS:
        path = project_root / relative
        if path.exists() or path.is_symlink():
            raise ValueError(f"Walnut protocol freeze is too late: {relative}")


def write_tsv(
    path: Path, header: tuple[str, ...], rows: Iterable[Iterable[Any]]
) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def chmod_read_only(root: Path) -> None:
    if os.name == "nt":
        return
    for path in sorted(root.rglob("*"), reverse=True):
        path.chmod(0o550 if path.is_dir() else 0o440)
    root.chmod(0o550)


def freeze_protocol(
    *,
    project_root: Path,
    code_root: Path,
    contract_path: Path,
    staged_inputs: Path,
    code_commit: str,
    output: Path,
) -> Path:
    partial = Path(str(output) + ".partial")
    if output.exists() or output.is_symlink() or partial.exists() or partial.is_symlink():
        raise FileExistsError("Refusing to overwrite Walnut protocol freeze or partial")
    verify_git_state(code_root, code_commit)
    contract = load_holdout_contract(contract_path)
    if (
        contract.holdout_id != HOLDOUT_ID
        or contract.model_version != CORE_H1_MODEL_VERSION
        or not isinstance(contract.scientific_parameters, CoreH1ScientificParameters)
        or contract.scientific_parameters.h2_or_topology_ranking != "forbidden"
        or contract.scientific_parameters.all_arm_collateral_loss_maximum != 0
    ):
        raise ValueError("Freeze requires the exact Walnut core H1-only contract")
    canonical_contract = code_root / PROTOCOL_FILES[0]
    if contract_path.resolve() != canonical_contract.resolve():
        raise ValueError("Freeze requires the canonical committed Walnut contract")
    require_regular_files(code_root, (*PROTOCOL_FILES, *IMPLEMENTATION_FILES))
    reject_too_late(project_root)

    preflight_path = code_root / "scripts/preflight_walnut_external_inputs_v0.8.py"
    preflight_namespace = runpy.run_path(str(preflight_path), run_name="walnut_preflight")
    report = preflight_namespace["run_preflight"](
        project_root=code_root,
        contract_path=contract_path,
        staged_inputs=staged_inputs,
    )
    if report.get("staged_inputs_verified") is not True:
        raise ValueError("Walnut stage was not verified")
    role_contract = json.loads(
        (staged_inputs / "role_contract.json").read_text(encoding="utf-8")
    )
    if role_contract.get("code_commit") != code_commit:
        raise ValueError("Staged inputs were not bound to this code commit")
    staged_checksum_sha256 = sha256_file(staged_inputs / "SHA256SUMS")

    partial.mkdir(parents=True)
    try:
        shutil.copyfile(contract_path, partial / "contract.json")
        shutil.copyfile(
            staged_inputs / "role_contract.json", partial / "role_contract.json"
        )
        shutil.copyfile(
            staged_inputs / "role_manifest.tsv", partial / "role_manifest.tsv"
        )
        shutil.copyfile(
            staged_inputs / "SHA256SUMS", partial / "staged_input_SHA256SUMS"
        )
        verify_sha256sums(staged_inputs, ignore_checksum_file=True)
        if (
            sha256_file(staged_inputs / "SHA256SUMS") != staged_checksum_sha256
            or sha256_file(partial / "role_contract.json")
            != sha256_file(staged_inputs / "role_contract.json")
            or sha256_file(partial / "role_manifest.tsv")
            != sha256_file(staged_inputs / "role_manifest.tsv")
        ):
            raise ValueError("Staged inputs changed during Walnut protocol freeze")
        artifact_root = partial / "protocol_artifacts"
        for relative in PROTOCOL_FILES[1:]:
            destination = artifact_root / relative
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(code_root / relative, destination)

        frozen_files = sorted(set((*PROTOCOL_FILES, *IMPLEMENTATION_FILES)))
        write_tsv(
            partial / "implementation_manifest.tsv",
            ("relative_path", "bytes", "sha256"),
            (
                (relative, (code_root / relative).stat().st_size,
                 sha256_file(code_root / relative))
                for relative in frozen_files
            ),
        )
        write_tsv(
            partial / "forbidden_target_artifacts.tsv",
            ("relative_path",),
            ((relative,) for relative in FORBIDDEN_TARGET_ARTIFACTS),
        )
        write_tsv(
            partial / "run_contract.tsv",
            ("field", "value"),
            (
                ("schema_version", SCHEMA_VERSION),
                ("holdout_id", contract.holdout_id),
                ("policy_id", contract.policy_id),
                ("protocol_profile", "core_H1_only_no_ranker"),
                ("model_version", contract.model_version),
                ("code_commit", code_commit),
                ("freeze_stage", "post_metadata_pre_pair_pre_candidate_pre_label"),
                ("wgd_pairs_enumerated_before_freeze", "false"),
                ("candidate_counts_computed_before_freeze", "false"),
                ("truth_labels_accessed_before_freeze", "false"),
                ("H2_or_topology_ranking", "forbidden"),
                ("bootstrap_replicates", "20000"),
                ("all_arm_collateral_loss_maximum", "0"),
            ),
        )
        no_ranker_root = code_root / "config/holdouts/walnut_walnut2_v0.8/no_ranker"
        verify_sha256sums(no_ranker_root, ignore_checksum_file=True)
        manifest = {
            "schema_version": SCHEMA_VERSION,
            "holdout_id": contract.holdout_id,
            "policy_id": contract.policy_id,
            "test_role": contract.test_role,
            "protocol_profile": "core_H1_only_no_ranker",
            "model_version": contract.model_version,
            "ranker_enabled": False,
            "h2_or_topology_ranking_enabled": False,
            "code_commit": code_commit,
            "freeze_stage": "post_metadata_pre_pair_pre_candidate_pre_label",
            "truth_access": False,
            "wgd_pairs_enumerated": False,
            "candidate_counts_computed": False,
            "truth_labels_accessed": False,
            "reference_roles": {"target": 1, "candidate_reference": 2,
                                "evaluator_reference": 2},
            "contract_sha256": sha256_file(partial / "contract.json"),
            "staged_input_SHA256SUMS_sha256": sha256_file(
                staged_inputs / "SHA256SUMS"
            ),
            "staged_role_contract_sha256": sha256_file(
                staged_inputs / "role_contract.json"
            ),
            "no_ranker_SHA256SUMS_sha256": sha256_file(
                no_ranker_root / "SHA256SUMS"
            ),
            "implementation_file_count": len(frozen_files),
            "all_arm_collateral_loss_maximum": 0,
            "formal_runner_frozen": False,
            "formal_runner_status": "deferred_not_part_of_core_protocol_skeleton",
            "forbidden_target_artifacts": list(FORBIDDEN_TARGET_ARTIFACTS),
        }
        (partial / "protocol_manifest.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        write_sha256sums(partial)
        verify_sha256sums(partial, ignore_checksum_file=True)
        chmod_read_only(partial)
        os.replace(partial, output)
    except BaseException:
        if partial.exists() and not partial.is_symlink():
            shutil.rmtree(partial)
        raise
    return output


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", required=True)
    parser.add_argument("--code-root", required=True)
    parser.add_argument("--contract", required=True)
    parser.add_argument("--staged-inputs", required=True)
    parser.add_argument("--code-commit", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args(argv)
    freeze_protocol(
        project_root=Path(args.project_root).resolve(),
        code_root=Path(args.code_root).resolve(),
        contract_path=Path(args.contract).resolve(),
        staged_inputs=Path(args.staged_inputs).resolve(),
        code_commit=args.code_commit,
        output=Path(args.output_dir).resolve(),
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
