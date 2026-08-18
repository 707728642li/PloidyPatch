#!/usr/bin/env python3
"""Freeze no-ranker Walnut core-H1 implementation and environment lineage."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import shutil
import subprocess
import tarfile
from typing import Any, Iterable

from ploidypatch.artifact_manifest import sha256_file, verify_sha256sums, write_sha256sums
from ploidypatch.holdout_contract import safe_relative_path
from ploidypatch.walnut_h1_framework import (
    BLIND_OUTPUTS,
    EXECUTION_SCHEMA,
    FORBIDDEN_PIPELINE_PATH_TOKENS,
    PIPELINE_ENTRIES,
    REQUIRED_ENVIRONMENTS,
    verify_protocol,
)


_FULL_SHA = re.compile(r"[0-9a-f]{40}")
_ENVIRONMENT = re.compile(r"[a-z0-9][a-z0-9_.-]*")
PATCH_FREEZE_STAGE = (
    "post_evaluator_truth_failed_blind_pre_formal_outputs_pre_label_execution_patch"
)
FAILED_STAGES = frozenset(
    {
        "structure_holdout_sentinel",
        "blind_normalization",
        "candidate_methods",
        "candidate_pools",
        "candidate_output_publication",
    }
)
_CANONICAL_BLIND_ROOT = PurePosixPath(
    "results/copy_collapse/external/walnut_v0.8_h1"
)
_FAILED_PUBLICATION_ROOT = PurePosixPath(
    "project/results/copy_collapse/external/walnut/v0.8_h1"
)
ENGINEERING_FILES = frozenset(
    {
        "src/ploidypatch/artifact_manifest.py",
        "src/ploidypatch/holdout_contract.py",
        "src/ploidypatch/walnut_h1_framework.py",
        "src/ploidypatch/published_output.py",
        "scripts/freeze_walnut_external_execution_v0.8.py",
        "scripts/build_walnut_blind_role_root_v0.8.py",
        "scripts/run_walnut_blind_isolated_v0.8.sh",
        "scripts/finalize_walnut_blind_custody_v0.8.py",
        "scripts/run_walnut_external_reveal_v0.8.py",
    }
)
SCIENTIFIC_FILES = frozenset(
    {
        "src/ploidypatch/audit.py",
        "src/ploidypatch/baseline.py",
        "src/ploidypatch/bootstrap.py",
        "src/ploidypatch/cli.py",
        "src/ploidypatch/consensus.py",
        "src/ploidypatch/copy_pair_sampling.py",
        "src/ploidypatch/gff.py",
        "src/ploidypatch/gff_compat.py",
        "src/ploidypatch/homeolog_pairs.py",
        "src/ploidypatch/io.py",
        "src/ploidypatch/normalize.py",
        "src/ploidypatch/perturb.py",
        "src/ploidypatch/score.py",
        "src/ploidypatch/self_wgd_pairs.py",
        "src/ploidypatch/structure_perturb.py",
        "src/ploidypatch/synteny_io.py",
        "src/ploidypatch/walnut_h1.py",
        "src/ploidypatch/wgdi_summary.py",
        "scripts/build_walnut_h1_candidate_pools_v0.8.py",
        "scripts/build_walnut_structure_holdout_v0.8.py",
        "scripts/finalize_walnut_evaluator_only_inputs_v0.8.py",
        "scripts/infer_walnut_external_pairs_v0.8.py",
        "scripts/prepare_walnut_blind_candidate_inputs_v0.8.py",
        "scripts/prepare_walnut_evaluator_wgdi_inputs_v0.8.py",
        "scripts/prepare_walnut_external_normalized_inputs_v0.8.py",
        "scripts/prepare_walnut_target_ks_inputs_v0.8.py",
        "scripts/run_walnut_candidate_methods_v0.8.sh",
        "scripts/verify_published_method_output_v0.8.py",
        "scripts/run_walnut_evaluator_wgdi_v0.8.sh",
        "scripts/build_wgdi_source_alias_gff.py",
        "scripts/run_gemoma_homology.sh",
        "scripts/publish_gemoma_working.sh",
        "scripts/run_lifton_transfer.sh",
    }
)


def parse_environment(value: str) -> tuple[str, Path]:
    name, separator, raw = value.partition("=")
    if not separator or _ENVIRONMENT.fullmatch(name) is None or not raw:
        raise argparse.ArgumentTypeError("environment must be NAME=PREFIX")
    return name, Path(raw).resolve()


def verify_git_state(code_root: Path, code_commit: str) -> None:
    if _FULL_SHA.fullmatch(code_commit) is None:
        raise ValueError("code_commit must be a full lowercase Git SHA")
    head = subprocess.run(
        ["git", "-C", str(code_root), "rev-parse", "HEAD"],
        check=True, capture_output=True, text=True,
    ).stdout.strip()
    if head != code_commit:
        raise ValueError(f"code_commit differs from HEAD: {head}")
    status = subprocess.run(
        ["git", "-C", str(code_root), "status", "--porcelain"],
        check=True, capture_output=True,
    ).stdout
    if status:
        raise ValueError("Code root must be completely clean before execution freeze")


def git_changed_files(
    code_root: Path, base_commit: str, patch_commit: str
) -> list[tuple[str, str]]:
    raw = subprocess.run(
        [
            "git", "-C", str(code_root), "diff", "--name-status",
            "--no-renames", base_commit, patch_commit, "--",
        ],
        check=True, capture_output=True, text=True,
    ).stdout
    changed: list[tuple[str, str]] = []
    for line in raw.splitlines():
        status, separator, relative = line.partition("\t")
        if not separator or status not in {"A", "M"}:
            raise ValueError(f"Walnut execution patch may only add or modify files: {line}")
        safe_relative_path(relative, "Walnut execution patch file")
        changed.append((status, relative))
    if not changed or len({relative for _status, relative in changed}) != len(changed):
        raise ValueError("Walnut execution patch diff is empty or duplicate")
    return changed


def validate_execution_patch(
    *,
    code_root: Path,
    code_commit: str,
    protocol: Path,
    superseded_execution: Path,
    failed_attempt_log: Path,
    failed_stage: str,
    patch_reason: Path,
    allowed_changed_files: list[str],
    failed_attempt_root: Path | None = None,
) -> dict[str, Any]:
    verify_sha256sums(superseded_execution, ignore_checksum_file=True)
    prior_manifest_path = superseded_execution / "execution_manifest.json"
    if (
        not prior_manifest_path.is_file()
        or prior_manifest_path.is_symlink()
        or prior_manifest_path.stat().st_size == 0
    ):
        raise ValueError("Superseded Walnut execution manifest is missing")
    prior = json.loads(prior_manifest_path.read_text(encoding="utf-8"))
    base_commit = prior.get("code_commit")
    if prior.get("schema_version") != EXECUTION_SCHEMA or not isinstance(
        base_commit, str
    ) or _FULL_SHA.fullmatch(base_commit) is None:
        raise ValueError("Superseded Walnut execution is not a valid base freeze")
    prior_patch = prior.get("execution_patch")
    previous_patch_depth = 0
    if prior_patch is not None:
        if (
            not isinstance(prior_patch, dict)
            or prior_patch.get("schema_version")
            != "ploidypatch.walnut_execution_patch.v0.8"
            or prior_patch.get("patch_code_commit") != base_commit
            or prior_patch.get("blind_candidate_generation_completed_before_patch")
            is not False
            or prior_patch.get("formal_scores_generated_before_patch") is not False
            or prior_patch.get("truth_labels_accessed_before_patch") is not False
        ):
            raise ValueError("Superseded Walnut execution patch lineage is unsafe")
        previous_patch_depth = prior_patch.get("chain_depth", 1)
        if not isinstance(previous_patch_depth, int) or previous_patch_depth < 1:
            raise ValueError("Superseded Walnut execution patch depth is invalid")
    if (
        prior.get("protocol_SHA256SUMS_sha256")
        != sha256_file(protocol / "SHA256SUMS")
        or prior.get("contract_sha256") != sha256_file(protocol / "contract.json")
        or prior.get("ranker_or_model_execution") is not False
    ):
        raise ValueError("Superseded Walnut execution lineage differs")
    if (
        not failed_attempt_log.is_file()
        or failed_attempt_log.is_symlink()
        or failed_attempt_log.stat().st_size == 0
    ):
        raise ValueError("Failed Walnut attempt log is missing")
    if failed_stage not in FAILED_STAGES:
        raise ValueError("Failed Walnut attempt stage is invalid")
    failed_attempt_audit: dict[str, Any] | None = None
    if failed_stage == "candidate_output_publication":
        if failed_attempt_root is None:
            raise ValueError(
                "Candidate-output publication patch requires the full failed attempt"
            )
        failed_attempt_audit = audit_failed_candidate_publication(
            failed_attempt_root, prior["blind_outputs"]
        )
        if sha256_file(failed_attempt_log) != failed_attempt_audit["stderr_sha256"]:
            raise ValueError("Failed-attempt log differs from the bound attempt stderr")
    elif failed_attempt_root is not None:
        raise ValueError(
            "Full failed attempt is only valid for candidate-output publication"
        )
    if (
        not patch_reason.is_file()
        or patch_reason.is_symlink()
        or patch_reason.stat().st_size == 0
    ):
        raise ValueError("Walnut execution patch reason is missing")
    try:
        reason_relative = patch_reason.resolve().relative_to(code_root).as_posix()
    except ValueError as error:
        raise ValueError("Walnut patch reason must be inside the code root") from error
    allowed = {
        safe_relative_path(value, "allowed Walnut patch file").as_posix()
        for value in allowed_changed_files
    }
    changed = git_changed_files(code_root, base_commit, code_commit)
    changed_paths = {relative for _status, relative in changed}
    if changed_paths != allowed or reason_relative not in allowed:
        raise ValueError(
            "Walnut patch Git diff must equal the exact changed-file whitelist "
            "and include the patch-reason file"
        )
    for _status, relative in changed:
        if relative != reason_relative and not relative.startswith(
            ("scripts/", "src/ploidypatch/", "tests/")
        ):
            raise ValueError(f"Forbidden non-implementation Walnut patch file: {relative}")
    return {
        "prior": prior,
        "base_commit": base_commit,
        "reason_relative": reason_relative,
        "changed": changed,
        "previous_patch_depth": previous_patch_depth,
        "failed_stage": failed_stage,
        "failed_attempt_audit": failed_attempt_audit,
    }


def audit_failed_candidate_publication(
    attempt: Path, blind_outputs: dict[str, str]
) -> dict[str, Any]:
    """Bind an isolated post-pool/pre-custody failure without reusing its outputs."""
    if not attempt.is_dir() or attempt.is_symlink():
        raise ValueError("Failed Walnut blind attempt must be a real directory")
    required = {
        "exit_status.txt",
        "stdout.log",
        "stderr.log",
        "bwrap_command.txt",
        "mount_manifest.json",
        "namespace_role_validation.json",
        "project/pipeline_commands.tsv",
    }
    rows: list[dict[str, Any]] = []
    for path in sorted(attempt.rglob("*")):
        if path.is_symlink():
            raise ValueError("Failed Walnut blind attempt may not contain symlinks")
        if path.is_file():
            relative = path.relative_to(attempt).as_posix()
            folded = relative.casefold()
            if any(
                token in folded
                for token in (
                    "/candidate_labels.tsv",
                    "/truth/",
                    "/labels/",
                    "/evaluator/",
                    "/evaluation/",
                    "custody_manifest.json",
                    "reveal_authorization",
                )
            ):
                raise ValueError(
                    f"Failed Walnut attempt contains forbidden reveal artifact: {relative}"
                )
            rows.append(
                {
                    "relative_path": relative,
                    "bytes": path.stat().st_size,
                    "sha256": sha256_file(path),
                }
            )
    paths = {row["relative_path"] for row in rows}
    if not required <= paths:
        raise ValueError("Failed Walnut publication attempt lacks isolation/log evidence")
    status = (attempt / "exit_status.txt").read_text(encoding="utf-8").strip()
    stderr = (attempt / "stderr.log").read_text(encoding="utf-8")
    if status != "72" or stderr.strip() != (
        "missing exact Walnut blind output: raw_predictions.manifest.json"
    ):
        raise ValueError("Failed Walnut publication signature differs from preregistered fault")
    command = (attempt / "bwrap_command.txt").read_text(encoding="utf-8")
    if any(flag not in command for flag in ("--unshare-all", "--unshare-net", "--clearenv")):
        raise ValueError("Failed Walnut publication attempt lacks isolation flags")
    mount = json.loads((attempt / "mount_manifest.json").read_text(encoding="utf-8"))
    namespace = json.loads(
        (attempt / "namespace_role_validation.json").read_text(encoding="utf-8")
    )
    roles = {
        item.get("role")
        for item in mount.get("mounts", [])
        if isinstance(item, dict)
        and item.get("role") in {"shared_target", "candidate_only", "blind_benchmark"}
    }
    if (
        roles != {"shared_target", "candidate_only", "blind_benchmark"}
        or namespace.get("shared_target_visible") is not True
        or namespace.get("candidate_only_visible") is not True
        or namespace.get("blind_benchmark_visible") is not True
        or namespace.get("evaluator_only_visible") is not False
        or namespace.get("truth_visible") is not False
        or namespace.get("complete_target_annotation_visible") is not False
        or namespace.get("nas_data_visible") is not False
    ):
        raise ValueError("Failed Walnut publication isolation evidence is incomplete")
    for item in mount.get("mounts", []):
        text = f"{item.get('host_path', '')}\n{item.get('namespace_path', '')}".casefold()
        if any(
            token in text
            for token in ("/nas_data", "evaluator_only", "/truth", "/labels", "target_complete")
        ):
            raise ValueError("Failed Walnut publication mounted forbidden data")
    for name, raw in blind_outputs.items():
        relative = safe_relative_path(raw, f"Walnut blind output {name}").as_posix()
        if name != "command_log" and f"project/{relative}" in paths:
            raise ValueError(
                f"Failed Walnut publication already contains canonical output: {name}"
            )
    wrong_root = attempt.joinpath(*_FAILED_PUBLICATION_ROOT.parts)
    if not wrong_root.is_dir() or wrong_root.is_symlink():
        raise ValueError("Failed Walnut publication lacks the noncanonical pool tree")
    verify_sha256sums(wrong_root, ignore_checksum_file=True)
    canonical_prefix = _CANONICAL_BLIND_ROOT.as_posix() + "/"
    for name, raw in blind_outputs.items():
        if name == "command_log":
            continue
        if not raw.startswith(canonical_prefix):
            raise ValueError(f"Unexpected canonical Walnut blind output: {name}")
        suffix = PurePosixPath(raw[len(canonical_prefix):])
        path = wrong_root.joinpath(*suffix.parts)
        if not path.is_file() or path.is_symlink() or path.stat().st_size == 0:
            raise ValueError(f"Failed Walnut publication lacks completed pool output: {name}")
    digest = hashlib.sha256()
    for row in rows:
        digest.update(
            f"{row['relative_path']}\0{row['bytes']}\0{row['sha256']}\n".encode()
        )
    return {
        "rows": rows,
        "tree_sha256": digest.hexdigest(),
        "files": len(rows),
        "bytes": sum(int(row["bytes"]) for row in rows),
        "exit_status": 72,
        "stderr_sha256": sha256_file(attempt / "stderr.log"),
        "canonical_outputs_present": False,
        "noncanonical_pool_checksum_verified": True,
        "truth_labels_accessed": False,
        "reuse_permitted": False,
    }


def validate_pipeline_files(code_root: Path, values: Iterable[str]) -> list[str]:
    paths = sorted(
        set(
            (
                *ENGINEERING_FILES,
                *SCIENTIFIC_FILES,
                *PIPELINE_ENTRIES.values(),
                *values,
            )
        )
    )
    if set(PIPELINE_ENTRIES.values()) - set(paths):
        raise ValueError("Execution implementation lacks a pipeline entry")
    for relative in paths:
        safe_relative_path(relative, "pipeline implementation path")
        folded = relative.casefold()
        hits = [token for token in FORBIDDEN_PIPELINE_PATH_TOKENS if token in folded]
        if hits:
            raise ValueError(f"Ranker/H2 implementation is forbidden: {relative}")
        path = code_root / relative
        if not path.is_file() or path.is_symlink() or path.stat().st_size == 0:
            raise ValueError(f"Missing, empty or symlinked implementation: {relative}")
    return paths


def safe_extract_archive(archive_path: Path, destination: Path) -> None:
    with tarfile.open(archive_path, "r:") as archive:
        members = archive.getmembers()
        for member in members:
            raw = member.name.rstrip("/")
            if not raw:
                continue
            try:
                safe_relative_path(raw, "source archive member")
            except ValueError as error:
                raise ValueError(f"Unsafe source archive member: {member.name}") from error
            if not (member.isdir() or member.isreg()) or member.issym() or member.islnk():
                raise ValueError(f"Unsafe source archive member type: {member.name}")
        archive.extractall(destination, members=members, filter="data")
    for path in destination.rglob("*"):
        if path.is_symlink():
            raise ValueError(f"Symlink appeared in frozen source: {path}")


def environment_locks(prefix: Path) -> tuple[bytes, bytes | None]:
    if (
        not prefix.is_dir()
        or prefix.is_symlink()
        or not (prefix / "conda-meta/history").is_file()
        or str(prefix).startswith("/nas_data/")
    ):
        raise ValueError(f"Unsafe or non-conda environment prefix: {prefix}")
    explicit = subprocess.run(
        ["conda", "list", "--explicit", "-p", str(prefix)],
        check=True, capture_output=True,
    ).stdout
    if not explicit.strip():
        raise ValueError(f"Empty environment lock: {prefix}")
    python = prefix / ("python.exe" if os.name == "nt" else "bin/python")
    pip: bytes | None = None
    if python.is_file() and not python.is_symlink():
        pip = subprocess.run(
            [str(python), "-m", "pip", "freeze", "--all"],
            check=True, capture_output=True,
        ).stdout
        if not pip.strip():
            raise ValueError(f"Empty pip lock for Python environment: {prefix}")
    return explicit, pip


def chmod_read_only(root: Path) -> None:
    if os.name == "nt":
        return
    for path in sorted(root.rglob("*"), reverse=True):
        path.chmod(0o550 if path.is_dir() else 0o440)
    root.chmod(0o550)


def freeze_execution(
    *,
    code_root: Path,
    protocol: Path,
    code_commit: str,
    environments: dict[str, Path],
    pipeline_files: list[str],
    output: Path,
    superseded_execution: Path | None = None,
    failed_attempt_log: Path | None = None,
    failed_attempt_root: Path | None = None,
    failed_stage: str | None = None,
    patch_reason: Path | None = None,
    allowed_changed_files: list[str] | None = None,
) -> Path:
    partial = Path(str(output) + ".partial")
    if output.exists() or output.is_symlink() or partial.exists() or partial.is_symlink():
        raise FileExistsError("Refusing to overwrite execution freeze or partial")
    protocol_manifest, contract = verify_protocol(protocol)
    if protocol_manifest.get("formal_runner_frozen") is not False:
        raise ValueError("Walnut protocol is not the pre-runner skeleton")
    if set(environments) != REQUIRED_ENVIRONMENTS:
        raise ValueError(f"Environment names must be exactly {sorted(REQUIRED_ENVIRONMENTS)}")
    if any("model" in name or "rank" in name for name in environments):
        raise ValueError("Model/ranker environments are forbidden")
    verify_git_state(code_root, code_commit)
    implementation = validate_pipeline_files(code_root, pipeline_files)
    patch_values = (
        superseded_execution,
        failed_attempt_log,
        failed_stage,
        patch_reason,
        allowed_changed_files,
    )
    patch_mode = all(value is not None and value != [] for value in patch_values)
    if any(value is not None and value != [] for value in patch_values) and not patch_mode:
        raise ValueError(
            "Walnut execution patch requires superseded execution, failed-attempt "
            "log, patch reason and exact changed-file whitelist"
        )
    if failed_attempt_root is not None and not patch_mode:
        raise ValueError("Failed-attempt root is only valid in execution patch mode")
    patch_audit: dict[str, Any] | None = None
    if patch_mode:
        assert superseded_execution is not None
        assert failed_attempt_log is not None
        assert failed_stage is not None
        assert patch_reason is not None
        assert allowed_changed_files is not None
        patch_audit = validate_execution_patch(
            code_root=code_root,
            code_commit=code_commit,
            protocol=protocol,
            superseded_execution=superseded_execution,
            failed_attempt_log=failed_attempt_log,
            failed_stage=failed_stage,
            patch_reason=patch_reason,
            allowed_changed_files=allowed_changed_files,
            failed_attempt_root=failed_attempt_root,
        )

    partial.mkdir(parents=True)
    try:
        archive = partial / "source.tar"
        subprocess.run(
            ["git", "-C", str(code_root), "archive", "--format=tar", "-o", str(archive), code_commit],
            check=True, capture_output=True,
        )
        source = partial / "source"
        source.mkdir()
        safe_extract_archive(archive, source)
        rows: list[dict[str, Any]] = []
        for relative in implementation:
            live = code_root / relative
            frozen = source / relative
            if (
                not frozen.is_file() or frozen.is_symlink()
                or frozen.stat().st_size != live.stat().st_size
                or sha256_file(frozen) != sha256_file(live)
            ):
                raise ValueError(f"Archived implementation differs: {relative}")
            rows.append({"relative_path": relative, "bytes": live.stat().st_size,
                         "sha256": sha256_file(live)})
        with (partial / "implementation_manifest.tsv").open(
            "x", encoding="utf-8", newline=""
        ) as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t",
                                    lineterminator="\n")
            writer.writeheader(); writer.writerows(rows)

        locks = partial / "environment_locks"
        locks.mkdir()
        environment_rows: list[dict[str, Any]] = []
        for name, prefix in sorted(environments.items()):
            explicit, pip = environment_locks(prefix)
            explicit_path = locks / f"{name}.explicit.txt"
            pip_path = locks / f"{name}.pip-freeze.txt"
            explicit_path.write_bytes(explicit)
            if pip is not None:
                pip_path.write_bytes(pip)
            environment_rows.append(
                {"name": name, "host_prefix": str(prefix),
                 "python_executable_present": pip is not None,
                 "explicit_lock": explicit_path.relative_to(partial).as_posix(),
                 "explicit_sha256": sha256_file(explicit_path),
                 "pip_lock": (
                     pip_path.relative_to(partial).as_posix() if pip is not None else ""
                 ),
                 "pip_sha256": sha256_file(pip_path) if pip is not None else ""}
            )
        with (partial / "environment_bindings.tsv").open(
            "x", encoding="utf-8", newline=""
        ) as handle:
            writer = csv.DictWriter(handle, fieldnames=list(environment_rows[0]),
                                    delimiter="\t", lineterminator="\n")
            writer.writeheader(); writer.writerows(environment_rows)

        if patch_audit is not None:
            prior_environment_rows = {
                row["name"]: row for row in patch_audit["prior"].get("environments", [])
            }
            for row in environment_rows:
                prior_row = prior_environment_rows.get(row["name"])
                fields = ("host_prefix", "explicit_sha256", "pip_sha256")
                if prior_row is None or any(
                    prior_row.get(field) != row.get(field) for field in fields
                ):
                    raise ValueError(
                        f"Walnut execution patch changed environment lineage: {row['name']}"
                    )

            assert failed_attempt_log is not None
            assert patch_reason is not None
            shutil.copyfile(
                failed_attempt_log, partial / "superseded_failed_attempt.log"
            )
            shutil.copyfile(patch_reason, partial / "patch_reason.md")
            failed_attempt_audit = patch_audit.get("failed_attempt_audit")
            if failed_attempt_audit is not None:
                with (partial / "superseded_failed_attempt_manifest.tsv").open(
                    "x", encoding="utf-8", newline=""
                ) as handle:
                    writer = csv.DictWriter(
                        handle,
                        fieldnames=["relative_path", "bytes", "sha256"],
                        delimiter="\t",
                        lineterminator="\n",
                    )
                    writer.writeheader()
                    writer.writerows(failed_attempt_audit["rows"])

        manifest = {
            "schema_version": EXECUTION_SCHEMA,
            "holdout_id": contract.holdout_id,
            "policy_id": contract.policy_id,
            "model_version": contract.model_version,
            "protocol_profile": "core_H1_only_no_ranker",
            "code_commit": code_commit,
            "protocol_code_commit": protocol_manifest["code_commit"],
            "freeze_stage": (
                PATCH_FREEZE_STAGE
                if patch_audit is not None
                else "post_protocol_pre_pair_pre_candidate_pre_label"
            ),
            "protocol_SHA256SUMS_sha256": sha256_file(protocol / "SHA256SUMS"),
            "contract_sha256": sha256_file(protocol / "contract.json"),
            "ranker_or_model_execution": False,
            "h2_or_topology_ranking_enabled": False,
            "network_access_in_blind_runner": False,
            "nas_data_mount_in_blind_runner": False,
            "complete_target_annotation_mount_in_blind_runner": False,
            "evaluator_only_mount_in_blind_runner": False,
            "truth_or_label_mount_in_blind_runner": False,
            "pipeline_entries": PIPELINE_ENTRIES,
            "blind_outputs": BLIND_OUTPUTS,
            "implementation_file_count": len(rows),
            "source_archive": {"path": "source.tar", "sha256": sha256_file(archive)},
            "environments": environment_rows,
        }
        if patch_audit is not None:
            assert superseded_execution is not None
            execution_patch = {
                "schema_version": "ploidypatch.walnut_execution_patch.v0.8",
                "chain_depth": patch_audit["previous_patch_depth"] + 1,
                "base_code_commit": patch_audit["base_commit"],
                "patch_code_commit": code_commit,
                "superseded_execution_SHA256SUMS_sha256": sha256_file(
                    superseded_execution / "SHA256SUMS"
                ),
                "failed_attempt_log": "superseded_failed_attempt.log",
                "failed_attempt_log_sha256": sha256_file(
                    partial / "superseded_failed_attempt.log"
                ),
                "failed_attempt_stage": patch_audit["failed_stage"],
                "patch_reason": "patch_reason.md",
                "patch_reason_sha256": sha256_file(partial / "patch_reason.md"),
                "changed_files": [
                    {"status": status, "relative_path": relative}
                    for status, relative in patch_audit["changed"]
                ],
                "evaluator_truth_construction_completed_before_patch": True,
                "structure_holdout_sentinel_failed_before_patch": (
                    patch_audit["failed_stage"] == "structure_holdout_sentinel"
                ),
                "blind_normalization_failed_before_patch": (
                    patch_audit["failed_stage"] == "blind_normalization"
                ),
                "candidate_methods_started_before_patch": (
                    patch_audit["failed_stage"]
                    in {"candidate_methods", "candidate_pools", "candidate_output_publication"}
                ),
                "candidate_methods_failed_before_patch": (
                    patch_audit["failed_stage"] == "candidate_methods"
                ),
                "candidate_pools_started_before_patch": (
                    patch_audit["failed_stage"]
                    in {"candidate_pools", "candidate_output_publication"}
                ),
                "candidate_pools_failed_before_patch": (
                    patch_audit["failed_stage"] == "candidate_pools"
                ),
                "candidate_output_publication_failed_before_patch": (
                    patch_audit["failed_stage"] == "candidate_output_publication"
                ),
                "blind_candidate_generation_completed_before_patch": (
                    patch_audit["failed_stage"] == "candidate_output_publication"
                ),
                "formal_scores_generated_before_patch": False,
                "truth_labels_accessed_before_patch": False,
                "scientific_thresholds_or_references_changed": False,
            }
            failed_attempt_audit = patch_audit.get("failed_attempt_audit")
            if failed_attempt_audit is not None:
                execution_patch["failed_attempt_tree"] = {
                    key: value
                    for key, value in failed_attempt_audit.items()
                    if key != "rows"
                }
                execution_patch["failed_attempt_tree_manifest"] = (
                    "superseded_failed_attempt_manifest.tsv"
                )
            manifest["execution_patch"] = execution_patch
        (partial / "execution_manifest.json").write_text(
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
    parser.add_argument("--code-root", required=True)
    parser.add_argument("--protocol-freeze", required=True)
    parser.add_argument("--code-commit", required=True)
    parser.add_argument("--environment", action="append", type=parse_environment, required=True)
    parser.add_argument("--pipeline-file", action="append", default=[])
    parser.add_argument("--superseded-execution")
    parser.add_argument("--failed-attempt-log")
    parser.add_argument("--failed-attempt-root")
    parser.add_argument("--failed-stage", choices=sorted(FAILED_STAGES))
    parser.add_argument("--patch-reason-file")
    parser.add_argument("--allow-changed-file", action="append")
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args(argv)
    environments = dict(args.environment)
    if len(environments) != len(args.environment):
        raise ValueError("Duplicate environment name")
    freeze_execution(
        code_root=Path(args.code_root).resolve(),
        protocol=Path(args.protocol_freeze).resolve(),
        code_commit=args.code_commit,
        environments=environments,
        pipeline_files=args.pipeline_file,
        output=Path(args.output_dir).resolve(),
        superseded_execution=(
            Path(args.superseded_execution).resolve()
            if args.superseded_execution else None
        ),
        failed_attempt_log=(
            Path(args.failed_attempt_log).resolve() if args.failed_attempt_log else None
        ),
        failed_attempt_root=(
            Path(args.failed_attempt_root).resolve() if args.failed_attempt_root else None
        ),
        failed_stage=args.failed_stage,
        patch_reason=(
            Path(args.patch_reason_file).resolve() if args.patch_reason_file else None
        ),
        allowed_changed_files=args.allow_changed_file,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
