#!/usr/bin/env python3
"""Freeze committed code and seven exact environments for a v0.5 holdout."""
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

from ploidypatch.artifact_manifest import (
    sha256_file,
    verify_sha256sums,
    write_sha256sums,
)
from ploidypatch.holdout_contract import load_holdout_contract, safe_relative_path


SCHEMA_VERSION = "ploidypatch.external_holdout_execution_freeze.v0.5"
PROTOCOL_SCHEMA = "ploidypatch.external_holdout_protocol_freeze.v0.5"
MODEL_SCHEMA = "ploidypatch.composite_ranker.v0.4"
PATCH_STAGE = (
    "post_evaluator_truth_failed_blind_pre_candidate_pre_score_pre_label_execution_patch"
)
REQUIRED_FAILED_ATTEMPT_FILES = frozenset(
    {
        "exit_status.txt",
        "stdout.log",
        "stderr.log",
        "bwrap_command.txt",
        "mount_manifest.json",
        "namespace_role_validation.json",
    }
)
REQUIRED_ENVIRONMENTS = frozenset(
    {
        "ploidypatch-dev",
        "ploidypatch-model",
        "ploidypatch-baseline",
        "ploidypatch-synteny",
        "ploidypatch-syngap",
        "ploidypatch-gemoma",
        "ploidypatch-lifton",
    }
)
_FULL_SHA = re.compile(r"[0-9a-f]{40}")
_ENVIRONMENT_NAME = re.compile(r"[a-z0-9][a-z0-9_.-]*")


def _reject_duplicate_keys(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    value: dict[str, Any] = {}
    for key, item in pairs:
        if key in value:
            raise ValueError(f"Duplicate JSON key: {key}")
        value[key] = item
    return value


def load_json(path: Path) -> dict[str, Any]:
    if not path.is_file() or path.is_symlink() or path.stat().st_size == 0:
        raise ValueError(f"Missing, empty or symlinked JSON input: {path}")
    value = json.loads(
        path.read_text(encoding="utf-8"), object_pairs_hook=_reject_duplicate_keys
    )
    if not isinstance(value, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return value


def parse_environment(value: str) -> tuple[str, Path]:
    name, separator, prefix = value.partition("=")
    if not separator or _ENVIRONMENT_NAME.fullmatch(name) is None or not prefix:
        raise argparse.ArgumentTypeError("--environment must be NAME=PREFIX")
    return name, Path(prefix).resolve()


def run_checked(command: list[str], *, cwd: Path | None = None) -> bytes:
    return subprocess.run(
        command,
        cwd=cwd,
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    ).stdout


def verify_git_state(code_root: Path, code_commit: str) -> None:
    if _FULL_SHA.fullmatch(code_commit) is None:
        raise ValueError("code_commit must be a full lowercase Git SHA")
    head = run_checked(["git", "-C", str(code_root), "rev-parse", "HEAD"]).decode().strip()
    if head != code_commit:
        raise ValueError(f"code_commit differs from code-root HEAD: {head}")
    if run_checked(["git", "-C", str(code_root), "status", "--porcelain"]):
        raise ValueError("Code root must be completely clean, including untracked files")


def create_git_archive(code_root: Path, code_commit: str, output: Path) -> None:
    with output.open("xb") as handle:
        subprocess.run(
            ["git", "-C", str(code_root), "archive", "--format=tar", code_commit],
            check=True,
            stdout=handle,
        )


def archive_names(archive: Path) -> set[str]:
    with tarfile.open(archive, "r:*") as handle:
        return {member.name.rstrip("/") for member in handle.getmembers()}


def extract_archive_safely(archive: Path, destination: Path) -> None:
    with tarfile.open(archive, "r:*") as handle:
        members = handle.getmembers()
        for member in members:
            relative = PurePosixPath(member.name)
            if (
                relative.is_absolute()
                or ".." in relative.parts
                or member.issym()
                or member.islnk()
                or not (member.isdir() or member.isfile())
            ):
                raise ValueError(f"Unsafe member in source archive: {member.name}")
        handle.extractall(destination, members=members)


def environment_lock(prefix: Path) -> bytes:
    if prefix.is_symlink() or not (prefix / "conda-meta/history").is_file():
        raise ValueError(f"Not a real conda environment prefix: {prefix}")
    return run_checked(["conda", "list", "--explicit", "-p", str(prefix)])


def pip_lock(prefix: Path) -> bytes:
    python = prefix / "bin/python"
    if not python.is_file():
        return b"# python unavailable; pip lock not applicable\n"
    try:
        return run_checked([str(python), "-m", "pip", "freeze", "--all"])
    except subprocess.CalledProcessError:
        return b"# pip unavailable; explicit conda lock is authoritative\n"


def read_manifest_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ["relative_path", "bytes", "sha256"]:
            raise ValueError(f"Implementation manifest fields differ: {path}")
        rows = list(reader)
    if not rows or len({row["relative_path"] for row in rows}) != len(rows):
        raise ValueError("Implementation manifest is empty or has duplicate paths")
    return rows


def verify_implementation(code_root: Path, rows: list[dict[str, str]]) -> None:
    for row in rows:
        relative = safe_relative_path(row["relative_path"], "implementation path")
        path = code_root.joinpath(*relative.parts)
        if (
            not path.is_file()
            or path.is_symlink()
            or str(path.stat().st_size) != row["bytes"]
            or sha256_file(path) != row["sha256"]
        ):
            raise ValueError(f"Implementation differs from protocol: {relative}")


def git_changed_files(
    code_root: Path, base_commit: str, patch_commit: str
) -> list[tuple[str, str]]:
    raw = run_checked(
        [
            "git", "-C", str(code_root), "diff", "--name-status",
            "--no-renames", base_commit, patch_commit, "--",
        ]
    ).decode("utf-8")
    rows: list[tuple[str, str]] = []
    for line in raw.splitlines():
        status, separator, path = line.partition("\t")
        if not separator or status not in {"A", "M"}:
            raise ValueError(f"Execution patch may only add or modify files: {line}")
        relative = safe_relative_path(path, "changed implementation path").as_posix()
        rows.append((status, relative))
    if not rows or len({path for _, path in rows}) != len(rows):
        raise ValueError("Execution patch diff is empty or duplicate")
    return rows


def verify_changed_whitelist(
    changed: list[tuple[str, str]], allowed: Iterable[str], patch_reason: str
) -> None:
    allowed_paths = {
        safe_relative_path(path, "allowed changed file").as_posix() for path in allowed
    }
    changed_paths = {path for _, path in changed}
    if changed_paths != allowed_paths or patch_reason not in allowed_paths:
        raise ValueError(
            "Patch changed files must exactly equal the explicit whitelist and include "
            "the patch-reason file"
        )
    for path in changed_paths:
        if path == patch_reason:
            if not path.startswith("docs/"):
                raise ValueError("Patch reason must be one committed docs/ file")
        elif not path.startswith(("scripts/", "src/ploidypatch/", "tests/")):
            raise ValueError(f"Non-implementation patch file is forbidden: {path}")


def _tree_digest(rows: list[dict[str, Any]]) -> str:
    digest = hashlib.sha256()
    for row in rows:
        digest.update(
            f"{row['relative_path']}\0{row['bytes']}\0{row['sha256']}\n".encode()
        )
    return digest.hexdigest()


def freeze_failed_attempt_manifest(
    attempt: Path, output: Path, blind_outputs: dict[str, str]
) -> dict[str, Any]:
    if not attempt.is_dir() or attempt.is_symlink():
        raise ValueError("Failed blind attempt must be a real directory")
    rows: list[dict[str, Any]] = []
    for path in sorted(attempt.rglob("*")):
        if path.is_symlink():
            raise ValueError("Failed blind attempt may not contain symlinks")
        if path.is_file():
            rows.append(
                {
                    "relative_path": path.relative_to(attempt).as_posix(),
                    "bytes": path.stat().st_size,
                    "sha256": sha256_file(path),
                }
            )
    paths = {row["relative_path"] for row in rows}
    if not REQUIRED_FAILED_ATTEMPT_FILES <= paths:
        raise ValueError("Failed blind attempt lacks required log/isolation evidence")
    raw_status = (attempt / "exit_status.txt").read_text(encoding="utf-8").strip()
    if not raw_status.isdigit() or int(raw_status) == 0:
        raise ValueError("Execution patch requires a retained nonzero blind attempt")
    command = (attempt / "bwrap_command.txt").read_text(encoding="utf-8")
    if any(flag not in command for flag in ("--unshare-all", "--unshare-net", "--clearenv")):
        raise ValueError("Failed attempt lacks frozen bubblewrap flags")
    mount = load_json(attempt / "mount_manifest.json")
    namespace = load_json(attempt / "namespace_role_validation.json")
    data_roles = {
        item.get("role")
        for item in mount.get("mounts", [])
        if isinstance(item, dict)
        and item.get("role") in {"shared_target", "candidate_only", "blind_benchmark"}
    }
    if (
        mount.get("schema_version") != "ploidypatch.blind_mount_manifest.v0.5"
        or data_roles != {"shared_target", "candidate_only", "blind_benchmark"}
        or namespace.get("schema_version")
        != "ploidypatch.blind_namespace_validation.v0.5"
        or namespace.get("shared_target_visible") is not True
        or namespace.get("candidate_only_visible") is not True
        or namespace.get("blind_benchmark_visible") is not True
        or namespace.get("evaluator_only_visible") is not False
        or namespace.get("truth_visible") is not False
        or namespace.get("complete_target_annotation_visible") is not False
        or namespace.get("nas_data_visible") is not False
    ):
        raise ValueError("Failed attempt isolation evidence is incomplete or violated")
    for item in mount["mounts"]:
        text = f"{item.get('host_path', '')}\n{item.get('namespace_path', '')}".lower()
        if any(
            token in text
            for token in ("/nas_data", "evaluator_only", "/truth", "/labels", "target_complete")
        ):
            raise ValueError("Failed attempt mounted evaluator/truth/labels/NAS")
    for name, raw in blind_outputs.items():
        relative = safe_relative_path(raw, f"blind output {name}").as_posix()
        if name != "command_log" and f"project/{relative}" in paths:
            raise ValueError(f"Failed attempt already contains formal blind output: {name}")
    forbidden: list[str] = []
    for relative in paths:
        folded = relative.casefold()
        name = PurePosixPath(relative).name.casefold()
        if not folded.startswith("project/"):
            continue
        if (
            name in {"candidate.gff3", "decisions.tsv", "candidate_labels.tsv"}
            or "/scores/" in folded
            or "/labels/" in folded
            or "/truth/" in folded
            or name == "custody_manifest.json"
        ):
            forbidden.append(relative)
    if forbidden:
        raise ValueError(
            "Failed attempt advanced beyond pre-candidate/pre-score/pre-label: "
            + ", ".join(sorted(forbidden))
        )
    with output.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["relative_path", "bytes", "sha256"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)
    return {
        "exit_status": int(raw_status),
        "tree_sha256": _tree_digest(rows),
        "files": len(rows),
        "bytes": sum(int(row["bytes"]) for row in rows),
    }


def patch_implementation(
    code_root: Path,
    base_rows: list[dict[str, str]],
    changed: list[tuple[str, str]],
) -> list[dict[str, str]]:
    changed_paths = {path for _, path in changed}
    base_paths = {row["relative_path"] for row in base_rows}
    rows: list[dict[str, str]] = []
    for relative in sorted(base_paths | changed_paths):
        path = code_root / relative
        if not path.is_file() or path.is_symlink() or path.stat().st_size == 0:
            raise ValueError(f"Patched implementation is missing: {relative}")
        if relative in base_paths and relative not in changed_paths:
            base = next(row for row in base_rows if row["relative_path"] == relative)
            if str(path.stat().st_size) != base["bytes"] or sha256_file(path) != base["sha256"]:
                raise ValueError(f"Non-whitelisted implementation changed: {relative}")
        rows.append(
            {
                "relative_path": relative,
                "bytes": str(path.stat().st_size),
                "sha256": sha256_file(path),
            }
        )
    return rows


def reject_target_artifacts(project_root: Path, manifest: dict[str, Any]) -> None:
    relatives = manifest.get("forbidden_target_artifacts")
    if not isinstance(relatives, list) or not relatives:
        raise ValueError("Protocol lacks forbidden target-artifact paths")
    for raw in relatives:
        relative = safe_relative_path(raw, "forbidden target-artifact path")
        path = project_root.joinpath(*relative.parts)
        if path.exists() or path.is_symlink():
            raise ValueError(
                "Execution freeze is too late; target-derived artifact exists: "
                f"{relative.as_posix()}"
            )


def chmod_read_only(root: Path) -> None:
    for path in sorted(root.rglob("*"), reverse=True):
        path.chmod(0o550 if path.is_dir() else 0o440)
    root.chmod(0o550)


def freeze_execution(
    *,
    project_root: Path,
    code_root: Path,
    protocol_root: Path,
    model_root: Path,
    environments: dict[str, Path],
    code_commit: str,
    output: Path,
    supersedes_execution: Path | None = None,
    failed_attempt: Path | None = None,
    patch_reason: Path | None = None,
    staged_inputs: Path | None = None,
    allowed_changed_files: list[str] | None = None,
) -> Path:
    partial = Path(str(output) + ".partial")
    if output.exists() or output.is_symlink() or partial.exists() or partial.is_symlink():
        raise FileExistsError("Refusing to overwrite execution freeze or partial output")
    if set(environments) != REQUIRED_ENVIRONMENTS:
        raise ValueError(f"Environment names must be exactly {sorted(REQUIRED_ENVIRONMENTS)}")
    patch_values = (
        supersedes_execution,
        failed_attempt,
        patch_reason,
        staged_inputs,
        allowed_changed_files,
    )
    patch_mode = all(value is not None and value != [] for value in patch_values)
    if any(value is not None and value != [] for value in patch_values) and not patch_mode:
        raise ValueError(
            "Execution patch requires superseded execution, failed attempt, patch "
            "reason, staged inputs and changed-file whitelist together"
        )
    verify_git_state(code_root, code_commit)
    verify_sha256sums(protocol_root, ignore_checksum_file=True)
    verify_sha256sums(model_root, ignore_checksum_file=True)
    protocol = load_json(protocol_root / "protocol_manifest.json")
    prior: dict[str, Any] | None = None
    changed: list[tuple[str, str]] = []
    patch_reason_relative: str | None = None
    if patch_mode:
        assert supersedes_execution is not None
        assert failed_attempt is not None
        assert patch_reason is not None
        assert staged_inputs is not None
        assert allowed_changed_files is not None
        verify_sha256sums(supersedes_execution, ignore_checksum_file=True)
        verify_sha256sums(staged_inputs, ignore_checksum_file=True)
        prior = load_json(supersedes_execution / "execution_manifest.json")
        if (
            prior.get("schema_version") != SCHEMA_VERSION
            or prior.get("execution_patch") is not None
            or prior.get("protocol_SHA256SUMS_sha256")
            != sha256_file(protocol_root / "SHA256SUMS")
            or prior.get("contract_sha256") != sha256_file(protocol_root / "contract.json")
            or prior.get("composite_model_SHA256SUMS_sha256")
            != sha256_file(model_root / "SHA256SUMS")
        ):
            raise ValueError("Superseded execution is not the exact base freeze")
        if (
            protocol.get("staged_input_SHA256SUMS_sha256")
            != sha256_file(staged_inputs / "SHA256SUMS")
            or (protocol_root / "role_manifest.tsv").read_bytes()
            != (staged_inputs / "role_manifest.tsv").read_bytes()
            or (protocol_root / "role_contract.json").read_bytes()
            != (staged_inputs / "role_contract.json").read_bytes()
        ):
            raise ValueError("Patch staged inputs differ from the base protocol")
        try:
            patch_reason_relative = patch_reason.resolve().relative_to(
                code_root.resolve()
            ).as_posix()
        except ValueError:
            raise ValueError("Patch reason must be committed inside code root") from None
        if not patch_reason.is_file() or patch_reason.is_symlink() or patch_reason.stat().st_size == 0:
            raise ValueError("Patch reason is missing, empty or symlinked")
        changed = git_changed_files(code_root, str(prior["code_commit"]), code_commit)
        verify_changed_whitelist(changed, allowed_changed_files, patch_reason_relative)
    if (
        protocol.get("schema_version") != PROTOCOL_SCHEMA
        or protocol.get("code_commit")
        != (prior.get("code_commit") if patch_mode and prior is not None else code_commit)
        or protocol.get("freeze_stage")
        != "post_metadata_pre_pair_pre_candidate_pre_label"
        or protocol.get("truth_access") is not False
        or protocol.get("wgd_pairs_enumerated") is not False
        or protocol.get("candidate_counts_computed") is not False
        or protocol.get("truth_labels_accessed") is not False
        or protocol.get("composite_model_SHA256SUMS_sha256")
        != sha256_file(model_root / "SHA256SUMS")
    ):
        raise ValueError("Protocol freeze fails pre-enumeration or model binding")
    contract = load_holdout_contract(protocol_root / "contract.json")
    if (
        contract.holdout_id != protocol.get("holdout_id")
        or contract.policy_id != protocol.get("policy_id")
        or contract.model_version != protocol.get("model_version")
    ):
        raise ValueError("Protocol manifest and frozen contract differ")
    model_manifest = load_json(model_root / "composite_manifest.json")
    if (
        model_manifest.get("schema_version") != MODEL_SCHEMA
        or model_manifest.get("automatic_approval") is not False
    ):
        raise ValueError("Execution requires the exact safe v0.4 composite model")
    if not patch_mode:
        reject_target_artifacts(project_root, protocol)
    implementation = read_manifest_tsv(protocol_root / "implementation_manifest.tsv")
    if patch_mode:
        implementation = patch_implementation(code_root, implementation, changed)
    else:
        verify_implementation(code_root, implementation)

    partial.mkdir(parents=True)
    try:
        archive = partial / "source.tar"
        create_git_archive(code_root, code_commit, archive)
        names = archive_names(archive)
        absent = [row["relative_path"] for row in implementation if row["relative_path"] not in names]
        if absent:
            raise ValueError(
                "Frozen implementation is absent from committed archive: "
                + ", ".join(absent)
            )
        source = partial / "source"
        source.mkdir()
        extract_archive_safely(archive, source)
        for row in implementation:
            path = source / row["relative_path"]
            if (
                not path.is_file()
                or str(path.stat().st_size) != row["bytes"]
                or sha256_file(path) != row["sha256"]
            ):
                raise ValueError(
                    "Archived implementation differs: " + row["relative_path"]
                )

        locks = partial / "environment_locks"
        locks.mkdir()
        environment_rows: list[dict[str, Any]] = []
        for name, prefix in sorted(environments.items()):
            explicit_path = locks / f"{name}.explicit.txt"
            pip_path = locks / f"{name}.pip-freeze.txt"
            explicit_path.write_bytes(environment_lock(prefix))
            pip_path.write_bytes(pip_lock(prefix))
            environment_rows.append(
                {
                    "name": name,
                    "host_prefix": str(prefix),
                    "explicit_lock": explicit_path.relative_to(partial).as_posix(),
                    "explicit_sha256": sha256_file(explicit_path),
                    "pip_lock": pip_path.relative_to(partial).as_posix(),
                    "pip_sha256": sha256_file(pip_path),
                }
            )
        if patch_mode:
            assert prior is not None
            if environment_rows != prior.get("environments"):
                raise ValueError(
                    "Execution patch changed an environment prefix or conda/pip lock"
                )
        with (partial / "environment_bindings.tsv").open(
            "x", encoding="utf-8", newline=""
        ) as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=list(environment_rows[0]),
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerows(environment_rows)
        with (partial / "implementation_manifest.tsv").open(
            "x", encoding="utf-8", newline=""
        ) as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=["relative_path", "bytes", "sha256"],
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerows(implementation)
        for name in ("pipeline_entries.tsv", "pipeline_files.tsv", "blind_outputs.tsv"):
            shutil.copyfile(protocol_root / name, partial / name)

        failed_summary: dict[str, Any] | None = None
        changed_rows: list[dict[str, str]] = []
        if patch_mode:
            assert failed_attempt is not None
            assert prior is not None
            assert patch_reason is not None
            assert patch_reason_relative is not None
            failed_summary = freeze_failed_attempt_manifest(
                failed_attempt,
                partial / "superseded_failed_attempt_manifest.tsv",
                prior["blind_outputs"],
            )
            evidence = partial / "superseded_failure_evidence"
            evidence.mkdir()
            for relative in sorted(REQUIRED_FAILED_ATTEMPT_FILES):
                shutil.copyfile(failed_attempt / relative, evidence / relative)
            for status, relative in changed:
                try:
                    base_blob = run_checked(
                        [
                            "git", "-C", str(code_root), "show",
                            f"{prior['code_commit']}:{relative}",
                        ]
                    )
                    base_sha = hashlib.sha256(base_blob).hexdigest()
                except subprocess.CalledProcessError:
                    base_sha = ""
                changed_rows.append(
                    {
                        "status": status,
                        "relative_path": relative,
                        "base_sha256": base_sha,
                        "patch_sha256": sha256_file(code_root / relative),
                    }
                )
            with (partial / "changed_files.tsv").open(
                "x", encoding="utf-8", newline=""
            ) as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=list(changed_rows[0]),
                    delimiter="\t",
                    lineterminator="\n",
                )
                writer.writeheader()
                writer.writerows(changed_rows)
            shutil.copyfile(patch_reason, partial / "patch_reason.md")

        execution = {
            "schema_version": SCHEMA_VERSION,
            "holdout_id": contract.holdout_id,
            "policy_id": contract.policy_id,
            "test_role": contract.test_role,
            "model_version": contract.model_version,
            "code_commit": code_commit,
            "freeze_stage": PATCH_STAGE if patch_mode else "post_metadata_pre_pair_pre_candidate_pre_label",
            "created_before": {
                "wgd_pair_enumeration": not patch_mode,
                "candidate_generation": True,
                "candidate_labels": True,
                "candidate_scores": True,
            },
            "network_access_in_blind_runner": False,
            "nas_data_mount_in_blind_runner": False,
            "complete_target_annotation_mount_in_blind_runner": False,
            "evaluator_only_mount_in_blind_runner": False,
            "truth_or_label_mount_in_blind_runner": False,
            "protocol_SHA256SUMS_sha256": sha256_file(protocol_root / "SHA256SUMS"),
            "contract_sha256": sha256_file(protocol_root / "contract.json"),
            "composite_model_SHA256SUMS_sha256": sha256_file(
                model_root / "SHA256SUMS"
            ),
            "source_archive": {
                "format": "git_archive_tar",
                "path": "source.tar",
                "sha256": sha256_file(archive),
            },
            "implementation_file_count": len(implementation),
            "pipeline_entries": protocol["pipeline_entries"],
            "blind_outputs": protocol["blind_outputs"],
            "forbidden_target_artifacts": protocol["forbidden_target_artifacts"],
            "environments": environment_rows,
        }
        if patch_mode:
            assert prior is not None
            assert failed_summary is not None
            assert patch_reason_relative is not None
            execution["execution_patch"] = {
                "schema_version": "ploidypatch.external_holdout_execution_patch.v0.5",
                "freeze_stage": PATCH_STAGE,
                "base_code_commit": prior["code_commit"],
                "patch_code_commit": code_commit,
                "superseded_execution_SHA256SUMS_sha256": sha256_file(
                    supersedes_execution / "SHA256SUMS"  # type: ignore[operator]
                ),
                "base_protocol_SHA256SUMS_sha256": sha256_file(
                    protocol_root / "SHA256SUMS"
                ),
                "contract_sha256": sha256_file(protocol_root / "contract.json"),
                "staged_input_SHA256SUMS_sha256": sha256_file(
                    staged_inputs / "SHA256SUMS"  # type: ignore[operator]
                ),
                "composite_model_SHA256SUMS_sha256": sha256_file(
                    model_root / "SHA256SUMS"
                ),
                "failed_attempt": failed_summary,
                "failed_attempt_manifest": "superseded_failed_attempt_manifest.tsv",
                "patch_reason_source_path": patch_reason_relative,
                "patch_reason_sha256": sha256_file(partial / "patch_reason.md"),
                "changed_files": changed_rows,
                "scientific_protocol_changed": False,
                "contract_or_policy_changed": False,
                "model_or_threshold_changed": False,
                "staged_inputs_changed": False,
                "truth_or_benchmark_regenerated": False,
                "evaluator_truth_construction_completed_before_patch": True,
                "blind_candidate_wgd_completed_before_patch": False,
                "candidate_generation_completed_before_patch": False,
                "formal_scores_generated_before_patch": False,
                "truth_labels_accessed_before_patch": False,
                "automatic_approval": False,
            }
        (partial / "execution_manifest.json").write_text(
            json.dumps(execution, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
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
    parser.add_argument("--protocol-freeze", required=True)
    parser.add_argument("--composite-model-freeze", required=True)
    parser.add_argument("--environment", action="append", type=parse_environment, required=True)
    parser.add_argument("--code-commit", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--supersedes-execution-freeze")
    parser.add_argument("--failed-blind-attempt")
    parser.add_argument("--patch-reason-file")
    parser.add_argument("--staged-inputs")
    parser.add_argument("--allow-changed-file", action="append")
    args = parser.parse_args(argv)
    environments = dict(args.environment)
    if len(environments) != len(args.environment):
        raise ValueError("Duplicate environment name")
    freeze_execution(
        project_root=Path(args.project_root).resolve(),
        code_root=Path(args.code_root).resolve(),
        protocol_root=Path(args.protocol_freeze).resolve(),
        model_root=Path(args.composite_model_freeze).resolve(),
        environments=environments,
        code_commit=args.code_commit,
        output=Path(args.output_dir).resolve(),
        supersedes_execution=(
            Path(args.supersedes_execution_freeze).resolve()
            if args.supersedes_execution_freeze
            else None
        ),
        failed_attempt=(
            Path(args.failed_blind_attempt).resolve()
            if args.failed_blind_attempt
            else None
        ),
        patch_reason=(
            Path(args.patch_reason_file).resolve() if args.patch_reason_file else None
        ),
        staged_inputs=(
            Path(args.staged_inputs).resolve() if args.staged_inputs else None
        ),
        allowed_changed_files=args.allow_changed_file,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
