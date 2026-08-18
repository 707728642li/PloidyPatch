#!/usr/bin/env python3
"""Create the immutable, pre-enumeration Populus v0.4 execution freeze.

The protocol freeze fixes the scientific question.  This second freeze fixes the
bytes that are allowed to implement it.  It deliberately refuses to run after
any Populus pair, benchmark, candidate, score, or label output exists.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import shutil
import subprocess
import tarfile
from pathlib import Path, PurePosixPath
from typing import Any, Iterable


SCHEMA_VERSION = "ploidypatch.populus_execution_freeze.v0.4"
PROTOCOL_SCHEMA = "ploidypatch.populus_external_protocol_freeze.v0.4"

# These names are part of the execution contract.  A freeze may not silently
# omit a stage merely because its implementation has not been written yet.
REQUIRED_PIPELINE_FILES = (
    "scripts/freeze_populus_execution_implementation_v0.4.py",
    "scripts/run_populus_blind_isolated_v0.4.sh",
    "scripts/run_populus_blind_pipeline_v0.4.sh",
    "scripts/run_populus_external_reveal_v0.4.sh",
    "scripts/finalize_populus_blind_custody_v0.4.py",
    "scripts/freeze_populus_external_protocol_v0.4.py",
    "scripts/preflight_external_inputs_v0.4.py",
    "scripts/stage_populus_external_inputs_v0.4.py",
    "scripts/audit_copy_pair_selection_truth.py",
    "scripts/build_wgdi_source_alias_gff.py",
    "scripts/prepare_populus_external_normalized_inputs_v0.4.sh",
    "scripts/prepare_populus_evaluator_wgdi_inputs_v0.4.sh",
    "scripts/run_populus_evaluator_wgdi_v0.4.sh",
    "scripts/infer_populus_external_pairs_v0.4.sh",
    "scripts/run_populus_copy_collapse_benchmark_v0.4.sh",
    "scripts/run_populus_miniprot_upstream_v0.4.sh",
    "scripts/synthesize_missing_transcript_exons.py",
    "scripts/run_gemoma_homology.sh",
    "scripts/publish_gemoma_working.sh",
    "scripts/run_lifton_transfer.sh",
    "scripts/build_populus_method_trio_candidate_pools_v0.4.sh",
    "scripts/run_populus_blind_union_self_wgd_v0.4.sh",
    "scripts/score_populus_candidates_blind_v0.4.sh",
    "scripts/build_populus_complete_control_reveal_inputs_v0.4.sh",
    "scripts/evaluate_external_v0.4.py",
    "scripts/evaluate_apple_external_v0.3.py",
    "docs/POPULUS_EXECUTION_PATCH_v0.4.1.md",
    "docs/POPULUS_EXECUTION_PATCH_v0.4.2.md",
)

# Presence of any of these outputs proves that execution-relevant target
# structure has already been enumerated.  Staged role-separated input bytes are
# intentionally allowed; they were selected by metadata-only preflight.
FORBIDDEN_PRE_FREEZE_PATHS = (
    "data/derived/external_inputs/populus_v0.4",
    "data/derived/external_evaluator/populus_v0.4_wgdi_inputs",
    "results/evaluator/populus/v0.4/wgdi",
    "results/evaluator/populus/v0.4/truth_pairs",
    "benchmark/structure/copy_collapse_v0.4/ptr_v4/annotation_copy_collapse_seed20260930",
    "results/baselines/populus_v0.4",
    "results/copy_collapse/external/populus_v0.4_method_trio",
    "results/copy_collapse/external/populus_v0.4_blind_self_wgd",
    "results/copy_collapse/external/populus_v0.4_blind_rankings",
    "results/copy_collapse/external/populus_v0.4_reveal",
    "results/blind_runs/populus_external_v0.4",
    "work/populus_external_v0.4",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def verify_sha256sums(root: Path) -> None:
    checksum = root / "SHA256SUMS"
    if not checksum.is_file() or checksum.stat().st_size == 0:
        raise ValueError(f"Missing SHA256SUMS: {root}")
    for number, line in enumerate(checksum.read_text(encoding="utf-8").splitlines(), 1):
        digest, separator, relative = line.partition("  ")
        path = root / relative
        if (
            not separator
            or not re.fullmatch(r"[0-9a-f]{64}", digest)
            or not path.is_file()
            or sha256(path) != digest
        ):
            raise ValueError(f"Checksum failure at line {number}: {root}")


def parse_environment(value: str) -> tuple[str, Path]:
    name, separator, prefix = value.partition("=")
    if (
        not separator
        or not re.fullmatch(r"[a-z0-9][a-z0-9_.-]*", name)
        or not prefix
    ):
        raise argparse.ArgumentTypeError("--environment must be NAME=PREFIX")
    return name, Path(prefix).resolve()


def run_checked(command: list[str], *, cwd: Path | None = None) -> bytes:
    return subprocess.run(
        command, cwd=cwd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    ).stdout


def verify_git_state(code_root: Path, commit: str) -> None:
    head = run_checked(["git", "-C", str(code_root), "rev-parse", "HEAD"]).decode().strip()
    if head != commit:
        raise ValueError(f"--code-commit {commit} is not code-root HEAD {head}")
    dirty = run_checked(["git", "-C", str(code_root), "status", "--porcelain"])
    if dirty:
        raise ValueError("Code root has tracked or untracked changes; commit before freeze")


def create_git_archive(code_root: Path, commit: str, output: Path) -> None:
    with output.open("xb") as handle:
        subprocess.run(
            ["git", "-C", str(code_root), "archive", "--format=tar", commit],
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
            pure = PurePosixPath(member.name)
            if (
                pure.is_absolute()
                or ".." in pure.parts
                or member.issym()
                or member.islnk()
                or not (member.isdir() or member.isfile())
            ):
                raise ValueError(f"Unsafe member in source archive: {member.name}")
        handle.extractall(destination, members=members)


def environment_lock(prefix: Path) -> bytes:
    if not (prefix / "conda-meta/history").is_file():
        raise ValueError(f"Not a conda environment prefix: {prefix}")
    return run_checked(["conda", "list", "--explicit", "-p", str(prefix)])


def pip_lock(prefix: Path) -> bytes:
    python = prefix / "bin/python"
    if not python.is_file():
        return b"# python unavailable; pip lock not applicable\n"
    try:
        return run_checked([str(python), "-m", "pip", "freeze", "--all"])
    except subprocess.CalledProcessError:
        return b"# pip unavailable; explicit conda lock is authoritative\n"


def implementation_files(code_root: Path) -> list[str]:
    core = [
        path.relative_to(code_root).as_posix()
        for path in sorted((code_root / "src/ploidypatch").glob("*.py"))
    ]
    configured = [
        path.relative_to(code_root).as_posix()
        for pattern in (
            "config/*populus*v0.4*",
            "config/primary_seqids/*.tsv",
            "docs/POPULUS_EXTERNAL_VALIDATION_PROTOCOL_v0.4.md",
        )
        for path in sorted(code_root.glob(pattern))
        if path.is_file()
    ]
    return sorted(set((*REQUIRED_PIPELINE_FILES, *core, *configured)))


def make_checksums(root: Path) -> None:
    checksum = root / "SHA256SUMS"
    with checksum.open("x", encoding="utf-8", newline="") as handle:
        for path in sorted(item for item in root.rglob("*") if item.is_file()):
            if path != checksum:
                handle.write(f"{sha256(path)}  {path.relative_to(root).as_posix()}\n")


def chmod_read_only(root: Path) -> None:
    for path in sorted(root.rglob("*"), reverse=True):
        path.chmod(0o550 if path.is_dir() else 0o440)
    root.chmod(0o550)


def freeze_failed_attempt_manifest(attempt: Path, output: Path) -> int:
    if not attempt.is_dir() or attempt.is_symlink():
        raise ValueError("Failed blind attempt must be a real directory")
    attempt_resolved = attempt.resolve()
    symlinks = sorted(path for path in attempt.rglob("*") if path.is_symlink())
    for path in symlinks:
        target = os.readlink(path)
        if os.path.isabs(target) or not path.resolve().is_relative_to(attempt_resolved):
            raise ValueError("Failed blind attempt contains an external symlink")
    exit_status = attempt / "exit_status.txt"
    if not exit_status.is_file():
        raise ValueError("Failed blind attempt lacks exit_status.txt")
    status_text = exit_status.read_text(encoding="utf-8").strip()
    if not status_text.isdigit() or int(status_text) == 0:
        raise ValueError("Patch freeze requires a retained nonzero blind attempt")
    forbidden = (
        "custody_manifest.json",
        "reveal_authorization.json",
        "v04.tsv",
        "candidate_labels.tsv",
    )
    files = sorted(
        path for path in attempt.rglob("*") if path.is_file() and not path.is_symlink()
    )
    if any(path.name in forbidden for path in files):
        raise ValueError("Failed attempt advanced beyond the pre-score/pre-reveal barrier")
    with output.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("relative_path", "kind", "bytes", "sha256"))
        for path in files:
            writer.writerow(
                (
                    path.relative_to(attempt).as_posix(),
                    "file",
                    path.stat().st_size,
                    sha256(path),
                )
            )
        for path in symlinks:
            target = os.readlink(path)
            encoded = f"symlink\0{target}".encode("utf-8")
            writer.writerow(
                (
                    path.relative_to(attempt).as_posix(),
                    "internal_relative_symlink",
                    len(target.encode("utf-8")),
                    hashlib.sha256(encoded).hexdigest(),
                )
            )
    return int(status_text)


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
    args = parser.parse_args(argv)
    if not re.fullmatch(r"[0-9a-f]{40}", args.code_commit):
        raise ValueError("--code-commit must be a full lowercase Git SHA")

    project_root = Path(args.project_root).resolve()
    code_root = Path(args.code_root).resolve()
    protocol_root = Path(args.protocol_freeze).resolve()
    model_root = Path(args.composite_model_freeze).resolve()
    output = Path(args.output_dir).resolve()
    partial = Path(str(output) + ".partial")
    if output.exists() or partial.exists():
        raise FileExistsError("Refusing to overwrite Populus execution freeze")
    patch_values = (
        args.supersedes_execution_freeze,
        args.failed_blind_attempt,
        args.patch_reason_file,
    )
    patch_mode = all(patch_values)
    if any(patch_values) and not patch_mode:
        raise ValueError(
            "Execution patch requires --supersedes-execution-freeze, "
            "--failed-blind-attempt and --patch-reason-file together"
        )
    prior_execution: Path | None = None
    failed_attempt: Path | None = None
    patch_reason: Path | None = None
    if patch_mode:
        prior_execution = Path(args.supersedes_execution_freeze).resolve()
        failed_attempt = Path(args.failed_blind_attempt).resolve()
        patch_reason = Path(args.patch_reason_file).resolve()
        verify_sha256sums(prior_execution)
        prior_manifest = json.loads(
            (prior_execution / "execution_manifest.json").read_text(encoding="utf-8")
        )
        if prior_manifest.get("schema_version") != SCHEMA_VERSION:
            raise ValueError("Superseded freeze has the wrong schema")
        if not patch_reason.is_file() or patch_reason.stat().st_size == 0:
            raise ValueError("Execution patch reason file is missing")
        if not patch_reason.is_relative_to(code_root):
            raise ValueError("Execution patch reason must be committed inside code root")
        for forbidden in (
            "results/blind_runs/populus_external_v0.4",
            "results/evaluator/populus/v0.4/reveal_inputs",
            "results/copy_collapse/external/populus_v0.4_reveal",
        ):
            if (project_root / forbidden).exists():
                raise ValueError(f"Execution patch is too late; revealed artifact exists: {forbidden}")
    else:
        for relative in FORBIDDEN_PRE_FREEZE_PATHS:
            if (project_root / relative).exists():
                raise ValueError(
                    f"Execution freeze is too late; target artifact exists: {relative}"
                )
    verify_sha256sums(protocol_root)
    verify_sha256sums(model_root)
    protocol_manifest = json.loads(
        (protocol_root / "protocol_manifest.json").read_text(encoding="utf-8")
    )
    if protocol_manifest.get("schema_version") != PROTOCOL_SCHEMA:
        raise ValueError("Wrong Populus protocol freeze schema")
    environments = dict(args.environment)
    required_environments = {
        "ploidypatch-dev",
        "ploidypatch-model",
        "ploidypatch-baseline",
        "ploidypatch-synteny",
        "ploidypatch-syngap",
        "ploidypatch-gemoma",
        "ploidypatch-lifton",
    }
    if len(environments) != len(args.environment) or set(environments) != required_environments:
        raise ValueError(
            f"Environment names must be exactly: {sorted(required_environments)}"
        )
    files = implementation_files(code_root)
    missing = [relative for relative in files if not (code_root / relative).is_file()]
    if missing:
        raise ValueError(f"Missing execution implementation: {', '.join(missing)}")
    verify_git_state(code_root, args.code_commit)

    partial.mkdir(parents=True)
    try:
        archive = partial / "source.tar"
        create_git_archive(code_root, args.code_commit, archive)
        names = archive_names(archive)
        absent = [relative for relative in files if relative not in names]
        if absent:
            raise ValueError(
                "Execution implementation is not in committed source archive: "
                + ", ".join(absent)
            )
        source = partial / "source"
        source.mkdir()
        extract_archive_safely(archive, source)

        locks = partial / "environment_locks"
        locks.mkdir()
        environment_rows: list[dict[str, Any]] = []
        for name, prefix in sorted(environments.items()):
            explicit = environment_lock(prefix)
            pip = pip_lock(prefix)
            explicit_path = locks / f"{name}.explicit.txt"
            pip_path = locks / f"{name}.pip-freeze.txt"
            explicit_path.write_bytes(explicit)
            pip_path.write_bytes(pip)
            environment_rows.append(
                {
                    "name": name,
                    "host_prefix": str(prefix),
                    "explicit_lock": explicit_path.relative_to(partial).as_posix(),
                    "explicit_sha256": sha256(explicit_path),
                    "pip_lock": pip_path.relative_to(partial).as_posix(),
                    "pip_sha256": sha256(pip_path),
                }
            )

        with (partial / "implementation_manifest.tsv").open(
            "x", encoding="utf-8", newline=""
        ) as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(("relative_path", "bytes", "sha256"))
            for relative in files:
                path = source / relative
                writer.writerow((relative, path.stat().st_size, sha256(path)))
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

        failed_exit_status: int | None = None
        if patch_mode:
            assert failed_attempt is not None
            failed_exit_status = freeze_failed_attempt_manifest(
                failed_attempt, partial / "superseded_failed_attempt_manifest.tsv"
            )

        contract = {
            "schema_version": SCHEMA_VERSION,
            "code_commit": args.code_commit,
            "freeze_stage": (
                "post_failed_blind_pre_score_pre_reveal_execution_patch"
                if patch_mode
                else "post_metadata_pre_pair_pre_candidate_pre_label"
            ),
            "created_before": {
                "wgd_pair_enumeration": not patch_mode,
                "candidate_generation": not patch_mode,
                "candidate_labels": True,
                "candidate_scores": True,
            },
            "network_access_in_blind_runner": False,
            "nas_data_mount_in_blind_runner": False,
            "complete_target_annotation_mount_in_blind_runner": False,
            "evaluator_only_mount_in_blind_runner": False,
            "truth_or_label_mount_in_blind_runner": False,
            "protocol_SHA256SUMS_sha256": sha256(protocol_root / "SHA256SUMS"),
            "composite_model_SHA256SUMS_sha256": sha256(model_root / "SHA256SUMS"),
            "source_archive": {
                "format": "git_archive_tar",
                "path": "source.tar",
                "sha256": sha256(archive),
            },
            "implementation_file_count": len(files),
            "environments": environment_rows,
        }
        if patch_mode:
            assert prior_execution is not None
            assert failed_attempt is not None
            assert patch_reason is not None
            contract["execution_patch"] = {
                "pre_reveal": True,
                "scientific_protocol_changed": False,
                "model_or_threshold_changed": False,
                "supersedes_execution_freeze": str(prior_execution),
                "superseded_SHA256SUMS_sha256": sha256(
                    prior_execution / "SHA256SUMS"
                ),
                "failed_attempt": str(failed_attempt),
                "failed_attempt_exit_status": failed_exit_status,
                "failed_attempt_manifest": "superseded_failed_attempt_manifest.tsv",
                "patch_reason_file": patch_reason.relative_to(code_root).as_posix(),
                "patch_reason_sha256": sha256(patch_reason),
            }
        (partial / "execution_manifest.json").write_text(
            json.dumps(contract, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        make_checksums(partial)
        verify_sha256sums(partial)
        chmod_read_only(partial)
        os.replace(partial, output)
    except BaseException:
        if partial.exists():
            shutil.rmtree(partial)
        raise
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
