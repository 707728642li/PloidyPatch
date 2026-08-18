#!/usr/bin/env python3
"""Build and smoke-test a checksum-bound PloidyPatch release candidate."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import shutil
import subprocess
import sys
import tarfile
import tomllib
from typing import Any, Sequence
import venv
from zipfile import ZipFile


SCHEMA = "ploidypatch.release_candidate_evidence.v2"
PACKAGE = "ploidypatch"
ROOT = Path(__file__).resolve().parents[1]
VERSION = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))["project"]["version"]
WHEEL_NAME = f"{PACKAGE}-{VERSION}-py3-none-any.whl"
SDIST_NAME = f"{PACKAGE}-{VERSION}.tar.gz"
COMMAND_FAMILIES = (
    "audit",
    "baseline",
    "benchmark",
    "evidence",
    "graph",
    "normalize",
    "patch",
    "report",
)
EXAMPLE_FILES = (
    "candidate.gff3",
    "candidate.gff3.manifest.json",
    "decisions.tsv",
    "input_manifest.tsv",
    "review_decisions.tsv",
    "report_copy_features.tsv",
    "report_scores.tsv",
    "report_topology.tsv",
    "run_example.py",
    "source.gff3",
)
REQUIRED_SDIST_SUFFIXES = (
    "README.md",
    "pyproject.toml",
    "docs/USER_GUIDE.md",
    "docs/REPRODUCIBILITY_GUIDE.md",
    "docs/CLI_COMMAND_INVENTORY_v0.1.json",
    "examples/minimal_reviewed_patch/run_example.py",
    "examples/minimal_reviewed_patch/input_manifest.tsv",
    "containers/Dockerfile",
    "scripts/smoke_container_v0.1.sh",
)


class ReleaseCandidateError(RuntimeError):
    """Raised when a release-candidate invariant is violated."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ReleaseCandidateError(message)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(8 * 1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def command(
    argv: Sequence[str | os.PathLike[str]],
    *,
    cwd: Path,
    records: list[dict[str, Any]],
    label: str,
    environment: dict[str, str] | None = None,
) -> subprocess.CompletedProcess[str]:
    normalized = [os.fspath(item) for item in argv]
    completed = subprocess.run(
        normalized,
        cwd=cwd,
        env=environment,
        text=True,
        capture_output=True,
        check=False,
    )
    records.append(
        {
            "argument_count": len(normalized),
            "argv_sha256": sha256_bytes("\0".join(normalized).encode("utf-8")),
            "executable": Path(normalized[0]).name,
            "label": label,
            "returncode": completed.returncode,
            "stderr_sha256": sha256_bytes(completed.stderr.encode("utf-8")),
            "stdout_sha256": sha256_bytes(completed.stdout.encode("utf-8")),
        }
    )
    if completed.returncode != 0:
        raise ReleaseCandidateError(
            f"{label} failed with exit code {completed.returncode}:\n{completed.stderr[-4000:]}"
        )
    return completed


def git_state(project_root: Path) -> tuple[str, str]:
    head = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=project_root,
        text=True,
        capture_output=True,
        check=True,
    ).stdout.strip()
    status = subprocess.run(
        ["git", "status", "--porcelain=v1", "--untracked-files=all"],
        cwd=project_root,
        text=True,
        capture_output=True,
        check=True,
    ).stdout
    require(len(head) == 40 and all(character in "0123456789abcdef" for character in head), "invalid git HEAD")
    return head, status


def safe_archive_name(raw: str) -> PurePosixPath:
    require(raw != "" and "\\" not in raw and "\x00" not in raw, f"unsafe archive path: {raw!r}")
    path = PurePosixPath(raw)
    require(
        not path.is_absolute() and all(part not in {"", ".", ".."} for part in path.parts),
        f"unsafe archive path: {raw!r}",
    )
    return path


def audit_wheel(path: Path) -> dict[str, Any]:
    with ZipFile(path) as archive:
        names = archive.namelist()
        require(names and len(names) == len(set(name.casefold() for name in names)), "wheel path collision")
        for name in names:
            safe_archive_name(name)
        metadata_names = [name for name in names if name.endswith(".dist-info/METADATA")]
        entry_names = [name for name in names if name.endswith(".dist-info/entry_points.txt")]
        require(len(metadata_names) == 1 and len(entry_names) == 1, "wheel metadata universe differs")
        metadata = archive.read(metadata_names[0]).decode("utf-8").replace("\r\n", "\n")
        entries = archive.read(entry_names[0]).decode("utf-8").replace("\r\n", "\n")
    require("Name: ploidypatch\n" in metadata, "wheel project name differs")
    require(f"Version: {VERSION}\n" in metadata, "wheel version differs")
    require("Requires-Python: >=3.11\n" in metadata, "wheel Python requirement differs")
    mandatory = [
        line
        for line in metadata.splitlines()
        if line.startswith("Requires-Dist:") and "; extra == " not in line
    ]
    require(not mandatory, f"wheel has unexpected mandatory runtime dependencies: {mandatory}")
    require("ploidypatch = ploidypatch.cli:main" in entries, "console entry point differs")
    return {
        "bytes": path.stat().st_size,
        "members": len(names),
        "name": path.name,
        "no_mandatory_runtime_dependencies": True,
        "sha256": sha256_file(path),
    }


def audit_sdist(path: Path) -> tuple[dict[str, Any], list[str]]:
    with tarfile.open(path, "r:gz") as archive:
        members = archive.getmembers()
        names = [member.name for member in members]
        require(names and len(names) == len(set(name.casefold() for name in names)), "sdist path collision")
        for member in members:
            safe_archive_name(member.name)
            require(member.isfile() or member.isdir(), f"unsupported sdist member type: {member.name}")
    root = f"{PACKAGE}-{VERSION}/"
    require(all(name == root.rstrip("/") or name.startswith(root) for name in names), "sdist root differs")
    for suffix in REQUIRED_SDIST_SUFFIXES:
        require(f"{root}{suffix}" in names, f"required sdist member is missing: {suffix}")
    return (
        {
            "bytes": path.stat().st_size,
            "members": len(names),
            "name": path.name,
            "required_public_assets": list(REQUIRED_SDIST_SUFFIXES),
            "sha256": sha256_file(path),
        },
        names,
    )


def write_sha256sums(root: Path) -> None:
    destination = root / "SHA256SUMS"
    require(not destination.exists(), f"refusing to overwrite checksum manifest: {destination}")
    paths = sorted(
        (path for path in root.rglob("*") if path.is_file()),
        key=lambda item: item.relative_to(root).as_posix().encode("utf-8"),
    )
    with destination.open("x", encoding="utf-8", newline="\n") as handle:
        for path in paths:
            handle.write(f"{sha256_file(path)}  {path.relative_to(root).as_posix()}\n")


def build(project_root: Path, work_root: Path, evidence_root: Path) -> None:
    project_root = project_root.resolve(strict=True)
    work_root = work_root.resolve()
    evidence_root = evidence_root.resolve()
    work_working = Path(f"{work_root}.working")
    evidence_working = Path(f"{evidence_root}.working")
    for path in (work_root, work_working, evidence_root, evidence_working):
        require(not path.exists() and not path.is_symlink(), f"refusing to overwrite: {path}")
    head, status = git_state(project_root)
    require(status == "", "release candidate requires a clean git worktree")

    records: list[dict[str, Any]] = []
    work_working.mkdir(parents=True)
    dist = work_working / "dist"
    command(
        [sys.executable, "-m", "build", "--outdir", dist, project_root],
        cwd=project_root,
        records=records,
        label="pep517_build",
    )
    wheel = dist / WHEEL_NAME
    sdist = dist / SDIST_NAME
    observed = sorted(path.name for path in dist.iterdir() if path.is_file())
    require(observed == sorted((WHEEL_NAME, SDIST_NAME)), f"distribution universe differs: {observed}")
    wheel_record = audit_wheel(wheel)
    sdist_record, sdist_members = audit_sdist(sdist)

    smoke_env = work_working / "smoke_env"
    venv.EnvBuilder(with_pip=True, clear=False, symlinks=False).create(smoke_env)
    smoke_python = smoke_env / ("Scripts/python.exe" if os.name == "nt" else "bin/python")
    require(smoke_python.is_file(), "clean smoke Python is missing")
    smoke_environment = os.environ.copy()
    smoke_environment.pop("PYTHONHOME", None)
    smoke_environment.pop("PYTHONPATH", None)
    command(
        [smoke_python, "-m", "pip", "install", "--no-deps", wheel],
        cwd=work_working,
        records=records,
        label="install_wheel_no_deps",
        environment=smoke_environment,
    )
    pip_check = command(
        [smoke_python, "-m", "pip", "check"],
        cwd=work_working,
        records=records,
        label="pip_check",
        environment=smoke_environment,
    )
    require("No broken requirements found" in pip_check.stdout, "pip check did not report a clean environment")
    version = command(
        [smoke_python, "-m", "ploidypatch.cli", "--version"],
        cwd=work_working,
        records=records,
        label="installed_version",
        environment=smoke_environment,
    ).stdout.strip()
    require(version == VERSION, f"installed version differs: {version!r}")
    help_text = command(
        [smoke_python, "-m", "ploidypatch.cli", "--help"],
        cwd=work_working,
        records=records,
        label="installed_help",
        environment=smoke_environment,
    ).stdout
    require(all(family in help_text for family in COMMAND_FAMILIES), "installed help omits a command family")

    source_example = project_root / "examples/minimal_reviewed_patch"
    copied_example = work_working / "example/input"
    copied_example.mkdir(parents=True)
    for name in EXAMPLE_FILES:
        source = source_example / name
        require(source.is_file() and not source.is_symlink(), f"missing example input: {source}")
        shutil.copyfile(source, copied_example / name)
    example_output = work_working / "example/output"
    example_run = command(
        [
            smoke_python,
            copied_example / "run_example.py",
            "--input-dir",
            copied_example,
            "--output-dir",
            example_output,
        ],
        cwd=work_working,
        records=records,
        label="installed_wheel_reviewed_patch_example",
        environment=smoke_environment,
    )
    example_stdout = json.loads(example_run.stdout)
    summary = json.loads((example_output / "run_summary.json").read_text(encoding="utf-8"))
    require(example_stdout["automatic_approval"] is False, "example enabled automatic approval")
    require(summary["accepted_additions"] == 2, "example accepted-addition count differs")
    require(summary["byte_identical_reversion"] is True, "example reversion is not byte-identical")
    require(summary["source_sha256"] == summary["reverted_sha256"], "example source/revert hashes differ")

    evidence_working.mkdir(parents=True)
    with ZipFile(wheel) as archive:
        wheel_members = archive.namelist()
    (evidence_working / "wheel_members.txt").write_text(
        "\n".join(sorted(wheel_members, key=lambda value: value.encode("utf-8"))) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    (evidence_working / "sdist_members.txt").write_text(
        "\n".join(sorted(sdist_members, key=lambda value: value.encode("utf-8"))) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    shutil.copyfile(example_output / "run_summary.json", evidence_working / "example_run_summary.json")
    shutil.copyfile(example_output / "SHA256SUMS", evidence_working / "example_output_SHA256SUMS")
    with (evidence_working / "command_audit.jsonl").open("x", encoding="utf-8", newline="\n") as handle:
        for record in records:
            handle.write(json.dumps(record, sort_keys=True) + "\n")
    manifest = {
        "schema_version": SCHEMA,
        "source_commit": head,
        "source_worktree_clean": True,
        "build": {
            "backend": "setuptools.build_meta",
            "frontend": "PyPA build",
            "isolation": True,
            "python": sys.version.split()[0],
        },
        "distributions": {"sdist": sdist_record, "wheel": wheel_record},
        "installed_smoke": {
            "command_families": list(COMMAND_FAMILIES),
            "example_accepted_additions": 2,
            "example_automatic_approval": False,
            "example_byte_identical_reversion": True,
            "pip_check": "passed",
            "source_tree_pythonpath_used": False,
            "version": version,
        },
        "release_boundary": {
            "distribution_metadata_complete": True,
            "github_publication_external_state_not_inferred": True,
            "archive_doi": None,
            "archive_doi_required_for_github_release": False,
        },
    }
    (evidence_working / "release_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    write_sha256sums(evidence_working)
    os.replace(work_working, work_root)
    os.replace(evidence_working, evidence_root)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--work-root", type=Path, required=True)
    parser.add_argument("--evidence-root", type=Path, required=True)
    args = parser.parse_args()
    build(args.project_root, args.work_root, args.evidence_root)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
