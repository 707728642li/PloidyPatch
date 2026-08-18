#!/usr/bin/env python3
"""Build, reproduce, install, and smoke-test PloidyPatch distributions twice."""

from __future__ import annotations

import argparse
from datetime import date, datetime, timezone
import gzip
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
from typing import Any
import venv
from zipfile import ZipFile


ROOT = Path(__file__).resolve().parents[1]
BASE_SCRIPT = ROOT / "scripts/build_release_candidate_v0.2.py"
SCHEMA = "ploidypatch.release_candidate_evidence.v3"
CANONICAL_SDIST_POLICY = {
    "archive_format": "pax",
    "directory_mode": "0755",
    "file_mode": "0755_if_source_executable_else_0644",
    "gzip_filename": "",
    "gzip_level": 9,
    "gid": 0,
    "gname": "",
    "member_order": "ascending_utf8_path_bytes",
    "mtime": "citation_release_date_midnight_utc",
    "pax_headers": "generated_only_when_required_by_tarfile",
    "uid": 0,
    "uname": "",
}


def _load_base():
    spec = importlib.util.spec_from_file_location("ploidypatch_release_builder_v02", BASE_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load base release builder: {BASE_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


base = _load_base()
ReleaseCandidateError = base.ReleaseCandidateError
require = base.require


def source_date_epoch(project_root: Path) -> int:
    """Return a version-stable epoch from the validated citation release date."""

    citation = project_root / "CITATION.cff"
    require(citation.is_file() and not citation.is_symlink(), "CITATION.cff is missing")
    values = [
        line.split(":", 1)[1].strip().strip('"\'')
        for line in citation.read_text(encoding="utf-8").splitlines()
        if line.startswith("date-released:")
    ]
    require(len(values) == 1, "CITATION.cff release date is missing or ambiguous")
    try:
        released = date.fromisoformat(values[0])
    except ValueError as error:
        raise ReleaseCandidateError("CITATION.cff release date is invalid") from error
    epoch = int(datetime.combine(released, datetime.min.time(), timezone.utc).timestamp())
    require(epoch > 0, "citation release-date timestamp must be positive")
    return epoch


def archive_doi(project_root: Path) -> str:
    citation = project_root / "CITATION.cff"
    require(citation.is_file() and not citation.is_symlink(), "CITATION.cff is missing")
    values = [
        line.split(":", 1)[1].strip().strip('"\'')
        for line in citation.read_text(encoding="utf-8").splitlines()
        if line.startswith("doi:")
    ]
    require(len(values) == 1 and values[0].startswith("10."), "CITATION.cff DOI is invalid")
    return values[0]


def _normalized_member(member: tarfile.TarInfo, epoch: int) -> tarfile.TarInfo:
    normalized = tarfile.TarInfo(member.name)
    normalized.type = member.type
    normalized.mode = 0o755 if member.isdir() or member.mode & 0o111 else 0o644
    normalized.mtime = epoch
    normalized.uid = 0
    normalized.gid = 0
    normalized.uname = ""
    normalized.gname = ""
    normalized.linkname = ""
    normalized.pax_headers = {}
    normalized.size = member.size if member.isfile() else 0
    return normalized


def canonicalize_sdist(source: Path, destination: Path, epoch: int) -> list[str]:
    """Rewrite a safe sdist with deterministic tar and gzip metadata."""

    source = source.resolve(strict=True)
    destination = destination.resolve()
    require(source != destination, "canonical sdist destination must differ from source")
    require(not destination.exists() and not destination.is_symlink(), f"refusing to overwrite: {destination}")
    require(epoch > 0, "canonical sdist epoch must be positive")
    try:
        with tarfile.open(source, "r:gz") as incoming:
            members = incoming.getmembers()
            names = [member.name for member in members]
            require(names, "sdist is empty")
            require(
                len(names) == len(set(name.casefold() for name in names)),
                "sdist path collision",
            )
            for member in members:
                base.safe_archive_name(member.name)
                require(
                    member.isfile() or member.isdir(),
                    f"unsupported sdist member type: {member.name}",
                )
            ordered = sorted(members, key=lambda item: item.name.encode("utf-8"))
            with destination.open("xb") as raw:
                with gzip.GzipFile(
                    filename="",
                    mode="wb",
                    compresslevel=9,
                    fileobj=raw,
                    mtime=epoch,
                ) as compressed:
                    with tarfile.open(
                        fileobj=compressed,
                        mode="w",
                        format=tarfile.PAX_FORMAT,
                    ) as outgoing:
                        for member in ordered:
                            normalized = _normalized_member(member, epoch)
                            if member.isfile():
                                payload = incoming.extractfile(member)
                                require(payload is not None, f"cannot read sdist member: {member.name}")
                                outgoing.addfile(normalized, payload)
                            else:
                                outgoing.addfile(normalized)
                raw.flush()
                os.fsync(raw.fileno())
        return [member.name for member in ordered]
    except BaseException:
        destination.unlink(missing_ok=True)
        raise


def _build_once(
    project_root: Path,
    root: Path,
    *,
    epoch: int,
    label: str,
    records: list[dict[str, Any]],
) -> tuple[dict[str, Any], dict[str, Any], list[str]]:
    root.mkdir(parents=True)
    dist = root / "dist"
    environment = os.environ.copy()
    environment["SOURCE_DATE_EPOCH"] = str(epoch)
    base.command(
        [sys.executable, "-m", "build", "--outdir", dist, project_root],
        cwd=project_root,
        records=records,
        label=f"{label}_pep517_build",
        environment=environment,
    )
    wheel = dist / base.WHEEL_NAME
    sdist = dist / base.SDIST_NAME
    observed = sorted(path.name for path in dist.iterdir() if path.is_file())
    require(
        observed == sorted((base.WHEEL_NAME, base.SDIST_NAME)),
        f"{label} distribution universe differs: {observed}",
    )
    canonical = dist / f"{base.SDIST_NAME}.canonical"
    canonical_members = canonicalize_sdist(sdist, canonical, epoch)
    os.replace(canonical, sdist)
    return base.audit_wheel(wheel), base.audit_sdist(sdist)[0], canonical_members


def _copy_example(project_root: Path, destination: Path) -> Path:
    destination.mkdir(parents=True)
    source_root = project_root / "examples/minimal_reviewed_patch"
    for name in base.EXAMPLE_FILES:
        source = source_root / name
        require(source.is_file() and not source.is_symlink(), f"missing example input: {source}")
        shutil.copyfile(source, destination / name)
    return destination


def _install_and_smoke(
    artifact: Path,
    *,
    kind: str,
    work_root: Path,
    project_root: Path,
    records: list[dict[str, Any]],
    epoch: int,
) -> dict[str, Any]:
    environment_root = work_root / f"{kind}_smoke_env"
    venv.EnvBuilder(with_pip=True, clear=False, symlinks=False).create(environment_root)
    python = environment_root / ("Scripts/python.exe" if os.name == "nt" else "bin/python")
    require(python.is_file(), f"{kind} smoke Python is missing")
    environment = os.environ.copy()
    environment.pop("PYTHONHOME", None)
    environment.pop("PYTHONPATH", None)
    environment["SOURCE_DATE_EPOCH"] = str(epoch)
    base.command(
        [python, "-m", "pip", "install", "--no-deps", artifact],
        cwd=work_root,
        records=records,
        label=f"install_{kind}_no_deps",
        environment=environment,
    )
    pip_check = base.command(
        [python, "-m", "pip", "check"],
        cwd=work_root,
        records=records,
        label=f"{kind}_pip_check",
        environment=environment,
    )
    require("No broken requirements found" in pip_check.stdout, f"{kind} pip check failed")
    version = base.command(
        [python, "-m", "ploidypatch.cli", "--version"],
        cwd=work_root,
        records=records,
        label=f"{kind}_installed_version",
        environment=environment,
    ).stdout.strip()
    require(version == base.VERSION, f"{kind} installed version differs: {version!r}")
    help_text = base.command(
        [python, "-m", "ploidypatch.cli", "--help"],
        cwd=work_root,
        records=records,
        label=f"{kind}_installed_help",
        environment=environment,
    ).stdout
    require(
        all(family in help_text for family in base.COMMAND_FAMILIES),
        f"{kind} installed help omits a command family",
    )
    example_input = _copy_example(project_root, work_root / f"{kind}_example/input")
    example_output = work_root / f"{kind}_example/output"
    completed = base.command(
        [
            python,
            example_input / "run_example.py",
            "--input-dir",
            example_input,
            "--output-dir",
            example_output,
        ],
        cwd=work_root,
        records=records,
        label=f"installed_{kind}_reviewed_patch_example",
        environment=environment,
    )
    stdout = json.loads(completed.stdout)
    summary_path = example_output / "run_summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    require(stdout["automatic_approval"] is False, f"{kind} example enabled automatic approval")
    require(summary["accepted_additions"] == 2, f"{kind} example accepted-addition count differs")
    require(summary["byte_identical_reversion"] is True, f"{kind} example did not revert exactly")
    require(summary["source_sha256"] == summary["reverted_sha256"], f"{kind} source/revert hashes differ")
    report = json.loads((example_output / "report/report.json").read_text(encoding="utf-8"))
    require(
        report["summary"]["state"] == "validated_reversible_run",
        f"{kind} report did not reach validated_reversible_run",
    )
    require(
        report["counts"]["automatic_approval_checked"] is True,
        f"{kind} report did not audit automatic-approval fields",
    )
    return {
        "example_output_checksum_manifest_sha256": base.sha256_file(example_output / "SHA256SUMS"),
        "example_report_checksum_manifest_sha256": base.sha256_file(
            example_output / "report/SHA256SUMS"
        ),
        "example_report_state": report["summary"]["state"],
        "example_run_summary_sha256": base.sha256_file(summary_path),
        "pip_check": "passed",
        "source_tree_pythonpath_used": False,
        "version": version,
    }


def _normalize_path_strings(value: Any, root: str) -> Any:
    if isinstance(value, str):
        return value.replace(root, "$EXAMPLE_ROOT")
    if isinstance(value, list):
        return [_normalize_path_strings(item, root) for item in value]
    if isinstance(value, dict):
        return {key: _normalize_path_strings(item, root) for key, item in value.items()}
    return value


def normalized_example_command_log(path: Path, example_root: Path) -> bytes:
    records = [
        json.loads(line)
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    require(records, f"example command log is empty: {path}")
    root = str(example_root.resolve(strict=True))
    normalized = [_normalize_path_strings(record, root) for record in records]
    return b"".join(
        (json.dumps(record, sort_keys=True) + "\n").encode("utf-8")
        for record in normalized
    )


def build(project_root: Path, work_root: Path, evidence_root: Path) -> None:
    project_root = project_root.resolve(strict=True)
    work_root = work_root.resolve()
    evidence_root = evidence_root.resolve()
    work_working = Path(f"{work_root}.working")
    evidence_working = Path(f"{evidence_root}.working")
    for path in (work_root, work_working, evidence_root, evidence_working):
        require(not path.exists() and not path.is_symlink(), f"refusing to overwrite: {path}")
    head, status = base.git_state(project_root)
    require(status == "", "release candidate requires a clean git worktree")
    epoch = source_date_epoch(project_root)
    records: list[dict[str, Any]] = []
    try:
        work_working.mkdir(parents=True)
        wheel_a, sdist_a, sdist_members = _build_once(
            project_root,
            work_working / "build_a",
            epoch=epoch,
            label="build_a",
            records=records,
        )
        wheel_b, sdist_b, members_b = _build_once(
            project_root,
            work_working / "build_b",
            epoch=epoch,
            label="build_b",
            records=records,
        )
        require(sdist_members == members_b, "canonical sdist member order differs between builds")
        require(wheel_a == wheel_b, "wheel records differ between independent builds")
        require(sdist_a == sdist_b, "sdist records differ between independent builds")
        for name in (base.WHEEL_NAME, base.SDIST_NAME):
            first = work_working / "build_a/dist" / name
            second = work_working / "build_b/dist" / name
            require(first.read_bytes() == second.read_bytes(), f"{name} is not byte-reproducible")

        wheel = work_working / "build_a/dist" / base.WHEEL_NAME
        sdist = work_working / "build_a/dist" / base.SDIST_NAME
        wheel_smoke = _install_and_smoke(
            wheel,
            kind="wheel",
            work_root=work_working,
            project_root=project_root,
            records=records,
            epoch=epoch,
        )
        sdist_smoke = _install_and_smoke(
            sdist,
            kind="sdist",
            work_root=work_working,
            project_root=project_root,
            records=records,
            epoch=epoch,
        )
        wheel_summary = work_working / "wheel_example/output/run_summary.json"
        sdist_summary = work_working / "sdist_example/output/run_summary.json"
        require(wheel_summary.read_bytes() == sdist_summary.read_bytes(), "wheel/sdist example summaries differ")
        wheel_output = work_working / "wheel_example/output"
        sdist_output = work_working / "sdist_example/output"
        wheel_files = {
            path.relative_to(wheel_output).as_posix(): path
            for path in wheel_output.rglob("*")
            if path.is_file()
        }
        sdist_files = {
            path.relative_to(sdist_output).as_posix(): path
            for path in sdist_output.rglob("*")
            if path.is_file()
        }
        require(wheel_files.keys() == sdist_files.keys(), "wheel/sdist example file universes differ")
        path_bound = {"SHA256SUMS", "command_log.jsonl"}
        for name in sorted(wheel_files.keys() - path_bound):
            require(
                wheel_files[name].read_bytes() == sdist_files[name].read_bytes(),
                f"wheel/sdist stable example output differs: {name}",
            )
        wheel_log = normalized_example_command_log(
            wheel_output / "command_log.jsonl", work_working / "wheel_example"
        )
        sdist_log = normalized_example_command_log(
            sdist_output / "command_log.jsonl", work_working / "sdist_example"
        )
        require(wheel_log == sdist_log, "wheel/sdist normalized command logs differ")

        evidence_working.mkdir(parents=True)
        with ZipFile(wheel) as archive:
            wheel_members = sorted(archive.namelist(), key=lambda value: value.encode("utf-8"))
        (evidence_working / "wheel_members.txt").write_text(
            "\n".join(wheel_members) + "\n", encoding="utf-8", newline="\n"
        )
        (evidence_working / "sdist_members.txt").write_text(
            "\n".join(sdist_members) + "\n", encoding="utf-8", newline="\n"
        )
        shutil.copyfile(wheel_summary, evidence_working / "example_run_summary.json")
        shutil.copyfile(
            work_working / "wheel_example/output/SHA256SUMS",
            evidence_working / "example_output_SHA256SUMS",
        )
        with (evidence_working / "command_audit.jsonl").open(
            "x", encoding="utf-8", newline="\n"
        ) as handle:
            for record in records:
                handle.write(json.dumps(record, sort_keys=True) + "\n")
        reproducibility = {
            "canonical_sdist_policy": CANONICAL_SDIST_POLICY,
            "independent_builds": 2,
            "normalized_example_command_logs_identical": True,
            "sdist_byte_identical": True,
            "source_date_epoch": epoch,
            "stable_example_outputs_identical": True,
            "wheel_byte_identical": True,
        }
        (evidence_working / "reproducibility_report.json").write_text(
            json.dumps(reproducibility, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
            newline="\n",
        )
        manifest = {
            "schema_version": SCHEMA,
            "source_commit": head,
            "source_worktree_clean": True,
            "build": {
                "backend": "setuptools.build_meta",
                "frontend": "PyPA build",
                "independent_builds": 2,
                "isolation": True,
                "python": sys.version.split()[0],
                "source_date_epoch": epoch,
            },
            "distributions": {"sdist": sdist_a, "wheel": wheel_a},
            "installed_smoke": {
                "command_families": list(base.COMMAND_FAMILIES),
                "example_accepted_additions": 2,
                "example_automatic_approval": False,
                "example_byte_identical_reversion": True,
                "sdist": sdist_smoke,
                "wheel": wheel_smoke,
            },
            "reproducibility": reproducibility,
            "release_boundary": {
                "artifact_status": "verified_release_candidate",
                "formal_publication_performed_by_builder": False,
                "distribution_metadata_complete": True,
                "github_publication_external_state_not_inferred": True,
                "archive_doi": archive_doi(project_root),
                "archive_doi_required_for_github_release": False,
            },
        }
        (evidence_working / "release_manifest.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
            newline="\n",
        )
        base.write_sha256sums(evidence_working)
        os.replace(work_working, work_root)
        os.replace(evidence_working, evidence_root)
    except BaseException:
        shutil.rmtree(evidence_working, ignore_errors=True)
        raise


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=ROOT)
    parser.add_argument("--work-root", type=Path, required=True)
    parser.add_argument("--evidence-root", type=Path, required=True)
    args = parser.parse_args()
    build(args.project_root, args.work_root, args.evidence_root)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
