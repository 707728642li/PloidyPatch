#!/usr/bin/env python3
"""Run the distributable PloidyPatch reviewed-copy example end to end."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path
from typing import Any


INPUT_FILES = frozenset(
    {
        "source.gff3",
        "candidate.gff3",
        "decisions.tsv",
        "candidate.gff3.manifest.json",
        "review_decisions.tsv",
        "report_copy_features.tsv",
        "report_scores.tsv",
        "report_topology.tsv",
    }
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_inputs(input_dir: Path) -> dict[str, dict[str, Any]]:
    manifest_path = input_dir / "input_manifest.tsv"
    if not manifest_path.is_file() or manifest_path.is_symlink():
        raise ValueError(f"missing regular input manifest: {manifest_path}")
    with manifest_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != ("file_name", "bytes", "sha256"):
            raise ValueError("input manifest fields differ")
        rows = list(reader)
    indexed = {row["file_name"]: row for row in rows}
    if len(indexed) != len(rows) or set(indexed) != INPUT_FILES:
        raise ValueError("input manifest file universe differs")
    verified: dict[str, dict[str, Any]] = {}
    for file_name in sorted(INPUT_FILES):
        if Path(file_name).name != file_name:
            raise ValueError(f"unsafe input file name: {file_name}")
        path = input_dir / file_name
        if not path.is_file() or path.is_symlink():
            raise ValueError(f"missing regular example input: {file_name}")
        expected_bytes = int(indexed[file_name]["bytes"])
        expected_sha256 = indexed[file_name]["sha256"]
        actual_sha256 = sha256_file(path)
        if path.stat().st_size != expected_bytes or actual_sha256 != expected_sha256:
            raise ValueError(f"input bytes or SHA-256 differ: {file_name}")
        verified[file_name] = {
            "bytes": expected_bytes,
            "sha256": actual_sha256,
        }
    return verified


def run_cli(args: list[str], *, repo_root: Path) -> dict[str, Any]:
    environment = os.environ.copy()
    source_root = repo_root / "src"
    if source_root.is_dir():
        prior = environment.get("PYTHONPATH", "")
        environment["PYTHONPATH"] = os.fspath(source_root) + (
            os.pathsep + prior if prior else ""
        )
    completed = subprocess.run(
        [sys.executable, "-m", "ploidypatch.cli", *args],
        cwd=repo_root,
        env=environment,
        text=True,
        capture_output=True,
        check=False,
    )
    if completed.returncode != 0:
        raise RuntimeError(
            "PloidyPatch command failed: "
            + " ".join(args)
            + "\n"
            + completed.stderr.strip()
        )
    try:
        report = json.loads(completed.stdout)
    except json.JSONDecodeError as error:
        raise RuntimeError("PloidyPatch command did not emit one JSON report") from error
    if not isinstance(report, dict):
        raise RuntimeError("PloidyPatch command report is not a JSON object")
    return report


def write_sha256sums(root: Path) -> None:
    paths = sorted(
        path
        for path in root.rglob("*")
        if path.is_file() and path != root / "SHA256SUMS"
    )
    with (root / "SHA256SUMS").open("x", encoding="utf-8", newline="") as handle:
        for path in paths:
            handle.write(f"{sha256_file(path)}  {path.relative_to(root).as_posix()}\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--output-dir", type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    input_dir = args.input_dir.resolve(strict=True)
    output_dir = (args.output_dir or (input_dir / "output")).resolve()
    working = Path(f"{output_dir}.working")
    if output_dir.exists() or working.exists():
        raise FileExistsError(f"refusing to overwrite example output: {output_dir}")
    verified_inputs = verify_inputs(input_dir)

    repo_root = Path(__file__).resolve().parents[2]
    working.mkdir(parents=True)
    command_log: list[dict[str, Any]] = []

    edits = working / "reviewed_copy_additions.edits.json"
    patch = working / "reviewed_copy_additions.patch.json"
    patched = working / "annotation.with_reviewed_additions.gff3"
    reverted = working / "annotation.reverted.gff3"
    commands = [
        [
            "patch",
            "compile-reviewed-copy-additions",
            "--annotation-gff",
            os.fspath(input_dir / "source.gff3"),
            "--candidate-gff",
            os.fspath(input_dir / "candidate.gff3"),
            "--pool-decisions",
            os.fspath(input_dir / "decisions.tsv"),
            "--pool-manifest",
            os.fspath(input_dir / "candidate.gff3.manifest.json"),
            "--review-decisions",
            os.fspath(input_dir / "review_decisions.tsv"),
            "--output-edits-json",
            os.fspath(edits),
        ],
        [
            "patch",
            "create",
            "--source-gff",
            os.fspath(input_dir / "source.gff3"),
            "--edits-json",
            os.fspath(edits),
            "--output-patch",
            os.fspath(patch),
        ],
        [
            "patch",
            "apply",
            "--source-gff",
            os.fspath(input_dir / "source.gff3"),
            "--patch",
            os.fspath(patch),
            "--output-gff",
            os.fspath(patched),
        ],
        [
            "patch",
            "revert",
            "--patched-gff",
            os.fspath(patched),
            "--patch",
            os.fspath(patch),
            "--output-gff",
            os.fspath(reverted),
        ],
    ]
    for command in commands:
        command_log.append({"argv": command, "report": run_cli(command, repo_root=repo_root)})

    source_bytes = (input_dir / "source.gff3").read_bytes()
    reverted_bytes = reverted.read_bytes()
    if source_bytes != reverted_bytes:
        raise RuntimeError("reverted annotation is not byte-identical to the source")
    patched_text = patched.read_text(encoding="utf-8")
    if not (
        "PPCONS_gene_a" in patched_text
        and "PPCONS_gene_c" in patched_text
        and "PPCONS_gene_b" not in patched_text
    ):
        raise RuntimeError("review decisions were not applied exactly")

    edit_report = command_log[0]["report"]
    summary = {
        "schema_version": "ploidypatch.minimal_reviewed_patch_example.v1",
        "automatic_approval": edit_report["policy"]["automatic_approval"],
        "accepted_additions": edit_report["counts"]["copy_addition_events"],
        "selected_candidate_digests": sorted(
            event["consensus_digest"] for event in edit_report["events"]
        ),
        "byte_identical_reversion": True,
        "source_sha256": sha256_file(input_dir / "source.gff3"),
        "reverted_sha256": sha256_file(reverted),
        "verified_inputs": verified_inputs,
    }
    (working / "run_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline=""
    )
    report_command = [
        "report",
        "--candidate-gff",
        os.fspath(input_dir / "candidate.gff3"),
        "--pool-decisions",
        os.fspath(input_dir / "decisions.tsv"),
        "--pool-manifest",
        os.fspath(input_dir / "candidate.gff3.manifest.json"),
        "--review-decisions",
        os.fspath(input_dir / "review_decisions.tsv"),
        "--scores",
        os.fspath(input_dir / "report_scores.tsv"),
        "--copy-features",
        os.fspath(input_dir / "report_copy_features.tsv"),
        "--topology-features",
        os.fspath(input_dir / "report_topology.tsv"),
        "--patch-edits",
        os.fspath(edits),
        "--run-summary",
        os.fspath(working / "run_summary.json"),
        "--title",
        "Minimal reviewed patch",
        "--output-dir",
        os.fspath(working / "report"),
        "--fail-on-attention",
    ]
    report_result = run_cli(report_command, repo_root=repo_root)
    if report_result["state"] != "validated_reversible_run":
        raise RuntimeError("example report did not reach validated_reversible_run")
    command_log.append({"argv": report_command, "report": report_result})
    with (working / "command_log.jsonl").open("x", encoding="utf-8", newline="") as handle:
        for item in command_log:
            handle.write(json.dumps(item, sort_keys=True) + "\n")
    write_sha256sums(working)
    working.rename(output_dir)
    print(json.dumps({"output_dir": os.fspath(output_dir), **summary}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
