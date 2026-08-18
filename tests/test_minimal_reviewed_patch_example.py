from __future__ import annotations

import hashlib
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).parents[1]
EXAMPLE = ROOT / "examples" / "minimal_reviewed_patch"
RUNNER = EXAMPLE / "run_example.py"
INPUT_FILES = (
    "input_manifest.tsv",
    "source.gff3",
    "candidate.gff3",
    "decisions.tsv",
    "candidate.gff3.manifest.json",
    "review_decisions.tsv",
    "report_copy_features.tsv",
    "report_scores.tsv",
    "report_topology.tsv",
)


def run_example(input_dir: Path, output_dir: Path) -> subprocess.CompletedProcess[str]:
    environment = os.environ.copy()
    environment["PYTHONPATH"] = os.fspath(ROOT / "src")
    return subprocess.run(
        [
            sys.executable,
            os.fspath(RUNNER),
            "--input-dir",
            os.fspath(input_dir),
            "--output-dir",
            os.fspath(output_dir),
        ],
        cwd=ROOT,
        env=environment,
        text=True,
        capture_output=True,
        check=False,
    )


def copy_inputs(destination: Path) -> None:
    destination.mkdir()
    for name in INPUT_FILES:
        shutil.copyfile(EXAMPLE / name, destination / name)


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_minimal_reviewed_patch_example_runs_end_to_end(tmp_path: Path) -> None:
    inputs = tmp_path / "inputs"
    output = tmp_path / "output"
    copy_inputs(inputs)
    completed = run_example(inputs, output)
    assert completed.returncode == 0, completed.stderr
    report = json.loads(completed.stdout)
    assert report["accepted_additions"] == 2
    assert report["automatic_approval"] is False
    assert report["byte_identical_reversion"] is True
    assert report["source_sha256"] == report["reverted_sha256"]
    assert report["selected_candidate_digests"] == ["a" * 64, "c" * 64]
    patched = (output / "annotation.with_reviewed_additions.gff3").read_text(encoding="utf-8")
    assert "PPCONS_gene_a" in patched
    assert "PPCONS_gene_b" not in patched
    assert "PPCONS_gene_c" in patched
    assert (output / "annotation.reverted.gff3").read_bytes() == (inputs / "source.gff3").read_bytes()
    assert len((output / "command_log.jsonl").read_text(encoding="utf-8").splitlines()) == 5
    assert (output / "report" / "index.html").is_file()
    assert json.loads((output / "report" / "report.json").read_text(encoding="utf-8"))[
        "summary"
    ]["state"] == "validated_reversible_run"
    checksum_lines = (output / "SHA256SUMS").read_text(encoding="utf-8").splitlines()
    recorded: dict[str, str] = {}
    for line in checksum_lines:
        digest, file_name = line.split("  ", maxsplit=1)
        assert file_name not in recorded
        recorded[file_name] = digest
    expected_files = {
        path.relative_to(output).as_posix()
        for path in output.rglob("*")
        if path.is_file() and path != output / "SHA256SUMS"
    }
    assert set(recorded) == expected_files
    for file_name, digest in recorded.items():
        assert sha256_file(output / Path(file_name)) == digest


def test_minimal_reviewed_patch_example_refuses_overwrite(tmp_path: Path) -> None:
    inputs = tmp_path / "inputs"
    output = tmp_path / "output"
    copy_inputs(inputs)
    assert run_example(inputs, output).returncode == 0
    repeated = run_example(inputs, output)
    assert repeated.returncode != 0
    assert "refusing to overwrite example output" in repeated.stderr


def test_minimal_reviewed_patch_example_rejects_changed_input(tmp_path: Path) -> None:
    inputs = tmp_path / "inputs"
    output = tmp_path / "output"
    copy_inputs(inputs)
    with (inputs / "candidate.gff3").open("a", encoding="utf-8", newline="") as handle:
        handle.write("# changed\n")
    completed = run_example(inputs, output)
    assert completed.returncode != 0
    assert "input bytes or SHA-256 differ: candidate.gff3" in completed.stderr
    assert not output.exists()
    assert not Path(f"{output}.working").exists()
