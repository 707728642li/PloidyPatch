from __future__ import annotations

import gzip
import json
import subprocess
import sys
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "scripts" / "filter_fastq_length_quality.py"


def test_streaming_fastq_filter_applies_length_and_mean_quality(tmp_path: Path) -> None:
    input_path = tmp_path / "reads.fastq.gz"
    with gzip.open(input_path, "wt", encoding="ascii", newline="") as handle:
        handle.write(
            "@pass\nACGTA\n+\nIIIII\n"
            "@too_short\nACG\n+\nIII\n"
            "@low_quality\nACGTA\n+\n#####\n"
        )
    summary = tmp_path / "summary.json"
    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            str(input_path),
            "--minimum-length",
            "5",
            "--minimum-mean-quality",
            "7",
            "--summary-json",
            str(summary),
        ],
        check=True,
        stdout=subprocess.PIPE,
    )

    assert result.stdout == b"@pass\nACGTA\n+\nIIIII\n"
    report = json.loads(summary.read_text(encoding="utf-8"))
    assert report == {
        "input_bases": 13,
        "input_records": 3,
        "length_failed_records": 1,
        "minimum_length": 5,
        "minimum_mean_quality": 7.0,
        "passed_bases": 5,
        "passed_records": 1,
        "phred_offset": 33,
        "quality_failed_records": 1,
    }


def test_streaming_fastq_filter_rejects_malformed_records(tmp_path: Path) -> None:
    input_path = tmp_path / "bad.fastq"
    input_path.write_bytes(b"@bad\nACGT\n+\nIII\n")
    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            str(input_path),
            "--minimum-length",
            "1",
            "--minimum-mean-quality",
            "0",
            "--summary-json",
            str(tmp_path / "summary.json"),
        ],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    assert result.returncode != 0
    assert b"Sequence/quality length mismatch" in result.stderr


def test_streaming_fastq_filter_uses_error_probability_mean_q(
    tmp_path: Path,
) -> None:
    input_path = tmp_path / "mixed_quality.fastq"
    # PHRED scores 0 and 20 have arithmetic mean Q=10, but their mean error
    # probability corresponds to Q~=2.97.  Pychopper therefore rejects this
    # record at Q7 and the preprocessing layer must do the same.
    input_path.write_bytes(b"@mixed\nAC\n+\n!5\n")
    summary = tmp_path / "summary.json"
    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            str(input_path),
            "--minimum-length",
            "1",
            "--minimum-mean-quality",
            "7",
            "--summary-json",
            str(summary),
        ],
        check=True,
        stdout=subprocess.PIPE,
    )

    assert result.stdout == b""
    report = json.loads(summary.read_text(encoding="utf-8"))
    assert report["quality_failed_records"] == 1
    assert report["passed_records"] == 0


def test_streaming_fastq_filter_can_defer_quality_gate(tmp_path: Path) -> None:
    input_path = tmp_path / "length_only.fastq"
    input_path.write_bytes(b"@kept\nACGTA\n+\n#####\n@short\nAC\n+\nII\n")
    summary = tmp_path / "summary.json"
    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            str(input_path),
            "--minimum-length",
            "5",
            "--summary-json",
            str(summary),
        ],
        check=True,
        stdout=subprocess.PIPE,
    )

    assert result.stdout == b"@kept\nACGTA\n+\n#####\n"
    report = json.loads(summary.read_text(encoding="utf-8"))
    assert report["minimum_mean_quality"] is None
    assert report["quality_failed_records"] == 0
    assert report["length_failed_records"] == 1
