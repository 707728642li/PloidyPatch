#!/usr/bin/env python3
"""Strict streaming FASTQ length/mean-quality filter.

The filtered records are written to stdout so large ONT inputs can be passed
directly to a full-length cDNA classifier without materializing another copy.
"""

from __future__ import annotations

import argparse
import gzip
import json
import math
import sys
from pathlib import Path
from typing import BinaryIO


def filter_fastq_stream(
    source: BinaryIO,
    destination: BinaryIO,
    *,
    minimum_length: int,
    minimum_mean_quality: float | None,
    phred_offset: int = 33,
) -> dict[str, int | float]:
    if minimum_length < 1:
        raise ValueError("minimum_length must be positive")
    if minimum_mean_quality is not None and minimum_mean_quality < 0:
        raise ValueError("minimum_mean_quality must be nonnegative")
    if phred_offset < 0:
        raise ValueError("phred_offset must be nonnegative")
    input_records = input_bases = passed_records = passed_bases = 0
    length_failed = quality_failed = 0
    while True:
        header = source.readline()
        if not header:
            break
        sequence_line = source.readline()
        plus = source.readline()
        quality_line = source.readline()
        input_records += 1
        if not sequence_line or not plus or not quality_line:
            raise ValueError(f"Truncated FASTQ record {input_records}")
        if not header.startswith(b"@") or not plus.startswith(b"+"):
            raise ValueError(f"Malformed FASTQ record {input_records}")
        sequence = sequence_line.rstrip(b"\r\n")
        quality = quality_line.rstrip(b"\r\n")
        if len(sequence) != len(quality):
            raise ValueError(
                f"Sequence/quality length mismatch in FASTQ record {input_records}"
            )
        if minimum_mean_quality is not None and quality and min(quality) < phred_offset:
            raise ValueError(f"Quality below PHRED offset in record {input_records}")
        length = len(sequence)
        input_bases += length
        if length < minimum_length:
            length_failed += 1
            continue
        # ONT/Pychopper report the mean read Q score by averaging base-error
        # probabilities and converting the mean back to PHRED space.  An
        # arithmetic mean of per-base Q values can materially overestimate
        # mixed-quality reads and would not reproduce Pychopper's -Q filter.
        if minimum_mean_quality is not None:
            mean_error_probability = (
                sum(
                    10.0 ** (-(value - phred_offset) / 10.0)
                    for value in quality
                )
                / length
                if length
                else 1.0
            )
            mean_quality = -10.0 * math.log10(mean_error_probability)
            if mean_quality < minimum_mean_quality:
                quality_failed += 1
                continue
        destination.write(header)
        destination.write(sequence_line)
        destination.write(plus)
        destination.write(quality_line)
        passed_records += 1
        passed_bases += length
    return {
        "input_records": input_records,
        "input_bases": input_bases,
        "passed_records": passed_records,
        "passed_bases": passed_bases,
        "length_failed_records": length_failed,
        "quality_failed_records": quality_failed,
        "minimum_length": minimum_length,
        "minimum_mean_quality": minimum_mean_quality,
        "phred_offset": phred_offset,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("input_fastq")
    parser.add_argument("--minimum-length", type=int, required=True)
    parser.add_argument("--minimum-mean-quality", type=float)
    parser.add_argument("--phred-offset", type=int, default=33)
    parser.add_argument("--summary-json", required=True)
    args = parser.parse_args()
    input_path = Path(args.input_fastq)
    summary_path = Path(args.summary_json)
    if not input_path.is_file() or input_path.stat().st_size == 0:
        raise ValueError(f"Missing or empty FASTQ input: {input_path}")
    if summary_path.exists():
        raise FileExistsError(f"Refusing to overwrite summary: {summary_path}")
    opener = gzip.open if input_path.suffix == ".gz" else open
    with opener(input_path, "rb") as source:
        report = filter_fastq_stream(
            source,
            sys.stdout.buffer,
            minimum_length=args.minimum_length,
            minimum_mean_quality=args.minimum_mean_quality,
            phred_offset=args.phred_offset,
        )
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    with summary_path.open("x", encoding="utf-8", newline="") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
