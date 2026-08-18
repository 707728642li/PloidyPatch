from __future__ import annotations

import csv
import hashlib
import json
import re
import subprocess
import tempfile
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable

from . import __version__
from .perturb import _file_sha256


JUNCTION_MANIFEST_SCHEMA = "ploidypatch.rna_junction_extract.v1"
JUNCTION_AGGREGATE_SCHEMA = "ploidypatch.rna_junction_aggregate.v1"
JUNCTION_GROUP_AGGREGATE_SCHEMA = "ploidypatch.rna_junction_group_aggregate.v1"
CIGAR_PATTERN = re.compile(r"(\d+)([MIDNSHP=X])")
REFERENCE_CONSUMING = frozenset({"M", "D", "N", "=", "X"})
ALIGNED_REFERENCE_CONSUMING = frozenset({"M", "=", "X"})


def junctions_from_cigar(
    seqid: str, position: int, cigar: str
) -> tuple[tuple[str, int, int], ...]:
    """Return genomic exon-end/right-exon-start pairs from one SAM CIGAR."""

    if position < 1 or cigar == "*":
        return ()
    operations = CIGAR_PATTERN.findall(cigar)
    if not operations or "".join(length + code for length, code in operations) != cigar:
        raise ValueError(f"Malformed SAM CIGAR: {cigar}")
    reference_position = position
    junctions: list[tuple[str, int, int]] = []
    aligned_before = 0
    aligned_after = sum(
        int(length_text)
        for length_text, code in operations
        if code in ALIGNED_REFERENCE_CONSUMING
    )
    for length_text, code in operations:
        length = int(length_text)
        if length < 1:
            raise ValueError(f"Malformed SAM CIGAR operation: {cigar}")
        if code in ALIGNED_REFERENCE_CONSUMING:
            aligned_before += length
            aligned_after -= length
        if code == "N":
            if aligned_before > 0 and aligned_after > 0:
                junctions.append(
                    (seqid, reference_position - 1, reference_position + length)
                )
        if code in REFERENCE_CONSUMING:
            reference_position += length
    return tuple(junctions)


def _write_json_exclusive(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8", newline="") as handle:
        json.dump(value, handle, indent=2, sort_keys=True, ensure_ascii=False)
        handle.write("\n")


def extract_bam_junctions(
    bam_path: str | Path,
    output_tsv_path: str | Path,
    *,
    sample: str,
    samtools_path: str | Path,
    threads: int = 4,
    min_mapq: int = 20,
    excluded_flag_mask: int = 2308,
) -> dict[str, Any]:
    """Extract unstranded splice-junction read counts from a coordinate BAM."""

    bam = Path(bam_path)
    output = Path(output_tsv_path)
    manifest_path = Path(str(output) + ".manifest.json")
    samtools = Path(samtools_path)
    if output.exists() or manifest_path.exists():
        raise FileExistsError(f"Refusing to overwrite junction output: {output}")
    if not bam.is_file() or not Path(str(bam) + ".bai").is_file():
        raise FileNotFoundError(f"BAM or adjacent BAI is absent: {bam}")
    if not samtools.is_file():
        raise FileNotFoundError(f"samtools executable is absent: {samtools}")
    if not sample or any(character in sample for character in "\t\r\n"):
        raise ValueError("sample must be a non-empty single-field identifier")
    if threads < 1 or not 0 <= min_mapq <= 255 or excluded_flag_mask < 0:
        raise ValueError("Invalid threads, MAPQ, or flag-mask setting")

    header = subprocess.run(
        [str(samtools), "view", "-H", str(bam)],
        check=True,
        capture_output=True,
        text=True,
        encoding="utf-8",
    ).stdout
    header_sha = hashlib.sha256(header.encode("utf-8")).hexdigest()
    command = [
        str(samtools),
        "view",
        "-@",
        str(threads),
        "-q",
        str(min_mapq),
        "-F",
        str(excluded_flag_mask),
        str(bam),
    ]
    started = time.monotonic()
    counts: Counter[tuple[str, int, int]] = Counter()
    alignments = 0
    split_alignments = 0
    junction_observations = 0
    with tempfile.TemporaryFile(mode="w+t", encoding="utf-8") as stderr_handle:
        process = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=stderr_handle,
            text=True,
            encoding="utf-8",
            bufsize=1024 * 1024,
        )
        assert process.stdout is not None
        try:
            for line_number, raw_line in enumerate(process.stdout, start=1):
                fields = raw_line.rstrip("\r\n").split("\t")
                if len(fields) < 6:
                    raise ValueError(
                        f"Malformed SAM alignment at stream line {line_number}"
                    )
                try:
                    position = int(fields[3])
                except ValueError as error:
                    raise ValueError(
                        f"Invalid SAM position at stream line {line_number}"
                    ) from error
                alignments += 1
                junctions = junctions_from_cigar(fields[2], position, fields[5])
                if junctions:
                    split_alignments += 1
                    counts.update(junctions)
                    junction_observations += len(junctions)
        except Exception:
            process.kill()
            process.wait()
            raise
        return_code = process.wait()
        stderr_handle.seek(0)
        stderr = stderr_handle.read()
    if return_code != 0:
        raise RuntimeError(
            f"samtools view failed with exit {return_code}: {stderr.strip()}"
        )
    elapsed = time.monotonic() - started

    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=(
                "sample",
                "seqid",
                "left_exon_end",
                "right_exon_start",
                "read_count",
            ),
            delimiter="\t",
        )
        writer.writeheader()
        for (seqid, left, right), read_count in sorted(counts.items()):
            writer.writerow(
                {
                    "sample": sample,
                    "seqid": seqid,
                    "left_exon_end": left,
                    "right_exon_start": right,
                    "read_count": read_count,
                }
            )
    manifest: dict[str, Any] = {
        "schema_version": JUNCTION_MANIFEST_SCHEMA,
        "generator": {"name": "PloidyPatch", "version": __version__},
        "input": {
            "bam_path": str(bam.resolve()),
            "bam_bytes": bam.stat().st_size,
            "bai_bytes": Path(str(bam) + ".bai").stat().st_size,
            "header_sha256": header_sha,
        },
        "parameters": {
            "sample": sample,
            "threads": threads,
            "min_mapq": min_mapq,
            "excluded_flag_mask": excluded_flag_mask,
            "strand_policy": "unstranded_library_do_not_infer_from_read_orientation",
        },
        "command": command,
        "counts": {
            "alignments_passing_samtools_filters": alignments,
            "split_alignments": split_alignments,
            "junction_observations": junction_observations,
            "unique_junctions": len(counts),
        },
        "elapsed_seconds": elapsed,
        "samtools_stderr": stderr,
        "output": {
            "file_name": output.name,
            "rows": len(counts),
            "sha256": _file_sha256(output),
        },
    }
    _write_json_exclusive(manifest_path, manifest)
    return manifest


def _read_junction_file(path: Path) -> tuple[str, dict[tuple[str, int, int], int]]:
    if not path.is_file():
        raise FileNotFoundError(f"Junction TSV is absent: {path}")
    manifest_path = Path(str(path) + ".manifest.json")
    if not manifest_path.is_file():
        raise FileNotFoundError(f"Junction manifest is absent: {manifest_path}")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("schema_version") != JUNCTION_MANIFEST_SCHEMA:
        raise ValueError(f"Unsupported junction manifest: {manifest_path}")
    if manifest.get("output", {}).get("sha256") != _file_sha256(path):
        raise ValueError(f"Junction TSV checksum mismatch: {path}")
    rows: dict[tuple[str, int, int], int] = {}
    sample = ""
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        expected_fields = {
            "sample",
            "seqid",
            "left_exon_end",
            "right_exon_start",
            "read_count",
        }
        if set(reader.fieldnames or ()) != expected_fields:
            raise ValueError(f"Unexpected junction TSV columns: {path}")
        for line_number, row in enumerate(reader, start=2):
            if sample and row["sample"] != sample:
                raise ValueError(f"Multiple samples in {path} at line {line_number}")
            sample = row["sample"]
            key = (
                row["seqid"],
                int(row["left_exon_end"]),
                int(row["right_exon_start"]),
            )
            if not key[0] or key[1] < 1 or key[2] <= key[1] + 1:
                raise ValueError(
                    f"Invalid junction coordinates in {path} at line {line_number}"
                )
            if key in rows:
                raise ValueError(f"Duplicate junction in {path} at line {line_number}")
            read_count = int(row["read_count"])
            if read_count < 1:
                raise ValueError(f"Invalid read count in {path} at line {line_number}")
            rows[key] = read_count
    if not sample:
        sample = str(manifest["parameters"]["sample"])
    if sample != str(manifest["parameters"]["sample"]):
        raise ValueError(f"Sample identity differs from manifest: {path}")
    return sample, rows


def aggregate_junctions(
    junction_tsv_paths: Iterable[str | Path],
    output_tsv_path: str | Path,
    *,
    primary_samples: frozenset[str],
    min_reads_per_sample: int = 2,
    min_supporting_samples: int = 2,
    sample_group_tsv_path: str | Path | None = None,
    min_samples_per_group: int = 2,
    min_secondary_groups: int = 2,
) -> dict[str, Any]:
    """Aggregate held-back primary and secondary hawthorn RNA junction support."""

    output = Path(output_tsv_path)
    manifest_path = Path(str(output) + ".manifest.json")
    if output.exists() or manifest_path.exists():
        raise FileExistsError(f"Refusing to overwrite junction aggregate: {output}")
    if min_reads_per_sample < 1 or min_supporting_samples < 1:
        raise ValueError("Junction support thresholds must be positive")
    if min_samples_per_group < 1 or min_secondary_groups < 1:
        raise ValueError("Junction group thresholds must be positive")
    if not primary_samples:
        raise ValueError("At least one primary RNA sample is required")
    paths = sorted(Path(path) for path in junction_tsv_paths)
    if not paths:
        raise ValueError("No junction TSV inputs were supplied")

    sample_groups: dict[str, str] | None = None
    if sample_group_tsv_path is not None:
        sample_groups = {}
        with Path(sample_group_tsv_path).open(
            "r", encoding="utf-8", newline=""
        ) as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if set(reader.fieldnames or ()) != {"sample", "filename_stem_group"}:
                raise ValueError("Sample-group TSV must contain sample and filename_stem_group")
            for line_number, row in enumerate(reader, start=2):
                sample = row["sample"]
                group = row["filename_stem_group"]
                if not sample or not group or sample in sample_groups:
                    raise ValueError(
                        f"Invalid or duplicate sample-group row at line {line_number}"
                    )
                sample_groups[sample] = group

    samples: set[str] = set()
    support: dict[tuple[str, int, int], dict[str, int]] = defaultdict(dict)
    sources: list[dict[str, Any]] = []
    for path in paths:
        sample, rows = _read_junction_file(path)
        if sample in samples:
            raise ValueError(f"Duplicate junction sample: {sample}")
        samples.add(sample)
        for key, count in rows.items():
            support[key][sample] = count
        sources.append(
            {
                "sample": sample,
                "file_name": path.name,
                "sha256": _file_sha256(path),
                "manifest_sha256": _file_sha256(str(path) + ".manifest.json"),
            }
        )
    missing_primary = primary_samples - samples
    if missing_primary:
        raise ValueError(
            "Primary RNA sample(s) absent from junction inputs: "
            + ", ".join(sorted(missing_primary))
        )
    primary_groups: frozenset[str] = frozenset()
    if sample_groups is not None:
        missing_group_samples = samples - set(sample_groups)
        extra_group_samples = set(sample_groups) - samples
        if missing_group_samples or extra_group_samples:
            raise ValueError(
                "Sample-group/input mismatch; missing="
                + ",".join(sorted(missing_group_samples))
                + "; extra="
                + ",".join(sorted(extra_group_samples))
            )
        primary_groups = frozenset(sample_groups[sample] for sample in primary_samples)
        if len(primary_groups) != 1:
            raise ValueError("Primary samples must share exactly one filename-stem group")

    output.parent.mkdir(parents=True, exist_ok=True)
    fields = (
        "seqid",
        "left_exon_end",
        "right_exon_start",
        "primary_samples_ge_threshold",
        "primary_read_count",
        "secondary_samples_ge_threshold",
        "secondary_read_count",
        "total_samples_ge_threshold",
        "total_read_count",
        "primary_support",
        "secondary_support",
    )
    if sample_groups is not None:
        fields += (
            "primary_groups_ge_threshold",
            "secondary_groups_ge_threshold",
            "total_groups_ge_threshold",
            "primary_group_support",
            "secondary_group_support",
        )
    primary_supported = 0
    secondary_supported = 0
    primary_group_supported = 0
    secondary_group_supported = 0
    with output.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for (seqid, left, right), per_sample in sorted(support.items()):
            qualifying_primary = [
                sample
                for sample, count in per_sample.items()
                if sample in primary_samples and count >= min_reads_per_sample
            ]
            qualifying_secondary = [
                sample
                for sample, count in per_sample.items()
                if sample not in primary_samples and count >= min_reads_per_sample
            ]
            primary_state = len(qualifying_primary) >= min_supporting_samples
            secondary_state = len(qualifying_secondary) >= min_supporting_samples
            primary_supported += int(primary_state)
            secondary_supported += int(secondary_state)
            output_row: dict[str, Any] = {
                "seqid": seqid,
                "left_exon_end": left,
                "right_exon_start": right,
                "primary_samples_ge_threshold": len(qualifying_primary),
                "primary_read_count": sum(
                    count
                    for sample, count in per_sample.items()
                    if sample in primary_samples
                ),
                "secondary_samples_ge_threshold": len(qualifying_secondary),
                "secondary_read_count": sum(
                    count
                    for sample, count in per_sample.items()
                    if sample not in primary_samples
                ),
                "total_samples_ge_threshold": (
                    len(qualifying_primary) + len(qualifying_secondary)
                ),
                "total_read_count": sum(per_sample.values()),
                "primary_support": str(primary_state).lower(),
                "secondary_support": str(secondary_state).lower(),
            }
            if sample_groups is not None:
                qualifying_by_group: Counter[str] = Counter(
                    sample_groups[sample]
                    for sample, count in per_sample.items()
                    if count >= min_reads_per_sample
                )
                qualifying_groups = {
                    group
                    for group, group_samples in qualifying_by_group.items()
                    if group_samples >= min_samples_per_group
                }
                qualifying_primary_groups = qualifying_groups & primary_groups
                qualifying_secondary_groups = qualifying_groups - primary_groups
                primary_group_state = bool(qualifying_primary_groups)
                secondary_group_state = (
                    len(qualifying_secondary_groups) >= min_secondary_groups
                )
                primary_group_supported += int(primary_group_state)
                secondary_group_supported += int(secondary_group_state)
                output_row.update(
                    {
                        "primary_groups_ge_threshold": len(
                            qualifying_primary_groups
                        ),
                        "secondary_groups_ge_threshold": len(
                            qualifying_secondary_groups
                        ),
                        "total_groups_ge_threshold": len(qualifying_groups),
                        "primary_group_support": str(primary_group_state).lower(),
                        "secondary_group_support": str(
                            secondary_group_state
                        ).lower(),
                    }
                )
            writer.writerow(output_row)
    parameters: dict[str, Any] = {
        "primary_samples": sorted(primary_samples),
        "secondary_samples": sorted(samples - primary_samples),
        "min_reads_per_sample": min_reads_per_sample,
        "min_supporting_samples": min_supporting_samples,
        "strand_policy": "unstranded",
        "negative_evidence_policy": "absence_is_missing_not_contradiction",
    }
    counts: dict[str, Any] = {
        "samples": len(samples),
        "primary_samples": len(primary_samples),
        "secondary_samples": len(samples - primary_samples),
        "unique_junctions": len(support),
        "primary_supported_junctions": primary_supported,
        "secondary_supported_junctions": secondary_supported,
    }
    if sample_groups is not None:
        parameters.update(
            {
                "sample_group_tsv_sha256": _file_sha256(sample_group_tsv_path),
                "sample_group_interpretation": (
                    "filename_stem_only_not_biological_metadata"
                ),
                "min_samples_per_group": min_samples_per_group,
                "min_secondary_groups": min_secondary_groups,
            }
        )
        counts.update(
            {
                "filename_stem_groups": len(set(sample_groups.values())),
                "primary_group_supported_junctions": primary_group_supported,
                "secondary_group_supported_junctions": secondary_group_supported,
            }
        )
    manifest = {
        "schema_version": (
            JUNCTION_GROUP_AGGREGATE_SCHEMA
            if sample_groups is not None
            else JUNCTION_AGGREGATE_SCHEMA
        ),
        "generator": {"name": "PloidyPatch", "version": __version__},
        "inputs": sources,
        "parameters": parameters,
        "counts": counts,
        "output": {
            "file_name": output.name,
            "sha256": _file_sha256(output),
        },
    }
    _write_json_exclusive(manifest_path, manifest)
    return manifest
