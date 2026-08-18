from __future__ import annotations

import csv
import io
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from ploidypatch.cli import build_parser
from ploidypatch.perturb import _file_sha256
from ploidypatch.rna import (
    JUNCTION_MANIFEST_SCHEMA,
    aggregate_junctions,
    extract_bam_junctions,
    junctions_from_cigar,
)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_junction_file(
    directory: Path,
    sample: str,
    rows: list[tuple[str, int, int, int]],
) -> Path:
    path = directory / f"{sample}.junctions.tsv"
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            ["sample", "seqid", "left_exon_end", "right_exon_start", "read_count"]
        )
        for seqid, left, right, count in rows:
            writer.writerow([sample, seqid, left, right, count])
    manifest = {
        "schema_version": JUNCTION_MANIFEST_SCHEMA,
        "parameters": {"sample": sample},
        "output": {"sha256": _file_sha256(path)},
    }
    Path(str(path) + ".manifest.json").write_text(
        json.dumps(manifest), encoding="utf-8", newline=""
    )
    return path


def test_junctions_from_cigar_uses_one_based_exon_boundaries() -> None:
    assert junctions_from_cigar("chr1", 100, "10M90N20M") == (
        ("chr1", 109, 200),
    )
    assert junctions_from_cigar("chr2", 1, "5S10M2I5M100N10M3D5M5S") == (
        ("chr2", 15, 116),
    )
    assert junctions_from_cigar("chr1", 100, "90N20M") == ()
    assert junctions_from_cigar("chr1", 100, "20M90N") == ()
    assert junctions_from_cigar("chr1", 0, "10M90N20M") == ()
    with pytest.raises(ValueError, match="Malformed SAM CIGAR"):
        junctions_from_cigar("chr1", 100, "10Mbroken")


def test_extract_bam_junctions_records_filters_and_manifest(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    bam = tmp_path / "Black-1.sorted.bam"
    bai = tmp_path / "Black-1.sorted.bam.bai"
    samtools = tmp_path / "samtools"
    for path in (bam, bai, samtools):
        path.write_bytes(b"fixture")

    def fake_run(*args: object, **kwargs: object) -> SimpleNamespace:
        return SimpleNamespace(stdout="@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n")

    class FakeProcess:
        def __init__(self, command: list[str], **kwargs: object) -> None:
            self.command = command
            self.stdout = io.StringIO(
                "r1\t0\tchr1\t100\t60\t10M90N20M\t*\t0\t0\t*\t*\n"
                "r2\t0\tchr1\t100\t60\t10M90N20M\t*\t0\t0\t*\t*\n"
                "r3\t0\tchr1\t300\t60\t20M\t*\t0\t0\t*\t*\n"
            )
            self.return_code = 0

        def wait(self) -> int:
            return self.return_code

        def kill(self) -> None:
            self.return_code = -9

    monkeypatch.setattr("ploidypatch.rna.subprocess.run", fake_run)
    monkeypatch.setattr("ploidypatch.rna.subprocess.Popen", FakeProcess)
    output = tmp_path / "Black-1.junctions.tsv"
    manifest = extract_bam_junctions(
        bam,
        output,
        sample="Black-1",
        samtools_path=samtools,
        threads=3,
        min_mapq=25,
    )

    assert read_tsv(output) == [
        {
            "sample": "Black-1",
            "seqid": "chr1",
            "left_exon_end": "109",
            "right_exon_start": "200",
            "read_count": "2",
        }
    ]
    assert manifest["counts"] == {
        "alignments_passing_samtools_filters": 3,
        "split_alignments": 2,
        "junction_observations": 2,
        "unique_junctions": 1,
    }
    assert manifest["command"][2:8] == ["-@", "3", "-q", "25", "-F", "2308"]
    assert manifest["parameters"]["strand_policy"].startswith("unstranded")
    with pytest.raises(FileExistsError, match="Refusing to overwrite"):
        extract_bam_junctions(
            bam,
            output,
            sample="Black-1",
            samtools_path=samtools,
        )


def test_aggregate_junctions_separates_primary_and_secondary_tiers(
    tmp_path: Path,
) -> None:
    paths = [
        write_junction_file(
            tmp_path, "Black-1", [("chr1", 109, 200, 3), ("chr1", 309, 400, 1)]
        ),
        write_junction_file(
            tmp_path, "Black-2", [("chr1", 109, 200, 2), ("chr1", 309, 400, 2)]
        ),
        write_junction_file(
            tmp_path, "Black-3", [("chr1", 109, 200, 4), ("chr1", 309, 400, 2)]
        ),
        write_junction_file(
            tmp_path, "Red-1", [("chr1", 309, 400, 5), ("chr2", 9, 100, 3)]
        ),
        write_junction_file(
            tmp_path, "Red-2", [("chr1", 309, 400, 3), ("chr2", 9, 100, 1)]
        ),
    ]
    output = tmp_path / "aggregate.tsv"
    manifest = aggregate_junctions(
        paths,
        output,
        primary_samples=frozenset({"Black-1", "Black-2", "Black-3"}),
    )
    rows = {(row["seqid"], row["left_exon_end"]): row for row in read_tsv(output)}

    assert rows[("chr1", "109")]["primary_support"] == "true"
    assert rows[("chr1", "109")]["secondary_support"] == "false"
    assert rows[("chr1", "309")]["primary_samples_ge_threshold"] == "2"
    assert rows[("chr1", "309")]["secondary_samples_ge_threshold"] == "2"
    assert rows[("chr1", "309")]["primary_read_count"] == "5"
    assert rows[("chr2", "9")]["secondary_support"] == "false"
    assert manifest["counts"] == {
        "samples": 5,
        "primary_samples": 3,
        "secondary_samples": 2,
        "unique_junctions": 3,
        "primary_supported_junctions": 2,
        "secondary_supported_junctions": 1,
    }
    assert manifest["parameters"]["negative_evidence_policy"] == (
        "absence_is_missing_not_contradiction"
    )

    groups = tmp_path / "groups.tsv"
    groups.write_text(
        "sample\tfilename_stem_group\n"
        "Black-1\tBlack\nBlack-2\tBlack\nBlack-3\tBlack\n"
        "Red-1\tRed\nRed-2\tRed\n",
        encoding="utf-8",
        newline="",
    )
    grouped_output = tmp_path / "grouped.tsv"
    grouped_manifest = aggregate_junctions(
        paths,
        grouped_output,
        primary_samples=frozenset({"Black-1", "Black-2", "Black-3"}),
        sample_group_tsv_path=groups,
        min_secondary_groups=1,
    )
    grouped_rows = {
        (row["seqid"], row["left_exon_end"]): row
        for row in read_tsv(grouped_output)
    }
    assert grouped_rows[("chr1", "309")]["primary_group_support"] == "true"
    assert grouped_rows[("chr1", "309")]["secondary_group_support"] == "true"
    assert grouped_rows[("chr1", "309")]["secondary_groups_ge_threshold"] == "1"
    assert grouped_manifest["counts"]["filename_stem_groups"] == 2
    assert grouped_manifest["parameters"]["sample_group_interpretation"] == (
        "filename_stem_only_not_biological_metadata"
    )


def test_aggregate_checks_primary_presence_hashes_and_cli() -> None:
    parser = build_parser()
    parsed = parser.parse_args(
        [
            "evidence",
            "aggregate-junctions",
            "--input-dir",
            "junctions",
            "--primary-sample",
            "Black-1",
            "--output-tsv",
            "aggregate.tsv",
        ]
    )
    assert parsed.min_reads_per_sample == 2
    assert parsed.min_supporting_samples == 2
    assert parsed.primary_sample == ["Black-1"]
    assert parsed.input_dir == ["junctions"]
    assert parsed.sample_groups is None
    assert parsed.min_samples_per_group == 2
    assert parsed.min_secondary_groups == 2
