from __future__ import annotations

import csv
import gzip
import importlib.util
import stat
from pathlib import Path

import pytest


SCRIPT = (
    Path(__file__).parents[1]
    / "scripts"
    / "verify_apple_flnc_result_v0.1.py"
)


def load_module():
    spec = importlib.util.spec_from_file_location("verify_apple_flnc", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_fastq_counter_is_strict(tmp_path: Path) -> None:
    module = load_module()
    valid = tmp_path / "valid.fastq.gz"
    with gzip.open(valid, "wt") as handle:
        handle.write("@a\nACGT\n+\nIIII\n@b\nA\n+\nI\n")
    assert module.count_fastq_records(valid) == 2
    truncated = tmp_path / "truncated.fastq.gz"
    with gzip.open(truncated, "wt") as handle:
        handle.write("@a\nACGT\n+\n")
    with pytest.raises(module.VerificationError, match="truncated"):
        module.count_fastq_records(truncated)


def test_tsv_reader_accepts_frozen_long_read_support_fields(tmp_path: Path) -> None:
    module = load_module()
    path = tmp_path / "evidence.tsv"
    oversized = "read," * 30_000
    path.write_text(f"candidate_digest\tsupporting_reads\nabc\t{oversized}\n", encoding="utf-8")
    previous_limit = csv.field_size_limit()
    try:
        csv.field_size_limit(131_072)
        rows = list(module.read_tsv(path))
    finally:
        csv.field_size_limit(previous_limit)
    assert rows == [{"candidate_digest": "abc", "supporting_reads": oversized}]
    assert len(oversized) > 131_072


def test_checksum_manifest_requires_exact_safe_file_universe(tmp_path: Path) -> None:
    module = load_module()
    payload = tmp_path / "payload.txt"
    payload.write_text("payload\n", encoding="utf-8")
    digest = module.sha256_file(payload)
    (tmp_path / "SHA256SUMS").write_text(
        f"{digest}  ./payload.txt\n", encoding="utf-8"
    )
    assert module.verify_root_manifest(tmp_path, verify_hashes=True)[0] == 1
    extra = tmp_path / "extra.txt"
    extra.write_text("extra\n", encoding="utf-8")
    with pytest.raises(module.VerificationError, match="file universe"):
        module.verify_root_manifest(tmp_path, verify_hashes=True)
    extra.unlink()
    (tmp_path / "SHA256SUMS").write_text(
        f"{digest}  ../payload.txt\n", encoding="utf-8"
    )
    with pytest.raises(module.VerificationError, match="traversal"):
        module.verify_root_manifest(tmp_path, verify_hashes=False)


def test_read_only_verifier_checks_every_path(tmp_path: Path) -> None:
    module = load_module()
    child = tmp_path / "child.txt"
    child.write_text("x", encoding="utf-8")
    child.chmod(stat.S_IRUSR | stat.S_IRGRP | stat.S_IROTH)
    tmp_path.chmod(stat.S_IRUSR | stat.S_IXUSR | stat.S_IRGRP | stat.S_IXGRP | stat.S_IROTH | stat.S_IXOTH)
    module.verify_read_only_tree(tmp_path)
    tmp_path.chmod(stat.S_IRWXU)
    child.chmod(stat.S_IRUSR | stat.S_IWUSR)
    with pytest.raises(module.VerificationError, match="write bits"):
        module.verify_read_only_tree(tmp_path)


def test_verifier_freezes_expected_scientific_contract() -> None:
    text = SCRIPT.read_text(encoding="utf-8")
    for literal in (
        "29144",
        "20261004",
        "20261005",
        "query_orientation",
        "full_chain_supported",
        "case_study_ready",
        "automatic_annotation_patch",
        "SHA256SUMS file universe",
        '"verifier"',
    ):
        assert literal in text
