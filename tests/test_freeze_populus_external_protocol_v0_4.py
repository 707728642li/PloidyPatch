from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest


MODULE_PATH = (
    Path(__file__).parents[1]
    / "scripts"
    / "freeze_populus_external_protocol_v0.4.py"
)
SPEC = importlib.util.spec_from_file_location("populus_protocol_freeze", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_policy_reader_rejects_duplicate_fields(tmp_path: Path) -> None:
    policy = tmp_path / "policy.tsv"
    policy.write_text("field\tvalue\na\t1\na\t2\n", encoding="utf-8")
    with pytest.raises(ValueError, match="duplicate"):
        MODULE.read_policy(policy)


def test_sha256_verifier_detects_mutation(tmp_path: Path) -> None:
    artifact = tmp_path / "artifact.txt"
    artifact.write_text("frozen\n", encoding="utf-8")
    (tmp_path / "SHA256SUMS").write_text(
        f"{MODULE.sha256(artifact)}  artifact.txt\n", encoding="utf-8"
    )
    MODULE.verify_sha256sums(tmp_path)
    artifact.write_text("changed\n", encoding="utf-8")
    with pytest.raises(ValueError, match="Checksum failure"):
        MODULE.verify_sha256sums(tmp_path)


def test_freeze_schema_declares_pre_pair_stage() -> None:
    assert MODULE.SCHEMA_VERSION.endswith("v0.4")
    assert MODULE.POLICY_ID == "ploidypatch_populus_external_validation_v0.4"
