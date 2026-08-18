from __future__ import annotations

import json
from pathlib import Path

from ploidypatch.artifact_manifest import sha256_file, verify_sha256sums


ROOT = Path(__file__).resolve().parents[1]
EVIDENCE = ROOT / "docs/release_evidence/wheel_core_v0.2"


def test_current_release_candidate_is_exact_and_installed_wheel_smoke_passed() -> None:
    assert verify_sha256sums(EVIDENCE, ignore_checksum_file=True)
    manifest = json.loads((EVIDENCE / "release_manifest.json").read_text(encoding="utf-8"))
    assert manifest["schema_version"] == "ploidypatch.release_candidate_evidence.v2"
    assert manifest["source_commit"] == "34988bdfa15a5f3bd107cf11329f26e842042e66"
    assert manifest["source_worktree_clean"] is True
    assert manifest["build"] == {
        "backend": "setuptools.build_meta",
        "frontend": "PyPA build",
        "isolation": True,
        "python": "3.11.15",
    }
    assert manifest["distributions"]["wheel"]["sha256"] == (
        "a5287972f6bf4b60d044b466b565cab3ad0b30063189180dbbb607c6fb027f53"
    )
    assert manifest["distributions"]["sdist"]["sha256"] == (
        "3ad3e3a2ef65e5ed33c120113d56d27564790129c6df124c4292e702b4cfd916"
    )
    assert manifest["installed_smoke"] == {
        "command_families": [
            "audit",
            "baseline",
            "benchmark",
            "evidence",
            "graph",
            "normalize",
            "patch",
        ],
        "example_accepted_additions": 2,
        "example_automatic_approval": False,
        "example_byte_identical_reversion": True,
        "pip_check": "passed",
        "source_tree_pythonpath_used": False,
        "version": "0.1.0a0",
    }
    assert manifest["release_boundary"]["formal_public_release"] is False


def test_release_evidence_has_sanitized_commands_and_exact_example_safety() -> None:
    commands = [
        json.loads(line)
        for line in (EVIDENCE / "command_audit.jsonl").read_text(encoding="utf-8").splitlines()
    ]
    assert [row["label"] for row in commands] == [
        "pep517_build",
        "install_wheel_no_deps",
        "pip_check",
        "installed_version",
        "installed_help",
        "installed_wheel_reviewed_patch_example",
    ]
    assert all(row["returncode"] == 0 for row in commands)
    assert all("argv" not in row for row in commands)
    evidence_text = "\n".join(
        path.read_text(encoding="utf-8")
        for path in EVIDENCE.iterdir()
        if path.is_file() and path.suffix != ".jsonl"
    )
    assert "D:\\" not in evidence_text
    assert "/data/codexli" not in evidence_text
    summary = json.loads((EVIDENCE / "example_run_summary.json").read_text(encoding="utf-8"))
    assert summary["accepted_additions"] == 2
    assert summary["automatic_approval"] is False
    assert summary["byte_identical_reversion"] is True
    assert summary["source_sha256"] == summary["reverted_sha256"]
    assert sha256_file(EVIDENCE / "SHA256SUMS") == (
        "8160101578dd8fa96863b33e13b626ebc2e89e40f523f4aa6923f41d3c3b1e97"
    )


def test_sdist_inventory_contains_the_release_builder_and_public_assets() -> None:
    members = set((EVIDENCE / "sdist_members.txt").read_text(encoding="utf-8").splitlines())
    prefix = "ploidypatch-0.1.0a0/"
    for relative in (
        "README.md",
        "docs/USER_GUIDE.md",
        "docs/REPRODUCIBILITY_GUIDE.md",
        "docs/CLI_COMMAND_INVENTORY_v0.1.json",
        "examples/minimal_reviewed_patch/run_example.py",
        "scripts/build_release_candidate_v0.2.py",
    ):
        assert prefix + relative in members
