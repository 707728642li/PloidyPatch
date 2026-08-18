from __future__ import annotations

import hashlib
import json
from pathlib import Path, PurePosixPath

from ploidypatch.artifact_manifest import verify_sha256sums


ROOT = Path(__file__).resolve().parents[1]
EVIDENCE = ROOT / "docs/release_evidence/wheel_core_v0.3"


def test_reproducible_release_evidence_is_exact_and_keeps_owner_boundary() -> None:
    assert verify_sha256sums(EVIDENCE, ignore_checksum_file=True)
    assert hashlib.sha256((EVIDENCE / "SHA256SUMS").read_bytes()).hexdigest() == (
        "45d53a8d0fe6b803e52918430c2dfeacc187d8545fd81c8423eacde430b5e4f0"
    )
    manifest = json.loads((EVIDENCE / "release_manifest.json").read_text(encoding="utf-8"))
    report = json.loads(
        (EVIDENCE / "reproducibility_report.json").read_text(encoding="utf-8")
    )
    assert manifest["schema_version"] == "ploidypatch.release_candidate_evidence.v3"
    assert manifest["source_commit"] == "e4d32e60a5ea16b09e565480cfb25a32cea6bde2"
    assert manifest["source_worktree_clean"] is True
    assert manifest["build"]["independent_builds"] == 2
    assert manifest["build"]["source_date_epoch"] == 1_786_321_076
    assert manifest["reproducibility"] == report
    assert report["wheel_byte_identical"] is True
    assert report["sdist_byte_identical"] is True
    assert report["stable_example_outputs_identical"] is True
    assert report["normalized_example_command_logs_identical"] is True
    assert manifest["release_boundary"] == {
        "formal_public_release": False,
        "missing_owner_fields_unchanged": [
            "authors",
            "license_spdx",
            "repository_url",
            "archive_doi",
        ],
    }


def test_reproducible_distributions_and_both_install_smokes_are_bound() -> None:
    manifest = json.loads((EVIDENCE / "release_manifest.json").read_text(encoding="utf-8"))
    wheel = manifest["distributions"]["wheel"]
    sdist = manifest["distributions"]["sdist"]
    assert wheel == {
        "bytes": 360_677,
        "members": 78,
        "name": "ploidypatch-0.1.0a0-py3-none-any.whl",
        "no_mandatory_runtime_dependencies": True,
        "sha256": "b14592beea9735091e26eb7184c8a38a846616398ac329eb2b20da92eb66a246",
    }
    assert sdist["bytes"] == 1_291_167
    assert sdist["members"] == 782
    assert sdist["sha256"] == (
        "d4be89e48d7768a9449fe14aa0cf0ac0be32856dd5184dff7b26e35d32c0cf01"
    )
    assert len((EVIDENCE / "wheel_members.txt").read_text(encoding="utf-8").splitlines()) == 78
    assert len((EVIDENCE / "sdist_members.txt").read_text(encoding="utf-8").splitlines()) == 782
    for kind in ("wheel", "sdist"):
        smoke = manifest["installed_smoke"][kind]
        assert smoke["pip_check"] == "passed"
        assert smoke["source_tree_pythonpath_used"] is False
        assert smoke["version"] == "0.1.0a0"
        assert smoke["example_run_summary_sha256"] == (
            "d36f2c5ff0a983b7780749281d5166062f1706ab0111d5094c90e75dae61f795"
        )
    assert manifest["installed_smoke"]["example_accepted_additions"] == 2
    assert manifest["installed_smoke"]["example_automatic_approval"] is False
    assert manifest["installed_smoke"]["example_byte_identical_reversion"] is True


def test_reproducible_release_command_audit_and_member_paths_are_safe() -> None:
    records = [
        json.loads(line)
        for line in (EVIDENCE / "command_audit.jsonl").read_text(encoding="utf-8").splitlines()
    ]
    assert [record["label"] for record in records] == [
        "build_a_pep517_build",
        "build_b_pep517_build",
        "install_wheel_no_deps",
        "wheel_pip_check",
        "wheel_installed_version",
        "wheel_installed_help",
        "installed_wheel_reviewed_patch_example",
        "install_sdist_no_deps",
        "sdist_pip_check",
        "sdist_installed_version",
        "sdist_installed_help",
        "installed_sdist_reviewed_patch_example",
    ]
    assert all(record["returncode"] == 0 for record in records)
    for file_name in ("wheel_members.txt", "sdist_members.txt"):
        names = (EVIDENCE / file_name).read_text(encoding="utf-8").splitlines()
        assert names == sorted(names, key=lambda value: value.encode("utf-8"))
        assert len(names) == len(set(name.casefold() for name in names))
        for name in names:
            path = PurePosixPath(name)
            assert not path.is_absolute()
            assert all(part not in {"", ".", ".."} for part in path.parts)
