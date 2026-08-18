from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
STATUS = ROOT / "docs/release_evidence/public_v1.0_external_status_2026-08-11.json"
ZENODO = ROOT / "docs/release_evidence/zenodo_v1.0_record_2026-08-11.json"


def test_public_v1_external_status_is_complete_and_claim_bounded() -> None:
    value = json.loads(STATUS.read_text(encoding="utf-8"))
    assert value["schema_version"] == "ploidypatch.public_release_external_status.v1"
    assert value["repository"] == {
        "full_name": "707728642li/PloidyPatch",
        "private": False,
        "default_branch": "main",
        "default_branch_head_sha": "729de93a639a905a6de3888a47c1a327dbee9544",
        "license_spdx": "BSD-3-Clause",
        "url": "https://github.com/707728642li/PloidyPatch",
    }
    release = value["release"]
    assert release["tag"] == "v1.0.0"
    assert release["draft"] is False
    assert release["prerelease"] is False
    assert {item["name"] for item in release["assets"]} == {
        "ploidypatch-1.0.0-py3-none-any.whl",
        "ploidypatch-1.0.0.tar.gz",
        "SHA256SUMS",
    }
    assert all(len(item["sha256"]) == 64 for item in release["assets"])
    zenodo = value["zenodo"]
    assert zenodo == {
        "record_id": 21875561,
        "doi": "10.5281/zenodo.21875561",
        "status": "published",
        "state": "done",
        "version": "1.0.0",
        "resource_type": "software",
        "access_right": "open",
        "doi_http_status": 200,
        "doi_final_url": "https://zenodo.org/records/21875561",
        "record_evidence": "docs/release_evidence/zenodo_v1.0_record_2026-08-11.json",
    }
    assert {item["name"]: item["conclusion"] for item in value["upstream_ci"]["checks"]} == {
        "test": "success",
        "codeql": "success",
    }
    bioconda = value["bioconda"]
    assert bioconda["state"] == "open"
    assert bioconda["mergeable"] is True
    assert bioconda["mergeable_state"] == "clean"
    assert bioconda["label"] == "please review & merge"
    assert bioconda["review_count"] == 0
    assert bioconda["actionable_review_comments"] == 0
    assert bioconda["installability_claim"] == "pending_merge"
    assert {item["conclusion"] for item in bioconda["checks"]} == {
        "success",
        "neutral",
    }
    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    installation = (ROOT / "docs/INSTALLATION.md").read_text(encoding="utf-8")
    pyproject = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
    assert "Production annotation and validation workflows are supported on **Linux**" in readme
    assert "Native Windows is not an officially supported bioinformatics platform" in readme
    assert "native Windows is not supported" in installation
    assert "Use WSL2, a Linux server, or the container" in installation
    assert '"Operating System :: POSIX :: Linux"' in pyproject
    assert '"Operating System :: MacOS"' in pyproject
    assert "Operating System :: OS Independent" not in pyproject


def test_zenodo_record_is_public_and_byte_identical_to_release_assets() -> None:
    evidence = json.loads(ZENODO.read_text(encoding="utf-8"))
    record = evidence["record"]
    assert record["doi"] == "10.5281/zenodo.21875561"
    assert record["status"] == "published"
    assert record["state"] == "done"
    assert record["submitted"] is True
    assert record["version"] == "1.0.0"
    assert record["resource_type"] == "software"
    assert record["access_right"] == "open"
    assert evidence["doi_resolution"] == {
        "http_status": 200,
        "final_url": "https://zenodo.org/records/21875561",
    }
    status = json.loads(STATUS.read_text(encoding="utf-8"))
    github_assets = {
        item["name"]: (item["bytes"], item["sha256"])
        for item in status["release"]["assets"]
    }
    zenodo_assets = {
        item["name"]: (item["bytes"], item["sha256"])
        for item in evidence["files"]
    }
    assert zenodo_assets == github_assets
    assert all(len(item["zenodo_md5"]) == 32 for item in evidence["files"])
    assert 'doi: "10.5281/zenodo.21875561"' in (
        ROOT / "CITATION.cff"
    ).read_text(encoding="utf-8")
