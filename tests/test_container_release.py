from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_container_uses_digest_pinned_two_stage_offline_build() -> None:
    dockerfile = (ROOT / "containers" / "Dockerfile").read_text(encoding="utf-8")
    digest = "sha256:d29f48a31a8b408ed19272ca1e7b10ebae13b240a27e862d3d4217c528e2e0c3"
    assert f"python:3.11.15-slim-bookworm@{digest}" in dockerfile
    assert "AS wheel-builder" in dockerfile
    assert "AS runtime" in dockerfile
    assert dockerfile.count("RUN --network=none") == 2
    assert "--no-build-isolation" in dockerfile
    assert dockerfile.count("--no-deps") == 2
    assert dockerfile.count("ploidypatch-1.0.1-py3-none-any.whl") == 3
    assert 'org.opencontainers.image.licenses="BSD-3-Clause"' in dockerfile
    assert "/tmp/ploidypatch.whl" not in dockerfile
    assert "USER 10001:10001" in dockerfile
    assert "mkdir -p /work/examples" in dockerfile
    assert 'ENTRYPOINT ["ploidypatch"]' in dockerfile
    assert "apt-get" not in dockerfile


def test_container_build_context_is_an_allowlist() -> None:
    lines = (ROOT / ".dockerignore").read_text(encoding="utf-8").splitlines()
    assert lines[0] == "**"
    assert set(lines[1:]) == {
        "!pyproject.toml",
        "!README.md",
        "!LICENSE",
        "!src/",
        "!src/**",
        "!containers/",
        "!containers/Dockerfile",
    }
    for forbidden in ("results", "data", "manuscript", "envs", ".git"):
        assert not any(line == f"!{forbidden}/" for line in lines)


def test_container_ci_runs_nonroot_offline_readonly_example() -> None:
    workflow = (ROOT / ".github" / "workflows" / "test.yml").read_text(
        encoding="utf-8"
    )
    smoke = (ROOT / "scripts" / "smoke_container_v0.1.sh").read_text(
        encoding="utf-8"
    )
    assert "Build digest-pinned non-root core image" in workflow
    assert "bash scripts/smoke_container_v0.1.sh" in workflow
    assert 'assert report["runtime_user"] == "10001:10001"' in workflow
    assert "--read-only" in smoke
    assert "--network none" in smoke
    assert "--tmpfs /tmp:rw,noexec,nosuid,size=16m" in smoke
    assert "/work/examples/minimal_reviewed_patch/run_example.py" in smoke
    assert "dst=/work/examples/minimal_reviewed_patch,readonly" in smoke
    assert 'example.get("automatic_approval") is False' in smoke
    assert 'example.get("source_sha256") == example.get("reverted_sha256")' in smoke
    assert "refusing to overwrite container smoke output" in smoke
    assert "find . -type f ! -name SHA256SUMS -printf '%P\\0'" in smoke


def test_container_documentation_preserves_scope_and_license_boundary() -> None:
    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    guide = (ROOT / "docs" / "CONTAINER_GUIDE.md").read_text(encoding="utf-8")
    assert "docs/CONTAINER_GUIDE.md" in readme
    assert "does not bundle miniprot" in guide
    assert "not an automatic annotation service" in guide
    assert "`BSD-3-Clause` license" in guide
    assert "10001:10001" in guide
    assert "--network none" in guide
    assert "automatic_approval=false" in guide


def test_container_smoke_evidence_matches_upstream_manifest() -> None:
    root = ROOT / "docs" / "release_evidence" / "container_core_v0.1"
    rows = []
    for line in (root / "upstream_SHA256SUMS").read_text(encoding="utf-8").splitlines():
        digest, name = line.split("  ", 1)
        rows.append((digest, name))
    assert {name for _, name in rows} == {
        "cli_help.txt",
        "container_smoke.json",
        "example_stdout.json",
        "image_inspect.json",
        "version.txt",
    }
    for expected, name in rows:
        assert hashlib.sha256((root / name).read_bytes()).hexdigest() == expected

    report = json.loads((root / "container_smoke.json").read_text(encoding="utf-8"))
    assert report["schema_version"] == "ploidypatch.container_smoke.v1"
    assert report["image_id"] == (
        "sha256:d4457eb7c89e25d674bb53f26b78f483ad072029e26c1dc47431737ab811bbe7"
    )
    assert report["image_size_bytes"] == 133_761_776
    assert report["runtime_user"] == "10001:10001"
    assert report["runtime_network"] == "none"
    assert report["runtime_read_only"] is True
    assert all(report["checks"].values())


def test_container_evidence_package_is_exact_and_distributable() -> None:
    root = ROOT / "docs" / "release_evidence" / "container_core_v0.1"
    manifest = root / "package_SHA256SUMS"
    rows = []
    for line in manifest.read_text(encoding="utf-8").splitlines():
        digest, name = line.split("  ", 1)
        rows.append((digest, name))
    expected_names = {
        path.name for path in root.iterdir() if path.is_file() and path != manifest
    }
    assert {name for _, name in rows} == expected_names
    for expected, name in rows:
        assert hashlib.sha256((root / name).read_bytes()).hexdigest() == expected

    readme = (root / "README.md").read_text(encoding="utf-8")
    assert "147eb7873ae25ea094c5c9c1ee70d9a3a4ff3a17860f612298b1f8310a3a175c" in readme
    assert "b8400509a622dccde6ed23c9fe7c31d94886eea64154231a466143ff88bbcd7a" in readme
    assert "registry" in readme
    source_manifest = (ROOT / "MANIFEST.in").read_text(encoding="utf-8")
    assert "include .dockerignore" in source_manifest
    assert "recursive-include containers *" in source_manifest
    assert "recursive-include docs/release_evidence *" in source_manifest
    workflow = (ROOT / ".github" / "workflows" / "test.yml").read_text(
        encoding="utf-8"
    )
    for path in (
        "/containers/Dockerfile$",
        "/scripts/smoke_container_v0.1.sh$",
        "/docs/release_evidence/container_core_v0.1/container_smoke.json$",
    ):
        assert path in workflow


def test_container_smoke_audit_binds_implementation_and_image() -> None:
    audit = (ROOT / "docs" / "CONTAINER_SMOKE_2026-08-10.md").read_text(
        encoding="utf-8"
    )
    commit = "f75c35a343ad09418ac75948e9402f0d49a166e8"
    image = "sha256:d4457eb7c89e25d674bb53f26b78f483ad072029e26c1dc47431737ab811bbe7"
    package = "240aa19b2f4a1deaad700600e5b22c17063dd93269bec141dd7d521bfb52e21a"
    assert "10001:10001" in audit
    assert "read-only" in audit
    assert "automatic" in audit
    assert commit in audit
    assert image in audit
    assert package in audit
    assert "Build attempt 1" in audit
    assert "smoke attempt 1" in audit
