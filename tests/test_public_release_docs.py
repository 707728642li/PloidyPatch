from __future__ import annotations

import json
import hashlib
import os
import re
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
INVENTORY = ROOT / "docs" / "CLI_COMMAND_INVENTORY_v0.1.json"
EXPORTER = ROOT / "scripts" / "export_cli_inventory_v0.1.py"


def _environment() -> dict[str, str]:
    environment = os.environ.copy()
    source_root = os.fspath(ROOT / "src")
    prior = environment.get("PYTHONPATH", "")
    environment["PYTHONPATH"] = source_root + (os.pathsep + prior if prior else "")
    return environment


def test_tracked_cli_inventory_matches_production_parser() -> None:
    completed = subprocess.run(
        [sys.executable, os.fspath(EXPORTER), "--output", os.fspath(INVENTORY), "--check"],
        cwd=ROOT,
        env=_environment(),
        text=True,
        capture_output=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr
    inventory = json.loads(INVENTORY.read_text(encoding="utf-8"))
    assert inventory["schema_version"] == "ploidypatch.cli_command_inventory.v1"
    assert inventory["package_version"] == "1.0.1"
    assert inventory["top_level_families"] == [
        "audit",
        "baseline",
        "benchmark",
        "evidence",
        "graph",
        "normalize",
        "patch",
        "report",
    ]
    assert inventory["leaf_command_count"] == 87
    assert len(inventory["leaf_commands"]) == len(set(inventory["leaf_commands"]))
    for command in (
        "audit",
        "patch compile-reviewed-copy-additions",
        "patch create",
        "patch apply",
        "patch revert",
        "report",
    ):
        assert command in inventory["leaf_commands"]


def test_exporter_writes_deterministic_inventory(tmp_path: Path) -> None:
    generated = tmp_path / "inventory.json"
    completed = subprocess.run(
        [sys.executable, os.fspath(EXPORTER), "--output", os.fspath(generated)],
        cwd=ROOT,
        env=_environment(),
        text=True,
        capture_output=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr
    assert generated.read_bytes() == INVENTORY.read_bytes()


def test_public_docs_define_review_and_reproducibility_boundaries() -> None:
    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    user_guide = (ROOT / "docs" / "USER_GUIDE.md").read_text(encoding="utf-8")
    cli_reference = (ROOT / "docs" / "CLI_REFERENCE.md").read_text(encoding="utf-8")
    reproducibility = (ROOT / "docs" / "REPRODUCIBILITY_GUIDE.md").read_text(
        encoding="utf-8"
    )
    for link in (
        "docs/INSTALLATION.md",
        "docs/TUTORIAL.md",
        "docs/USER_GUIDE.md",
        "docs/CLI_REFERENCE.md",
        "docs/CLI_PARAMETERS.md",
        "docs/REPRODUCIBILITY_GUIDE.md",
        "docs/CLI_COMMAND_INVENTORY_v0.1.json",
    ):
        assert link in readme
    assert "automatic_approval=false" in readme
    assert "byte-identical" in readme
    assert "compile-reviewed-copy-additions" in user_guide
    assert "missing review decision: no automatic acceptance" in user_guide
    assert "87 executable leaf" in cli_reference
    assert "invalid_run" in reproducibility
    assert "not_evaluable" in reproducibility
    assert "formal negative" in reproducibility
    assert "formal positive" in reproducibility
    assert "87 executable leaf commands" in readme


def test_readme_is_a_concise_public_entrypoint() -> None:
    readme_path = ROOT / "README.md"
    readme = readme_path.read_text(encoding="utf-8")
    assert readme_path.stat().st_size < 20_000
    assert "/nas_data/" not in readme
    assert "/data/codexli/" not in readme
    assert "PloidyPatch contributors" not in readme
    assert "docs/examples/report_preview/index.html" in readme
    assert ".github/assets/report-workbench.svg" in readme


def test_github_actions_are_immutable_and_cover_report_quality() -> None:
    workflow_root = ROOT / ".github" / "workflows"
    workflows = {
        path.name: path.read_text(encoding="utf-8")
        for path in sorted(workflow_root.glob("*.yml"))
    }
    action_references = []
    for text in workflows.values():
        action_references.extend(re.findall(r"^\s*- uses:\s+([^\s#]+)", text, re.MULTILINE))
    assert action_references
    assert all(re.fullmatch(r"[^@\s]+@[0-9a-f]{40}", item) for item in action_references)
    assert "javascript-typescript" in workflows["codeql.yml"]
    test_workflow = workflows["test.yml"]
    for token in (
        "runs-on: macos-latest",
        "node --check src/ploidypatch/report_assets/report.js",
        "python -m ruff check",
        "python -m mypy",
        "--cov-fail-under=80",
    ):
        assert token in test_workflow


def test_public_report_preview_is_complete_and_checksum_bound() -> None:
    preview = ROOT / "docs/examples/report_preview"
    expected = {
        "index.html",
        "report.json",
        "candidates.tsv",
        "review_ledger_template.tsv",
    }
    recorded: dict[str, str] = {}
    for line in (preview / "SHA256SUMS").read_text(encoding="utf-8").splitlines():
        digest, relative = line.split("  ", 1)
        assert relative not in recorded
        recorded[relative] = digest
    assert set(recorded) == expected
    for relative, digest in recorded.items():
        assert hashlib.sha256((preview / relative).read_bytes()).hexdigest() == digest
    report = json.loads((preview / "report.json").read_text(encoding="utf-8"))
    assert report["summary"]["state"] == "validated_reversible_run"
    assert report["truth_access"] is False
    assert report["counts"]["automatic_approval_checked"] is True


def test_source_distribution_manifest_includes_public_docs() -> None:
    manifest = (ROOT / "MANIFEST.in").read_text(encoding="utf-8")
    workflow = (ROOT / ".github" / "workflows" / "test.yml").read_text(
        encoding="utf-8"
    )
    assert "recursive-include docs *.md *.json *.html *.tsv *.png *.svg *.pdf SHA256SUMS" in manifest
    for path in (
        "/docs/USER_GUIDE.md$",
        "/docs/REPRODUCIBILITY_GUIDE.md$",
        "/docs/CLI_COMMAND_INVENTORY_v0.1.json$",
    ):
        assert path in workflow


def test_historical_release_evidence_is_labeled_and_internally_bound() -> None:
    release = (ROOT / "docs" / "RELEASE_WHEEL_SMOKE_2026-08-08.md").read_text(
        encoding="utf-8"
    )
    commit = "e4d32e60a5ea16b09e565480cfb25a32cea6bde2"
    wheel_sha = "b14592beea9735091e26eb7184c8a38a846616398ac329eb2b20da92eb66a246"
    sdist_sha = "d4be89e48d7768a9449fe14aa0cf0ac0be32856dd5184dff7b26e35d32c0cf01"
    assert commit in release
    assert wheel_sha in release
    assert "no source-tree `PYTHONPATH`" in release
    assert sdist_sha in release
    timeline = (ROOT / "docs/release_evidence/README.md").read_text(encoding="utf-8")
    assert "development `0.1.0a0`" in timeline
    assert "deliberately not rewritten" in timeline
    assert "eight command families" in " ".join(timeline.split())
