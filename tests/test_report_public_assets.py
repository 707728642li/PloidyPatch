from __future__ import annotations

from pathlib import Path
import subprocess


ROOT = Path(__file__).resolve().parents[1]
REQUIRED = (
    "docs/RESULT_REPORT.md",
    "docs/examples/report_preview/SHA256SUMS",
    "docs/examples/report_preview/candidates.tsv",
    "docs/examples/report_preview/index.html",
    "docs/examples/report_preview/report.json",
    "docs/examples/report_preview/review_ledger_template.tsv",
    "examples/minimal_reviewed_patch/report_copy_features.tsv",
    "examples/minimal_reviewed_patch/report_scores.tsv",
    "examples/minimal_reviewed_patch/report_topology.tsv",
)


def test_documented_report_assets_exist_and_are_git_tracked() -> None:
    for relative in REQUIRED:
        path = ROOT / relative
        assert path.is_file() and not path.is_symlink(), relative
    completed = subprocess.run(
        ["git", "ls-files", "--error-unmatch", "--", *REQUIRED],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr


def test_result_report_documentation_uses_the_verified_public_paths() -> None:
    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    guide = (ROOT / "docs/RESULT_REPORT.md").read_text(encoding="utf-8")
    assert "docs/examples/report_preview/index.html" in readme
    for relative in REQUIRED[6:]:
        assert relative in guide
