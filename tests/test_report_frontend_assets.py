from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
ASSETS = ROOT / "src" / "ploidypatch" / "report_assets"


def _relative_luminance(hex_color: str) -> float:
    channels = [int(hex_color[index : index + 2], 16) / 255 for index in (1, 3, 5)]
    linear = [
        value / 12.92
        if value <= 0.04045
        else ((value + 0.055) / 1.055) ** 2.4
        for value in channels
    ]
    return 0.2126 * linear[0] + 0.7152 * linear[1] + 0.0722 * linear[2]


def _contrast_ratio(first: str, second: str) -> float:
    high, low = sorted(
        (_relative_luminance(first), _relative_luminance(second)), reverse=True
    )
    return (high + 0.05) / (low + 0.05)


def test_frontend_is_maintainable_and_packaged_as_three_assets() -> None:
    template = (ASSETS / "report.html").read_text(encoding="utf-8")
    css = (ASSETS / "report.css").read_text(encoding="utf-8")
    javascript = (ASSETS / "report.js").read_text(encoding="utf-8")
    source = (ROOT / "src" / "ploidypatch" / "review_report.py").read_text(
        encoding="utf-8"
    )
    assert template.count("__REPORT_TITLE__") == 2
    assert template.count("__REPORT_CSS__") == 1
    assert template.count("__REPORT_JS__") == 1
    assert template.count("__REPORT_PAYLOAD__") == 1
    assert "const RAW=" in javascript
    assert "const RAW=" not in source
    assert "template.replace(\"{{\"" not in source
    assert "prefers-color-scheme:dark" in css
    assert "prefers-reduced-motion:reduce" in css
    assert ":focus-visible" in css


def test_frontend_javascript_parses_when_node_is_available() -> None:
    node = shutil.which("node")
    if node is None:
        pytest.skip("Node.js is verified in GitHub Actions")
    completed = subprocess.run(
        [node, "--check", str(ASSETS / "report.js")],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr


def test_light_and_dark_report_tokens_meet_wcag_contrast_targets() -> None:
    text_pairs = (
        ("#132338", "#ffffff"),
        ("#607086", "#ffffff"),
        ("#2c78b8", "#ffffff"),
        ("#edf4fb", "#142235"),
        ("#b4c2d2", "#142235"),
        ("#75b9ef", "#142235"),
    )
    focus_pairs = (("#8b5a00", "#ffffff"), ("#f2c66d", "#142235"))
    assert all(_contrast_ratio(*pair) >= 4.5 for pair in text_pairs)
    assert all(_contrast_ratio(*pair) >= 3.0 for pair in focus_pairs)

    css = (ASSETS / "report.css").read_text(encoding="utf-8")
    assert "--focus:#8b5a00" in css
    assert "--focus:#f2c66d" in css
    assert "outline-color:var(--focus)" in css


def test_source_distribution_manifest_and_package_data_include_frontend() -> None:
    manifest = (ROOT / "MANIFEST.in").read_text(encoding="utf-8")
    pyproject = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
    assert "recursive-include src/ploidypatch/report_assets *.html *.css *.js" in manifest
    assert 'ploidypatch = ["report_assets/*.html", "report_assets/*.css", "report_assets/*.js"]' in pyproject
