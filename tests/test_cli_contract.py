from __future__ import annotations

import argparse
import tomllib
from pathlib import Path

import pytest

from ploidypatch import __version__
from ploidypatch.cli import build_parser


ROOT = Path(__file__).parents[1]


def top_level_commands(parser: argparse.ArgumentParser) -> set[str]:
    actions = [
        action
        for action in parser._actions
        if isinstance(action, argparse._SubParsersAction)
    ]
    assert len(actions) == 1
    return set(actions[0].choices)


def test_distribution_and_runtime_versions_are_identical() -> None:
    metadata = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    assert metadata["project"]["version"] == __version__


def test_public_cli_keeps_documented_command_families() -> None:
    parser = build_parser()
    assert top_level_commands(parser) == {
        "audit",
        "baseline",
        "benchmark",
        "evidence",
        "graph",
        "normalize",
        "patch",
        "report",
    }
    help_text = parser.format_help()
    assert "usage: ploidypatch" in help_text
    for command in top_level_commands(parser):
        assert command in help_text


def test_cli_version_action_reports_runtime_version(capsys: pytest.CaptureFixture[str]) -> None:
    with pytest.raises(SystemExit) as caught:
        build_parser().parse_args(["--version"])
    assert caught.value.code == 0
    assert capsys.readouterr().out.strip() == __version__
