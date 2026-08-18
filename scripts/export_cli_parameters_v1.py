#!/usr/bin/env python3
"""Export deterministic, human-readable parameter tables from the production CLI."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from ploidypatch import __version__
from ploidypatch.cli import build_parser


def _leaf_parsers(
    parser: argparse.ArgumentParser,
    prefix: tuple[str, ...] = (),
) -> list[tuple[tuple[str, ...], argparse.ArgumentParser]]:
    children: list[tuple[tuple[str, ...], argparse.ArgumentParser]] = []
    subparser_actions = [
        action
        for action in parser._actions
        if isinstance(action, argparse._SubParsersAction)
    ]
    if not subparser_actions:
        return [(prefix, parser)]
    for action in subparser_actions:
        for name in sorted(action.choices):
            children.extend(_leaf_parsers(action.choices[name], (*prefix, name)))
    return children


def _display_default(action: argparse.Action) -> str:
    if action.default is argparse.SUPPRESS:
        return "—"
    if action.default is None:
        return "None"
    if isinstance(action.default, bool):
        return str(action.default).lower()
    if isinstance(action.default, (str, int, float)):
        return str(action.default)
    return json.dumps(action.default, sort_keys=True, default=str)


def _display_choices(action: argparse.Action) -> str:
    if action.choices is None:
        return "—"
    return ", ".join(f"`{choice}`" for choice in action.choices)


def _escape(value: Any) -> str:
    return str(value).replace("|", "\\|").replace("\n", " ").strip()


def _argument_label(action: argparse.Action) -> str:
    if action.option_strings:
        return "<br>".join(f"`{value}`" for value in action.option_strings)
    return f"`{action.dest}`"


def _argument_type(action: argparse.Action) -> str:
    if action.type is not None:
        return getattr(action.type, "__name__", str(action.type))
    if isinstance(action, (argparse._StoreTrueAction, argparse._StoreFalseAction)):
        return "flag"
    return "string"


def reference_bytes() -> bytes:
    leaves = _leaf_parsers(build_parser())
    lines = [
        "# PloidyPatch v1.0 command and parameter reference",
        "",
        "This file is generated from `ploidypatch.cli.build_parser` for package",
        f"version `{__version__}`. It documents {len(leaves)} executable leaf commands.",
        "Do not edit it by hand; run `scripts/export_cli_parameters_v1.py`.",
        "",
        "Every output path is non-overwriting unless a command explicitly documents",
        "otherwise. Built-in `--help` remains the authoritative runtime view.",
        "",
        "## Contents",
        "",
    ]
    for path, _ in leaves:
        command = " ".join(path)
        anchor = "ploidypatch-" + "-".join(path)
        lines.append(f"- [`ploidypatch {command}`](#{anchor})")
    lines.append("")
    for path, parser in leaves:
        command = " ".join(path)
        lines.extend(
            [
                f"## `ploidypatch {command}`",
                "",
                _escape(parser.description or parser.prog),
                "",
                "```text",
                parser.format_usage().strip(),
                "```",
                "",
                "| Argument | Required | Type | Default | Choices | Description |",
                "|---|---:|---|---|---|---|",
            ]
        )
        actions = [
            action
            for action in parser._actions
            if not isinstance(action, argparse._SubParsersAction)
        ]
        for action in actions:
            required = "yes" if getattr(action, "required", False) else "no"
            help_text = action.help or "—"
            lines.append(
                "| "
                + " | ".join(
                    (
                        _argument_label(action),
                        required,
                        f"`{_escape(_argument_type(action))}`",
                        f"`{_escape(_display_default(action))}`",
                        _display_choices(action),
                        _escape(help_text),
                    )
                )
                + " |"
            )
        lines.append("")
    return ("\n".join(lines).rstrip() + "\n").encode("utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--check", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    expected = reference_bytes()
    if args.check:
        if not args.output.is_file() or args.output.is_symlink():
            raise FileNotFoundError(f"missing regular parameter reference: {args.output}")
        if args.output.read_bytes() != expected:
            raise ValueError(f"CLI parameter reference differs: {args.output}")
        print(f"CLI parameter reference is current: {args.output}")
        return 0
    if args.output.exists() or args.output.is_symlink():
        raise FileExistsError(f"refusing to overwrite parameter reference: {args.output}")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_bytes(expected)
    print(f"wrote CLI parameter reference: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
