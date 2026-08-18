#!/usr/bin/env python3
"""Run every test shipped in a clean public PloidyPatch checkout."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def command(extra: list[str]) -> list[str]:
    return [sys.executable, "-m", "pytest", "-q", *extra]


def main() -> int:
    completed = subprocess.run(command(sys.argv[1:]), cwd=ROOT, check=False)
    return completed.returncode


if __name__ == "__main__":
    raise SystemExit(main())
