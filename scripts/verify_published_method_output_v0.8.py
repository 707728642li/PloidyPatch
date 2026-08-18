#!/usr/bin/env python3
"""Verify a published GeMoMa/LiftOn output manifest after atomic rename."""
from __future__ import annotations

from pathlib import Path
import sys

from ploidypatch.published_output import verify_published_method_output as verify


def main(argv: list[str] | None = None) -> int:
    values = sys.argv[1:] if argv is None else argv
    if len(values) < 2:
        raise SystemExit("usage: verify_published_method_output_v0.8.py ROOT ROLE...")
    verify(Path(values[0]), set(values[1:]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
