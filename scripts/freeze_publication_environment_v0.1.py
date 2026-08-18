#!/usr/bin/env python3
"""Freeze an isolated plotting environment for publication artifacts."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import stat
import subprocess
import sys
from pathlib import Path


class FreezeError(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise FreezeError(message)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def capture(argv: list[str]) -> str:
    completed = subprocess.run(
        argv,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    require(
        completed.returncode == 0,
        f"command failed ({completed.returncode}): {' '.join(argv)}: {completed.stderr}",
    )
    return completed.stdout


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--environment-prefix", required=True, type=Path)
    parser.add_argument("--conda", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    prefix = args.environment_prefix.resolve(strict=True)
    conda = args.conda.resolve(strict=True)
    python = prefix / "bin/python"
    require(python.is_file(), f"environment Python is missing: {python}")
    output = args.output_dir.resolve()
    working = Path(f"{output}.working")
    require(not output.exists() and not working.exists(), f"refusing to overwrite: {output}")
    working.mkdir(parents=True)
    try:
        (working / "conda-explicit.txt").write_text(
            capture([os.fspath(conda), "list", "--explicit", "-p", os.fspath(prefix)]),
            encoding="utf-8",
        )
        (working / "pip-freeze.txt").write_text(
            capture([os.fspath(python), "-m", "pip", "freeze", "--all"]),
            encoding="utf-8",
        )
        software = json.loads(
            capture(
                [
                    os.fspath(python),
                    "-c",
                    (
                        "import json,matplotlib,numpy,platform,sys;"
                        "print(json.dumps({'python':platform.python_version(),"
                        "'implementation':platform.python_implementation(),"
                        "'matplotlib':matplotlib.__version__,'numpy':numpy.__version__,"
                        "'platform':platform.platform(),'executable':sys.executable},"
                        "sort_keys=True))"
                    ),
                ]
            )
        )
        software.update(
            {
                "schema_version": "ploidypatch.publication_environment.v1",
                "environment_prefix": os.fspath(prefix),
                "conda": {
                    "path": os.fspath(conda),
                    "sha256": sha256_file(conda),
                    "version": capture([os.fspath(conda), "--version"]).strip(),
                },
                "pip": capture(
                    [os.fspath(python), "-m", "pip", "--version"]
                ).strip(),
                "generator": {
                    "path": os.fspath(Path(__file__).resolve()),
                    "sha256": sha256_file(Path(__file__).resolve()),
                },
            }
        )
        (working / "software.json").write_text(
            json.dumps(software, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        files = sorted(path for path in working.iterdir() if path.is_file())
        with (working / "SHA256SUMS").open(
            "x", encoding="utf-8", newline=""
        ) as handle:
            for path in files:
                handle.write(f"{sha256_file(path)}  ./{path.name}\n")
        for path in working.iterdir():
            path.chmod(stat.S_IMODE(path.stat().st_mode) & ~0o222)
        working.chmod(stat.S_IMODE(working.stat().st_mode) & ~0o222)
        working.replace(output)
    except Exception:
        if working.exists():
            shutil.rmtree(working)
        raise
    print(
        json.dumps(
            {
                "environment_freeze": os.fspath(output),
                "sha256sums_sha256": sha256_file(output / "SHA256SUMS"),
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except FreezeError as error:
        print(f"publication environment freeze failed: {error}", file=sys.stderr)
        raise SystemExit(1) from error
