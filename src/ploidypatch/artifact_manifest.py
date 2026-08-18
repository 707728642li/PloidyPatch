"""Strict SHA-256 artifact manifests with exact file-universe validation."""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
from pathlib import Path, PurePosixPath
import re
from typing import Iterable


_DIGEST = re.compile(r"[0-9a-f]{64}")


@dataclass(frozen=True)
class ManifestEntry:
    relative_path: PurePosixPath
    sha256: str


def sha256_file(path: str | Path) -> str:
    target = Path(path)
    digest = hashlib.sha256()
    with target.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def _safe_manifest_path(raw: str, *, line_number: int) -> PurePosixPath:
    if (
        not raw
        or "\\" in raw
        or "\x00" in raw
        or raw.startswith("/")
        or raw.endswith("/")
        or "//" in raw
        or ":" in raw
    ):
        raise ValueError(f"Unsafe manifest path at line {line_number}: {raw!r}")
    parts = raw.split("/")
    if any(part in {"", ".", ".."} for part in parts):
        raise ValueError(f"Unsafe manifest path at line {line_number}: {raw!r}")
    path = PurePosixPath(*parts)
    if path.is_absolute() or path.as_posix() != raw:
        raise ValueError(f"Non-canonical manifest path at line {line_number}: {raw!r}")
    return path


def read_sha256sums(path: str | Path) -> tuple[ManifestEntry, ...]:
    manifest = Path(path)
    if not manifest.is_file() or manifest.is_symlink() or manifest.stat().st_size == 0:
        raise ValueError(f"Missing, empty or symlinked SHA256SUMS manifest: {manifest}")
    entries: list[ManifestEntry] = []
    seen: set[str] = set()
    seen_folded: set[str] = set()
    for line_number, line in enumerate(
        manifest.read_text(encoding="utf-8").splitlines(), start=1
    ):
        digest, separator, raw = line.partition("  ")
        if not separator or not _DIGEST.fullmatch(digest):
            raise ValueError(f"Malformed SHA256SUMS line {line_number}: {manifest}")
        relative = _safe_manifest_path(raw, line_number=line_number)
        name = relative.as_posix()
        folded = name.casefold()
        if name in seen or folded in seen_folded:
            raise ValueError(f"Duplicate SHA256SUMS path at line {line_number}: {name}")
        seen.add(name)
        seen_folded.add(folded)
        entries.append(ManifestEntry(relative_path=relative, sha256=digest))
    if not entries:
        raise ValueError(f"SHA256SUMS manifest has no entries: {manifest}")
    return tuple(entries)


def _reject_symlinks(root: Path) -> None:
    if root.is_symlink():
        raise ValueError(f"Artifact root may not be a symlink: {root}")
    for path in root.rglob("*"):
        if path.is_symlink():
            raise ValueError(f"Symlink is forbidden in an artifact tree: {path}")


def _relative_files(root: Path) -> dict[str, Path]:
    files: dict[str, Path] = {}
    folded: set[str] = set()
    for path in sorted(root.rglob("*")):
        if not path.is_file():
            continue
        relative = path.relative_to(root).as_posix()
        key = relative.casefold()
        if key in folded:
            raise ValueError(f"Case-colliding artifact paths are forbidden: {relative}")
        folded.add(key)
        files[relative] = path
    return files


def verify_sha256sums(
    root: str | Path,
    manifest_path: str | Path | None = None,
    *,
    ignore_checksum_file: bool = False,
) -> dict[str, str]:
    """Verify digests and require the manifest to cover the exact file universe.

    ``ignore_checksum_file`` may be set only when the manifest itself is inside
    ``root``.  In that mode the manifest must not list itself and it is excluded
    from the actual file universe.
    """

    artifact_root = Path(root)
    if not artifact_root.is_dir():
        raise ValueError(f"Artifact root is not a directory: {artifact_root}")
    manifest = Path(manifest_path) if manifest_path is not None else artifact_root / "SHA256SUMS"
    _reject_symlinks(artifact_root)
    entries = read_sha256sums(manifest)
    actual = _relative_files(artifact_root)

    manifest_relative: str | None = None
    try:
        manifest_relative = manifest.resolve().relative_to(artifact_root.resolve()).as_posix()
    except ValueError:
        if ignore_checksum_file:
            raise ValueError(
                "ignore_checksum_file requires a manifest located inside artifact root"
            ) from None
    if ignore_checksum_file:
        assert manifest_relative is not None
        actual.pop(manifest_relative, None)
        if any(entry.relative_path.as_posix() == manifest_relative for entry in entries):
            raise ValueError("An ignored SHA256SUMS file must not list itself")

    expected = {entry.relative_path.as_posix(): entry.sha256 for entry in entries}
    missing = sorted(set(expected) - set(actual))
    extra = sorted(set(actual) - set(expected))
    if missing or extra:
        raise ValueError(
            f"SHA256SUMS file universe differs; missing={missing}, extra={extra}"
        )
    for relative, expected_digest in expected.items():
        path = actual[relative]
        if not path.is_file() or path.is_symlink():
            raise ValueError(f"Manifest entry is not a regular file: {relative}")
        observed = sha256_file(path)
        if observed != expected_digest:
            raise ValueError(
                f"SHA256 mismatch for {relative}: expected={expected_digest}, observed={observed}"
            )
    return expected


def write_sha256sums(
    root: str | Path,
    manifest_path: str | Path | None = None,
    *,
    excluded_relative_paths: Iterable[str] = (),
) -> Path:
    """Write a canonical non-overwriting manifest for regular, non-symlink files."""

    artifact_root = Path(root)
    if not artifact_root.is_dir():
        raise ValueError(f"Artifact root is not a directory: {artifact_root}")
    manifest = Path(manifest_path) if manifest_path is not None else artifact_root / "SHA256SUMS"
    if manifest.exists() or manifest.is_symlink():
        raise FileExistsError(f"Refusing to overwrite SHA256SUMS manifest: {manifest}")
    _reject_symlinks(artifact_root)
    excluded = {
        _safe_manifest_path(raw, line_number=0).as_posix()
        for raw in excluded_relative_paths
    }
    try:
        relative_manifest = manifest.resolve().relative_to(artifact_root.resolve()).as_posix()
    except ValueError:
        relative_manifest = None
    if relative_manifest is not None:
        excluded.add(relative_manifest)
    files = _relative_files(artifact_root)
    selected = [relative for relative in sorted(files) if relative not in excluded]
    if not selected:
        raise ValueError("Refusing to write an empty SHA256SUMS manifest")
    manifest.parent.mkdir(parents=True, exist_ok=True)
    with manifest.open("x", encoding="utf-8", newline="") as handle:
        for relative in selected:
            handle.write(f"{sha256_file(files[relative])}  {relative}\n")
    return manifest
