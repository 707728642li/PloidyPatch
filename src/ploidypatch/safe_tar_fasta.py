"""Fail-closed extraction of one immutable FASTA from a tar-gzip container."""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import os
from pathlib import Path, PurePosixPath
import re
import stat
import tarfile
import uuid
from typing import BinaryIO


_SHA256 = re.compile(r"[0-9a-f]{64}")
_FASTA_ALPHABET = frozenset(b"ACGTRYSWKMBDHVN")
_MAX_FASTA_LINE_BYTES = 1_000_000
_PATH_CHANGING_PAX_KEYS = frozenset({"path", "linkpath"})


class TarFastaValidationError(ValueError):
    """The container, member, FASTA, or immutable binding is unsafe."""


@dataclass(frozen=True)
class TarFastaExtractionResult:
    source: Path
    destination: Path
    outer_bytes: int
    outer_sha256: str
    member_name: str
    member_bytes: int
    member_sha256: str
    fasta_records: int
    fasta_bases: int


def _require_digest(value: str, label: str) -> str:
    if not _SHA256.fullmatch(value):
        raise TarFastaValidationError(
            f"{label} must be an exact lowercase SHA-256 digest"
        )
    return value


def _require_positive_size(value: int, label: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise TarFastaValidationError(f"{label} must be a positive integer")
    return value


def _safe_member_name(value: str) -> str:
    if (
        not value
        or "\\" in value
        or "\x00" in value
        or value.startswith("/")
        or value.endswith("/")
        or "//" in value
        or ":" in value
    ):
        raise TarFastaValidationError(f"Unsafe tar member path: {value!r}")
    parts = value.split("/")
    if any(part in {"", ".", ".."} for part in parts):
        raise TarFastaValidationError(f"Unsafe tar member path: {value!r}")
    path = PurePosixPath(*parts)
    if path.is_absolute() or path.as_posix() != value:
        raise TarFastaValidationError(f"Non-canonical tar member path: {value!r}")
    return value


def _absolute_without_resolution(path: Path) -> Path:
    return Path(os.path.abspath(os.fspath(path)))


def _reject_symlink_components(path: Path, *, include_leaf: bool) -> None:
    absolute = _absolute_without_resolution(path)
    current = Path(absolute.anchor)
    components = absolute.parts[1:] if absolute.anchor else absolute.parts
    if not include_leaf:
        components = components[:-1]
    for component in components:
        current /= component
        try:
            mode = current.lstat().st_mode
        except FileNotFoundError:
            continue
        if stat.S_ISLNK(mode):
            raise TarFastaValidationError(
                f"Symlinked path components are forbidden: {current}"
            )


def _sha256_stream(handle: BinaryIO) -> str:
    digest = hashlib.sha256()
    while chunk := handle.read(8 * 1024 * 1024):
        digest.update(chunk)
    return digest.hexdigest()


class _FastaValidator:
    def __init__(self) -> None:
        self.identifiers: set[str] = set()
        self.records = 0
        self.total_bases = 0
        self._current_identifier: str | None = None
        self._current_bases = 0

    def consume(self, raw_line: bytes) -> None:
        if len(raw_line) > _MAX_FASTA_LINE_BYTES:
            raise TarFastaValidationError("FASTA line exceeds the safety limit")
        line = raw_line[:-1] if raw_line.endswith(b"\n") else raw_line
        if line.endswith(b"\r"):
            line = line[:-1]
        if not line:
            raise TarFastaValidationError("Blank FASTA lines are forbidden")
        if b"\r" in line or b"\x00" in line:
            raise TarFastaValidationError("FASTA contains forbidden control bytes")

        if line.startswith(b">"):
            if self._current_identifier is not None and self._current_bases == 0:
                raise TarFastaValidationError(
                    f"FASTA record has no sequence: {self._current_identifier}"
                )
            body = line[1:]
            try:
                decoded = body.decode("ascii")
            except UnicodeDecodeError as error:
                raise TarFastaValidationError("FASTA headers must be ASCII") from error
            if any(ord(character) < 32 and character != "\t" for character in decoded):
                raise TarFastaValidationError("FASTA header contains control characters")
            fields = decoded.split()
            if not fields:
                raise TarFastaValidationError("FASTA header has no identifier")
            identifier = fields[0]
            if identifier in self.identifiers:
                raise TarFastaValidationError(
                    f"Duplicate FASTA identifier: {identifier}"
                )
            self.identifiers.add(identifier)
            self.records += 1
            self._current_identifier = identifier
            self._current_bases = 0
            return

        if self._current_identifier is None:
            raise TarFastaValidationError("FASTA sequence appears before the first header")
        upper = line.upper()
        invalid = set(upper) - _FASTA_ALPHABET
        if invalid:
            rendered = ",".join(f"0x{value:02x}" for value in sorted(invalid))
            raise TarFastaValidationError(
                f"FASTA sequence contains invalid symbols: {rendered}"
            )
        self._current_bases += len(line)
        self.total_bases += len(line)

    def finish(self) -> tuple[int, int]:
        if self.records == 0:
            raise TarFastaValidationError("FASTA has no records")
        if self._current_bases == 0:
            raise TarFastaValidationError(
                f"FASTA record has no sequence: {self._current_identifier}"
            )
        return self.records, self.total_bases


def _validate_member(member: tarfile.TarInfo, expected_name: str) -> None:
    _safe_member_name(member.name)
    if member.name != expected_name:
        raise TarFastaValidationError(
            f"Unexpected tar member name: expected={expected_name!r}, "
            f"observed={member.name!r}"
        )
    if not member.isreg() or member.type not in {tarfile.REGTYPE, tarfile.AREGTYPE}:
        raise TarFastaValidationError("The sole tar member must be a regular file")
    if getattr(member, "sparse", None) is not None:
        raise TarFastaValidationError("Sparse tar members are forbidden")
    pax_keys = set(member.pax_headers)
    if pax_keys & _PATH_CHANGING_PAX_KEYS or any(
        key.startswith("GNU.sparse") for key in pax_keys
    ):
        raise TarFastaValidationError("Path-changing or sparse PAX metadata is forbidden")


def _temporary_path(destination: Path) -> Path:
    return destination.with_name(f".{destination.name}.tmp.{uuid.uuid4().hex}")


def _fsync_directory_best_effort(directory: Path) -> None:
    if os.name == "nt":
        return
    descriptor = os.open(directory, os.O_RDONLY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def extract_single_member_tar_fasta(
    source: str | Path,
    destination: str | Path,
    *,
    expected_outer_bytes: int,
    expected_outer_sha256: str,
    expected_member_name: str,
    expected_member_bytes: int,
    expected_member_sha256: str,
) -> TarFastaExtractionResult:
    """Extract one hash-bound regular FASTA without traversal or overwrite.

    The source is opened once without following its leaf symlink. Both the
    container and member require exact byte counts and SHA-256 digests. The
    destination is published through a same-directory hard link, which is
    atomic and fails if any filesystem object already occupies that name.
    """

    outer_bytes = _require_positive_size(expected_outer_bytes, "outer bytes")
    member_bytes = _require_positive_size(expected_member_bytes, "member bytes")
    outer_sha256 = _require_digest(expected_outer_sha256, "outer SHA-256")
    member_sha256 = _require_digest(expected_member_sha256, "member SHA-256")
    member_name = _safe_member_name(expected_member_name)

    source_path = Path(source)
    destination_path = Path(destination)
    _reject_symlink_components(source_path, include_leaf=True)
    _reject_symlink_components(destination_path, include_leaf=False)
    parent = destination_path.parent
    if not parent.is_dir() or parent.is_symlink():
        raise TarFastaValidationError(
            f"Destination parent must be an existing non-symlink directory: {parent}"
        )
    if os.path.lexists(destination_path):
        raise FileExistsError(f"Refusing to overwrite extraction target: {destination_path}")

    flags = os.O_RDONLY | getattr(os, "O_BINARY", 0) | getattr(os, "O_NOFOLLOW", 0)
    try:
        descriptor = os.open(source_path, flags)
    except OSError as error:
        raise TarFastaValidationError(
            f"Cannot safely open tar source: {source_path}"
        ) from error

    temporary: Path | None = None
    try:
        with os.fdopen(descriptor, "rb", closefd=True) as source_handle:
            opened = os.fstat(source_handle.fileno())
            if not stat.S_ISREG(opened.st_mode):
                raise TarFastaValidationError("Tar source must be a regular file")
            if opened.st_size != outer_bytes:
                raise TarFastaValidationError(
                    f"Outer byte count differs: expected={outer_bytes}, "
                    f"observed={opened.st_size}"
                )
            observed_outer = _sha256_stream(source_handle)
            if observed_outer != outer_sha256:
                raise TarFastaValidationError(
                    f"Outer SHA-256 differs: expected={outer_sha256}, "
                    f"observed={observed_outer}"
                )
            source_handle.seek(0)
            try:
                archive_context = tarfile.open(fileobj=source_handle, mode="r:gz")
            except (OSError, EOFError, tarfile.TarError) as error:
                raise TarFastaValidationError("Invalid tar-gzip container") from error
            with archive_context as archive:
                try:
                    members = archive.getmembers()
                except (OSError, EOFError, tarfile.TarError) as error:
                    raise TarFastaValidationError("Invalid tar-gzip member table") from error
                global_pax_keys = set(archive.pax_headers)
                if global_pax_keys & _PATH_CHANGING_PAX_KEYS or any(
                    key.startswith("GNU.sparse") for key in global_pax_keys
                ):
                    raise TarFastaValidationError(
                        "Path-changing or sparse global PAX metadata is forbidden"
                    )
                if len(members) != 1:
                    raise TarFastaValidationError(
                        f"Tar archive must contain exactly one member; observed={len(members)}"
                    )
                member = members[0]
                _validate_member(member, member_name)
                if member.size != member_bytes:
                    raise TarFastaValidationError(
                        f"Member byte count differs: expected={member_bytes}, "
                        f"observed={member.size}"
                    )
                extracted = archive.extractfile(member)
                if extracted is None:
                    raise TarFastaValidationError("Cannot stream the sole tar member")

                temporary = _temporary_path(destination_path)
                write_flags = (
                    os.O_WRONLY
                    | os.O_CREAT
                    | os.O_EXCL
                    | getattr(os, "O_BINARY", 0)
                    | getattr(os, "O_NOFOLLOW", 0)
                )
                output_descriptor = os.open(temporary, write_flags, 0o600)
                digest = hashlib.sha256()
                observed_member_bytes = 0
                validator = _FastaValidator()
                with extracted, os.fdopen(
                    output_descriptor, "wb", closefd=True
                ) as output_handle:
                    for raw_line in extracted:
                        validator.consume(raw_line)
                        output_handle.write(raw_line)
                        digest.update(raw_line)
                        observed_member_bytes += len(raw_line)
                        if observed_member_bytes > member_bytes:
                            raise TarFastaValidationError(
                                "Extracted member exceeds its frozen byte count"
                            )
                    output_handle.flush()
                    os.fsync(output_handle.fileno())
                records, bases = validator.finish()
                if observed_member_bytes != member_bytes:
                    raise TarFastaValidationError(
                        f"Extracted byte count differs: expected={member_bytes}, "
                        f"observed={observed_member_bytes}"
                    )
                observed_member_sha256 = digest.hexdigest()
                if observed_member_sha256 != member_sha256:
                    raise TarFastaValidationError(
                        f"Member SHA-256 differs: expected={member_sha256}, "
                        f"observed={observed_member_sha256}"
                    )

            source_handle.seek(0)
            repeated_outer = _sha256_stream(source_handle)
            closed = os.fstat(source_handle.fileno())
            if repeated_outer != outer_sha256 or (
                opened.st_dev,
                opened.st_ino,
                opened.st_size,
                opened.st_mtime_ns,
            ) != (
                closed.st_dev,
                closed.st_ino,
                closed.st_size,
                closed.st_mtime_ns,
            ):
                raise TarFastaValidationError("Tar source changed during extraction")

        assert temporary is not None
        # POSIX permits unlinking a read-only file; Windows does not. The
        # immutable guarantee here is hash binding plus atomic non-overwrite,
        # not a platform-dependent permission bit.
        if os.name != "nt":
            os.chmod(temporary, 0o444)
        _reject_symlink_components(destination_path, include_leaf=False)
        try:
            os.link(temporary, destination_path, follow_symlinks=False)
        except FileExistsError as error:
            raise FileExistsError(
                f"Refusing to overwrite extraction target: {destination_path}"
            ) from error
        except OSError as error:
            raise TarFastaValidationError(
                "Atomic non-overwriting publication by hard link failed"
            ) from error
        temporary.unlink()
        temporary = None
        _fsync_directory_best_effort(parent)
        return TarFastaExtractionResult(
            source=source_path,
            destination=destination_path,
            outer_bytes=outer_bytes,
            outer_sha256=outer_sha256,
            member_name=member_name,
            member_bytes=member_bytes,
            member_sha256=member_sha256,
            fasta_records=records,
            fasta_bases=bases,
        )
    except (TarFastaValidationError, FileExistsError):
        raise
    except (OSError, EOFError, tarfile.TarError) as error:
        raise TarFastaValidationError("Safe tar FASTA extraction failed") from error
    finally:
        if temporary is not None and os.path.lexists(temporary):
            temporary.unlink()
