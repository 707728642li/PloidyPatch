from __future__ import annotations

from io import BytesIO
import hashlib
import os
from pathlib import Path
import tarfile

import pytest

from ploidypatch.safe_tar_fasta import (
    TarFastaValidationError,
    extract_single_member_tar_fasta,
)


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def write_tar(
    path: Path,
    members: list[tuple[tarfile.TarInfo, bytes | None]],
    *,
    tar_format: int = tarfile.GNU_FORMAT,
) -> None:
    with tarfile.open(path, "w:gz", format=tar_format) as archive:
        for member, payload in members:
            if payload is not None:
                member.size = len(payload)
                archive.addfile(member, BytesIO(payload))
            else:
                archive.addfile(member)


def regular_member(name: str, content: bytes) -> tuple[tarfile.TarInfo, bytes]:
    member = tarfile.TarInfo(name)
    member.mode = 0o644
    return member, content


def extract(
    source: Path,
    destination: Path,
    content: bytes,
    *,
    member_name: str = "genome.fa",
    outer_bytes: int | None = None,
    outer_sha256: str | None = None,
    member_bytes: int | None = None,
    member_sha256: str | None = None,
):
    return extract_single_member_tar_fasta(
        source,
        destination,
        expected_outer_bytes=(
            source.stat().st_size if outer_bytes is None else outer_bytes
        ),
        expected_outer_sha256=(
            sha256_bytes(source.read_bytes())
            if outer_sha256 is None
            else outer_sha256
        ),
        expected_member_name=member_name,
        expected_member_bytes=len(content) if member_bytes is None else member_bytes,
        expected_member_sha256=(
            sha256_bytes(content) if member_sha256 is None else member_sha256
        ),
    )


def make_valid_tar(tmp_path: Path) -> tuple[Path, bytes]:
    content = b">chr1 description\nACGTN\n>chr2\nryswkmbdhv\n"
    source = tmp_path / "source.tar.gz"
    write_tar(source, [regular_member("genome.fa", content)])
    return source, content


def temporary_outputs(tmp_path: Path, destination: Path) -> list[Path]:
    return list(tmp_path.glob(f".{destination.name}.tmp.*"))


def test_extracts_one_hash_bound_regular_fasta_atomically(tmp_path: Path) -> None:
    source, content = make_valid_tar(tmp_path)
    destination = tmp_path / "genome.fa"
    result = extract(source, destination, content)

    assert destination.read_bytes() == content
    assert result.outer_bytes == source.stat().st_size
    assert result.outer_sha256 == sha256_bytes(source.read_bytes())
    assert result.member_name == "genome.fa"
    assert result.member_bytes == len(content)
    assert result.member_sha256 == sha256_bytes(content)
    assert result.fasta_records == 2
    assert result.fasta_bases == 15
    assert temporary_outputs(tmp_path, destination) == []


@pytest.mark.parametrize(
    ("keyword", "value", "message"),
    (
        ("expected_outer_bytes", 0, "positive integer"),
        ("expected_member_bytes", -1, "positive integer"),
        ("expected_outer_sha256", "A" * 64, "lowercase SHA-256"),
        ("expected_member_sha256", "not-a-digest", "lowercase SHA-256"),
        ("expected_member_name", "../genome.fa", "Unsafe tar member"),
        ("expected_member_name", "/genome.fa", "Unsafe tar member"),
    ),
)
def test_rejects_invalid_frozen_parameters(
    tmp_path: Path, keyword: str, value: object, message: str
) -> None:
    source, content = make_valid_tar(tmp_path)
    arguments: dict[str, object] = {
        "expected_outer_bytes": source.stat().st_size,
        "expected_outer_sha256": sha256_bytes(source.read_bytes()),
        "expected_member_name": "genome.fa",
        "expected_member_bytes": len(content),
        "expected_member_sha256": sha256_bytes(content),
    }
    arguments[keyword] = value
    with pytest.raises(TarFastaValidationError, match=message):
        extract_single_member_tar_fasta(
            source, tmp_path / "out.fa", **arguments  # type: ignore[arg-type]
        )


def test_rejects_empty_archive(tmp_path: Path) -> None:
    source = tmp_path / "empty.tar.gz"
    write_tar(source, [])
    with pytest.raises(TarFastaValidationError, match="exactly one member"):
        extract(source, tmp_path / "out.fa", b">x\nA\n")


def test_rejects_multiple_members(tmp_path: Path) -> None:
    content = b">x\nA\n"
    source = tmp_path / "multiple.tar.gz"
    write_tar(
        source,
        [
            regular_member("genome.fa", content),
            regular_member("second.fa", content),
        ],
    )
    with pytest.raises(TarFastaValidationError, match="exactly one member"):
        extract(source, tmp_path / "out.fa", content)


@pytest.mark.parametrize("member_name", ("../escape.fa", "/escape.fa", "a\\b.fa"))
def test_rejects_unsafe_archive_member_paths(
    tmp_path: Path, member_name: str
) -> None:
    content = b">x\nA\n"
    source = tmp_path / "unsafe.tar.gz"
    write_tar(source, [regular_member(member_name, content)])
    with pytest.raises(TarFastaValidationError, match="Unsafe tar member"):
        extract(source, tmp_path / "out.fa", content)
    assert not (tmp_path / "escape.fa").exists()


@pytest.mark.parametrize(
    "member_type",
    (
        tarfile.DIRTYPE,
        tarfile.SYMTYPE,
        tarfile.LNKTYPE,
        tarfile.CHRTYPE,
        tarfile.FIFOTYPE,
        tarfile.GNUTYPE_SPARSE,
    ),
)
def test_rejects_non_regular_and_sparse_members(
    tmp_path: Path, member_type: bytes
) -> None:
    source = tmp_path / "nonregular.tar.gz"
    member = tarfile.TarInfo("genome.fa")
    member.type = member_type
    if member_type in {tarfile.SYMTYPE, tarfile.LNKTYPE}:
        member.linkname = "elsewhere.fa"
    write_tar(source, [(member, None)])
    with pytest.raises(TarFastaValidationError, match="regular file|Sparse"):
        extract(source, tmp_path / "out.fa", b">x\nA\n")


def test_rejects_path_changing_pax_metadata(tmp_path: Path) -> None:
    content = b">x\nA\n"
    source = tmp_path / "pax.tar.gz"
    member = tarfile.TarInfo("genome.fa")
    member.pax_headers = {"path": "genome.fa"}
    write_tar(source, [(member, content)], tar_format=tarfile.PAX_FORMAT)
    with pytest.raises(TarFastaValidationError, match="PAX"):
        extract(source, tmp_path / "out.fa", content)


def test_rejects_unexpected_member_name(tmp_path: Path) -> None:
    content = b">x\nA\n"
    source = tmp_path / "wrong-name.tar.gz"
    write_tar(source, [regular_member("other.fa", content)])
    with pytest.raises(TarFastaValidationError, match="Unexpected tar member"):
        extract(source, tmp_path / "out.fa", content)


@pytest.mark.parametrize("field", ("outer_bytes", "outer_sha", "member_bytes", "member_sha"))
def test_rejects_byte_or_digest_mismatch(tmp_path: Path, field: str) -> None:
    source, content = make_valid_tar(tmp_path)
    kwargs: dict[str, object] = {}
    if field == "outer_bytes":
        kwargs["outer_bytes"] = source.stat().st_size + 1
    elif field == "outer_sha":
        kwargs["outer_sha256"] = "0" * 64
    elif field == "member_bytes":
        kwargs["member_bytes"] = len(content) + 1
    else:
        kwargs["member_sha256"] = "0" * 64
    with pytest.raises(TarFastaValidationError, match="byte count|SHA-256"):
        extract(source, tmp_path / "out.fa", content, **kwargs)
    assert not (tmp_path / "out.fa").exists()
    assert temporary_outputs(tmp_path, tmp_path / "out.fa") == []


def test_rejects_truncated_tar_gzip(tmp_path: Path) -> None:
    source, _ = make_valid_tar(tmp_path)
    truncated = tmp_path / "truncated.tar.gz"
    value = source.read_bytes()
    truncated.write_bytes(value[: len(value) // 2])
    with pytest.raises(TarFastaValidationError, match="Invalid|failed"):
        extract(truncated, tmp_path / "out.fa", b">x\nA\n")


@pytest.mark.parametrize(
    ("content", "message"),
    (
        (b"ACGT\n", "before the first header"),
        (b">\nACGT\n", "no identifier"),
        (b">a\n>b\nACGT\n", "has no sequence"),
        (b">a\nACGT\n>a\nACGT\n", "Duplicate FASTA"),
        (b">a\nACGT!\n", "invalid symbols"),
        (b">a\n\nACGT\n", "Blank FASTA"),
        (b">a\xff\nACGT\n", "headers must be ASCII"),
    ),
)
def test_rejects_malformed_fasta_and_cleans_temporary_output(
    tmp_path: Path, content: bytes, message: str
) -> None:
    source = tmp_path / "malformed.tar.gz"
    write_tar(source, [regular_member("genome.fa", content)])
    destination = tmp_path / "out.fa"
    with pytest.raises(TarFastaValidationError, match=message):
        extract(source, destination, content)
    assert not destination.exists()
    assert temporary_outputs(tmp_path, destination) == []


def test_refuses_existing_destination(tmp_path: Path) -> None:
    source, content = make_valid_tar(tmp_path)
    destination = tmp_path / "out.fa"
    destination.write_bytes(b"keep\n")
    with pytest.raises(FileExistsError, match="overwrite"):
        extract(source, destination, content)
    assert destination.read_bytes() == b"keep\n"


def test_refuses_symlinked_source(tmp_path: Path) -> None:
    source, content = make_valid_tar(tmp_path)
    link = tmp_path / "source-link.tar.gz"
    try:
        link.symlink_to(source.name)
    except OSError:
        pytest.skip("Host does not permit symlink creation")
    with pytest.raises(TarFastaValidationError, match="Symlink"):
        extract(link, tmp_path / "out.fa", content)


def test_refuses_symlinked_destination_parent(tmp_path: Path) -> None:
    source, content = make_valid_tar(tmp_path)
    real = tmp_path / "real"
    real.mkdir()
    link = tmp_path / "linked"
    try:
        link.symlink_to(real, target_is_directory=True)
    except OSError:
        pytest.skip("Host does not permit directory symlink creation")
    with pytest.raises(TarFastaValidationError, match="Symlink"):
        extract(source, link / "out.fa", content)
    assert not (real / "out.fa").exists()


def test_refuses_symlink_at_destination(tmp_path: Path) -> None:
    source, content = make_valid_tar(tmp_path)
    real = tmp_path / "real.fa"
    real.write_bytes(b"keep\n")
    destination = tmp_path / "out.fa"
    try:
        destination.symlink_to(real.name)
    except OSError:
        pytest.skip("Host does not permit symlink creation")
    with pytest.raises(FileExistsError, match="overwrite"):
        extract(source, destination, content)
    assert real.read_bytes() == b"keep\n"


def test_publication_race_does_not_overwrite(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    source, content = make_valid_tar(tmp_path)
    destination = tmp_path / "out.fa"
    real_link = os.link

    def competing_link(
        source_path: str | os.PathLike[str],
        destination_path: str | os.PathLike[str],
        *,
        follow_symlinks: bool = True,
    ) -> None:
        Path(destination_path).write_bytes(b"competitor\n")
        real_link(
            source_path,
            destination_path,
            follow_symlinks=follow_symlinks,
        )

    monkeypatch.setattr(os, "link", competing_link)
    with pytest.raises(FileExistsError, match="overwrite"):
        extract(source, destination, content)
    assert destination.read_bytes() == b"competitor\n"
    assert temporary_outputs(tmp_path, destination) == []
