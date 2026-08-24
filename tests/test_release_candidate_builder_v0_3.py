from __future__ import annotations

import gzip
import importlib.util
import io
import json
from pathlib import Path
import tarfile

import pytest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/build_release_candidate_v0.3.py"


def load_module():
    spec = importlib.util.spec_from_file_location("release_candidate_v03", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_sdist(
    path: Path,
    *,
    root_mtime: float,
    file_mtime: float,
    gzip_mtime: int,
    executable: bool = False,
) -> None:
    with path.open("xb") as raw:
        with gzip.GzipFile(
            filename=path.name,
            mode="wb",
            fileobj=raw,
            mtime=gzip_mtime,
        ) as compressed:
            with tarfile.open(fileobj=compressed, mode="w", format=tarfile.PAX_FORMAT) as archive:
                root = tarfile.TarInfo("ploidypatch-0.1.0a0")
                root.type = tarfile.DIRTYPE
                root.mode = 0o775
                root.mtime = root_mtime
                archive.addfile(root)
                member = tarfile.TarInfo("ploidypatch-0.1.0a0/README.md")
                payload = b"reproducible\n"
                member.size = len(payload)
                member.mode = 0o775 if executable else 0o664
                member.mtime = file_mtime
                member.uid = 1000
                member.gid = 1001
                member.uname = "builder"
                member.gname = "group"
                archive.addfile(member, io.BytesIO(payload))


def test_canonical_sdist_is_byte_identical_across_metadata_variation(tmp_path: Path) -> None:
    module = load_module()
    first = tmp_path / "first.tar.gz"
    second = tmp_path / "second.tar.gz"
    write_sdist(first, root_mtime=100.25, file_mtime=200.5, gzip_mtime=300)
    write_sdist(second, root_mtime=400.75, file_mtime=500.125, gzip_mtime=600)
    canonical_first = tmp_path / "canonical-first.tar.gz"
    canonical_second = tmp_path / "canonical-second.tar.gz"
    names_first = module.canonicalize_sdist(first, canonical_first, 1_700_000_000)
    names_second = module.canonicalize_sdist(second, canonical_second, 1_700_000_000)
    assert names_first == names_second
    assert canonical_first.read_bytes() == canonical_second.read_bytes()
    with tarfile.open(canonical_first, "r:gz") as archive:
        members = archive.getmembers()
        assert [member.name for member in members] == sorted(names_first)
        assert all(member.mtime == 1_700_000_000 for member in members)
        assert all(member.uid == member.gid == 0 for member in members)
        assert all(member.uname == member.gname == "" for member in members)
        assert members[0].mode == 0o755
        assert members[1].mode == 0o644
        assert archive.extractfile(members[1]).read() == b"reproducible\n"


def test_canonical_sdist_preserves_executable_semantics_and_refuses_overwrite(
    tmp_path: Path,
) -> None:
    module = load_module()
    source = tmp_path / "source.tar.gz"
    write_sdist(source, root_mtime=1, file_mtime=2, gzip_mtime=3, executable=True)
    destination = tmp_path / "canonical.tar.gz"
    module.canonicalize_sdist(source, destination, 1_700_000_000)
    with tarfile.open(destination, "r:gz") as archive:
        assert archive.getmember("ploidypatch-0.1.0a0/README.md").mode == 0o755
    with pytest.raises(module.ReleaseCandidateError, match="refusing to overwrite"):
        module.canonicalize_sdist(source, destination, 1_700_000_000)


def test_canonical_sdist_rejects_unsupported_members(tmp_path: Path) -> None:
    module = load_module()
    source = tmp_path / "source.tar.gz"
    with tarfile.open(source, "w:gz") as archive:
        member = tarfile.TarInfo("ploidypatch-0.1.0a0/link")
        member.type = tarfile.SYMTYPE
        member.linkname = "../../outside"
        archive.addfile(member)
    destination = tmp_path / "canonical.tar.gz"
    with pytest.raises(module.ReleaseCandidateError, match="unsupported sdist member type"):
        module.canonicalize_sdist(source, destination, 1_700_000_000)
    assert not destination.exists()


def test_v03_policy_is_explicit_and_owner_boundary_is_unchanged() -> None:
    module = load_module()
    assert module.SCHEMA == "ploidypatch.release_candidate_evidence.v3"
    assert module.CANONICAL_SDIST_POLICY["mtime"] == "citation_release_date_midnight_utc"
    assert module.CANONICAL_SDIST_POLICY["uid"] == 0
    assert module.CANONICAL_SDIST_POLICY["gzip_filename"] == ""
    assert module.archive_doi(ROOT) == "10.5281/zenodo.21875560"
    assert module.source_date_epoch(ROOT) == 1_787_529_600
    source = SCRIPT.read_text(encoding="utf-8")
    assert '"formal_public_release": False' not in source
    assert '"missing_owner_fields_unchanged"' not in source
    assert '"artifact_status": "verified_release_candidate"' in source


def test_example_command_log_normalization_removes_only_run_root(tmp_path: Path) -> None:
    module = load_module()
    wheel_root = tmp_path / "wheel_example"
    sdist_root = tmp_path / "sdist_example"
    wheel_root.mkdir()
    sdist_root.mkdir()
    wheel_log = wheel_root / "command_log.jsonl"
    sdist_log = sdist_root / "command_log.jsonl"
    wheel_log.write_text(
        json.dumps(
            {
                "argv": ["patch", "apply", str(wheel_root / "input.gff3")],
                "report": {"output": str(wheel_root / "output.gff3")},
            }
        )
        + "\n",
        encoding="utf-8",
    )
    sdist_log.write_text(
        json.dumps(
            {
                "argv": ["patch", "apply", str(sdist_root / "input.gff3")],
                "report": {"output": str(sdist_root / "output.gff3")},
            }
        )
        + "\n",
        encoding="utf-8",
    )
    assert module.normalized_example_command_log(
        wheel_log, wheel_root
    ) == module.normalized_example_command_log(sdist_log, sdist_root)
