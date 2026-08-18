from __future__ import annotations

import importlib.util
from pathlib import Path
import tarfile
from zipfile import ZipFile

import pytest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/build_release_candidate_v0.2.py"


def load_module():
    spec = importlib.util.spec_from_file_location("release_candidate_v02", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_archive_path_validation_rejects_unsafe_names() -> None:
    module = load_module()
    for value in ("", "/absolute", "../escape", "a/../../escape", r"a\b"):
        with pytest.raises(module.ReleaseCandidateError, match="unsafe archive path"):
            module.safe_archive_name(value)
    assert module.safe_archive_name("ploidypatch-0.1.0a0/src/ploidypatch/cli.py").as_posix().endswith(
        "src/ploidypatch/cli.py"
    )


def test_wheel_audit_rejects_a_mandatory_runtime_dependency(tmp_path: Path) -> None:
    module = load_module()
    wheel = tmp_path / module.WHEEL_NAME
    root = f"{module.PACKAGE}-{module.VERSION}.dist-info"
    with ZipFile(wheel, "w") as archive:
        archive.writestr(
            f"{root}/METADATA",
            f"Metadata-Version: 2.4\r\nName: ploidypatch\r\nVersion: {module.VERSION}\r\n"
            "Requires-Python: >=3.11\r\nRequires-Dist: numpy>=1.26\r\n",
        )
        archive.writestr(
            f"{root}/entry_points.txt",
            "[console_scripts]\nploidypatch = ploidypatch.cli:main\n",
        )
    with pytest.raises(module.ReleaseCandidateError, match="mandatory runtime dependencies"):
        module.audit_wheel(wheel)


def test_wheel_audit_accepts_optional_extra_dependencies(tmp_path: Path) -> None:
    module = load_module()
    wheel = tmp_path / module.WHEEL_NAME
    root = f"{module.PACKAGE}-{module.VERSION}.dist-info"
    with ZipFile(wheel, "w") as archive:
        archive.writestr(
            f"{root}/METADATA",
            f"Metadata-Version: 2.4\r\nName: ploidypatch\r\nVersion: {module.VERSION}\r\n"
            "Requires-Python: >=3.11\r\n"
            'Requires-Dist: numpy>=1.26; extra == "evaluation"\r\n'
            'Requires-Dist: matplotlib>=3.8; extra == "visualization"\r\n',
        )
        archive.writestr(
            f"{root}/entry_points.txt",
            "[console_scripts]\nploidypatch = ploidypatch.cli:main\n",
        )
    record = module.audit_wheel(wheel)
    assert record["no_mandatory_runtime_dependencies"] is True


def test_sdist_audit_rejects_links(tmp_path: Path) -> None:
    module = load_module()
    sdist = tmp_path / module.SDIST_NAME
    with tarfile.open(sdist, "w:gz") as archive:
        member = tarfile.TarInfo(f"{module.PACKAGE}-{module.VERSION}/unsafe")
        member.type = tarfile.SYMTYPE
        member.linkname = "../../outside"
        archive.addfile(member)
    with pytest.raises(module.ReleaseCandidateError, match="unsupported sdist member type"):
        module.audit_sdist(sdist)


def test_release_builder_does_not_hardcode_obsolete_owner_blockers() -> None:
    source = SCRIPT.read_text(encoding="utf-8")
    assert '"distribution_metadata_complete": True' in source
    assert '"archive_doi_required_for_github_release": False' in source
    assert '"formal_public_release": False' not in source
    assert '"missing_owner_fields_unchanged"' not in source
