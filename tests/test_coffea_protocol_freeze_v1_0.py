from __future__ import annotations

import hashlib
from pathlib import Path
import tarfile

import pytest

from ploidypatch.artifact_manifest import sha256_file


ROOT = Path(__file__).resolve().parents[1]


def _namespace() -> dict[str, object]:
    value: dict[str, object] = {"__name__": "coffea_freeze_test"}
    exec(
        compile(
            (ROOT / "scripts/freeze_coffea_external_protocol_v1.0.py").read_text(
                encoding="utf-8"
            ),
            "freeze_coffea_external_protocol_v1.0.py",
            "exec",
        ),
        value,
    )
    return value


def _source_freeze(tmp_path: Path, *, drift: bool = False) -> tuple[Path, str]:
    commit = "a" * 40
    root = tmp_path / "freeze"
    source = root / "source"
    (source / "src").mkdir(parents=True)
    (source / "src/example.py").write_text("value = 1\n", encoding="utf-8")
    archive = root / "source.tar"
    with tarfile.open(archive, "w") as handle:
        handle.add(source / "src", arcname="src")
    (root / "source.tar.sha256").write_text(
        f"{sha256_file(archive)}  source.tar\n", encoding="utf-8"
    )
    (root / "source_commit.txt").write_text(commit + "\n", encoding="utf-8")
    if drift:
        (source / "src/example.py").write_text("value = 2\n", encoding="utf-8")
    return root, commit


def test_source_freeze_binds_archive_commit_and_extracted_bytes(tmp_path: Path) -> None:
    root, commit = _source_freeze(tmp_path)
    code_root = _namespace()["verify_source_freeze"](root, commit)
    assert code_root == root / "source"


def test_source_freeze_rejects_extracted_drift(tmp_path: Path) -> None:
    root, commit = _source_freeze(tmp_path, drift=True)
    with pytest.raises(ValueError, match="differs from archive"):
        _namespace()["verify_source_freeze"](root, commit)


def test_coffea_freeze_forbids_results_but_allows_metadata_inputs() -> None:
    namespace = _namespace()
    forbidden = set(namespace["FORBIDDEN_TARGET_ARTIFACTS"])
    assert "data/derived/external_evaluator/coffea_et39_v1.0" in forbidden
    assert "results/evaluator/coffea_et39_v1.0" in forbidden
    assert "data/derived/external_inputs/coffea/v1.0" not in forbidden
    assert all("rank" not in path for path in namespace["IMPLEMENTATION_FILES"])

