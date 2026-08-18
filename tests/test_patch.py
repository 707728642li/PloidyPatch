from __future__ import annotations

import json
from pathlib import Path

import pytest

from ploidypatch.patch import (
    apply_annotation_patch,
    create_annotation_patch,
    revert_annotation_patch,
)


SOURCE = (
    "##gff-version 3\n"
    "chr1\ttest\tgene\t1\t30\t.\t+\t.\tID=G1\n"
    "chr1\ttest\tmRNA\t1\t30\t.\t+\t.\tID=T1;Parent=G1\n"
    "chr1\ttest\tCDS\t1\t30\t.\t+\t0\tParent=T1\n"
)


def _write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8", newline="")
    return path


def test_annotation_patch_apply_and_revert_are_exact(tmp_path: Path) -> None:
    source = _write(tmp_path / "source.gff3", SOURCE)
    new_gene = "chr1\tpatch\tgene\t1\t30\t.\t+\t.\tID=G1.fixed\n"
    new_transcript = (
        "chr1\tpatch\tmRNA\t1\t30\t.\t+\t.\tID=T1.fixed;Parent=G1.fixed\n"
    )
    edits = tmp_path / "edits.json"
    edits.write_text(
        json.dumps(
            {
                "event_ids": ["EV1"],
                "operations": [
                    {
                        "source_line_number": 2,
                        "replacement_lines": [new_gene, new_transcript],
                    },
                    {"source_line_number": 3, "replacement_lines": []},
                ],
            }
        ),
        encoding="utf-8",
    )
    patch_path = tmp_path / "patch.json"
    patch = create_annotation_patch(source, edits, patch_path)
    assert patch["schema_version"] == "ploidypatch.annotation_patch.v1"
    assert patch["event_ids"] == ["EV1"]
    assert len(patch["operations"]) == 2

    patched = tmp_path / "patched.gff3"
    apply_report = apply_annotation_patch(source, patch_path, patched)
    assert apply_report["operations"] == 2
    assert patched.read_text(encoding="utf-8") == (
        "##gff-version 3\n" + new_gene + new_transcript
        + "chr1\ttest\tCDS\t1\t30\t.\t+\t0\tParent=T1\n"
    )

    reverted = tmp_path / "reverted.gff3"
    revert_report = revert_annotation_patch(patched, patch_path, reverted)
    assert revert_report["text_sha256"] == patch["source"]["text_sha256"]
    assert reverted.read_bytes() == source.read_bytes()


def test_patch_accepts_hidden_truth_line_edits(tmp_path: Path) -> None:
    source = _write(tmp_path / "source.gff3", SOURCE)
    replacement = "chr1\ttest\tgene\t2\t30\t.\t+\t.\tID=G1\n"
    truth_like = tmp_path / "truth.json"
    truth_like.write_text(
        json.dumps(
            {
                "events": [
                    {
                        "event_id": "EV2",
                        "line_edits": [
                            {
                                "source_line_number": 2,
                                "source_raw_line": SOURCE.splitlines(keepends=True)[1],
                                "perturbed_lines": [{"raw_line": replacement}],
                            }
                        ],
                    }
                ]
            }
        ),
        encoding="utf-8",
    )
    patch_path = tmp_path / "patch.json"
    patch = create_annotation_patch(source, truth_like, patch_path)
    assert patch["event_ids"] == ["EV2"]
    patched = tmp_path / "patched.gff3"
    apply_annotation_patch(source, patch_path, patched)
    assert replacement in patched.read_text(encoding="utf-8")


def test_patch_rejects_tampering_and_overwrite(tmp_path: Path) -> None:
    source = _write(tmp_path / "source.gff3", SOURCE)
    edits = tmp_path / "edits.json"
    edits.write_text(
        json.dumps(
            {
                "operations": [
                    {"source_line_number": 2, "replacement_lines": []}
                ]
            }
        ),
        encoding="utf-8",
    )
    patch_path = tmp_path / "patch.json"
    create_annotation_patch(source, edits, patch_path)
    with pytest.raises(FileExistsError, match="Refusing to overwrite"):
        create_annotation_patch(source, edits, patch_path)

    tampered = _write(tmp_path / "tampered.gff3", SOURCE + "# changed\n")
    with pytest.raises(ValueError, match="line count|checksum"):
        apply_annotation_patch(tampered, patch_path, tmp_path / "bad.gff3")

    patched = tmp_path / "patched.gff3"
    apply_annotation_patch(source, patch_path, patched)
    patched.write_text(
        patched.read_text(encoding="utf-8") + "# changed\n",
        encoding="utf-8",
        newline="",
    )
    with pytest.raises(ValueError, match="line count|checksum"):
        revert_annotation_patch(patched, patch_path, tmp_path / "bad_revert.gff3")
