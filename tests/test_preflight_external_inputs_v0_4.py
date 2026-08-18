from __future__ import annotations

import csv
import gzip
import importlib.util
import json
from pathlib import Path
import sys

import pytest


MODULE_PATH = Path(__file__).parents[1] / "scripts" / "preflight_external_inputs_v0.4.py"
SPEC = importlib.util.spec_from_file_location("external_preflight", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def write_gzip(path: Path, text: str) -> None:
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        handle.write(text)


def test_preflight_reports_metadata_without_pairs_or_labels(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    genome = tmp_path / "genome.fa.gz"
    gff = tmp_path / "genes.gff3.gz"
    protein = tmp_path / "protein.fa.gz"
    write_gzip(genome, ">Chr01\nATGAAATAG\n>scaffold1\nNNNN\n")
    write_gzip(
        gff,
        "Chr01\tx\tgene\t1\t9\t.\t+\t.\tID=g1;Name=gene1\n"
        "Chr01\tx\tmRNA\t1\t9\t.\t+\t.\tID=t1;Parent=g1\n"
        "Chr01\tx\tCDS\t1\t9\t.\t+\t0\tParent=t1;protein_id=p1\n",
    )
    write_gzip(protein, ">p1 gene=g1\nMK\n")
    sources = tmp_path / "sources.tsv"
    with sources.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(MODULE.SOURCE_FIELDS)
        for artifact, path, pattern in (
            ("genome", genome, "^Chr01$"),
            ("gff3", gff, "^Chr01$"),
            ("protein", protein, ""),
        ):
            writer.writerow(("target", "toy", "v1", artifact, pattern, path))
    output = tmp_path / "preflight"
    monkeypatch.setattr(
        sys,
        "argv",
        [str(MODULE_PATH), "--sources", str(sources), "--output-dir", str(output)],
    )
    assert MODULE.main() == 0
    report = json.loads((output / "metadata.json").read_text(encoding="utf-8"))
    assert report["wgd_pairs_enumerated"] is False
    assert report["truth_labels_accessed"] is False
    assert report["species"]["toy"]["genome"]["primary_records"] == 1
    assert report["species"]["toy"]["protein"]["exact_alias_fraction"] == 1.0
    assert (output / "SHA256SUMS").is_file()


def test_preflight_rejects_overwrite(tmp_path: Path) -> None:
    output = tmp_path / "already"
    output.mkdir()
    with pytest.raises(FileExistsError, match="overwrite"):
        monkeypatch = pytest.MonkeyPatch()
        try:
            monkeypatch.setattr(
                sys,
                "argv",
                [
                    str(MODULE_PATH),
                    "--sources",
                    str(tmp_path / "missing.tsv"),
                    "--output-dir",
                    str(output),
                ],
            )
            MODULE.main()
        finally:
            monkeypatch.undo()
