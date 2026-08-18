from __future__ import annotations

import gzip
import hashlib
import importlib.util
from pathlib import Path

import pytest


SCRIPT = Path(__file__).parents[1] / "scripts/normalize_maker_transcript_hierarchy_v0.5.py"
SPEC = importlib.util.spec_from_file_location("maker_hierarchy_v05", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def _write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8", newline="")
    return path


def _cds_digest(text: str) -> str:
    digest = hashlib.sha256()
    for line in text.splitlines():
        if not line.startswith("#") and line.split("\t")[2] == "CDS":
            digest.update((line + "\n").encode())
    return digest.hexdigest()


def test_shared_gene_transcript_ids_are_namespaced_without_cds_change(tmp_path: Path) -> None:
    source_text = (
        "##gff-version 3\n"
        "Chr01\tmaker\tgene\t1\t30\t.\t+\t.\tID=Ae01g1\n"
        "Chr01\tmaker\tmRNA\t1\t30\t.\t+\t.\tID=Ae01g1;Parent=Ae01g1\n"
        "Chr01\tmaker\tCDS\t1\t9\t.\t+\t0\tID=Ae01g1:CDS;Parent=Ae01g1\n"
        "Chr01\tmaker\tCDS\t20\t30\t.\t+\t0\tID=Ae01g1:CDS;Parent=Ae01g1\n"
    )
    source = _write(tmp_path / "source.gff3", source_text)
    output = tmp_path / "normalized.gff3"
    report = MODULE.normalize_hierarchy(
        input_gff=source, output_gff=output, mode="shared_gene_transcript_id"
    )
    observed = output.read_text(encoding="utf-8")
    assert "ID=gene:Ae01g1;ploidypatch_original_gene_id=Ae01g1" in observed
    assert "ID=Ae01g1;Parent=gene:Ae01g1" in observed
    assert observed.count("\tCDS\t") == 2
    assert report["coordinate_or_cds_changes"] is False
    assert report["cds_rows_sha256"] == {
        "input": _cds_digest(source_text),
        "output": _cds_digest(source_text),
    }
    assert Path(str(output) + ".manifest.json").is_file()


def test_direct_gene_children_receive_one_synthetic_transcript(tmp_path: Path) -> None:
    source_text = (
        "##gff-version 3\n"
        "ARU1.0ch01\tmaker\tgene\t10\t80\t.\t-\t.\tID=Aru01g1\n"
        "ARU1.0ch01\tmaker\texon\t60\t80\t.\t-\t.\tID=e1;Parent=Aru01g1\n"
        "ARU1.0ch01\tmaker\tCDS\t60\t80\t.\t-\t0\tID=c1;Parent=Aru01g1\n"
        "ARU1.0ch01\tmaker\tCDS\t10\t30\t.\t-\t0\tID=c2;Parent=Aru01g1\n"
    )
    source = _write(tmp_path / "source.gff3", source_text)
    output = tmp_path / "normalized.gff3"
    report = MODULE.normalize_hierarchy(
        input_gff=source, output_gff=output, mode="gene_with_direct_children"
    )
    observed = output.read_text(encoding="utf-8")
    assert observed.count("\tmRNA\t") == 1
    assert (
        "ARU1.0ch01\tmaker\tmRNA\t10\t80\t.\t-\t.\t"
        "ID=Aru01g1;Parent=gene:Aru01g1;ploidypatch_synthetic_transcript=true\n"
    ) in observed
    assert "ID=c1;Parent=Aru01g1" in observed
    assert report["source_genes"] == report["output_transcripts"] == 1
    assert report["cds_rows_sha256"] == {
        "input": _cds_digest(source_text),
        "output": _cds_digest(source_text),
    }


def test_gzip_provider_input_is_supported_and_cds_rows_remain_exact(tmp_path: Path) -> None:
    source_text = (
        "##gff-version 3\n"
        "Chr29\tmaker\tgene\t4\t42\t.\t+\t.\tID=Ae29g1\n"
        "Chr29\tmaker\tmRNA\t4\t42\t.\t+\t.\tID=Ae29g1;Parent=Ae29g1\n"
        "Chr29\tmaker\tCDS\t4\t12\t7.5\t+\t1\tID=cds1;Parent=Ae29g1\n"
    )
    source = tmp_path / "source.gff3.gz"
    with gzip.open(source, "wt", encoding="utf-8", newline="") as handle:
        handle.write(source_text)
    output = tmp_path / "normalized.gff3"
    report = MODULE.normalize_hierarchy(
        input_gff=source, output_gff=output, mode="shared_gene_transcript_id"
    )
    input_cds = [line for line in source_text.splitlines() if "\tCDS\t" in line]
    output_cds = [
        line for line in output.read_text(encoding="utf-8").splitlines() if "\tCDS\t" in line
    ]
    assert output_cds == input_cds
    assert report["source"]["path"] == "source.gff3.gz"
    assert report["cds_rows_sha256"]["input"] == report["cds_rows_sha256"]["output"]


@pytest.mark.parametrize(
    ("mode", "body", "message"),
    [
        (
            "shared_gene_transcript_id",
            "1\tx\tgene\t1\t9\t.\t+\t.\tID=g1\n"
            "1\tx\tmRNA\t1\t9\t.\t+\t.\tID=t1;Parent=g1\n"
            "1\tx\tCDS\t1\t9\t.\t+\t0\tID=c1;Parent=t1\n",
            "Parent to equal its ID",
        ),
        (
            "gene_with_direct_children",
            "1\tx\tgene\t1\t9\t.\t+\t.\tID=g1\n"
            "1\tx\tmRNA\t1\t9\t.\t+\t.\tID=t1;Parent=g1\n",
            "forbids source mRNA",
        ),
    ],
)
def test_adapter_fails_closed_on_wrong_provider_shape(
    tmp_path: Path, mode: str, body: str, message: str
) -> None:
    with pytest.raises(ValueError, match=message):
        MODULE.normalize_hierarchy(
            input_gff=_write(tmp_path / "source.gff3", body),
            output_gff=tmp_path / "normalized.gff3",
            mode=mode,
        )


def test_adapter_refuses_overwrite(tmp_path: Path) -> None:
    source = _write(
        tmp_path / "source.gff3",
        "1\tx\tgene\t1\t9\t.\t+\t.\tID=g1\n"
        "1\tx\tCDS\t1\t9\t.\t+\t0\tID=c1;Parent=g1\n",
    )
    output = _write(tmp_path / "normalized.gff3", "occupied\n")
    with pytest.raises(FileExistsError, match="Refusing to overwrite"):
        MODULE.normalize_hierarchy(
            input_gff=source,
            output_gff=output,
            mode="gene_with_direct_children",
        )


def test_adapter_removes_partial_output_after_unresolved_parent(tmp_path: Path) -> None:
    source = _write(
        tmp_path / "source.gff3",
        "1\tx\tgene\t1\t9\t.\t+\t.\tID=g1\n"
        "1\tx\tCDS\t1\t9\t.\t+\t0\tID=c1;Parent=unknown\n",
    )
    output = tmp_path / "normalized.gff3"
    with pytest.raises(ValueError, match="unresolved"):
        MODULE.normalize_hierarchy(
            input_gff=source,
            output_gff=output,
            mode="gene_with_direct_children",
        )
    assert not output.exists()
    assert not Path(str(output) + ".working").exists()
    assert not Path(str(output) + ".manifest.json").exists()
