from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest


SCRIPT = Path(__file__).parents[1] / "scripts" / "build_wgdi_source_alias_gff.py"
SPEC = importlib.util.spec_from_file_location("build_wgdi_source_alias_gff", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def test_builds_exact_name_alias_without_coordinate_change(tmp_path: Path) -> None:
    source = tmp_path / "source.gff3"
    source.write_text(
        "##gff-version 3\n"
        "Chr01\tx\tgene\t10\t20\t.\t+\t.\tID=gene:MD01;Name=MD01\n"
        "Chr01\tx\tmRNA\t10\t20\t.\t+\t.\tID=mRNA:MD01;Parent=gene:MD01\n"
        "Chr01\tx\tgene\t30\t40\t.\t-\t.\tID=gene:MD02;Name=MD02\n"
        "Chr02\tx\tgene\t50\t60\t.\t+\t.\tID=PPCONS_gene_x\n",
        encoding="utf-8",
    )
    representatives = tmp_path / "representatives.tsv"
    representatives.write_text(
        "gene_id\tprotein_id\nMD02\tp2\nMD01\tp1\nPPCONS_gene_x\tpx\n",
        encoding="utf-8",
    )
    output = tmp_path / "aliases.gff3"
    manifest = MODULE.build_alias_gff(
        source_gff=source,
        representatives_tsv=representatives,
        output_gff=output,
    )
    lines = output.read_text(encoding="utf-8").splitlines()
    assert lines[1].split("\t")[3:5] == ["30", "40"]
    assert lines[2].split("\t")[3:5] == ["10", "20"]
    assert lines[3].split("\t")[3:5] == ["50", "60"]
    assert "ID=gene:MD02;Name=MD02;gene_id=MD02" in lines[1]
    assert manifest["truth_access"] is False
    assert manifest["coordinate_change"] is False
    assert "ID=PPCONS_gene_x;gene_id=PPCONS_gene_x" in lines[3]
    assert manifest["counts"]["mapped_gene_ids"] == 3


def test_rejects_missing_or_nonunique_name_alias(tmp_path: Path) -> None:
    source = tmp_path / "source.gff3"
    source.write_text(
        "Chr01\tx\tgene\t1\t2\t.\t+\t.\tID=g1;Name=A\n"
        "Chr02\tx\tgene\t3\t4\t.\t+\t.\tID=g2;Name=A\n",
        encoding="utf-8",
    )
    representatives = tmp_path / "representatives.tsv"
    representatives.write_text("gene_id\nA\n", encoding="utf-8")
    with pytest.raises(ValueError, match="multiple genes"):
        MODULE.build_alias_gff(
            source_gff=source,
            representatives_tsv=representatives,
            output_gff=tmp_path / "out.gff3",
        )
