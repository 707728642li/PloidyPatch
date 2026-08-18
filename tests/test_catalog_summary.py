from __future__ import annotations

import json
from pathlib import Path

import pytest

from ploidypatch.catalog_summary import (
    MISSING_VALUE,
    summarize_candidate_catalog,
    write_candidate_catalog_summary,
)


CATALOG = (
    "candidate_id\tgene_id\tsubgenome\tsynteny_stratum\n"
    "PPC-1\tgene:G1\tA\texpected_only\n"
    "PPC-2\tgene:G2\tA\texpected_and_cross\n"
    "PPC-3\tgene:G3\tC\texpected_only\n"
    "PPC-4\tgene:G4\t\t\n"
)


def test_catalog_summary_counts_missing_and_joint_strata(tmp_path: Path) -> None:
    catalog = tmp_path / "catalog.tsv"
    catalog.write_text(CATALOG, encoding="utf-8", newline="")

    report = summarize_candidate_catalog(
        catalog,
        columns=["subgenome", "synteny_stratum"],
        crossings=["subgenome,synteny_stratum"],
    )

    assert report["catalog"]["rows"] == 4
    assert report["one_way_counts"]["subgenome"] == {
        MISSING_VALUE: 1,
        "A": 2,
        "C": 1,
    }
    assert {
        (row["subgenome"], row["synteny_stratum"]): row["count"]
        for row in report["joint_counts"]["subgenome__x__synteny_stratum"]
    } == {
        (MISSING_VALUE, MISSING_VALUE): 1,
        ("A", "expected_and_cross"): 1,
        ("A", "expected_only"): 1,
        ("C", "expected_only"): 1,
    }


def test_catalog_summary_writes_json_and_refuses_overwrite(tmp_path: Path) -> None:
    catalog = tmp_path / "catalog.tsv"
    catalog.write_text(CATALOG, encoding="utf-8", newline="")
    output = tmp_path / "summary.json"

    report = write_candidate_catalog_summary(
        catalog, output, columns=["subgenome"]
    )
    assert json.loads(output.read_text(encoding="utf-8")) == report
    with pytest.raises(FileExistsError, match="Refusing to overwrite"):
        write_candidate_catalog_summary(catalog, output, columns=["subgenome"])


def test_catalog_summary_validates_columns_and_crossings(tmp_path: Path) -> None:
    catalog = tmp_path / "catalog.tsv"
    catalog.write_text(CATALOG, encoding="utf-8", newline="")

    with pytest.raises(ValueError, match="missing required"):
        summarize_candidate_catalog(catalog, columns=["unknown"])
    with pytest.raises(ValueError, match="expected COLUMN_A,COLUMN_B"):
        summarize_candidate_catalog(catalog, columns=["subgenome"], crossings=["x"])
