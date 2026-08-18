from __future__ import annotations

import csv
from pathlib import Path, PurePosixPath
import re


ROOT = Path(__file__).parents[1]
MANIFEST = ROOT / "config" / "walnut_external_input_sources_v0.8.tsv"
RATIONALE = ROOT / "docs" / "WALNUT_EXTERNAL_SELECTION_AND_CONTAMINATION_RATIONALE_v0.8.md"


EXPECTED_COLUMNS = (
    "role",
    "species_id",
    "release",
    "artifact",
    "primary_seqid_regex",
    "primary_seqid_table",
    "source_path",
    "bytes",
    "sha256",
    "source_url",
    "container_format",
    "member_name",
    "member_bytes",
    "member_sha256",
)


EXPECTED_SPECIES = {
    "target": {"Juglans_regia"},
    "candidate_reference": {"Juglans_mandshurica", "Carya_illinoinensis"},
    "evaluator_reference": {"Corylus_avellana", "Castanea_mollissima"},
}


EXPECTED_PRIMARY_COUNTS = {
    "Juglans_regia": 16,
    "Juglans_mandshurica": 16,
    "Carya_illinoinensis": 16,
    "Corylus_avellana": 11,
    "Castanea_mollissima": 12,
}


def read_manifest() -> tuple[tuple[str, ...], list[dict[str, str]]]:
    with MANIFEST.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return tuple(reader.fieldnames or ()), list(reader)


def read_primary(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        assert tuple(reader.fieldnames or ()) == ("seqid", "chromosome_label")
        return list(reader)


def test_walnut_manifest_has_exact_disjoint_one_two_two_roles() -> None:
    columns, rows = read_manifest()
    assert columns == EXPECTED_COLUMNS
    assert len(rows) == 15
    assert {row["role"] for row in rows} == set(EXPECTED_SPECIES)
    for role, expected in EXPECTED_SPECIES.items():
        assert {row["species_id"] for row in rows if row["role"] == role} == expected
    role_sets = list(EXPECTED_SPECIES.values())
    assert all(
        role_sets[left].isdisjoint(role_sets[right])
        for left in range(len(role_sets))
        for right in range(left + 1, len(role_sets))
    )
    for species in set().union(*role_sets):
        species_rows = [row for row in rows if row["species_id"] == species]
        assert {row["artifact"] for row in species_rows} == {
            "genome",
            "gff3",
            "protein",
        }
        assert len({row["release"] for row in species_rows}) == 1
        assert len({row["primary_seqid_table"] for row in species_rows}) == 1


def test_walnut_manifest_binds_bytes_digests_urls_and_tar_member() -> None:
    _, rows = read_manifest()
    for row in rows:
        assert int(row["bytes"]) > 0
        assert re.fullmatch(r"[0-9a-f]{64}", row["sha256"])
        assert PurePosixPath(row["source_path"]).is_absolute()
        if row["species_id"] in {"Carya_illinoinensis", "Castanea_mollissima"}:
            assert row["source_url"].startswith("https://ftp.ncbi.nlm.nih.gov/")
        else:
            assert row["source_url"] == ""

    tar_rows = [row for row in rows if row["container_format"] == "tar.gz"]
    assert len(tar_rows) == 1
    tar_row = tar_rows[0]
    assert (tar_row["species_id"], tar_row["artifact"]) == (
        "Juglans_mandshurica",
        "genome",
    )
    assert tar_row["member_name"] == "Juglans_mandshurica.genome.fa"
    assert tar_row["member_bytes"] == "554212158"
    assert tar_row["member_sha256"] == (
        "7b1be8bfb34096526ac94bc6f1806f0257c7b31293b81590225d283f8949679f"
    )
    direct = [row for row in rows if row["container_format"] == "direct"]
    assert len(direct) == 14
    assert all(
        not row["member_name"]
        and not row["member_bytes"]
        and not row["member_sha256"]
        for row in direct
    )


def test_primary_tables_are_exact_unique_and_match_frozen_regexes() -> None:
    _, manifest_rows = read_manifest()
    for species, expected_count in EXPECTED_PRIMARY_COUNTS.items():
        species_rows = [row for row in manifest_rows if row["species_id"] == species]
        table = ROOT / species_rows[0]["primary_seqid_table"]
        rows = read_primary(table)
        assert len(rows) == expected_count
        assert len({row["seqid"] for row in rows}) == expected_count
        assert len({row["chromosome_label"] for row in rows}) == expected_count
        pattern = re.compile(species_rows[0]["primary_seqid_regex"])
        assert all(pattern.fullmatch(row["seqid"]) for row in rows)
        assert all(
            row["primary_seqid_regex"]
            in {"", species_rows[0]["primary_seqid_regex"]}
            for row in species_rows
        )


def test_rationale_keeps_metadata_boundary_and_rejected_roles_explicit() -> None:
    text = RATIONALE.read_text(encoding="utf-8")
    for required in (
        "does not",
        "No WGD pair yield",
        "0.9528485844348267",
        "0.9997615493865762",
        "0.9965136505080359",
        "J. sigillata",
        "Quercus lobata",
        "Cucurbita",
        "reference-anchored topology/backbone construction",
        "final holdout contract",
    ):
        assert required in text
