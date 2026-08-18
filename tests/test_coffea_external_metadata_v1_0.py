from __future__ import annotations

from dataclasses import asdict
import csv
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

from ploidypatch.holdout_contract import (
    FIXED_KNOWN_SUBGENOME_CORE_H1_SCIENTIFIC_PARAMETERS,
    KNOWN_SUBGENOME_CORE_H1_MODEL_VERSION,
    load_holdout_contract,
)


ROOT = Path(__file__).parents[1]
POLICY = ROOT / "config/coffea_external_validation_policy_v1.0.tsv"
EVENT = ROOT / "config/coffea_external_event_definition_v1.0.tsv"
HOMOEOLOGS = ROOT / "config/coffea_et39_homoeolog_groups_v1.0.tsv"
GARDENIA_ALIASES = ROOT / "config/coffea_gardenia_seqid_aliases_v1.0.tsv"
SOURCE_MANIFEST = ROOT / "config/coffea_external_input_sources_v1.0.tsv"
CONTRACT = ROOT / "config/holdouts/coffea_et39_v1.0/contract.json"
PRIMARY_TABLES = {
    "target": ("coffea_arabica_et39_hifi.tsv", 22),
    "eugenioides": ("coffea_eugenioides_bua_v1.tsv", 11),
    "mauritiana": ("coffea_mauritiana_v1.tsv", 11),
    "gardenia": ("gardenia_jasminoides_asm1310374v1.tsv", 11),
    "ophiorrhiza": ("ophiorrhiza_pumila_v1.tsv", 11),
}


def load_script(relative: str, name: str):
    path = ROOT / relative
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


PREFLIGHT = load_script(
    "scripts/preflight_coffea_external_inputs_v1.0.py", "coffea_preflight_v10"
)


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_field_table(path: Path) -> dict[str, str]:
    rows = read_rows(path)
    assert tuple(rows[0]) == ("field", "value")
    values = {row["field"]: row["value"] for row in rows}
    assert len(values) == len(rows)
    assert all(values)
    assert all(values.values())
    return values


def test_coffea_primary_tables_are_exact_and_role_separated() -> None:
    observed_sets: list[set[str]] = []
    for filename, expected_count in PRIMARY_TABLES.values():
        rows = read_rows(ROOT / "config/primary_seqids" / filename)
        assert tuple(rows[0]) == ("seqid", "chromosome_label")
        assert len(rows) == expected_count
        seqids = {row["seqid"] for row in rows}
        labels = {row["chromosome_label"] for row in rows}
        assert len(seqids) == len(rows)
        assert len(labels) == len(rows)
        observed_sets.append(seqids)
    assert sum(len(values) for values in observed_sets) == len(
        set().union(*observed_sets)
    )


def test_et39_homoeolog_groups_partition_all_target_chromosomes() -> None:
    target_rows = read_rows(
        ROOT / "config/primary_seqids/coffea_arabica_et39_hifi.tsv"
    )
    target = {row["seqid"]: row["chromosome_label"] for row in target_rows}
    groups = read_rows(HOMOEOLOGS)
    assert tuple(groups[0]) == (
        "homoeolog_group",
        "c_subgenome_seqid",
        "e_subgenome_seqid",
        "c_label",
        "e_label",
    )
    assert [row["homoeolog_group"] for row in groups] == [
        str(index) for index in range(1, 12)
    ]
    members: list[str] = []
    for row in groups:
        group = row["homoeolog_group"]
        assert row["c_label"] == f"{group}c"
        assert row["e_label"] == f"{group}e"
        assert target[row["c_subgenome_seqid"]] == row["c_label"]
        assert target[row["e_subgenome_seqid"]] == row["e_label"]
        members.extend((row["c_subgenome_seqid"], row["e_subgenome_seqid"]))
    assert len(members) == len(set(members)) == 22
    assert set(members) == set(target)


def test_gardenia_aliases_are_bijective_and_match_primary_table() -> None:
    primary = read_rows(
        ROOT / "config/primary_seqids/gardenia_jasminoides_asm1310374v1.tsv"
    )
    aliases = read_rows(GARDENIA_ALIASES)
    assert tuple(aliases[0]) == (
        "gff_seqid",
        "genome_seqid",
        "chromosome_label",
    )
    assert len(aliases) == 11
    assert len({row["gff_seqid"] for row in aliases}) == 11
    assert len({row["genome_seqid"] for row in aliases}) == 11
    assert {(row["genome_seqid"], row["chromosome_label"]) for row in aliases} == {
        (row["seqid"], row["chromosome_label"]) for row in primary
    }


def test_coffea_policy_and_event_bind_known_subgenome_profile() -> None:
    policy = read_field_table(POLICY)
    event = read_field_table(EVENT)
    frozen = asdict(FIXED_KNOWN_SUBGENOME_CORE_H1_SCIENTIFIC_PARAMETERS)
    assert policy["model_version"] == KNOWN_SUBGENOME_CORE_H1_MODEL_VERSION
    assert policy["protocol_profile"] == frozen["protocol_profile"]
    assert policy["truth_pair_event_discriminator"] == frozen[
        "truth_pair_event_discriminator"
    ]
    assert policy["truth_pair_yn00_ks_policy"] == frozen[
        "truth_pair_yn00_ks_policy"
    ]
    assert policy["truth_sampler_seed"] == event["event_sampling_seed"] == "20260911"
    assert policy["minimum_target_chromosomes"] == event[
        "minimum_target_chromosomes"
    ] == "17"
    assert policy["minimum_homoeolog_groups"] == event[
        "minimum_homoeolog_groups"
    ] == "9"
    assert policy["H1_bootstrap_seed"] == "20260912"
    assert policy["bootstrap_replicates"] == str(frozen["bootstrap_replicates"])
    assert policy["H2_or_topology_ranking"] == "forbidden"
    assert event["self_stream_yn00_ks_policy"] == (
        "descriptive_only_not_used_for_selection"
    )
    assert event["formal_support_groups_required"] == (
        "self_WGDI,Gardenia_jasminoides,Ophiorrhiza_pumila"
    )


def test_selection_record_preserves_metadata_only_boundary() -> None:
    text = (
        ROOT / "docs/COFFEA_EXTERNAL_SELECTION_AND_CONTAMINATION_RATIONALE_v1.0.md"
    ).read_text(encoding="utf-8")
    required = (
        "No Coffea\nWGDI pair, candidate pool, truth label, score, event count, or performance",
        "*C. canephora* DH200-94 v2 was rejected at 87.09%",
        "same-species cultivar projection would make the primary\n  test easier",
        "untouched prospective\nexternal H1 replication",
    )
    for phrase in required:
        assert phrase in text


def test_coffea_no_ranker_sentinel_is_exact() -> None:
    root = ROOT / "config/holdouts/coffea_et39_v1.0/no_ranker"
    manifest_path = root / "composite_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    assert manifest == {
        "schema_version": "ploidypatch.no_ranker_known_subgenome_h1.v1.0",
        "model_version": KNOWN_SUBGENOME_CORE_H1_MODEL_VERSION,
        "ranker_enabled": False,
        "h2_or_topology_ranking_enabled": False,
        "primary_hypothesis": "H1_retain_distinct_vs_suppress_overlap_only",
        "event_discriminator": (
            "predeclared_homoeolog_group_and_subgenome_labels"
        ),
        "yn00_ks_policy": "descriptive_only_not_used_for_selection",
        "retired_ranker": "ploidypatch.stable_reference_ranker.v0.9",
        "retirement_status": "retire_v09_ranker_keep_chain_workflow",
        "automatic_approval": False,
        "truth_access": False,
    }
    expected = hashlib.sha256(manifest_path.read_bytes()).hexdigest()
    assert (root / "SHA256SUMS").read_text(encoding="utf-8") == (
        f"{expected}  composite_manifest.json\n"
    )


def test_coffea_source_manifest_is_exact_complete_and_role_separated() -> None:
    rows = read_rows(SOURCE_MANIFEST)
    assert tuple(rows[0]) == (
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
    assert len(rows) == 15
    assert len(
        {(row["role"], row["species_id"], row["artifact"]) for row in rows}
    ) == 15
    assert {row["artifact"] for row in rows} == {"genome", "gff3", "protein"}
    assert {
        role: len({row["species_id"] for row in rows if row["role"] == role})
        for role in ("target", "candidate_reference", "evaluator_reference")
    } == {"target": 1, "candidate_reference": 2, "evaluator_reference": 2}
    assert all(int(row["bytes"]) > 0 for row in rows)
    assert all(len(row["sha256"]) == 64 for row in rows)
    assert all(row["source_path"].startswith("/") for row in rows)
    assert all(row["primary_seqid_table"] for row in rows)

    containers = [row for row in rows if row["container_format"] != "direct"]
    assert len(containers) == 1
    assert containers[0]["species_id"] == "Ophiorrhiza_pumila"
    assert containers[0]["artifact"] == "genome"
    assert containers[0]["container_format"] == "tar.gz"
    assert containers[0]["member_name"] == "Ophiorrhiza_pumila.genome.fa"
    assert containers[0]["member_bytes"] == "447657994"
    assert containers[0]["member_sha256"] == (
        "3c73323c673b100b1d9c4fc4bbbfb5dba3b00d9d7b7cb42da3d6e1496a68b3e7"
    )


def test_coffea_contract_binds_every_manifest_artifact_byte_for_byte() -> None:
    contract = load_holdout_contract(CONTRACT)
    rows = read_rows(SOURCE_MANIFEST)
    manifest = {
        (row["role"], row["species_id"], row["artifact"]): row for row in rows
    }
    assert contract.holdout_id == "coffea_et39_v1.0"
    assert contract.policy_id == "ploidypatch_coffea_external_core_h1_v1.0"
    assert contract.model_version == KNOWN_SUBGENOME_CORE_H1_MODEL_VERSION
    assert contract.seeds.truth_sampler == 20260911
    assert contract.seeds.h1_bootstrap == 20260912
    assert contract.target_resolved_parameters.primary_chromosome_count == 22
    assert contract.target_resolved_parameters.minimum_target_chromosomes == 17

    observed = set()
    for reference in contract.references:
        for artifact_name, artifact in reference.artifact_items():
            key = (reference.role, reference.species_id, artifact_name)
            row = manifest[key]
            observed.add(key)
            assert str(artifact.source_relative_path) == row["source_path"].lstrip("/")
            assert artifact.bytes == int(row["bytes"])
            assert artifact.sha256 == row["sha256"]
            assert str(reference.primary_seqid_table) == row["primary_seqid_table"]
            if row["container_format"] == "tar.gz":
                assert artifact.container is not None
                assert artifact.container.format == "tar.gz"
                assert str(artifact.container.member_name) == row["member_name"]
                assert artifact.container.member_bytes == int(row["member_bytes"])
                assert artifact.container.member_sha256 == row["member_sha256"]
            else:
                assert artifact.container is None
    assert observed == set(manifest)


def test_coffea_metadata_preflight_is_truth_blind_and_no_ranker() -> None:
    report = PREFLIGHT.run_preflight(project_root=ROOT, contract_path=CONTRACT)
    assert report == {
        "schema_version": PREFLIGHT.SCHEMA_VERSION,
        "holdout_id": "coffea_et39_v1.0",
        "policy_id": "ploidypatch_coffea_external_core_h1_v1.0",
        "model_version": KNOWN_SUBGENOME_CORE_H1_MODEL_VERSION,
        "protocol_profile": "core_H1_known_subgenome_no_ranker",
        "references": 5,
        "artifacts": 15,
        "target_primary_chromosomes": 22,
        "target_homoeolog_groups": 11,
        "minimum_target_chromosomes": 17,
        "minimum_homoeolog_groups": 9,
        "truth_access": False,
        "wgd_pairs_enumerated": False,
        "candidate_counts_computed": False,
        "truth_labels_accessed": False,
        "staged_inputs_verified": False,
        "staged_code_commit": None,
    }
