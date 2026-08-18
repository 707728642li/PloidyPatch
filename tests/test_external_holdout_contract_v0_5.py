from __future__ import annotations

from dataclasses import asdict
import csv
import json
from pathlib import Path

import pytest

from ploidypatch import artifact_manifest
from ploidypatch.holdout_contract import (
    FIXED_KNOWN_SUBGENOME_CORE_H1_SCIENTIFIC_PARAMETERS,
    FIXED_SCIENTIFIC_PARAMETERS,
    KNOWN_SUBGENOME_CORE_H1_MODEL_VERSION,
    KnownSubgenomeCoreH1ScientificParameters,
    SCHEMA_VERSION,
    TRUTH_BLIND_DECLARATIONS,
    load_holdout_contract,
    staged_relative_path,
)


def reference(
    role: str,
    species: str,
    release: str,
    bundle: str,
    prefix: str,
) -> dict[str, object]:
    return {
        "role": role,
        "species_id": species,
        "release": release,
        "bundle_id": bundle,
        "wgdi_prefix": prefix,
        "primary_seqid_table": f"config/primary_seqids/{bundle}.tsv",
        "artifacts": {
            artifact: {
                "source_relative_path": f"raw/{species}/{artifact}.{suffix}",
                "bytes": index + 10,
                "sha256": f"{index + 1:x}" * 64,
            }
            for index, (artifact, suffix) in enumerate(
                (("genome", "fa.gz"), ("gff3", "gff3.gz"), ("protein", "faa.gz"))
            )
        },
    }


def valid_payload() -> dict[str, object]:
    return {
        "schema_version": SCHEMA_VERSION,
        "holdout_id": "actinidia_red5_v0.5",
        "policy_id": "ploidypatch_actinidia_external_validation_v0.5",
        "test_role": "target_level_predeclared_untouched_secondary_replication",
        "model_version": "PloidyPatch_ranker_v0.4",
        "references": [
            reference(
                "target", "Actinidia_chinensis", "Red5_v1.0", "target_red5", "ach"
            ),
            reference(
                "candidate_reference",
                "Actinidia_eriantha",
                "release_v1",
                "candidate_aeriantha",
                "aer",
            ),
            reference(
                "candidate_reference",
                "Actinidia_rufa",
                "release_v1",
                "candidate_arufa",
                "aru",
            ),
            reference(
                "evaluator_reference",
                "Rhododendron_simsii",
                "release_v1",
                "evaluator_rsimsii",
                "rsi",
            ),
            reference(
                "evaluator_reference",
                "Diospyros_oleifera",
                "release_v1",
                "evaluator_doleifera",
                "dol",
            ),
        ],
        "seeds": {
            "truth_sampler": 20261101,
            "h1_bootstrap": 20261102,
            "h2_bootstrap": 20261103,
            "guard_v03_bootstrap": 20261104,
        },
        "target_resolved_parameters": {
            "primary_chromosome_count": 29,
            "minimum_target_chromosomes_fraction": 0.75,
            "minimum_target_chromosomes": 22,
        },
        "scientific_parameters": asdict(FIXED_SCIENTIFIC_PARAMETERS),
        "truth_blind": dict(TRUTH_BLIND_DECLARATIONS),
    }


def write_contract(path: Path, payload: dict[str, object] | None = None) -> Path:
    path.write_text(
        json.dumps(valid_payload() if payload is None else payload, indent=2) + "\n",
        encoding="utf-8",
    )
    return path


def known_subgenome_h1_payload() -> dict[str, object]:
    payload = valid_payload()
    payload["holdout_id"] = "coffea_et39_v1.0"
    payload["policy_id"] = "ploidypatch_coffea_external_core_h1_v1.0"
    payload["test_role"] = "untouched_confirmatory_external_species"
    payload["model_version"] = KNOWN_SUBGENOME_CORE_H1_MODEL_VERSION
    payload["seeds"] = {
        "truth_sampler": 20260831,
        "h1_bootstrap": 20260901,
    }
    payload["target_resolved_parameters"] = {
        "primary_chromosome_count": 22,
        "minimum_target_chromosomes_fraction": 0.75,
        "minimum_target_chromosomes": 17,
    }
    payload["scientific_parameters"] = asdict(
        FIXED_KNOWN_SUBGENOME_CORE_H1_SCIENTIFIC_PARAMETERS
    )
    return payload


def test_valid_contract_is_strongly_typed_and_role_separated(tmp_path: Path) -> None:
    contract = load_holdout_contract(write_contract(tmp_path / "contract.json"))
    assert contract.target.species_id == "Actinidia_chinensis"
    assert len(contract.references_for_role("candidate_reference")) == 2
    assert len(contract.references_for_role("evaluator_reference")) == 2
    assert contract.target_resolved_parameters.primary_chromosome_count == 29
    assert contract.target_resolved_parameters.minimum_target_chromosomes == 22
    assert (
        staged_relative_path(contract.target, "genome").parts[0]
        == "shared_target"
    )
    assert staged_relative_path(contract.target, "gff3").as_posix().startswith(
        "evaluator_only/target_complete/"
    )
    assert staged_relative_path(
        contract.references_for_role("candidate_reference")[0], "protein"
    ).as_posix().startswith("candidate_only/")


def test_known_subgenome_core_h1_profile_is_strongly_typed(tmp_path: Path) -> None:
    payload = known_subgenome_h1_payload()
    contract = load_holdout_contract(
        write_contract(tmp_path / "coffea_contract.json", payload)
    )
    assert contract.model_version == KNOWN_SUBGENOME_CORE_H1_MODEL_VERSION
    assert contract.seeds.truth_sampler == 20260831
    assert contract.target_resolved_parameters.minimum_target_chromosomes == 17
    assert isinstance(
        contract.scientific_parameters,
        KnownSubgenomeCoreH1ScientificParameters,
    )
    assert contract.scientific_parameters.truth_pair_yn00_ks_policy == (
        "descriptive_only_not_used_for_selection"
    )
    assert contract.scientific_parameters.truth_pair_target_subgenome_pairing == (
        "same_group_exactly_one_member_from_each_subgenome"
    )


def test_known_subgenome_core_h1_profile_rejects_event_rule_drift(
    tmp_path: Path,
) -> None:
    payload = known_subgenome_h1_payload()
    payload["scientific_parameters"][  # type: ignore[index]
        "truth_pair_yn00_ks_policy"
    ] = "selection_threshold"
    with pytest.raises(ValueError, match="truth_pair_yn00_ks_policy"):
        load_holdout_contract(
            write_contract(tmp_path / "coffea_contract.json", payload)
        )


@pytest.mark.parametrize(
    ("role_index", "replacement_role", "message"),
    (
        (4, "candidate_reference", "exactly one target"),
        (2, "evaluator_reference", "exactly one target"),
        (0, "candidate_reference", "exactly one target"),
    ),
)
def test_contract_requires_exactly_one_two_two_roles(
    tmp_path: Path, role_index: int, replacement_role: str, message: str
) -> None:
    payload = valid_payload()
    payload["references"][role_index]["role"] = replacement_role  # type: ignore[index]
    with pytest.raises(ValueError, match=message):
        load_holdout_contract(write_contract(tmp_path / "contract.json", payload))


@pytest.mark.parametrize("field", ("species_id", "bundle_id", "wgdi_prefix"))
def test_contract_rejects_cross_role_identifier_reuse(
    tmp_path: Path, field: str
) -> None:
    payload = valid_payload()
    references = payload["references"]  # type: ignore[assignment]
    references[4][field] = references[0][field]  # type: ignore[index]
    with pytest.raises(ValueError, match=field):
        load_holdout_contract(write_contract(tmp_path / "contract.json", payload))


@pytest.mark.parametrize(
    "unsafe",
    (
        "/absolute/genome.fa",
        "../escape/genome.fa",
        "raw/../escape.fa",
        "raw\\windows\\genome.fa",
        "C:/drive/genome.fa",
        "raw//genome.fa",
        "./raw/genome.fa",
    ),
)
def test_contract_rejects_unsafe_artifact_paths(tmp_path: Path, unsafe: str) -> None:
    payload = valid_payload()
    payload["references"][0]["artifacts"]["genome"][  # type: ignore[index]
        "source_relative_path"
    ] = unsafe
    with pytest.raises(ValueError, match="Unsafe|canonical"):
        load_holdout_contract(write_contract(tmp_path / "contract.json", payload))


def test_artifact_schema_structurally_prevents_role_or_release_drift(
    tmp_path: Path,
) -> None:
    payload = valid_payload()
    artifact = payload["references"][1]["artifacts"]["genome"]  # type: ignore[index]
    artifact["role"] = "evaluator_reference"
    artifact["release"] = "different"
    with pytest.raises(ValueError, match="fields differ"):
        load_holdout_contract(write_contract(tmp_path / "contract.json", payload))


def test_contract_rejects_scientific_parameter_drift(tmp_path: Path) -> None:
    payload = valid_payload()
    payload["scientific_parameters"]["adapter_min_identity"] = 0.49  # type: ignore[index]
    with pytest.raises(ValueError, match="adapter_min_identity"):
        load_holdout_contract(write_contract(tmp_path / "contract.json", payload))


def test_contract_rejects_truth_blind_relaxation(tmp_path: Path) -> None:
    payload = valid_payload()
    payload["truth_blind"]["blind_evaluator_reference_mount"] = True  # type: ignore[index]
    with pytest.raises(ValueError, match="blind_evaluator_reference_mount"):
        load_holdout_contract(write_contract(tmp_path / "contract.json", payload))


def test_target_chromosome_gate_is_resolved_as_ceil_75_percent(
    tmp_path: Path,
) -> None:
    payload = valid_payload()
    payload["target_resolved_parameters"]["minimum_target_chromosomes"] = 15  # type: ignore[index]
    with pytest.raises(ValueError, match=r"ceil\(0.75 \* 29\) = 22"):
        load_holdout_contract(write_contract(tmp_path / "contract.json", payload))

    payload = valid_payload()
    payload["target_resolved_parameters"][  # type: ignore[index]
        "minimum_target_chromosomes_fraction"
    ] = 0.5
    with pytest.raises(ValueError, match="frozen 0.75"):
        load_holdout_contract(write_contract(tmp_path / "contract2.json", payload))


def test_contract_requires_unique_randomization_seeds(tmp_path: Path) -> None:
    payload = valid_payload()
    payload["seeds"]["h2_bootstrap"] = payload["seeds"]["h1_bootstrap"]  # type: ignore[index]
    with pytest.raises(ValueError, match="unique seed"):
        load_holdout_contract(write_contract(tmp_path / "contract.json", payload))


def test_contract_rejects_duplicate_json_keys(tmp_path: Path) -> None:
    path = tmp_path / "duplicate.json"
    path.write_text('{"schema_version":"x","schema_version":"y"}\n', encoding="utf-8")
    with pytest.raises(ValueError, match="Duplicate JSON object key"):
        load_holdout_contract(path)


def test_real_actinidia_contract_loads_and_matches_frozen_sources_and_policy() -> None:
    root = Path(__file__).parents[1]
    contract_path = (
        root / "config/holdouts/actinidia_red5_v0.5/contract.json"
    )
    contract = load_holdout_contract(contract_path)
    policy_path = root / "config/actinidia_external_validation_policy_v0.5.tsv"
    with policy_path.open(encoding="utf-8", newline="") as handle:
        policy_rows = list(csv.DictReader(handle, delimiter="\t"))
    policy = {row["field"]: row["value"] for row in policy_rows}
    assert len(policy) == len(policy_rows)
    assert contract.policy_id == policy["policy_id"]
    assert contract.model_version == policy["model_version"] == policy["estimator_version"]
    assert contract.test_role == "target_level_predeclared_untouched_secondary_replication"
    assert policy["test_role"] == (
        "target_level_predeclared_untouched_secondary_replication"
    )
    assert contract.seeds.truth_sampler == int(policy["truth_sampler_seed"])
    assert contract.seeds.h1_bootstrap == int(policy["H1_bootstrap_seed"])
    assert contract.seeds.h2_bootstrap == int(policy["H2_bootstrap_seed"])
    assert contract.seeds.guard_v03_bootstrap == int(
        policy["guard_v03_bootstrap_seed"]
    )
    target = contract.target_resolved_parameters
    assert target.primary_chromosome_count == int(
        policy["target_primary_chromosome_count"]
    ) == 29
    assert target.minimum_target_chromosomes_fraction == float(
        policy["minimum_target_chromosome_fraction"]
    ) == 0.75
    assert target.minimum_target_chromosomes == int(
        policy["minimum_target_chromosomes"]
    ) == 22

    science = contract.scientific_parameters
    assert ",".join(science.candidate_method_families) == policy[
        "candidate_method_families"
    ]
    for name in (
        "multiple_references_per_method_vote",
        "candidate_topology_identity",
        "primary_candidate_policy",
        "legacy_candidate_comparator",
        "truth_removal_policy",
        "truth_sampler_balance",
    ):
        assert str(getattr(science, name)) == policy[name]
    assert science.truth_pair_self_wgdi_require_cross_seqid is True
    assert policy["truth_pair_self_wgdi_require_cross_primary_chromosome"] == "true"
    assert policy["truth_pair_outgroup_counterpart_multiplicity"].startswith(
        science.truth_pair_outgroup_counterpart_multiplicity
    )
    assert science.truth_pair_outgroup_require_reciprocal_unique is True
    assert policy["truth_pair_outgroup_require_reciprocal_unique"].startswith("true")
    assert science.truth_pair_final_rule == (
        "exact_unordered_pair_intersection_of_self_wgdi_and_two_outgroup_support"
    )
    assert policy["truth_pair_final_rule"] == (
        "exact_unordered_pair_intersection_of_target_self_wgdi_and_two_evaluator_groups"
    )
    for name in (
        "adapter_min_identity",
        "adapter_min_query_coverage",
        "adapter_max_existing_cds_overlap",
        "adapter_max_redundancy_overlap",
        "minimum_topology_coverage_among_positive_candidates",
        "minimum_v03_AP_gain_retained_fraction",
    ):
        assert float(getattr(science, name)) == float(policy[name])
    assert science.truth_pair_self_wgdi_min_block_pairs == int(
        policy["truth_pair_self_wgdi_min_block_pairs"]
    )
    assert science.truth_pair_outgroup_min_block_pairs == int(
        policy["truth_pair_outgroup_min_block_pairs"]
    )
    assert science.truth_pair_outgroup_min_support_groups == int(
        policy["truth_pair_outgroup_min_support_groups"]
    )
    assert science.truth_event_count == 800
    assert policy["truth_event_count"].startswith("800_")
    assert science.minimum_formal_event_count == int(
        policy["minimum_formal_event_count"]
    )
    assert science.minimum_events_per_complexity_bin == int(
        policy["minimum_events_per_complexity_bin"]
    )
    assert science.bootstrap_replicates == int(policy["bootstrap_replicates"])
    assert science.minimum_chromosome_bootstrap_valid_replicates == int(
        policy["minimum_chromosome_bootstrap_valid_replicates"]
    )
    assert science.automatic_copy_addition_approval is False
    assert policy["automatic_copy_addition_approval"] == "false"
    assert {
        reference.bundle_id: reference.wgdi_prefix
        for reference in contract.references
    } == {
        "target_red5": "red5",
        "candidate_eriantha": "aer",
        "candidate_rufa": "aru",
        "evaluator_rhododendron": "rhs",
        "evaluator_diospyros": "dol",
    }
    truth_blind = contract.truth_blind
    assert (
        truth_blind[
            "selection_by_pair_yield_candidate_count_or_model_performance"
        ]
        is False
    )
    assert policy[
        "selection_by_pair_yield_candidate_count_label_or_model_performance"
    ] == "forbidden"
    assert truth_blind["reference_role_change_after_freeze"] is False
    assert policy["reference_role_change_after_freeze"] == "forbidden"
    assert truth_blind["wgd_pairs_enumerated_before_protocol_freeze"] is False
    assert policy["preflight_wgd_pairs_enumerated"] == "false"
    assert truth_blind["candidate_counts_computed_before_protocol_freeze"] is False
    assert policy["preflight_candidate_counts_computed"] == "false"
    assert truth_blind["truth_labels_accessed_before_protocol_freeze"] is False
    assert policy["preflight_truth_labels_accessed"] == "false"
    for contract_name, policy_name in (
        ("blind_truth_mount", "blind_runner_truth_mount"),
        (
            "blind_complete_target_annotation_mount",
            "blind_runner_complete_target_annotation_mount",
        ),
        (
            "blind_evaluator_reference_mount",
            "blind_runner_evaluator_reference_mount",
        ),
        ("blind_nas_data_mount", "blind_runner_nas_data_mount"),
        ("blind_network_access", "blind_runner_network_access"),
    ):
        assert truth_blind[contract_name] is False
        assert policy[policy_name] == "false"
    assert truth_blind[
        "complete_control_generated_after_blind_raw_prediction_freeze"
    ] is True
    assert policy["paired_complete_annotation_control"] == (
        "required_evaluator_generated_after_blind_raw_prediction_freeze"
    )

    source_table = root / "config/actinidia_external_input_sources_v0.5.tsv"
    with source_table.open(encoding="utf-8", newline="") as handle:
        source_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(source_rows) == 15
    source_root = "/nas_data/NFS/Public_genome_data/"
    observed = {
        (row["role"], row["species_id"], row["release"], row["artifact"]): row
        for row in source_rows
    }
    assert len(observed) == 15
    for reference in contract.references:
        primary = root / reference.primary_seqid_table.as_posix()
        assert primary.is_file() and not primary.is_symlink()
        for artifact_name, artifact in reference.artifact_items():
            row = observed[(
                reference.role,
                reference.species_id,
                reference.release,
                artifact_name,
            )]
            assert row["source_path"].startswith(source_root)
            assert row["source_path"][len(source_root) :] == (
                artifact.source_relative_path.as_posix()
            )
            assert int(row["bytes"]) == artifact.bytes
            assert row["sha256"] == artifact.sha256
    target_primary = root / contract.target.primary_seqid_table.as_posix()
    with target_primary.open(encoding="utf-8", newline="") as handle:
        primary_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(primary_rows) == target.primary_chromosome_count == 29


def test_strict_manifest_verifies_exact_nested_file_universe(tmp_path: Path) -> None:
    root = tmp_path / "artifacts"
    (root / "nested").mkdir(parents=True)
    (root / "a.txt").write_text("alpha\n", encoding="utf-8")
    (root / "nested/b.txt").write_text("beta\n", encoding="utf-8")
    manifest = artifact_manifest.write_sha256sums(root)
    entries = artifact_manifest.verify_sha256sums(
        root, manifest, ignore_checksum_file=True
    )
    assert set(entries) == {"a.txt", "nested/b.txt"}


def test_strict_manifest_rejects_missing_extra_and_mutated_files(
    tmp_path: Path,
) -> None:
    root = tmp_path / "artifacts"
    root.mkdir()
    item = root / "item.txt"
    item.write_text("original\n", encoding="utf-8")
    manifest = artifact_manifest.write_sha256sums(root)

    item.write_text("mutated!\n", encoding="utf-8")
    with pytest.raises(ValueError, match="SHA256 mismatch"):
        artifact_manifest.verify_sha256sums(root, manifest, ignore_checksum_file=True)

    item.write_text("original\n", encoding="utf-8")
    (root / "extra.txt").write_text("extra\n", encoding="utf-8")
    with pytest.raises(ValueError, match="extra=.*extra.txt"):
        artifact_manifest.verify_sha256sums(root, manifest, ignore_checksum_file=True)

    (root / "extra.txt").unlink()
    item.unlink()
    with pytest.raises(ValueError, match="missing=.*item.txt"):
        artifact_manifest.verify_sha256sums(root, manifest, ignore_checksum_file=True)


@pytest.mark.parametrize(
    "line",
    (
        f"{'a' * 64}  /absolute.txt\n",
        f"{'a' * 64}  ../escape.txt\n",
        f"{'a' * 64}  dir\\windows.txt\n",
        f"{'a' * 64}  ./noncanonical.txt\n",
    ),
)
def test_strict_manifest_rejects_unsafe_paths(tmp_path: Path, line: str) -> None:
    manifest = tmp_path / "SHA256SUMS"
    manifest.write_text(line, encoding="utf-8")
    with pytest.raises(ValueError, match="Unsafe|canonical"):
        artifact_manifest.read_sha256sums(manifest)


def test_strict_manifest_rejects_duplicate_and_case_colliding_paths(
    tmp_path: Path,
) -> None:
    manifest = tmp_path / "SHA256SUMS"
    manifest.write_text(
        f"{'a' * 64}  item.txt\n{'b' * 64}  ITEM.txt\n", encoding="utf-8"
    )
    with pytest.raises(ValueError, match="Duplicate"):
        artifact_manifest.read_sha256sums(manifest)


def test_ignored_checksum_file_must_not_list_itself(tmp_path: Path) -> None:
    root = tmp_path / "root"
    root.mkdir()
    (root / "item").write_text("x", encoding="utf-8")
    manifest = root / "SHA256SUMS"
    manifest.write_text(
        f"{artifact_manifest.sha256_file(root / 'item')}  item\n"
        f"{'0' * 64}  SHA256SUMS\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="must not list itself"):
        artifact_manifest.verify_sha256sums(
            root, manifest, ignore_checksum_file=True
        )


def test_strict_manifest_rejects_symlinks(tmp_path: Path) -> None:
    root = tmp_path / "root"
    root.mkdir()
    target = root / "target"
    target.write_text("x", encoding="utf-8")
    link = root / "link"
    try:
        link.symlink_to(target.name)
    except OSError:
        pytest.skip("Host does not permit symlink creation")
    manifest = root / "SHA256SUMS"
    manifest.write_text(
        f"{artifact_manifest.sha256_file(target)}  target\n", encoding="utf-8"
    )
    with pytest.raises(ValueError, match="Symlink"):
        artifact_manifest.verify_sha256sums(
            root, manifest, ignore_checksum_file=True
        )
