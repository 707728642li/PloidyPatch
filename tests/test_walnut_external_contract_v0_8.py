from __future__ import annotations

import csv
import importlib.util
import json
from pathlib import Path
import shutil
import sys

import pytest

from ploidypatch.artifact_manifest import (
    sha256_file,
    verify_sha256sums,
    write_sha256sums,
)
from ploidypatch.holdout_contract import (
    CORE_H1_MODEL_VERSION,
    CoreH1ScientificParameters,
    CoreH1Seeds,
    load_holdout_contract,
    staged_relative_path,
)


ROOT = Path(__file__).parents[1]
CONTRACT = ROOT / "config/holdouts/walnut_walnut2_v0.8/contract.json"
NO_RANKER = ROOT / "config/holdouts/walnut_walnut2_v0.8/no_ranker"


def load_script(relative: str, name: str):
    path = ROOT / relative
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


PREFLIGHT = load_script(
    "scripts/preflight_walnut_external_inputs_v0.8.py", "walnut_preflight_v08"
)
FREEZE = load_script(
    "scripts/freeze_walnut_external_protocol_v0.8.py", "walnut_freeze_v08"
)
GENERIC_FREEZE = load_script(
    "scripts/freeze_external_holdout_protocol_v0.5.py", "generic_freeze_for_walnut"
)


def read_manifest() -> list[dict[str, str]]:
    with (ROOT / "config/walnut_external_input_sources_v0.8.tsv").open(
        encoding="utf-8", newline=""
    ) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def mutate_contract(tmp_path: Path, dotted: tuple[str | int, ...], value: object) -> Path:
    payload = json.loads(CONTRACT.read_text(encoding="utf-8"))
    current: object = payload
    for key in dotted[:-1]:
        current = current[key]  # type: ignore[index]
    current[dotted[-1]] = value  # type: ignore[index]
    path = tmp_path / "contract.json"
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return path


def test_real_walnut_contract_is_exact_h1_only_and_role_separated() -> None:
    contract = load_holdout_contract(CONTRACT)
    assert contract.holdout_id == "walnut_walnut2_v0.8"
    assert contract.model_version == CORE_H1_MODEL_VERSION
    assert isinstance(contract.seeds, CoreH1Seeds)
    assert set(vars(contract.seeds)) == {"truth_sampler", "h1_bootstrap"}
    assert isinstance(contract.scientific_parameters, CoreH1ScientificParameters)
    science = contract.scientific_parameters
    assert science.protocol_profile == "core_H1_only_no_ranker"
    assert science.h2_or_topology_ranking == "forbidden"
    assert science.retired_ranker == "ploidypatch.stable_reference_ranker.v0.9"
    assert science.bootstrap_replicates == 20_000
    assert science.bootstrap_unit == "paired_event"
    assert science.truth_pair_yn00_ks_minimum == 0.10
    assert science.truth_pair_yn00_ks_maximum == 0.75
    assert science.truth_pair_missing_or_out_of_range_ks == "abstain"
    assert science.all_arm_collateral_loss_maximum == 0
    assert contract.target_resolved_parameters.primary_chromosome_count == 16
    assert contract.target_resolved_parameters.minimum_target_chromosomes == 12
    assert [len(contract.references_for_role(role)) for role in (
        "target", "candidate_reference", "evaluator_reference"
    )] == [1, 2, 2]
    assert {reference.species_id for reference in contract.references_for_role(
        "candidate_reference"
    )} == {"Juglans_mandshurica", "Carya_illinoinensis"}
    assert {reference.species_id for reference in contract.references_for_role(
        "evaluator_reference"
    )} == {"Corylus_avellana", "Castanea_mollissima"}


def test_real_contract_matches_all_metadata_sources_and_tar_member() -> None:
    contract = load_holdout_contract(CONTRACT)
    rows = {(row["species_id"], row["artifact"]): row for row in read_manifest()}
    assert len(rows) == 15
    for reference in contract.references:
        for artifact_name, artifact in reference.artifact_items():
            row = rows[(reference.species_id, artifact_name)]
            assert row["source_path"] == "/" + artifact.source_relative_path.as_posix()
            assert row["bytes"] == str(artifact.bytes)
            assert row["sha256"] == artifact.sha256
            assert row["role"] == reference.role
            assert row["release"] == reference.release
    mandshurica = next(
        reference for reference in contract.references
        if reference.species_id == "Juglans_mandshurica"
    )
    genome = mandshurica.genome
    assert genome.container is not None
    assert genome.container.format == "tar.gz"
    assert genome.staged_filename == "Juglans_mandshurica.genome.fa"
    assert genome.staged_bytes == 554_212_158
    assert genome.staged_sha256 == (
        "7b1be8bfb34096526ac94bc6f1806f0257c7b31293b81590225d283f8949679f"
    )
    assert staged_relative_path(mandshurica, "genome").name == genome.staged_filename


@pytest.mark.parametrize(
    ("path", "value", "message"),
    (
        (("scientific_parameters", "h2_or_topology_ranking"), "enabled", "H1-only"),
        (("scientific_parameters", "truth_pair_yn00_ks_maximum"), 0.8, "H1-only"),
        (("scientific_parameters", "all_arm_collateral_loss_maximum"), 1, "H1-only"),
        (("seeds", "h2_bootstrap"), 20260823, "fields differ"),
        (("references", 1, "artifacts", "genome", "container", "member_name"),
         "../escape.fa", "Unsafe"),
        (("references", 1, "artifacts", "gff3", "container"),
         {
             "format": "tar.gz",
             "member_name": "annotation.fa",
             "member_bytes": 1,
             "member_sha256": "0" * 64,
         }, "Only a genome"),
    ),
)
def test_walnut_contract_rejects_h2_scientific_and_tar_drift(
    tmp_path: Path, path: tuple[str | int, ...], value: object, message: str
) -> None:
    with pytest.raises(ValueError, match=message):
        load_holdout_contract(mutate_contract(tmp_path, path, value))


def test_walnut_metadata_preflight_is_truth_blind_and_no_ranker() -> None:
    report = PREFLIGHT.run_preflight(
        project_root=ROOT,
        contract_path=CONTRACT,
    )
    assert report == {
        "schema_version": PREFLIGHT.SCHEMA_VERSION,
        "holdout_id": "walnut_walnut2_v0.8",
        "policy_id": "ploidypatch_walnut_external_core_h1_v0.8",
        "model_version": CORE_H1_MODEL_VERSION,
        "protocol_profile": "core_H1_only_no_ranker",
        "references": 5,
        "artifacts": 15,
        "target_primary_chromosomes": 16,
        "minimum_target_chromosomes": 12,
        "truth_access": False,
        "wgd_pairs_enumerated": False,
        "candidate_counts_computed": False,
        "truth_labels_accessed": False,
        "staged_inputs_verified": False,
        "staged_code_commit": None,
    }
    verify_sha256sums(NO_RANKER, ignore_checksum_file=True)
    manifest = json.loads((NO_RANKER / "composite_manifest.json").read_text())
    assert manifest["ranker_enabled"] is False
    assert manifest["h2_or_topology_ranking_enabled"] is False
    with pytest.raises(ValueError, match="v0.4 ranker"):
        GENERIC_FREEZE.verify_model(NO_RANKER, load_holdout_contract(CONTRACT))


def test_generic_freezer_rejects_ranker_for_walnut_core_h1(tmp_path: Path) -> None:
    model = tmp_path / "model"
    model.mkdir()
    (model / "composite_manifest.json").write_text(
        json.dumps(
            {
                "schema_version": "ploidypatch.composite_ranker.v0.4",
                "automatic_approval": False,
            }
        ) + "\n",
        encoding="utf-8",
    )
    write_sha256sums(model)
    with pytest.raises(ValueError, match="v0.4 ranker"):
        GENERIC_FREEZE.verify_model(model, load_holdout_contract(CONTRACT))


def test_walnut_freeze_skeleton_binds_h1_protocol_inputs_and_code(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    code_root = tmp_path / "code"
    for relative in (*FREEZE.PROTOCOL_FILES, *FREEZE.IMPLEMENTATION_FILES):
        source = ROOT / relative
        destination = code_root / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, destination)
    stage = tmp_path / "stage"
    stage.mkdir()
    commit = "a" * 40
    (stage / "role_contract.json").write_text(
        json.dumps({"code_commit": commit}) + "\n", encoding="utf-8"
    )
    (stage / "role_manifest.tsv").write_text("role\n", encoding="utf-8")
    write_sha256sums(stage)
    monkeypatch.setattr(FREEZE, "verify_git_state", lambda *_: None)
    monkeypatch.setattr(
        FREEZE.runpy,
        "run_path",
        lambda *_args, **_kwargs: {
            "run_preflight": lambda **_values: {"staged_inputs_verified": True}
        },
    )
    output = tmp_path / "freeze"
    FREEZE.freeze_protocol(
        project_root=tmp_path / "empty-project",
        code_root=code_root,
        contract_path=code_root / FREEZE.PROTOCOL_FILES[0],
        staged_inputs=stage,
        code_commit=commit,
        output=output,
    )
    verify_sha256sums(output, ignore_checksum_file=True)
    manifest = json.loads((output / "protocol_manifest.json").read_text())
    assert manifest["protocol_profile"] == "core_H1_only_no_ranker"
    assert manifest["ranker_enabled"] is False
    assert manifest["h2_or_topology_ranking_enabled"] is False
    assert manifest["formal_runner_frozen"] is False
    assert manifest["truth_access"] is False
    assert manifest["wgd_pairs_enumerated"] is False
    assert manifest["all_arm_collateral_loss_maximum"] == 0


def test_walnut_freeze_refuses_existing_target_derived_artifact(tmp_path: Path) -> None:
    forbidden = tmp_path / FREEZE.FORBIDDEN_TARGET_ARTIFACTS[0]
    forbidden.mkdir(parents=True)
    with pytest.raises(ValueError, match="too late"):
        FREEZE.reject_too_late(tmp_path)
