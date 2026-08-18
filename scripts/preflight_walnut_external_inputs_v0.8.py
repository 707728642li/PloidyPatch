#!/usr/bin/env python3
"""Truth-blind metadata and staged-input preflight for Walnut core H1."""
from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path, PurePosixPath
from typing import Any, Iterable

from ploidypatch.artifact_manifest import sha256_file, verify_sha256sums
from ploidypatch.holdout_contract import (
    CORE_H1_MODEL_VERSION,
    CoreH1ScientificParameters,
    HoldoutContract,
    load_holdout_contract,
    staged_relative_path,
)


SCHEMA_VERSION = "ploidypatch.walnut_external_preflight.v0.8"
HOLDOUT_ID = "walnut_walnut2_v0.8"
POLICY_ID = "ploidypatch_walnut_external_core_h1_v0.8"
STAGE_SCHEMA = "ploidypatch.external_holdout_input_stage.v0.5"
ROLE_MANIFEST_FIELDS = (
    "role",
    "species_id",
    "release",
    "bundle_id",
    "wgdi_prefix",
    "artifact",
    "bytes",
    "sha256",
    "source_relative_path",
    "staged_relative_path",
    "staged_sha256",
)


def _reject_duplicate_json_keys(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"Duplicate JSON key: {key}")
        result[key] = value
    return result


def load_json_object(path: Path) -> dict[str, Any]:
    if not path.is_file() or path.is_symlink() or path.stat().st_size == 0:
        raise ValueError(f"Missing, empty or symlinked JSON: {path}")
    value = json.loads(
        path.read_text(encoding="utf-8"),
        object_pairs_hook=_reject_duplicate_json_keys,
    )
    if not isinstance(value, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return value


def read_tsv(path: Path, expected_fields: tuple[str, ...]) -> list[dict[str, str]]:
    if not path.is_file() or path.is_symlink() or path.stat().st_size == 0:
        raise ValueError(f"Missing, empty or symlinked TSV: {path}")
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != expected_fields:
            raise ValueError(f"Unexpected TSV fields in {path}")
        rows = list(reader)
    if not rows:
        raise ValueError(f"TSV has no data rows: {path}")
    return rows


def read_field_table(path: Path) -> dict[str, str]:
    rows = read_tsv(path, ("field", "value"))
    values = {row["field"]: row["value"] for row in rows}
    if len(values) != len(rows) or any(not key or not value for key, value in values.items()):
        raise ValueError(f"Duplicate or empty field/value entry: {path}")
    return values


def verify_policy_and_event(project_root: Path, contract: HoldoutContract) -> None:
    policy = read_field_table(
        project_root / "config/walnut_external_validation_policy_v0.8.tsv"
    )
    event = read_field_table(
        project_root / "config/walnut_external_event_definition_v0.8.tsv"
    )
    required_policy = {
        "policy_id": contract.policy_id,
        "protocol_profile": "core_H1_only_no_ranker",
        "model_version": CORE_H1_MODEL_VERSION,
        "stable_reference_ranker_v0.9_status": (
            "retire_v09_ranker_keep_chain_workflow"
        ),
        "H2_or_topology_ranking": "forbidden",
        "primary_hypothesis_order": "H1_only",
        "truth_pair_yn00_ks_interval": "closed_0.10_to_0.75",
        "truth_pair_missing_or_out_of_range_ks": "abstain",
        "minimum_formal_event_count": "500",
        "target_primary_chromosome_count": "16",
        "minimum_target_chromosomes": "12",
        "minimum_events_per_complexity_bin": "20",
        "bootstrap_replicates": "20000",
        "collateral_loss_gate": (
            "zero_baseline_transcript_structure_loss_in_every_arm"
        ),
    }
    for name, expected in required_policy.items():
        if policy.get(name) != expected:
            raise ValueError(f"Walnut policy differs at {name}")
    required_event = {
        "policy_id": contract.policy_id,
        "target_event": "Juglandaceae_ancestral_juglandoid_WGD",
        "target_event_paleoallopolyploid_risk": "declared",
        "evaluator_groups_outside_target_event": "true",
        "self_stream_yn00_ks_interval": "closed_0.10_to_0.75",
        "self_stream_missing_or_out_of_range_ks": "abstain",
        "formal_support_groups_required": (
            "self_WGDI,Corylus_avellana,Castanea_mollissima"
        ),
        "minimum_formal_events": "500",
        "minimum_target_chromosomes": "12",
        "minimum_events_each_complexity_bin": "20",
    }
    for name, expected in required_event.items():
        if event.get(name) != expected:
            raise ValueError(f"Walnut event definition differs at {name}")


def verify_source_manifest(project_root: Path, contract: HoldoutContract) -> None:
    fields = (
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
    rows = read_tsv(project_root / "config/walnut_external_input_sources_v0.8.tsv", fields)
    observed = {(row["species_id"], row["artifact"]): row for row in rows}
    if len(rows) != 15 or len(observed) != 15:
        raise ValueError("Walnut source manifest must contain 15 unique artifacts")
    for reference in contract.references:
        for artifact_name, artifact in reference.artifact_items():
            key = (reference.species_id, artifact_name)
            row = observed.get(key)
            if row is None:
                raise ValueError(f"Source manifest lacks {key}")
            expected_source = "/" + artifact.source_relative_path.as_posix()
            if (
                row["role"] != reference.role
                or row["release"] != reference.release
                or row["primary_seqid_table"]
                != reference.primary_seqid_table.as_posix()
                or row["source_path"] != expected_source
                or row["bytes"] != str(artifact.bytes)
                or row["sha256"] != artifact.sha256
            ):
                raise ValueError(f"Source manifest differs from contract at {key}")
            if artifact.container is None:
                container_values = (
                    row["container_format"], row["member_name"],
                    row["member_bytes"], row["member_sha256"],
                )
                if container_values != ("direct", "", "", ""):
                    raise ValueError(f"Unexpected container metadata at {key}")
            elif (
                row["container_format"] != artifact.container.format
                or row["member_name"] != artifact.container.member_name.as_posix()
                or row["member_bytes"] != str(artifact.container.member_bytes)
                or row["member_sha256"] != artifact.container.member_sha256
            ):
                raise ValueError(f"Tar member binding differs at {key}")


def verify_primary_tables(project_root: Path, contract: HoldoutContract) -> None:
    for reference in contract.references:
        rows = read_tsv(
            project_root.joinpath(*reference.primary_seqid_table.parts),
            ("seqid", "chromosome_label"),
        )
        if len({row["seqid"] for row in rows}) != len(rows) or len(
            {row["chromosome_label"] for row in rows}
        ) != len(rows):
            raise ValueError(f"Duplicate primary seqid or label for {reference.species_id}")
        if reference is contract.target and len(rows) != 16:
            raise ValueError("Walnut target primary-seqid table must contain 16 rows")


def verify_staged_inputs(
    stage_root: Path, contract: HoldoutContract, contract_path: Path
) -> dict[str, Any]:
    verify_sha256sums(stage_root, ignore_checksum_file=True)
    role_contract_path = stage_root / "role_contract.json"
    if not role_contract_path.is_file() or role_contract_path.is_symlink():
        raise ValueError("Missing or symlinked staged role contract")
    role_contract = load_json_object(role_contract_path)
    if (
        not isinstance(role_contract, dict)
        or role_contract.get("schema_version") != STAGE_SCHEMA
        or role_contract.get("holdout_id") != contract.holdout_id
        or role_contract.get("policy_id") != contract.policy_id
        or role_contract.get("model_version") != contract.model_version
        or role_contract.get("contract", {}).get("sha256")
        != sha256_file(contract_path)
        or role_contract.get("truth_blind") != dict(contract.truth_blind)
        or role_contract.get("role_boundaries", {}).get(
            "candidate_evaluator_species_overlap"
        )
        is not False
    ):
        raise ValueError("Staged role contract differs from Walnut contract")
    rows = read_tsv(stage_root / "role_manifest.tsv", ROLE_MANIFEST_FIELDS)
    observed = {(row["species_id"], row["artifact"]): row for row in rows}
    if len(rows) != 15 or len(observed) != 15:
        raise ValueError("Staged role manifest must contain 15 unique artifacts")
    for reference in contract.references:
        for artifact_name, artifact in reference.artifact_items():
            key = (reference.species_id, artifact_name)
            row = observed.get(key)
            relative = staged_relative_path(reference, artifact_name)
            path = stage_root.joinpath(*relative.parts)
            if row is None or (
                row["role"] != reference.role
                or row["release"] != reference.release
                or row["bundle_id"] != reference.bundle_id
                or row["wgdi_prefix"] != reference.wgdi_prefix
                or row["source_relative_path"]
                != artifact.source_relative_path.as_posix()
                or row["staged_relative_path"] != relative.as_posix()
                or row["bytes"] != str(artifact.staged_bytes)
                or row["sha256"] != artifact.staged_sha256
                or row["staged_sha256"] != artifact.staged_sha256
                or not path.is_file()
                or path.is_symlink()
                or path.stat().st_size != artifact.staged_bytes
                or sha256_file(path) != artifact.staged_sha256
            ):
                raise ValueError(f"Staged artifact differs at {key}")
    return role_contract


def run_preflight(
    *,
    project_root: Path,
    contract_path: Path,
    staged_inputs: Path | None = None,
    output: Path | None = None,
) -> dict[str, Any]:
    contract = load_holdout_contract(contract_path)
    if contract.holdout_id != HOLDOUT_ID or contract.policy_id != POLICY_ID:
        raise ValueError("Preflight received a non-Walnut contract")
    if contract.model_version != CORE_H1_MODEL_VERSION or not isinstance(
        contract.scientific_parameters, CoreH1ScientificParameters
    ):
        raise ValueError("Walnut preflight requires the exact core H1-only profile")
    verify_policy_and_event(project_root, contract)
    verify_source_manifest(project_root, contract)
    verify_primary_tables(project_root, contract)
    no_ranker = project_root / "config/holdouts/walnut_walnut2_v0.8/no_ranker"
    verify_sha256sums(no_ranker, ignore_checksum_file=True)
    no_ranker_manifest = load_json_object(no_ranker / "composite_manifest.json")
    expected_no_ranker = {
        "schema_version": "ploidypatch.no_ranker_h1.v0.8",
        "model_version": CORE_H1_MODEL_VERSION,
        "ranker_enabled": False,
        "h2_or_topology_ranking_enabled": False,
        "primary_hypothesis": "H1_retain_distinct_vs_suppress_overlap_only",
        "retired_ranker": "ploidypatch.stable_reference_ranker.v0.9",
        "retirement_status": "retire_v09_ranker_keep_chain_workflow",
        "automatic_approval": False,
        "truth_access": False,
    }
    if no_ranker_manifest != expected_no_ranker:
        raise ValueError("Walnut no-ranker sentinel differs from the exact H1 profile")
    role_contract = (
        verify_staged_inputs(staged_inputs, contract, contract_path)
        if staged_inputs is not None
        else None
    )
    report = {
        "schema_version": SCHEMA_VERSION,
        "holdout_id": contract.holdout_id,
        "policy_id": contract.policy_id,
        "model_version": contract.model_version,
        "protocol_profile": "core_H1_only_no_ranker",
        "references": 5,
        "artifacts": 15,
        "target_primary_chromosomes": 16,
        "minimum_target_chromosomes": 12,
        "truth_access": False,
        "wgd_pairs_enumerated": False,
        "candidate_counts_computed": False,
        "truth_labels_accessed": False,
        "staged_inputs_verified": staged_inputs is not None,
        "staged_code_commit": (
            role_contract.get("code_commit") if role_contract is not None else None
        ),
    }
    if output is not None:
        if output.exists() or output.is_symlink():
            raise FileExistsError(f"Refusing to overwrite preflight report: {output}")
        output.parent.mkdir(parents=True, exist_ok=True)
        temporary = output.with_name(f".{output.name}.partial.{os.getpid()}")
        if temporary.exists() or temporary.is_symlink():
            raise FileExistsError(f"Preflight temporary path already exists: {temporary}")
        try:
            temporary.write_text(
                json.dumps(report, indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
            )
            os.replace(temporary, output)
        finally:
            if temporary.exists() and not temporary.is_symlink():
                temporary.unlink()
    return report


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", required=True)
    parser.add_argument("--contract", required=True)
    parser.add_argument("--staged-inputs")
    parser.add_argument("--output-json")
    args = parser.parse_args(argv)
    run_preflight(
        project_root=Path(args.project_root).resolve(),
        contract_path=Path(args.contract).resolve(),
        staged_inputs=(Path(args.staged_inputs).resolve() if args.staged_inputs else None),
        output=(Path(args.output_json).resolve() if args.output_json else None),
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
