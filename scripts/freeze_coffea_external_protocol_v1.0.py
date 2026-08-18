#!/usr/bin/env python3
"""Freeze Coffea H1 science, normalized inputs and exact source archive."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import shutil
import tarfile
from typing import Any, Iterable

from ploidypatch.artifact_manifest import (
    sha256_file,
    verify_sha256sums,
    write_sha256sums,
)
from ploidypatch.holdout_contract import (
    KNOWN_SUBGENOME_CORE_H1_MODEL_VERSION,
    KnownSubgenomeCoreH1ScientificParameters,
    load_holdout_contract,
)


SCHEMA_VERSION = "ploidypatch.coffea_core_h1_protocol_freeze.v1.0"
HOLDOUT_ID = "coffea_et39_v1.0"
POLICY_ID = "ploidypatch_coffea_external_core_h1_v1.0"
FULL_SHA = re.compile(r"[0-9a-f]{40}")
PROTOCOL_FILES = (
    "config/holdouts/coffea_et39_v1.0/contract.json",
    "config/holdouts/coffea_et39_v1.0/no_ranker/composite_manifest.json",
    "config/holdouts/coffea_et39_v1.0/no_ranker/SHA256SUMS",
    "config/coffea_external_input_sources_v1.0.tsv",
    "config/coffea_external_validation_policy_v1.0.tsv",
    "config/coffea_external_event_definition_v1.0.tsv",
    "config/coffea_et39_homoeolog_groups_v1.0.tsv",
    "config/coffea_gardenia_seqid_aliases_v1.0.tsv",
    "config/primary_seqids/coffea_arabica_et39_hifi.tsv",
    "config/primary_seqids/coffea_eugenioides_bua_v1.tsv",
    "config/primary_seqids/coffea_mauritiana_v1.tsv",
    "config/primary_seqids/gardenia_jasminoides_asm1310374v1.tsv",
    "config/primary_seqids/ophiorrhiza_pumila_v1.tsv",
    "docs/COFFEA_EXTERNAL_SELECTION_AND_CONTAMINATION_RATIONALE_v1.0.md",
    "docs/COFFEA_EXTERNAL_CORE_H1_PROTOCOL_v1.0.md",
)
IMPLEMENTATION_FILES = (
    "src/ploidypatch/artifact_manifest.py",
    "src/ploidypatch/audit.py",
    "src/ploidypatch/gff.py",
    "src/ploidypatch/holdout_contract.py",
    "src/ploidypatch/io.py",
    "src/ploidypatch/known_subgenome_h1.py",
    "src/ploidypatch/normalize.py",
    "src/ploidypatch/perturb.py",
    "src/ploidypatch/protein_universe.py",
    "src/ploidypatch/safe_tar_fasta.py",
    "src/ploidypatch/seqid_alias.py",
    "src/ploidypatch/synteny_io.py",
    "scripts/stage_external_holdout_inputs_v0.5.py",
    "scripts/preflight_coffea_external_inputs_v1.0.py",
    "scripts/prepare_coffea_external_normalized_inputs_v1.0.py",
    "scripts/build_coffea_protein_supported_universes_v1.0.py",
    "scripts/freeze_coffea_external_protocol_v1.0.py",
)
FORBIDDEN_TARGET_ARTIFACTS = (
    "data/derived/external_evaluator/coffea_et39_v1.0",
    "benchmark/structure/coffea_et39_v1.0",
    "results/evaluator/coffea_et39_v1.0",
    "results/copy_collapse/external/coffea_et39_v1.0",
    "results/blind_runs/coffea_et39_v1.0",
    "results/holdouts/coffea_et39_v1.0",
    "work/coffea_et39_v1.0",
)


def _read_json(path: Path) -> dict[str, Any]:
    if not path.is_file() or path.is_symlink() or path.stat().st_size == 0:
        raise ValueError(f"Missing, empty or symlinked JSON: {path}")
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return value


def _write_tsv(
    path: Path, header: tuple[str, ...], rows: Iterable[Iterable[Any]]
) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def _safe_tar_name(name: str) -> PurePosixPath:
    path = PurePosixPath(name)
    if (
        not name
        or name.startswith("/")
        or "\\" in name
        or path.is_absolute()
        or any(part in {"", ".", ".."} for part in path.parts)
    ):
        raise ValueError(f"Unsafe frozen source archive path: {name!r}")
    return path


def verify_source_freeze(source_freeze: Path, code_commit: str) -> Path:
    """Verify archive bytes and every extracted regular file without Git."""

    if FULL_SHA.fullmatch(code_commit) is None:
        raise ValueError("code_commit must be one full lowercase Git SHA")
    archive = source_freeze / "source.tar"
    checksum = source_freeze / "source.tar.sha256"
    commit_file = source_freeze / "source_commit.txt"
    code_root = source_freeze / "source"
    for path in (archive, checksum, commit_file):
        if not path.is_file() or path.is_symlink() or path.stat().st_size == 0:
            raise ValueError(f"Frozen source metadata is missing: {path}")
    if commit_file.read_text(encoding="utf-8") != code_commit + "\n":
        raise ValueError("Frozen source commit differs")
    expected_line = f"{sha256_file(archive)}  source.tar\n"
    if checksum.read_text(encoding="utf-8") != expected_line:
        raise ValueError("Frozen source archive checksum differs")
    if not code_root.is_dir() or code_root.is_symlink():
        raise ValueError("Frozen extracted source root is missing or symlinked")

    archived: dict[str, str] = {}
    with tarfile.open(archive, mode="r:") as handle:
        for member in handle.getmembers():
            path = _safe_tar_name(member.name.rstrip("/") if member.isdir() else member.name)
            if member.isdir():
                continue
            if not member.isfile():
                raise ValueError(f"Non-regular frozen source member: {member.name}")
            stream = handle.extractfile(member)
            if stream is None:
                raise ValueError(f"Unreadable frozen source member: {member.name}")
            digest = hashlib.sha256()
            while chunk := stream.read(8 * 1024 * 1024):
                digest.update(chunk)
            relative = path.as_posix()
            if relative in archived:
                raise ValueError(f"Duplicate frozen source member: {relative}")
            archived[relative] = digest.hexdigest()
    extracted = {
        path.relative_to(code_root).as_posix(): path
        for path in code_root.rglob("*")
        if path.is_file()
    }
    if set(extracted) != set(archived):
        raise ValueError("Extracted source file universe differs from Git archive")
    for relative, expected in archived.items():
        path = extracted[relative]
        if path.is_symlink() or sha256_file(path) != expected:
            raise ValueError(f"Extracted source differs from archive: {relative}")
    return code_root


def _verify_stage(stage: Path, contract_sha256: str, code_commit: str) -> dict[str, Any]:
    verify_sha256sums(stage, ignore_checksum_file=True)
    manifest = _read_json(stage / "role_contract.json")
    if (
        manifest.get("schema_version")
        != "ploidypatch.external_holdout_input_stage.v0.5"
        or manifest.get("holdout_id") != HOLDOUT_ID
        or manifest.get("policy_id") != POLICY_ID
        or manifest.get("code_commit") != code_commit
        or manifest.get("contract", {}).get("sha256") != contract_sha256
    ):
        raise ValueError("Coffea staged inputs differ from the frozen code/contract")
    return manifest


def _verify_derived_inputs(normalized: Path, universes: Path) -> tuple[dict[str, Any], dict[str, Any]]:
    verify_sha256sums(normalized, ignore_checksum_file=True)
    verify_sha256sums(universes, ignore_checksum_file=True)
    normalized_manifest = _read_json(normalized / "manifest.json")
    universe_manifest = _read_json(universes / "manifest.json")
    if (
        normalized_manifest.get("schema_version")
        != "ploidypatch.coffea_external_normalized_inputs.v1.0"
        or normalized_manifest.get("holdout_id") != HOLDOUT_ID
        or any(normalized_manifest.get("firewall", {}).values())
    ):
        raise ValueError("Coffea normalized input firewall differs")
    expected_species = {
        "Coffea_arabica_ET39",
        "Coffea_eugenioides_BuA",
        "Coffea_mauritiana",
        "Gardenia_jasminoides",
        "Ophiorrhiza_pumila",
    }
    if (
        universe_manifest.get("schema_version")
        != "ploidypatch.coffea_protein_universes.v1.0"
        or universe_manifest.get("holdout_id") != HOLDOUT_ID
        or universe_manifest.get("policy_id") != POLICY_ID
        or universe_manifest.get("labels_used") is not False
        or universe_manifest.get("truth_pairs_used") is not False
        or universe_manifest.get("candidate_predictions_used") is not False
        or set(universe_manifest.get("species", {})) != expected_species
    ):
        raise ValueError("Coffea protein-supported universes differ")
    for species_id in expected_species:
        species_manifest = _read_json(universes / species_id / "manifest.json")
        counts = species_manifest.get("counts", {})
        if (
            species_manifest.get("truth_access") is not False
            or species_manifest.get("candidate_access") is not False
            or species_manifest.get("fuzzy_mapping_used") is not False
            or counts.get("genes_with_exact_provider_protein", 0) < 1
            or counts.get("unresolved_CDS_rows_with_Parent") != 0
        ):
            raise ValueError(f"Unsafe Coffea protein universe: {species_id}")
    return normalized_manifest, universe_manifest


def _chmod_read_only(root: Path) -> None:
    if os.name == "nt":
        return
    for path in sorted(root.rglob("*"), reverse=True):
        path.chmod(0o550 if path.is_dir() else 0o440)
    root.chmod(0o550)


def freeze_protocol(
    *,
    project_root: Path,
    source_freeze: Path,
    staged_inputs: Path,
    normalized_inputs: Path,
    protein_universes: Path,
    code_commit: str,
    output: Path,
) -> Path:
    partial = Path(str(output) + ".partial")
    if output.exists() or output.is_symlink() or partial.exists() or partial.is_symlink():
        raise FileExistsError("Refusing to overwrite Coffea protocol freeze or partial")
    code_root = verify_source_freeze(source_freeze, code_commit)
    contract_path = code_root / PROTOCOL_FILES[0]
    contract = load_holdout_contract(contract_path)
    if (
        contract.holdout_id != HOLDOUT_ID
        or contract.policy_id != POLICY_ID
        or contract.model_version != KNOWN_SUBGENOME_CORE_H1_MODEL_VERSION
        or not isinstance(
            contract.scientific_parameters, KnownSubgenomeCoreH1ScientificParameters
        )
        or contract.scientific_parameters.h2_or_topology_ranking != "forbidden"
        or contract.scientific_parameters.all_arm_collateral_loss_maximum != 0
    ):
        raise ValueError("Freeze requires exact Coffea known-subgenome H1 contract")
    for relative in (*PROTOCOL_FILES, *IMPLEMENTATION_FILES):
        path = code_root / relative
        if not path.is_file() or path.is_symlink() or path.stat().st_size == 0:
            raise ValueError(f"Missing frozen Coffea implementation: {relative}")
    for relative in FORBIDDEN_TARGET_ARTIFACTS:
        path = project_root / relative
        if path.exists() or path.is_symlink():
            raise ValueError(f"Coffea protocol freeze is too late: {relative}")
    stage_manifest = _verify_stage(
        staged_inputs, sha256_file(contract_path), code_commit
    )
    normalized_manifest, universe_manifest = _verify_derived_inputs(
        normalized_inputs, protein_universes
    )
    no_ranker = code_root / "config/holdouts/coffea_et39_v1.0/no_ranker"
    verify_sha256sums(no_ranker, ignore_checksum_file=True)

    partial.mkdir(parents=True)
    try:
        shutil.copyfile(contract_path, partial / "contract.json")
        shutil.copyfile(staged_inputs / "role_contract.json", partial / "role_contract.json")
        shutil.copyfile(staged_inputs / "role_manifest.tsv", partial / "role_manifest.tsv")
        shutil.copyfile(staged_inputs / "SHA256SUMS", partial / "staged_input_SHA256SUMS")
        shutil.copyfile(normalized_inputs / "manifest.json", partial / "normalized_manifest.json")
        shutil.copyfile(normalized_inputs / "SHA256SUMS", partial / "normalized_SHA256SUMS")
        shutil.copyfile(protein_universes / "manifest.json", partial / "protein_universes_manifest.json")
        shutil.copyfile(protein_universes / "SHA256SUMS", partial / "protein_universes_SHA256SUMS")
        artifact_root = partial / "protocol_artifacts"
        for relative in PROTOCOL_FILES[1:]:
            destination = artifact_root / relative
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(code_root / relative, destination)

        frozen_files = sorted(set((*PROTOCOL_FILES, *IMPLEMENTATION_FILES)))
        _write_tsv(
            partial / "implementation_manifest.tsv",
            ("relative_path", "bytes", "sha256"),
            (
                (relative, (code_root / relative).stat().st_size, sha256_file(code_root / relative))
                for relative in frozen_files
            ),
        )
        _write_tsv(
            partial / "forbidden_target_artifacts.tsv",
            ("relative_path",),
            ((relative,) for relative in FORBIDDEN_TARGET_ARTIFACTS),
        )
        _write_tsv(
            partial / "run_contract.tsv",
            ("field", "value"),
            (
                ("schema_version", SCHEMA_VERSION),
                ("holdout_id", HOLDOUT_ID),
                ("policy_id", POLICY_ID),
                ("protocol_profile", "core_H1_known_subgenome_no_ranker"),
                ("code_commit", code_commit),
                ("freeze_stage", "post_metadata_compatibility_pre_pair_pre_candidate_pre_label"),
                ("wgd_pairs_enumerated_before_freeze", "false"),
                ("candidate_counts_computed_before_freeze", "false"),
                ("truth_labels_accessed_before_freeze", "false"),
                ("H2_or_topology_ranking", "forbidden"),
                ("bootstrap_replicates", "20000"),
                ("all_arm_collateral_loss_maximum", "0"),
            ),
        )
        manifest = {
            "schema_version": SCHEMA_VERSION,
            "holdout_id": HOLDOUT_ID,
            "policy_id": POLICY_ID,
            "test_role": contract.test_role,
            "protocol_profile": "core_H1_known_subgenome_no_ranker",
            "model_version": contract.model_version,
            "ranker_enabled": False,
            "h2_or_topology_ranking_enabled": False,
            "formal_runner_frozen": False,
            "formal_runner_status": "execution_freeze_required_before_any_pair_or_candidate",
            "code_commit": code_commit,
            "source_archive_sha256": sha256_file(source_freeze / "source.tar"),
            "freeze_stage": "post_metadata_compatibility_pre_pair_pre_candidate_pre_label",
            "truth_access": False,
            "wgd_pairs_enumerated": False,
            "candidate_counts_computed": False,
            "truth_labels_accessed": False,
            "contract_sha256": sha256_file(partial / "contract.json"),
            "staged_input_SHA256SUMS_sha256": sha256_file(staged_inputs / "SHA256SUMS"),
            "normalized_input_SHA256SUMS_sha256": sha256_file(normalized_inputs / "SHA256SUMS"),
            "protein_universes_SHA256SUMS_sha256": sha256_file(protein_universes / "SHA256SUMS"),
            "no_ranker_SHA256SUMS_sha256": sha256_file(no_ranker / "SHA256SUMS"),
            "implementation_file_count": len(frozen_files),
            "reference_roles": {"target": 1, "candidate_reference": 2, "evaluator_reference": 2},
            "stage_role_boundaries": stage_manifest["role_boundaries"],
            "normalized_firewall": normalized_manifest["firewall"],
            "protein_universe_species": sorted(universe_manifest["species"]),
            "forbidden_target_artifacts": list(FORBIDDEN_TARGET_ARTIFACTS),
        }
        (partial / "protocol_manifest.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        write_sha256sums(partial)
        verify_sha256sums(partial, ignore_checksum_file=True)
        _chmod_read_only(partial)
        os.replace(partial, output)
    except BaseException:
        if partial.exists() and not partial.is_symlink():
            shutil.rmtree(partial)
        raise
    return output


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", required=True)
    parser.add_argument("--source-freeze", required=True)
    parser.add_argument("--staged-inputs", required=True)
    parser.add_argument("--normalized-inputs", required=True)
    parser.add_argument("--protein-universes", required=True)
    parser.add_argument("--code-commit", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args(argv)
    freeze_protocol(
        project_root=Path(args.project_root).resolve(),
        source_freeze=Path(args.source_freeze).resolve(),
        staged_inputs=Path(args.staged_inputs).resolve(),
        normalized_inputs=Path(args.normalized_inputs).resolve(),
        protein_universes=Path(args.protein_universes).resolve(),
        code_commit=args.code_commit,
        output=Path(args.output_dir).resolve(),
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

