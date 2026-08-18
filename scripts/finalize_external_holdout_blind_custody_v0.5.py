#!/usr/bin/env python3
"""Seal generic blind outputs and negative-access evidence before reveal."""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
from typing import Any, Iterable

from ploidypatch.artifact_manifest import sha256_file, verify_sha256sums
from ploidypatch.holdout_contract import load_holdout_contract, safe_relative_path


SCHEMA_VERSION = "ploidypatch.external_holdout_blind_custody.v0.5"
FRAMEWORK_VERSION = "PloidyPatch_external_validation_framework_v0.5"
EXECUTION_SCHEMA = "ploidypatch.external_holdout_execution_freeze.v0.5"
PROTOCOL_SCHEMA = "ploidypatch.external_holdout_protocol_freeze.v0.5"
ROLE_SCHEMA = "ploidypatch.blind_role_manifest.v0.5"
MOUNT_SCHEMA = "ploidypatch.blind_mount_manifest.v0.5"
NAMESPACE_SCHEMA = "ploidypatch.blind_namespace_validation.v0.5"
SCORE_SCHEMA = "ploidypatch.conflict_winner_guard_scores.v1"
POOL_SCHEMA = "ploidypatch.method_candidate_pool.v2"


def _reject_duplicate_keys(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    value: dict[str, Any] = {}
    for key, item in pairs:
        if key in value:
            raise ValueError(f"Duplicate JSON key: {key}")
        value[key] = item
    return value


def load_json(path: Path) -> dict[str, Any]:
    if not path.is_file() or path.is_symlink() or path.stat().st_size == 0:
        raise ValueError(f"Missing, empty or symlinked JSON input: {path}")
    value = json.loads(
        path.read_text(encoding="utf-8"), object_pairs_hook=_reject_duplicate_keys
    )
    if not isinstance(value, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return value


def safe_join(root: Path, raw: str, context: str) -> Path:
    relative = safe_relative_path(raw, context)
    path = root.joinpath(*relative.parts)
    try:
        path.resolve().relative_to(root.resolve())
    except ValueError:
        raise ValueError(f"{context} escapes blind project root") from None
    return path


def failed_attempt_manifest_summary(path: Path) -> dict[str, Any]:
    fields, rows = read_tsv(path)
    if fields != ["relative_path", "bytes", "sha256"] or not rows:
        raise ValueError("Failed-attempt manifest fields differ or are empty")
    digest = hashlib.sha256()
    total_bytes = 0
    seen: set[str] = set()
    for row in rows:
        relative = safe_relative_path(
            row["relative_path"], "failed-attempt manifest path"
        ).as_posix()
        if relative in seen or not row["bytes"].isdigit():
            raise ValueError("Failed-attempt manifest path/bytes are malformed")
        seen.add(relative)
        item_sha = row["sha256"]
        if len(item_sha) != 64 or any(c not in "0123456789abcdef" for c in item_sha):
            raise ValueError("Failed-attempt manifest SHA-256 is malformed")
        byte_count = int(row["bytes"])
        total_bytes += byte_count
        digest.update(f"{relative}\0{byte_count}\0{item_sha}\n".encode())
    return {
        "tree_sha256": digest.hexdigest(),
        "files": len(rows),
        "bytes": total_bytes,
    }


def validate_execution_patch(
    execution_manifest: dict[str, Any],
    protocol_manifest: dict[str, Any],
    execution: Path,
    protocol: Path,
    model: Path,
) -> None:
    patch = execution_manifest.get("execution_patch")
    if patch is None:
        if (
            execution_manifest.get("freeze_stage")
            != "post_metadata_pre_pair_pre_candidate_pre_label"
            or execution_manifest.get("code_commit")
            != protocol_manifest.get("code_commit")
        ):
            raise ValueError("Base execution freeze stage or commit differs")
        return
    changed = patch.get("changed_files") if isinstance(patch, dict) else None
    failed = patch.get("failed_attempt") if isinstance(patch, dict) else None
    failed_manifest_path = execution / "superseded_failed_attempt_manifest.tsv"
    failed_summary = (
        failed_attempt_manifest_summary(failed_manifest_path)
        if failed_manifest_path.is_file()
        else None
    )
    if (
        not isinstance(patch, dict)
        or patch.get("schema_version")
        != "ploidypatch.external_holdout_execution_patch.v0.5"
        or patch.get("freeze_stage")
        != "post_evaluator_truth_failed_blind_pre_candidate_pre_score_pre_label_execution_patch"
        or execution_manifest.get("freeze_stage") != patch.get("freeze_stage")
        or patch.get("base_code_commit") != protocol_manifest.get("code_commit")
        or patch.get("patch_code_commit") != execution_manifest.get("code_commit")
        or patch.get("base_protocol_SHA256SUMS_sha256")
        != sha256_file(protocol / "SHA256SUMS")
        or not isinstance(
            patch.get("superseded_execution_SHA256SUMS_sha256"), str
        )
        or len(patch["superseded_execution_SHA256SUMS_sha256"]) != 64
        or patch.get("contract_sha256") != sha256_file(protocol / "contract.json")
        or patch.get("composite_model_SHA256SUMS_sha256")
        != sha256_file(model / "SHA256SUMS")
        or patch.get("scientific_protocol_changed") is not False
        or patch.get("contract_or_policy_changed") is not False
        or patch.get("model_or_threshold_changed") is not False
        or patch.get("staged_inputs_changed") is not False
        or patch.get("truth_or_benchmark_regenerated") is not False
        or patch.get("evaluator_truth_construction_completed_before_patch") is not True
        or patch.get("blind_candidate_wgd_completed_before_patch") is not False
        or patch.get("candidate_generation_completed_before_patch") is not False
        or patch.get("formal_scores_generated_before_patch") is not False
        or patch.get("truth_labels_accessed_before_patch") is not False
        or patch.get("automatic_approval") is not False
        or execution_manifest.get("created_before")
        != {
            "wgd_pair_enumeration": False,
            "candidate_generation": True,
            "candidate_labels": True,
            "candidate_scores": True,
        }
        or not isinstance(changed, list)
        or not changed
        or not all(
            isinstance(row, dict)
            and row.get("status") in {"A", "M"}
            and isinstance(row.get("relative_path"), str)
            and isinstance(row.get("patch_sha256"), str)
            and len(row["patch_sha256"]) == 64
            for row in changed
        )
        or not isinstance(failed, dict)
        or not isinstance(failed.get("exit_status"), int)
        or failed["exit_status"] == 0
        or not isinstance(failed.get("tree_sha256"), str)
        or len(failed["tree_sha256"]) != 64
        or failed_summary is None
        or failed_summary.get("tree_sha256") != failed.get("tree_sha256")
        or failed_summary.get("files") != failed.get("files")
        or failed_summary.get("bytes") != failed.get("bytes")
        or not (execution / "patch_reason.md").is_file()
    ):
        raise ValueError("Execution patch provenance or scientific firewall differs")


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or len(reader.fieldnames) != len(set(reader.fieldnames)):
            raise ValueError(f"Malformed TSV header: {path}")
        return list(reader.fieldnames), list(reader)


def accepted_decision_digests(path: Path) -> set[str]:
    fields, rows = read_tsv(path)
    if not {"consensus_digest", "status"} <= set(fields):
        raise ValueError("Pool decisions lack consensus_digest/status")
    accepted = [row["consensus_digest"] for row in rows if row["status"] == "accepted"]
    if not accepted or len(accepted) != len(set(accepted)):
        raise ValueError("Accepted pool decisions are empty or duplicate")
    return set(accepted)


def score_digests(path: Path) -> set[str]:
    fields, rows = read_tsv(path)
    if "candidate_digest" not in fields:
        raise ValueError("Score table lacks candidate_digest")
    forbidden = [field for field in fields if "truth" in field.lower() or "label" in field.lower()]
    if forbidden:
        raise ValueError(f"Blind scores contain truth/label fields: {forbidden}")
    digests = [row["candidate_digest"] for row in rows]
    if not digests or len(digests) != len(set(digests)):
        raise ValueError("Blind score candidates are empty or duplicate")
    return set(digests)


def tree_digest(root: Path) -> tuple[str, int, int]:
    digest = hashlib.sha256()
    count = 0
    byte_count = 0
    for path in sorted(root.rglob("*")):
        if path.is_symlink():
            raise ValueError(f"Symlink is forbidden in blind output: {path}")
        if not path.is_file():
            continue
        relative = path.relative_to(root).as_posix()
        item_digest = sha256_file(path)
        size = path.stat().st_size
        digest.update(f"{relative}\0{size}\0{item_digest}\n".encode("utf-8"))
        count += 1
        byte_count += size
    if not count:
        raise ValueError("Blind project output is empty")
    return digest.hexdigest(), count, byte_count


def validate_mount_manifest(
    path: Path, holdout_id: str, role_root: Path, execution: Path,
    protocol: Path, model: Path, output_root: Path
) -> dict[str, Any]:
    manifest = load_json(path)
    execution_manifest = load_json(execution / "execution_manifest.json")
    environment_prefixes = {
        row["name"]: Path(row["host_prefix"]).resolve()
        for row in execution_manifest.get("environments", [])
        if isinstance(row, dict)
        and isinstance(row.get("name"), str)
        and isinstance(row.get("host_prefix"), str)
    }
    if len(environment_prefixes) != 7:
        raise ValueError("Execution freeze does not bind seven environments")
    mounts = manifest.get("mounts")
    if (
        manifest.get("schema_version") != MOUNT_SCHEMA
        or manifest.get("holdout_id") != holdout_id
        or not isinstance(mounts, list)
        or not mounts
    ):
        raise ValueError("Blind mount manifest schema or holdout differs")
    data_roles: set[str] = set()
    source_namespaces: set[str] = set()
    environment_namespaces: dict[str, set[str]] = {
        name: set() for name in environment_prefixes
    }
    metadata_namespaces: set[str] = set()
    singleton_roles: set[str] = set()
    seen_mounts: set[tuple[str, str, str]] = set()
    for item in mounts:
        if not isinstance(item, dict) or set(item) != {
            "role", "host_path", "namespace_path", "read_only"
        }:
            raise ValueError("Malformed blind mount entry")
        role = item["role"]
        host = Path(str(item["host_path"])).resolve()
        namespace = str(item["namespace_path"])
        identity = (str(role), str(host), namespace)
        if identity in seen_mounts:
            raise ValueError(f"Duplicate blind mount entry: {identity}")
        seen_mounts.add(identity)
        if (
            not isinstance(role, str)
            or not isinstance(namespace, str)
            or not namespace.startswith("/")
            or host == Path("/nas_data")
            or str(host).startswith("/nas_data/")
            or {part.lower() for part in host.parts}
            & {"evaluator_only", "target_complete", "truth", "labels"}
        ):
            raise ValueError(f"Unsafe blind mount entry: {item}")
        if role in {"shared_target", "candidate_only", "blind_benchmark"}:
            data_roles.add(role)
            expected = role_root / role
            if (
                host != expected.resolve()
                or namespace != f"/holdout/{role}"
                or item["read_only"] is not True
            ):
                raise ValueError(f"Blind data mount differs from sealed role: {role}")
        elif role == "frozen_source":
            source_namespaces.add(namespace)
            if host != (execution / "source").resolve() or item["read_only"] is not True:
                raise ValueError("Frozen source mount differs")
        elif role == "frozen_execution":
            singleton_roles.add(role)
            if (
                host != execution.resolve()
                or namespace != "/frozen/execution"
                or item["read_only"] is not True
            ):
                raise ValueError("Frozen execution mount differs")
        elif role == "frozen_protocol":
            singleton_roles.add(role)
            if (
                host != protocol.resolve()
                or namespace != "/frozen/protocol"
                or item["read_only"] is not True
            ):
                raise ValueError("Frozen protocol mount differs")
        elif role == "frozen_model":
            singleton_roles.add(role)
            if (
                host != model.resolve()
                or namespace != "/frozen/model"
                or item["read_only"] is not True
            ):
                raise ValueError("Frozen model mount differs")
        elif role == "blind_output":
            singleton_roles.add(role)
            if (
                host != output_root.resolve()
                or namespace != "/run/blind-run"
                or item["read_only"] is not False
            ):
                raise ValueError("Blind output must be the only writable project mount")
        elif role.startswith("frozen_environment:"):
            name = role.split(":", 1)[1]
            if name not in environment_prefixes or host != environment_prefixes[name]:
                raise ValueError(f"Frozen environment host differs: {role}")
            environment_namespaces[name].add(namespace)
            if item["read_only"] is not True:
                raise ValueError(f"Frozen environment mount is writable: {role}")
        elif role == "system_role_metadata":
            metadata_namespaces.add(namespace)
            expected_metadata = {
                "/holdout/blind_role_manifest.json": (
                    role_root / "role_manifest.json"
                ).resolve(),
                "/holdout/role_manifest.tsv": (
                    protocol / "role_manifest.tsv"
                ).resolve(),
                "/holdout/role_contract.json": (
                    protocol / "role_contract.json"
                ).resolve(),
            }
            if (
                expected_metadata.get(namespace) != host
                or item["read_only"] is not True
            ):
                raise ValueError("Role metadata mount differs or is writable")
        elif role in {"system_usr", "system_etc", "system_bin", "system_lib", "system_lib64"}:
            expected_namespace = {
                "system_usr": "/usr",
                "system_etc": "/etc",
                "system_bin": "/bin",
                "system_lib": "/lib",
                "system_lib64": "/lib64",
            }[role]
            if (
                str(host) != str(Path(expected_namespace).resolve())
                or namespace != expected_namespace
                or item["read_only"] is not True
            ):
                raise ValueError(f"System mount differs or is writable: {role}")
        elif role == "system_conda":
            if (
                str(host) != namespace
                or not (host / "conda-meta/history").is_file()
                or not (host / "condabin/conda").is_file()
                or item["read_only"] is not True
            ):
                raise ValueError("System conda mount differs or is writable")
        else:
            raise ValueError(f"Unexpected blind mount role: {role}")
    if data_roles != {"shared_target", "candidate_only", "blind_benchmark"}:
        raise ValueError("Blind namespace must mount exactly three data roles")
    if source_namespaces != {
        "/frozen/source",
        "/run/blind-run/project/code",
    }:
        raise ValueError("Frozen source aliases differ from the isolated project contract")
    for name, namespaces in environment_namespaces.items():
        expected = {
            str(environment_prefixes[name]),
            f"/frozen/envs/{name}",
            f"/run/blind-run/project/envs/{name}",
        }
        if namespaces != expected:
            raise ValueError(f"Frozen environment aliases differ: {name}")
    if metadata_namespaces != {
        "/holdout/blind_role_manifest.json",
        "/holdout/role_manifest.tsv",
        "/holdout/role_contract.json",
    }:
        raise ValueError("Role metadata namespace mounts differ")
    if singleton_roles != {
        "frozen_execution", "frozen_protocol", "frozen_model", "blind_output"
    }:
        raise ValueError("A singleton frozen/output mount is missing")
    return manifest


def reject_forbidden_command_content(path: Path) -> None:
    text = path.read_text(encoding="utf-8", errors="strict").lower()
    if not text.strip():
        raise ValueError(f"Empty blind command record: {path}")
    forbidden = ("/nas_data", "evaluator_only", "target_complete", "/truth", "/labels")
    if any(token in text for token in forbidden):
        raise ValueError(f"Blind command record names forbidden data: {path}")


def finalize_custody(
    *,
    blind_project: Path,
    blind_output_root: Path,
    execution: Path,
    protocol: Path,
    model: Path,
    role_root: Path,
    runner_command: Path,
    mount_manifest_path: Path,
    namespace_validation_path: Path,
    bwrap_version: str,
    output: Path,
) -> Path:
    for root in (execution, protocol, model, role_root):
        verify_sha256sums(root, ignore_checksum_file=True)
    execution_manifest = load_json(execution / "execution_manifest.json")
    protocol_manifest = load_json(protocol / "protocol_manifest.json")
    contract = load_holdout_contract(protocol / "contract.json")
    if (
        execution_manifest.get("schema_version") != EXECUTION_SCHEMA
        or protocol_manifest.get("schema_version") != PROTOCOL_SCHEMA
        or execution_manifest.get("holdout_id") != contract.holdout_id
        or protocol_manifest.get("holdout_id") != contract.holdout_id
        or execution_manifest.get("contract_sha256") != sha256_file(protocol / "contract.json")
        or execution_manifest.get("protocol_SHA256SUMS_sha256")
        != sha256_file(protocol / "SHA256SUMS")
        or execution_manifest.get("composite_model_SHA256SUMS_sha256")
        != sha256_file(model / "SHA256SUMS")
        or any(
            execution_manifest.get(field) is not False
            for field in (
                "network_access_in_blind_runner",
                "nas_data_mount_in_blind_runner",
                "complete_target_annotation_mount_in_blind_runner",
                "evaluator_only_mount_in_blind_runner",
                "truth_or_label_mount_in_blind_runner",
            )
        )
    ):
        raise ValueError("Execution/protocol/model holdout bindings differ")
    validate_execution_patch(
        execution_manifest, protocol_manifest, execution, protocol, model
    )

    role_manifest = load_json(role_root / "role_manifest.json")
    if (
        role_manifest.get("schema_version") != ROLE_SCHEMA
        or role_manifest.get("holdout_id") != contract.holdout_id
        or role_manifest.get("contract_sha256") != sha256_file(protocol / "contract.json")
        or role_manifest.get("truth_access") is not False
        or role_manifest.get("complete_target_annotation_present") is not False
        or role_manifest.get("evaluator_references_present") is not False
        or role_manifest.get("network_access") is not False
        or role_manifest.get("roles")
        != ["shared_target", "candidate_only", "blind_benchmark"]
        or role_manifest.get("protocol_SHA256SUMS_sha256")
        != sha256_file(protocol / "SHA256SUMS")
        or role_manifest.get("blind_benchmark_SHA256SUMS_sha256")
        != sha256_file(role_root / "blind_benchmark/SHA256SUMS")
        or role_manifest.get("blind_benchmark_manifest_sha256")
        != sha256_file(role_root / "blind_benchmark/blind_manifest.json")
    ):
        raise ValueError("Blind role manifest violates the role firewall")
    if any(
        not (role_root / role).is_dir()
        for role in ("shared_target", "candidate_only", "blind_benchmark")
    ):
        raise ValueError("Blind role root lacks a required truth-free data role")
    verify_sha256sums(role_root / "blind_benchmark", ignore_checksum_file=True)
    benchmark_manifest = load_json(role_root / "blind_benchmark/blind_manifest.json")
    if (
        benchmark_manifest.get("schema_version")
        != "ploidypatch.blind_benchmark_input.v0.5"
        or benchmark_manifest.get("truth_access") is not False
        or benchmark_manifest.get("complete_target_annotation_access") is not False
        or benchmark_manifest.get("perturbed_annotation", {}).get("file_name")
        != "perturbed.gff3"
        or benchmark_manifest.get("perturbed_annotation", {}).get("sha256")
        != sha256_file(role_root / "blind_benchmark/perturbed.gff3")
    ):
        raise ValueError("Blind benchmark role is not a sealed truth-free perturbation")

    mount_manifest = validate_mount_manifest(
        mount_manifest_path,
        contract.holdout_id,
        role_root,
        execution,
        protocol,
        model,
        blind_output_root,
    )
    namespace = load_json(namespace_validation_path)
    if (
        namespace.get("schema_version") != NAMESPACE_SCHEMA
        or namespace.get("holdout_id") != contract.holdout_id
        or namespace.get("mount_manifest_sha256") != sha256_file(mount_manifest_path)
        or namespace.get("host_role_manifest_sha256")
        != sha256_file(role_root / "role_manifest.json")
        or namespace.get("shared_target_visible") is not True
        or namespace.get("candidate_only_visible") is not True
        or namespace.get("blind_benchmark_visible") is not True
        or namespace.get("blind_benchmark_manifest_sha256")
        != sha256_file(role_root / "blind_benchmark/blind_manifest.json")
        or namespace.get("evaluator_only_visible") is not False
        or namespace.get("truth_visible") is not False
        or namespace.get("complete_target_annotation_visible") is not False
        or namespace.get("nas_data_visible") is not False
    ):
        raise ValueError("Namespace role validation is incomplete or violated")

    outputs = execution_manifest.get("blind_outputs")
    if not isinstance(outputs, dict) or set(outputs) != {
        "scores", "score_manifest", "pool_decisions", "pool_manifest", "command_log"
    }:
        raise ValueError("Execution freeze lacks exact blind output paths")
    paths = {
        name: safe_join(blind_project, relative, f"blind output {name}")
        for name, relative in outputs.items()
    }
    for name, path in paths.items():
        if not path.is_file() or path.is_symlink() or path.stat().st_size == 0:
            raise ValueError(f"Missing, empty or symlinked blind output {name}: {path}")
    pool_manifest = load_json(paths["pool_manifest"])
    score_manifest = load_json(paths["score_manifest"])
    if (
        pool_manifest.get("schema_version") != POOL_SCHEMA
        or pool_manifest.get("outputs", {}).get("decisions", {}).get("sha256")
        != sha256_file(paths["pool_decisions"])
    ):
        raise ValueError("Blind pool manifest does not bind decisions")
    if (
        score_manifest.get("schema_version") != SCORE_SCHEMA
        or score_manifest.get("truth_access") is not False
        or score_manifest.get("inputs", {}).get("pool_decisions")
        != sha256_file(paths["pool_decisions"])
        or score_manifest.get("inputs", {}).get("pool_manifest")
        != sha256_file(paths["pool_manifest"])
        or score_manifest.get("outputs", {}).get("scores", {}).get("sha256")
        != sha256_file(paths["scores"])
    ):
        raise ValueError("Blind score manifest does not bind pool and scores")
    if score_digests(paths["scores"]) != accepted_decision_digests(paths["pool_decisions"]):
        raise ValueError("Blind score and accepted-pool candidate universes differ")
    for path in (runner_command, paths["command_log"]):
        if not path.is_file() or path.is_symlink():
            raise ValueError(f"Missing blind command record: {path}")
        reject_forbidden_command_content(path)
    runner_text = runner_command.read_text(encoding="utf-8")
    for required in (
        "--unshare-all",
        "--unshare-net",
        "--clearenv",
        "PLOIDYPATCH_NETWORK_ACCESS",
        "none",
    ):
        if required not in runner_text:
            raise ValueError(f"Blind runner command lacks {required}")

    tree_hash, file_count, byte_count = tree_digest(blind_project)
    custody = {
        "schema_version": SCHEMA_VERSION,
        "framework_version": FRAMEWORK_VERSION,
        "holdout_id": contract.holdout_id,
        "policy_id": contract.policy_id,
        "model_version": contract.model_version,
        "runner_identity": f"bubblewrap_external_holdout_v0.5:{contract.holdout_id}",
        "frozen_before_truth_reveal_at": datetime.now(timezone.utc)
        .replace(microsecond=0)
        .isoformat()
        .replace("+00:00", "Z"),
        "bubblewrap": {
            "version": bwrap_version,
            "required_flags": ["--unshare-all", "--unshare-net", "--clearenv"],
            "command_sha256": sha256_file(runner_command),
        },
        "truth_mounted": False,
        "complete_target_annotation_mounted": False,
        "evaluator_references_mounted": False,
        "nas_data_mounted": False,
        "network_access": False,
        "mount_manifest": {
            "sha256": sha256_file(mount_manifest_path),
            "mounts": mount_manifest["mounts"],
        },
        "namespace_validation_sha256": sha256_file(namespace_validation_path),
        "host_role_manifest_sha256": sha256_file(role_root / "role_manifest.json"),
        "commands": {
            "runner_relative_path": runner_command.name,
            "runner_sha256": sha256_file(runner_command),
            "pipeline_relative_path": paths["command_log"].relative_to(blind_project).as_posix(),
            "pipeline_sha256": sha256_file(paths["command_log"]),
        },
        "frozen_inputs": {
            "execution_SHA256SUMS_sha256": sha256_file(execution / "SHA256SUMS"),
            "protocol_SHA256SUMS_sha256": sha256_file(protocol / "SHA256SUMS"),
            "composite_model_SHA256SUMS_sha256": sha256_file(model / "SHA256SUMS"),
            "blind_role_SHA256SUMS_sha256": sha256_file(role_root / "SHA256SUMS"),
        },
        "execution_patch": (
            {
                "active": True,
                "freeze_stage": execution_manifest["execution_patch"]["freeze_stage"],
                "base_code_commit": execution_manifest["execution_patch"]["base_code_commit"],
                "patch_code_commit": execution_manifest["execution_patch"]["patch_code_commit"],
                "superseded_execution_SHA256SUMS_sha256": execution_manifest[
                    "execution_patch"
                ]["superseded_execution_SHA256SUMS_sha256"],
                "failed_attempt_tree_sha256": execution_manifest["execution_patch"][
                    "failed_attempt"
                ]["tree_sha256"],
            }
            if execution_manifest.get("execution_patch") is not None
            else {"active": False}
        ),
        "blind_outputs": {
            name: {
                "relative_path": path.relative_to(blind_project).as_posix(),
                "sha256": sha256_file(path),
            }
            for name, path in paths.items()
        }
        | {"output_tree_sha256": tree_hash, "file_count": file_count, "bytes": byte_count},
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    descriptor = os.open(output, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o440)
    with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as handle:
        json.dump(custody, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return output


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--blind-project-root", required=True)
    parser.add_argument("--blind-output-root", required=True)
    parser.add_argument("--execution-freeze", required=True)
    parser.add_argument("--protocol-freeze", required=True)
    parser.add_argument("--composite-model-freeze", required=True)
    parser.add_argument("--blind-role-root", required=True)
    parser.add_argument("--runner-command", required=True)
    parser.add_argument("--mount-manifest", required=True)
    parser.add_argument("--namespace-validation", required=True)
    parser.add_argument("--bwrap-version", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args(argv)
    finalize_custody(
        blind_project=Path(args.blind_project_root).resolve(),
        blind_output_root=Path(args.blind_output_root).resolve(),
        execution=Path(args.execution_freeze).resolve(),
        protocol=Path(args.protocol_freeze).resolve(),
        model=Path(args.composite_model_freeze).resolve(),
        role_root=Path(args.blind_role_root).resolve(),
        runner_command=Path(args.runner_command).resolve(),
        mount_manifest_path=Path(args.mount_manifest).resolve(),
        namespace_validation_path=Path(args.namespace_validation).resolve(),
        bwrap_version=args.bwrap_version,
        output=Path(args.output).resolve(),
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
