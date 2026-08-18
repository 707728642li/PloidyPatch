#!/usr/bin/env python3
"""Seal Walnut core-H1 blind outputs and negative-access evidence before reveal."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
from typing import Any, Iterable

from ploidypatch.artifact_manifest import sha256_file, verify_sha256sums, write_sha256sums
from ploidypatch.walnut_h1 import POOL_SCHEMAS
from ploidypatch.walnut_h1_framework import (
    BLIND_OUTPUTS,
    CUSTODY_SCHEMA,
    HOLDOUT_ID,
    MOUNT_SCHEMA,
    NAMESPACE_SCHEMA,
    POLICY_ID,
    RAW_PREDICTION_TREE_KEYS,
    REQUIRED_ENVIRONMENTS,
    ROLE_SCHEMA,
    load_json,
    read_tsv,
    reject_forbidden_text,
    safe_join,
    verify_execution,
)


RAW_SCHEMA = "ploidypatch.walnut_h1_raw_predictions.v0.8"


def tree_digest(root: Path) -> tuple[str, int, int]:
    digest = hashlib.sha256()
    count = 0
    byte_count = 0
    for path in sorted(root.rglob("*")):
        if path.is_symlink():
            raise ValueError(f"Symlink is forbidden in blind tree: {path}")
        if not path.is_file():
            continue
        relative = path.relative_to(root).as_posix()
        size = path.stat().st_size
        item_sha = sha256_file(path)
        digest.update(f"{item_sha}  {relative}\n".encode())
        count += 1; byte_count += size
    if count == 0:
        raise ValueError(f"Blind tree is empty: {root}")
    return digest.hexdigest(), count, byte_count


def validate_mount_manifest(
    path: Path, *, role_root: Path, execution: Path, protocol: Path,
    blind_output_root: Path
) -> dict[str, Any]:
    manifest = load_json(path)
    execution_manifest = load_json(execution / "execution_manifest.json")
    environments = {
        row["name"]: Path(row["host_prefix"]).resolve()
        for row in execution_manifest["environments"]
    }
    if set(environments) != REQUIRED_ENVIRONMENTS:
        raise ValueError("Execution environment bindings differ")
    mounts = manifest.get("mounts")
    if manifest.get("schema_version") != MOUNT_SCHEMA or not isinstance(mounts, list):
        raise ValueError("Mount manifest schema differs")
    required = {
        ("shared_target", str((role_root / "shared_target").resolve()),
         "/holdout/shared_target", True),
        ("candidate_only", str((role_root / "candidate_only").resolve()),
         "/holdout/candidate_only", True),
        ("blind_benchmark", str((role_root / "blind_benchmark").resolve()),
         "/holdout/blind_benchmark", True),
        ("frozen_execution", str(execution.resolve()), "/frozen/execution", True),
        ("frozen_protocol", str(protocol.resolve()), "/frozen/protocol", True),
        ("blind_output", str(blind_output_root.resolve()), "/run/blind-run", False),
    }
    observed: set[tuple[str, str, str, bool]] = set()
    data_roles: set[str] = set()
    source_namespaces: set[str] = set()
    environment_namespaces = {name: set() for name in environments}
    metadata_namespaces: set[str] = set()
    system_roles: set[str] = set()
    for item in mounts:
        if not isinstance(item, dict) or set(item) != {
            "role", "host_path", "namespace_path", "read_only"
        }:
            raise ValueError("Malformed mount entry")
        role = item["role"]
        host = Path(str(item["host_path"])).resolve()
        namespace = str(item["namespace_path"])
        read_only = item["read_only"]
        if (
            not isinstance(role, str) or not isinstance(read_only, bool)
            or not namespace.startswith("/")
            or host == Path("/nas_data") or str(host).startswith("/nas_data/")
            or set(part.casefold() for part in host.parts)
            & {"evaluator_only", "target_complete", "truth", "labels", "model", "ranker"}
        ):
            raise ValueError(f"Unsafe blind mount: {item}")
        record = (role, str(host), namespace, read_only)
        if record in observed:
            raise ValueError("Duplicate mount entry")
        observed.add(record)
        if role in {"shared_target", "candidate_only", "blind_benchmark"}:
            data_roles.add(role)
        elif role.startswith("frozen_environment:"):
            name = role.split(":", 1)[1]
            if name not in environments or host != environments[name] or read_only is not True:
                raise ValueError("Frozen environment mount differs")
            environment_namespaces[name].add(namespace)
        elif role == "frozen_source":
            if read_only is not True:
                raise ValueError("Frozen source mount is writable")
            if host != (execution / "source").resolve():
                raise ValueError("Frozen source host differs")
            source_namespaces.add(namespace)
        elif role == "system_role_metadata":
            expected_metadata = {
                "/holdout/blind_role_manifest.json": (role_root / "role_manifest.json").resolve(),
                "/holdout/role_manifest.tsv": (protocol / "role_manifest.tsv").resolve(),
                "/holdout/role_contract.json": (protocol / "role_contract.json").resolve(),
            }
            if expected_metadata.get(namespace) != host or read_only is not True:
                raise ValueError("Role metadata mount differs")
            metadata_namespaces.add(namespace)
        elif role in {
            "system_usr", "system_etc", "system_bin", "system_lib",
            "system_lib64", "system_conda",
        }:
            if read_only is not True:
                raise ValueError(f"System mount is writable: {role}")
            if role == "system_conda":
                if str(host) != namespace or not (host / "conda-meta/history").is_file():
                    raise ValueError("System conda mount differs")
            else:
                expected_namespace = {
                    "system_usr": "/usr", "system_etc": "/etc", "system_bin": "/bin",
                    "system_lib": "/lib", "system_lib64": "/lib64",
                }[role]
                if namespace != expected_namespace or host != Path(expected_namespace).resolve():
                    raise ValueError(f"System mount host differs: {role}")
            system_roles.add(role)
        elif record not in required:
            raise ValueError(f"Unexpected blind mount role: {role}")
    if not required <= observed or data_roles != {
        "shared_target", "candidate_only", "blind_benchmark"
    }:
        raise ValueError("Blind namespace lacks exact three roles/frozen lineage")
    if source_namespaces != {"/frozen/source", "/run/blind-run/project/code"}:
        raise ValueError("Frozen source namespace aliases differ")
    for name, namespaces in environment_namespaces.items():
        if namespaces != {
            str(environments[name]), f"/frozen/envs/{name}",
            f"/run/blind-run/project/envs/{name}",
        }:
            raise ValueError(f"Frozen environment aliases differ: {name}")
    if metadata_namespaces != {
        "/holdout/blind_role_manifest.json", "/holdout/role_manifest.tsv",
        "/holdout/role_contract.json",
    }:
        raise ValueError("Role metadata namespaces differ")
    if not {"system_usr", "system_etc", "system_conda"} <= system_roles:
        raise ValueError("Required read-only system mounts are missing")
    if any(role in {"frozen_model", "evaluator_only", "target_complete", "truth", "labels"}
           for role, *_ in observed):
        raise ValueError("Forbidden model/evaluator/truth role is mounted")
    return manifest


def expected_raw_input_hashes(
    *, role_root: Path, protocol: Path, execution: Path
) -> dict[str, str]:
    benchmark = role_root / "blind_benchmark"
    benchmark_manifest = load_json(benchmark / "blind_manifest.json")
    fields, rows = read_tsv(protocol / "role_manifest.tsv")
    required_fields = {"role", "bundle_id", "artifact", "sha256"}
    if not required_fields <= set(fields):
        raise ValueError("Frozen role manifest lacks raw-lineage fields")
    candidate_hashes: dict[str, str] = {}
    for row in rows:
        if row["role"] != "candidate_reference":
            continue
        bundle = row["bundle_id"]
        artifact = row["artifact"]
        if bundle not in {"candidate_mandshurica", "candidate_carya"}:
            raise ValueError("Unexpected candidate reference in frozen role manifest")
        if artifact not in {"genome", "gff3", "protein"}:
            raise ValueError("Unexpected candidate artifact in frozen role manifest")
        key = f"{bundle}_{artifact}_sha256"
        if key in candidate_hashes:
            raise ValueError("Duplicate candidate artifact in frozen role manifest")
        candidate_hashes[key] = row["sha256"]
    if len(candidate_hashes) != 6:
        raise ValueError("Frozen role manifest lacks exact two candidate triplets")
    role_manifest = load_json(role_root / "role_manifest.json")
    expected = {
        "staged_input_SHA256SUMS_sha256": role_manifest[
            "staged_input_SHA256SUMS_sha256"
        ],
        "blind_benchmark_SHA256SUMS_sha256": sha256_file(
            benchmark / "SHA256SUMS"
        ),
        "protocol_SHA256SUMS_sha256": sha256_file(protocol / "SHA256SUMS"),
        "execution_SHA256SUMS_sha256": sha256_file(execution / "SHA256SUMS"),
        "target_genome_sha256": benchmark_manifest.get("target_genome", {}).get(
            "sha256"
        ),
        "perturbed_gff3_sha256": sha256_file(benchmark / "perturbed.gff3"),
        **candidate_hashes,
    }
    if any(
        not isinstance(value, str)
        or len(value) != 64
        or any(character not in "0123456789abcdef" for character in value)
        for value in expected.values()
    ):
        raise ValueError("Raw prediction input lineage contains a malformed SHA-256")
    return expected


def validate_raw_manifest(
    path: Path, blind_project: Path, *, expected_input_hashes: dict[str, str]
) -> dict[str, Any]:
    manifest = load_json(path)
    trees = manifest.get("raw_prediction_trees")
    if (
        manifest.get("schema_version") != RAW_SCHEMA
        or manifest.get("holdout_id") != HOLDOUT_ID
        or manifest.get("policy_id") != POLICY_ID
        or manifest.get("truth_access") is not False
        or manifest.get("ranker_access") is not False
        or manifest.get("method_families") != ["miniprot", "gemoma", "lifton"]
        or manifest.get("candidate_references")
        != ["candidate_mandshurica", "candidate_carya"]
        or manifest.get("within_method_reference_vote_count") != 1
        or manifest.get("input_hashes") != expected_input_hashes
        or manifest.get("tree_hash_algorithm")
        != "sha256_of_sorted_sha256_two_space_relative_path_newline"
        or not isinstance(trees, dict) or set(trees) != RAW_PREDICTION_TREE_KEYS
    ):
        raise ValueError("Raw prediction manifest schema/firewall differs")
    relative_paths: set[str] = set()
    for name, item in trees.items():
        if not isinstance(name, str) or not isinstance(item, dict) or set(item) != {
            "relative_path", "file_count", "bytes", "tree_sha256"
        }:
            raise ValueError("Malformed raw prediction tree binding")
        if item["relative_path"] in relative_paths:
            raise ValueError("Raw prediction tree paths must be one-to-one")
        relative_paths.add(item["relative_path"])
        root = safe_join(blind_project, item["relative_path"], "raw prediction tree")
        observed = tree_digest(root)
        if observed != (item["tree_sha256"], item["file_count"], item["bytes"]):
            raise ValueError(f"Raw prediction tree differs: {name}")
    return manifest


def validate_pool_manifest(
    path: Path, *, pool: Path, decisions: Path, raw_sha: str, arm: str
) -> None:
    manifest = load_json(path)
    if (
        arm not in POOL_SCHEMAS
        or manifest.get("schema_version") != POOL_SCHEMAS[arm]
        or manifest.get("holdout_id") != HOLDOUT_ID
        or manifest.get("policy_id") != POLICY_ID
        or manifest.get("truth_access") is not False
        or manifest.get("ranker_access") is not False
        or manifest.get("automatic_approval") is not False
        or manifest.get("within_method_reference_vote_count") != 1
        or manifest.get("candidate_references_per_method") != 2
        or manifest.get("policy_arm") != arm
        or manifest.get("inputs", {}).get("raw_predictions_manifest")
        != {"file_name": "raw_predictions.manifest.json", "sha256": raw_sha}
        or manifest.get("outputs", {}).get("candidate_gff", {}).get("sha256")
        != sha256_file(pool)
        or manifest.get("outputs", {}).get("decisions", {}).get("sha256")
        != sha256_file(decisions)
    ):
        raise ValueError(f"Blind {arm} pool manifest lineage differs")


def finalize_custody(
    *, blind_project: Path, blind_output_root: Path, execution: Path, protocol: Path,
    role_root: Path, runner_command: Path, mount_manifest_path: Path,
    namespace_validation_path: Path, bwrap_version: str, output: Path
) -> Path:
    execution_manifest, protocol_manifest, contract = verify_execution(execution, protocol)
    verify_sha256sums(role_root, ignore_checksum_file=True)
    role_manifest = load_json(role_root / "role_manifest.json")
    if (
        role_manifest.get("schema_version") != ROLE_SCHEMA
        or role_manifest.get("holdout_id") != contract.holdout_id
        or role_manifest.get("protocol_SHA256SUMS_sha256")
        != sha256_file(protocol / "SHA256SUMS")
        or role_manifest.get("roles")
        != ["shared_target", "candidate_only", "blind_benchmark"]
        or any(role_manifest.get(key) is not False for key in (
            "truth_access", "complete_target_annotation_present",
            "evaluator_references_present", "nas_data_present", "network_access",
            "ranker_or_model_present", "h2_or_topology_ranking_present",
        ))
    ):
        raise ValueError("Blind role root violates Walnut firewall")
    mount_manifest = validate_mount_manifest(
        mount_manifest_path, role_root=role_root, execution=execution,
        protocol=protocol, blind_output_root=blind_output_root,
    )
    namespace = load_json(namespace_validation_path)
    if (
        namespace.get("schema_version") != NAMESPACE_SCHEMA
        or namespace.get("holdout_id") != contract.holdout_id
        or namespace.get("mount_manifest_sha256") != sha256_file(mount_manifest_path)
        or namespace.get("shared_target_visible") is not True
        or namespace.get("candidate_only_visible") is not True
        or namespace.get("blind_benchmark_visible") is not True
        or any(namespace.get(key) is not False for key in (
            "evaluator_only_visible", "complete_target_annotation_visible",
            "truth_visible", "labels_visible", "nas_data_visible", "model_visible",
            "ranker_visible",
        ))
    ):
        raise ValueError("Blind namespace negative-access validation differs")

    paths = {name: safe_join(blind_project, relative, f"blind output {name}")
             for name, relative in BLIND_OUTPUTS.items()}
    for name, path in paths.items():
        if not path.is_file() or path.is_symlink() or path.stat().st_size == 0:
            raise ValueError(f"Missing, empty or symlinked blind output: {name}")
    raw_manifest = validate_raw_manifest(
        paths["raw_predictions_manifest"],
        blind_project,
        expected_input_hashes=expected_raw_input_hashes(
            role_root=role_root, protocol=protocol, execution=execution
        ),
    )
    raw_sha = sha256_file(paths["raw_predictions_manifest"])
    validate_pool_manifest(
        paths["retain_manifest"], pool=paths["retain_pool"],
        decisions=paths["retain_decisions"], raw_sha=raw_sha,
        arm="retain_distinct_phased_CDS_chains",
    )
    validate_pool_manifest(
        paths["suppress_manifest"], pool=paths["suppress_pool"],
        decisions=paths["suppress_decisions"], raw_sha=raw_sha,
        arm="suppress_strongly_overlapping_alternative_chains",
    )
    for decision in (paths["retain_decisions"], paths["suppress_decisions"]):
        fields, rows = read_tsv(decision)
        if not rows or not {"consensus_digest", "status"} <= set(fields):
            raise ValueError("Blind decisions are empty or malformed")
    for command in (runner_command, paths["command_log"]):
        if not command.is_file() or command.is_symlink():
            raise ValueError("Missing blind command audit")
        reject_forbidden_text(command)
    runner_text = runner_command.read_text(encoding="utf-8")
    for token in ("--unshare-all", "--unshare-net", "--clearenv",
                  "PLOIDYPATCH_NETWORK_ACCESS", "none"):
        if token not in runner_text:
            raise ValueError(f"Blind launcher command lacks {token}")

    tree_sha, file_count, byte_count = tree_digest(blind_project)
    custody = {
        "schema_version": CUSTODY_SCHEMA,
        "holdout_id": contract.holdout_id,
        "policy_id": contract.policy_id,
        "protocol_profile": "core_H1_only_no_ranker",
        "frozen_before_reveal_at": datetime.now(timezone.utc).replace(
            microsecond=0
        ).isoformat().replace("+00:00", "Z"),
        "protocol_SHA256SUMS_sha256": sha256_file(protocol / "SHA256SUMS"),
        "execution_SHA256SUMS_sha256": sha256_file(execution / "SHA256SUMS"),
        "role_root_SHA256SUMS_sha256": sha256_file(role_root / "SHA256SUMS"),
        "ranker_or_model_executed": False,
        "h2_or_topology_ranking_executed": False,
        "truth_mounted": False,
        "complete_target_annotation_mounted": False,
        "evaluator_references_mounted": False,
        "nas_data_mounted": False,
        "network_access": False,
        "bubblewrap": {"version": bwrap_version, "command_sha256": sha256_file(runner_command)},
        "mount_manifest_sha256": sha256_file(mount_manifest_path),
        "namespace_validation_sha256": sha256_file(namespace_validation_path),
        "raw_predictions_manifest_sha256": sha256_file(paths["raw_predictions_manifest"]),
        "raw_prediction_tree_count": len(raw_manifest["raw_prediction_trees"]),
        "blind_project": {"tree_sha256": tree_sha, "file_count": file_count,
                          "bytes": byte_count},
        "blind_outputs": {
            name: {"relative_path": path.relative_to(blind_project).as_posix(),
                   "sha256": sha256_file(path), "bytes": path.stat().st_size}
            for name, path in paths.items()
        },
    }
    if output != blind_output_root / "custody_manifest.json":
        raise ValueError("Custody output must be blind root/custody_manifest.json")
    if output.exists() or output.is_symlink() or (blind_output_root / "SHA256SUMS").exists():
        raise FileExistsError("Refusing to overwrite blind custody")
    output.parent.mkdir(parents=True, exist_ok=True)
    descriptor = os.open(output, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o440)
    with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as handle:
        json.dump(custody, handle, indent=2, sort_keys=True); handle.write("\n")
    write_sha256sums(blind_output_root)
    verify_sha256sums(blind_output_root, ignore_checksum_file=True)
    return output


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--blind-project-root", required=True)
    parser.add_argument("--blind-output-root", required=True)
    parser.add_argument("--execution-freeze", required=True)
    parser.add_argument("--protocol-freeze", required=True)
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
