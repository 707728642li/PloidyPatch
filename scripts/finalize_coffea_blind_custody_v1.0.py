#!/usr/bin/env python3
"""Seal Coffea blind outputs and negative-access evidence before reveal."""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path, PurePosixPath
from typing import Any, Iterable

from ploidypatch.artifact_manifest import sha256_file, verify_sha256sums
from ploidypatch.coffea_h1_framework import (
    BLIND_OUTPUTS,
    CUSTODY_SCHEMA,
    HOLDOUT_ID,
    MOUNT_SCHEMA,
    NAMESPACE_SCHEMA,
    PATCH2_STAGE,
    PATCH3_STAGE,
    POLICY_ID,
    RAW_PREDICTION_TREE_KEYS,
    RAW_SCHEMA,
    ROLE_SCHEMA,
    load_json,
    read_tsv,
    verify_execution,
)


def safe_join(root: Path, raw: str, context: str) -> Path:
    relative = PurePosixPath(raw)
    if relative.is_absolute() or any(part in {"", ".", ".."} for part in relative.parts):
        raise ValueError(f"Unsafe relative path for {context}: {raw}")
    path = root.joinpath(*relative.parts)
    try:
        path.resolve().relative_to(root.resolve())
    except ValueError:
        raise ValueError(f"{context} escapes blind project") from None
    return path


def tree_digest(root: Path) -> tuple[str, int, int]:
    digest = hashlib.sha256(); count = 0; byte_count = 0
    for path in sorted(root.rglob("*")):
        if path.is_symlink():
            raise ValueError(f"Symlink is forbidden in Coffea blind tree: {path}")
        if not path.is_file():
            continue
        relative = path.relative_to(root).as_posix()
        size = path.stat().st_size
        digest.update(f"{sha256_file(path)}  {relative}\n".encode())
        count += 1; byte_count += size
    if count == 0:
        raise ValueError(f"Empty Coffea tree: {root}")
    return digest.hexdigest(), count, byte_count


def expected_raw_hashes(role_root: Path, protocol: Path, execution: Path) -> dict[str, str]:
    fields, rows = read_tsv(role_root / "role_manifest.tsv")
    if not {"role", "bundle_id", "artifact", "sha256"} <= set(fields):
        raise ValueError("Coffea blind role manifest lacks lineage fields")
    candidate: dict[str, str] = {}
    for row in rows:
        if row["role"] != "candidate_reference":
            continue
        key = f"{row['bundle_id']}_{row['artifact']}_sha256"
        if row["bundle_id"] not in {"candidate_bua", "candidate_mauritiana"}:
            raise ValueError("Unexpected Coffea candidate bundle")
        if row["artifact"] not in {"genome", "gff3", "protein"} or key in candidate:
            raise ValueError("Malformed Coffea candidate artifact binding")
        candidate[key] = row["sha256"]
    if len(candidate) != 6:
        raise ValueError("Coffea blind role lacks exact two candidate triplets")
    role = load_json(role_root / "role_manifest.json")
    benchmark = role_root / "blind_benchmark"
    benchmark_manifest = load_json(benchmark / "blind_manifest.json")
    expected = {
        "staged_input_SHA256SUMS_sha256": role["staged_input_SHA256SUMS_sha256"],
        "blind_benchmark_SHA256SUMS_sha256": sha256_file(benchmark / "SHA256SUMS"),
        "protocol_SHA256SUMS_sha256": sha256_file(protocol / "SHA256SUMS"),
        "execution_SHA256SUMS_sha256": sha256_file(execution / "SHA256SUMS"),
        "target_genome_sha256": benchmark_manifest["target_genome"]["sha256"],
        "perturbed_gff3_sha256": sha256_file(benchmark / "perturbed.gff3"),
        **candidate,
    }
    if any(not isinstance(value, str) or len(value) != 64 for value in expected.values()):
        raise ValueError("Malformed Coffea raw input hash")
    return expected


def validate_raw_manifest(path: Path, project: Path, expected_hashes: dict[str, str]) -> None:
    raw = load_json(path)
    trees = raw.get("raw_prediction_trees")
    if (
        raw.get("schema_version") != RAW_SCHEMA
        or raw.get("holdout_id") != HOLDOUT_ID
        or raw.get("policy_id") != POLICY_ID
        or raw.get("truth_access") is not False
        or raw.get("ranker_access") is not False
        or raw.get("candidate_references")
        != ["candidate_bua", "candidate_mauritiana"]
        or raw.get("within_method_reference_vote_count") != 1
        or raw.get("input_hashes") != expected_hashes
        or not isinstance(trees, dict)
        or set(trees) != RAW_PREDICTION_TREE_KEYS
    ):
        raise ValueError("Coffea raw prediction manifest differs")
    for label, item in trees.items():
        if not isinstance(item, dict) or set(item) != {
            "relative_path", "file_count", "bytes", "tree_sha256"
        }:
            raise ValueError(f"Malformed Coffea raw tree binding: {label}")
        observed = tree_digest(safe_join(project, item["relative_path"], label))
        if observed != (item["tree_sha256"], item["file_count"], item["bytes"]):
            raise ValueError(f"Coffea raw prediction tree differs: {label}")


def validate_pool(
    *, manifest: Path, pool: Path, decisions: Path, raw_sha: str,
    scope: str, arm: str
) -> None:
    value = load_json(manifest)
    schema = (
        "ploidypatch.method_candidate_pool.v2"
        if arm == "retain_distinct"
        else "ploidypatch.method_consensus.v1"
    )
    policy_arm = (
        "retain_distinct_phased_CDS_chains"
        if arm == "retain_distinct"
        else "suppress_strongly_overlapping_alternative_chains"
    )
    expected_count = 2 if scope == "combined" else 1
    if (
        value.get("schema_version") != schema
        or value.get("holdout_id") != HOLDOUT_ID
        or value.get("policy_id") != POLICY_ID
        or value.get("policy_arm") != policy_arm
        or value.get("reference_scope") != scope
        or value.get("candidate_references_per_method") != expected_count
        or value.get("within_method_reference_vote_count") != 1
        or value.get("truth_access") is not False
        or value.get("ranker_access") is not False
        or value.get("automatic_approval") is not False
        or value.get("inputs", {}).get("raw_predictions_manifest")
        != {"file_name": "raw_predictions.manifest.json", "sha256": raw_sha}
        or value.get("outputs", {}).get("candidate_gff", {}).get("sha256")
        != sha256_file(pool)
        or value.get("outputs", {}).get("decisions", {}).get("sha256")
        != sha256_file(decisions)
    ):
        raise ValueError(f"Coffea blind pool lineage differs: {scope}/{arm}")
    fields, rows = read_tsv(decisions)
    forbidden = {field.casefold() for field in fields} & {
        "label", "truth", "is_positive", "target_label"
    }
    if not rows or not {"consensus_digest", "status"} <= set(fields) or forbidden:
        raise ValueError(f"Coffea decisions malformed or truth-bearing: {scope}/{arm}")


def validate_patch2_primary_replay(
    execution_manifest: dict[str, Any], paths: dict[str, Path]
) -> None:
    if execution_manifest.get("freeze_stage") != PATCH2_STAGE:
        return
    patch = execution_manifest.get("execution_patch")
    expected = (
        patch.get("primary_combined_replay_outputs")
        if isinstance(patch, dict)
        else None
    )
    keys = {
        "combined_retain_distinct_pool",
        "combined_retain_distinct_decisions",
        "combined_suppress_overlap_pool",
        "combined_suppress_overlap_decisions",
    }
    if not isinstance(expected, dict) or set(expected) != keys:
        raise ValueError("Coffea patch-2 primary replay binding differs")
    observed = {name: sha256_file(paths[name]) for name in sorted(keys)}
    if observed != expected:
        raise ValueError("Coffea patch-2 changed a primary combined blind output")


def validate_patch3_reconciliation(
    execution_manifest: dict[str, Any], path: Path, raw_manifest: Path
) -> dict[str, Any] | None:
    if execution_manifest.get("freeze_stage") != PATCH3_STAGE:
        return None
    value = load_json(path)
    patch = execution_manifest.get("execution_patch")
    raw = load_json(raw_manifest)
    lineage = raw.get("reproducibility_reconciliation")
    arms = value.get("arms")
    expected_arms = {
        f"{method}__{bundle}"
        for method in ("miniprot", "gemoma", "lifton")
        for bundle in ("candidate_bua", "candidate_mauritiana")
    }
    if (
        not isinstance(patch, dict)
        or value.get("schema_version")
        != "ploidypatch.coffea_blind_reproducibility_reconciliation.v1.1"
        or value.get("holdout_id") != HOLDOUT_ID
        or value.get("policy_id") != POLICY_ID
        or value.get("truth_access") is not False
        or value.get("label_access") is not False
        or value.get("ranker_access") is not False
        or value.get("automatic_approval") is not False
        or value.get("biological_rules_or_thresholds_changed") is not False
        or value.get("unstable_projection_policy") != "abstain"
        or value.get("selection_rule")
        != "exact_complete_adapter_decision_row_intersection_across_two_independent_blind_runs"
        or not isinstance(arms, dict)
        or set(arms) != expected_arms
        or not isinstance(lineage, dict)
        or lineage.get("execution_patch_sequence") != 3
        or lineage.get("selection_rule")
        != "same_model_id_and_exact_complete_decision_row_in_both_runs"
        or lineage.get("run_a_manifest_sha256")
        != patch.get("reproducibility_run_a_manifest_sha256")
        or lineage.get("run_b_manifest_sha256")
        != patch.get("reproducibility_run_b_manifest_sha256")
    ):
        raise ValueError("Coffea patch-3 reconciliation lineage differs")
    for key, arm in arms.items():
        if (
            not isinstance(arm, dict)
            or arm.get("label_access") is not False
            or arm.get("truth_access") is not False
            or arm.get("selection_rule")
            != "same_model_id_and_exact_complete_decision_row_in_both_runs"
            or not isinstance(arm.get("unstable_models"), int)
            or arm["unstable_models"] < 0
            or arm.get("stable_accepted_models") != arm.get("kept_models")
        ):
            raise ValueError(f"Malformed Coffea reconciliation arm: {key}")
    return value


def validate_isolation(
    mount_path: Path, namespace_path: Path, role_root: Path, *, patch3: bool
) -> None:
    mount = load_json(mount_path); namespace = load_json(namespace_path)
    mounts = mount.get("mounts")
    if mount.get("schema_version") != MOUNT_SCHEMA or not isinstance(mounts, list):
        raise ValueError("Coffea mount manifest differs")
    observed = {item.get("role") for item in mounts if isinstance(item, dict)}
    if not {"shared_target", "candidate_only", "blind_benchmark"} <= observed:
        raise ValueError("Coffea blind namespace lacks exact safe roles")
    if patch3 and not {"reproducibility_run_a", "reproducibility_run_b"} <= observed:
        raise ValueError("Coffea reconciliation namespace lacks both blind runs")
    for item in mounts:
        text = f"{item.get('host_path', '')}\n{item.get('namespace_path', '')}".casefold()
        if any(token in text for token in ("/nas_data", "evaluator_only", "/truth", "/labels", "target_complete")):
            raise ValueError("Coffea blind mount exposes forbidden data")
    if (
        namespace.get("schema_version") != NAMESPACE_SCHEMA
        or namespace.get("shared_target_visible") is not True
        or namespace.get("candidate_only_visible") is not True
        or namespace.get("blind_benchmark_visible") is not True
        or any(namespace.get(key) is not False for key in (
            "evaluator_only_visible", "complete_target_annotation_visible",
            "truth_visible", "labels_visible", "nas_data_visible",
            "model_visible", "ranker_visible",
        ))
        or namespace.get("role_root_SHA256SUMS_sha256")
        != sha256_file(role_root / "SHA256SUMS")
        or (
            patch3
            and (
                namespace.get("reproducibility_run_a_visible") is not True
                or namespace.get("reproducibility_run_b_visible") is not True
            )
        )
    ):
        raise ValueError("Coffea namespace negative-access evidence differs")


def finalize(
    *, blind_project: Path, blind_run: Path, execution: Path, protocol: Path,
    role_root: Path, runner_command: Path, mount_manifest: Path,
    namespace_validation: Path, output: Path
) -> Path:
    execution_manifest, _, _ = verify_execution(execution, protocol)
    verify_sha256sums(role_root, ignore_checksum_file=True)
    role = load_json(role_root / "role_manifest.json")
    if (
        role.get("schema_version") != ROLE_SCHEMA
        or role.get("roles") != ["shared_target", "candidate_only", "blind_benchmark"]
        or any(role.get(key) is not False for key in (
            "truth_access", "complete_target_annotation_present",
            "evaluator_references_present", "nas_data_present", "network_access",
            "ranker_or_model_present", "h2_or_topology_ranking_present",
        ))
    ):
        raise ValueError("Coffea blind role root violates firewall")
    patch3 = execution_manifest.get("freeze_stage") == PATCH3_STAGE
    validate_isolation(
        mount_manifest, namespace_validation, role_root, patch3=patch3
    )
    paths = {
        name: safe_join(blind_project, relative, name)
        for name, relative in BLIND_OUTPUTS.items()
    }
    for name, path in paths.items():
        if not path.is_file() or path.is_symlink() or path.stat().st_size == 0:
            raise ValueError(f"Missing Coffea blind output: {name}")
    validate_patch2_primary_replay(execution_manifest, paths)
    validate_raw_manifest(
        paths["raw_predictions_manifest"], blind_project,
        expected_raw_hashes(role_root, protocol, execution),
    )
    raw_sha = sha256_file(paths["raw_predictions_manifest"])
    reconciliation_path = (
        blind_project
        / "results/copy_collapse/external/coffea_v1.0_h1/reconciliation_manifest.json"
    )
    reconciliation = validate_patch3_reconciliation(
        execution_manifest, reconciliation_path, paths["raw_predictions_manifest"]
    )
    for scope in ("combined", "bua_only", "mauritiana_only"):
        for arm in ("retain_distinct", "suppress_overlap"):
            prefix = f"{scope}_{arm}"
            validate_pool(
                manifest=paths[f"{prefix}_manifest"],
                pool=paths[f"{prefix}_pool"],
                decisions=paths[f"{prefix}_decisions"],
                raw_sha=raw_sha, scope=scope, arm=arm,
            )
    for command in (runner_command, paths["command_log"]):
        text = command.read_text(encoding="utf-8").casefold()
        if any(token in text for token in ("/nas_data", "evaluator_only", "/truth", "/labels", "ranker")):
            raise ValueError(f"Forbidden Coffea command text: {command}")
    if output.exists() or output.is_symlink():
        raise FileExistsError(f"Refusing to overwrite Coffea custody: {output}")
    tree_sha, files, byte_count = tree_digest(blind_project)
    payload = {
        "schema_version": CUSTODY_SCHEMA,
        "holdout_id": HOLDOUT_ID,
        "policy_id": POLICY_ID,
        "frozen_before_reveal_at": datetime.now(timezone.utc).replace(
            microsecond=0
        ).isoformat().replace("+00:00", "Z"),
        "protocol_SHA256SUMS_sha256": sha256_file(protocol / "SHA256SUMS"),
        "execution_SHA256SUMS_sha256": sha256_file(execution / "SHA256SUMS"),
        "role_root_SHA256SUMS_sha256": sha256_file(role_root / "SHA256SUMS"),
        "blind_project_tree_sha256": tree_sha,
        "blind_project_files": files,
        "blind_project_bytes": byte_count,
        "blind_outputs": {
            name: {
                "relative_path": BLIND_OUTPUTS[name],
                "sha256": sha256_file(path),
                "bytes": path.stat().st_size,
            }
            for name, path in sorted(paths.items())
        },
        "mount_manifest_sha256": sha256_file(mount_manifest),
        "namespace_validation_sha256": sha256_file(namespace_validation),
        "ranker_or_model_executed": False,
        "h2_or_topology_ranking_executed": False,
        "truth_mounted": False,
        "complete_target_annotation_mounted": False,
        "evaluator_references_mounted": False,
        "nas_data_mounted": False,
        "network_access": False,
        "automatic_approval": False,
    }
    if reconciliation is not None:
        payload["reproducibility_reconciliation"] = {
            "relative_path": reconciliation_path.relative_to(blind_project).as_posix(),
            "sha256": sha256_file(reconciliation_path),
            "unstable_models": sum(
                int(arm["unstable_models"])
                for arm in reconciliation["arms"].values()
            ),
        }
    output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return output


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--blind-project", required=True)
    parser.add_argument("--blind-run", required=True)
    parser.add_argument("--execution-freeze", required=True)
    parser.add_argument("--protocol-freeze", required=True)
    parser.add_argument("--blind-role-root", required=True)
    parser.add_argument("--runner-command", required=True)
    parser.add_argument("--mount-manifest", required=True)
    parser.add_argument("--namespace-validation", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args(argv)
    finalize(
        blind_project=Path(args.blind_project).resolve(),
        blind_run=Path(args.blind_run).resolve(),
        execution=Path(args.execution_freeze).resolve(),
        protocol=Path(args.protocol_freeze).resolve(),
        role_root=Path(args.blind_role_root).resolve(),
        runner_command=Path(args.runner_command).resolve(),
        mount_manifest=Path(args.mount_manifest).resolve(),
        namespace_validation=Path(args.namespace_validation).resolve(),
        output=Path(args.output).resolve(),
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
