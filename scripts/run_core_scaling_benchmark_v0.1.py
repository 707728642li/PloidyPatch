#!/usr/bin/env python3
"""Run the frozen multi-species PloidyPatch core-scaling benchmark."""

from __future__ import annotations

import argparse
import os
import random
import shutil
import subprocess
import sys
import time
from pathlib import Path

import core_scaling_common_v0_1 as common
from core_scaling_common_v0_1 import (
    COMPONENTS,
    HEX40,
    REFERENCE_SCHEMA,
    RESULT_SCHEMA,
    ScalingError,
    canonical_json,
    command_specs,
    environment_manifest,
    freeze_read_only,
    load_config,
    load_json,
    require,
    resolve_project_file,
    run_specs,
    safely_extract_archive,
    sha256_file,
    verify_artifacts,
    verify_sha256sums,
    write_sha256sums,
    write_tsv,
)


SEED = 20261006


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", required=True, type=Path)
    parser.add_argument("--python-bin", required=True, type=Path)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--reference-root", required=True, type=Path)
    parser.add_argument("--code-commit", required=True)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--warm-repeats", type=int, default=3)
    parser.add_argument("--seed", type=int, default=SEED)
    parser.add_argument("--conda", default="conda")
    return parser.parse_args()


def current_environment_locks(conda: str, python_bin: Path) -> dict[str, str]:
    env_prefix = python_bin.resolve().parent.parent
    commands = {
        "environment.conda-explicit.txt": [conda, "list", "--explicit", "-p", str(env_prefix)],
        "environment.pip-freeze.txt": [str(python_bin), "-m", "pip", "freeze", "--all"],
    }
    output: dict[str, str] = {}
    for name, argv in commands.items():
        completed = subprocess.run(argv, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        require(completed.returncode == 0, f"environment query failed: {completed.stderr}")
        output[name] = completed.stdout
    return output


def replicate_order(workloads: list[dict], phase: str, repeat: int, seed: int) -> list[dict]:
    ordered = list(workloads)
    if phase == "warm_repeats":
        random.Random(seed + repeat).shuffle(ordered)
    return ordered


def main() -> int:
    args = parse_args()
    require(args.warm_repeats == 3, "v0.1 protocol requires exactly three warm repeats")
    require(args.seed == SEED, f"v0.1 protocol requires seed {SEED}")
    require(HEX40.fullmatch(args.code_commit) is not None, "code commit must be a full lowercase SHA")
    project_root = args.project_root.resolve(strict=True)
    python_bin = args.python_bin.resolve(strict=True)
    config_path = args.config.resolve(strict=True)
    reference_root = args.reference_root.resolve(strict=True)
    output = args.output_dir.resolve()
    working = Path(str(output) + ".working")
    require(not output.exists() and not working.exists(), f"refusing to overwrite: {output}")
    verify_sha256sums(reference_root)
    config = load_config(config_path, project_root)
    registry = load_json(reference_root / "reference_registry.json")
    require(registry.get("schema_version") == REFERENCE_SCHEMA, "invalid reference registry schema")
    code = registry.get("code", {})
    require(code.get("code_commit") == args.code_commit, "code commit differs from reference freeze")
    require(code.get("workload_config_sha256") == sha256_file(config_path), "config differs from reference freeze")
    model_path = resolve_project_file(project_root, config["model_v03"])
    require(code.get("model_v03_sha256") == sha256_file(model_path), "model differs from reference freeze")
    archive = reference_root / "code_source.tar.gz"
    require(code.get("source_archive_sha256") == sha256_file(archive), "source archive differs")
    locks = current_environment_locks(args.conda, python_bin)
    for name, text in locks.items():
        require((reference_root / name).read_text(encoding="utf-8") == text, f"environment lock differs: {name}")
    working.mkdir(parents=True)
    any_execution = False
    try:
        frozen_code = working / "frozen_code"
        safely_extract_archive(archive, frozen_code)
        for relative, live in (
            ("scripts/run_core_scaling_benchmark_v0.1.py", Path(__file__).resolve()),
            ("scripts/core_scaling_common_v0_1.py", Path(common.__file__).resolve()),
        ):
            require(sha256_file(frozen_code / relative) == sha256_file(live), f"orchestrator differs from frozen code: {relative}")
        (working / "environment.json").write_text(canonical_json(environment_manifest()), encoding="utf-8")
        (working / "run_contract.json").write_text(
            canonical_json(
                {
                    "schema_version": RESULT_SCHEMA,
                    "code_commit": args.code_commit,
                    "reference_registry_sha256": sha256_file(reference_root / "reference_registry.json"),
                    "reference_sha256sums_sha256": sha256_file(reference_root / "SHA256SUMS"),
                    "seed": args.seed,
                    "cold_executions": 1,
                    "warm_repeats": args.warm_repeats,
                    "cache_policy": "no_privileged_cache_drop; first execution labelled cold_order_first",
                    "thread_scaling": "not_applicable_core_components_expose_no_thread_parameter",
                    "upstream_tools_included": False,
                    "truth_access": False,
                    "automatic_approval": False,
                }
            ),
            encoding="utf-8",
        )
        workload_metadata = {item["workload_id"]: item for item in registry["workloads"]}
        require(set(workload_metadata) == {item["workload_id"] for item in config["workloads"]}, "reference/config workload universe differs")
        rows: list[dict[str, object]] = []
        env = os.environ.copy()
        env["PYTHONPATH"] = str(frozen_code / "src")
        env["PYTHONHASHSEED"] = "0"
        env["PYTHONDONTWRITEBYTECODE"] = "1"
        phases = [("cold_order_first", 0)] + [("warm_repeats", repeat) for repeat in range(1, args.warm_repeats + 1)]
        for phase, repeat in phases:
            ordered = replicate_order(config["workloads"], phase, repeat, args.seed)
            for order_index, workload in enumerate(ordered, start=1):
                workload_id = workload["workload_id"]
                metadata = workload_metadata[workload_id]
                for component in COMPONENTS:
                    run_root = working / "replicates" / phase / f"repeat_{repeat}" / workload_id / component
                    artifact_root = run_root / "artifacts"
                    log_root = run_root / "logs"
                    status = "pass"
                    reason = ""
                    metrics: dict[str, object] = {}
                    output_bytes = 0
                    output_digest = ""
                    started = time.monotonic()
                    try:
                        specs = command_specs(
                            component=component,
                            workload=workload,
                            project_root=project_root,
                            python_bin=python_bin,
                            artifact_root=artifact_root,
                            reference_root=reference_root,
                            model_path=model_path,
                        )
                        any_execution = True
                        metrics = run_specs(specs, code_root=frozen_code, log_root=log_root, env=env)
                        output_bytes, output_digest = verify_artifacts(
                            component=component,
                            artifact_root=artifact_root,
                            reference_workload_root=reference_root / "workloads" / workload_id,
                        )
                    except Exception as exc:  # preserve and report every failed replicate
                        status = "failed"
                        reason = f"{type(exc).__name__}: {exc}"
                        metrics.setdefault("wall_seconds", time.monotonic() - started)
                    candidate_count = int(metadata["pool"]["candidate_count"])
                    conflict_sets = int(metadata["pool"]["conflict_sets"])
                    wall_seconds = float(metrics["wall_seconds"])
                    rows.append(
                        {
                            "phase": phase,
                            "repeat": repeat,
                            "order_index": order_index,
                            "workload_id": workload_id,
                            "species": workload["species"],
                            "component": component,
                            "status": status,
                            "failure_reason": reason,
                            "candidate_count": candidate_count,
                            "conflict_sets": conflict_sets,
                            "wall_seconds": format(wall_seconds, ".9f"),
                            "user_seconds": "" if "user_seconds" not in metrics else format(float(metrics["user_seconds"]), ".9f"),
                            "system_seconds": "" if "system_seconds" not in metrics else format(float(metrics["system_seconds"]), ".9f"),
                            "peak_rss_kb": metrics.get("peak_rss_kb", ""),
                            "filesystem_inputs": metrics.get("filesystem_inputs", ""),
                            "filesystem_outputs": metrics.get("filesystem_outputs", ""),
                            "output_bytes": output_bytes,
                            "candidates_per_second": format(candidate_count / wall_seconds, ".9f") if wall_seconds else "",
                            "conflict_sets_per_second": format(conflict_sets / wall_seconds, ".9f") if wall_seconds else "",
                            "output_digest": output_digest,
                        }
                    )
        fields = [
            "phase", "repeat", "order_index", "workload_id", "species", "component", "status",
            "failure_reason", "candidate_count", "conflict_sets", "wall_seconds", "user_seconds",
            "system_seconds", "peak_rss_kb", "filesystem_inputs", "filesystem_outputs", "output_bytes",
            "candidates_per_second", "conflict_sets_per_second", "output_digest",
        ]
        write_tsv(working / "replicates.tsv", fields, rows)
        failed = [row for row in rows if row["status"] != "pass"]
        summary = {
            "schema_version": RESULT_SCHEMA,
            "replicates": len(rows),
            "passed": len(rows) - len(failed),
            "failed": len(failed),
            "all_outputs_exact": not failed,
            "largest_verified_workload": max(
                (
                    {
                        "workload_id": item["workload_id"],
                        "candidate_count": int(item["pool"]["candidate_count"]),
                    }
                    for item in registry["workloads"]
                    if not any(row["workload_id"] == item["workload_id"] and row["status"] != "pass" for row in rows)
                ),
                key=lambda item: item["candidate_count"],
                default=None,
            ),
            "thread_scaling": "not_applicable_core_components_expose_no_thread_parameter",
        }
        (working / "summary.json").write_text(canonical_json(summary), encoding="utf-8")
        write_sha256sums(working)
        working.rename(output)
        freeze_read_only(output)
    except Exception:
        if not any_execution and working.is_dir():
            shutil.rmtree(working)
        else:
            print(f"benchmark failed; retained attempt at {working}", file=sys.stderr)
        raise
    print(canonical_json({"result_root": str(output), **summary}), end="")
    return 0 if summary["failed"] == 0 else 1


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except ScalingError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
