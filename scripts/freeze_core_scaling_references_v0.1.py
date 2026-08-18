#!/usr/bin/env python3
"""Build immutable real-workload references for the PloidyPatch scaling run."""

from __future__ import annotations

import argparse
import csv
import os
import shutil
import subprocess
import sys
from collections import Counter
from pathlib import Path

from core_scaling_common_v0_1 import (
    ARTIFACTS,
    HEX40,
    REFERENCE_SCHEMA,
    ScalingError,
    canonical_json,
    command_specs,
    compare_canonical,
    environment_manifest,
    fasta_stats,
    freeze_read_only,
    gff_stats,
    load_config,
    pool_stats,
    require,
    resolve_project_file,
    run_specs,
    safely_extract_archive,
    sha256_file,
    write_sha256sums,
    write_tsv,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", required=True, type=Path)
    parser.add_argument("--python-bin", required=True, type=Path)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--source-archive", required=True, type=Path)
    parser.add_argument("--code-commit", required=True)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--conda", default="conda")
    return parser.parse_args()


def decision_stats(path: Path) -> dict[str, object]:
    statuses: Counter[str] = Counter()
    conflict_members: dict[str, list[str]] = {}
    declared_conflict_sizes: dict[str, int] = {}
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require(reader.fieldnames is not None, f"missing decisions header: {path}")
        require("status" in reader.fieldnames, f"missing status: {path}")
        for row_number, row in enumerate(reader, start=2):
            statuses[row["status"]] += 1
            if row["status"] == "accepted" and row.get("conflict_set_digest"):
                conflict = row["conflict_set_digest"]
                digest = row.get("consensus_digest") or row.get("candidate_digest")
                require(bool(digest), f"missing candidate digest at {path}:{row_number}")
                declared = int(row["conflict_member_count"])
                if conflict in declared_conflict_sizes:
                    require(
                        declared_conflict_sizes[conflict] == declared,
                        f"inconsistent conflict size in {path}: {conflict}",
                    )
                declared_conflict_sizes[conflict] = declared
                conflict_members.setdefault(conflict, []).append(digest)
    conflict_sets_by_size: Counter[int] = Counter()
    for conflict, members in conflict_members.items():
        require(
            len(members) == len(set(members)) == declared_conflict_sizes[conflict],
            f"conflict membership differs in {path}: {conflict}",
        )
        conflict_sets_by_size[len(members)] += 1
    return {
        "rows": sum(statuses.values()),
        "status_counts": dict(sorted(statuses.items())),
        "conflict_sets_by_size": {
            str(size): count for size, count in sorted(conflict_sets_by_size.items())
        },
    }


def capture_environment_locks(conda: str, python_bin: Path, root: Path) -> None:
    env_prefix = python_bin.resolve().parent.parent
    commands = {
        "environment.conda-explicit.txt": [conda, "list", "--explicit", "-p", str(env_prefix)],
        "environment.pip-freeze.txt": [str(python_bin), "-m", "pip", "freeze", "--all"],
    }
    for name, argv in commands.items():
        completed = subprocess.run(argv, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        require(completed.returncode == 0, f"environment lock failed: {' '.join(argv)}: {completed.stderr}")
        (root / name).write_text(completed.stdout, encoding="utf-8")


def main() -> int:
    args = parse_args()
    project_root = args.project_root.resolve(strict=True)
    python_bin = args.python_bin.resolve(strict=True)
    config_path = args.config.resolve(strict=True)
    source_archive = args.source_archive.resolve(strict=True)
    require(HEX40.fullmatch(args.code_commit) is not None, "code commit must be a full lowercase SHA")
    output = args.output_dir.resolve()
    working = Path(str(output) + ".working")
    require(not output.exists() and not working.exists(), f"refusing to overwrite: {output}")
    config = load_config(config_path, project_root)
    model_path = resolve_project_file(project_root, config["model_v03"])
    working.mkdir(parents=True)
    try:
        shutil.copyfile(config_path, working / "workload_config.json")
        shutil.copyfile(source_archive, working / "code_source.tar.gz")
        frozen_code = working / "frozen_code"
        safely_extract_archive(working / "code_source.tar.gz", frozen_code)
        for relative, live in (
            ("scripts/freeze_core_scaling_references_v0.1.py", Path(__file__).resolve()),
            ("scripts/core_scaling_common_v0_1.py", Path(sys.modules["core_scaling_common_v0_1"].__file__).resolve()),
        ):
            require(
                sha256_file(frozen_code / relative) == sha256_file(live),
                f"reference orchestrator differs from frozen code: {relative}",
            )
        (working / "environment.json").write_text(canonical_json(environment_manifest()), encoding="utf-8")
        capture_environment_locks(args.conda, python_bin, working)
        code_manifest = {
            "code_commit": args.code_commit,
            "source_archive_sha256": sha256_file(working / "code_source.tar.gz"),
            "workload_config_sha256": sha256_file(working / "workload_config.json"),
            "model_v03_sha256": sha256_file(model_path),
            "python_sha256": sha256_file(python_bin),
        }
        (working / "code_manifest.json").write_text(canonical_json(code_manifest), encoding="utf-8")
        registry_workloads: list[dict[str, object]] = []
        input_rows: list[dict[str, object]] = []
        timing_rows: list[dict[str, object]] = []
        env = os.environ.copy()
        env["PYTHONPATH"] = str(frozen_code / "src")
        env["PYTHONHASHSEED"] = "0"
        env["PYTHONDONTWRITEBYTECODE"] = "1"
        for workload in config["workloads"]:
            workload_id = workload["workload_id"]
            reference = working / "workloads" / workload_id
            logs = working / "reference_logs" / workload_id
            specs = command_specs(
                component="end_to_end",
                workload=workload,
                project_root=project_root,
                python_bin=python_bin,
                artifact_root=reference,
                reference_root=working,
                model_path=model_path,
            )
            metrics = run_specs(specs, code_root=frozen_code, log_root=logs, env=env)
            canonical_pool = (project_root / workload["canonical_pool_root"]).resolve(strict=True)
            canonical_features = (project_root / workload["canonical_feature_root"]).resolve(strict=True)
            compare_canonical(reference, canonical_pool, ARTIFACTS["consensus"])
            compare_canonical(reference, canonical_features, ARTIFACTS["feature_construction"])
            if "canonical_v03_root" in workload:
                compare_canonical(
                    reference,
                    (project_root / workload["canonical_v03_root"]).resolve(strict=True),
                    ARTIFACTS["score_v03"],
                )
            if "canonical_v04_root" in workload:
                compare_canonical(
                    reference,
                    (project_root / workload["canonical_v04_root"]).resolve(strict=True),
                    ARTIFACTS["guard_v04"],
                )
            canonical_roots = {
                "consensus": (
                    canonical_pool,
                    ARTIFACTS["consensus"],
                ),
                "feature_construction": (
                    canonical_features,
                    ARTIFACTS["feature_construction"],
                ),
            }
            if "canonical_v03_root" in workload:
                canonical_roots["score_v03"] = (
                    (project_root / workload["canonical_v03_root"]).resolve(strict=True),
                    ARTIFACTS["score_v03"],
                )
            if "canonical_v04_root" in workload:
                canonical_roots["guard_v04"] = (
                    (project_root / workload["canonical_v04_root"]).resolve(strict=True),
                    ARTIFACTS["guard_v04"],
                )
            canonical_provenance: dict[str, dict[str, dict[str, object]]] = {}
            for source_component, (source_root, names) in canonical_roots.items():
                canonical_provenance[source_component] = {}
                for relative in names:
                    name = Path(relative).name
                    source = source_root / name
                    project_relative = source.relative_to(project_root).as_posix()
                    canonical_provenance[source_component][name] = {
                        "path": project_relative,
                        "bytes": source.stat().st_size,
                        "sha256": sha256_file(source),
                    }
                    input_rows.append(
                        {
                            "workload_id": workload_id,
                            "role": f"canonical_{source_component}_{name}",
                            "bytes": source.stat().st_size,
                            "sha256": sha256_file(source),
                            "path": project_relative,
                        }
                    )
            pool = pool_stats(reference / "consensus/candidate.gff3.manifest.json")
            genome = fasta_stats(resolve_project_file(project_root, workload["genome_fasta"]))
            annotation = gff_stats(resolve_project_file(project_root, workload["annotation_gff"]))
            decisions = decision_stats(reference / "consensus/decisions.tsv")
            method_inputs = {
                method: decision_stats(resolve_project_file(project_root, relative))
                for method, relative in sorted(workload["method_decisions"].items())
            }
            artifacts = {
                relative: {
                    "bytes": (reference / relative).stat().st_size,
                    "sha256": sha256_file(reference / relative),
                }
                for relative in ARTIFACTS["end_to_end"]
            }
            input_paths: dict[str, str] = {
                "genome_fasta": workload["genome_fasta"],
                "annotation_gff": workload["annotation_gff"],
                "base_gff": workload["base_gff"],
                "prior_wgd_selection": workload["prior_wgd_selection"],
                **{f"method_candidate_{key}": value for key, value in workload["method_candidate_gffs"].items()},
                **{f"method_decisions_{key}": value for key, value in workload["method_decisions"].items()},
            }
            inputs: dict[str, dict[str, object]] = {}
            for role, relative in sorted(input_paths.items()):
                path = resolve_project_file(project_root, relative)
                inputs[role] = {"path": relative, "bytes": path.stat().st_size, "sha256": sha256_file(path)}
                input_rows.append(
                    {
                        "workload_id": workload_id,
                        "role": role,
                        "bytes": path.stat().st_size,
                        "sha256": sha256_file(path),
                        "path": relative,
                    }
                )
            registry_workloads.append(
                {
                    "workload_id": workload_id,
                    "species": workload["species"],
                    "assembly": workload["assembly"],
                    "annotation_release": workload["annotation_release"],
                    "split_role": workload["split_role"],
                    "genome": genome,
                    "annotation": annotation,
                    "pool": pool,
                    "decisions": decisions,
                    "method_inputs": method_inputs,
                    "inputs": inputs,
                    "canonical_provenance": canonical_provenance,
                    "reference_artifacts": artifacts,
                }
            )
            timing_rows.append(
                {
                    "workload_id": workload_id,
                    "wall_seconds": format(metrics["wall_seconds"], ".9f"),
                    "user_seconds": format(metrics["user_seconds"], ".9f"),
                    "system_seconds": format(metrics["system_seconds"], ".9f"),
                    "peak_rss_kb": metrics["peak_rss_kb"],
                }
            )
        registry = {
            "schema_version": REFERENCE_SCHEMA,
            "truth_access": False,
            "automatic_approval": False,
            "code": code_manifest,
            "components": list(ARTIFACTS),
            "workloads": registry_workloads,
        }
        (working / "reference_registry.json").write_text(canonical_json(registry), encoding="utf-8")
        write_tsv(
            working / "input_manifest.tsv",
            ["workload_id", "role", "bytes", "sha256", "path"],
            input_rows,
        )
        write_tsv(
            working / "reference_timing.tsv",
            ["workload_id", "wall_seconds", "user_seconds", "system_seconds", "peak_rss_kb"],
            timing_rows,
        )
        write_sha256sums(working)
        working.rename(output)
        freeze_read_only(output)
    except Exception:
        print(f"reference freeze failed; retained attempt at {working}", file=sys.stderr)
        raise
    print(canonical_json({"reference_root": str(output), "schema_version": REFERENCE_SCHEMA}), end="")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except ScalingError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
