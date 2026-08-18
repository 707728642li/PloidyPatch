#!/usr/bin/env python3
"""Shared primitives for the frozen PloidyPatch core-scaling benchmark."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import re
import shutil
import stat
import subprocess
import sys
import tarfile
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable


CONFIG_SCHEMA = "ploidypatch.core_scaling_config.v1"
REFERENCE_SCHEMA = "ploidypatch.core_scaling_references.v1"
RESULT_SCHEMA = "ploidypatch.core_scaling_benchmark.v1"
COMPONENTS = (
    "consensus",
    "feature_construction",
    "score_v03",
    "guard_v04",
    "end_to_end",
)
ARTIFACTS = {
    "consensus": (
        "consensus/candidate.gff3",
        "consensus/decisions.tsv",
        "consensus/candidate.gff3.manifest.json",
    ),
    "feature_construction": (
        "features/wgd_selection.tsv",
        "features/wgd_selection.tsv.manifest.json",
        "features/copy_features.tsv",
        "features/copy_features.tsv.manifest.json",
        "features/topology_features.tsv",
        "features/topology_features.tsv.manifest.json",
    ),
    "score_v03": (
        "scores/v03.tsv",
        "scores/v03.tsv.manifest.json",
    ),
    "guard_v04": (
        "scores/v04.tsv",
        "scores/v04.tsv.manifest.json",
    ),
}
ARTIFACTS["end_to_end"] = tuple(
    item
    for component in COMPONENTS[:-1]
    for item in ARTIFACTS[component]
)
TIME_KEYS = {
    "User time (seconds)": "user_seconds",
    "System time (seconds)": "system_seconds",
    "Maximum resident set size (kbytes)": "peak_rss_kb",
    "File system inputs": "filesystem_inputs",
    "File system outputs": "filesystem_outputs",
}
SAFE_ID = re.compile(r"^[a-z0-9][a-z0-9_.-]*$")
HEX40 = re.compile(r"^[0-9a-f]{40}$")


class ScalingError(RuntimeError):
    """Raised for an invalid benchmark contract or divergent execution."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ScalingError(message)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical_json(value: Any) -> str:
    return json.dumps(value, indent=2, sort_keys=True) + "\n"


def load_json(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        value = json.load(handle)
    require(isinstance(value, dict), f"JSON root is not an object: {path}")
    return value


def resolve_project_file(project_root: Path, relative: str) -> Path:
    require(relative and not Path(relative).is_absolute(), f"non-relative path: {relative}")
    root = project_root.resolve()
    path = (root / relative).resolve(strict=True)
    try:
        path.relative_to(root)
    except ValueError as exc:
        raise ScalingError(f"path escapes project root: {relative}") from exc
    require(path.is_file() and path.stat().st_size > 0, f"missing/empty input: {relative}")
    return path


def load_config(config_path: Path, project_root: Path) -> dict[str, Any]:
    config = load_json(config_path)
    require(config.get("schema_version") == CONFIG_SCHEMA, "invalid scaling config schema")
    model = config.get("model_v03")
    require(isinstance(model, str), "model_v03 must be a project-relative path")
    resolve_project_file(project_root, model)
    workloads = config.get("workloads")
    require(isinstance(workloads, list) and workloads, "scaling config has no workloads")
    observed: set[str] = set()
    for workload in workloads:
        require(isinstance(workload, dict), "workload is not an object")
        workload_id = workload.get("workload_id")
        require(
            isinstance(workload_id, str)
            and SAFE_ID.fullmatch(workload_id) is not None
            and workload_id not in observed,
            f"invalid/duplicate workload_id: {workload_id}",
        )
        observed.add(workload_id)
        for key in (
            "species",
            "assembly",
            "annotation_release",
            "split_role",
            "genome_fasta",
            "annotation_gff",
            "base_gff",
            "prior_wgd_selection",
            "canonical_pool_root",
            "canonical_feature_root",
        ):
            require(isinstance(workload.get(key), str) and workload[key], f"{workload_id}: missing {key}")
        for key in (
            "genome_fasta",
            "annotation_gff",
            "base_gff",
            "prior_wgd_selection",
        ):
            resolve_project_file(project_root, workload[key])
        for key, required_names in (
            ("method_candidate_gffs", {"miniprot", "gemoma", "lifton"}),
            ("method_decisions", {"miniprot", "gemoma", "lifton"}),
        ):
            mapping = workload.get(key)
            require(isinstance(mapping, dict) and set(mapping) == required_names, f"{workload_id}: invalid {key}")
            for relative in mapping.values():
                require(isinstance(relative, str), f"{workload_id}: invalid {key} path")
                resolve_project_file(project_root, relative)
        pool_root = Path(workload["canonical_pool_root"])
        feature_root = Path(workload["canonical_feature_root"])
        require(not pool_root.is_absolute(), f"{workload_id}: pool root must be relative")
        require(not feature_root.is_absolute(), f"{workload_id}: feature root must be relative")
        for name in ARTIFACTS["consensus"]:
            resolve_project_file(project_root, str(pool_root / Path(name).name))
        for name in ARTIFACTS["feature_construction"]:
            resolve_project_file(project_root, str(feature_root / Path(name).name))
        for optional in ("canonical_v03_root", "canonical_v04_root"):
            if optional in workload:
                root = Path(workload[optional])
                require(not root.is_absolute(), f"{workload_id}: {optional} must be relative")
                names = ARTIFACTS["score_v03" if optional == "canonical_v03_root" else "guard_v04"]
                for name in names:
                    resolve_project_file(project_root, str(root / Path(name).name))
    return config


def parse_time_v(path: Path) -> dict[str, int | float]:
    values: dict[str, int | float] = {}
    with path.open(encoding="utf-8", errors="strict") as handle:
        for raw in handle:
            line = raw.strip()
            for source, target in TIME_KEYS.items():
                prefix = source + ":"
                if line.startswith(prefix):
                    text = line[len(prefix) :].strip()
                    values[target] = float(text) if "seconds" in target else int(text)
    missing = set(TIME_KEYS.values()) - set(values)
    require(not missing, f"time -v output lacks fields {sorted(missing)}: {path}")
    return values


def fasta_stats(path: Path) -> dict[str, int]:
    opener = gzip.open if path.suffix == ".gz" else open
    sequences = 0
    bases = 0
    with opener(path, "rt", encoding="utf-8") as handle:
        for raw in handle:
            if raw.startswith(">"):
                sequences += 1
            else:
                bases += len(raw.strip())
    require(sequences > 0 and bases > 0, f"empty FASTA: {path}")
    return {"sequences": sequences, "bases": bases}


def gff_stats(path: Path) -> dict[str, int]:
    counts = {"features": 0, "genes": 0, "transcripts": 0}
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as handle:
        for raw in handle:
            if not raw or raw.startswith("#"):
                continue
            fields = raw.rstrip("\r\n").split("\t")
            require(len(fields) == 9, f"malformed GFF row: {path}")
            counts["features"] += 1
            if fields[2] == "gene":
                counts["genes"] += 1
            if fields[2] in {"mRNA", "transcript"}:
                counts["transcripts"] += 1
    require(counts["features"] > 0, f"empty GFF: {path}")
    return counts


def pool_stats(manifest_path: Path) -> dict[str, int]:
    manifest = load_json(manifest_path)
    require(manifest.get("schema_version") == "ploidypatch.method_candidate_pool.v2", f"invalid pool manifest: {manifest_path}")
    counts = manifest.get("counts", {})
    output = {
        "candidate_count": int(counts["accepted_models"]),
        "conflict_sets": int(counts["conflict_sets"]),
        "conflicted_chains": int(counts["conflicted_chains"]),
    }
    for key, value in output.items():
        require(value >= 0, f"negative pool count {key}: {manifest_path}")
    return output


@dataclass(frozen=True)
class CommandSpec:
    label: str
    argv: tuple[str, ...]


def _cli(python_bin: Path, *args: str) -> tuple[str, ...]:
    return (str(python_bin), "-m", "ploidypatch.cli", *args)


def command_specs(
    *,
    component: str,
    workload: dict[str, Any],
    project_root: Path,
    python_bin: Path,
    artifact_root: Path,
    reference_root: Path,
    model_path: Path,
) -> list[CommandSpec]:
    require(component in COMPONENTS, f"unsupported component: {component}")
    workload_id = workload["workload_id"]
    reference = reference_root / "workloads" / workload_id
    artifact_root.mkdir(parents=True, exist_ok=True)
    for directory in ("consensus", "features", "scores"):
        (artifact_root / directory).mkdir(exist_ok=True)
    use_generated_consensus = component in {"consensus", "end_to_end"}
    use_generated_features = component in {"feature_construction", "end_to_end"}
    use_generated_v03 = component in {"score_v03", "end_to_end"}
    pool_candidate = (
        artifact_root / "consensus/candidate.gff3"
        if use_generated_consensus
        else reference / "consensus/candidate.gff3"
    )
    pool_decisions = (
        artifact_root / "consensus/decisions.tsv"
        if use_generated_consensus
        else reference / "consensus/decisions.tsv"
    )
    pool_manifest = (
        artifact_root / "consensus/candidate.gff3.manifest.json"
        if use_generated_consensus
        else reference / "consensus/candidate.gff3.manifest.json"
    )
    feature_root = artifact_root / "features" if use_generated_features else reference / "features"
    v03_path = artifact_root / "scores/v03.tsv" if use_generated_v03 else reference / "scores/v03.tsv"
    base_gff = resolve_project_file(project_root, workload["base_gff"])
    prior_wgd = resolve_project_file(project_root, workload["prior_wgd_selection"])
    candidate_gffs = {
        key: resolve_project_file(project_root, value)
        for key, value in workload["method_candidate_gffs"].items()
    }
    decisions = {
        key: resolve_project_file(project_root, value)
        for key, value in workload["method_decisions"].items()
    }
    specs: list[CommandSpec] = []
    if component in {"consensus", "end_to_end"}:
        specs.append(
            CommandSpec(
                "consensus",
                _cli(
                    python_bin,
                    "baseline",
                    "select-method-consensus",
                    "--base-gff",
                    str(base_gff),
                    *(item for method in ("miniprot", "gemoma", "lifton") for item in ("--candidate", f"{method}={candidate_gffs[method]}")),
                    "--min-method-support",
                    "1",
                    "--max-redundancy-overlap",
                    "0.5",
                    "--redundancy-policy",
                    "retain_distinct_chains",
                    "--output-gff",
                    str(pool_candidate),
                    "--decisions-tsv",
                    str(pool_decisions),
                ),
            )
        )
    if component in {"feature_construction", "end_to_end"}:
        specs.extend(
            [
                CommandSpec(
                    "wgd_propagation",
                    _cli(
                        python_bin,
                        "evidence",
                        "propagate-wgd-conflict-partners",
                        "--base-gff",
                        str(base_gff),
                        "--candidate-gff",
                        str(pool_candidate),
                        "--pool-decisions",
                        str(pool_decisions),
                        "--prior-wgd-selection",
                        str(prior_wgd),
                        "--output-selection",
                        str(artifact_root / "features/wgd_selection.tsv"),
                    ),
                ),
                CommandSpec(
                    "copy_features",
                    _cli(
                        python_bin,
                        "evidence",
                        "build-copy-features",
                        "--consensus-decisions",
                        str(pool_decisions),
                        *(item for method in ("miniprot", "gemoma", "lifton") for item in ("--method-decisions", f"{method}={decisions[method]}")),
                        "--wgd-selection",
                        str(artifact_root / "features/wgd_selection.tsv"),
                        "--output-tsv",
                        str(artifact_root / "features/copy_features.tsv"),
                    ),
                ),
                CommandSpec(
                    "topology_features",
                    _cli(
                        python_bin,
                        "evidence",
                        "build-homeolog-topology-features",
                        "--copy-features",
                        str(artifact_root / "features/copy_features.tsv"),
                        "--wgd-selection",
                        str(artifact_root / "features/wgd_selection.tsv"),
                        "--candidate-gff",
                        str(pool_candidate),
                        "--base-gff",
                        str(base_gff),
                        "--output-tsv",
                        str(artifact_root / "features/topology_features.tsv"),
                    ),
                ),
            ]
        )
    if component in {"score_v03", "end_to_end"}:
        specs.append(
            CommandSpec(
                "score_v03",
                _cli(
                    python_bin,
                    "evidence",
                    "score-support-conditioned-candidates",
                    "--copy-features",
                    str(feature_root / "copy_features.tsv"),
                    "--topology-features",
                    str(feature_root / "topology_features.tsv"),
                    "--model-json",
                    str(model_path),
                    "--output-tsv",
                    str(artifact_root / "scores/v03.tsv"),
                ),
            )
        )
    if component in {"guard_v04", "end_to_end"}:
        specs.append(
            CommandSpec(
                "guard_v04",
                _cli(
                    python_bin,
                    "evidence",
                    "apply-conflict-winner-guard",
                    "--v03-scores",
                    str(v03_path),
                    "--pool-decisions",
                    str(pool_decisions),
                    "--pool-manifest",
                    str(pool_manifest),
                    "--output-tsv",
                    str(artifact_root / "scores/v04.tsv"),
                ),
            )
        )
    return specs


def run_specs(
    specs: Iterable[CommandSpec],
    *,
    code_root: Path,
    log_root: Path,
    env: dict[str, str],
) -> dict[str, Any]:
    require(Path("/usr/bin/time").is_file(), "/usr/bin/time is required")
    log_root.mkdir(parents=True, exist_ok=True)
    commands: list[dict[str, Any]] = []
    overall_start = time.monotonic()
    for index, spec in enumerate(specs, start=1):
        prefix = log_root / f"{index:02d}_{spec.label}"
        time_path = Path(str(prefix) + ".time.txt")
        stdout_path = Path(str(prefix) + ".stdout.json")
        stderr_path = Path(str(prefix) + ".stderr.log")
        command = ("/usr/bin/time", "-v", "-o", str(time_path), *spec.argv)
        started = time.monotonic()
        with stdout_path.open("xb") as stdout_handle, stderr_path.open("xb") as stderr_handle:
            completed = subprocess.run(
                command,
                cwd=code_root,
                env=env,
                stdout=stdout_handle,
                stderr=stderr_handle,
                check=False,
            )
        wall = time.monotonic() - started
        metrics = parse_time_v(time_path)
        commands.append(
            {
                "label": spec.label,
                "argv": list(spec.argv),
                "exit_status": completed.returncode,
                "wall_seconds": wall,
                **metrics,
            }
        )
        if completed.returncode != 0:
            raise ScalingError(f"command failed ({completed.returncode}): {spec.label}; see {stderr_path}")
    return {
        "wall_seconds": time.monotonic() - overall_start,
        "user_seconds": sum(float(item["user_seconds"]) for item in commands),
        "system_seconds": sum(float(item["system_seconds"]) for item in commands),
        "peak_rss_kb": max(int(item["peak_rss_kb"]) for item in commands),
        "filesystem_inputs": sum(int(item["filesystem_inputs"]) for item in commands),
        "filesystem_outputs": sum(int(item["filesystem_outputs"]) for item in commands),
        "commands": commands,
    }


def verify_artifacts(
    *, component: str, artifact_root: Path, reference_workload_root: Path
) -> tuple[int, str]:
    lines: list[str] = []
    total_bytes = 0
    for relative in ARTIFACTS[component]:
        observed = artifact_root / relative
        reference = reference_workload_root / relative
        require(observed.is_file(), f"missing benchmark artifact: {observed}")
        require(reference.is_file(), f"missing reference artifact: {reference}")
        observed_sha = sha256_file(observed)
        reference_sha = sha256_file(reference)
        require(observed_sha == reference_sha, f"output divergence for {relative}")
        total_bytes += observed.stat().st_size
        lines.append(f"{relative}\t{observed_sha}\n")
    digest = hashlib.sha256("".join(lines).encode("utf-8")).hexdigest()
    return total_bytes, digest


def compare_canonical(
    generated_root: Path,
    canonical_root: Path,
    names: Iterable[str],
) -> None:
    for relative in names:
        name = Path(relative).name
        generated = generated_root / relative
        canonical = canonical_root / name
        require(canonical.is_file(), f"missing canonical artifact: {canonical}")
        require(sha256_file(generated) == sha256_file(canonical), f"reference reconstruction differs: {canonical}")


def write_sha256sums(root: Path) -> None:
    output = root / "SHA256SUMS"
    require(not output.exists(), f"checksum manifest already exists: {output}")
    files = sorted(
        path
        for path in root.rglob("*")
        if path.is_file() and path != output
    )
    with output.open("x", encoding="utf-8", newline="") as handle:
        for path in files:
            relative = path.relative_to(root).as_posix()
            handle.write(f"{sha256_file(path)}  {relative}\n")


def verify_sha256sums(root: Path) -> int:
    """Verify both hashes and the exact file universe of a frozen result."""

    manifest = root / "SHA256SUMS"
    require(manifest.is_file(), f"missing SHA256SUMS: {root}")
    expected: dict[str, str] = {}
    with manifest.open(encoding="utf-8") as handle:
        for number, raw in enumerate(handle, start=1):
            fields = raw.rstrip("\n").split("  ", 1)
            require(
                len(fields) == 2 and len(fields[0]) == 64,
                f"malformed SHA256SUMS line {number}",
            )
            relative = fields[1]
            path = Path(relative)
            require(
                not path.is_absolute() and ".." not in path.parts,
                f"unsafe SHA path: {relative}",
            )
            require(relative not in expected, f"duplicate SHA path: {relative}")
            expected[relative] = fields[0]
    observed = {
        path.relative_to(root).as_posix()
        for path in root.rglob("*")
        if path.is_file() and path != manifest
    }
    require(observed == set(expected), "SHA256SUMS file universe differs")
    for relative, digest in expected.items():
        require(
            sha256_file(root / relative) == digest,
            f"checksum differs: {relative}",
        )
    return len(expected)


def freeze_read_only(root: Path) -> None:
    for path in sorted(root.rglob("*"), reverse=True):
        mode = path.stat().st_mode
        if path.is_dir():
            path.chmod((mode & ~0o222) | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
        else:
            path.chmod(mode & ~0o222)
    root.chmod((root.stat().st_mode & ~0o222) | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


def environment_manifest() -> dict[str, Any]:
    def capture(*argv: str) -> str:
        completed = subprocess.run(argv, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False)
        return completed.stdout.strip()

    return {
        "hostname": capture("hostname"),
        "kernel": capture("uname", "-srvmo"),
        "lscpu": capture("lscpu"),
        "memory": capture("free", "-b"),
        "filesystem": capture("df", "-T", "."),
        "python": capture(os.environ.get("PYTHON", "python3"), "--version"),
        "time": capture("/usr/bin/time", "--version").splitlines()[0],
    }


def safely_extract_archive(archive: Path, target: Path) -> None:
    """Extract a regular-file-only tar.gz after rejecting traversal and links."""

    target.mkdir(parents=True)
    with tarfile.open(archive, "r:gz") as handle:
        members = handle.getmembers()
        for member in members:
            path = Path(member.name)
            require(
                not path.is_absolute() and ".." not in path.parts,
                f"unsafe code archive member: {member.name}",
            )
            require(
                not member.issym() and not member.islnk(),
                f"links forbidden in code archive: {member.name}",
            )
            require(
                member.isfile() or member.isdir(),
                f"special archive member forbidden: {member.name}",
            )
        for member in members:
            destination = target / member.name
            if member.isdir():
                destination.mkdir(parents=True, exist_ok=True)
                continue
            destination.parent.mkdir(parents=True, exist_ok=True)
            source = handle.extractfile(member)
            require(source is not None, f"cannot read archive member: {member.name}")
            with source, destination.open("xb") as output:
                shutil.copyfileobj(source, output)


def write_tsv(path: Path, fieldnames: list[str], rows: Iterable[dict[str, Any]]) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def remove_working_on_preexecution_failure(path: Path) -> None:
    """Remove only a new, empty pre-execution directory.

    Once any command log exists, the failed attempt is evidence and is retained.
    """

    if path.is_dir() and not any(path.rglob("*")):
        path.rmdir()


def _validation_main() -> int:
    parser = argparse.ArgumentParser(description="Validate a core-scaling workload config")
    parser.add_argument("--project-root", required=True, type=Path)
    parser.add_argument("--config", required=True, type=Path)
    args = parser.parse_args()
    config = load_config(args.config.resolve(strict=True), args.project_root.resolve(strict=True))
    print(
        canonical_json(
            {
                "schema_version": config["schema_version"],
                "workloads": [item["workload_id"] for item in config["workloads"]],
            }
        ),
        end="",
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(_validation_main())
    except ScalingError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
