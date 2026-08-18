#!/usr/bin/env python3
"""Freeze committed Coffea H1 code and six exact conda environments."""
from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path, PurePosixPath
import re
import shutil
import subprocess
import tarfile
import tempfile
from typing import Iterable

from ploidypatch.artifact_manifest import sha256_file, verify_sha256sums, write_sha256sums
from ploidypatch.coffea_h1_framework import (
    BLIND_OUTPUTS,
    EXECUTION_SCHEMA,
    HOLDOUT_ID,
    PIPELINE_ENTRIES,
    POLICY_ID,
    REQUIRED_ENVIRONMENTS,
    verify_protocol,
)
from ploidypatch.reproducible_projection import verify_tree_manifest


FULL_SHA = re.compile(r"[0-9a-f]{40}")
ENVIRONMENT = re.compile(r"[a-z0-9][a-z0-9_.-]*")
IMPLEMENTATION_PREFIXES = (
    "src/ploidypatch/",
    "scripts/",
    "config/primary_seqids/",
)
IMPLEMENTATION_EXACT = {
    "config/coffea_external_validation_policy_v1.0.tsv",
    "config/coffea_external_event_definition_v1.0.tsv",
    "config/coffea_et39_homoeolog_groups_v1.0.tsv",
    "config/coffea_gardenia_seqid_aliases_v1.0.tsv",
    "config/holdouts/coffea_et39_v1.0/contract.json",
}
PATCH_STAGE = "post_evaluator_truth_failed_blind_pre_candidate_pre_label_execution_patch"
PATCH2_STAGE = (
    "post_evaluator_truth_failed_blind_partial_candidate_pre_label_execution_patch_2"
)
PATCH3_STAGE = (
    "post_evaluator_truth_two_complete_blind_runs_pre_label_"
    "reproducibility_patch_3"
)
PATCH4_STAGE = (
    "post_blind_custody_reveal_authorized_pre_evaluator_environment_patch_4"
)
PATCH5_STAGE = (
    "post_blind_custody_pre_truth_authorization_custody_lineage_patch_5"
)
PATCH1_ALLOWED_PATHS = {
    "docs/COFFEA_EXECUTION_PATCH_MANIFEST_ALIAS_v1.0.md",
    "scripts/freeze_coffea_external_execution_v1.0.py",
    "scripts/run_coffea_blind_isolated_v1.0.sh",
    "src/ploidypatch/coffea_h1_framework.py",
    "tests/test_coffea_h1_execution_v1_0.py",
}
PATCH2_ALLOWED_PATHS = {
    "docs/COFFEA_EXECUTION_PATCH_SINGLE_REFERENCE_SCOPE_v1.0.md",
    "scripts/finalize_coffea_blind_custody_v1.0.py",
    "scripts/freeze_coffea_external_execution_v1.0.py",
    "src/ploidypatch/candidate_merge.py",
    "src/ploidypatch/coffea_h1_framework.py",
    "tests/test_candidate_merge.py",
    "tests/test_coffea_h1_execution_v1_0.py",
}
PATCH3_ALLOWED_PATHS = {
    "docs/COFFEA_BLIND_REPRODUCIBILITY_AMENDMENT_v1.1.md",
    "scripts/finalize_coffea_blind_custody_v1.0.py",
    "scripts/freeze_coffea_external_execution_v1.0.py",
    "scripts/reconcile_coffea_blind_replicates_v1.1.py",
    "scripts/run_coffea_reconciliation_isolated_v1.1.sh",
    "src/ploidypatch/coffea_h1_framework.py",
    "src/ploidypatch/reproducible_projection.py",
    "tests/test_coffea_h1_execution_v1_0.py",
    "tests/test_coffea_reproducibility_amendment_v1_1.py",
    "tests/test_reproducible_projection.py",
}
PATCH4_ALLOWED_PATHS = {
    "docs/COFFEA_REVEAL_EXECUTION_PATCH_v1.2.md",
    "scripts/freeze_coffea_external_execution_v1.0.py",
    "scripts/run_coffea_blind_isolated_v1.0.sh",
    "scripts/run_coffea_evaluator_pipeline_v1.0.sh",
    "scripts/run_coffea_external_reveal_v1.0.py",
    "scripts/run_coffea_reconciliation_isolated_v1.1.sh",
    "src/ploidypatch/coffea_h1_framework.py",
    "tests/test_coffea_h1_execution_v1_0.py",
}
PATCH5_ALLOWED_PATHS = {
    "docs/COFFEA_CUSTODY_LINEAGE_PATCH_v1.3.md",
    "scripts/freeze_coffea_external_execution_v1.0.py",
    "scripts/run_coffea_external_reveal_v1.0.py",
    "src/ploidypatch/coffea_h1_framework.py",
    "tests/test_coffea_h1_execution_v1_0.py",
}


def run(command: list[str], *, cwd: Path | None = None) -> bytes:
    return subprocess.run(
        command, cwd=cwd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    ).stdout


def parse_environment(value: str) -> tuple[str, Path]:
    name, separator, raw = value.partition("=")
    if not separator or ENVIRONMENT.fullmatch(name) is None or not raw:
        raise argparse.ArgumentTypeError("environment must be NAME=PREFIX")
    return name, Path(raw).resolve()


def verify_git(code_root: Path, commit: str) -> None:
    if FULL_SHA.fullmatch(commit) is None:
        raise ValueError("Coffea execution commit must be a full lowercase SHA")
    head = run(["git", "-C", str(code_root), "rev-parse", "HEAD"]).decode().strip()
    if head != commit:
        raise ValueError(f"Coffea execution commit differs from HEAD: {head}")
    if run(["git", "-C", str(code_root), "status", "--porcelain"]):
        raise ValueError("Coffea execution code root must be completely clean")


def safe_extract(archive: Path, destination: Path) -> None:
    with tarfile.open(archive, "r:") as handle:
        members = handle.getmembers()
        for member in members:
            relative = PurePosixPath(member.name)
            if (
                relative.is_absolute()
                or not relative.parts
                or any(part in {"", ".", ".."} for part in relative.parts)
                or member.issym()
                or member.islnk()
                or not (member.isfile() or member.isdir())
            ):
                raise ValueError(f"Unsafe Coffea source archive member: {member.name}")
        handle.extractall(destination, members=members, filter="data")


def environment_locks(prefix: Path) -> tuple[bytes, bytes]:
    if prefix.is_symlink() or not (prefix / "conda-meta/history").is_file():
        raise ValueError(f"Not a real Coffea conda environment: {prefix}")
    explicit = run(["conda", "list", "--explicit", "-p", str(prefix)])
    python = prefix / "bin/python"
    if not python.is_file():
        pip = b"# python unavailable; pip lock not applicable\n"
    else:
        try:
            pip = run([str(python), "-m", "pip", "freeze", "--all"])
        except subprocess.CalledProcessError:
            pip = b"# pip unavailable; explicit conda lock is authoritative\n"
    return explicit, pip


def implementation_paths(code_root: Path, commit: str) -> list[str]:
    paths = run(
        ["git", "-C", str(code_root), "ls-tree", "-r", "--name-only", commit]
    ).decode().splitlines()
    selected = sorted(
        path
        for path in paths
        if path in IMPLEMENTATION_EXACT or path.startswith(IMPLEMENTATION_PREFIXES)
    )
    required = {
        *PIPELINE_ENTRIES.values(),
        "scripts/build_coffea_blind_role_root_v1.0.py",
        "scripts/build_coffea_evaluator_role_root_v1.0.py",
        "scripts/run_coffea_blind_isolated_v1.0.sh",
        "scripts/finalize_coffea_blind_custody_v1.0.py",
        "scripts/freeze_coffea_external_execution_v1.0.py",
        "scripts/run_coffea_external_reveal_v1.0.py",
        "scripts/prepare_coffea_evaluator_wgdi_inputs_v1.0.py",
        "scripts/run_coffea_evaluator_wgdi_v1.0.sh",
        "scripts/infer_coffea_external_pairs_v1.0.py",
        "scripts/build_coffea_structure_holdout_v1.0.py",
        "scripts/run_coffea_evaluator_pipeline_v1.0.sh",
        "scripts/prepare_coffea_blind_candidate_inputs_v1.0.py",
        "scripts/run_coffea_candidate_methods_v1.0.sh",
        "scripts/build_coffea_h1_candidate_pools_v1.0.py",
        "src/ploidypatch/coffea_h1_framework.py",
        "src/ploidypatch/core_h1_evaluation.py",
        "src/ploidypatch/core_h1_pools.py",
        "src/ploidypatch/known_subgenome_h1.py",
        "src/ploidypatch/reproducible_projection.py",
        "scripts/reconcile_coffea_blind_replicates_v1.1.py",
        "scripts/run_coffea_reconciliation_isolated_v1.1.sh",
    }
    missing = required - set(selected)
    if missing:
        raise ValueError(f"Coffea execution archive lacks implementation: {sorted(missing)}")
    if any("ranker" in path.casefold() or "topology" in path.casefold() for path in required):
        raise ValueError("Coffea core H1 cannot freeze ranker/topology entries")
    return selected


def patch_files(
    code_root: Path,
    base_commit: str,
    patch_commit: str,
    *,
    allowed_paths: set[str],
) -> list[dict[str, object]]:
    lines = run(
        [
            "git", "-C", str(code_root), "diff", "--name-status", "--no-renames",
            base_commit, patch_commit, "--",
        ]
    ).decode().splitlines()
    changed: dict[str, str] = {}
    for line in lines:
        status, separator, relative = line.partition("\t")
        if not separator or status not in {"A", "M"} or not relative:
            raise ValueError(f"Unsafe Coffea execution patch change: {line}")
        changed[relative] = status
    if set(changed) != allowed_paths:
        raise ValueError(
            f"Coffea execution patch differs from exact whitelist: {sorted(changed)}"
        )
    rows: list[dict[str, object]] = []
    for relative, status in sorted(changed.items()):
        path = code_root / relative
        if not path.is_file() or path.is_symlink() or path.stat().st_size == 0:
            raise ValueError(f"Missing Coffea execution patch file: {relative}")
        rows.append(
            {
                "status": status,
                "relative_path": relative,
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
        )
    return rows


def failed_attempt_files(
    root: Path, *, patch_sequence: int
) -> tuple[
    int,
    list[dict[str, object]],
    list[dict[str, object]],
    dict[str, str] | None,
]:
    root_required = {
        "bwrap_command.txt",
        "exit_status.txt",
        "mount_manifest.json",
        "namespace_role_validation.json",
        "project/pipeline_commands.tsv",
        "stderr.log",
        "stdout.log",
    }
    failure_root = (
        "project/results/copy_collapse/external/"
        "coffea_v1.0_h1.invalid_run"
    )
    patch2_small_evidence = {
        f"{failure_root}/raw_predictions.manifest.json",
        f"{failure_root}/combined/retain_distinct/blind/candidate.gff3.manifest.json",
        f"{failure_root}/combined/suppress_overlap/blind/candidate.gff3.manifest.json",
        f"{failure_root}/scopes/bua_only/methods/miniprot/stderr.log",
        f"{failure_root}/scopes/bua_only/methods/miniprot/stdout.json",
    }
    replay_paths = {
        "combined_retain_distinct_pool": (
            f"{failure_root}/combined/retain_distinct/blind/candidate.gff3"
        ),
        "combined_retain_distinct_decisions": (
            f"{failure_root}/combined/retain_distinct/blind/decisions.tsv"
        ),
        "combined_suppress_overlap_pool": (
            f"{failure_root}/combined/suppress_overlap/blind/candidate.gff3"
        ),
        "combined_suppress_overlap_decisions": (
            f"{failure_root}/combined/suppress_overlap/blind/decisions.tsv"
        ),
    }
    completed_root = (
        "project/results/copy_collapse/external/coffea_v1.0_h1"
    )
    patch3_outputs = {
        f"{completed_root}/raw_predictions.manifest.json",
        *{
            f"{completed_root}/{scope}/{arm}/blind/{name}"
            for scope in ("combined", "bua_only", "mauritiana_only")
            for arm in ("retain_distinct", "suppress_overlap")
            for name in (
                "candidate.gff3",
                "decisions.tsv",
                "candidate.gff3.manifest.json",
            )
        },
    }
    if not root.is_dir() or root.is_symlink():
        raise ValueError("Coffea execution patch requires a retained failed attempt")
    files = sorted(path for path in root.rglob("*") if path.is_file())
    if any(path.is_symlink() for path in root.rglob("*")):
        raise ValueError("Symlink in retained Coffea failed attempt")
    relatives = {path.relative_to(root).as_posix() for path in files}
    if patch_sequence == 4:
        if relatives != {"reveal_authorization.json"}:
            raise ValueError("Coffea patch-4 reveal failure universe differs")
        authorization = json.loads(
            (root / "reveal_authorization.json").read_text(encoding="utf-8")
        )
        if (
            authorization.get("schema_version")
            != "ploidypatch.coffea_h1_reveal_authorization.v1.0"
            or authorization.get("holdout_id") != HOLDOUT_ID
            or authorization.get("truth_reveal_authorized") is not True
            or authorization.get("ranker_or_model_authorized") is not False
            or authorization.get("h2_or_topology_ranking_authorized") is not False
        ):
            raise ValueError("Coffea patch-4 reveal authorization differs")
        rows = [
            {
                "relative_path": "reveal_authorization.json",
                "bytes": files[0].stat().st_size,
                "sha256": sha256_file(files[0]),
            }
        ]
        return 1, rows, rows, None
    if patch_sequence == 5:
        if relatives != {"error.log"}:
            raise ValueError("Coffea patch-5 reveal failure universe differs")
        error_log = root / "error.log"
        if "Coffea custody failed pre-truth authorization" not in error_log.read_text(
            encoding="utf-8"
        ):
            raise ValueError("Coffea patch-5 reveal failure differs")
        rows = [
            {
                "relative_path": "error.log",
                "bytes": error_log.stat().st_size,
                "sha256": sha256_file(error_log),
            }
        ]
        return 1, rows, rows, None
    if patch_sequence == 1:
        if relatives != root_required or any(
            relative.startswith("project/results/") for relative in relatives
        ):
            raise ValueError(
                f"Unexpected Coffea failed-attempt universe: {sorted(relatives)}"
            )
        evidence_relatives = root_required
    elif patch_sequence == 2:
        required = root_required | patch2_small_evidence | set(replay_paths.values())
        if not required <= relatives:
            raise ValueError(
                "Coffea patch-2 failed attempt lacks exact evidence: "
                + ", ".join(sorted(required - relatives))
            )
        forbidden = ("/truth", "/labels", "evaluator_only", "target_complete")
        if any(any(token in relative.casefold() for token in forbidden) for relative in relatives):
            raise ValueError("Coffea patch-2 failed tree contains forbidden-role paths")
        evidence_relatives = root_required | patch2_small_evidence
    elif patch_sequence == 3:
        required = root_required | patch3_outputs
        if not required <= relatives:
            raise ValueError(
                "Coffea patch-3 blind run lacks complete outputs: "
                + ", ".join(sorted(required - relatives))
            )
        if "custody_manifest.json" in relatives:
            raise ValueError("Coffea patch-3 input unexpectedly completed custody")
        forbidden = ("/truth", "/labels", "evaluator_only", "target_complete")
        if any(
            any(token in relative.casefold() for token in forbidden)
            for relative in relatives
        ):
            raise ValueError("Coffea patch-3 blind tree contains forbidden-role paths")
        evidence_relatives = root_required | {
            f"{completed_root}/raw_predictions.manifest.json",
            *{
                f"{completed_root}/{scope}/{arm}/blind/candidate.gff3.manifest.json"
                for scope in ("combined", "bua_only", "mauritiana_only")
                for arm in ("retain_distinct", "suppress_overlap")
            },
        }
    else:
        raise ValueError(f"Unsupported Coffea execution patch sequence: {patch_sequence}")
    try:
        exit_status = int((root / "exit_status.txt").read_text(encoding="utf-8").strip())
    except ValueError as error:
        raise ValueError("Malformed Coffea failed-attempt exit status") from error
    if exit_status == 0 and patch_sequence != 3:
        raise ValueError("Coffea execution patch requires a nonzero failed attempt")
    if exit_status != 0 and patch_sequence == 3:
        raise ValueError("Coffea patch-3 requires a completed blind pipeline")
    namespace = json.loads((root / "namespace_role_validation.json").read_text())
    mount = json.loads((root / "mount_manifest.json").read_text())
    if (
        namespace.get("schema_version")
        != "ploidypatch.coffea_core_h1_namespace_validation.v1.0"
        or mount.get("schema_version")
        != "ploidypatch.coffea_core_h1_mount_manifest.v1.0"
        or any(
            namespace.get(key) is not False
            for key in (
                "complete_target_annotation_visible", "evaluator_only_visible",
                "truth_visible", "labels_visible", "nas_data_visible",
                "model_visible", "ranker_visible",
            )
        )
        or mount.get("network_access") is not False
        or mount.get("ranker_or_model_mounted") is not False
    ):
        raise ValueError("Failed Coffea attempt did not preserve blind isolation")
    rows = [
        {
            "relative_path": path.relative_to(root).as_posix(),
            "bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        }
        for path in files
    ]
    indexed = {str(row["relative_path"]): row for row in rows}
    evidence_rows = [indexed[relative] for relative in sorted(evidence_relatives)]
    replay = (
        {
            name: str(indexed[relative]["sha256"])
            for name, relative in sorted(replay_paths.items())
        }
        if patch_sequence == 2
        else None
    )
    return exit_status, rows, evidence_rows, replay


def freeze(
    *, code_root: Path, code_commit: str, protocol: Path,
    environments: dict[str, Path], absence_paths: list[Path], output: Path,
    superseded_execution: Path | None = None,
    failed_attempt: Path | None = None,
    patch_reason: Path | None = None,
    reproducibility_run_a: Path | None = None,
    failed_attempt_log: Path | None = None,
    completed_blind_run: Path | None = None,
) -> Path:
    verify_git(code_root, code_commit)
    protocol_manifest, contract = verify_protocol(protocol)
    if set(environments) != REQUIRED_ENVIRONMENTS:
        raise ValueError(
            f"Coffea execution environments differ: {sorted(environments)}"
        )
    patch_values = (superseded_execution, failed_attempt, patch_reason)
    patch_mode = all(value is not None for value in patch_values)
    if any(value is not None for value in patch_values) and not patch_mode:
        raise ValueError(
            "Coffea execution patch requires superseded execution, failed attempt, and reason"
        )
    if not patch_mode and (
        reproducibility_run_a is not None
        or failed_attempt_log is not None
        or completed_blind_run is not None
    ):
        raise ValueError("Coffea run-A/log bindings require execution patch mode")
    prior_manifest: dict[str, object] | None = None
    changed_rows: list[dict[str, object]] = []
    failed_status: int | None = None
    failed_rows: list[dict[str, object]] = []
    failed_evidence_rows: list[dict[str, object]] = []
    primary_replay: dict[str, str] | None = None
    patch_sequence = 0
    patch_stage = "post_metadata_pre_pair_pre_candidate_pre_label"
    if patch_mode:
        assert superseded_execution is not None and failed_attempt is not None
        assert patch_reason is not None
        from ploidypatch.coffea_h1_framework import verify_execution

        prior_manifest, _, _ = verify_execution(superseded_execution, protocol)
        prior_patch = prior_manifest.get("execution_patch")
        if prior_patch is None:
            patch_sequence = 1
            patch_stage = PATCH_STAGE
            allowed_paths = PATCH1_ALLOWED_PATHS
        elif (
            isinstance(prior_patch, dict)
            and prior_patch.get("schema_version")
            == "ploidypatch.coffea_core_h1_execution_patch.v1.0"
            and prior_patch.get("patch_sequence", 1) == 1
            and prior_patch.get("freeze_stage") == PATCH_STAGE
        ):
            patch_sequence = 2
            patch_stage = PATCH2_STAGE
            allowed_paths = PATCH2_ALLOWED_PATHS
        elif (
            isinstance(prior_patch, dict)
            and prior_patch.get("schema_version")
            == "ploidypatch.coffea_core_h1_execution_patch.v1.0"
            and prior_patch.get("patch_sequence") == 2
            and prior_patch.get("freeze_stage") == PATCH2_STAGE
        ):
            patch_sequence = 3
            patch_stage = PATCH3_STAGE
            allowed_paths = PATCH3_ALLOWED_PATHS
        elif (
            isinstance(prior_patch, dict)
            and prior_patch.get("schema_version")
            == "ploidypatch.coffea_core_h1_execution_patch.v1.0"
            and prior_patch.get("patch_sequence") == 3
            and prior_patch.get("freeze_stage") == PATCH3_STAGE
        ):
            patch_sequence = 4
            patch_stage = PATCH4_STAGE
            allowed_paths = PATCH4_ALLOWED_PATHS
        elif (
            isinstance(prior_patch, dict)
            and prior_patch.get("schema_version")
            == "ploidypatch.coffea_core_h1_execution_patch.v1.0"
            and prior_patch.get("patch_sequence") == 4
            and prior_patch.get("freeze_stage") == PATCH4_STAGE
        ):
            patch_sequence = 5
            patch_stage = PATCH5_STAGE
            allowed_paths = PATCH5_ALLOWED_PATHS
        else:
            raise ValueError(
                "Coffea execution patch chain is not original/patch1/patch2/patch3/patch4"
            )
        base_commit = str(prior_manifest["code_commit"])
        changed_rows = patch_files(
            code_root, base_commit, code_commit, allowed_paths=allowed_paths
        )
        (
            failed_status,
            failed_rows,
            failed_evidence_rows,
            primary_replay,
        ) = failed_attempt_files(failed_attempt, patch_sequence=patch_sequence)
        if not patch_reason.is_file() or patch_reason.is_symlink() or patch_reason.stat().st_size == 0:
            raise ValueError("Missing Coffea execution patch reason")
        if patch_sequence == 3:
            if (
                reproducibility_run_a is None
                or failed_attempt_log is None
                or completed_blind_run is not None
            ):
                raise ValueError(
                    "Coffea patch-3 requires run A and the custody-failure log"
                )
            verify_tree_manifest(
                root=reproducibility_run_a,
                manifest=(
                    superseded_execution / "superseded_failed_attempt_manifest.tsv"
                ),
            )
            if (
                not failed_attempt_log.is_file()
                or failed_attempt_log.is_symlink()
                or failed_attempt_log.stat().st_size == 0
                or "Coffea patch-2 changed a primary combined blind output"
                not in failed_attempt_log.read_text(encoding="utf-8")
            ):
                raise ValueError("Coffea patch-3 custody-failure log differs")
        elif patch_sequence == 4:
            if reproducibility_run_a is not None:
                raise ValueError("Coffea run-A binding is patch-3 only")
            if failed_attempt_log is None or completed_blind_run is None:
                raise ValueError(
                    "Coffea patch-4 requires the reveal failure log and blind custody"
                )
            if (
                not failed_attempt_log.is_file()
                or failed_attempt_log.is_symlink()
                or "Coffea frozen dev Python is missing"
                not in failed_attempt_log.read_text(encoding="utf-8")
            ):
                raise ValueError("Coffea patch-4 reveal failure log differs")
            verify_sha256sums(completed_blind_run, ignore_checksum_file=True)
            custody = json.loads(
                (completed_blind_run / "custody_manifest.json").read_text(
                    encoding="utf-8"
                )
            )
            if (
                custody.get("schema_version")
                != "ploidypatch.coffea_core_h1_blind_custody.v1.0"
                or custody.get("execution_SHA256SUMS_sha256")
                != sha256_file(superseded_execution / "SHA256SUMS")
                or custody.get("truth_mounted") is not False
                or custody.get("automatic_approval") is not False
            ):
                raise ValueError("Coffea patch-4 blind custody differs")
        elif patch_sequence == 5:
            if reproducibility_run_a is not None or failed_attempt_log is not None:
                raise ValueError("Coffea patch-5 accepts only completed blind custody")
            if completed_blind_run is None:
                raise ValueError("Coffea patch-5 requires completed blind custody")
            verify_sha256sums(completed_blind_run, ignore_checksum_file=True)
            custody = json.loads(
                (completed_blind_run / "custody_manifest.json").read_text(
                    encoding="utf-8"
                )
            )
            prior_patch = prior_manifest.get("execution_patch")
            if (
                not isinstance(prior_patch, dict)
                or custody.get("schema_version")
                != "ploidypatch.coffea_core_h1_blind_custody.v1.0"
                or custody.get("execution_SHA256SUMS_sha256")
                != prior_patch.get("superseded_execution_SHA256SUMS_sha256")
                or custody.get("truth_mounted") is not False
                or custody.get("automatic_approval") is not False
            ):
                raise ValueError("Coffea patch-5 custody lineage differs")
        elif (
            reproducibility_run_a is not None
            or failed_attempt_log is not None
            or completed_blind_run is not None
        ):
            raise ValueError("Coffea special execution bindings have the wrong patch")
    for path in absence_paths:
        if path.exists() or path.is_symlink():
            raise ValueError(f"Coffea freeze is too late; artifact exists: {path}")
    if output.exists() or output.is_symlink():
        raise FileExistsError(f"Refusing to overwrite Coffea execution freeze: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    working = Path(tempfile.mkdtemp(prefix=f".{output.name}.working.", dir=output.parent))
    try:
        archive = working / "source.tar"
        with archive.open("xb") as handle:
            subprocess.run(
                ["git", "-C", str(code_root), "archive", "--format=tar", code_commit],
                check=True, stdout=handle,
            )
        source = working / "source"; source.mkdir()
        safe_extract(archive, source)
        paths = implementation_paths(code_root, code_commit)
        with (working / "implementation_manifest.tsv").open(
            "x", encoding="utf-8", newline=""
        ) as handle:
            writer = csv.DictWriter(
                handle, fieldnames=("relative_path", "bytes", "sha256"),
                delimiter="\t", lineterminator="\n",
            )
            writer.writeheader()
            for relative in paths:
                path = source / relative
                if not path.is_file() or path.is_symlink() or path.stat().st_size == 0:
                    raise ValueError(f"Coffea frozen implementation differs: {relative}")
                writer.writerow(
                    {
                        "relative_path": relative,
                        "bytes": path.stat().st_size,
                        "sha256": sha256_file(path),
                    }
                )
        locks = working / "environment_locks"; locks.mkdir()
        bindings: list[dict[str, object]] = []
        for name, prefix in sorted(environments.items()):
            explicit, pip = environment_locks(prefix)
            explicit_path = locks / f"{name}.explicit.txt"
            pip_path = locks / f"{name}.pip-freeze.txt"
            explicit_path.write_bytes(explicit); pip_path.write_bytes(pip)
            bindings.append(
                {
                    "name": name,
                    "host_prefix": str(prefix),
                    "explicit_relative_path": explicit_path.relative_to(working).as_posix(),
                    "explicit_sha256": sha256_file(explicit_path),
                    "pip_relative_path": pip_path.relative_to(working).as_posix(),
                    "pip_sha256": sha256_file(pip_path),
                }
            )
        with (working / "environment_bindings.tsv").open(
            "x", encoding="utf-8", newline=""
        ) as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=(
                    "name", "host_prefix", "explicit_relative_path", "explicit_sha256",
                    "pip_relative_path", "pip_sha256",
                ),
                delimiter="\t", lineterminator="\n",
            )
            writer.writeheader(); writer.writerows(bindings)
        if patch_mode and prior_manifest is not None and bindings != prior_manifest.get("environments"):
            raise ValueError("Coffea execution patch changed an environment binding or lock")
        if patch_mode:
            assert superseded_execution is not None and failed_attempt is not None
            assert patch_reason is not None and failed_status is not None
            with (working / "patch_changed_files.tsv").open(
                "x", encoding="utf-8", newline=""
            ) as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=("status", "relative_path", "bytes", "sha256"),
                    delimiter="\t", lineterminator="\n",
                )
                writer.writeheader(); writer.writerows(changed_rows)
            with (working / "superseded_failed_attempt_manifest.tsv").open(
                "x", encoding="utf-8", newline=""
            ) as handle:
                writer = csv.DictWriter(
                    handle, fieldnames=("relative_path", "bytes", "sha256"),
                    delimiter="\t", lineterminator="\n",
                )
                writer.writeheader(); writer.writerows(failed_rows)
            evidence = working / "failed_attempt_evidence"; evidence.mkdir()
            for row in failed_evidence_rows:
                relative = str(row["relative_path"])
                destination = evidence / relative
                destination.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(failed_attempt / relative, destination)
            shutil.copyfile(patch_reason, working / "patch_reason.md")
            if patch_sequence == 3:
                assert failed_attempt_log is not None
                shutil.copyfile(
                    superseded_execution / "superseded_failed_attempt_manifest.tsv",
                    working / "reproducibility_run_a_manifest.tsv",
                )
                shutil.copyfile(
                    failed_attempt_log, working / "failed_attempt_custody_error.log"
                )
            elif patch_sequence == 4:
                assert failed_attempt_log is not None
                assert completed_blind_run is not None
                shutil.copyfile(
                    failed_attempt_log, working / "failed_reveal_error.log"
                )
                shutil.copyfile(
                    completed_blind_run / "custody_manifest.json",
                    working / "completed_blind_custody_manifest.json",
                )
                shutil.copyfile(
                    completed_blind_run / "SHA256SUMS",
                    working / "completed_blind_root_SHA256SUMS",
                )
        manifest = {
            "schema_version": EXECUTION_SCHEMA,
            "holdout_id": HOLDOUT_ID,
            "policy_id": POLICY_ID,
            "protocol_profile": "core_H1_known_subgenome_no_ranker",
            "freeze_stage": patch_stage,
            "code_commit": code_commit,
            "protocol_code_commit": protocol_manifest["code_commit"],
            "protocol_SHA256SUMS_sha256": sha256_file(protocol / "SHA256SUMS"),
            "contract_sha256": sha256_file(protocol / "contract.json"),
            "source_archive_sha256": sha256_file(archive),
            "implementation_manifest_sha256": sha256_file(
                working / "implementation_manifest.tsv"
            ),
            "pipeline_entries": PIPELINE_ENTRIES,
            "blind_outputs": BLIND_OUTPUTS,
            "environments": bindings,
            "truth_pairs_enumerated_before_execution_freeze": patch_mode,
            "candidate_counts_computed_before_execution_freeze": patch_sequence >= 2,
            "truth_labels_accessed_before_execution_freeze": False,
            "ranker_or_model_execution": False,
            "h2_or_topology_ranking_enabled": False,
            "network_access_in_blind_runner": False,
            "nas_data_mount_in_blind_runner": False,
            "complete_target_annotation_mount_in_blind_runner": False,
            "evaluator_only_mount_in_blind_runner": False,
            "truth_or_label_mount_in_blind_runner": False,
            "all_arm_collateral_loss_maximum": 0,
            "automatic_approval": False,
        }
        if patch_mode:
            assert prior_manifest is not None and superseded_execution is not None
            assert failed_status is not None
            execution_patch = {
                "schema_version": "ploidypatch.coffea_core_h1_execution_patch.v1.0",
                "patch_sequence": patch_sequence,
                "freeze_stage": patch_stage,
                "base_code_commit": prior_manifest["code_commit"],
                "patch_code_commit": code_commit,
                "superseded_execution_SHA256SUMS_sha256": sha256_file(
                    superseded_execution / "SHA256SUMS"
                ),
                "changed_files_manifest_sha256": sha256_file(
                    working / "patch_changed_files.tsv"
                ),
                "failed_attempt_exit_status": failed_status,
                "failed_attempt_manifest_sha256": sha256_file(
                    working / "superseded_failed_attempt_manifest.tsv"
                ),
                "failed_attempt_files": len(failed_rows),
                "failed_attempt_evidence_files": len(failed_evidence_rows),
                "patch_reason_sha256": sha256_file(working / "patch_reason.md"),
                "evaluator_truth_construction_completed_before_patch": True,
                "candidate_generation_completed_before_patch": patch_sequence >= 3,
                "partial_candidate_generation_before_patch": patch_sequence == 2,
                "formal_blind_outputs_frozen_before_patch": patch_sequence >= 4,
                "blind_custody_completed_before_patch": patch_sequence >= 4,
                "formal_scores_generated_before_patch": False,
                "truth_labels_accessed_before_patch": False,
                "scientific_rules_or_thresholds_changed": patch_sequence == 3,
            }
            if patch_sequence == 2:
                if primary_replay is None:
                    raise AssertionError("Coffea patch-2 replay binding is missing")
                execution_patch["primary_combined_replay_outputs"] = primary_replay
                prior_patch = prior_manifest.get("execution_patch")
                execution_patch["superseded_execution_patch"] = prior_patch
            elif patch_sequence == 3:
                execution_patch.update(
                    {
                        "two_blind_candidate_executions_completed_before_patch": True,
                        "projection_reproducibility_abstention_rule_added": True,
                        "biological_rules_or_performance_thresholds_changed": False,
                        "label_informed_selection": False,
                        "unstable_projection_policy": "abstain",
                        "reproducibility_run_a_manifest_relative_path": (
                            "reproducibility_run_a_manifest.tsv"
                        ),
                        "reproducibility_run_a_manifest_sha256": sha256_file(
                            working / "reproducibility_run_a_manifest.tsv"
                        ),
                        "reproducibility_run_b_manifest_relative_path": (
                            "superseded_failed_attempt_manifest.tsv"
                        ),
                        "reproducibility_run_b_manifest_sha256": sha256_file(
                            working / "superseded_failed_attempt_manifest.tsv"
                        ),
                        "failed_attempt_custody_error_log_sha256": sha256_file(
                            working / "failed_attempt_custody_error.log"
                        ),
                        "superseded_execution_patch": prior_manifest.get(
                            "execution_patch"
                        ),
                    }
                )
            elif patch_sequence == 4:
                assert completed_blind_run is not None
                execution_patch.update(
                    {
                        "truth_reveal_authorized_before_patch": True,
                        "evaluator_truth_bytes_hashed_before_patch": True,
                        "evaluator_invoked_before_patch": False,
                        "performance_metrics_computed_before_patch": False,
                        "environment_interpreter_symlink_validation_fixed": True,
                        "canonical_nested_manifest_writer_fixed": True,
                        "blind_custody_manifest_sha256": sha256_file(
                            working / "completed_blind_custody_manifest.json"
                        ),
                        "blind_root_SHA256SUMS_sha256": sha256_file(
                            working / "completed_blind_root_SHA256SUMS"
                        ),
                        "failed_reveal_error_log_sha256": sha256_file(
                            working / "failed_reveal_error.log"
                        ),
                        "superseded_execution_patch": prior_manifest.get(
                            "execution_patch"
                        ),
                    }
                )
            elif patch_sequence == 5:
                assert completed_blind_run is not None
                custody = json.loads(
                    (completed_blind_run / "custody_manifest.json").read_text(
                        encoding="utf-8"
                    )
                )
                execution_patch.update(
                    {
                        "truth_reveal_authorized_before_patch": False,
                        "evaluator_truth_bytes_hashed_before_patch": False,
                        "evaluator_invoked_before_patch": False,
                        "performance_metrics_computed_before_patch": False,
                        "custody_execution_chain_validation_fixed": True,
                        "blind_custody_manifest_sha256": sha256_file(
                            completed_blind_run / "custody_manifest.json"
                        ),
                        "blind_custody_execution_SHA256SUMS_sha256": custody[
                            "execution_SHA256SUMS_sha256"
                        ],
                        "blind_root_SHA256SUMS_sha256": sha256_file(
                            completed_blind_run / "SHA256SUMS"
                        ),
                        "failed_reveal_error_log_sha256": sha256_file(
                            working / "failed_attempt_evidence/error.log"
                        ),
                        "superseded_execution_patch": prior_manifest.get(
                            "execution_patch"
                        ),
                    }
                )
            manifest["execution_patch"] = execution_patch
        (working / "execution_manifest.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        write_sha256sums(working)
        verify_sha256sums(working, ignore_checksum_file=True)
        if os.name != "nt":
            for path in sorted(working.rglob("*"), reverse=True):
                path.chmod(0o550 if path.is_dir() else 0o440)
            working.chmod(0o550)
        os.replace(working, output)
    except BaseException:
        shutil.rmtree(working, ignore_errors=True)
        raise
    return output


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--code-root", required=True)
    parser.add_argument("--code-commit", required=True)
    parser.add_argument("--protocol-freeze", required=True)
    parser.add_argument("--environment", action="append", type=parse_environment, default=[])
    parser.add_argument("--absence-path", action="append", default=[])
    parser.add_argument("--superseded-execution")
    parser.add_argument("--failed-attempt")
    parser.add_argument("--patch-reason-file")
    parser.add_argument("--reproducibility-run-a")
    parser.add_argument("--failed-attempt-log")
    parser.add_argument("--completed-blind-run")
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args(argv)
    environments = dict(args.environment)
    if len(environments) != len(args.environment):
        raise ValueError("Duplicate Coffea environment binding")
    freeze(
        code_root=Path(args.code_root).resolve(),
        code_commit=args.code_commit,
        protocol=Path(args.protocol_freeze).resolve(),
        environments=environments,
        absence_paths=[Path(value).resolve() for value in args.absence_path],
        output=Path(args.output_dir).resolve(),
        superseded_execution=(
            Path(args.superseded_execution).resolve() if args.superseded_execution else None
        ),
        failed_attempt=(Path(args.failed_attempt).resolve() if args.failed_attempt else None),
        patch_reason=(Path(args.patch_reason_file).resolve() if args.patch_reason_file else None),
        reproducibility_run_a=(
            Path(args.reproducibility_run_a).resolve()
            if args.reproducibility_run_a
            else None
        ),
        failed_attempt_log=(
            Path(args.failed_attempt_log).resolve() if args.failed_attempt_log else None
        ),
        completed_blind_run=(
            Path(args.completed_blind_run).resolve()
            if args.completed_blind_run
            else None
        ),
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
