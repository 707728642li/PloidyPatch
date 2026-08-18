#!/usr/bin/env python3
"""Freeze one untouched external-holdout contract before target enumeration."""
from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path, PurePosixPath
import re
import shutil
import subprocess
from typing import Any, Iterable

from ploidypatch.artifact_manifest import (
    sha256_file,
    verify_sha256sums,
    write_sha256sums,
)
from ploidypatch.holdout_contract import (
    HoldoutContract,
    load_holdout_contract,
    safe_relative_path,
    staged_relative_path,
)


SCHEMA_VERSION = "ploidypatch.external_holdout_protocol_freeze.v0.5"
STAGE_SCHEMA = "ploidypatch.external_holdout_input_stage.v0.5"
MODEL_SCHEMA = "ploidypatch.composite_ranker.v0.4"
REQUIRED_PIPELINE_ENTRIES = frozenset(
    {"blind_pipeline", "reveal_input_builder", "evaluator"}
)
REQUIRED_BLIND_OUTPUTS = frozenset(
    {"scores", "score_manifest", "pool_decisions", "pool_manifest", "command_log"}
)
GENERIC_IMPLEMENTATION_FILES = (
    "scripts/freeze_external_holdout_protocol_v0.5.py",
    "scripts/freeze_external_holdout_execution_v0.5.py",
    "scripts/run_external_holdout_blind_isolated_v0.5.sh",
    "scripts/finalize_external_holdout_blind_custody_v0.5.py",
    "scripts/run_external_holdout_reveal_v0.5.sh",
    "scripts/build_external_holdout_blind_role_root_v0.5.py",
    "scripts/stage_external_holdout_inputs_v0.5.py",
)
DEFAULT_FORBIDDEN_TEMPLATES = (
    "data/derived/external_evaluator/{holdout_id}",
    "results/evaluator/{holdout_id}",
    "benchmark/structure/{holdout_id}",
    "results/baselines/{holdout_id}",
    "results/copy_collapse/external/{holdout_id}",
    "results/blind_runs/{holdout_id}",
    "results/holdouts/{holdout_id}",
    "work/{holdout_id}",
)
HOLDOUT_FORBIDDEN_PATHS = {
    "actinidia_red5_v0.5": (
        "data/derived/external_evaluator/actinidia_v0.5_wgdi_inputs",
        "data/derived/external_inputs/actinidia/v0.5",
        "results/evaluator/actinidia/v0.5",
        "benchmark/structure/copy_collapse_v0.5/red5_ps1_1.69.0/annotation_copy_collapse_seed20261010",
        "results/baselines/actinidia_v0.5",
        "results/copy_collapse/external/actinidia_v0.5_method_trio",
        "results/copy_collapse/external/actinidia_v0.5_blind_self_wgd",
        "results/copy_collapse/external/actinidia_v0.5_blind_rankings",
    ),
}
HOLDOUT_REQUIRED_PIPELINE_FILES = {
    "actinidia_red5_v0.5": frozenset(
        {
            "scripts/prepare_actinidia_external_normalized_inputs_v0.5.sh",
            "scripts/prepare_actinidia_evaluator_wgdi_inputs_v0.5.sh",
            "scripts/run_actinidia_evaluator_wgdi_v0.5.sh",
            "scripts/infer_actinidia_external_pairs_v0.5.sh",
            "scripts/run_actinidia_copy_collapse_benchmark_v0.5.sh",
            "scripts/build_actinidia_complete_control_reveal_inputs_v0.5.sh",
            "scripts/evaluate_actinidia_external_v0.5.py",
            "scripts/evaluate_external_v0.5.py",
            "scripts/run_actinidia_blind_pipeline_v0.5.sh",
            "scripts/run_actinidia_miniprot_upstream_v0.5.sh",
            "scripts/build_actinidia_method_trio_candidate_pools_v0.5.sh",
            "scripts/run_actinidia_blind_union_self_wgd_v0.5.sh",
            "scripts/score_actinidia_candidates_blind_v0.5.sh",
            "scripts/verify_external_holdout_blind_context_v0.5.py",
            "scripts/run_gemoma_homology.sh",
            "scripts/publish_gemoma_working.sh",
            "scripts/run_lifton_transfer.sh",
            "scripts/normalize_maker_transcript_hierarchy_v0.5.py",
            "scripts/synthesize_missing_transcript_exons.py",
            "scripts/build_wgdi_source_alias_gff.py",
            "scripts/audit_copy_pair_selection_truth.py",
        }
    )
}
HOLDOUT_REQUIRED_PIPELINE_ENTRIES = {
    "actinidia_red5_v0.5": {
        "blind_pipeline": "scripts/run_actinidia_blind_pipeline_v0.5.sh",
        "reveal_input_builder": "scripts/build_actinidia_complete_control_reveal_inputs_v0.5.sh",
        "evaluator": "scripts/evaluate_actinidia_external_v0.5.py",
    }
}
HOLDOUT_REQUIRED_BLIND_OUTPUTS = {
    "actinidia_red5_v0.5": {
        "scores": "results/copy_collapse/external/actinidia_v0.5_blind_rankings/scores/v04.tsv",
        "score_manifest": "results/copy_collapse/external/actinidia_v0.5_blind_rankings/scores/v04.tsv.manifest.json",
        "pool_decisions": "results/copy_collapse/external/actinidia_v0.5_method_trio/consensus/primary_union/blind/decisions.tsv",
        "pool_manifest": "results/copy_collapse/external/actinidia_v0.5_method_trio/consensus/primary_union/blind/candidate.gff3.manifest.json",
        "command_log": "pipeline_commands.tsv",
    }
}
_MAPPING_NAME = re.compile(r"[a-z][a-z0-9_]{1,63}")
_FULL_SHA = re.compile(r"[0-9a-f]{40}")


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
        raise ValueError(f"Expected a JSON object: {path}")
    return value


def parse_mapping(value: str) -> tuple[str, PurePosixPath]:
    name, separator, raw_path = value.partition("=")
    if not separator or _MAPPING_NAME.fullmatch(name) is None:
        raise argparse.ArgumentTypeError("mapping must be NAME=SAFE/RELATIVE/PATH")
    try:
        relative = safe_relative_path(raw_path, f"mapping {name}")
    except ValueError as error:
        raise argparse.ArgumentTypeError(str(error)) from error
    return name, relative


def safe_code_relative(code_root: Path, value: str | Path, context: str) -> str:
    raw = Path(value)
    path = raw.resolve() if raw.is_absolute() else (code_root / raw).resolve()
    try:
        relative = path.relative_to(code_root).as_posix()
    except ValueError:
        raise ValueError(f"{context} must be inside code root: {value}") from None
    safe_relative_path(relative, context)
    if not path.is_file() or path.is_symlink() or path.stat().st_size == 0:
        raise ValueError(f"Missing, empty or symlinked {context}: {path}")
    return relative


def run_checked(command: list[str]) -> bytes:
    return subprocess.run(
        command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    ).stdout


def verify_git_state(code_root: Path, code_commit: str) -> None:
    if _FULL_SHA.fullmatch(code_commit) is None:
        raise ValueError("code_commit must be a full lowercase Git SHA")
    head = run_checked(["git", "-C", str(code_root), "rev-parse", "HEAD"]).decode().strip()
    if head != code_commit:
        raise ValueError(f"code_commit differs from code-root HEAD: {head}")
    if run_checked(["git", "-C", str(code_root), "status", "--porcelain"]):
        raise ValueError("Code root must be completely clean, including untracked files")


def default_forbidden_paths(contract: HoldoutContract) -> tuple[str, ...]:
    generic = tuple(
        template.format(holdout_id=contract.holdout_id)
        for template in DEFAULT_FORBIDDEN_TEMPLATES
    )
    return generic + HOLDOUT_FORBIDDEN_PATHS.get(contract.holdout_id, ())


def reject_target_artifacts(project_root: Path, relatives: Iterable[str]) -> None:
    for raw in sorted(set(relatives)):
        relative = safe_relative_path(raw, "forbidden target-artifact path")
        path = project_root.joinpath(*relative.parts)
        if path.exists() or path.is_symlink():
            raise ValueError(
                "Protocol freeze is too late; target-derived artifact exists: "
                f"{relative.as_posix()}"
            )


def verify_model(model_root: Path, contract: HoldoutContract) -> dict[str, Any]:
    verify_sha256sums(model_root, ignore_checksum_file=True)
    manifest = load_json(model_root / "composite_manifest.json")
    if (
        manifest.get("schema_version") != MODEL_SCHEMA
        or manifest.get("automatic_approval") is not False
        or contract.model_version != "PloidyPatch_ranker_v0.4"
    ):
        raise ValueError("Composite model is not the exact safe v0.4 ranker")
    return manifest


def verify_staged_inputs(
    stage_root: Path, contract: HoldoutContract, contract_path: Path, code_commit: str
) -> dict[str, Any]:
    verify_sha256sums(stage_root, ignore_checksum_file=True)
    role_contract = load_json(stage_root / "role_contract.json")
    if (
        role_contract.get("schema_version") != STAGE_SCHEMA
        or role_contract.get("holdout_id") != contract.holdout_id
        or role_contract.get("policy_id") != contract.policy_id
        or role_contract.get("model_version") != contract.model_version
        or role_contract.get("code_commit") != code_commit
        or role_contract.get("contract", {}).get("sha256") != sha256_file(contract_path)
        or role_contract.get("truth_blind") != dict(contract.truth_blind)
        or role_contract.get("role_boundaries", {}).get(
            "candidate_evaluator_species_overlap"
        )
        is not False
    ):
        raise ValueError("Staged role contract differs from the holdout contract")

    manifest_path = stage_root / "role_manifest.tsv"
    with manifest_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        expected_fields = [
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
        ]
        if reader.fieldnames != expected_fields:
            raise ValueError("Staged role manifest fields differ")
        rows = list(reader)

    expected: dict[tuple[str, str], tuple[str, int, str]] = {}
    for reference in contract.references:
        for artifact_name, artifact in reference.artifact_items():
            expected[(reference.species_id, artifact_name)] = (
                staged_relative_path(reference, artifact_name).as_posix(),
                artifact.staged_bytes,
                artifact.staged_sha256,
            )
    if len(rows) != len(expected):
        raise ValueError("Staged role manifest has the wrong artifact count")
    seen: set[tuple[str, str]] = set()
    for row in rows:
        key = (row["species_id"], row["artifact"])
        if key in seen or key not in expected:
            raise ValueError(f"Unexpected or duplicate staged role row: {key}")
        seen.add(key)
        relative, byte_count, digest = expected[key]
        staged_path = stage_root.joinpath(*PurePosixPath(relative).parts)
        if (
            row["staged_relative_path"] != relative
            or row["bytes"] != str(byte_count)
            or row["sha256"] != digest
            or row["staged_sha256"] != digest
            or not staged_path.is_file()
            or staged_path.is_symlink()
            or staged_path.stat().st_size != byte_count
            or sha256_file(staged_path) != digest
        ):
            raise ValueError(f"Staged role row fails exact binding: {key}")
    return role_contract


def implementation_files(
    code_root: Path,
    contract_path: Path,
    contract: HoldoutContract,
    protocol_artifacts: Iterable[str],
    pipeline_files: Iterable[str],
) -> list[str]:
    core = [
        path.relative_to(code_root).as_posix()
        for path in sorted((code_root / "src/ploidypatch").glob("*.py"))
        if path.is_file()
    ]
    contract_relative = contract_path.relative_to(code_root).as_posix()
    primary_seqids = [
        reference.primary_seqid_table.as_posix() for reference in contract.references
    ]
    return sorted(
        set(
            (*GENERIC_IMPLEMENTATION_FILES, *core, contract_relative,
             *primary_seqids, *protocol_artifacts, *pipeline_files)
        )
    )


def write_tsv(path: Path, header: tuple[str, ...], rows: Iterable[Iterable[Any]]) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def chmod_read_only(root: Path) -> None:
    for path in sorted(root.rglob("*"), reverse=True):
        path.chmod(0o550 if path.is_dir() else 0o440)
    root.chmod(0o550)


def freeze_protocol(
    *,
    project_root: Path,
    code_root: Path,
    contract_path: Path,
    staged_inputs: Path,
    model_root: Path,
    protocol_artifacts: list[str],
    pipeline_entries: dict[str, str],
    pipeline_files: list[str],
    blind_outputs: dict[str, str],
    forbidden_paths: list[str],
    code_commit: str,
    output: Path,
) -> Path:
    partial = Path(str(output) + ".partial")
    if output.exists() or output.is_symlink() or partial.exists() or partial.is_symlink():
        raise FileExistsError("Refusing to overwrite protocol freeze or partial output")
    verify_git_state(code_root, code_commit)
    contract = load_holdout_contract(contract_path)
    if set(pipeline_entries) != REQUIRED_PIPELINE_ENTRIES:
        raise ValueError(
            f"pipeline entries must be exactly {sorted(REQUIRED_PIPELINE_ENTRIES)}"
        )
    expected_entries = HOLDOUT_REQUIRED_PIPELINE_ENTRIES.get(contract.holdout_id)
    if expected_entries is not None and pipeline_entries != expected_entries:
        raise ValueError(
            "Holdout pipeline entries differ from the frozen species entry points"
        )
    if set(blind_outputs) != REQUIRED_BLIND_OUTPUTS:
        raise ValueError(
            f"blind outputs must be exactly {sorted(REQUIRED_BLIND_OUTPUTS)}"
        )
    if len(set(blind_outputs.values())) != len(blind_outputs):
        raise ValueError("Blind output relative paths must be unique")
    expected_outputs = HOLDOUT_REQUIRED_BLIND_OUTPUTS.get(contract.holdout_id)
    if expected_outputs is not None and blind_outputs != expected_outputs:
        raise ValueError(
            "Holdout blind output paths differ from the species pipeline namespace"
        )
    all_pipeline_files = sorted(set((*pipeline_files, *pipeline_entries.values())))
    if not all_pipeline_files:
        raise ValueError("At least one species pipeline file is required")
    absent_required = sorted(
        HOLDOUT_REQUIRED_PIPELINE_FILES.get(contract.holdout_id, frozenset())
        - set(all_pipeline_files)
    )
    if absent_required:
        raise ValueError(
            "Holdout pipeline dependency closure is incomplete: "
            + ", ".join(absent_required)
        )
    for relative in (*protocol_artifacts, *all_pipeline_files):
        safe_code_relative(code_root, relative, "frozen code artifact")
    for relative in blind_outputs.values():
        safe_relative_path(relative, "blind output relative path")

    resolved_forbidden = sorted(
        set((*default_forbidden_paths(contract), *forbidden_paths))
    )
    reject_target_artifacts(project_root, resolved_forbidden)
    verify_model(model_root, contract)
    role_contract = verify_staged_inputs(
        staged_inputs, contract, contract_path, code_commit
    )
    files = implementation_files(
        code_root,
        contract_path,
        contract,
        protocol_artifacts,
        all_pipeline_files,
    )
    missing = [relative for relative in files if not (code_root / relative).is_file()]
    if missing:
        raise ValueError("Missing frozen implementation: " + ", ".join(missing))

    partial.mkdir(parents=True)
    try:
        shutil.copyfile(contract_path, partial / "contract.json")
        shutil.copyfile(staged_inputs / "role_contract.json", partial / "role_contract.json")
        shutil.copyfile(staged_inputs / "role_manifest.tsv", partial / "role_manifest.tsv")
        shutil.copyfile(
            staged_inputs / "role_manifest.tsv",
            partial / "preflight_input_manifest.tsv",
        )
        shutil.copyfile(staged_inputs / "SHA256SUMS", partial / "staged_input_SHA256SUMS")
        protocol_copy_root = partial / "protocol_artifacts"
        for relative in protocol_artifacts:
            destination = protocol_copy_root / relative
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(code_root / relative, destination)

        write_tsv(
            partial / "pipeline_entries.tsv",
            ("entry", "relative_path", "sha256"),
            (
                (name, relative, sha256_file(code_root / relative))
                for name, relative in sorted(pipeline_entries.items())
            ),
        )
        write_tsv(
            partial / "pipeline_files.tsv",
            ("relative_path", "bytes", "sha256"),
            (
                (relative, (code_root / relative).stat().st_size,
                 sha256_file(code_root / relative))
                for relative in all_pipeline_files
            ),
        )
        write_tsv(
            partial / "blind_outputs.tsv",
            ("name", "relative_path"),
            sorted(blind_outputs.items()),
        )
        write_tsv(
            partial / "forbidden_target_artifacts.tsv",
            ("relative_path",),
            ((relative,) for relative in resolved_forbidden),
        )
        write_tsv(
            partial / "implementation_manifest.tsv",
            ("relative_path", "bytes", "sha256"),
            (
                (relative, (code_root / relative).stat().st_size,
                 sha256_file(code_root / relative))
                for relative in files
            ),
        )
        write_tsv(
            partial / "run_contract.tsv",
            ("field", "value"),
            (
                ("schema_version", SCHEMA_VERSION),
                ("holdout_id", contract.holdout_id),
                ("policy_id", contract.policy_id),
                ("test_role", contract.test_role),
                ("model_version", contract.model_version),
                ("code_commit", code_commit),
                ("freeze_stage", "post_metadata_pre_pair_pre_candidate_pre_label"),
                ("wgd_pairs_enumerated_before_freeze", "false"),
                ("candidate_counts_computed_before_freeze", "false"),
                ("truth_labels_accessed_before_freeze", "false"),
                ("automatic_copy_addition_approval", "false"),
            ),
        )
        manifest = {
            "schema_version": SCHEMA_VERSION,
            "holdout_id": contract.holdout_id,
            "policy_id": contract.policy_id,
            "test_role": contract.test_role,
            "model_version": contract.model_version,
            "code_commit": code_commit,
            "freeze_stage": "post_metadata_pre_pair_pre_candidate_pre_label",
            "truth_access": False,
            "wgd_pairs_enumerated": False,
            "candidate_counts_computed": False,
            "truth_labels_accessed": False,
            "contract_sha256": sha256_file(partial / "contract.json"),
            "staged_input_SHA256SUMS_sha256": sha256_file(
                staged_inputs / "SHA256SUMS"
            ),
            "staged_role_contract_sha256": sha256_file(
                staged_inputs / "role_contract.json"
            ),
            "composite_model_SHA256SUMS_sha256": sha256_file(
                model_root / "SHA256SUMS"
            ),
            "pipeline_entries": pipeline_entries,
            "blind_outputs": blind_outputs,
            "forbidden_target_artifacts": resolved_forbidden,
            "implementation_file_count": len(files),
            "staged_role_boundaries": role_contract["role_boundaries"],
        }
        (partial / "protocol_manifest.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        write_sha256sums(partial)
        verify_sha256sums(partial, ignore_checksum_file=True)
        chmod_read_only(partial)
        os.replace(partial, output)
    except BaseException:
        if partial.exists() and not partial.is_symlink():
            shutil.rmtree(partial)
        raise
    return output


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", required=True)
    parser.add_argument("--code-root", required=True)
    parser.add_argument("--contract", required=True)
    parser.add_argument("--staged-inputs", required=True)
    parser.add_argument("--composite-model-freeze", required=True)
    parser.add_argument("--protocol-artifact", action="append", required=True)
    parser.add_argument("--pipeline-entry", action="append", type=parse_mapping, required=True)
    parser.add_argument("--pipeline-file", action="append", default=[])
    parser.add_argument("--blind-output", action="append", type=parse_mapping, required=True)
    parser.add_argument("--forbidden-target-artifact", action="append", default=[])
    parser.add_argument("--code-commit", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args(argv)
    pipeline_entries = dict(args.pipeline_entry)
    blind_outputs = dict(args.blind_output)
    if len(pipeline_entries) != len(args.pipeline_entry):
        raise ValueError("Duplicate pipeline-entry name")
    if len(blind_outputs) != len(args.blind_output):
        raise ValueError("Duplicate blind-output name")
    project_root = Path(args.project_root).resolve()
    code_root = Path(args.code_root).resolve()
    contract_path = Path(args.contract).resolve()
    if not contract_path.is_relative_to(code_root):
        raise ValueError("Contract must be a committed file inside code root")
    freeze_protocol(
        project_root=project_root,
        code_root=code_root,
        contract_path=contract_path,
        staged_inputs=Path(args.staged_inputs).resolve(),
        model_root=Path(args.composite_model_freeze).resolve(),
        protocol_artifacts=[
            safe_code_relative(code_root, value, "protocol artifact")
            for value in args.protocol_artifact
        ],
        pipeline_entries={
            name: safe_code_relative(code_root, relative.as_posix(), "pipeline entry")
            for name, relative in args.pipeline_entry
        },
        pipeline_files=[
            safe_code_relative(code_root, value, "pipeline file")
            for value in args.pipeline_file
        ],
        blind_outputs={name: relative.as_posix() for name, relative in args.blind_output},
        forbidden_paths=[
            safe_relative_path(value, "forbidden target artifact").as_posix()
            for value in args.forbidden_target_artifact
        ],
        code_commit=args.code_commit,
        output=Path(args.output_dir).resolve(),
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
