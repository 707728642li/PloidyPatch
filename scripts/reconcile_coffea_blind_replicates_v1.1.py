#!/usr/bin/env python3
"""Build truth-blind Coffea pools from projections reproduced in two blind runs."""
from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path
import shutil
import sys
import tempfile
from typing import Iterable

from ploidypatch.artifact_manifest import (
    sha256_file,
    verify_sha256sums,
    write_sha256sums,
)
from ploidypatch.core_h1_pools import (
    assemble_method_candidate_gff,
    build_raw_predictions_manifest,
    collapse_duplicate_exact_chains,
    run_command,
    seal_blind_pool_manifest,
)
from ploidypatch.coffea_h1_framework import verify_execution
from ploidypatch.known_subgenome_h1 import HOLDOUT_ID, POLICY_ID
from ploidypatch.reproducible_projection import (
    compare_decision_tables,
    filter_candidate_gff_by_upstream_models,
    write_comparison_audit,
)


METHODS = ("miniprot", "gemoma", "lifton")
BUNDLES = ("candidate_bua", "candidate_mauritiana")
SCOPES = {
    "combined": BUNDLES,
    "bua_only": ("candidate_bua",),
    "mauritiana_only": ("candidate_mauritiana",),
}
ARMS = {
    "retain_distinct": (
        "retain_distinct_chains",
        "retain_distinct_phased_CDS_chains",
    ),
    "suppress_overlap": (
        "suppress_overlapping",
        "suppress_strongly_overlapping_alternative_chains",
    ),
}


def _formal_input_hashes(
    *, role_root: Path, benchmark: Path, protocol: Path, execution: Path,
    project_root: Path
) -> dict[str, str]:
    with (role_root / "role_manifest.tsv").open(
        encoding="utf-8", newline=""
    ) as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    species = {
        "candidate_bua": "Coffea_eugenioides_BuA",
        "candidate_mauritiana": "Coffea_mauritiana",
    }
    candidate: dict[str, str] = {}
    for bundle, species_id in species.items():
        for artifact in ("genome", "gff3", "protein"):
            matches = [
                row
                for row in rows
                if row.get("role") == "candidate_reference"
                and row.get("species_id") == species_id
                and row.get("artifact") == artifact
            ]
            if len(matches) != 1:
                raise ValueError(f"Non-unique formal role: {bundle}/{artifact}")
            candidate[f"{bundle}_{artifact}_sha256"] = matches[0]["sha256"]
    role = json.loads((role_root / "role_manifest.json").read_text(encoding="utf-8"))
    normalized = project_root / "results/baselines/coffea/v1.0/normalized"
    return {
        "staged_input_SHA256SUMS_sha256": role[
            "staged_input_SHA256SUMS_sha256"
        ],
        "blind_benchmark_SHA256SUMS_sha256": sha256_file(
            benchmark / "SHA256SUMS"
        ),
        "protocol_SHA256SUMS_sha256": sha256_file(protocol / "SHA256SUMS"),
        "execution_SHA256SUMS_sha256": sha256_file(execution / "SHA256SUMS"),
        "target_genome_sha256": sha256_file(
            normalized / "target_et39/primary_chromosomes.genome.fa"
        ),
        "perturbed_gff3_sha256": sha256_file(benchmark / "perturbed.gff3"),
        **candidate,
    }


def _build_formal_raw_manifest(
    *, project_root: Path, role_root: Path, benchmark: Path,
    protocol: Path, execution: Path, output: Path
) -> None:
    execution_manifest, _, _ = verify_execution(execution, protocol)
    if execution_manifest.get("freeze_stage") != (
        "post_evaluator_truth_two_complete_blind_runs_pre_label_"
        "reproducibility_patch_3"
    ):
        raise ValueError("Formal reconciliation requires the Coffea patch-3 freeze")
    raw_root = project_root / "results/baselines/coffea/v1.0"
    trees = {
        f"{method}__{bundle}": raw_root / method / bundle
        for method in METHODS
        for bundle in BUNDLES
    }
    payload = build_raw_predictions_manifest(
        holdout_id=HOLDOUT_ID,
        policy_id=POLICY_ID,
        project_root=project_root,
        raw_prediction_trees=trees,
        candidate_references=BUNDLES,
        input_hashes=_formal_input_hashes(
            role_root=role_root,
            benchmark=benchmark,
            protocol=protocol,
            execution=execution,
            project_root=project_root,
        ),
        output=output,
    )
    patch = execution_manifest.get("execution_patch", {})
    payload["reproducibility_reconciliation"] = {
        "execution_patch_sequence": patch.get("patch_sequence"),
        "selection_rule": (
            "same_model_id_and_exact_complete_decision_row_in_both_runs"
        ),
        "run_a_manifest_sha256": patch.get(
            "reproducibility_run_a_manifest_sha256"
        ),
        "run_b_manifest_sha256": patch.get(
            "reproducibility_run_b_manifest_sha256"
        ),
    }
    temporary = output.with_name(output.name + ".working")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    os.replace(temporary, output)


def _manifest_rows(root: Path, relatives: Iterable[Path]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for relative in sorted(relatives, key=lambda value: value.as_posix()):
        path = root / relative
        if not path.is_file() or path.is_symlink() or path.stat().st_size == 0:
            raise ValueError(f"Missing reconciliation input: {path}")
        rows.append(
            {
                "relative_path": relative.as_posix(),
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
        )
    return rows


def _write_manifest_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("relative_path", "bytes", "sha256"),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def _merge_stable_candidates(
    *, stable: dict[tuple[str, str], Path], base_gff: Path, working: Path
) -> None:
    for scope, bundles in SCOPES.items():
        method_candidates: dict[str, Path] = {}
        for method in METHODS:
            merged = working / "scopes" / scope / "methods" / method
            merged.mkdir(parents=True)
            merged_candidate_only = merged / "candidate_only.gff3"
            command = [
                sys.executable,
                "-m",
                "ploidypatch.cli",
                "baseline",
                "merge-candidate-gffs",
            ]
            for bundle in bundles:
                command += ["--candidate", f"{bundle}={stable[(method, bundle)]}"]
            command += [
                "--output-gff",
                str(merged_candidate_only),
                "--provenance-tsv",
                str(merged / "provenance.tsv"),
            ]
            run_command(
                command,
                stdout=merged / "stdout.json",
                stderr=merged / "stderr.log",
            )
            deduplicated = merged / "candidate_deduplicated.gff3"
            collapse_duplicate_exact_chains(
                candidate_gff=merged_candidate_only, output=deduplicated
            )
            assemble_method_candidate_gff(
                base_gff=base_gff,
                merged_candidate_gff=deduplicated,
                output=merged / "candidate.gff3",
            )
            method_candidates[method] = merged / "candidate.gff3"

        for arm, (redundancy_policy, policy_arm) in ARMS.items():
            destination = working / scope / arm / "blind"
            destination.mkdir(parents=True)
            command = [
                sys.executable,
                "-m",
                "ploidypatch.cli",
                "baseline",
                "select-method-consensus",
                "--base-gff",
                str(base_gff),
            ]
            for method in METHODS:
                command += ["--candidate", f"{method}={method_candidates[method]}"]
            command += [
                "--min-method-support",
                "1",
                "--max-redundancy-overlap",
                "0.5",
                "--redundancy-policy",
                redundancy_policy,
                "--output-gff",
                str(destination / "candidate.gff3"),
                "--decisions-tsv",
                str(destination / "decisions.tsv"),
            ]
            run_command(
                command,
                stdout=destination / "stdout.json",
                stderr=destination / "stderr.log",
            )
            seal_blind_pool_manifest(
                holdout_id=HOLDOUT_ID,
                policy_id=POLICY_ID,
                manifest_path=destination / "candidate.gff3.manifest.json",
                candidate_gff_path=destination / "candidate.gff3",
                decisions_path=destination / "decisions.tsv",
                raw_predictions_manifest_path=working / "raw_predictions.manifest.json",
                policy_arm=policy_arm,
                reference_scope=scope,
                candidate_reference_count=len(bundles),
            )


def reconcile(
    *, pools_a: Path, pools_b: Path, base_gff: Path,
    raw_predictions_manifest: Path | None, output: Path,
    formal_project_root: Path | None = None,
    role_root: Path | None = None,
    protocol: Path | None = None,
    execution: Path | None = None,
) -> Path:
    if output.exists() or output.is_symlink():
        raise FileExistsError(f"Refusing to overwrite reconciliation output: {output}")
    for source in (pools_a, pools_b):
        if not source.is_dir() or source.is_symlink():
            raise ValueError(f"Missing or symlinked blind pool tree: {source}")
    if not base_gff.is_file() or base_gff.is_symlink():
        raise ValueError("Reconciliation requires a truth-free perturbed GFF")
    formal_values = (formal_project_root, role_root, protocol, execution)
    formal_mode = all(value is not None for value in formal_values)
    if any(value is not None for value in formal_values) and not formal_mode:
        raise ValueError("Formal reconciliation context is incomplete")
    if formal_mode == (raw_predictions_manifest is not None):
        raise ValueError("Choose exactly one raw-manifest construction mode")

    output.parent.mkdir(parents=True, exist_ok=True)
    working = Path(tempfile.mkdtemp(prefix=f".{output.name}.working.", dir=output.parent))
    try:
        relatives = [
            Path("adapted") / method / bundle / name
            for method in METHODS
            for bundle in BUNDLES
            for name in ("decisions.tsv", "candidate_only.gff3")
        ]
        rows_a = _manifest_rows(pools_a, relatives)
        rows_b = _manifest_rows(pools_b, relatives)
        _write_manifest_tsv(working / "blind_run_a_inputs.tsv", rows_a)
        _write_manifest_tsv(working / "blind_run_b_inputs.tsv", rows_b)
        raw_output = working / "raw_predictions.manifest.json"
        if formal_mode:
            assert formal_project_root is not None and role_root is not None
            assert protocol is not None and execution is not None
            _build_formal_raw_manifest(
                project_root=formal_project_root,
                role_root=role_root,
                benchmark=base_gff.parent,
                protocol=protocol,
                execution=execution,
                output=raw_output,
            )
        else:
            assert raw_predictions_manifest is not None
            shutil.copyfile(raw_predictions_manifest, raw_output)
        raw_value = json.loads(raw_output.read_text(encoding="utf-8"))
        if (
            raw_value.get("schema_version")
            != "ploidypatch.core_h1_raw_predictions.v1"
            or raw_value.get("truth_access") is not False
            or raw_value.get("ranker_access") is not False
        ):
            raise ValueError("Reconciliation raw manifest is not truth blind")

        stable: dict[tuple[str, str], Path] = {}
        arm_summaries: dict[str, object] = {}
        for method in METHODS:
            for bundle in BUNDLES:
                key = (method, bundle)
                decisions_a = pools_a / "adapted" / method / bundle / "decisions.tsv"
                decisions_b = pools_b / "adapted" / method / bundle / "decisions.tsv"
                comparison = compare_decision_tables(decisions_a, decisions_b)
                audit_root = working / "reproducibility" / method / bundle
                audit_root.mkdir(parents=True)
                audit = write_comparison_audit(
                    comparison=comparison,
                    decisions_a=decisions_a,
                    decisions_b=decisions_b,
                    output_tsv=audit_root / "model_stability.tsv",
                    output_json=audit_root / "model_stability.json",
                    method=method,
                    reference=bundle,
                )
                filtered = working / "stable_adapted" / method / bundle / "candidate_only.gff3"
                counts = filter_candidate_gff_by_upstream_models(
                    source=pools_b / "adapted" / method / bundle / "candidate_only.gff3",
                    allowed_models=comparison.stable_accepted_models,
                    output=filtered,
                )
                stable[key] = filtered
                arm_summaries[f"{method}__{bundle}"] = {**audit, **counts}

        _merge_stable_candidates(stable=stable, base_gff=base_gff, working=working)
        manifest = {
            "schema_version": "ploidypatch.coffea_blind_reproducibility_reconciliation.v1.1",
            "holdout_id": HOLDOUT_ID,
            "policy_id": POLICY_ID,
            "truth_access": False,
            "label_access": False,
            "ranker_access": False,
            "automatic_approval": False,
            "selection_rule": "exact_complete_adapter_decision_row_intersection_across_two_independent_blind_runs",
            "selection_run": "run_b_after_exact_intersection",
            "biological_rules_or_thresholds_changed": False,
            "unstable_projection_policy": "abstain",
            "base_gff_sha256": sha256_file(base_gff),
            "raw_predictions_manifest_sha256": sha256_file(
                working / "raw_predictions.manifest.json"
            ),
            "blind_run_a_inputs_manifest_sha256": sha256_file(
                working / "blind_run_a_inputs.tsv"
            ),
            "blind_run_b_inputs_manifest_sha256": sha256_file(
                working / "blind_run_b_inputs.tsv"
            ),
            "arms": arm_summaries,
        }
        (working / "reconciliation_manifest.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        write_sha256sums(working)
        verify_sha256sums(working, ignore_checksum_file=True)
        os.replace(working, output)
    except BaseException:
        failed = output.with_name(output.name + ".invalid_run")
        if working.exists() and not working.is_symlink():
            if failed.exists() or failed.is_symlink():
                raise RuntimeError(f"Reconciliation failure tree exists: {failed}")
            os.replace(working, failed)
        raise
    return output


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pools-a", required=True)
    parser.add_argument("--pools-b", required=True)
    parser.add_argument("--base-gff", required=True)
    parser.add_argument("--raw-predictions-manifest")
    parser.add_argument("--formal-project-root")
    parser.add_argument("--role-root")
    parser.add_argument("--protocol-freeze")
    parser.add_argument("--execution-freeze")
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args(argv)
    reconcile(
        pools_a=Path(args.pools_a).resolve(),
        pools_b=Path(args.pools_b).resolve(),
        base_gff=Path(args.base_gff).resolve(),
        raw_predictions_manifest=(
            Path(args.raw_predictions_manifest).resolve()
            if args.raw_predictions_manifest
            else None
        ),
        output=Path(args.output_dir).resolve(),
        formal_project_root=(
            Path(args.formal_project_root).resolve()
            if args.formal_project_root
            else None
        ),
        role_root=Path(args.role_root).resolve() if args.role_root else None,
        protocol=(
            Path(args.protocol_freeze).resolve() if args.protocol_freeze else None
        ),
        execution=(
            Path(args.execution_freeze).resolve() if args.execution_freeze else None
        ),
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
