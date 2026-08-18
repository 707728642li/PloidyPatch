#!/usr/bin/env python3
"""Adapt six Coffea projections into six frozen core-H1 pool arms."""
from __future__ import annotations

import csv
import os
from pathlib import Path
import shutil
import sys
import tempfile

from ploidypatch.artifact_manifest import sha256_file, verify_sha256sums, write_sha256sums
from ploidypatch.core_h1_pools import (
    assemble_method_candidate_gff,
    build_raw_predictions_manifest,
    collapse_duplicate_exact_chains,
    extract_candidate_suffix,
    load_json_object,
    run_command,
    seal_blind_pool_manifest,
)
from ploidypatch.known_subgenome_h1 import HOLDOUT_ID, POLICY_ID
from ploidypatch.published_output import PUBLISHED_ROLES, verify_published_method_output


BUNDLES = ("candidate_bua", "candidate_mauritiana")
METHODS = ("miniprot", "gemoma", "lifton")
SCOPES = {
    "combined": BUNDLES,
    "bua_only": ("candidate_bua",),
    "mauritiana_only": ("candidate_mauritiana",),
}
RESULT_RELATIVE = Path("results/copy_collapse/external/coffea_v1.0_h1")


def role_hashes(staged: Path) -> dict[str, str]:
    with (staged / "role_manifest.tsv").open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    result: dict[str, str] = {}
    species = {
        "candidate_bua": "Coffea_eugenioides_BuA",
        "candidate_mauritiana": "Coffea_mauritiana",
    }
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
                raise ValueError(f"Non-unique Coffea staged role: {bundle}/{artifact}")
            result[f"{bundle}_{artifact}_sha256"] = matches[0]["sha256"]
    return result


def build_pools(
    *,
    project_root: Path,
    raw_baselines: Path,
    base_gff: Path,
    output: Path,
    include_raw_manifest: bool,
    raw_manifest_path: Path | None = None,
    seal_manifests: bool = True,
) -> None:
    if output.exists() or output.is_symlink():
        raise FileExistsError(f"Refusing to overwrite Coffea H1 pools: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    working = Path(tempfile.mkdtemp(prefix=f".{output.name}.working.", dir=output.parent))
    try:
        raw_trees = {
            f"{method}__{bundle}": raw_baselines / method / bundle
            for method in METHODS
            for bundle in BUNDLES
        }
        for label, tree in raw_trees.items():
            method = label.split("__", 1)[0]
            if method == "miniprot":
                verify_sha256sums(tree, ignore_checksum_file=True)
            else:
                verify_published_method_output(tree, PUBLISHED_ROLES[method])
        if include_raw_manifest:
            staged = Path(os.environ["PLOIDYPATCH_STAGED_INPUT_ROOT"]).resolve()
            benchmark = Path(os.environ["PLOIDYPATCH_BLIND_BENCHMARK_ROOT"]).resolve()
            protocol = Path(os.environ["PLOIDYPATCH_PROTOCOL_FREEZE"]).resolve()
            execution = Path(os.environ["PLOIDYPATCH_EXECUTION_FREEZE"]).resolve()
            blind_role = load_json_object(staged / "role_manifest.json")
            normalized = raw_baselines / "normalized"
            hashes = {
                "staged_input_SHA256SUMS_sha256": blind_role[
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
                "perturbed_gff3_sha256": sha256_file(base_gff),
                **role_hashes(staged),
            }
            raw_manifest_path = working / "raw_predictions.manifest.json"
            build_raw_predictions_manifest(
                holdout_id=HOLDOUT_ID,
                policy_id=POLICY_ID,
                project_root=project_root,
                raw_prediction_trees=raw_trees,
                candidate_references=BUNDLES,
                input_hashes=hashes,
                output=raw_manifest_path,
            )
        if raw_manifest_path is None:
            raise ValueError("Coffea pool construction requires frozen raw manifest")

        adapted: dict[tuple[str, str], Path] = {}
        for method in METHODS:
            for bundle in BUNDLES:
                destination = working / "adapted" / method / bundle
                destination.mkdir(parents=True)
                common = [sys.executable, "-m", "ploidypatch.cli", "baseline"]
                if method == "miniprot":
                    command = [
                        *common,
                        "adapt-miniprot",
                        "--perturbed-gff",
                        str(base_gff),
                        "--miniprot-gff",
                        str(raw_baselines / method / bundle / "raw/miniprot.gff3"),
                        "--protein-map",
                        str(raw_baselines / method / bundle / "reference/protein.map.tsv"),
                        "--min-identity",
                        "0.5",
                        "--min-query-coverage",
                        "0.5",
                        "--max-existing-cds-overlap",
                        "0.2",
                        "--max-redundancy-overlap",
                        "0.5",
                    ]
                else:
                    raw = raw_baselines / method / bundle / "upstream" / (
                        "final_annotation.gff" if method == "gemoma" else "lifton.gff3"
                    )
                    command = [
                        *common,
                        "adapt-gff",
                        "--perturbed-gff",
                        str(base_gff),
                        "--candidate-gff",
                        str(raw),
                        "--source",
                        method,
                        "--max-existing-cds-overlap",
                        "0.2",
                        "--max-redundancy-overlap",
                        "0.5",
                    ]
                command += [
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
                candidate_only = destination / "candidate_only.gff3"
                extract_candidate_suffix(
                    base_gff=base_gff,
                    adapted_gff=destination / "candidate.gff3",
                    output=candidate_only,
                )
                adapted[(method, bundle)] = candidate_only

        for scope, scope_bundles in SCOPES.items():
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
                for bundle in scope_bundles:
                    command += ["--candidate", f"{bundle}={adapted[(method, bundle)]}"]
                command += [
                    "--output-gff",
                    str(merged_candidate_only),
                    "--provenance-tsv",
                    str(merged / "provenance.tsv"),
                ]
                run_command(
                    command, stdout=merged / "stdout.json", stderr=merged / "stderr.log"
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

            arms = {
                "retain_distinct": (
                    "retain_distinct_chains",
                    "retain_distinct_phased_CDS_chains",
                ),
                "suppress_overlap": (
                    "suppress_overlapping",
                    "suppress_strongly_overlapping_alternative_chains",
                ),
            }
            for arm, (redundancy_policy, policy_arm) in arms.items():
                destination = working / scope / arm / (
                    "blind" if seal_manifests else "complete_control"
                )
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
                if seal_manifests:
                    seal_blind_pool_manifest(
                        holdout_id=HOLDOUT_ID,
                        policy_id=POLICY_ID,
                        manifest_path=destination / "candidate.gff3.manifest.json",
                        candidate_gff_path=destination / "candidate.gff3",
                        decisions_path=destination / "decisions.tsv",
                        raw_predictions_manifest_path=raw_manifest_path,
                        policy_arm=policy_arm,
                        reference_scope=scope,
                        candidate_reference_count=len(scope_bundles),
                    )
        write_sha256sums(working)
        verify_sha256sums(working, ignore_checksum_file=True)
        os.replace(working, output)
    except BaseException:
        failed = output.with_name(output.name + ".invalid_run")
        if working.exists() and not working.is_symlink():
            if failed.exists() or failed.is_symlink():
                raise RuntimeError(f"Coffea candidate-pool failure tree exists: {failed}")
            os.replace(working, failed)
        raise


def main(argv: list[str] | None = None) -> int:
    values = sys.argv[1:] if argv is None else argv
    if len(values) != 1:
        raise SystemExit("usage: build_coffea_h1_candidate_pools_v1.0.py PROJECT_ROOT")
    project_root = Path(values[0]).resolve()
    if os.environ.get("PLOIDYPATCH_BLIND_RUNNER") != "1":
        raise ValueError("Coffea H1 pools require the blind runner")
    if os.environ.get("PLOIDYPATCH_NETWORK_ACCESS") != "none":
        raise ValueError("Coffea H1 pools require network=none")
    if any(
        (Path("/holdout") / name).exists()
        for name in ("evaluator_only", "truth", "labels", "target_complete")
    ):
        raise ValueError("Forbidden evaluator role is visible to Coffea pool builder")
    build_pools(
        project_root=project_root,
        raw_baselines=project_root / "results/baselines/coffea/v1.0",
        base_gff=Path(os.environ["PLOIDYPATCH_BLIND_BENCHMARK_ROOT"]).resolve()
        / "perturbed.gff3",
        output=project_root / RESULT_RELATIVE,
        include_raw_manifest=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

