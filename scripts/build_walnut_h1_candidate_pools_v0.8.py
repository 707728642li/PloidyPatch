#!/usr/bin/env python3
"""Adapt six frozen raw projections into the two Walnut core-H1 pools."""
from __future__ import annotations

import csv
from collections import defaultdict
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
from typing import Iterable

from ploidypatch.artifact_manifest import sha256_file, verify_sha256sums, write_sha256sums
from ploidypatch.gff import parse_attributes
from ploidypatch.published_output import PUBLISHED_ROLES, verify_published_method_output
from ploidypatch.walnut_h1 import (
    build_raw_predictions_manifest,
    load_json_object,
    seal_blind_pool_manifest,
)


BUNDLES = ("candidate_mandshurica", "candidate_carya")
METHODS = ("miniprot", "gemoma", "lifton")
RESULT_RELATIVE = Path("results/copy_collapse/external/walnut_v0.8_h1")


def run(command: Iterable[str], *, stdout: Path, stderr: Path) -> None:
    stdout.parent.mkdir(parents=True, exist_ok=True)
    with stdout.open("x", encoding="utf-8") as out, stderr.open(
        "x", encoding="utf-8"
    ) as err:
        completed = subprocess.run(list(command), stdout=out, stderr=err, check=False)
    if completed.returncode != 0:
        raise RuntimeError(f"Walnut candidate-pool command exited {completed.returncode}")


def role_hashes(staged: Path) -> dict[str, str]:
    with (staged / "role_manifest.tsv").open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    result: dict[str, str] = {}
    species = {
        "candidate_mandshurica": "Juglans_mandshurica",
        "candidate_carya": "Carya_illinoinensis",
    }
    for bundle, species_id in species.items():
        for artifact in ("genome", "gff3", "protein"):
            matches = [
                row for row in rows
                if row.get("role") == "candidate_reference"
                and row.get("species_id") == species_id
                and row.get("artifact") == artifact
            ]
            if len(matches) != 1:
                raise ValueError(f"Non-unique Walnut staged role: {bundle}/{artifact}")
            result[f"{bundle}_{artifact}_sha256"] = matches[0]["sha256"]
    return result


def extract_candidate_suffix(*, base_gff: Path, adapted_gff: Path, output: Path) -> None:
    """Remove one byte-identical base prefix before cross-reference merging."""
    if output.exists() or output.is_symlink():
        raise FileExistsError(f"Refusing to overwrite Walnut candidate suffix: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    with base_gff.open("rb") as base, adapted_gff.open("rb") as adapted, output.open(
        "xb"
    ) as destination:
        last_byte = b""
        for chunk in iter(lambda: base.read(1024 * 1024), b""):
            if adapted.read(len(chunk)) != chunk:
                raise ValueError(
                    f"Walnut adapted candidate does not preserve base prefix: {adapted_gff}"
                )
            last_byte = chunk[-1:]
        if last_byte not in {b"\n", b"\r"} and adapted.read(1) != b"\n":
            raise ValueError(f"Walnut adapted candidate lacks base separator: {adapted_gff}")
        if adapted.readline() not in {b"###\n", b"###\r\n"}:
            raise ValueError(f"Walnut adapted candidate lacks boundary marker: {adapted_gff}")
        destination.write(b"##gff-version 3\n")
        shutil.copyfileobj(adapted, destination, length=1024 * 1024)


def assemble_method_candidate_gff(
    *, base_gff: Path, merged_candidate_gff: Path, output: Path
) -> None:
    """Restore the exact base once, followed by merged candidate-only features."""
    if output.exists() or output.is_symlink():
        raise FileExistsError(f"Refusing to overwrite Walnut method candidate GFF: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("xb") as destination, base_gff.open("rb") as base:
        last_byte = b""
        for chunk in iter(lambda: base.read(1024 * 1024), b""):
            destination.write(chunk)
            last_byte = chunk[-1:]
        if last_byte not in {b"\n", b"\r"}:
            destination.write(b"\n")
        destination.write(b"###\n")
        with merged_candidate_gff.open("rb") as merged:
            for raw in merged:
                if raw.strip() and not raw.startswith(b"#"):
                    if len(raw.rstrip(b"\r\n").split(b"\t")) != 9:
                        raise ValueError("Malformed merged Walnut candidate feature")
                    destination.write(raw)


def collapse_duplicate_exact_chains(*, candidate_gff: Path, output: Path) -> dict[str, int]:
    """Give one within-method vote to an exact phased-CDS chain."""
    manifest_path = Path(str(output) + ".manifest.json")
    if any(path.exists() or path.is_symlink() for path in (output, manifest_path)):
        raise FileExistsError(f"Refusing to overwrite Walnut exact-chain collapse: {output}")
    transcripts: dict[str, tuple[str, str, str]] = {}
    cds: dict[str, list[tuple[int, int, str]]] = defaultdict(list)
    gene_transcripts: dict[str, set[str]] = defaultdict(set)
    with candidate_gff.open(encoding="utf-8", newline="") as handle:
        for line_number, raw in enumerate(handle, start=1):
            stripped = raw.rstrip("\r\n")
            if not stripped or stripped.startswith("#"):
                continue
            fields = stripped.split("\t")
            if len(fields) != 9 or fields[2] not in {"gene", "mRNA", "transcript", "exon", "CDS"}:
                raise ValueError(f"Unsupported merged Walnut feature at line {line_number}")
            attributes, malformed = parse_attributes(fields[8])
            if malformed:
                raise ValueError(f"Malformed merged Walnut attributes at line {line_number}")
            if fields[2] in {"mRNA", "transcript"}:
                transcript_id = attributes.get("ID", "")
                parents = tuple(value for value in attributes.get("Parent", "").split(",") if value)
                if not transcript_id or transcript_id in transcripts or len(parents) != 1:
                    raise ValueError(f"Invalid merged Walnut transcript at line {line_number}")
                transcripts[transcript_id] = (fields[0], fields[6], parents[0])
                gene_transcripts[parents[0]].add(transcript_id)
            elif fields[2] == "CDS":
                parents = tuple(value for value in attributes.get("Parent", "").split(",") if value)
                if len(parents) != 1:
                    raise ValueError(f"Invalid merged Walnut CDS parent at line {line_number}")
                try:
                    start, end = int(fields[3]), int(fields[4])
                except ValueError as error:
                    raise ValueError(f"Invalid merged Walnut CDS coordinate at line {line_number}") from error
                cds[parents[0]].append((start, end, fields[7]))
    if set(cds) - set(transcripts):
        raise ValueError("Merged Walnut CDS references an unknown transcript")
    grouped: dict[tuple[str, str, tuple[tuple[int, int, str], ...]], list[str]] = defaultdict(list)
    for transcript_id, (seqid, strand, _gene_id) in transcripts.items():
        chain = tuple(sorted(cds.get(transcript_id, [])))
        if not chain or any(phase not in {"0", "1", "2"} for _start, _end, phase in chain):
            raise ValueError(f"Merged Walnut transcript lacks a phased CDS: {transcript_id}")
        grouped[(seqid, strand, chain)].append(transcript_id)
    dropped_transcripts: set[str] = set()
    duplicate_groups = 0
    for values in grouped.values():
        if len(values) > 1:
            duplicate_groups += 1
            dropped_transcripts.update(sorted(values)[1:])
    dropped_genes = {
        gene_id
        for gene_id, values in gene_transcripts.items()
        if values and values <= dropped_transcripts
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    with candidate_gff.open(encoding="utf-8", newline="") as source, output.open(
        "x", encoding="utf-8", newline=""
    ) as destination:
        for raw in source:
            stripped = raw.rstrip("\r\n")
            if not stripped or stripped.startswith("#"):
                destination.write(raw)
                continue
            fields = stripped.split("\t")
            attributes, _malformed = parse_attributes(fields[8])
            identifier = attributes.get("ID", "")
            parents = {value for value in attributes.get("Parent", "").split(",") if value}
            if identifier in dropped_genes or identifier in dropped_transcripts:
                continue
            if parents & dropped_transcripts:
                continue
            destination.write(raw)
    audit = {
        "schema_version": "ploidypatch.walnut_within_method_exact_chain_collapse.v0.8",
        "input_sha256": sha256_file(candidate_gff),
        "output_sha256": sha256_file(output),
        "input_transcripts": len(transcripts),
        "retained_transcripts": len(transcripts) - len(dropped_transcripts),
        "collapsed_duplicate_transcripts": len(dropped_transcripts),
        "duplicate_chain_groups": duplicate_groups,
        "selection_rule": "lexicographically_smallest_namespaced_transcript_id",
        "truth_access": False,
    }
    manifest_path.write_text(
        json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return {
        key: int(value)
        for key, value in audit.items()
        if key in {
            "input_transcripts",
            "retained_transcripts",
            "collapsed_duplicate_transcripts",
            "duplicate_chain_groups",
        }
    }


def build_pools(
    *,
    project_root: Path,
    base_gff: Path,
    output: Path,
    include_raw_manifest: bool,
    raw_manifest_path: Path | None = None,
    seal_manifests: bool = True,
) -> None:
    code_root = project_root / "code"
    normalized = project_root / "results/baselines/walnut/v0.8/normalized"
    baselines = project_root / "results/baselines/walnut/v0.8"
    if output.exists() or output.is_symlink():
        raise FileExistsError(f"Refusing to overwrite Walnut H1 pools: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    working = Path(tempfile.mkdtemp(prefix=f".{output.name}.working.", dir=output.parent))
    try:
        raw_trees = {
            f"{method}__{bundle}": baselines / method / bundle
            for method in METHODS for bundle in BUNDLES
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
            blind_role = load_json_object(staged / "blind_role_manifest.json")
            hashes = {
                "staged_input_SHA256SUMS_sha256": blind_role[
                    "staged_input_SHA256SUMS_sha256"
                ],
                "blind_benchmark_SHA256SUMS_sha256": sha256_file(benchmark / "SHA256SUMS"),
                "protocol_SHA256SUMS_sha256": sha256_file(protocol / "SHA256SUMS"),
                "execution_SHA256SUMS_sha256": sha256_file(execution / "SHA256SUMS"),
                "target_genome_sha256": sha256_file(
                    normalized / "target_walnut2/primary_chromosomes.genome.fa"
                ),
                "perturbed_gff3_sha256": sha256_file(base_gff),
                **role_hashes(staged),
            }
            raw_manifest_path = working / "raw_predictions.manifest.json"
            build_raw_predictions_manifest(
                project_root=project_root,
                raw_prediction_trees=raw_trees,
                input_hashes=hashes,
                output=raw_manifest_path,
            )
        if raw_manifest_path is None:
            raise ValueError("Walnut pool construction requires frozen raw manifest")

        method_candidates: dict[str, Path] = {}
        for method in METHODS:
            adapted: list[tuple[str, Path]] = []
            for bundle in BUNDLES:
                destination = working / "adapted" / method / bundle
                destination.mkdir(parents=True)
                common = [
                    sys.executable, "-m", "ploidypatch.cli", "baseline",
                ]
                if method == "miniprot":
                    command = [
                        *common, "adapt-miniprot", "--perturbed-gff", str(base_gff),
                        "--miniprot-gff", str(baselines / method / bundle / "raw/miniprot.gff3"),
                        "--protein-map", str(baselines / method / bundle / "reference/protein.map.tsv"),
                        "--min-identity", "0.5", "--min-query-coverage", "0.5",
                        "--max-existing-cds-overlap", "0.2",
                        "--max-redundancy-overlap", "0.5",
                    ]
                else:
                    raw = (
                        baselines / method / bundle / "upstream"
                        / ("final_annotation.gff" if method == "gemoma" else "lifton.gff3")
                    )
                    command = [
                        *common, "adapt-gff", "--perturbed-gff", str(base_gff),
                        "--candidate-gff", str(raw), "--source", method,
                        "--max-existing-cds-overlap", "0.2",
                        "--max-redundancy-overlap", "0.5",
                    ]
                command += [
                    "--output-gff", str(destination / "candidate.gff3"),
                    "--decisions-tsv", str(destination / "decisions.tsv"),
                ]
                run(
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
                adapted.append((bundle, candidate_only))
            merged = working / "methods" / method
            merged.mkdir(parents=True)
            merged_candidate_only = merged / "candidate_only.gff3"
            command = [
                sys.executable, "-m", "ploidypatch.cli", "baseline",
                "merge-candidate-gffs",
            ]
            for bundle, path in adapted:
                command += ["--candidate", f"{bundle}={path}"]
            command += [
                "--output-gff", str(merged_candidate_only),
                "--provenance-tsv", str(merged / "provenance.tsv"),
            ]
            run(command, stdout=merged / "stdout.json", stderr=merged / "stderr.log")
            deduplicated = merged / "candidate_deduplicated.gff3"
            collapse_duplicate_exact_chains(
                candidate_gff=merged_candidate_only,
                output=deduplicated,
            )
            assemble_method_candidate_gff(
                base_gff=base_gff,
                merged_candidate_gff=deduplicated,
                output=merged / "candidate.gff3",
            )
            method_candidates[method] = merged / "candidate.gff3"

        arms = {
            "retain_distinct": (
                "retain_distinct_chains", "retain_distinct_phased_CDS_chains"
            ),
            "suppress_overlap": (
                "suppress_overlapping",
                "suppress_strongly_overlapping_alternative_chains",
            ),
        }
        for arm, (redundancy_policy, policy_arm) in arms.items():
            destination = working / arm / ("blind" if seal_manifests else "complete_control")
            destination.mkdir(parents=True)
            command = [
                sys.executable, "-m", "ploidypatch.cli", "baseline",
                "select-method-consensus", "--base-gff", str(base_gff),
            ]
            for method in METHODS:
                command += ["--candidate", f"{method}={method_candidates[method]}"]
            command += [
                "--min-method-support", "1", "--max-redundancy-overlap", "0.5",
                "--redundancy-policy", redundancy_policy,
                "--output-gff", str(destination / "candidate.gff3"),
                "--decisions-tsv", str(destination / "decisions.tsv"),
            ]
            run(
                command,
                stdout=destination / "stdout.json",
                stderr=destination / "stderr.log",
            )
            if seal_manifests:
                seal_blind_pool_manifest(
                    manifest_path=destination / "candidate.gff3.manifest.json",
                    candidate_gff_path=destination / "candidate.gff3",
                    decisions_path=destination / "decisions.tsv",
                    raw_predictions_manifest_path=raw_manifest_path,
                    policy_arm=policy_arm,
                )
        write_sha256sums(working)
        verify_sha256sums(working, ignore_checksum_file=True)
        os.replace(working, output)
    except BaseException:
        failed = output.with_name(output.name + ".invalid_run")
        if working.exists() and not working.is_symlink():
            if failed.exists() or failed.is_symlink():
                raise RuntimeError(
                    f"Walnut candidate-pool failure tree already exists: {failed}"
                )
            os.replace(working, failed)
        raise


def main(argv: list[str] | None = None) -> int:
    values = sys.argv[1:] if argv is None else argv
    if len(values) != 1:
        raise SystemExit("usage: build_walnut_h1_candidate_pools_v0.8.py PROJECT_ROOT")
    project_root = Path(values[0]).resolve()
    if os.environ.get("PLOIDYPATCH_BLIND_RUNNER") != "1":
        raise ValueError("Walnut H1 pools require the blind runner")
    if os.environ.get("PLOIDYPATCH_NETWORK_ACCESS") != "none":
        raise ValueError("Walnut H1 pools require network=none")
    if any((Path("/holdout") / name).exists() for name in ("evaluator_only", "truth", "labels")):
        raise ValueError("Forbidden evaluator role is visible to Walnut H1 pool builder")
    build_pools(
        project_root=project_root,
        base_gff=Path(os.environ["PLOIDYPATCH_BLIND_BENCHMARK_ROOT"]).resolve()
        / "perturbed.gff3",
        output=project_root / RESULT_RELATIVE,
        include_raw_manifest=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
