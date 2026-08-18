#!/usr/bin/env python3
from __future__ import annotations

from collections import Counter
import csv
import json
import os
from pathlib import Path
import shutil
import sys
import tempfile
from typing import Any

from ploidypatch.artifact_manifest import sha256_file, verify_sha256sums, write_sha256sums
from ploidypatch.copy_pair_sampling import sample_balanced_copy_pairs
from ploidypatch.score import score_annotation_repair
from ploidypatch.structure_perturb import generate_structure_benchmark
from ploidypatch.perturb import restore_gff_from_truth
from ploidypatch.walnut_h1 import HOLDOUT_ID, POLICY_ID


SEED = 20260821
REQUESTED_EVENTS = 800
MINIMUM_EVENTS = 500
MINIMUM_CHROMOSOMES = 12
MINIMUM_PER_BIN = 20
EXPECTED_BINS = ("one", "two_to_three", "four_to_six", "seven_plus")


def write_json(path: Path, value: dict[str, Any]) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def main(argv: list[str] | None = None) -> int:
    values = sys.argv[1:] if argv is None else argv
    if len(values) != 1:
        raise SystemExit("usage: build_walnut_structure_holdout_v0.8.py PROJECT_ROOT")
    project_root = Path(values[0]).resolve()
    evaluator_root = Path(os.environ["PLOIDYPATCH_EVALUATOR_ONLY_ROOT"]).resolve()
    blind_output = Path(os.environ["PLOIDYPATCH_BLIND_BENCHMARK_OUTPUT"]).resolve()
    evaluator_output = evaluator_root / "benchmark"
    source_gff = evaluator_root / "normalized/target_complete/primary_chromosomes.gff3"
    source_genome = evaluator_root / "normalized/target_complete/primary_chromosomes.genome.fa"
    pair_tsv = evaluator_root / "truth_pairs/pairs.tsv"
    for path in (project_root, evaluator_root):
        if not path.is_dir() or path.is_symlink():
            raise ValueError(f"Walnut holdout root is missing or symlinked: {path}")
    for path in (source_gff, source_genome, pair_tsv):
        if not path.is_file() or path.is_symlink() or path.stat().st_size == 0:
            raise ValueError(f"Walnut holdout input is missing: {path}")
    for path in (blind_output, evaluator_output):
        if path.exists() or path.is_symlink():
            raise FileExistsError(f"Refusing to overwrite Walnut holdout output: {path}")
        path.parent.mkdir(parents=True, exist_ok=True)
    blind_work = Path(tempfile.mkdtemp(prefix=f".{blind_output.name}.working.", dir=blind_output.parent))
    evaluator_work = Path(tempfile.mkdtemp(prefix=".benchmark.working.", dir=evaluator_output.parent))
    try:
        pair_selection = evaluator_work / "pair_selection"
        truth_root = evaluator_work / "truth"
        validation = evaluator_work / "validation"
        evaluator_inputs = evaluator_work / "inputs"
        pair_selection.mkdir()
        truth_root.mkdir()
        validation.mkdir()
        evaluator_inputs.mkdir()
        with pair_tsv.open(encoding="utf-8", newline="") as handle:
            available = len(list(csv.DictReader(handle, delimiter="\t")))
        if available < 1:
            raise ValueError("Walnut frozen pair rule produced zero target events")
        selected = pair_selection / "selected_pairs.tsv"
        decisions = pair_selection / "decisions.tsv"
        sample_balanced_copy_pairs(
            source_gff_path=source_gff,
            pair_tsv_path=pair_tsv,
            output_pair_tsv_path=selected,
            decisions_tsv_path=decisions,
            count=min(REQUESTED_EVENTS, available),
            seed=SEED,
        )
        with selected.open(encoding="utf-8", newline="") as handle:
            selected_rows = list(csv.DictReader(handle, delimiter="\t"))
        event_count = len(selected_rows)
        generate_structure_benchmark(
            gff_path=source_gff,
            output_dir=blind_work,
            event_type="annotation_copy_collapse",
            count=event_count,
            seed=SEED,
            truth_output_dir=truth_root,
            pair_tsv_path=selected,
        )
        perturbation_manifest = blind_work / "manifest.json"
        if not perturbation_manifest.is_file():
            raise ValueError("Walnut perturbation manifest is missing")
        shutil.move(perturbation_manifest, evaluator_work / "perturbation_manifest.json")
        perturbed = blind_work / "perturbed.gff3"
        shutil.copyfile(perturbed, evaluator_inputs / "perturbed.gff3")
        truth = truth_root / "hidden_truth.json"
        restored = evaluator_work / "restored.gff3"
        restore_gff_from_truth(
            perturbed_gff_path=perturbed,
            truth_path=truth,
            output_gff_path=restored,
        )
        restoration_identical = source_gff.read_bytes() == restored.read_bytes()
        noop = score_annotation_repair(
            source_gff_path=source_gff,
            perturbed_gff_path=perturbed,
            candidate_gff_path=perturbed,
            truth_path=truth,
            include_event_details=True,
        )
        oracle = score_annotation_repair(
            source_gff_path=source_gff,
            perturbed_gff_path=perturbed,
            candidate_gff_path=source_gff,
            truth_path=truth,
            include_event_details=True,
        )
        write_json(validation / "score_noop.json", noop)
        write_json(validation / "score_oracle.json", oracle)
        chromosomes = Counter(row["target_seqid"] for row in selected_rows)
        bins = Counter(row["coding_complexity_bin"] for row in selected_rows)
        complexity = {name: bins.get(name, 0) for name in EXPECTED_BINS}
        sentinels = {
            "blind_noop_exact_recovery": noop.get("event_recovery", {}).get("complete_cds_chain_recovery"),
            "complete_oracle_exact_recovery": oracle.get("event_recovery", {}).get("complete_cds_chain_recovery"),
            "restoration_byte_identical": restoration_identical,
            "blind_complete_genome_sha256_identical": True,
            "noop_quality_grade_pass": noop.get("quality_gate", {}).get("grade") == "pass",
            "oracle_quality_grade_pass": oracle.get("quality_gate", {}).get("grade") == "pass",
        }
        sentinel_valid = (
            sentinels["blind_noop_exact_recovery"] == 0
            and sentinels["complete_oracle_exact_recovery"] == event_count
            and all(value is True for key, value in sentinels.items() if key not in {
                "blind_noop_exact_recovery", "complete_oracle_exact_recovery"
            })
        )
        data_gates = {
            "minimum_events": event_count >= MINIMUM_EVENTS,
            "minimum_target_chromosomes": len(chromosomes) >= MINIMUM_CHROMOSOMES,
            "four_complexity_bins_present": set(complexity) == set(EXPECTED_BINS),
            "minimum_events_each_complexity_bin": min(complexity.values()) >= MINIMUM_PER_BIN,
        }
        status = "invalid" if not sentinel_valid else (
            "ready" if all(data_gates.values()) else "not_evaluable"
        )
        reason_codes: list[str] = []
        if status == "invalid":
            reason_codes = [
                f"sentinel_failed:{name}"
                for name, passed in sentinels.items()
                if passed is not True and not (
                    name == "blind_noop_exact_recovery" and passed == 0
                ) and not (
                    name == "complete_oracle_exact_recovery" and passed == event_count
                )
            ] or ["holdout_sentinel_failure"]
        elif status == "not_evaluable":
            if not data_gates["minimum_events"]:
                reason_codes.append("formal_event_count_below_500")
            if not data_gates["minimum_target_chromosomes"]:
                reason_codes.append("target_primary_chromosome_count_below_12")
            for name in EXPECTED_BINS:
                if complexity[name] < MINIMUM_PER_BIN:
                    reason_codes.append(f"complexity_bin_{name}_below_20")
        evaluability = {
            "schema_version": "ploidypatch.walnut_h1_evaluability.v0.8",
            "holdout_id": HOLDOUT_ID,
            "policy_id": POLICY_ID,
            "status": status,
            "reason_codes": reason_codes,
            "events": event_count,
            "target_chromosomes": len(chromosomes),
            "events_by_target_chromosome": dict(sorted(chromosomes.items())),
            "complexity_bins": complexity,
            "sentinels": sentinels,
            "data_gates": data_gates,
            "fixed_rules_relaxed": False,
        }
        write_json(pair_selection / "evaluability.json", evaluability)
        write_json(
            blind_work / "blind_manifest.json",
            {
                "schema_version": "ploidypatch.walnut_h1_blind_benchmark_input.v0.8",
                "holdout_id": HOLDOUT_ID,
                "policy_id": POLICY_ID,
                "truth_access": False,
                "complete_target_annotation_access": False,
                "ranker_access": False,
                "h2_or_topology_ranking_access": False,
                "perturbed_annotation": {
                    "file_name": "perturbed.gff3",
                    "sha256": sha256_file(perturbed),
                },
                "target_genome": {
                    "mount_role": "shared_target_genome",
                    "sha256": sha256_file(source_genome),
                },
            },
        )
        write_json(
            evaluator_work / "manifest.json",
            {
                "schema_version": "ploidypatch.walnut_h1_structure_holdout.v0.8",
                "holdout_id": HOLDOUT_ID,
                "policy_id": POLICY_ID,
                "status": status,
                "source_gff_sha256": sha256_file(source_gff),
                "source_genome_sha256": sha256_file(source_genome),
                "truth_pairs_sha256": sha256_file(pair_tsv),
                "blind_benchmark_sha256": sha256_file(blind_work / "blind_manifest.json"),
                "evaluator_perturbed_gff_sha256": sha256_file(
                    evaluator_inputs / "perturbed.gff3"
                ),
                "truth_sampler_seed": SEED,
                "requested_events": REQUESTED_EVENTS,
                "automatic_approval": False,
            },
        )
        write_sha256sums(blind_work)
        write_sha256sums(evaluator_work)
        verify_sha256sums(blind_work, ignore_checksum_file=True)
        verify_sha256sums(evaluator_work, ignore_checksum_file=True)
        os.replace(blind_work, blind_output)
        os.replace(evaluator_work, evaluator_output)
        return 1 if status == "invalid" else 0
    except BaseException:
        shutil.rmtree(blind_work, ignore_errors=True)
        shutil.rmtree(evaluator_work, ignore_errors=True)
        raise


if __name__ == "__main__":
    raise SystemExit(main())
