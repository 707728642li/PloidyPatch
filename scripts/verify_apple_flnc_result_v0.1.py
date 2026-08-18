#!/usr/bin/env python3
"""Independently verify the frozen apple Pychopper FLNC result."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import re
import stat
import sys
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path


EXPECTED_TISSUES = {
    "CRR1429911": "root",
    "CRR1429912": "stem",
    "CRR1429913": "leaf",
    "CRR1429914": "flower",
    "CRR1429915": "young_fruit_45DAB",
    "CRR1429916": "expanding_fruit_95DAB",
    "CRR1429917": "mature_fruit_145DAB",
}
EXPECTED_BUDGETS = [100, 146, 250, 292, 500, 583]
EXPECTED_ESTIMATORS = ["baseline", "v03_primary", "v04_guard"]
EVIDENCE_STATES = {
    "not_assessable",
    "no_qualifying_observation",
    "single_exon_span_support",
    "partial_junction_support",
    "all_junctions_supported",
    "full_chain_supported",
}
SHA_LINE = re.compile(r"^([0-9a-f]{64}) [ *](.+)$")
TSV_FIELD_SIZE_LIMIT = 1 << 30


class VerificationError(RuntimeError):
    """Raised when a frozen-result invariant fails."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise VerificationError(message)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path) -> dict:
    with path.open(encoding="utf-8") as handle:
        value = json.load(handle)
    require(isinstance(value, dict), f"{path}: JSON root is not an object")
    return value


def read_unique_key_value_tsv(path: Path) -> dict[str, str]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require(reader.fieldnames == ["field", "value"], f"{path}: invalid header")
        output: dict[str, str] = {}
        for row in reader:
            key = row["field"]
            require(key not in output, f"{path}: duplicate field {key}")
            output[key] = row["value"]
    return output


def read_tsv(path: Path):
    # Full read-support ID lists can legitimately exceed the stdlib default
    # of 131,072 characters.  Keep the limit explicit and finite so that the
    # independent verifier accepts the production schema without silently
    # truncating or special-casing any evidence row.
    csv.field_size_limit(TSV_FIELD_SIZE_LIMIT)
    with path.open(encoding="utf-8", newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def count_fastq_records(path: Path) -> int:
    count = 0
    with gzip.open(path, "rt", encoding="ascii", newline="") as handle:
        while True:
            name = handle.readline()
            if not name:
                break
            sequence = handle.readline()
            plus = handle.readline()
            quality = handle.readline()
            require(bool(sequence and plus and quality), f"{path}: truncated record")
            sequence = sequence.rstrip("\r\n")
            quality = quality.rstrip("\r\n")
            require(name.startswith("@") and plus.startswith("+"), f"{path}: malformed FASTQ")
            require(len(sequence) == len(quality), f"{path}: sequence/quality mismatch")
            count += 1
    return count


def verify_root_manifest(root: Path, *, verify_hashes: bool) -> tuple[int, str]:
    manifest = root / "SHA256SUMS"
    require(manifest.is_file(), f"missing {manifest}")
    expected: dict[str, str] = {}
    with manifest.open(encoding="utf-8", newline="") as handle:
        for line_number, raw in enumerate(handle, start=1):
            match = SHA_LINE.fullmatch(raw.rstrip("\r\n"))
            require(match is not None, f"SHA256SUMS:{line_number}: malformed line")
            digest, relative = match.groups()
            relative = relative.removeprefix("./")
            candidate = Path(relative)
            require(not candidate.is_absolute(), f"absolute checksum path: {relative}")
            require(".." not in candidate.parts, f"traversal checksum path: {relative}")
            normalized = candidate.as_posix()
            require(normalized not in expected, f"duplicate checksum path: {normalized}")
            expected[normalized] = digest
    actual = {
        path.relative_to(root).as_posix()
        for path in root.rglob("*")
        if path.is_file() and path.name != "SHA256SUMS"
    }
    missing = sorted(actual - set(expected))
    extra = sorted(set(expected) - actual)
    require(
        not missing and not extra,
        f"SHA256SUMS file universe is incomplete or has extras; unlisted={missing}, absent={extra}",
    )
    if verify_hashes:
        for relative, digest in expected.items():
            require(sha256_file(root / relative) == digest, f"checksum mismatch: {relative}")
    return len(expected), sha256_file(manifest)


def verify_read_only_tree(root: Path) -> None:
    for path in [root, *root.rglob("*")]:
        require(not path.is_symlink(), f"frozen result contains symlink: {path}")
        mode = stat.S_IMODE(path.stat().st_mode)
        require(mode & 0o222 == 0, f"frozen path retains write bits: {path}")


def verify_referenced_manifests(root: Path, project_root: Path, *, verify_hashes: bool) -> None:
    input_manifest = root / "input_manifest.tsv"
    with input_manifest.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require(reader.fieldnames == ["role", "bytes", "sha256", "path"], "invalid input manifest")
        roles: set[str] = set()
        for row in reader:
            require(row["role"] not in roles, f"duplicate input role {row['role']}")
            roles.add(row["role"])
            path = Path(row["path"])
            require(path.is_file(), f"missing referenced input: {path}")
            require(path.stat().st_size == int(row["bytes"]), f"input size mismatch: {path}")
            if verify_hashes:
                require(sha256_file(path) == row["sha256"], f"input checksum mismatch: {path}")
    required_roles = {
        "genome",
        "genome_fai",
        "base_gff",
        "candidate_gff",
        "review_rankings",
        "repeat_gff",
        "minimap2_index",
        "raw_alignment_freeze",
        "raw_ts_patch2_freeze",
        "reused_self_map",
        "rank_freeze",
        "ont_freeze",
        "te_freeze",
    }
    require(roles == required_roles, "input manifest role set changed")

    with (root / "freeze/code_manifest.tsv").open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require(reader.fieldnames == ["path", "sha256"], "invalid code manifest")
        seen: set[str] = set()
        for row in reader:
            require(row["path"] not in seen, f"duplicate code path {row['path']}")
            seen.add(row["path"])
            path = project_root / "code" / row["path"]
            if row["path"] == "envs/ploidypatch-pychopper/bin/pychopper":
                path = project_root / row["path"]
            require(path.is_file(), f"missing referenced code: {path}")
            require(sha256_file(path) == row["sha256"], f"code checksum mismatch: {path}")


def verify_contract(root: Path, expected_commit: str) -> dict[str, str]:
    contract = read_unique_key_value_tsv(root / "run_contract.tsv")
    exact = {
        "code_commit": expected_commit,
        "target": "Malus_domestica_GDDH13_v1.1",
        "validation_source": "Golden_Delicious_GSA_CRA021523_raw_ONT",
        "analysis_role": "descriptive_posthoc_natural_validation",
        "labels_seen": "true",
        "candidate_coordinates_modified": "false",
        "candidate_ranks_modified": "false",
        "evidence_thresholds_modified": "false",
        "pychopper_version": "2.7.10",
        "pychopper_kit": "PCS109",
        "pychopper_method": "phmm",
        "pychopper_autotune_numpy_seed": "20261005",
        "source_date_epoch": "0",
        "pychopper_rescued_segments_used": "false",
        "alignment_strand_source": "query_orientation",
        "single_exon_is_not_full_chain_positive": "true",
        "review_budgets": "100,146,250,292,500,583",
        "bootstrap_replicates": "20000",
        "bootstrap_seed": "20261004",
        "automatic_annotation_patch": "false",
    }
    for key, value in exact.items():
        require(contract.get(key) == value, f"run contract mismatch: {key}")
    require(contract.get("pandas_version", "").split(".")[0] == "2", "pandas is not 2.x")
    require(contract.get("tqdm_runtime_policy") == "Bioconda_curated_modern_runtime_metadata_exception", "tqdm policy changed")
    return contract


def verify_classification(root: Path, *, count_fastq: bool) -> dict:
    report = load_json(root / "flnc/classification_summary.json")
    require(report.get("schema_version") == "ploidypatch.apple_pychopper_flnc_summary.v1", "classification schema changed")
    tissues = report.get("tissues")
    require(isinstance(tissues, list) and len(tissues) == 7, "classification must have seven tissues")
    observed = {row["accession"]: row["tissue"] for row in tissues}
    require(observed == EXPECTED_TISSUES, "classification tissue/accession set changed")
    actual_fastq_counts: dict[str, int] = {}
    if count_fastq:
        with ThreadPoolExecutor(max_workers=7) as executor:
            actual_fastq_counts = dict(
                zip(
                    (row["accession"] for row in tissues),
                    executor.map(
                        count_fastq_records,
                        (
                            root / "flnc/by_run" / f"{row['accession']}.flnc.fastq.gz"
                            for row in tissues
                        ),
                    ),
                    strict=True,
                )
            )
    numeric_keys = [
        "input_reads",
        "prefilter_pass_reads",
        "prefilter_length_failed",
        "pychopper_quality_failed",
        "pychopper_primary_classified",
        "pychopper_primary_flnc",
        "pychopper_output_segments_below_50",
        "pychopper_rescued_reads_excluded",
        "pychopper_rescued_segments_excluded",
        "pychopper_unusable_reads",
    ]
    for row in tissues:
        for key in numeric_keys:
            require(isinstance(row[key], int) and row[key] >= 0, f"invalid {key}")
        require(row["input_reads"] == row["prefilter_pass_reads"] + row["prefilter_length_failed"], "prefilter count mismatch")
        require(
            row["prefilter_pass_reads"]
            == row["pychopper_primary_classified"]
            + row["pychopper_rescued_reads_excluded"]
            + row["pychopper_unusable_reads"]
            + row["pychopper_quality_failed"],
            "Pychopper read classification count mismatch",
        )
        require(row["pychopper_primary_flnc"] <= row["pychopper_primary_classified"], "FLNC output exceeds primary classification")
        require(
            row["pychopper_primary_classified"] - row["pychopper_primary_flnc"]
            <= row["pychopper_output_segments_below_50"],
            "unexplained classified/output FLNC discrepancy",
        )
        expected_fraction = row["pychopper_primary_flnc"] / (
            row["prefilter_pass_reads"] - row["pychopper_quality_failed"]
        )
        require(math.isclose(row["primary_flnc_fraction_of_q7_pass"], expected_fraction, rel_tol=1e-12), "FLNC fraction mismatch")
        if count_fastq:
            require(
                actual_fastq_counts[row["accession"]]
                == row["pychopper_primary_flnc"],
                f"actual FLNC count mismatch: {row['accession']}",
            )
    totals = report.get("totals")
    require(isinstance(totals, dict), "classification totals missing")
    for key in numeric_keys:
        require(totals.get(key) == sum(row[key] for row in tissues), f"classification total mismatch: {key}")
    require(report.get("reported_ranges_are_descriptive_not_validity_gates") is True, "publication range became a gate")
    return report


def verify_evidence(root: Path, expected_candidates: int) -> tuple[dict, dict[str, dict], Counter]:
    path = root / "evidence/candidate_flnc_evidence.tsv"
    manifest = load_json(Path(f"{path}.manifest.json"))
    require(manifest.get("schema_version") == "ploidypatch.isoseq_candidate_validation.v2", "evidence schema changed")
    parameters = manifest["parameters"]
    expected_parameters = {
        "alignment_strand_source": "query_orientation",
        "minimum_query_coverage": 0.85,
        "minimum_identity": 0.9,
        "minimum_mapq": 20,
        "maximum_secondary_score_fraction": 0.95,
        "minimum_candidate_cds_coverage": 0.9,
        "flank_bp": 5000,
        "single_exon_is_not_full_chain_positive": True,
    }
    for key, value in expected_parameters.items():
        require(parameters.get(key) == value, f"evidence parameter changed: {key}")
    require(manifest["output"]["rows"] == expected_candidates, "evidence manifest row count changed")
    require(manifest["output"]["sha256"] == sha256_file(path), "evidence output hash mismatch")
    evidence: dict[str, dict] = {}
    states: Counter = Counter()
    for row in read_tsv(path):
        digest = row["candidate_digest"]
        require(digest not in evidence, f"duplicate evidence candidate {digest}")
        require(row["evidence_state"] in EVIDENCE_STATES, f"unknown evidence state {row['evidence_state']}")
        segments = int(row["cds_segments"])
        reads = int(row["supporting_full_length_reads"])
        if row["evidence_state"] == "full_chain_supported":
            require(segments > 1 and reads > 0, f"invalid full-chain positive {digest}")
        if row["evidence_state"] == "single_exon_span_support":
            require(segments == 1, f"single-exon support assigned to multi-exon candidate {digest}")
        evidence[digest] = row
        states[row["evidence_state"]] += 1
    require(len(evidence) == expected_candidates, "evidence candidate universe changed")
    require(dict(states) == manifest["counts"]["evidence_states"], "evidence state counts disagree with manifest")
    require(manifest["counts"]["candidate_models"] == expected_candidates, "candidate count changed")
    return manifest, evidence, states


def verify_ranked_review(root: Path, evidence: dict[str, dict], states: Counter) -> tuple[dict, dict]:
    ranked_path = root / "evidence/ranked_flnc_evidence.tsv"
    by_estimator: dict[str, list[tuple[int, str, str]]] = defaultdict(list)
    for row in read_tsv(ranked_path):
        estimator = row["estimator"]
        require(estimator in EXPECTED_ESTIMATORS, f"unknown estimator {estimator}")
        digest = row["candidate_digest"]
        require(digest in evidence, f"ranked candidate absent from evidence {digest}")
        require(row["evidence_state"] == evidence[digest]["evidence_state"], f"ranked evidence mismatch {digest}")
        by_estimator[estimator].append((int(row["review_rank"]), digest, row["evidence_state"]))
    universe = set(evidence)
    for estimator in EXPECTED_ESTIMATORS:
        values = sorted(by_estimator[estimator])
        require([rank for rank, _, _ in values] == list(range(1, len(universe) + 1)), f"rank sequence changed: {estimator}")
        require({digest for _, digest, _ in values} == universe, f"rank universe changed: {estimator}")

    review = load_json(root / "evidence/review_yield.json")
    require(review.get("schema_version") == "ploidypatch.isoseq_review_yield.v1", "review schema changed")
    require(review["parameters"]["review_budgets"] == EXPECTED_BUDGETS, "review budgets changed")
    require(review["parameters"]["primary_estimator"] == "v04_guard", "primary estimator changed")
    require(review["parameters"]["comparator_estimator"] == "baseline", "comparator changed")
    observed_rows = {(row["estimator"], row["review_budget"]): row for row in review["yield"]}
    require(len(observed_rows) == len(EXPECTED_ESTIMATORS) * len(EXPECTED_BUDGETS), "review-yield grid incomplete")
    for estimator in EXPECTED_ESTIMATORS:
        values = sorted(by_estimator[estimator])
        for budget in EXPECTED_BUDGETS:
            counts = Counter(state for _, _, state in values[:budget])
            observed = observed_rows[(estimator, budget)]
            for state in EVIDENCE_STATES:
                require(observed.get(state, 0) == counts[state], f"review yield mismatch: {estimator}/{budget}/{state}")
    primary = review["primary"]["full_chain_supported"]
    for estimator in EXPECTED_ESTIMATORS:
        require(primary[estimator] == observed_rows[(estimator, 100)]["full_chain_supported"], f"primary count mismatch: {estimator}")
    expected_delta = primary["v04_guard"] - primary["baseline"]
    require(review["primary"]["estimator_delta"]["delta_full_chain_supported"] == expected_delta, "primary delta mismatch")
    baseline_top = {digest for _, digest, _ in sorted(by_estimator["baseline"])[:100]}
    primary_top = {digest for _, digest, _ in sorted(by_estimator["v04_guard"])[:100]}
    require(review["primary"]["top_set_overlap"] == len(baseline_top & primary_top), "top-set overlap mismatch")

    bootstrap = load_json(root / "evidence/bootstrap.json")
    require(bootstrap.get("schema_version") == "ploidypatch.isoseq_review_bootstrap.v2", "bootstrap schema changed")
    require(bootstrap["parameters"]["replicates"] == 20000, "bootstrap replicates changed")
    require(bootstrap["parameters"]["seed"] == 20261004, "bootstrap seed changed")
    require(bootstrap["parameters"]["review_budgets"] == EXPECTED_BUDGETS, "bootstrap budgets changed")
    require(bootstrap["counts"]["candidates"] == len(evidence), "bootstrap candidate count mismatch")
    require(bootstrap["counts"]["full_chain_supported"] == states["full_chain_supported"], "bootstrap positive count mismatch")
    delta = bootstrap["primary_estimator_delta"]
    require(delta["observed_delta_supported_per_100"] == expected_delta, "bootstrap observed delta mismatch")
    require(delta["ci_lower"] <= expected_delta <= delta["ci_upper"], "bootstrap interval excludes observation")
    random_rows = {(row["estimator"], row["review_budget"]): row for row in bootstrap["random_ranking_enrichment"]}
    require(set(random_rows) == set(observed_rows), "random-enrichment grid incomplete")
    for key, row in random_rows.items():
        require(row["observed_full_chain_supported"] == observed_rows[key]["full_chain_supported"], f"random observed count mismatch: {key}")
        require(0 < row["empirical_p_random_ge_observed"] <= 1, f"invalid empirical p-value: {key}")
        require(row["random_mean"] >= 0 and row["fold_enrichment_over_random_mean"] >= 0, f"invalid enrichment: {key}")
    return review, bootstrap


def verify_audit(root: Path, evidence: dict[str, dict], states: Counter) -> tuple[dict, set[str], dict[str, dict[str, int]]]:
    summary = load_json(root / "audit/summary.json")
    require(summary.get("schema_version") == "ploidypatch.natural_candidate_biological_audit.v1", "audit schema changed")
    require(summary["interpretation"]["automatic_approval"] is False, "audit enables automatic approval")
    require(summary["parameters"]["minimum_full_length_read_support"] == 2, "case read threshold changed")
    require(summary["parameters"]["review_budgets"] == EXPECTED_BUDGETS, "audit budgets changed")
    require(summary["counts"]["candidates"] == len(evidence), "audit candidate count mismatch")
    require(summary["counts"]["isoseq_evidence_state"] == dict(states), "audit evidence-state counts mismatch")

    signatures: dict[str, tuple[str, int, int, int, int, int, int, str]] = {}
    estimator_digests: dict[str, set[str]] = defaultdict(set)
    cases: set[str] = set()
    case_ranks: dict[str, dict[str, int]] = defaultdict(dict)
    chain_supported: set[str] = set()
    for row in read_tsv(root / "audit/candidate_biological_audit.tsv"):
        estimator = row["estimator"]
        require(estimator in EXPECTED_ESTIMATORS, f"unknown audit estimator {estimator}")
        digest = row["candidate_digest"]
        require(digest in evidence, f"audit candidate absent from evidence {digest}")
        require(digest not in estimator_digests[estimator], f"duplicate audit row {estimator}/{digest}")
        estimator_digests[estimator].add(digest)
        signature = (
            row["isoseq_evidence_state"],
            int(row["isoseq_supporting_full_length_reads"]),
            int(row["transcript_chain_supported"]),
            int(row["strict_collision_free"]),
            int(row["repeat_risk_flag"]),
            int(row["case_study_ready"]),
            int(row["complete_orf"]),
            row["mappability_class"],
        )
        if digest in signatures:
            require(signatures[digest] == signature, f"audit evidence differs by estimator: {digest}")
        else:
            signatures[digest] = signature
        require(row["isoseq_evidence_state"] == evidence[digest]["evidence_state"], f"audit/evidence mismatch {digest}")
        if int(row["transcript_chain_supported"]):
            chain_supported.add(digest)
            require(row["isoseq_evidence_state"] == "full_chain_supported", f"invalid chain support {digest}")
            require(int(row["isoseq_supporting_full_length_reads"]) >= 2, f"chain support below threshold {digest}")
        if int(row["case_study_ready"]):
            cases.add(digest)
            case_ranks[digest][estimator] = int(row["review_rank"])
            require(int(row["cds_segments"]) > 1, f"single-exon case {digest}")
            require(int(row["complete_orf"]) == 1, f"incomplete ORF case {digest}")
            require(int(row["strict_collision_free"]) == 1, f"collision-positive case {digest}")
            require(
                row["mappability_class"] in {"unique_locus", "duplicated_locus"},
                f"ambiguous-mapping case {digest}",
            )
            require(int(row["repeat_risk_flag"]) == 0, f"repeat-risk case {digest}")
            require(int(row["transcript_chain_supported"]) == 1, f"unsupported case {digest}")
    universe = set(evidence)
    for estimator in EXPECTED_ESTIMATORS:
        require(estimator_digests[estimator] == universe, f"audit universe changed: {estimator}")
    require(len(cases) == summary["counts"]["case_study_ready"], "case count mismatch")
    require(len(chain_supported) == summary["counts"]["transcript_chain_supported"], "chain-supported count mismatch")
    require(summary["counts"]["ranking_rows"] == len(evidence) * len(EXPECTED_ESTIMATORS), "audit row count mismatch")
    return summary, cases, case_ranks


def verify_comparison(root: Path, evidence: dict[str, dict], cases: set[str], review: dict) -> dict:
    comparison = load_json(root / "evidence/raw_ts_patch2_vs_flnc.json")
    counts = comparison["counts"]
    full = {digest for digest, row in evidence.items() if row["evidence_state"] == "full_chain_supported"}
    reads2 = {digest for digest in full if int(evidence[digest]["supporting_full_length_reads"]) >= 2}
    require(counts["candidate_universe"] == len(evidence), "comparison candidate count mismatch")
    require(counts["pychopper_flnc_full_chain"] == len(full), "comparison FLNC full-chain count mismatch")
    require(counts["pychopper_flnc_full_chain_reads_ge_2"] == len(reads2), "comparison FLNC reads>=2 mismatch")
    require(counts["pychopper_flnc_case_study_ready"] == len(cases), "comparison case count mismatch")
    require(sum(comparison["state_transitions"].values()) == len(evidence), "comparison transition count mismatch")
    require(comparison["review_primary"]["pychopper_flnc"] == review["primary"], "comparison review-primary mismatch")
    return comparison


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("result_root", type=Path)
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--expected-code-commit", required=True)
    parser.add_argument("--expected-candidates", type=int, default=29144)
    parser.add_argument("--fast", action="store_true", help="skip large file/input hashes and FASTQ recount")
    parser.add_argument("--output-json", type=Path)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    root = args.result_root.resolve(strict=True)
    project_root = args.project_root.resolve(strict=True)
    require(root.is_dir() and not root.name.endswith(".working"), "result root is not a finalized directory")
    require(re.fullmatch(r"[0-9a-f]{40}", args.expected_code_commit) is not None, "expected commit is not a full SHA")
    manifest_files, manifest_sha = verify_root_manifest(root, verify_hashes=not args.fast)
    verify_read_only_tree(root)
    verify_referenced_manifests(root, project_root, verify_hashes=not args.fast)
    contract = verify_contract(root, args.expected_code_commit)
    classification = verify_classification(root, count_fastq=not args.fast)
    evidence_manifest, evidence, states = verify_evidence(root, args.expected_candidates)
    review, bootstrap = verify_ranked_review(root, evidence, states)
    audit, cases, case_ranks = verify_audit(root, evidence, states)
    comparison = verify_comparison(root, evidence, cases, review)
    report = {
        "schema_version": "ploidypatch.apple_flnc_independent_verification.v1",
        "valid": True,
        "verification_mode": "fast" if args.fast else "full",
        "verifier": {
            "path": os.fspath(Path(__file__).resolve()),
            "sha256": sha256_file(Path(__file__).resolve()),
        },
        "result_root": os.fspath(root),
        "code_commit": contract["code_commit"],
        "sha256sums_sha256": manifest_sha,
        "manifest_files": manifest_files,
        "counts": {
            "input_reads": classification["totals"]["input_reads"],
            "primary_flnc_reads": classification["totals"]["pychopper_primary_flnc"],
            "candidates": len(evidence),
            "full_chain_supported": states["full_chain_supported"],
            "full_chain_reads_ge_2": audit["counts"]["transcript_chain_supported"],
            "case_study_ready": len(cases),
        },
        "primary_review": review["primary"],
        "primary_bootstrap": bootstrap["primary_estimator_delta"],
        "raw_ts_vs_flnc_counts": comparison["counts"],
        "case_ranks": {digest: case_ranks[digest] for digest in sorted(cases)},
        "evidence_manifest_sha256": sha256_file(root / "evidence/candidate_flnc_evidence.tsv.manifest.json"),
    }
    serialized = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output_json:
        args.output_json.parent.mkdir(parents=True, exist_ok=True)
        with args.output_json.open("x", encoding="utf-8", newline="") as handle:
            handle.write(serialized)
    else:
        sys.stdout.write(serialized)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except VerificationError as error:
        print(f"verification failed: {error}", file=sys.stderr)
        raise SystemExit(1) from error
