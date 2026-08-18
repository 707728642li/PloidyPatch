#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any


SCHEMA_VERSION = "ploidypatch.conflict_winner_guard_production_replay.v1"
GUARD_MANIFEST_SCHEMA = "ploidypatch.conflict_winner_guard_scores.v1"
VALUE_FIELDS = {
    "baseline": "v03_baseline_logit",
    "primary": "v03_primary_rank_score",
    "v04_guard": "v04_primary_rank_score",
}
FLAG_FIELDS = {
    "v04_guard_applied": "v04_conflict_guard_applied",
    "v04_topology_abstained": "v04_topology_abstained",
    "v04_automatic_approval": "v04_automatic_approval",
}
MAX_ABSOLUTE_DIFFERENCE = 1e-12
MAX_ULP_DIFFERENCE = 32.0


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def parse_production(value: str) -> tuple[str, Path]:
    name, separator, path_value = value.partition("=")
    path = Path(path_value)
    if not separator or not name or not path.is_file() or path.stat().st_size == 0:
        raise ValueError("--production requires NAME=NONEMPTY_V04_TSV")
    return name, path


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Missing TSV header: {path}")
        rows = list(reader)
    return list(reader.fieldnames), rows


def indexed(rows: list[dict[str, str]], path: Path) -> dict[str, dict[str, str]]:
    output = {row.get("candidate_digest", ""): row for row in rows}
    if "" in output or len(output) != len(rows):
        raise ValueError(f"Empty or duplicate candidate digest: {path}")
    return output


def parse_float(value: str, field: str, digest: str) -> float:
    parsed = float(value)
    if not math.isfinite(parsed):
        raise ValueError(f"Non-finite {field} for {digest}")
    return parsed


def ulp_difference(left: float, right: float) -> float:
    if left == right:
        return 0.0
    return abs(left - right) / max(math.ulp(left), math.ulp(right))


def ranking(rows: dict[str, dict[str, str]], field: str) -> list[str]:
    return sorted(
        rows,
        key=lambda digest: (
            -parse_float(rows[digest][field], field, digest),
            digest,
        ),
    )


def verify_dataset(
    name: str,
    production_path: Path,
    development_rows: list[dict[str, str]],
) -> dict[str, Any]:
    manifest_path = Path(str(production_path) + ".manifest.json")
    if not manifest_path.is_file() or manifest_path.stat().st_size == 0:
        raise ValueError(f"Missing production guard manifest: {manifest_path}")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if (
        manifest.get("schema_version") != GUARD_MANIFEST_SCHEMA
        or manifest.get("truth_access") is not False
        or manifest.get("outputs", {}).get("scores", {}).get("sha256")
        != sha256(production_path)
        or manifest.get("winner_audit", {}).get("mismatch_count") != 0
        or manifest.get("winner_audit", {}).get("baseline_mapping_sha256")
        != manifest.get("winner_audit", {}).get("v04_guard_mapping_sha256")
    ):
        raise ValueError(f"Production guard manifest fails safety gate: {name}")

    production_fields, production_rows = read_tsv(production_path)
    forbidden = [
        field
        for field in production_fields
        if "label" in field.lower() or "truth" in field.lower()
    ]
    required = {"candidate_digest", *VALUE_FIELDS.values(), *FLAG_FIELDS.values()}
    if forbidden or not required <= set(production_fields):
        raise ValueError(f"Production replay contains truth or lacks fields: {name}")
    development = indexed(development_rows, Path(f"development:{name}"))
    production = indexed(production_rows, production_path)
    if set(development) != set(production):
        raise ValueError(f"Candidate universe mismatch in production replay: {name}")

    tolerance_violation_count = 0
    flag_mismatch_count = 0
    bitwise_difference_count = 0
    maximum_absolute_difference = 0.0
    maximum_ulp_difference = 0.0
    canonical_lines: list[str] = []
    for digest in sorted(development):
        left = development[digest]
        right = production[digest]
        values: list[str] = []
        for development_field, production_field in VALUE_FIELDS.items():
            left_value = parse_float(
                left[development_field], development_field, digest
            )
            right_value = parse_float(
                right[production_field], production_field, digest
            )
            difference = abs(left_value - right_value)
            ulps = ulp_difference(left_value, right_value)
            bitwise_difference_count += int(left_value.hex() != right_value.hex())
            tolerance_violation_count += int(
                difference > MAX_ABSOLUTE_DIFFERENCE and ulps > MAX_ULP_DIFFERENCE
            )
            maximum_absolute_difference = max(maximum_absolute_difference, difference)
            maximum_ulp_difference = max(maximum_ulp_difference, ulps)
            values.extend((left_value.hex(), right_value.hex()))
        for development_field, production_field in FLAG_FIELDS.items():
            left_value = left[development_field]
            right_value = right[production_field]
            flag_mismatch_count += int(left_value != right_value)
            values.extend((left_value, right_value))
        canonical_lines.append("\t".join((digest, *values)) + "\n")
    ranking_mismatch_fields = []
    for development_field, production_field in VALUE_FIELDS.items():
        if ranking(development, development_field) != ranking(
            production, production_field
        ):
            ranking_mismatch_fields.append(development_field)
    mismatch_count = (
        tolerance_violation_count
        + flag_mismatch_count
        + len(ranking_mismatch_fields)
    )
    if mismatch_count:
        raise AssertionError(
            f"Production replay differs for {name}: tolerance="
            f"{tolerance_violation_count}, flags={flag_mismatch_count}, "
            f"ranking={ranking_mismatch_fields}, max_abs="
            f"{maximum_absolute_difference:.17g}, max_ulp="
            f"{maximum_ulp_difference:.17g}"
        )
    return {
        "candidates": len(production),
        "field_comparisons_per_candidate": len(VALUE_FIELDS) + len(FLAG_FIELDS),
        "mismatch_count": 0,
        "score_tolerance": {
            "acceptance_rule": (
                "absolute_difference_at_most_threshold_or_ulp_difference_"
                "at_most_threshold"
            ),
            "maximum_allowed_absolute_difference": MAX_ABSOLUTE_DIFFERENCE,
            "maximum_allowed_ulp_difference": MAX_ULP_DIFFERENCE,
            "maximum_observed_absolute_difference": maximum_absolute_difference,
            "maximum_observed_ulp_difference": maximum_ulp_difference,
            "bitwise_difference_count": bitwise_difference_count,
            "tolerance_violation_count": tolerance_violation_count,
        },
        "flag_mismatch_count": flag_mismatch_count,
        "ranking_mismatch_fields": ranking_mismatch_fields,
        "comparison_sha256": hashlib.sha256(
            "".join(canonical_lines).encode("utf-8")
        ).hexdigest(),
        "production_scores": {
            "path": str(production_path),
            "bytes": production_path.stat().st_size,
            "sha256": sha256(production_path),
        },
        "production_manifest": {
            "path": str(manifest_path),
            "sha256": sha256(manifest_path),
        },
        "winner_audit": manifest["winner_audit"],
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Verify production v0.4 scores against frozen development predictions"
    )
    parser.add_argument("--development-predictions", required=True)
    parser.add_argument("--production", action="append", required=True)
    parser.add_argument("--code-commit", required=True)
    parser.add_argument("--output-json", required=True)
    args = parser.parse_args()
    development_path = Path(args.development_predictions)
    if not development_path.is_file() or development_path.stat().st_size == 0:
        raise ValueError("Missing development predictions")
    _, development_rows = read_tsv(development_path)
    if not {"dataset", "candidate_digest", *VALUE_FIELDS, *FLAG_FIELDS} <= set(
        development_rows[0]
    ):
        raise ValueError("Development predictions lack replay fields")
    production_inputs = [parse_production(value) for value in args.production]
    names = [name for name, _ in production_inputs]
    if len(names) != len(set(names)):
        raise ValueError("Production dataset names must be unique")
    reports = {}
    for name, path in production_inputs:
        selected = [row for row in development_rows if row["dataset"] == name]
        if not selected:
            raise ValueError(f"No development predictions for {name}")
        reports[name] = verify_dataset(name, path, selected)
    output = Path(args.output_json)
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite replay audit: {output}")
    report = {
        "schema_version": SCHEMA_VERSION,
        "code_commit": args.code_commit,
        "truth_access": "seen_development_predictions_only",
        "labels_used": False,
        "development_predictions": {
            "path": str(development_path),
            "bytes": development_path.stat().st_size,
            "sha256": sha256(development_path),
        },
        "datasets": reports,
        "all_replay_gates_pass": all(
            item["mismatch_count"] == 0 for item in reports.values()
        ),
        "all_rankings_exact": all(
            not item["ranking_mismatch_fields"] for item in reports.values()
        ),
        "all_scores_bitwise_exact": all(
            item["score_tolerance"]["bitwise_difference_count"] == 0
            for item in reports.values()
        ),
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x", encoding="utf-8", newline="") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
