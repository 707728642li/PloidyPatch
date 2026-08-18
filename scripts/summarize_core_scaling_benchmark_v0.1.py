#!/usr/bin/env python3
"""Summarize a frozen PloidyPatch core-scaling benchmark without cherry-picking."""

from __future__ import annotations

import argparse
import csv
import math
import statistics
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable

from core_scaling_common_v0_1 import (
    COMPONENTS,
    REFERENCE_SCHEMA,
    RESULT_SCHEMA,
    ScalingError,
    canonical_json,
    freeze_read_only,
    load_json,
    require,
    sha256_file,
    verify_sha256sums,
    write_sha256sums,
    write_tsv,
)


SUMMARY_SCHEMA = "ploidypatch.core_scaling_summary.v1"
T975 = {
    1: 12.7062,
    2: 4.3027,
    3: 3.1824,
    4: 2.7764,
    5: 2.5706,
    6: 2.4469,
    7: 2.3646,
    8: 2.3060,
    9: 2.2622,
    10: 2.2281,
    11: 2.2010,
    12: 2.1788,
    13: 2.1604,
    14: 2.1448,
    15: 2.1314,
    16: 2.1199,
    17: 2.1098,
    18: 2.1009,
    19: 2.0930,
    20: 2.0860,
    21: 2.0796,
    22: 2.0739,
    23: 2.0687,
    24: 2.0639,
    25: 2.0595,
    26: 2.0555,
    27: 2.0518,
    28: 2.0484,
    29: 2.0452,
    30: 2.0423,
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--benchmark-root", required=True, type=Path)
    parser.add_argument("--reference-root", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require(reader.fieldnames is not None, f"missing TSV header: {path}")
        return list(reader)


def ols_loglog(points: list[tuple[str, float, float]]) -> tuple[dict[str, float], list[dict[str, Any]]]:
    require(len(points) >= 3, "log-log regression requires at least three workloads")
    names = [item[0] for item in points]
    xs = [math.log10(item[1]) for item in points]
    ys = [math.log10(item[2]) for item in points]
    x_mean = statistics.mean(xs)
    y_mean = statistics.mean(ys)
    sxx = sum((value - x_mean) ** 2 for value in xs)
    require(sxx > 0, "regression predictor has zero variance")
    slope = sum((x - x_mean) * (y - y_mean) for x, y in zip(xs, ys, strict=True)) / sxx
    intercept = y_mean - slope * x_mean
    predictions = [intercept + slope * value for value in xs]
    residuals = [observed - predicted for observed, predicted in zip(ys, predictions, strict=True)]
    residual_ss = sum(value * value for value in residuals)
    total_ss = sum((value - y_mean) ** 2 for value in ys)
    degrees = len(points) - 2
    residual_variance = residual_ss / degrees
    slope_se = math.sqrt(residual_variance / sxx)
    tcrit = T975[degrees] if degrees <= 30 else 1.96
    summary = {
        "n": float(len(points)),
        "slope": slope,
        "slope_se": slope_se,
        "slope_ci_lower": slope - tcrit * slope_se,
        "slope_ci_upper": slope + tcrit * slope_se,
        "intercept": intercept,
        "r_squared": 1.0 - residual_ss / total_ss if total_ss else 1.0,
        "residual_ss": residual_ss,
    }
    rows = [
        {
            "workload_id": name,
            "log10_candidate_count": x,
            "log10_observed": y,
            "log10_predicted": predicted,
            "log10_residual": residual,
        }
        for name, x, y, predicted, residual in zip(
            names, xs, ys, predictions, residuals, strict=True
        )
    ]
    return summary, rows


def leave_one_out_slopes(points: list[tuple[str, float, float]]) -> dict[str, float]:
    output: dict[str, float] = {}
    for omitted, _, _ in points:
        subset = [point for point in points if point[0] != omitted]
        output[omitted] = ols_loglog(subset)[0]["slope"]
    return output


def numeric(rows: Iterable[dict[str, str]], field: str) -> list[float]:
    return [float(row[field]) for row in rows]


def main() -> int:
    args = parse_args()
    benchmark = args.benchmark_root.resolve(strict=True)
    reference = args.reference_root.resolve(strict=True)
    output = args.output_dir.resolve()
    working = Path(str(output) + ".working")
    require(not output.exists() and not working.exists(), f"refusing to overwrite: {output}")
    verify_sha256sums(benchmark)
    verify_sha256sums(reference)
    contract = load_json(benchmark / "run_contract.json")
    benchmark_summary = load_json(benchmark / "summary.json")
    registry = load_json(reference / "reference_registry.json")
    require(contract.get("schema_version") == RESULT_SCHEMA, "invalid benchmark contract schema")
    require(benchmark_summary.get("schema_version") == RESULT_SCHEMA, "invalid benchmark summary schema")
    require(registry.get("schema_version") == REFERENCE_SCHEMA, "invalid reference registry schema")
    require(contract.get("warm_repeats") == 3, "benchmark did not use three warm repeats")
    require(contract.get("seed") == 20261006, "benchmark seed differs")
    require(
        contract.get("reference_registry_sha256") == sha256_file(reference / "reference_registry.json"),
        "benchmark/reference binding differs",
    )
    metadata = {item["workload_id"]: item for item in registry["workloads"]}
    raw_rows = read_tsv(benchmark / "replicates.tsv")
    expected_keys = {
        (phase, repeat, workload_id, component)
        for phase, repeats in (("cold_order_first", [0]), ("warm_repeats", [1, 2, 3]))
        for repeat in repeats
        for workload_id in metadata
        for component in COMPONENTS
    }
    observed_keys = {
        (row["phase"], int(row["repeat"]), row["workload_id"], row["component"])
        for row in raw_rows
    }
    exact_universe = len(observed_keys) == len(raw_rows) and observed_keys == expected_keys
    failures = [row for row in raw_rows if row["status"] != "pass"]
    grouped: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in raw_rows:
        require(row["workload_id"] in metadata, f"unknown workload row: {row['workload_id']}")
        require(row["component"] in COMPONENTS, f"unknown component row: {row['component']}")
        grouped[(row["workload_id"], row["component"])].append(row)
    summary_rows: list[dict[str, Any]] = []
    for workload_id in metadata:
        for component in COMPONENTS:
            rows = grouped[(workload_id, component)]
            cold = [row for row in rows if row["phase"] == "cold_order_first" and row["status"] == "pass"]
            warm = [row for row in rows if row["phase"] == "warm_repeats" and row["status"] == "pass"]
            candidate_count = int(metadata[workload_id]["pool"]["candidate_count"])
            conflict_sets = int(metadata[workload_id]["pool"]["conflict_sets"])
            warm_wall = numeric(warm, "wall_seconds")
            warm_rss = numeric(warm, "peak_rss_kb")
            warm_cpu = [
                float(row["user_seconds"]) + float(row["system_seconds"])
                for row in warm
            ]
            digests = {row["output_digest"] for row in rows if row["status"] == "pass"}
            summary_rows.append(
                {
                    "workload_id": workload_id,
                    "species": metadata[workload_id]["species"],
                    "genome_bp": metadata[workload_id]["genome"]["bases"],
                    "candidate_count": candidate_count,
                    "conflict_sets": conflict_sets,
                    "component": component,
                    "cold_pass": len(cold),
                    "cold_wall_seconds": "" if not cold else format(float(cold[0]["wall_seconds"]), ".9f"),
                    "cold_peak_rss_kb": "" if not cold else cold[0]["peak_rss_kb"],
                    "warm_pass": len(warm),
                    "warm_wall_median": "" if not warm else format(statistics.median(warm_wall), ".9f"),
                    "warm_wall_min": "" if not warm else format(min(warm_wall), ".9f"),
                    "warm_wall_max": "" if not warm else format(max(warm_wall), ".9f"),
                    "warm_cpu_median": "" if not warm else format(statistics.median(warm_cpu), ".9f"),
                    "warm_peak_rss_median_kb": "" if not warm else format(statistics.median(warm_rss), ".3f"),
                    "warm_peak_rss_min_kb": "" if not warm else format(min(warm_rss), ".3f"),
                    "warm_peak_rss_max_kb": "" if not warm else format(max(warm_rss), ".3f"),
                    "warm_candidates_per_second_median": "" if not warm else format(statistics.median(numeric(warm, "candidates_per_second")), ".9f"),
                    "output_bytes": "" if not warm else warm[0]["output_bytes"],
                    "output_digest_count": len(digests),
                }
            )
    lookup = {(row["workload_id"], row["component"]): row for row in summary_rows}
    for row in summary_rows:
        end_to_end = lookup[(row["workload_id"], "end_to_end")]
        if row["component"] == "end_to_end" or not row["warm_wall_median"] or not end_to_end["warm_wall_median"]:
            row["component_share_of_end_to_end"] = ""
        else:
            row["component_share_of_end_to_end"] = format(
                float(row["warm_wall_median"]) / float(end_to_end["warm_wall_median"]),
                ".9f",
            )
    regression_rows: list[dict[str, Any]] = []
    point_rows: list[dict[str, Any]] = []
    for component in COMPONENTS:
        eligible = [
            row
            for row in summary_rows
            if row["component"] == component and row["warm_pass"] == 3
        ]
        for endpoint, field in (
            ("wall_seconds", "warm_wall_median"),
            ("peak_rss_kb", "warm_peak_rss_median_kb"),
        ):
            points = [
                (row["workload_id"], float(row["candidate_count"]), float(row[field]))
                for row in eligible
                if row[field] and float(row[field]) > 0
            ]
            if len(points) < 3:
                continue
            regression, residuals = ols_loglog(points)
            loo = leave_one_out_slopes(points)
            regression_rows.append(
                {
                    "component": component,
                    "endpoint": endpoint,
                    **{key: format(value, ".12g") for key, value in regression.items()},
                    "leave_one_out_slope_min": format(min(loo.values()), ".12g"),
                    "leave_one_out_slope_max": format(max(loo.values()), ".12g"),
                }
            )
            for row in residuals:
                point_rows.append(
                    {
                        "component": component,
                        "endpoint": endpoint,
                        **{key: format(value, ".12g") if isinstance(value, float) else value for key, value in row.items()},
                        "leave_one_out_slope": format(loo[row["workload_id"]], ".12g"),
                    }
                )
    working.mkdir(parents=True)
    summary_fields = list(summary_rows[0])
    write_tsv(working / "scaling_summary.tsv", summary_fields, summary_rows)
    write_tsv(working / "figure_source_data.tsv", summary_fields, summary_rows)
    regression_fields = [
        "component", "endpoint", "n", "slope", "slope_se", "slope_ci_lower",
        "slope_ci_upper", "intercept", "r_squared", "residual_ss",
        "leave_one_out_slope_min", "leave_one_out_slope_max",
    ]
    write_tsv(working / "loglog_regressions.tsv", regression_fields, regression_rows)
    write_tsv(
        working / "loglog_regression_points.tsv",
        [
            "component", "endpoint", "workload_id", "log10_candidate_count",
            "log10_observed", "log10_predicted", "log10_residual", "leave_one_out_slope",
        ],
        point_rows,
    )
    write_tsv(
        working / "failed_replicates.tsv",
        list(raw_rows[0]),
        failures,
    )
    complete_groups = all(
        row["cold_pass"] == 1
        and row["warm_pass"] == 3
        and row["output_digest_count"] == 1
        for row in summary_rows
    )
    result = {
        "schema_version": SUMMARY_SCHEMA,
        "benchmark_sha256sums_sha256": sha256_file(benchmark / "SHA256SUMS"),
        "reference_sha256sums_sha256": sha256_file(reference / "SHA256SUMS"),
        "expected_replicates": len(expected_keys),
        "observed_replicates": len(raw_rows),
        "exact_replicate_universe": exact_universe,
        "failed_replicates": len(failures),
        "all_groups_complete_and_deterministic": complete_groups,
        "operational_gate_pass": exact_universe and not failures and complete_groups,
        "regressions": len(regression_rows),
        "thread_scaling": "not_applicable_core_components_expose_no_thread_parameter",
        "interpretation": "empirical_scaling_not_asymptotic_complexity",
        "selection_policy": "cold_first_plus_median_and_range_of_all_three_warm_repeats",
    }
    (working / "summary.json").write_text(canonical_json(result), encoding="utf-8")
    write_sha256sums(working)
    working.rename(output)
    freeze_read_only(output)
    print(canonical_json({"summary_root": str(output), **result}), end="")
    return 0 if result["operational_gate_pass"] else 1


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except ScalingError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
