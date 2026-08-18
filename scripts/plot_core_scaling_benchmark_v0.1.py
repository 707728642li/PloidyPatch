#!/usr/bin/env python3
"""Render the preregistered PloidyPatch core-scaling figure from frozen tables."""

from __future__ import annotations

import argparse
import csv
import datetime
import os
import platform
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

matplotlib.rcParams["svg.hashsalt"] = "ploidypatch-core-scaling-v0.1"

from core_scaling_common_v0_1 import (  # noqa: E402
    COMPONENTS,
    ScalingError,
    canonical_json,
    freeze_read_only,
    load_json,
    require,
    sha256_file,
    verify_sha256sums,
    write_sha256sums,
)


COLORS = {
    "consensus": "#0072B2",
    "feature_construction": "#D55E00",
    "score_v03": "#009E73",
    "guard_v04": "#CC79A7",
    "end_to_end": "#000000",
}
LABELS = {
    "consensus": "Chain-preserving consensus",
    "feature_construction": "Plant feature construction",
    "score_v03": "v0.3 ranking",
    "guard_v04": "v0.4 safety guard",
    "end_to_end": "End-to-end core",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary-root", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require(reader.fieldnames is not None, f"missing TSV header: {path}")
        return list(reader)


def _ordered(rows: list[dict[str, str]], component: str) -> list[dict[str, str]]:
    selected = [row for row in rows if row["component"] == component]
    return sorted(selected, key=lambda row: int(row["candidate_count"]))


def main() -> int:
    args = parse_args()
    summary_root = args.summary_root.resolve(strict=True)
    output = args.output_dir.resolve()
    working = Path(str(output) + ".working")
    require(not output.exists() and not working.exists(), f"refusing to overwrite: {output}")
    verify_sha256sums(summary_root)
    summary = load_json(summary_root / "summary.json")
    require(summary.get("operational_gate_pass") is True, "scaling gate did not pass")
    rows = read_rows(summary_root / "figure_source_data.tsv")
    require(
        len(rows) == 6 * len(COMPONENTS),
        "figure source does not contain six workloads and five components",
    )
    require(
        all(int(row["warm_pass"]) == 3 and row["warm_wall_median"] for row in rows),
        "figure source has incomplete warm repetitions",
    )
    species_order = [
        row["workload_id"] for row in _ordered(rows, "end_to_end")
    ]
    species_labels = {
        row["workload_id"]: row["species"].replace(" ", "\n", 1)
        for row in rows
    }
    fig, axes = plt.subplots(2, 2, figsize=(11.2, 8.2), constrained_layout=True)
    ax_wall, ax_rss, ax_share, ax_cold = axes.ravel()

    for component in COMPONENTS:
        selected = _ordered(rows, component)
        x = [int(row["candidate_count"]) for row in selected]
        wall = [float(row["warm_wall_median"]) for row in selected]
        wall_low = [float(row["warm_wall_min"]) for row in selected]
        wall_high = [float(row["warm_wall_max"]) for row in selected]
        rss = [float(row["warm_peak_rss_median_kb"]) / (1024 * 1024) for row in selected]
        rss_low = [float(row["warm_peak_rss_min_kb"]) / (1024 * 1024) for row in selected]
        rss_high = [float(row["warm_peak_rss_max_kb"]) / (1024 * 1024) for row in selected]
        linewidth = 2.2 if component == "end_to_end" else 1.3
        markersize = 6.0 if component == "end_to_end" else 4.5
        ax_wall.errorbar(
            x,
            wall,
            yerr=[
                [median - low for median, low in zip(wall, wall_low, strict=True)],
                [high - median for median, high in zip(wall, wall_high, strict=True)],
            ],
            color=COLORS[component],
            marker="o",
            markersize=markersize,
            linewidth=linewidth,
            capsize=2,
            label=LABELS[component],
        )
        ax_rss.errorbar(
            x,
            rss,
            yerr=[
                [median - low for median, low in zip(rss, rss_low, strict=True)],
                [high - median for median, high in zip(rss, rss_high, strict=True)],
            ],
            color=COLORS[component],
            marker="o",
            markersize=markersize,
            linewidth=linewidth,
            capsize=2,
        )

    for axis, ylabel in (
        (ax_wall, "Wall time (s)"),
        (ax_rss, "Peak RSS (GiB)"),
    ):
        axis.set_xscale("log")
        axis.set_yscale("log")
        axis.set_xlabel("Accepted candidate chains")
        axis.set_ylabel(ylabel)
        axis.grid(True, which="both", color="#dddddd", linewidth=0.6)
    ax_wall.legend(frameon=False, fontsize=8, ncol=2)

    largest = max(
        (row for row in rows if row["component"] == "end_to_end"),
        key=lambda row: int(row["candidate_count"]),
    )["workload_id"]
    isolated = [component for component in COMPONENTS if component != "end_to_end"]
    shares = {
        row["component"]: 100.0 * float(row["component_share_of_end_to_end"])
        for row in rows
        if row["workload_id"] == largest and row["component"] in isolated
    }
    y_positions = list(range(len(isolated)))
    ax_share.barh(
        y_positions,
        [shares[component] for component in isolated],
        color=[COLORS[component] for component in isolated],
    )
    ax_share.set_yticks(y_positions, [LABELS[component] for component in isolated])
    ax_share.invert_yaxis()
    ax_share.set_xlabel("Isolated component / end-to-end wall time (%)")
    ax_share.set_title(f"Largest workload: {species_labels[largest].replace(chr(10), ' ')}")
    ax_share.grid(True, axis="x", color="#dddddd", linewidth=0.6)

    end_rows = {row["workload_id"]: row for row in rows if row["component"] == "end_to_end"}
    ratios = [
        float(end_rows[workload]["cold_wall_seconds"])
        / float(end_rows[workload]["warm_wall_median"])
        for workload in species_order
    ]
    ax_cold.bar(
        list(range(len(species_order))),
        ratios,
        color=COLORS["end_to_end"],
        width=0.72,
    )
    ax_cold.axhline(1.0, color="#888888", linewidth=1.0, linestyle="--")
    ax_cold.set_xticks(
        list(range(len(species_order))),
        [species_labels[item] for item in species_order],
        rotation=25,
        ha="right",
    )
    ax_cold.set_ylabel("First run / warm median")
    ax_cold.set_title("Cache-order sensitivity")
    ax_cold.grid(True, axis="y", color="#dddddd", linewidth=0.6)

    for label, axis in zip(("a", "b", "c", "d"), axes.ravel(), strict=True):
        axis.text(-0.12, 1.06, label, transform=axis.transAxes, fontsize=14, fontweight="bold", va="top")
    fig.suptitle("PloidyPatch core scaling on six real plant workloads", fontsize=14)
    working.mkdir(parents=True)
    outputs = {
        "png": working / "core_scaling_v0.1.png",
        "pdf": working / "core_scaling_v0.1.pdf",
        "svg": working / "core_scaling_v0.1.svg",
    }
    os.environ["SOURCE_DATE_EPOCH"] = "0"
    fixed_date = datetime.datetime(1970, 1, 1, tzinfo=datetime.timezone.utc)
    fig.savefig(
        outputs["png"],
        dpi=300,
        facecolor="white",
        metadata={"Software": "PloidyPatch"},
    )
    fig.savefig(
        outputs["pdf"],
        facecolor="white",
        metadata={
            "Creator": "PloidyPatch",
            "CreationDate": fixed_date,
            "ModDate": fixed_date,
        },
    )
    fig.savefig(
        outputs["svg"],
        facecolor="white",
        metadata={"Creator": "PloidyPatch", "Date": "1970-01-01T00:00:00Z"},
    )
    plt.close(fig)
    manifest = {
        "schema_version": "ploidypatch.core_scaling_figure.v1",
        "summary_sha256sums_sha256": sha256_file(summary_root / "SHA256SUMS"),
        "source_data_sha256": sha256_file(summary_root / "figure_source_data.tsv"),
        "matplotlib_version": matplotlib.__version__,
        "python_version": platform.python_version(),
        "panels": {
            "a": "warm wall median and full three-repeat range versus candidate count",
            "b": "warm peak RSS median and full three-repeat range versus candidate count",
            "c": "isolated component time divided by end-to-end time for the largest workload",
            "d": "first-order end-to-end wall time divided by warm median",
        },
        "thread_scaling": "not_applicable_core_components_expose_no_thread_parameter",
        "outputs": {
            key: {"file_name": path.name, "bytes": path.stat().st_size, "sha256": sha256_file(path)}
            for key, path in outputs.items()
        },
    }
    (working / "figure_manifest.json").write_text(canonical_json(manifest), encoding="utf-8")
    write_sha256sums(working)
    working.rename(output)
    freeze_read_only(output)
    print(canonical_json({"figure_root": str(output), **manifest}), end="")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except ScalingError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
