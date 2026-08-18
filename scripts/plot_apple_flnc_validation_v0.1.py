#!/usr/bin/env python3
"""Render the checksum-bound Golden Delicious FLNC validation figure."""

from __future__ import annotations

import argparse
import csv
import datetime
import hashlib
import json
import os
import re
import shutil
import stat
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


matplotlib.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 8.5,
        "svg.hashsalt": "ploidypatch-apple-flnc-v0.1",
    }
)

BLUE = "#0072B2"
ORANGE = "#D55E00"
GREEN = "#009E73"
PURPLE = "#7A5195"
GRAY = "#666666"
LIGHT_GRAY = "#DDDDDD"
ESTIMATOR_COLORS = {
    "baseline": GRAY,
    "v03_primary": GREEN,
    "v04_guard": PURPLE,
}
ESTIMATOR_LABELS = {
    "baseline": "Frozen baseline",
    "v03_primary": "v0.3",
    "v04_guard": "v0.4 guard",
}
SHA_LINE = re.compile(r"^([0-9a-f]{64}) [ *](\./[^/].*)$")


class FigureError(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise FigureError(message)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        value = json.load(handle)
    require(isinstance(value, dict), f"JSON root is not an object: {path}")
    return value


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require(reader.fieldnames is not None, f"TSV header is missing: {path}")
        return list(reader)


def verify_source(source: Path, expected_sha256: str) -> dict[str, Any]:
    require(re.fullmatch(r"[0-9a-f]{64}", expected_sha256) is not None, "expected SHA-256 is invalid")
    checksum = source / "SHA256SUMS"
    require(checksum.is_file(), "source SHA256SUMS is missing")
    require(sha256_file(checksum) == expected_sha256, "source SHA256SUMS hash differs")
    declared: set[str] = set()
    for line_number, line in enumerate(checksum.read_text(encoding="utf-8").splitlines(), start=1):
        match = SHA_LINE.fullmatch(line)
        require(match is not None, f"malformed checksum line {line_number}")
        relative = match.group(2)[2:]
        require("/" not in relative and "\\" not in relative, "nested checksum path is forbidden")
        require(relative not in declared, f"duplicate checksum path: {relative}")
        declared.add(relative)
        path = source / relative
        require(path.is_file() and sha256_file(path) == match.group(1), f"source checksum mismatch: {relative}")
    actual = {
        path.name
        for path in source.iterdir()
        if path.is_file() and path.name != "SHA256SUMS"
    }
    require(actual == declared, "source file universe differs from SHA256SUMS")
    required = {
        "summary.json",
        "manifest.json",
        "table_pychopper_classification.tsv",
        "figure_review_yield.tsv",
        "figure_random_enrichment.tsv",
        "figure_raw_ts_to_flnc_transitions.tsv",
        "figure_primary_rank_delta.tsv",
        "table_strict_natural_cases.tsv",
    }
    require(required <= declared, "source checksum universe is incomplete")
    summary = load_json(source / "summary.json")
    require(
        summary.get("schema_version") == "ploidypatch.apple_flnc_manuscript_source.v1",
        "source summary schema differs",
    )
    require(summary.get("interpretation") == "descriptive_posthoc_natural_validation", "source role differs")
    require(summary.get("automatic_annotation_patch") is False, "natural validation cannot approve patches")
    return summary


def panel_label(ax, label: str) -> None:
    ax.text(-0.12, 1.07, label, transform=ax.transAxes, fontsize=14, fontweight="bold", va="top")


def draw_classification(ax, rows: list[dict[str, str]], plotted: list[dict[str, Any]]) -> None:
    require(len(rows) == 7, "classification must contain seven tissues")
    labels = [row["tissue"].replace("_", " ") for row in rows]
    fractions = [float(row["primary_flnc_fraction_of_q7_pass"]) for row in rows]
    colors = [BLUE if index % 2 == 0 else "#56B4E9" for index in range(len(rows))]
    positions = list(range(len(rows)))
    ax.bar(positions, fractions, color=colors, width=0.72)
    ax.set_xticks(positions, labels, rotation=32, ha="right")
    ax.set_ylim(0, max(0.05, max(fractions) * 1.22))
    ax.set_ylabel("Primary FLNC / Q7-pass reads")
    ax.set_title("Seven-tissue FLNC classification", fontweight="bold")
    ax.grid(axis="y", color=LIGHT_GRAY, linewidth=0.6)
    for position, row, fraction in zip(positions, rows, fractions, strict=True):
        ax.text(position, fraction + max(fractions) * 0.035, f"{100 * fraction:.1f}%", ha="center", fontsize=7.5)
        plotted.append(
            {
                "panel": "a",
                "metric": "primary_flnc_fraction_of_q7_pass",
                "group": row["accession"],
                "estimator": row["tissue"],
                "value": fraction,
                "numerator": row["pychopper_primary_flnc"],
                "denominator": row["prefilter_pass_reads"],
                "ci_lower": "",
                "ci_upper": "",
            }
        )


def draw_review_yield(
    ax,
    rows: list[dict[str, str]],
    primary_rows: list[dict[str, str]],
    summary: dict[str, Any],
    plotted: list[dict[str, Any]],
) -> None:
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[row["estimator"]].append(row)
    require(set(ESTIMATOR_COLORS) <= set(grouped), "review-yield estimator grid is incomplete")
    maximum_budget = 0
    maximum_supported = 0
    for estimator in ESTIMATOR_COLORS:
        ordered = sorted(grouped[estimator], key=lambda row: int(row["review_budget"]))
        budgets = [int(row["review_budget"]) for row in ordered]
        supported = [int(row["full_chain_supported"]) for row in ordered]
        maximum_budget = max(maximum_budget, max(budgets))
        maximum_supported = max(maximum_supported, max(supported))
        ax.plot(
            budgets,
            supported,
            marker="o",
            linewidth=1.7,
            markersize=4.5,
            color=ESTIMATOR_COLORS[estimator],
            label=ESTIMATOR_LABELS[estimator],
        )
        for row, budget, count in zip(ordered, budgets, supported, strict=True):
            plotted.append(
                {
                    "panel": "b",
                    "metric": "full_chain_supported_review_yield",
                    "group": budget,
                    "estimator": estimator,
                    "value": count,
                    "numerator": count,
                    "denominator": budget,
                    "ci_lower": "",
                    "ci_upper": "",
                }
            )
    ax.plot([0, maximum_budget], [0, maximum_budget], color=LIGHT_GRAY, linestyle="--", linewidth=1, zorder=0)
    ax.set_xlim(0, maximum_budget * 1.04)
    ax.set_ylim(0, max(1, maximum_supported * 1.17))
    ax.set_xlabel("Candidates reviewed")
    ax.set_ylabel("FLNC full-chain-supported candidates")
    ax.set_title("Natural evidence at fixed review budgets", fontweight="bold")
    ax.grid(color=LIGHT_GRAY, linewidth=0.6)
    ax.legend(frameon=False, loc="upper left")
    require(len(primary_rows) == 1, "primary rank delta must contain one row")
    primary = primary_rows[0]
    delta = float(primary["observed_delta_supported"])
    lower = float(primary["ci_lower"])
    upper = float(primary["ci_upper"])
    audit_text = (
        f"Primary v0.4 − baseline at top {int(primary['review_budget']):,}: "
        f"{delta:+.1f} supported\n95% CI {lower:+.1f} to {upper:+.1f}; "
        f"strict cases = {int(summary['counts']['strict_cases'])}"
    )
    ax.text(
        0.98,
        0.045,
        audit_text,
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=7.5,
        color=PURPLE,
        fontweight="bold",
        bbox={"boxstyle": "round,pad=0.35", "facecolor": "white", "edgecolor": PURPLE},
    )
    plotted.append(
        {
            "panel": "b",
            "metric": "primary_supported_delta",
            "group": primary["review_budget"],
            "estimator": f"{primary['primary_estimator']}_minus_{primary['comparator_estimator']}",
            "value": delta,
            "numerator": "",
            "denominator": "",
            "ci_lower": lower,
            "ci_upper": upper,
        }
    )


def draw_enrichment(ax, rows: list[dict[str, str]], plotted: list[dict[str, Any]]) -> None:
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[row["estimator"]].append(row)
    require(set(ESTIMATOR_COLORS) <= set(grouped), "enrichment estimator grid is incomplete")
    for estimator in ESTIMATOR_COLORS:
        ordered = sorted(grouped[estimator], key=lambda row: int(row["review_budget"]))
        budgets = [int(row["review_budget"]) for row in ordered]
        folds = [float(row["fold_enrichment_over_random_mean"]) for row in ordered]
        ax.plot(
            budgets,
            folds,
            marker="o",
            linewidth=1.7,
            markersize=4.5,
            color=ESTIMATOR_COLORS[estimator],
            label=ESTIMATOR_LABELS[estimator],
        )
        for row, budget, fold in zip(ordered, budgets, folds, strict=True):
            plotted.append(
                {
                    "panel": "c",
                    "metric": "fold_enrichment_over_random_mean",
                    "group": budget,
                    "estimator": estimator,
                    "value": fold,
                    "numerator": row["observed_full_chain_supported"],
                    "denominator": row["random_mean"],
                    "ci_lower": row["random_ci_lower"],
                    "ci_upper": row["random_ci_upper"],
                }
            )
    ax.axhline(1.0, color=ORANGE, linestyle="--", linewidth=1.2, label="Random expectation")
    ax.set_xlabel("Candidates reviewed")
    ax.set_ylabel("Fold enrichment over random mean")
    ax.set_title("Chromosome-matched random-rank enrichment", fontweight="bold")
    ax.grid(color=LIGHT_GRAY, linewidth=0.6)
    ax.legend(frameon=False)


def draw_transition_audit(
    ax,
    transition_rows: list[dict[str, str]],
    summary: dict[str, Any],
    plotted: list[dict[str, Any]],
) -> None:
    ordered = sorted(transition_rows, key=lambda row: int(row["candidates"]), reverse=True)
    displayed = ordered[:6]
    labels = [
        f"{row['raw_ts_state'].replace('_', ' ')}\n→ {row['pychopper_flnc_state'].replace('_', ' ')}"
        for row in displayed
    ]
    counts = [int(row["candidates"]) for row in displayed]
    positions = list(range(len(displayed)))
    ax.barh(
        positions,
        [max(1, count - 1) for count in counts],
        left=1,
        color=[PURPLE if index == 0 else "#B8A7CE" for index in positions],
    )
    ax.set_yticks(positions, labels, fontsize=7.2)
    ax.invert_yaxis()
    ax.set_xscale("log")
    ax.set_xlabel("Candidates (log scale)")
    ax.set_title("Raw-read to FLNC evidence-state audit", fontweight="bold")
    ax.grid(axis="x", color=LIGHT_GRAY, linewidth=0.6)
    for position, row, count in zip(positions, displayed, counts, strict=True):
        ax.text(count * 1.05, position, f"{count:,}", va="center", fontsize=7.5)
        plotted.append(
            {
                "panel": "d",
                "metric": "raw_ts_to_flnc_transition",
                "group": f"{row['raw_ts_state']}->{row['pychopper_flnc_state']}",
                "estimator": "sensitivity_audit",
                "value": count,
                "numerator": count,
                "denominator": summary["counts"]["candidate_universe"],
                "ci_lower": "",
                "ci_upper": "",
            }
        )


def write_source_data(path: Path, rows: list[dict[str, Any]]) -> None:
    fields = (
        "panel",
        "metric",
        "group",
        "estimator",
        "value",
        "numerator",
        "denominator",
        "ci_lower",
        "ci_upper",
    )
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", required=True, type=Path)
    parser.add_argument("--expected-sha256s-sha256", required=True)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    source = args.source_dir.resolve(strict=True)
    output = args.output_dir.resolve()
    working = Path(f"{output}.working")
    require(not output.exists() and not working.exists(), f"refusing to overwrite: {output}")
    summary = verify_source(source, args.expected_sha256s_sha256)
    classification = read_tsv(source / "table_pychopper_classification.tsv")
    review = read_tsv(source / "figure_review_yield.tsv")
    enrichment = read_tsv(source / "figure_random_enrichment.tsv")
    transitions = read_tsv(source / "figure_raw_ts_to_flnc_transitions.tsv")
    primary = read_tsv(source / "figure_primary_rank_delta.tsv")

    working.mkdir(parents=True)
    try:
        plotted: list[dict[str, Any]] = []
        fig, axes = plt.subplots(2, 2, figsize=(12.0, 8.7), layout="constrained")
        draw_classification(axes[0, 0], classification, plotted)
        draw_review_yield(axes[0, 1], review, primary, summary, plotted)
        draw_enrichment(axes[1, 0], enrichment, plotted)
        draw_transition_audit(axes[1, 1], transitions, summary, plotted)
        for label, ax in zip(("a", "b", "c", "d"), axes.ravel(), strict=True):
            panel_label(ax, label)
        fig.suptitle(
            "Golden Delicious long-read validation of frozen review prioritization",
            fontsize=14,
            fontweight="bold",
        )

        os.environ["SOURCE_DATE_EPOCH"] = "0"
        fixed_date = datetime.datetime(1970, 1, 1, tzinfo=datetime.timezone.utc)
        outputs = {
            "png": working / "figure5_apple_flnc_validation.png",
            "pdf": working / "figure5_apple_flnc_validation.pdf",
            "svg": working / "figure5_apple_flnc_validation.svg",
        }
        fig.savefig(outputs["png"], dpi=300, facecolor="white", metadata={"Software": "PloidyPatch"})
        fig.savefig(
            outputs["pdf"],
            facecolor="white",
            metadata={"Creator": "PloidyPatch", "CreationDate": fixed_date, "ModDate": fixed_date},
        )
        fig.savefig(outputs["svg"], facecolor="white", metadata={"Creator": "PloidyPatch", "Date": "1970-01-01T00:00:00Z"})
        plt.close(fig)

        source_data = working / "figure5_source_data.tsv"
        write_source_data(source_data, plotted)
        manifest = {
            "schema_version": "ploidypatch.apple_flnc_figure.v1",
            "interpretation": summary["interpretation"],
            "automatic_annotation_patch": False,
            "source": {
                "path": os.fspath(source),
                "sha256s_sha256": args.expected_sha256s_sha256,
            },
            "generator_sha256": sha256_file(Path(__file__).resolve()),
            "matplotlib_version": matplotlib.__version__,
            "source_data": {
                "file_name": source_data.name,
                "rows": len(plotted),
                "sha256": sha256_file(source_data),
            },
            "outputs": {
                key: {"file_name": path.name, "bytes": path.stat().st_size, "sha256": sha256_file(path)}
                for key, path in outputs.items()
            },
        }
        (working / "figure_manifest.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        with (working / "SHA256SUMS").open("x", encoding="utf-8", newline="") as handle:
            for path in sorted(item for item in working.iterdir() if item.is_file() and item.name != "SHA256SUMS"):
                handle.write(f"{sha256_file(path)}  ./{path.name}\n")
        working.replace(output)
        for path in output.iterdir():
            path.chmod(stat.S_IMODE(path.stat().st_mode) & ~0o222)
        output.chmod(stat.S_IMODE(output.stat().st_mode) & ~0o222)
    except Exception:
        if working.exists():
            shutil.rmtree(working)
        raise
    print(json.dumps({"figure_root": os.fspath(output), "source_rows": len(plotted)}, sort_keys=True))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except FigureError as error:
        print(f"Apple FLNC figure failed: {error}", file=sys.stderr)
        raise SystemExit(1) from error
