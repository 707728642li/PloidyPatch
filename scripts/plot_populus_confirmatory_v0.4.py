#!/usr/bin/env python3
"""Render the checksum-bound Populus confirmatory manuscript figure."""

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
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


matplotlib.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "svg.hashsalt": "ploidypatch-populus-confirmatory-v0.4",
    }
)

BLUE = "#0072B2"
ORANGE = "#D55E00"
GREEN = "#009E73"
PURPLE = "#7A5195"
GRAY = "#666666"
LIGHT_GRAY = "#DDDDDD"


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


def load_evaluation(path: Path, expected_sha256: str) -> dict[str, Any]:
    require(re.fullmatch(r"[0-9a-f]{64}", expected_sha256) is not None, "expected SHA-256 is invalid")
    require(path.is_file() and path.stat().st_size > 0, f"evaluation is missing: {path}")
    require(sha256_file(path) == expected_sha256, "evaluation SHA-256 mismatch")
    with path.open(encoding="utf-8") as handle:
        value = json.load(handle)
    require(isinstance(value, dict), "evaluation root is not an object")
    require(value.get("schema_version") == "ploidypatch.external_ranker_evaluation.v0.4", "evaluation schema differs")
    require(value.get("evaluation_role") == "untouched_confirmatory_external_species", "evaluation role differs")
    require(value.get("confirmatory_pass") is True, "confirmatory run did not pass")
    require(all(value.get("gates", {}).values()), "not every confirmatory gate passed")
    return value


def panel_label(ax, label: str) -> None:
    ax.text(-0.12, 1.07, label, transform=ax.transAxes, fontsize=14, fontweight="bold", va="top")


def draw_h1(ax, evaluation: dict[str, Any], rows: list[dict[str, Any]]) -> None:
    h1 = evaluation["H1_chain_preserving_candidate_ceiling"]
    events = h1["events"]
    labels = ["Legacy overlap\nsuppression", "Chain-preserving\npool"]
    recovered = [h1["legacy_recovered"], h1["primary_recovered"]]
    rates = [item / events for item in recovered]
    ax.bar([0, 1], rates, color=[GRAY, BLUE], width=0.62)
    ax.set_ylim(0, 1.03)
    ax.set_xticks([0, 1], labels)
    ax.set_ylabel("Exact phased-CDS event recall")
    ax.set_title("H1: candidate recovery ceiling", fontweight="bold")
    ax.grid(axis="y", color=LIGHT_GRAY, linewidth=0.6)
    for index, (count, rate) in enumerate(zip(recovered, rates, strict=True)):
        ax.text(index, rate + 0.025, f"{count}/{events}\n({rate:.3f})", ha="center", fontsize=8.5)
    ax.text(
        0.50,
        0.11,
        f"Δ = {h1['observed_delta']:+.5f}\n95% CI {h1['ci_lower']:+.5f} to {h1['ci_upper']:+.5f}",
        transform=ax.transAxes,
        ha="center",
        va="bottom",
        color=BLUE,
        fontweight="bold",
        bbox={"boxstyle": "round,pad=0.35", "facecolor": "white", "edgecolor": BLUE},
    )
    rows.extend(
        [
            {"panel": "a", "metric": "exact_event_recall", "estimator": "legacy", "value": rates[0], "numerator": recovered[0], "denominator": events, "ci_lower": "", "ci_upper": "", "gate": ""},
            {"panel": "a", "metric": "exact_event_recall", "estimator": "chain_preserving", "value": rates[1], "numerator": recovered[1], "denominator": events, "ci_lower": "", "ci_upper": "", "gate": ""},
            {"panel": "a", "metric": "paired_recall_delta", "estimator": "chain_preserving_minus_legacy", "value": h1["observed_delta"], "numerator": "", "denominator": "", "ci_lower": h1["ci_lower"], "ci_upper": h1["ci_upper"], "gate": "delta_and_lower_bound_positive"},
        ]
    )


def draw_h2(ax, evaluation: dict[str, Any], rows: list[dict[str, Any]]) -> None:
    metrics = evaluation["ranking_metrics"]
    h2 = evaluation["H2_v04_guard_minus_baseline"]
    baseline = metrics["baseline"]["average_precision"]
    guard = metrics["v04_guard"]["average_precision"]
    ax.plot([baseline, guard], [1, 0], color=LIGHT_GRAY, linewidth=2.0, zorder=1)
    ax.scatter([baseline], [1], color=GRAY, s=70, zorder=3, label="Frozen baseline")
    ax.scatter([guard], [0], color=PURPLE, s=70, zorder=3, label="v0.4 guard")
    ax.set_yticks([1, 0], ["Frozen baseline", "v0.4 guard"])
    lower = min(baseline, guard) - 0.002
    upper = max(baseline, guard) + 0.002
    ax.set_xlim(lower, upper)
    ax.set_ylim(-0.7, 1.7)
    ax.set_xlabel("Candidate average precision")
    ax.set_title("H2: untouched review ranking", fontweight="bold")
    ax.grid(axis="x", color=LIGHT_GRAY, linewidth=0.6)
    ax.text(baseline, 1.18, f"{baseline:.6f}", ha="center", color=GRAY)
    ax.text(guard, -0.28, f"{guard:.6f}", ha="center", color=PURPLE)
    ax.text(
        0.50,
        0.55,
        f"ΔAP = {h2['observed_delta']:+.6f}\n95% CI {h2['ci_lower']:+.6f} to {h2['ci_upper']:+.6f}",
        transform=ax.transAxes,
        ha="center",
        va="center",
        color=PURPLE,
        fontweight="bold",
        bbox={"boxstyle": "round,pad=0.35", "facecolor": "white", "edgecolor": PURPLE},
    )
    rows.extend(
        [
            {"panel": "b", "metric": "average_precision", "estimator": "baseline", "value": baseline, "numerator": "", "denominator": "", "ci_lower": "", "ci_upper": "", "gate": ""},
            {"panel": "b", "metric": "average_precision", "estimator": "v04_guard", "value": guard, "numerator": "", "denominator": "", "ci_lower": "", "ci_upper": "", "gate": ""},
            {"panel": "b", "metric": "chromosome_bootstrap_delta_AP", "estimator": "v04_guard_minus_baseline", "value": h2["observed_delta"], "numerator": h2["replicates_valid"], "denominator": h2["replicates_requested"], "ci_lower": h2["ci_lower"], "ci_upper": h2["ci_upper"], "gate": "delta_and_lower_bound_positive"},
        ]
    )


def draw_review(ax, evaluation: dict[str, Any], rows: list[dict[str, Any]]) -> None:
    budgets = evaluation["review_budgets"]
    estimators = (("baseline", "Frozen baseline", GRAY), ("v03_primary", "v0.3", GREEN), ("v04_guard", "v0.4 guard", PURPLE))
    for estimator, label, color in estimators:
        points = sorted(
            (
                (record["reviewed"], record["true_positive"], key, record["precision"])
                for key, record in budgets[estimator].items()
            ),
            key=lambda item: item[0],
        )
        ax.plot(
            [item[0] for item in points],
            [item[1] for item in points],
            marker="o",
            markersize=4.5,
            linewidth=1.6,
            color=color,
            label=label,
        )
        for reviewed, positives, key, precision in points:
            rows.append(
                {"panel": "c", "metric": "review_true_positives", "estimator": estimator, "value": positives, "numerator": positives, "denominator": reviewed, "ci_lower": "", "ci_upper": "", "gate": key, "precision": precision}
            )
    ax.plot([0, 520], [0, 520], color=LIGHT_GRAY, linewidth=1, linestyle="--", zorder=0)
    ax.set_xlim(0, 520)
    ax.set_ylim(0, 510)
    ax.set_xlabel("Candidates reviewed")
    ax.set_ylabel("Exact positives retrieved")
    ax.set_title("Fixed review-budget yield", fontweight="bold")
    ax.grid(color=LIGHT_GRAY, linewidth=0.6)
    ax.legend(frameon=False, loc="upper left")


def draw_safety(ax, evaluation: dict[str, Any], rows: list[dict[str, Any]]) -> None:
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.set_title("Predeclared safety gates", fontweight="bold")
    conflict_total = evaluation["conflict_sets"]["v04_guard"]["conflict_sets_total"]
    mismatch = evaluation["winner_audit"]["mismatch_count"]
    gain = evaluation["v03_AP_gain_retained_fraction"]
    coverage = evaluation["counts"]["topology_positive_coverage"]
    collateral = evaluation["collateral_gate_by_arm"]
    approvals = evaluation["counts"]["automatic_approved"]
    items = [
        ("Conflict winner identity", f"{conflict_total:,}/{conflict_total:,}; mismatches = {mismatch}", mismatch == 0),
        ("v0.3 AP gain retained", f"{100 * gain:.2f}%  (gate ≥90%)", gain >= 0.90),
        ("Positive topology coverage", f"{100 * coverage:.2f}%  (gate ≥70%)", coverage >= 0.70),
        ("Collateral transcript loss", f"0 in {len(collateral)}/{len(collateral)} arms", all(collateral.values())),
        ("Automatic approvals", str(approvals), approvals == 0),
    ]
    y = 0.83
    for label, value, passed in items:
        ax.text(0.06, y, "✓" if passed else "✕", color=GREEN if passed else ORANGE, fontsize=15, fontweight="bold", va="center")
        ax.text(0.14, y + 0.025, label, fontweight="bold", fontsize=8.5, va="center")
        ax.text(0.14, y - 0.035, value, color=GRAY, fontsize=8.2, va="center")
        rows.append({"panel": "d", "metric": label, "estimator": "v04_guard", "value": value, "numerator": "", "denominator": "", "ci_lower": "", "ci_upper": "", "gate": "pass" if passed else "fail"})
        y -= 0.16
    ax.text(
        0.50,
        0.055,
        "Review priority only — no calibrated probability or automatic patch",
        color=ORANGE,
        fontsize=8.2,
        fontweight="bold",
        ha="center",
    )


def write_source_data(path: Path, rows: list[dict[str, Any]]) -> None:
    fields = (
        "panel",
        "metric",
        "estimator",
        "value",
        "numerator",
        "denominator",
        "ci_lower",
        "ci_upper",
        "gate",
        "precision",
    )
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evaluation-json", required=True, type=Path)
    parser.add_argument("--expected-sha256", required=True)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    evaluation_path = args.evaluation_json.resolve(strict=True)
    evaluation = load_evaluation(evaluation_path, args.expected_sha256)
    output = args.output_dir.resolve()
    working = Path(f"{output}.working")
    require(not output.exists() and not working.exists(), f"refusing to overwrite: {output}")
    working.mkdir(parents=True)
    try:
        source_rows: list[dict[str, Any]] = []
        fig, axes = plt.subplots(2, 2, figsize=(11.5, 8.5), layout="constrained")
        draw_h1(axes[0, 0], evaluation, source_rows)
        draw_h2(axes[0, 1], evaluation, source_rows)
        draw_review(axes[1, 0], evaluation, source_rows)
        draw_safety(axes[1, 1], evaluation, source_rows)
        for label, ax in zip(("a", "b", "c", "d"), axes.ravel(), strict=True):
            panel_label(ax, label)
        fig.suptitle("Untouched Populus confirmation of PloidyPatch v0.4", fontsize=14, fontweight="bold")

        os.environ["SOURCE_DATE_EPOCH"] = "0"
        fixed_date = datetime.datetime(1970, 1, 1, tzinfo=datetime.timezone.utc)
        outputs = {
            "png": working / "figure4_populus_confirmatory.png",
            "pdf": working / "figure4_populus_confirmatory.pdf",
            "svg": working / "figure4_populus_confirmatory.svg",
        }
        fig.savefig(outputs["png"], dpi=300, facecolor="white", metadata={"Software": "PloidyPatch"})
        fig.savefig(outputs["pdf"], facecolor="white", metadata={"Creator": "PloidyPatch", "CreationDate": fixed_date, "ModDate": fixed_date})
        fig.savefig(outputs["svg"], facecolor="white", metadata={"Creator": "PloidyPatch", "Date": "1970-01-01T00:00:00Z"})
        plt.close(fig)
        source_data = working / "figure4_source_data.tsv"
        write_source_data(source_data, source_rows)
        manifest = {
            "schema_version": "ploidypatch.populus_confirmatory_figure.v1",
            "evaluation": {"path": os.fspath(evaluation_path), "sha256": args.expected_sha256},
            "confirmatory_pass": evaluation["confirmatory_pass"],
            "matplotlib_version": matplotlib.__version__,
            "generator_sha256": sha256_file(Path(__file__).resolve()),
            "source_data": {"file_name": source_data.name, "rows": len(source_rows), "sha256": sha256_file(source_data)},
            "outputs": {
                key: {"file_name": path.name, "bytes": path.stat().st_size, "sha256": sha256_file(path)}
                for key, path in outputs.items()
            },
        }
        (working / "figure_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        with (working / "SHA256SUMS").open("x", encoding="utf-8", newline="") as handle:
            for path in sorted(item for item in working.iterdir() if item.is_file() and item.name != "SHA256SUMS"):
                handle.write(f"{sha256_file(path)}  ./{path.name}\n")
        working.rename(output)
        for path in output.iterdir():
            path.chmod(stat.S_IMODE(path.stat().st_mode) & ~0o222)
        output.chmod(stat.S_IMODE(output.stat().st_mode) & ~0o222)
    except Exception:
        if working.exists():
            shutil.rmtree(working)
        raise
    print(json.dumps({"figure_root": os.fspath(output), "source_rows": len(source_rows)}, sort_keys=True))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except FigureError as error:
        print(f"Populus figure failed: {error}", file=sys.stderr)
        raise SystemExit(1) from error
