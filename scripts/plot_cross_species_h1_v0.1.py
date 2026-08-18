#!/usr/bin/env python3
"""Render the checksum-bound cross-species H1 manuscript figure."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from ploidypatch.artifact_manifest import write_sha256sums, verify_sha256sums


BLUE = "#0072B2"
ORANGE = "#D55E00"
GREEN = "#009E73"
GRAY = "#666666"
LIGHT_GRAY = "#DDDDDD"

matplotlib.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "svg.hashsalt": "ploidypatch-cross-species-h1-v0.1",
    }
)


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


def load_json(path: Path, expected_sha256: str) -> dict[str, Any]:
    require(path.is_file() and path.stat().st_size > 0, f"input is missing: {path}")
    require(sha256_file(path) == expected_sha256, f"input SHA-256 mismatch: {path}")
    value = json.loads(path.read_text(encoding="utf-8"))
    require(isinstance(value, dict), f"input root is not an object: {path}")
    return value


def extract(
    populus: dict[str, Any], actinidia: dict[str, Any], coffea: dict[str, Any]
) -> list[dict[str, Any]]:
    pop = populus["H1_chain_preserving_candidate_ceiling"]
    act = actinidia["H1_chain_preserving_candidate_ceiling"]
    cof_primary = coffea["arms"]["combined"]
    cof_delta = coffea["paired_event_bootstrap"]
    rows = [
        {
            "species": "Populus",
            "events": pop["events"],
            "retain": pop["primary_recovered"],
            "suppress": pop["legacy_recovered"],
            "delta": pop["observed_delta"],
            "ci_lower": pop["ci_lower"],
            "ci_upper": pop["ci_upper"],
            "formal_status": "confirmatory pass",
        },
        {
            "species": "Actinidia",
            "events": act["events"],
            "retain": act["primary_recovered"],
            "suppress": act["legacy_recovered"],
            "delta": act["observed_delta"],
            "ci_lower": act["ci_lower"],
            "ci_upper": act["ci_upper"],
            "formal_status": "H1 pass; composite negative",
        },
        {
            "species": "Coffea",
            "events": cof_primary["retain_distinct"]["events"],
            "retain": cof_primary["retain_distinct"][
                "exact_phased_cds_chain_recovered"
            ],
            "suppress": cof_primary["suppress_overlap"][
                "exact_phased_cds_chain_recovered"
            ],
            "delta": cof_delta["observed_delta"],
            "ci_lower": cof_delta["ci_lower"],
            "ci_upper": cof_delta["ci_upper"],
            "formal_status": "formal H1 positive",
        },
    ]
    require(all(row["events"] == 800 for row in rows), "H1 event universes differ")
    require(
        all(row["delta"] > 0 and row["ci_lower"] > 0 for row in rows),
        "a fixed H1 effect is not positive",
    )
    for row in rows:
        observed = (row["retain"] - row["suppress"]) / row["events"]
        require(abs(observed - row["delta"]) < 1e-12, "H1 count/delta mismatch")
    return rows


def draw(rows: list[dict[str, Any]], output: Path) -> list[dict[str, Any]]:
    fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.5), layout="constrained")
    ax = axes[0]
    x = list(range(len(rows)))
    suppress = [row["suppress"] / row["events"] for row in rows]
    retain = [row["retain"] / row["events"] for row in rows]
    for index, row in enumerate(rows):
        ax.plot(
            [index - 0.12, index + 0.12],
            [suppress[index], retain[index]],
            color=LIGHT_GRAY,
            linewidth=2.2,
            zorder=1,
        )
        ax.scatter(index - 0.12, suppress[index], color=GRAY, s=58, zorder=3)
        ax.scatter(index + 0.12, retain[index], color=BLUE, s=58, zorder=3)
        ax.text(
            index - 0.12,
            suppress[index] - 0.045,
            f"{row['suppress']}/800",
            ha="center",
            va="top",
            fontsize=8,
            color=GRAY,
        )
        ax.text(
            index + 0.12,
            retain[index] + 0.035,
            f"{row['retain']}/800",
            ha="center",
            va="bottom",
            fontsize=8,
            color=BLUE,
        )
    ax.set_xticks(x, [row["species"] for row in rows])
    ax.set_ylim(0.30, 1.00)
    ax.set_ylabel("Exact phased-CDS event recall")
    ax.set_title("A  Chain-preserving recovery", loc="left", fontweight="bold")
    ax.grid(axis="y", color=LIGHT_GRAY, linewidth=0.7)
    ax.scatter([], [], color=GRAY, s=50, label="Suppress overlap")
    ax.scatter([], [], color=BLUE, s=50, label="Retain distinct chains")
    ax.legend(frameon=False, loc="lower left")

    ax = axes[1]
    y = list(reversed(range(len(rows))))
    for position, row in zip(y, rows, strict=True):
        lower = row["delta"] - row["ci_lower"]
        upper = row["ci_upper"] - row["delta"]
        color = GREEN if row["species"] != "Actinidia" else ORANGE
        ax.errorbar(
            row["delta"],
            position,
            xerr=[[lower], [upper]],
            fmt="o",
            color=color,
            capsize=4,
            markersize=6,
            linewidth=1.8,
        )
        ax.text(
            0.117,
            position,
            f"{row['delta']:+.5f}  [{row['ci_lower']:+.5f}, {row['ci_upper']:+.5f}]",
            va="center",
            fontsize=8,
        )
        ax.text(
            0.002,
            position - 0.27,
            row["formal_status"],
            va="top",
            fontsize=7.5,
            color=GRAY,
        )
    ax.axvline(0, color=GRAY, linewidth=1, linestyle="--")
    ax.set_xlim(-0.005, 0.205)
    ax.set_ylim(-0.65, 2.55)
    ax.set_yticks(y, [row["species"] for row in rows])
    ax.set_xlabel("Recall difference: retain distinct - suppress overlap")
    ax.set_title("B  Paired-event effect (95% CI)", loc="left", fontweight="bold")
    ax.grid(axis="x", color=LIGHT_GRAY, linewidth=0.7)

    fig.suptitle(
        "Chain preservation reproduces across three prospective plant systems",
        fontsize=12,
        fontweight="bold",
    )
    for suffix in ("png", "pdf", "svg"):
        fig.savefig(output / f"figure_cross_species_h1_v0.1.{suffix}", dpi=300)
    plt.close(fig)

    source_rows: list[dict[str, Any]] = []
    for row in rows:
        for policy, recovered in (
            ("suppress_overlap", row["suppress"]),
            ("retain_distinct", row["retain"]),
        ):
            source_rows.append(
                {
                    "species": row["species"],
                    "metric": "exact_phased_cds_event_recall",
                    "policy": policy,
                    "recovered": recovered,
                    "events": row["events"],
                    "value": recovered / row["events"],
                    "ci_lower": "",
                    "ci_upper": "",
                    "formal_status": row["formal_status"],
                }
            )
        source_rows.append(
            {
                "species": row["species"],
                "metric": "paired_recall_delta",
                "policy": "retain_distinct_minus_suppress_overlap",
                "recovered": "",
                "events": row["events"],
                "value": row["delta"],
                "ci_lower": row["ci_lower"],
                "ci_upper": row["ci_upper"],
                "formal_status": row["formal_status"],
            }
        )
    return source_rows


def write_source_data(path: Path, rows: list[dict[str, Any]]) -> None:
    fields = (
        "species",
        "metric",
        "policy",
        "recovered",
        "events",
        "value",
        "ci_lower",
        "ci_upper",
        "formal_status",
    )
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--populus", required=True, type=Path)
    parser.add_argument("--populus-sha256", required=True)
    parser.add_argument("--actinidia", required=True, type=Path)
    parser.add_argument("--actinidia-sha256", required=True)
    parser.add_argument("--coffea", required=True, type=Path)
    parser.add_argument("--coffea-sha256", required=True)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args(argv)
    inputs = {
        "populus": (args.populus.resolve(strict=True), args.populus_sha256),
        "actinidia": (args.actinidia.resolve(strict=True), args.actinidia_sha256),
        "coffea": (args.coffea.resolve(strict=True), args.coffea_sha256),
    }
    values = {
        name: load_json(path, digest) for name, (path, digest) in inputs.items()
    }
    rows = extract(values["populus"], values["actinidia"], values["coffea"])
    output = args.output_dir.resolve()
    working = Path(f"{output}.working")
    require(not output.exists() and not working.exists(), f"refusing to overwrite: {output}")
    working.mkdir(parents=True)
    try:
        source_rows = draw(rows, working)
        source_path = working / "figure_cross_species_h1_v0.1_source_data.tsv"
        write_source_data(source_path, source_rows)
        outputs = {
            suffix: {
                "relative_path": f"figure_cross_species_h1_v0.1.{suffix}",
                "sha256": sha256_file(
                    working / f"figure_cross_species_h1_v0.1.{suffix}"
                ),
            }
            for suffix in ("png", "pdf", "svg")
        }
        manifest = {
            "schema_version": "ploidypatch.cross_species_h1_figure.v1",
            "inputs": {
                name: {
                    "project_relative_path": path.relative_to(Path.cwd()).as_posix(),
                    "sha256": digest,
                }
                for name, (path, digest) in inputs.items()
            },
            "species": [row["species"] for row in rows],
            "source_data": {
                "relative_path": source_path.name,
                "rows": len(source_rows),
                "sha256": sha256_file(source_path),
            },
            "outputs": outputs,
            "formal_interpretation": {
                "all_H1_point_estimates_positive": True,
                "all_H1_CI_lower_bounds_positive": True,
                "actinidia_composite_result_overridden": False,
                "coffea_ranker_claim": False,
            },
        }
        (working / "figure_manifest.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        write_sha256sums(working)
        verify_sha256sums(working, ignore_checksum_file=True)
        os.replace(working, output)
    except BaseException:
        shutil.rmtree(working, ignore_errors=True)
        raise
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
