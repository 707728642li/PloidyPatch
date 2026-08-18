#!/usr/bin/env python3
"""Render the PloidyPatch system and evidence-firewall manuscript figure."""

from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import os
import shutil
import stat
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle  # noqa: E402


matplotlib.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "svg.hashsalt": "ploidypatch-system-architecture-v0.1",
    }
)

BLUE = "#0072B2"
LIGHT_BLUE = "#DCEEF8"
ORANGE = "#D55E00"
LIGHT_ORANGE = "#FBE4D5"
GREEN = "#009E73"
LIGHT_GREEN = "#D9F2E9"
PURPLE = "#7A5195"
LIGHT_PURPLE = "#E9DFF0"
GRAY = "#5A5A5A"
LIGHT_GRAY = "#F2F2F2"
RED = "#B2182B"


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


def box(ax, x, y, width, height, text, *, face, edge, fontsize=8.5, linewidth=1.2):
    patch = FancyBboxPatch(
        (x, y),
        width,
        height,
        boxstyle="round,pad=0.012,rounding_size=0.018",
        facecolor=face,
        edgecolor=edge,
        linewidth=linewidth,
    )
    ax.add_patch(patch)
    ax.text(x + width / 2, y + height / 2, text, ha="center", va="center", fontsize=fontsize)
    return patch


def arrow(ax, start, end, *, color=GRAY, style="-|>", linewidth=1.3, linestyle="-"):
    patch = FancyArrowPatch(
        start,
        end,
        arrowstyle=style,
        mutation_scale=11,
        color=color,
        linewidth=linewidth,
        linestyle=linestyle,
        shrinkA=2,
        shrinkB=2,
    )
    ax.add_patch(patch)
    return patch


def panel_label(ax, label: str) -> None:
    ax.text(-0.02, 1.04, label, transform=ax.transAxes, fontsize=14, fontweight="bold", va="top")


def draw_workflow(ax) -> None:
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    panel_label(ax, "a")
    ax.text(0.02, 0.96, "Duplication-aware candidate-to-review workflow", fontsize=11, fontweight="bold", va="top")

    y = 0.54
    width = 0.125
    height = 0.25
    xs = [0.01, 0.17, 0.33, 0.49, 0.65, 0.81]
    labels = [
        "Target genome\n+ annotation\n\nCandidate-only\nrelatives",
        "Mature projection\nmethods\n\nminiprot\nGeMoMa\nLiftOn",
        "Exact phased-CDS\nchain consensus\n\nRetain distinct\nalternatives",
        "Optional plant evidence\n(applicability-gated)\n\nWGD partner\nchain topology",
        "Review ranking\n+ conflict guard\n\nBaseline preserved\nwhen abstaining",
        "Human review\n\nImmutable,\nreversible patch",
    ]
    faces = [LIGHT_BLUE, LIGHT_GRAY, LIGHT_BLUE, LIGHT_GREEN, LIGHT_PURPLE, LIGHT_ORANGE]
    edges = [BLUE, GRAY, BLUE, GREEN, PURPLE, ORANGE]
    for x, label, face, edge in zip(xs, labels, faces, edges, strict=True):
        box(ax, x, y, width, height, label, face=face, edge=edge, fontsize=7.8)
    for left, right in zip(xs[:-1], xs[1:], strict=True):
        arrow(ax, (left + width, y + height / 2), (right, y + height / 2))

    ax.text(0.337, 0.47, "conflict sets remain explicit", color=BLUE, fontsize=7.5, ha="center")
    ax.text(0.712, 0.47, "rank, not probability", color=PURPLE, fontsize=7.5, ha="center")
    ax.text(0.873, 0.47, "no automatic approval", color=ORANGE, fontsize=7.5, ha="center")

    box(
        ax,
        0.17,
        0.12,
        0.445,
        0.18,
        "Complete-annotation control (evaluator-owned)\nsubtracts reproducible background proposals",
        face="#FFF7E6",
        edge=ORANGE,
        fontsize=8.2,
    )
    arrow(ax, (0.235, y), (0.31, 0.30), color=ORANGE, linestyle="--")
    arrow(ax, (0.615, 0.21), (0.49, y), color=ORANGE, linestyle="--")


def draw_firewall(ax) -> None:
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    panel_label(ax, "b")
    ax.text(0.02, 0.96, "Confirmatory evidence firewall", fontsize=11, fontweight="bold", va="top")

    ax.add_patch(Rectangle((0.02, 0.10), 0.42, 0.76, facecolor="#F4FAFD", edgecolor=BLUE, linewidth=1.3))
    ax.add_patch(Rectangle((0.56, 0.10), 0.42, 0.76, facecolor="#FFF6F4", edgecolor=RED, linewidth=1.3))
    ax.plot([0.50, 0.50], [0.08, 0.88], color=GRAY, linewidth=1.5, linestyle=(0, (4, 3)))
    ax.text(0.50, 0.91, "truth boundary", ha="center", color=GRAY, fontsize=8, fontweight="bold")

    ax.text(0.23, 0.81, "BLIND CANDIDATE SIDE", color=BLUE, ha="center", fontweight="bold")
    ax.text(0.77, 0.81, "EVALUATOR SIDE", color=RED, ha="center", fontweight="bold")
    left = [
        "Perturbed target genome/GFF3",
        "Candidate-only relatives",
        "Frozen source, model and environments",
        "No network; no /nas_data mount",
        "Raw predictions, pools, features, scores",
    ]
    right = [
        "Complete target annotation",
        "Independent evaluator-only relatives",
        "WGD pairs and hidden event truth",
        "Complete-control adaptation",
        "Labels, sentinels and statistics",
    ]
    for index, text in enumerate(left):
        box(ax, 0.06, 0.69 - index * 0.11, 0.34, 0.075, text, face="white", edge=BLUE, fontsize=7.6)
    for index, text in enumerate(right):
        box(ax, 0.60, 0.69 - index * 0.11, 0.34, 0.075, text, face="white", edge=RED, fontsize=7.6)

    arrow(ax, (0.40, 0.145), (0.60, 0.145), color=GREEN, linewidth=2.0)
    ax.text(0.50, 0.17, "checksum freeze before reveal", color=GREEN, fontsize=7.5, ha="center", fontweight="bold")
    ax.text(0.50, 0.112, "no evaluator information flows back", color=RED, fontsize=7.3, ha="center")


def draw_safety(ax) -> None:
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    panel_label(ax, "c")
    ax.text(0.02, 0.96, "Review-only safety contract", fontsize=11, fontweight="bold", va="top")

    items = [
        ("1", "Candidate ceiling", "Keep distinct phased-CDS chains;\ndo not erase alternatives by overlap alone", BLUE, LIGHT_BLUE),
        (
            "2",
            "Applicability-gated\nplant evidence",
            "Topology augments review ranking only where valid;\nabstain and retain the baseline otherwise",
            GREEN,
            LIGHT_GREEN,
        ),
        ("3", "Conflict guard", "Baseline and v0.4 winner mappings must\nmatch exactly in every conflict set", PURPLE, LIGHT_PURPLE),
        ("4", "Human decision", "Reviewer, UTC time and rationale required;\nunreviewed candidates are excluded", ORANGE, LIGHT_ORANGE),
        ("5", "Reversible output", "Append complete hierarchies only;\napply and revert require exact bytes", GRAY, LIGHT_GRAY),
    ]
    y = 0.78
    for number, title, body, edge, face in items:
        circle = FancyBboxPatch(
            (0.04, y), 0.08, 0.085, boxstyle="round,pad=0.005,rounding_size=0.04",
            facecolor=edge, edgecolor=edge,
        )
        ax.add_patch(circle)
        ax.text(0.08, y + 0.043, number, color="white", fontsize=10, fontweight="bold", ha="center", va="center")
        box(ax, 0.15, y - 0.003, 0.80, 0.095, "", face=face, edge=edge, linewidth=1.0)
        title_is_multiline = "\n" in title
        ax.text(
            0.18,
            y + (0.044 if title_is_multiline else 0.058),
            title,
            color=edge,
            fontweight="bold",
            fontsize=7.8 if title_is_multiline else 8.3,
            linespacing=0.9,
            va="center",
        )
        ax.text(0.56, y + 0.044, body, color="#222222", fontsize=7.2, va="center")
        y -= 0.145
    ax.text(
        0.50,
        0.055,
        "Output: an auditable review queue — never a silent annotation rewrite",
        ha="center",
        color=RED,
        fontsize=8.2,
        fontweight="bold",
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    output = args.output_dir.resolve()
    working = Path(f"{output}.working")
    require(not output.exists() and not working.exists(), f"refusing to overwrite: {output}")
    working.mkdir(parents=True)
    try:
        fig = plt.figure(figsize=(14.0, 8.0), layout="constrained")
        grid = fig.add_gridspec(2, 2, height_ratios=(1.0, 1.18), width_ratios=(1.25, 0.75))
        draw_workflow(fig.add_subplot(grid[0, :]))
        draw_firewall(fig.add_subplot(grid[1, 0]))
        draw_safety(fig.add_subplot(grid[1, 1]))
        fig.suptitle("PloidyPatch: duplication-aware candidate prioritization with explicit evidence custody", fontsize=15, fontweight="bold")

        os.environ["SOURCE_DATE_EPOCH"] = "0"
        fixed_date = datetime.datetime(1970, 1, 1, tzinfo=datetime.timezone.utc)
        outputs = {
            "png": working / "figure1_system_architecture.png",
            "pdf": working / "figure1_system_architecture.pdf",
            "svg": working / "figure1_system_architecture.svg",
        }
        fig.savefig(outputs["png"], dpi=300, facecolor="white", metadata={"Software": "PloidyPatch"})
        fig.savefig(
            outputs["pdf"],
            facecolor="white",
            metadata={"Creator": "PloidyPatch", "CreationDate": fixed_date, "ModDate": fixed_date},
        )
        fig.savefig(
            outputs["svg"],
            facecolor="white",
            metadata={"Creator": "PloidyPatch", "Date": "1970-01-01T00:00:00Z"},
        )
        plt.close(fig)
        manifest = {
            "schema_version": "ploidypatch.system_architecture_figure.v1",
            "claim_boundary": {
                "automatic_approval": False,
                "calibrated_probability": False,
                "review_required": True,
            },
            "panels": {
                "a": "candidate-to-review workflow and complete-control subtraction",
                "b": "blind candidate and evaluator-only data custody",
                "c": "review-only safety and reversible patch contract",
            },
            "matplotlib_version": matplotlib.__version__,
            "generator_sha256": sha256_file(Path(__file__).resolve()),
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
        working.rename(output)
        for path in output.iterdir():
            path.chmod(stat.S_IMODE(path.stat().st_mode) & ~0o222)
        output.chmod(stat.S_IMODE(output.stat().st_mode) & ~0o222)
    except Exception:
        if working.exists():
            shutil.rmtree(working)
        raise
    print(json.dumps({"figure_root": os.fspath(output), "manifest_sha256": sha256_file(output / "figure_manifest.json")}, sort_keys=True))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except FigureError as error:
        print(f"system figure failed: {error}", file=sys.stderr)
        raise SystemExit(1) from error
