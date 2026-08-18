#!/usr/bin/env python3
"""Render the checksum-bound Actinidia external-validation figure."""

from __future__ import annotations

import argparse
import csv
import datetime
import hashlib
import json
import math
import os
import re
import shutil
import sys
from pathlib import Path
from typing import Any, Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


matplotlib.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "svg.hashsalt": "ploidypatch-actinidia-external-v0.5",
    }
)

BLUE = "#0072B2"
ORANGE = "#D55E00"
GREEN = "#009E73"
PURPLE = "#7A5195"
RED = "#B2182B"
GRAY = "#666666"
LIGHT_GRAY = "#D9D9D9"
PALE_RED = "#FBE9E7"
PALE_GREEN = "#E8F4EF"

SCHEMA_VERSION = "ploidypatch.actinidia_external_figure.v1"
EVALUATION_SCHEMA = "ploidypatch.external_ranker_evaluation.v0.5"
EVALUATION_ROLE = "target_level_predeclared_untouched_secondary_replication"
FORMAL_OUTCOME = "formal_negative_external_result"
TOPOLOGY_THRESHOLD = 0.70
EVIDENCE_SHA256 = {
    "bootstrap_deltas.tsv": "9dbb3a9d1f57633c7eda3fc3d2dfbfdb66c9e76f2301fa9a81cac1661188e953",
    "candidates.tsv": "09621368fd693dbf4dce74fa74e800eb2f3aedcc7cc4caa1a67bc0d64ac2a247",
    "evaluation.json": "f3c187bd5e86cf172343b7ac2b4ad10f331b0d4ab8a2273c5540c333109a0c3e",
}
GATE_ORDER = (
    "formal_evaluable",
    "H1_chain_preserving_ceiling",
    "H2_guard_AP_numerical",
    "H2_tested_in_fixed_sequence",
    "conflict_winner_mapping_identical_to_baseline",
    "top_1pct_review_noninferiority",
    "topology_positive_coverage",
    "retain_at_least_90pct_positive_v03_AP_gain",
    "zero_collateral_loss_all_arms",
    "automatic_approval_absent",
    "bootstrap_valid_replicates",
)
GATE_LABELS = {
    "formal_evaluable": "Formally evaluable",
    "H1_chain_preserving_ceiling": "H1 candidate ceiling",
    "H2_guard_AP_numerical": "H2 AP improvement",
    "H2_tested_in_fixed_sequence": "H2 fixed sequence",
    "conflict_winner_mapping_identical_to_baseline": "Conflict-winner identity",
    "top_1pct_review_noninferiority": "Top-1% review safety",
    "topology_positive_coverage": "Positive topology coverage",
    "retain_at_least_90pct_positive_v03_AP_gain": ">=90% v0.3 AP gain retained",
    "zero_collateral_loss_all_arms": "Zero collateral loss",
    "automatic_approval_absent": "No automatic approval",
    "bootstrap_valid_replicates": "Bootstrap validity",
}
SOURCE_FIELDS = (
    "panel",
    "metric",
    "estimator",
    "value",
    "numerator",
    "denominator",
    "ci_lower",
    "ci_upper",
    "threshold",
    "status",
    "note",
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


def parse_sha256sums(path: Path) -> dict[str, str]:
    require(path.is_file(), f"checksum manifest is missing: {path}")
    records: dict[str, str] = {}
    for line_number, raw in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        match = re.fullmatch(r"([0-9a-f]{64})  ([^/][^\r\n]*)", raw)
        require(match is not None, f"invalid SHA256SUMS line {line_number}")
        digest, relative = match.groups()
        require(relative not in records, f"duplicate checksum path: {relative}")
        records[relative] = digest
    return records


def verify_evidence(root: Path) -> dict[str, dict[str, Any]]:
    require(root.is_dir(), f"evidence directory is missing: {root}")
    records = parse_sha256sums(root / "SHA256SUMS")
    require(records == EVIDENCE_SHA256, "formal evidence SHA256SUMS differs")
    actual = {path.name for path in root.iterdir() if path.is_file()}
    require(actual == set(EVIDENCE_SHA256) | {"SHA256SUMS"}, "formal evidence file universe differs")
    output: dict[str, dict[str, Any]] = {}
    for name, expected in EVIDENCE_SHA256.items():
        path = root / name
        require(path.is_file() and path.stat().st_size > 0, f"evidence is missing: {name}")
        observed = sha256_file(path)
        require(observed == expected, f"evidence SHA-256 mismatch: {name}")
        output[name] = {"bytes": path.stat().st_size, "sha256": observed}
    output["SHA256SUMS"] = {
        "bytes": (root / "SHA256SUMS").stat().st_size,
        "sha256": sha256_file(root / "SHA256SUMS"),
    }
    return output


def linear_quantile(values: Iterable[float], probability: float) -> float:
    ordered = sorted(values)
    require(bool(ordered), "bootstrap endpoint is empty")
    index = (len(ordered) - 1) * probability
    lower = math.floor(index)
    upper = math.ceil(index)
    if lower == upper:
        return ordered[lower]
    fraction = index - lower
    return ordered[lower] * (1.0 - fraction) + ordered[upper] * fraction


def validate_bootstrap(path: Path, evaluation: dict[str, Any]) -> dict[str, Any]:
    endpoints: dict[str, list[tuple[int, float]]] = {
        "H1": [],
        "H2": [],
        "guard_v03_descriptive": [],
    }
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require(reader.fieldnames == ["endpoint", "replicate", "delta"], "bootstrap header differs")
        for row in reader:
            endpoint = row["endpoint"]
            require(endpoint in endpoints, f"unexpected bootstrap endpoint: {endpoint}")
            replicate = int(row["replicate"])
            value = float(row["delta"])
            require(math.isfinite(value), "non-finite bootstrap delta")
            endpoints[endpoint].append((replicate, value))
    reports = {
        "H1": evaluation["H1_chain_preserving_candidate_ceiling"],
        "H2": evaluation["H2_v04_guard_minus_baseline"],
        "guard_v03_descriptive": evaluation["descriptive_v04_guard_minus_v03"],
    }
    output: dict[str, Any] = {}
    for endpoint, records in endpoints.items():
        require(len(records) == 20_000, f"{endpoint} bootstrap count differs")
        require([item[0] for item in records] == list(range(1, 20_001)), f"{endpoint} replicate serial differs")
        values = [item[1] for item in records]
        lower = linear_quantile(values, 0.025)
        upper = linear_quantile(values, 0.975)
        require(math.isclose(lower, reports[endpoint]["ci_lower"], rel_tol=0.0, abs_tol=1e-15), f"{endpoint} lower CI differs")
        require(math.isclose(upper, reports[endpoint]["ci_upper"], rel_tol=0.0, abs_tol=1e-15), f"{endpoint} upper CI differs")
        output[endpoint] = {
            "replicates": len(records),
            "ci_lower": lower,
            "ci_upper": upper,
        }
    return output


def validate_candidates(path: Path, evaluation: dict[str, Any]) -> dict[str, Any]:
    expected_header = [
        "candidate_digest",
        "seqid",
        "label_exact_cds",
        "conflict_set_digest",
        "v03_baseline_logit",
        "v03_primary_rank_score",
        "v04_primary_rank_score",
        "v04_topology_abstained",
    ]
    rows = 0
    positives = 0
    digests: set[str] = set()
    seqids: set[str] = set()
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require(reader.fieldnames == expected_header, "candidate header differs")
        for row in reader:
            rows += 1
            digest = row["candidate_digest"]
            require(re.fullmatch(r"[0-9a-f]{64}", digest) is not None, "invalid candidate digest")
            require(digest not in digests, "duplicate candidate digest")
            digests.add(digest)
            label = int(row["label_exact_cds"])
            require(label in (0, 1), "candidate label is not binary")
            positives += label
            seqids.add(row["seqid"])
            require(int(row["v04_topology_abstained"]) in (0, 1), "topology abstention is not binary")
            for field in ("v03_baseline_logit", "v03_primary_rank_score", "v04_primary_rank_score"):
                require(math.isfinite(float(row[field])), f"non-finite candidate score: {field}")
    counts = evaluation["counts"]
    require(rows == counts["candidates"], "candidate row count differs")
    require(positives == counts["positive_exact_cds_candidates"], "positive candidate count differs")
    require(rows - positives == counts["negative_candidates"], "negative candidate count differs")
    require(len(seqids) == counts["target_chromosomes"], "candidate chromosome count differs")
    return {
        "rows": rows,
        "positive": positives,
        "negative": rows - positives,
        "target_chromosomes": len(seqids),
    }


def validate_reviews(evaluation: dict[str, Any]) -> None:
    counts = evaluation["counts"]
    candidates = counts["candidates"]
    positives = counts["positive_exact_cds_candidates"]
    expected_reviewed = {
        "top_0.5pct": math.ceil(candidates * 0.005),
        "top_1pct": math.ceil(candidates * 0.01),
        "top_2pct": math.ceil(candidates * 0.02),
        "top_100": 100,
        "top_250": 250,
        "top_500": 500,
    }
    for estimator in ("baseline", "v03_primary", "v04_guard"):
        records = evaluation["review_budgets"][estimator]
        require(set(records) == set(expected_reviewed), f"review budgets differ: {estimator}")
        for budget, expected in expected_reviewed.items():
            record = records[budget]
            reviewed = int(record["reviewed"])
            true_positive = int(record["true_positive"])
            require(reviewed == expected, f"review count differs: {estimator}/{budget}")
            require(math.isclose(record["precision"], true_positive / reviewed, abs_tol=1e-15), f"review precision differs: {estimator}/{budget}")
            require(math.isclose(record["positive_candidate_recall"], true_positive / positives, abs_tol=1e-15), f"review recall differs: {estimator}/{budget}")
            require(re.fullmatch(r"[0-9a-f]{64}", record["selected_digest_sha256"]) is not None, f"review digest differs: {estimator}/{budget}")


def load_and_validate_evaluation(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        evaluation = json.load(handle)
    require(isinstance(evaluation, dict), "evaluation root is not an object")
    require(evaluation.get("schema_version") == EVALUATION_SCHEMA, "evaluation schema differs")
    require(evaluation.get("evaluation_role") == EVALUATION_ROLE, "evaluation role differs")
    require(evaluation.get("formal_outcome") == FORMAL_OUTCOME, "formal outcome differs")
    require(evaluation.get("confirmatory_pass") is False, "Actinidia result is not formal negative")
    gates = evaluation.get("gates", {})
    require(set(gates) == set(GATE_ORDER), "gate universe differs")
    failed = [gate for gate in GATE_ORDER if gates[gate] is not True]
    require(failed == ["topology_positive_coverage"], "topology coverage is not the unique failed gate")
    counts = evaluation["counts"]
    numerator = int(counts["effective_topology_positive_candidates"])
    denominator = int(counts["positive_exact_cds_candidates"])
    coverage = float(counts["topology_positive_coverage"])
    require((numerator, denominator) == (252, 406), "topology coverage counts differ")
    require(math.isclose(coverage, numerator / denominator, abs_tol=1e-15), "topology coverage ratio differs")
    require(coverage < TOPOLOGY_THRESHOLD, "topology gate unexpectedly passes")
    h1 = evaluation["H1_chain_preserving_candidate_ceiling"]
    require((h1["primary_recovered"], h1["legacy_recovered"], h1["events"]) == (406, 333, 800), "H1 counts differ")
    require(h1["observed_delta"] > 0 and h1["ci_lower"] > 0, "H1 gate is not positive")
    h2 = evaluation["H2_v04_guard_minus_baseline"]
    require(h2["fixed_sequence_status"] == "confirmatory_tested", "H2 fixed-sequence status differs")
    require(h2["replicates_valid"] == h2["replicates_requested"] == 20_000, "H2 bootstrap count differs")
    require(h2["observed_delta"] > 0 and h2["ci_lower"] > 0, "H2 gate is not positive")
    require(evaluation["winner_audit"]["mismatch_count"] == 0, "winner audit differs")
    require(all(evaluation["collateral_gate_by_arm"].values()), "collateral gate differs")
    validate_reviews(evaluation)
    return evaluation


def panel_label(ax: Any, label: str) -> None:
    ax.text(-0.10, 1.07, label, transform=ax.transAxes, fontsize=14, fontweight="bold", va="top")


def source_row(**values: Any) -> dict[str, Any]:
    return {field: values.get(field, "") for field in SOURCE_FIELDS}


def draw_h1(ax: Any, evaluation: dict[str, Any], rows: list[dict[str, Any]]) -> None:
    h1 = evaluation["H1_chain_preserving_candidate_ceiling"]
    events = int(h1["events"])
    recovered = [int(h1["legacy_recovered"]), int(h1["primary_recovered"])]
    rates = [value / events for value in recovered]
    ax.bar([0, 1], rates, color=[GRAY, BLUE], width=0.62)
    ax.set_ylim(0, 0.63)
    ax.set_xticks([0, 1], ["Legacy overlap\nsuppression", "Chain-preserving\nprimary pool"])
    ax.set_ylabel("Exact phased-CDS event recall")
    ax.set_title("H1: candidate recovery ceiling", fontweight="bold")
    ax.grid(axis="y", color=LIGHT_GRAY, linewidth=0.6, zorder=0)
    ax.set_axisbelow(True)
    for index, (count, rate) in enumerate(zip(recovered, rates, strict=True)):
        ax.text(index, rate + 0.015, f"{count}/{events}\n{rate:.3f}", ha="center", fontsize=8.5)
    ax.text(
        0.50,
        0.08,
        f"Paired delta = {h1['observed_delta']:+.5f}\n95% CI [{h1['ci_lower']:+.5f}, {h1['ci_upper']:+.5f}]",
        transform=ax.transAxes,
        ha="center",
        va="bottom",
        color=BLUE,
        fontweight="bold",
        bbox={"boxstyle": "round,pad=0.35", "facecolor": "white", "edgecolor": BLUE},
    )
    rows.extend(
        [
            source_row(panel="a", metric="exact_event_recall", estimator="legacy_union", value=rates[0], numerator=recovered[0], denominator=events, status="descriptive"),
            source_row(panel="a", metric="exact_event_recall", estimator="primary_union", value=rates[1], numerator=recovered[1], denominator=events, status="descriptive"),
            source_row(panel="a", metric="paired_recall_delta", estimator="primary_minus_legacy", value=h1["observed_delta"], ci_lower=h1["ci_lower"], ci_upper=h1["ci_upper"], status="pass", note="20,000-replicate paired event bootstrap"),
        ]
    )


def draw_h2(ax: Any, evaluation: dict[str, Any], rows: list[dict[str, Any]]) -> None:
    metrics = evaluation["ranking_metrics"]
    h2 = evaluation["H2_v04_guard_minus_baseline"]
    names = ("baseline", "v03_primary", "v04_guard")
    labels = ("Frozen baseline", "v0.3 primary", "v0.4 guard")
    colors = (GRAY, GREEN, PURPLE)
    values = [float(metrics[name]["average_precision"]) for name in names]
    y = [2.25, 1.35, 0.45]
    ax.hlines(y, xmin=0.37, xmax=values, color=colors, linewidth=2.0, alpha=0.7)
    ax.scatter(values, y, color=colors, s=66, zorder=3)
    ax.set_xlim(0.37, 0.438)
    ax.set_ylim(-1.05, 2.8)
    ax.set_yticks(y, labels)
    ax.set_title("H2: untouched review ranking", fontweight="bold")
    ax.grid(axis="x", color=LIGHT_GRAY, linewidth=0.6)
    ax.tick_params(axis="x", bottom=False, labelbottom=False)
    ax.text(0.50, 0.315, "Candidate average precision", transform=ax.transAxes, ha="center", fontsize=8.2)
    for x, yy, color in zip(values, y, colors, strict=True):
        ax.text(x + 0.0012, yy, f"{x:.6f}", color=color, va="center", fontsize=8.3)

    inset = ax.inset_axes([0.18, 0.035, 0.72, 0.22])
    observed = float(h2["observed_delta"])
    lower = float(h2["ci_lower"])
    upper = float(h2["ci_upper"])
    inset.errorbar(
        observed,
        0,
        xerr=[[observed - lower], [upper - observed]],
        fmt="o",
        color=PURPLE,
        ecolor=PURPLE,
        capsize=4,
        markersize=5,
    )
    inset.axvline(0, color=GRAY, linestyle="--", linewidth=0.8)
    inset.set_xlim(-0.005, 0.065)
    inset.set_ylim(-0.45, 0.45)
    inset.set_yticks([])
    inset.set_xlabel("v0.4 - baseline AP (95% chromosome-bootstrap CI)", fontsize=7.5, labelpad=1)
    inset.tick_params(axis="x", labelsize=7)
    inset.spines["top"].set_visible(False)
    inset.spines["right"].set_visible(False)
    inset.spines["left"].set_visible(False)
    inset.text(observed, 0.27, f"{observed:+.4f}\n[{lower:+.4f}, {upper:+.4f}]", ha="center", va="bottom", fontsize=7.2, color=PURPLE)

    for name, value in zip(names, values, strict=True):
        rows.append(source_row(panel="b", metric="average_precision", estimator=name, value=value, status="descriptive" if name == "v03_primary" else "confirmatory_input"))
    rows.append(source_row(panel="b", metric="chromosome_bootstrap_delta_AP", estimator="v04_guard_minus_baseline", value=observed, numerator=h2["replicates_valid"], denominator=h2["replicates_requested"], ci_lower=lower, ci_upper=upper, status="pass", note="fixed-sequence H2"))


def draw_review(ax: Any, evaluation: dict[str, Any], rows: list[dict[str, Any]]) -> None:
    budgets = evaluation["review_budgets"]
    estimators = (
        ("baseline", "Frozen baseline", GRAY, "o", "-"),
        ("v03_primary", "v0.3 primary", GREEN, "s", "--"),
        ("v04_guard", "v0.4 guard", PURPLE, "^", "-"),
    )
    for estimator, label, color, marker, linestyle in estimators:
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
            marker=marker,
            markersize=5,
            linewidth=1.6,
            linestyle=linestyle,
            color=color,
            markerfacecolor="white" if estimator == "v04_guard" else color,
            label=label,
        )
        for reviewed, positives, key, precision in points:
            rows.append(source_row(panel="c", metric="review_true_positives", estimator=estimator, value=positives, numerator=positives, denominator=reviewed, status="fixed_budget", note=f"{key}; precision={precision:.12g}"))
    ax.set_xlim(0, 1460)
    ax.set_ylim(0, 285)
    ax.set_xlabel("Candidates reviewed")
    ax.set_ylabel("Exact positives retrieved")
    ax.set_title("Fixed review-budget yield", fontweight="bold")
    ax.grid(color=LIGHT_GRAY, linewidth=0.6)
    ax.legend(frameon=False, loc="lower right")
    ax.text(0.03, 0.91, "v0.3 and v0.4 retrieve identical TP counts\nat all six budgets", transform=ax.transAxes, color=PURPLE, fontsize=7.8, va="top")


def draw_gates(ax: Any, evaluation: dict[str, Any], rows: list[dict[str, Any]]) -> None:
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.set_title("Predeclared confirmatory gates", fontweight="bold")
    gates = evaluation["gates"]
    counts = evaluation["counts"]
    for index, gate in enumerate(GATE_ORDER):
        column = 0 if index < 6 else 1
        row_index = index if index < 6 else index - 6
        x = 0.02 + column * 0.49
        y = 0.88 - row_index * 0.115
        passed = gates[gate] is True
        color = GREEN if passed else RED
        ax.text(x, y, "PASS" if passed else "FAIL", color=color, fontsize=7.3, fontweight="bold", va="center")
        ax.text(x + 0.088, y, GATE_LABELS[gate], fontsize=7.6, fontweight="bold" if not passed else "normal", va="center", color=RED if not passed else "black")
        numerator: Any = ""
        denominator: Any = ""
        threshold: Any = ""
        value: Any = int(passed)
        note = ""
        if gate == "topology_positive_coverage":
            numerator = counts["effective_topology_positive_candidates"]
            denominator = counts["positive_exact_cds_candidates"]
            threshold = TOPOLOGY_THRESHOLD
            value = counts["topology_positive_coverage"]
            note = f"{numerator}/{denominator} = {100 * value:.2f}% < 70%"
            ax.text(x + 0.088, y - 0.037, note, fontsize=7.1, color=RED, va="center")
        rows.append(source_row(panel="d", metric=gate, estimator="predeclared_gate", value=value, numerator=numerator, denominator=denominator, threshold=threshold, status="pass" if passed else "fail", note=note))

    ax.text(
        0.5,
        0.105,
        "FORMAL NEGATIVE EXTERNAL RESULT",
        ha="center",
        va="center",
        fontsize=10.5,
        fontweight="bold",
        color=RED,
        bbox={"boxstyle": "round,pad=0.45", "facecolor": PALE_RED, "edgecolor": RED, "linewidth": 1.2},
    )
    ax.text(0.5, 0.035, "H1 and H2 pass; topology coverage is the only failed gate.", ha="center", va="center", fontsize=7.8, color=RED)
    rows.append(source_row(panel="d", metric="formal_outcome", estimator="actinidia_external_v0.5", value=FORMAL_OUTCOME, status="formal_negative", note="unique failed gate: topology_positive_coverage"))


def write_source_data(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=SOURCE_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evidence-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--replace-existing", action="store_true", help="atomically replace a complete prior figure bundle")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    evidence_root = args.evidence_dir.resolve(strict=True)
    evidence_manifest = verify_evidence(evidence_root)
    evaluation = load_and_validate_evaluation(evidence_root / "evaluation.json")
    bootstrap_audit = validate_bootstrap(evidence_root / "bootstrap_deltas.tsv", evaluation)
    candidate_audit = validate_candidates(evidence_root / "candidates.tsv", evaluation)
    output = args.output_dir.resolve()
    working = Path(f"{output}.working")
    output_names = {
        "figure_actinidia_external_v0.5.png",
        "figure_actinidia_external_v0.5.pdf",
        "figure_actinidia_external_v0.5.svg",
        "figure_actinidia_external_v0.5_source_data.tsv",
        "figure_manifest.json",
        "SHA256SUMS",
    }
    if output.exists():
        require(args.replace_existing, f"refusing to overwrite: {output}")
        require(output.is_dir(), f"existing output is not a directory: {output}")
        require({path.name for path in output.iterdir()} == output_names, "existing figure bundle file universe differs")
        prior_manifest = json.loads((output / "figure_manifest.json").read_text(encoding="utf-8"))
        require(prior_manifest.get("schema_version") == SCHEMA_VERSION, "existing figure bundle schema differs")
    require(not working.exists(), f"refusing existing working directory: {working}")
    working.mkdir(parents=True)
    try:
        source_rows: list[dict[str, Any]] = []
        fig, axes = plt.subplots(2, 2, figsize=(12.2, 8.7))
        fig.subplots_adjust(left=0.065, right=0.985, bottom=0.075, top=0.865, hspace=0.34, wspace=0.24)
        draw_h1(axes[0, 0], evaluation, source_rows)
        draw_h2(axes[0, 1], evaluation, source_rows)
        draw_review(axes[1, 0], evaluation, source_rows)
        draw_gates(axes[1, 1], evaluation, source_rows)
        for label, ax in zip(("a", "b", "c", "d"), axes.ravel(), strict=True):
            panel_label(ax, label)
        fig.suptitle("Actinidia external validation under the predeclared v0.5 protocol", y=0.975, fontsize=14, fontweight="bold")
        fig.text(0.5, 0.932, "Chain preservation and ranking improve, but topology coverage fails the formal gate.", ha="center", color=RED, fontsize=9)

        os.environ["SOURCE_DATE_EPOCH"] = "0"
        fixed_date = datetime.datetime(1970, 1, 1, tzinfo=datetime.timezone.utc)
        outputs = {
            "png": working / "figure_actinidia_external_v0.5.png",
            "pdf": working / "figure_actinidia_external_v0.5.pdf",
            "svg": working / "figure_actinidia_external_v0.5.svg",
        }
        fig.savefig(outputs["png"], dpi=300, facecolor="white", metadata={"Software": "PloidyPatch"})
        fig.savefig(outputs["pdf"], facecolor="white", metadata={"Creator": "PloidyPatch", "CreationDate": fixed_date, "ModDate": fixed_date})
        fig.savefig(outputs["svg"], facecolor="white", metadata={"Creator": "PloidyPatch", "Date": "1970-01-01T00:00:00Z"})
        plt.close(fig)

        source_data = working / "figure_actinidia_external_v0.5_source_data.tsv"
        write_source_data(source_data, source_rows)
        failed_gates = [gate for gate in GATE_ORDER if evaluation["gates"][gate] is not True]
        manifest = {
            "schema_version": SCHEMA_VERSION,
            "title": "Actinidia external validation under the predeclared v0.5 protocol",
            "evidence": evidence_manifest,
            "formal_outcome": evaluation["formal_outcome"],
            "confirmatory_pass": evaluation["confirmatory_pass"],
            "gate_count": len(GATE_ORDER),
            "passed_gate_count": len(GATE_ORDER) - len(failed_gates),
            "failed_gates": failed_gates,
            "topology_coverage": {
                "numerator": evaluation["counts"]["effective_topology_positive_candidates"],
                "denominator": evaluation["counts"]["positive_exact_cds_candidates"],
                "value": evaluation["counts"]["topology_positive_coverage"],
                "threshold": TOPOLOGY_THRESHOLD,
            },
            "bootstrap_audit": bootstrap_audit,
            "candidate_audit": candidate_audit,
            "matplotlib_version": matplotlib.__version__,
            "generator_sha256": sha256_file(Path(__file__).resolve()),
            "source_data": {
                "file_name": source_data.name,
                "rows": len(source_rows),
                "sha256": sha256_file(source_data),
            },
            "outputs": {
                key: {"file_name": path.name, "bytes": path.stat().st_size, "sha256": sha256_file(path)}
                for key, path in outputs.items()
            },
        }
        manifest_path = working / "figure_manifest.json"
        manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        with (working / "SHA256SUMS").open("x", encoding="utf-8", newline="") as handle:
            for path in sorted(item for item in working.iterdir() if item.is_file() and item.name != "SHA256SUMS"):
                handle.write(f"{sha256_file(path)}  {path.name}\n")
        if output.exists():
            for path in working.iterdir():
                os.replace(path, output / path.name)
            working.rmdir()
        else:
            working.rename(output)
    except Exception:
        if working.exists():
            shutil.rmtree(working)
        raise
    print(json.dumps({"figure_root": os.fspath(output), "formal_outcome": evaluation["formal_outcome"], "source_rows": len(source_rows)}, sort_keys=True))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except FigureError as error:
        print(f"Actinidia figure failed: {error}", file=sys.stderr)
        raise SystemExit(1) from error
