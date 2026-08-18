from __future__ import annotations

import json
import random
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable

from . import __version__
from .baseline import _file_sha256
from .bootstrap import _quantile
from .structure_hypothesis_score import (
    STRUCTURE_HYPOTHESIS_SCORE_SCHEMA_VERSION,
)


STRUCTURE_BOOTSTRAP_SCHEMA_VERSION = (
    "ploidypatch.structure_hypothesis_bootstrap.v1"
)


def paired_structure_hypothesis_bootstrap(
    *,
    score_inputs: Iterable[tuple[str, str | Path]],
    output_json_path: str | Path,
    replicates: int = 10_000,
    seed: int = 20260807,
    alpha: float = 0.05,
) -> dict[str, Any]:
    """Bootstrap complete audit-event detection, stratified by error type."""

    if replicates < 100:
        raise ValueError("At least 100 bootstrap replicates are required")
    if not 0 < alpha < 1:
        raise ValueError("alpha must be between zero and one")
    output = Path(output_json_path)
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite bootstrap output: {output}")

    grouped: dict[str, list[Path]] = defaultdict(list)
    label_order: list[str] = []
    for label, path in score_inputs:
        if not label:
            raise ValueError("Bootstrap labels must be non-empty")
        if label not in grouped:
            label_order.append(label)
        grouped[label].append(Path(path))
    if len(grouped) < 2:
        raise ValueError("At least two bootstrap labels are required")

    methods: dict[str, dict[str, tuple[str, bool]]] = {}
    input_manifests: list[dict[str, Any]] = []
    for label in label_order:
        events: dict[str, tuple[str, bool]] = {}
        for score_path in grouped[label]:
            report = json.loads(score_path.read_text(encoding="utf-8"))
            if (
                report.get("schema_version")
                != STRUCTURE_HYPOTHESIS_SCORE_SCHEMA_VERSION
            ):
                raise ValueError(
                    f"Unsupported structure-hypothesis score: {score_path}"
                )
            if report.get("quality_gate", {}).get("grade") != "pass":
                raise ValueError(f"Score quality gate failed: {score_path}")
            details = report.get("event_details")
            if not isinstance(details, list) or not details:
                raise ValueError(f"Score lacks event details: {score_path}")
            for row in details:
                event_id = row.get("event_id")
                event_type = row.get("event_type")
                recovered = row.get("recovered")
                if (
                    not isinstance(event_id, str)
                    or not event_id
                    or event_id in events
                ):
                    raise ValueError(
                        f"Invalid or duplicate event ID for {label}: {event_id}"
                    )
                if not isinstance(event_type, str) or not event_type:
                    raise ValueError(f"Missing event type for {event_id}")
                if not isinstance(recovered, bool):
                    raise ValueError(f"Non-Boolean recovery for {event_id}")
                events[event_id] = (event_type, recovered)
            input_manifests.append(
                {
                    "label": label,
                    "file_name": score_path.name,
                    "sha256": _file_sha256(score_path),
                    "events": len(details),
                }
            )
        methods[label] = events

    first_label = label_order[0]
    event_ids = set(methods[first_label])
    for label in label_order[1:]:
        if set(methods[label]) != event_ids:
            raise ValueError(f"Paired event IDs differ for {label}")
        for event_id in event_ids:
            if methods[label][event_id][0] != methods[first_label][event_id][0]:
                raise ValueError(f"Event type differs for {event_id} in {label}")

    strata: dict[str, list[str]] = defaultdict(list)
    for event_id, (event_type, _) in methods[first_label].items():
        strata[event_type].append(event_id)
    for values in strata.values():
        values.sort()

    rng = random.Random(seed)
    method_samples = {label: [] for label in label_order}
    pairs = [
        (left, right)
        for left_index, left in enumerate(label_order)
        for right in label_order[left_index + 1 :]
    ]
    delta_samples = {pair: [] for pair in pairs}
    total_events = len(event_ids)
    for _ in range(replicates):
        sampled_ids = [
            rng.choice(stratum_ids)
            for stratum_ids in strata.values()
            for _ in range(len(stratum_ids))
        ]
        rates = {
            label: sum(
                methods[label][event_id][1] for event_id in sampled_ids
            )
            / total_events
            for label in label_order
        }
        for label, rate in rates.items():
            method_samples[label].append(rate)
        for pair in pairs:
            delta_samples[pair].append(rates[pair[0]] - rates[pair[1]])

    low = alpha / 2
    high = 1 - alpha / 2
    method_reports = []
    for label in label_order:
        successes = sum(value for _, value in methods[label].values())
        samples = method_samples[label]
        method_reports.append(
            {
                "label": label,
                "events": total_events,
                "successes": successes,
                "observed_rate": successes / total_events,
                "bootstrap_mean": sum(samples) / replicates,
                "ci_lower": _quantile(samples, low),
                "ci_upper": _quantile(samples, high),
            }
        )
    observed = {row["label"]: row["observed_rate"] for row in method_reports}
    pair_reports = []
    for left, right in pairs:
        samples = delta_samples[(left, right)]
        pair_reports.append(
            {
                "left": left,
                "right": right,
                "delta_definition": "left_minus_right",
                "observed_delta": observed[left] - observed[right],
                "bootstrap_mean_delta": sum(samples) / replicates,
                "ci_lower": _quantile(samples, low),
                "ci_upper": _quantile(samples, high),
                "probability_delta_gt_zero": sum(x > 0 for x in samples)
                / replicates,
                "probability_delta_eq_zero": sum(x == 0 for x in samples)
                / replicates,
            }
        )

    report = {
        "schema_version": STRUCTURE_BOOTSTRAP_SCHEMA_VERSION,
        "generator": {"name": "PloidyPatch", "version": __version__},
        "metric": "complete_structure_hypothesis_event_recovery",
        "parameters": {
            "replicates": replicates,
            "seed": seed,
            "alpha": alpha,
            "interval": "percentile",
            "resampling_unit": "event",
            "stratified_by_event_type": True,
            "paired_method_resampling": True,
        },
        "counts": {
            "events": total_events,
            "event_types": {
                event_type: len(ids)
                for event_type, ids in sorted(strata.items())
            },
        },
        "inputs": input_manifests,
        "methods": method_reports,
        "paired_differences": pair_reports,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x", encoding="utf-8", newline="") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return report
