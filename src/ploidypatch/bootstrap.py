from __future__ import annotations

import hashlib
import json
import random
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable

from . import __version__


BOOTSTRAP_SCHEMA_VERSION = "ploidypatch.paired_event_bootstrap.v1"
INDEPENDENT_BOOTSTRAP_SCHEMA_VERSION = (
    "ploidypatch.independent_event_bootstrap.v1"
)
SCORE_SCHEMA_VERSION = "ploidypatch.annotation_repair_score.v5"
EVENT_BOOLEAN_METRICS = frozenset(
    {
        "complete_cds_chain_recovery",
        "complete_transcript_recovery",
        "complete_error_removal",
        "exact_cds_gene_grouping",
        "exact_gene_grouping",
    }
)


def _file_sha256(path: str | Path, block_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        while block := handle.read(block_size):
            digest.update(block)
    return digest.hexdigest()


def _quantile(values: list[float], probability: float) -> float:
    if not values:
        raise ValueError("Cannot take a quantile of an empty bootstrap sample")
    if not 0 <= probability <= 1:
        raise ValueError("Quantile probability must be within [0, 1]")
    ordered = sorted(values)
    position = probability * (len(ordered) - 1)
    lower = int(position)
    upper = min(lower + 1, len(ordered) - 1)
    fraction = position - lower
    return ordered[lower] * (1 - fraction) + ordered[upper] * fraction


def _read_score(
    label: str, path: str | Path, metric: str
) -> tuple[dict[str, tuple[str, bool]], dict[str, Any]]:
    score_path = Path(path)
    report = json.loads(score_path.read_text(encoding="utf-8"))
    if report.get("schema_version") != SCORE_SCHEMA_VERSION:
        raise ValueError(f"Unsupported score schema for {label}: {score_path}")
    if report.get("quality_gate", {}).get("grade") != "pass":
        raise ValueError(f"Score quality gate did not pass for {label}")
    details = report.get("event_details")
    if not isinstance(details, list) or not details:
        raise ValueError(f"Score lacks event details for {label}")
    events: dict[str, tuple[str, bool]] = {}
    for index, row in enumerate(details, start=1):
        event_id = row.get("event_id")
        event_type = row.get("event_type")
        value = row.get(metric)
        if not isinstance(event_id, str) or not event_id or event_id in events:
            raise ValueError(f"Invalid or duplicate event ID for {label} at row {index}")
        if not isinstance(event_type, str) or not event_type:
            raise ValueError(f"Missing event type for {label} at row {index}")
        if not isinstance(value, bool):
            raise ValueError(f"Metric {metric} is not Boolean for {label} at row {index}")
        events[event_id] = (event_type, value)
    return events, {
        "label": label,
        "file_name": score_path.name,
        "sha256": _file_sha256(score_path),
        "events": len(events),
    }


def paired_event_bootstrap(
    *,
    score_inputs: Iterable[tuple[str, str | Path]],
    output_json_path: str | Path,
    metric: str = "complete_cds_chain_recovery",
    replicates: int = 10_000,
    seed: int = 20260807,
    alpha: float = 0.05,
) -> dict[str, Any]:
    """Bootstrap event recovery and paired method deltas within event type."""

    if metric not in EVENT_BOOLEAN_METRICS:
        raise ValueError(f"Unsupported event bootstrap metric: {metric}")
    if replicates < 100:
        raise ValueError("At least 100 bootstrap replicates are required")
    if not 0 < alpha < 1:
        raise ValueError("alpha must be between zero and one")
    output = Path(output_json_path)
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite bootstrap output: {output}")
    inputs = list(score_inputs)
    if len(inputs) < 2:
        raise ValueError("At least two score inputs are required")
    labels = [label for label, _ in inputs]
    if any(not label for label in labels) or len(set(labels)) != len(labels):
        raise ValueError("Bootstrap labels must be non-empty and unique")

    methods: dict[str, dict[str, tuple[str, bool]]] = {}
    input_manifests: list[dict[str, Any]] = []
    for label, path in inputs:
        methods[label], manifest = _read_score(label, path, metric)
        input_manifests.append(manifest)
    first_label = labels[0]
    event_ids = set(methods[first_label])
    for label in labels[1:]:
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
    method_samples: dict[str, list[float]] = {label: [] for label in labels}
    delta_samples: dict[tuple[str, str], list[float]] = {
        (left, right): []
        for left_index, left in enumerate(labels)
        for right in labels[left_index + 1 :]
    }
    total_events = len(event_ids)
    for _ in range(replicates):
        sampled_ids = [
            rng.choice(stratum_ids)
            for stratum_ids in strata.values()
            for _ in range(len(stratum_ids))
        ]
        replicate_rates = {
            label: sum(methods[label][event_id][1] for event_id in sampled_ids)
            / total_events
            for label in labels
        }
        for label, rate in replicate_rates.items():
            method_samples[label].append(rate)
        for pair in delta_samples:
            delta_samples[pair].append(
                replicate_rates[pair[0]] - replicate_rates[pair[1]]
            )

    lower_probability = alpha / 2
    upper_probability = 1 - alpha / 2
    method_reports = []
    for label in labels:
        observed = sum(value for _, value in methods[label].values()) / total_events
        method_reports.append(
            {
                "label": label,
                "events": total_events,
                "successes": sum(value for _, value in methods[label].values()),
                "observed_rate": observed,
                "bootstrap_mean": sum(method_samples[label]) / replicates,
                "ci_lower": _quantile(method_samples[label], lower_probability),
                "ci_upper": _quantile(method_samples[label], upper_probability),
            }
        )
    pair_reports = []
    for (left, right), samples in delta_samples.items():
        observed_delta = next(
            row["observed_rate"] for row in method_reports if row["label"] == left
        ) - next(
            row["observed_rate"] for row in method_reports if row["label"] == right
        )
        pair_reports.append(
            {
                "left": left,
                "right": right,
                "delta_definition": "left_minus_right",
                "observed_delta": observed_delta,
                "bootstrap_mean_delta": sum(samples) / replicates,
                "ci_lower": _quantile(samples, lower_probability),
                "ci_upper": _quantile(samples, upper_probability),
                "probability_delta_gt_zero": sum(value > 0 for value in samples)
                / replicates,
                "probability_delta_eq_zero": sum(value == 0 for value in samples)
                / replicates,
            }
        )
    report = {
        "schema_version": BOOTSTRAP_SCHEMA_VERSION,
        "generator": {"name": "PloidyPatch", "version": __version__},
        "metric": metric,
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
                event_type: len(ids) for event_type, ids in sorted(strata.items())
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


def independent_event_bootstrap(
    *,
    score_inputs: Iterable[tuple[str, str | Path]],
    output_json_path: str | Path,
    metric: str = "complete_cds_chain_recovery",
    replicates: int = 10_000,
    seed: int = 20260807,
    alpha: float = 0.05,
) -> dict[str, Any]:
    """Bootstrap recovery rates from independent benchmark event sets.

    Each score is resampled within its own event-type strata.  Unlike
    :func:`paired_event_bootstrap`, event identifiers need not be shared and
    method deltas are therefore independent-sample contrasts.
    """

    if metric not in EVENT_BOOLEAN_METRICS:
        raise ValueError(f"Unsupported event bootstrap metric: {metric}")
    if replicates < 100:
        raise ValueError("At least 100 bootstrap replicates are required")
    if not 0 < alpha < 1:
        raise ValueError("alpha must be between zero and one")
    output = Path(output_json_path)
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite bootstrap output: {output}")
    inputs = list(score_inputs)
    if not inputs:
        raise ValueError("At least one score input is required")
    labels = [label for label, _ in inputs]
    if any(not label for label in labels) or len(set(labels)) != len(labels):
        raise ValueError("Bootstrap labels must be non-empty and unique")

    methods: dict[str, dict[str, tuple[str, bool]]] = {}
    strata: dict[str, dict[str, list[str]]] = {}
    input_manifests: list[dict[str, Any]] = []
    for label, path in inputs:
        methods[label], manifest = _read_score(label, path, metric)
        label_strata: dict[str, list[str]] = defaultdict(list)
        for event_id, (event_type, _) in methods[label].items():
            label_strata[event_type].append(event_id)
        for event_ids in label_strata.values():
            event_ids.sort()
        strata[label] = dict(label_strata)
        manifest["event_types"] = {
            event_type: len(event_ids)
            for event_type, event_ids in sorted(label_strata.items())
        }
        input_manifests.append(manifest)

    rng = random.Random(seed)
    method_samples: dict[str, list[float]] = {label: [] for label in labels}
    pairs = [
        (left, right)
        for left_index, left in enumerate(labels)
        for right in labels[left_index + 1 :]
    ]
    delta_samples: dict[tuple[str, str], list[float]] = {
        pair: [] for pair in pairs
    }
    for _ in range(replicates):
        rates: dict[str, float] = {}
        for label in labels:
            sampled_ids = [
                rng.choice(event_ids)
                for event_ids in strata[label].values()
                for _ in range(len(event_ids))
            ]
            rates[label] = sum(
                methods[label][event_id][1] for event_id in sampled_ids
            ) / len(sampled_ids)
            method_samples[label].append(rates[label])
        for pair in pairs:
            delta_samples[pair].append(rates[pair[0]] - rates[pair[1]])

    lower_probability = alpha / 2
    upper_probability = 1 - alpha / 2
    method_reports: list[dict[str, Any]] = []
    for label in labels:
        successes = sum(value for _, value in methods[label].values())
        events = len(methods[label])
        samples = method_samples[label]
        method_reports.append(
            {
                "label": label,
                "events": events,
                "successes": successes,
                "observed_rate": successes / events,
                "bootstrap_mean": sum(samples) / replicates,
                "ci_lower": _quantile(samples, lower_probability),
                "ci_upper": _quantile(samples, upper_probability),
            }
        )
    observed = {row["label"]: row["observed_rate"] for row in method_reports}
    difference_reports: list[dict[str, Any]] = []
    for left, right in pairs:
        samples = delta_samples[(left, right)]
        difference_reports.append(
            {
                "left": left,
                "right": right,
                "delta_definition": "left_minus_right",
                "observed_delta": observed[left] - observed[right],
                "bootstrap_mean_delta": sum(samples) / replicates,
                "ci_lower": _quantile(samples, lower_probability),
                "ci_upper": _quantile(samples, upper_probability),
                "probability_delta_gt_zero": sum(x > 0 for x in samples)
                / replicates,
                "probability_delta_eq_zero": sum(x == 0 for x in samples)
                / replicates,
            }
        )

    report = {
        "schema_version": INDEPENDENT_BOOTSTRAP_SCHEMA_VERSION,
        "generator": {"name": "PloidyPatch", "version": __version__},
        "metric": metric,
        "parameters": {
            "replicates": replicates,
            "seed": seed,
            "alpha": alpha,
            "interval": "percentile",
            "resampling_unit": "event",
            "stratified_by_event_type": True,
            "paired_method_resampling": False,
            "independent_score_resampling": True,
        },
        "inputs": input_manifests,
        "methods": method_reports,
        "independent_differences": difference_reports,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x", encoding="utf-8", newline="") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return report
