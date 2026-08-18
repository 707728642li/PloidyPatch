from __future__ import annotations

import csv
import hashlib
import json
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable

from . import __version__
from .baseline import _file_sha256


PAIR_INTERSECTION_SCHEMA_VERSION = "ploidypatch.copy_pair_intersection.v1"
PAIR_INTERSECTION_MANIFEST_SCHEMA_VERSION = (
    "ploidypatch.copy_pair_intersection_manifest.v1"
)
SAFE_LABEL = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")


def _parse_inputs(values: Iterable[str]) -> list[tuple[str, Path]]:
    parsed: list[tuple[str, Path]] = []
    labels: set[str] = set()
    for value in values:
        if "=" not in value:
            raise ValueError("Pair evidence must use LABEL=PATH")
        label, path_text = value.split("=", 1)
        path = Path(path_text)
        if not SAFE_LABEL.fullmatch(label) or label in labels:
            raise ValueError(f"Unsafe or duplicate pair-evidence label: {label!r}")
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f"Missing or empty pair-evidence input: {path}")
        labels.add(label)
        parsed.append((label, path))
    if len(parsed) < 2:
        raise ValueError("Pair intersection requires at least two evidence inputs")
    return parsed


def intersect_copy_pair_evidence(
    *,
    pair_inputs: Iterable[str],
    output_pair_tsv_path: str | Path,
    decisions_tsv_path: str | Path,
    pair_set_label: str,
    require_reciprocal_unique: bool = True,
) -> dict[str, Any]:
    """Intersect exact unordered gene pairs from independent evidence tables."""

    if not SAFE_LABEL.fullmatch(pair_set_label):
        raise ValueError(f"Unsafe pair-set label: {pair_set_label!r}")
    inputs = _parse_inputs(pair_inputs)
    output = Path(output_pair_tsv_path)
    decisions = Path(decisions_tsv_path)
    manifest_path = Path(str(output) + ".manifest.json")
    collisions = [path for path in (output, decisions, manifest_path) if path.exists()]
    if collisions:
        raise FileExistsError(
            "Refusing to overwrite pair-intersection artifact(s): "
            + ", ".join(str(path) for path in collisions)
        )

    evidence: dict[tuple[str, str], dict[str, str]] = defaultdict(dict)
    input_manifest: dict[str, dict[str, Any]] = {}
    for label, path in inputs:
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames is None or not {"gene_id_a", "gene_id_b"} <= set(
                reader.fieldnames
            ):
                raise ValueError(f"Pair evidence {label} lacks gene_id_a/gene_id_b")
            rows = list(reader)
        observed: set[tuple[str, str]] = set()
        for line_number, row in enumerate(rows, start=2):
            gene_a = row.get("gene_id_a", "")
            gene_b = row.get("gene_id_b", "")
            if not gene_a or not gene_b or gene_a == gene_b:
                raise ValueError(f"Invalid {label} pair on line {line_number}")
            pair = tuple(sorted((gene_a, gene_b)))
            if pair in observed:
                raise ValueError(f"Duplicate unordered pair in {label} on line {line_number}")
            observed.add(pair)
            evidence[pair][label] = row.get("pair_id", "")
        input_manifest[label] = {
            "file_name": path.name,
            "sha256": _file_sha256(path),
            "rows": len(rows),
        }

    labels = tuple(label for label, _ in inputs)
    all_labels = frozenset(labels)
    exact_intersection = {
        pair for pair, support in evidence.items() if frozenset(support) == all_labels
    }
    partners: dict[str, set[str]] = defaultdict(set)
    for gene_a, gene_b in exact_intersection:
        partners[gene_a].add(gene_b)
        partners[gene_b].add(gene_a)

    fields = (
        "pair_id",
        "gene_id_a",
        "gene_id_b",
        "evidence_count",
        "evidence_sources",
        "source_pair_ids",
        "relationship_scope",
    )
    decision_fields = (*fields, "status", "reason", "missing_evidence_sources")
    accepted: list[dict[str, str]] = []
    decision_rows: list[dict[str, str]] = []
    for gene_a, gene_b in sorted(evidence):
        pair = (gene_a, gene_b)
        support = evidence[pair]
        supported_labels = tuple(label for label in labels if label in support)
        missing = tuple(label for label in labels if label not in support)
        digest = hashlib.sha256(f"{gene_a}\0{gene_b}".encode()).hexdigest()[:24]
        row = {
            "pair_id": f"{pair_set_label}:{digest}",
            "gene_id_a": gene_a,
            "gene_id_b": gene_b,
            "evidence_count": str(len(supported_labels)),
            "evidence_sources": ",".join(supported_labels),
            "source_pair_ids": ";".join(
                f"{label}={support[label]}" for label in supported_labels
            ),
            "relationship_scope": "exact_unordered_pair_evidence_intersection",
        }
        if pair not in exact_intersection:
            status, reason = "rejected", "missing_required_evidence"
        elif require_reciprocal_unique and (
            len(partners[gene_a]) != 1 or len(partners[gene_b]) != 1
        ):
            status, reason = "rejected", "nonreciprocal_intersection_partner"
        else:
            status, reason = "accepted", "all_required_pair_evidence_exact_match"
            accepted.append(row)
        decision_rows.append(
            {
                **row,
                "status": status,
                "reason": reason,
                "missing_evidence_sources": ",".join(missing),
            }
        )

    output.parent.mkdir(parents=True, exist_ok=True)
    decisions.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(accepted)
    with decisions.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=decision_fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(decision_rows)

    reasons = Counter(row["reason"] for row in decision_rows)
    manifest: dict[str, Any] = {
        "schema_version": PAIR_INTERSECTION_MANIFEST_SCHEMA_VERSION,
        "pair_schema_version": PAIR_INTERSECTION_SCHEMA_VERSION,
        "generator": {"name": "PloidyPatch", "version": __version__},
        "access": "evaluator_only",
        "inputs": input_manifest,
        "parameters": {
            "pair_set_label": pair_set_label,
            "evidence_sources_in_order": list(labels),
            "pair_identity": "exact_unordered_gene_feature_ID_pair",
            "required_evidence_count": len(labels),
            "require_reciprocal_unique": require_reciprocal_unique,
        },
        "counts": {
            "union_pairs": len(evidence),
            "exact_intersection_pairs_before_reciprocal_gate": len(exact_intersection),
            "accepted_pairs": len(accepted),
            "decision_reason_counts": dict(sorted(reasons.items())),
        },
        "outputs": {
            "pairs": {
                "file_name": output.name,
                "sha256": _file_sha256(output),
                "rows": len(accepted),
            },
            "decisions": {
                "file_name": decisions.name,
                "sha256": _file_sha256(decisions),
                "rows": len(decision_rows),
            },
        },
        "claim_boundary": (
            "exact intersection of independently generated pair evidence; "
            "not gene-tree proof of homeology"
        ),
    }
    with manifest_path.open("x", encoding="utf-8", newline="") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return manifest
