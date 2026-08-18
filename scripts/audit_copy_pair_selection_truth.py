#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Audit selected copy-pair metadata against hidden perturbation truth"
    )
    parser.add_argument("--selected-pairs", required=True)
    parser.add_argument("--truth", required=True)
    parser.add_argument("--output-json", required=True)
    args = parser.parse_args()
    pair_path = Path(args.selected_pairs)
    truth_path = Path(args.truth)
    output_path = Path(args.output_json)
    if output_path.exists():
        raise FileExistsError(output_path)
    with pair_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "gene_id_a",
            "gene_id_b",
            "collapsed_gene_id",
            "retained_gene_id",
        }
        if reader.fieldnames is None or not required <= set(reader.fieldnames):
            raise ValueError("Selected pair table lacks audit fields")
        rows = list(reader)
    expected: dict[frozenset[str], tuple[str, str]] = {}
    for row in rows:
        pair = frozenset((row["gene_id_a"], row["gene_id_b"]))
        if len(pair) != 2 or pair in expected:
            raise ValueError("Selected pair table contains invalid or duplicate pairs")
        collapsed, retained = row["collapsed_gene_id"], row["retained_gene_id"]
        if frozenset((collapsed, retained)) != pair:
            raise ValueError("Sampler partner metadata disagrees with pair IDs")
        expected[pair] = (collapsed, retained)
    truth = json.loads(truth_path.read_text(encoding="utf-8"))
    events = truth.get("events")
    if not isinstance(events, list):
        raise ValueError("Hidden truth lacks events")
    observed: set[frozenset[str]] = set()
    for event in events:
        details = event.get("details", {})
        pair = frozenset(details.get("pair_gene_ids", ()))
        if pair not in expected or pair in observed:
            raise ValueError("Hidden truth contains an unknown or duplicate copy pair")
        collapsed, retained = expected[pair]
        if event.get("target", {}).get("gene_id") != collapsed:
            raise ValueError("Hidden truth collapsed partner disagrees with sampler")
        if details.get("retained_partner_gene_id") != retained:
            raise ValueError("Hidden truth retained partner disagrees with sampler")
        observed.add(pair)
    if observed != set(expected):
        raise ValueError("Hidden truth does not cover every selected pair")
    report = {
        "schema_version": "ploidypatch.copy_pair_truth_audit.v1",
        "grade": "pass",
        "selected_pair_rows": len(rows),
        "hidden_truth_events": len(events),
        "exact_pair_coverage": len(observed),
        "collapsed_partner_matches": len(observed),
        "retained_partner_matches": len(observed),
        "inputs": {
            "selected_pairs_sha256": _sha256(pair_path),
            "hidden_truth_sha256": _sha256(truth_path),
        },
    }
    with output_path.open("x", encoding="utf-8", newline="") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
