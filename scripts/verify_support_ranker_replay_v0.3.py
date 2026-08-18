#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from typing import Any


TOLERANCE = 1e-12


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Missing table header: {path}")
        return list(reader)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--canonical-predictions", required=True)
    parser.add_argument(
        "--score", action="append", required=True, help="SPECIES=portable_scores.tsv"
    )
    parser.add_argument("--model-json", required=True)
    parser.add_argument("--output-json", required=True)
    args = parser.parse_args()

    canonical_path = Path(args.canonical_predictions)
    model_path = Path(args.model_json)
    output_path = Path(args.output_json)
    if output_path.exists():
        raise FileExistsError(output_path)
    for path in (canonical_path, model_path):
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(path)

    canonical: dict[str, dict[str, dict[str, str]]] = {}
    for row in read_rows(canonical_path):
        species = row["dataset"]
        if species in {"cotton", "maize"}:
            canonical.setdefault(species, {})[row["candidate_digest"]] = row

    report: dict[str, Any] = {}
    for value in args.score:
        species, separator, path_value = value.partition("=")
        if not separator or not species:
            raise ValueError("--score requires SPECIES=PATH")
        path = Path(path_value)
        manifest_path = Path(str(path) + ".manifest.json")
        for required in (path, manifest_path):
            if not required.is_file() or required.stat().st_size == 0:
                raise FileNotFoundError(required)
        with manifest_path.open("r", encoding="utf-8") as handle:
            manifest = json.load(handle)
        if (
            manifest.get("truth_access") is not False
            or manifest.get("inputs", {}).get("model") != sha256(model_path)
            or manifest.get("outputs", {}).get("scores", {}).get("sha256")
            != sha256(path)
        ):
            raise ValueError(f"Portable score manifest failed for {species}")
        observed_rows = read_rows(path)
        observed = {row["candidate_digest"]: row for row in observed_rows}
        if len(observed) != len(observed_rows) or not observed:
            raise ValueError(f"Empty or duplicate portable scores for {species}")
        if species not in canonical or set(observed) != set(canonical[species]):
            raise ValueError(f"Canonical/portable candidate universes differ for {species}")
        baseline_delta = max(
            abs(
                float(observed[digest]["v03_baseline_logit"])
                - float(canonical[species][digest]["baseline"])
            )
            for digest in observed
        )
        primary_delta = max(
            abs(
                float(observed[digest]["v03_primary_rank_score"])
                - float(canonical[species][digest]["offset_support_conditioned"])
            )
            for digest in observed
        )
        observed_order = sorted(
            observed,
            key=lambda digest: (
                -float(observed[digest]["v03_primary_rank_score"]),
                digest,
            ),
        )
        canonical_order = sorted(
            canonical[species],
            key=lambda digest: (
                -float(
                    canonical[species][digest]["offset_support_conditioned"]
                ),
                digest,
            ),
        )
        rank_order_equal = observed_order == canonical_order
        if (
            baseline_delta > TOLERANCE
            or primary_delta > TOLERANCE
            or not rank_order_equal
        ):
            raise ValueError(f"Portable ranker replay failed for {species}")
        report[species] = {
            "candidates": len(observed),
            "max_absolute_baseline_delta": baseline_delta,
            "max_absolute_primary_delta": primary_delta,
            "rank_order_equal": rank_order_equal,
            "scores_sha256": sha256(path),
            "manifest_sha256": sha256(manifest_path),
        }

    if set(report) != {"cotton", "maize"}:
        raise ValueError("Replay audit requires exactly cotton and maize")
    output: dict[str, Any] = {
        "schema_version": "ploidypatch.support_ranker_replay_audit.v1",
        "status": "pass",
        "tolerance": TOLERANCE,
        "canonical_predictions": {
            "path": str(canonical_path),
            "sha256": sha256(canonical_path),
        },
        "model": {"path": str(model_path), "sha256": sha256(model_path)},
        "species": report,
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("x", encoding="utf-8", newline="") as handle:
        json.dump(output, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(json.dumps(output, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

