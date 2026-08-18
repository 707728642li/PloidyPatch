#!/usr/bin/env python3
from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
from sklearn.metrics import average_precision_score


def _load_evaluator():
    path = Path(__file__).with_name("evaluate_frozen_homeolog_ranker_v0.2.py")
    spec = importlib.util.spec_from_file_location("ploidypatch_rank_evaluator", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import evaluator: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main() -> int:
    evaluator = _load_evaluator()
    labels = np.asarray([1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0])
    scores = np.asarray([0.9, 0.9, 0.8, 0.7, 0.7, 0.6, 0.5, 0.5, 0.4, 0.3, 0.2, 0.1])
    group_codes = np.asarray([0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3])
    counts = np.asarray(
        [
            [1, 1, 1, 1],
            [2, 0, 1, 1],
            [0, 2, 2, 0],
            [3, 1, 0, 2],
        ]
    )
    contract = evaluator._weighted_ap_contract(labels, scores, group_codes, 4)
    observed = evaluator._weighted_ap_batch(counts, *contract)
    expected = np.asarray(
        [
            average_precision_score(labels, scores, sample_weight=row[group_codes])
            for row in counts
        ]
    )
    if not np.allclose(observed, expected, rtol=0, atol=1e-12):
        raise AssertionError(f"weighted AP mismatch: {observed=} {expected=}")
    report = evaluator._group_bootstrap_delta(
        labels,
        scores,
        scores - np.asarray([0.05, 0, 0.02, 0, 0.01, 0, 0, 0.03, 0, 0.02, 0, 0]),
        np.asarray(["chr1", "chr1", "chr1", "chr2", "chr2", "chr2", "chr3", "chr3", "chr3", "chr4", "chr4", "chr4"]),
        replicates=128,
        seed=20260829,
    )
    if report["replicates_requested"] != 128 or report["replicates_valid"] < 1:
        raise AssertionError(report)
    print("weighted AP and chromosome-group bootstrap numeric verification passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
