from __future__ import annotations

import importlib.util
import sys
import csv
from pathlib import Path

import pytest

np = pytest.importorskip("numpy")
pytest.importorskip("scipy")
pytest.importorskip("sklearn")


SCRIPT = Path(__file__).parents[1] / "scripts" / "evaluate_candidate_pool_rankers_v0.3.py"
SPEC = importlib.util.spec_from_file_location("candidate_pool_ranker_v0_3", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def test_read_index_accepts_legacy_consensus_digest(tmp_path: Path) -> None:
    path = tmp_path / "decisions.tsv"
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("consensus_digest", "conflict_set_digest"),
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerow({"consensus_digest": "candidate-a", "conflict_set_digest": "conflict-1"})

    fields, indexed = MODULE.read_index(path, ("candidate_digest", "consensus_digest"))

    assert fields == ["consensus_digest", "conflict_set_digest"]
    assert indexed["candidate-a"]["conflict_set_digest"] == "conflict-1"


def test_offset_fit_learns_correction_without_changing_unavailable_rows() -> None:
    labels = np.asarray([0, 0, 1, 1], dtype=float)
    offset = np.zeros(4)
    correction = np.asarray([[0.0], [-2.0], [0.0], [2.0]])

    beta = MODULE.fit_offset(offset, correction, labels)

    assert beta.shape == (1,)
    assert beta[0] > 0
    corrected = offset + correction @ beta
    assert np.allclose(offset[[0, 2]], corrected[[0, 2]])


def test_percentile_scores_average_ties() -> None:
    observed = MODULE.percentile_scores(np.asarray([3.0, 1.0, 3.0, 2.0]))
    assert np.allclose(observed, np.asarray([5 / 6, 0.0, 5 / 6, 1 / 3]))


def test_conflict_metric_uses_only_sets_with_exactly_one_positive() -> None:
    dataset = MODULE.Dataset(
        name="toy",
        role="development",
        paths=(Path("a"), Path("b"), Path("c")),
        digests=["a", "b", "c", "d", "e"],
        rows=[{} for _ in range(5)],
        labels=np.asarray([1, 0, 0, 1, 1]),
        groups=np.asarray(["1"] * 5),
        global_correction=np.empty((5, 6)),
        support_correction=np.empty((5, 19)),
        conflict_sets=np.asarray(["x", "x", "y", "y", "y"]),
        topology_available=np.zeros(5),
    )

    report = MODULE.conflict_metrics(dataset, np.asarray([2.0, 1.0, 0.0, 1.0, 2.0]))

    assert report["conflict_sets_total"] == 2
    assert report["evaluable_exactly_one_positive"] == 1
    assert report["top1_accuracy"] == 1.0


def test_grouped_oof_uses_distinct_chromosome_partitions(monkeypatch: pytest.MonkeyPatch) -> None:
    groups = np.repeat(np.asarray([f"chr{index}" for index in range(10)]), 4)
    labels = np.tile(np.asarray([0, 0, 1, 1], dtype=np.uint8), 10)
    dataset = MODULE.Dataset(
        name="toy",
        role="development",
        paths=(Path("a"), Path("b"), Path("c")),
        digests=[f"candidate-{index}" for index in range(len(labels))],
        rows=[{} for _ in range(len(labels))],
        labels=labels,
        groups=groups,
        global_correction=np.empty((len(labels), 6)),
        support_correction=np.empty((len(labels), 19)),
        conflict_sets=np.asarray([f"conflict-{index}" for index in range(len(labels))]),
        topology_available=np.zeros(len(labels)),
    )

    monkeypatch.setattr(
        MODULE,
        "fit_predict",
        lambda train, test: (
            {name: test.labels.astype(float) for name in MODULE.MODEL_NAMES},
            {},
        ),
    )

    _, reports = MODULE.grouped_oof(dataset)
    partitions = {
        tuple(sorted(tuple(model["test_groups"]) for model in report["fold_models"]))
        for report in reports
    }

    assert len(partitions) == len(MODULE.OOF_SEEDS)
    for report in reports:
        for model in report["fold_models"]:
            assert set(model["train_groups"]).isdisjoint(model["test_groups"])
