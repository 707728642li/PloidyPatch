from __future__ import annotations

import importlib.util
import csv
import hashlib
import json
import sys
from pathlib import Path

import pytest

np = pytest.importorskip("numpy")
pytest.importorskip("scipy")
pytest.importorskip("sklearn")


SCRIPT = Path(__file__).parents[1] / "scripts" / "evaluate_stable_reference_ranker_v0.9.py"
SPEC = importlib.util.spec_from_file_location("stable_reference_ranker_v0_9", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def _dataset(name: str, rows: int) -> object:
    labels = np.asarray(([0, 1] * (rows // 2))[:rows], dtype=np.uint8)
    return MODULE.Dataset(
        name=name,
        paths=(Path("a"), Path("b"), Path("c"), Path("d")),
        manifest_paths=(Path("e"), Path("f")),
        digests=tuple(f"{index:064x}" for index in range(rows)),
        copy_rows=tuple({} for _ in range(rows)),
        topology_rows=tuple({} for _ in range(rows)),
        labels=labels,
        groups=np.asarray([f"chr{index % 5}" for index in range(rows)]),
        conflict_sets=np.asarray(["" for _ in range(rows)]),
        correction=np.zeros((rows, len(MODULE.CORRECTION_FEATURE_NAMES))),
        topology_available=np.zeros(rows, dtype=np.uint8),
    )


def test_species_equal_weights_have_mean_one_and_equal_species_mass() -> None:
    first = _dataset("actinidia", 10)
    second = _dataset("populus", 30)

    weights = MODULE.species_equal_weights([first, second])

    assert weights.mean() == pytest.approx(1.0)
    assert weights[:10].sum() == pytest.approx(weights[10:].sum())


def test_weighted_offset_preserves_zero_correction_rows() -> None:
    labels = np.asarray([0, 0, 1, 1], dtype=float)
    offset = np.zeros(4)
    correction = np.asarray([[0.0], [-2.0], [0.0], [2.0]])
    weights = np.ones(4)

    coefficient = MODULE.fit_offset(offset, correction, labels, weights)

    assert coefficient[0] > 0
    corrected = offset + correction @ coefficient
    assert np.array_equal(corrected[[0, 2]], offset[[0, 2]])


def test_review_metrics_has_all_six_budgets_and_digest_tie_break() -> None:
    dataset = _dataset("actinidia", 600)
    scores = np.zeros(600)

    report = MODULE.review_metrics(dataset, scores)

    assert set(report) == {
        "top_0.5pct",
        "top_1pct",
        "top_2pct",
        "top_100",
        "top_250",
        "top_500",
    }
    expected = "".join(f"{index:064x}\n" for index in range(3))
    assert report["top_0.5pct"]["selection_digest_sha256"] == MODULE.sha256_text(expected)


def test_stable_constants_match_frozen_protocol() -> None:
    assert MODULE.FIXED_C == 1.0
    assert MODULE.FOLDS == 5
    assert len(MODULE.OOF_SEEDS) == 5
    assert len(set(MODULE.OOF_SEEDS)) == 5
    assert MODULE.BOOTSTRAP_REPLICATES == 20_000
    assert MODULE.FIXED_BUDGETS == (100, 250, 500)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _write_development_dataset(root: Path, name: str) -> tuple[Path, ...]:
    root.mkdir()
    copy_path = root / "copy.tsv"
    topology_path = root / "topology.tsv"
    labels_path = root / "labels.tsv"
    decisions_path = root / "decisions.tsv"
    digests = [hashlib.sha256(f"{name}:{index}".encode()).hexdigest() for index in range(60)]
    with copy_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=MODULE.FEATURE_COLUMNS, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for index, digest in enumerate(digests):
            row = {field: "0" for field in MODULE.FEATURE_COLUMNS}
            row.update(
                {
                    "candidate_digest": digest,
                    "seqid": f"chr{index % 10}",
                    "start": str(index * 100 + 1),
                    "end": str(index * 100 + 90),
                    "strand": "+",
                    "span_bp": "90",
                    "cds_segments": "2",
                    "cds_bp": "75",
                    "support_method_count": "1",
                    "support_methods": "miniprot",
                    "has_miniprot": "1",
                    "miniprot_identity": format(0.5 + 0.01 * (index % 5), ".2f"),
                    "miniprot_query_coverage": "0.8",
                }
            )
            writer.writerow(row)
    with topology_path.open("w", encoding="utf-8", newline="") as handle:
        fields = (
            "candidate_digest",
            "topology_available",
            "cds_bp_ratio",
            "cds_segment_count_ratio",
            "phase_lcs_similarity",
            "junction_fraction_similarity",
            "coding_span_ratio",
        )
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for index, digest in enumerate(digests):
            positive = index % 6 == 0
            value = "0.9" if positive else "0.2"
            writer.writerow(
                {
                    "candidate_digest": digest,
                    "topology_available": 1,
                    **{field: value for field in fields[2:]},
                }
            )
    with labels_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("candidate_digest", "label_exact_cds"),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for index, digest in enumerate(digests):
            writer.writerow(
                {"candidate_digest": digest, "label_exact_cds": int(index % 6 == 0)}
            )
    with decisions_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("candidate_digest", "status", "conflict_set_digest"),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for index, digest in enumerate(digests):
            writer.writerow(
                {
                    "candidate_digest": digest,
                    "status": "accepted",
                    "conflict_set_digest": f"conflict-{index // 2}",
                }
            )
    copy_manifest = {
        "truth_access": False,
        "outputs": {"features": {"sha256": _sha256(copy_path)}},
    }
    topology_manifest = {
        "truth_access": False,
        "inputs": {"copy_features": _sha256(copy_path)},
        "outputs": {"features": {"sha256": _sha256(topology_path)}},
    }
    Path(str(copy_path) + ".manifest.json").write_text(
        json.dumps(copy_manifest), encoding="utf-8"
    )
    Path(str(topology_path) + ".manifest.json").write_text(
        json.dumps(topology_manifest), encoding="utf-8"
    )
    return copy_path, topology_path, labels_path, decisions_path


def test_end_to_end_writes_exact_frozen_artifact_tree(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    actinidia = _write_development_dataset(tmp_path / "actinidia", "actinidia")
    populus = _write_development_dataset(tmp_path / "populus", "populus")
    protocol = tmp_path / "protocol.md"
    protocol.write_text("frozen before labels\n", encoding="utf-8")
    output = tmp_path / "evaluation"
    monkeypatch.setattr(MODULE, "BOOTSTRAP_REPLICATES", 64)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            str(MODULE.SCRIPT if hasattr(MODULE, "SCRIPT") else SCRIPT),
            "--dataset",
            "actinidia=" + ",".join(str(path) for path in actinidia),
            "--dataset",
            "populus=" + ",".join(str(path) for path in populus),
            "--protocol",
            str(protocol),
            "--output-dir",
            str(output),
        ],
    )

    assert MODULE.main() == 0

    assert not Path(str(output) + ".working").exists()
    expected = {
        "evaluation.json",
        "input_manifest.tsv",
        "partitions.json",
        "pooled_model.json",
        "predictions.tsv",
        "run_contract.tsv",
        "transfer_models.json",
    }
    manifest_paths = {
        line.split("  ", 1)[1]
        for line in (output / "SHA256SUMS").read_text(encoding="utf-8").splitlines()
    }
    assert manifest_paths == expected
    report = json.loads((output / "evaluation.json").read_text(encoding="utf-8"))
    assert set(report["oof"]) == {"actinidia", "populus"}
    assert set(report["cross_species_transfer"]) == {
        "actinidia_to_populus",
        "populus_to_actinidia",
    }
