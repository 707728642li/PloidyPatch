from __future__ import annotations

import csv
import hashlib
from pathlib import Path

from ploidypatch.calibration import (
    PROJECTION_RULE_GRID,
    RULE_GRID,
    calibrate_synteny_tiers,
)


def _seqid_for_split(validation: bool) -> str:
    for number in range(1000):
        seqid = f"chr{number}"
        held_out = (
            int(hashlib.sha256(seqid.encode("utf-8")).hexdigest()[:8], 16) % 10
            < 3
        )
        if held_out == validation:
            return seqid
    raise AssertionError("Could not construct deterministic split fixture")


def test_calibrates_rules_without_using_validation_for_selection(
    tmp_path: Path,
) -> None:
    fields = [field for field, _, _ in RULE_GRID]
    path = tmp_path / "labels.tsv"
    rows = []
    for index, (validation, label, strong) in enumerate(
        ((False, 1, True), (False, 0, False), (True, 1, True), (True, 0, False))
    ):
        values = {
            "model_rank": 1 if strong else 20,
            "model_identity": 0.99 if strong else 0.5,
            "model_query_coverage": 1.0 if strong else 0.5,
            "baseline_existing_cds_overlap_fraction": 0.0 if strong else 0.2,
            "supporting_block_count": 3 if strong else 1,
            "target_gap_genes": 1 if strong else 5,
            "locus_span_bp": 10_000 if strong else 500_000,
            "best_block_pairs": 500 if strong else 2,
            "best_block_pvalue": 0.001 if strong else 0.2,
        }
        rows.append(
            {
                "model_id": f"M{index}",
                "query_seqid": _seqid_for_split(validation),
                "label_exact_cds_chain": label,
                "eligible": 1,
                **values,
            }
        )
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=(
                "model_id",
                "query_seqid",
                "label_exact_cds_chain",
                "eligible",
                *fields,
            ),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)

    output = tmp_path / "calibration.json"
    report = calibrate_synteny_tiers(
        labeled_tsv_path=path,
        output_json_path=output,
        label_column="label_exact_cds_chain",
        high_precision_floor=0.5,
        high_precision_min_selected=1,
        eligible_column="eligible",
    )

    assert report["split"]["train_rows"] == 2
    assert report["split"]["validation_rows"] == 2
    assert report["input"]["eligible_column"] == "eligible"
    assert report["selected_policies"]["balanced_train_f1"]["train"]["f1"] == 1.0
    assert (
        report["selected_policies"]["balanced_train_f1"]["validation"]["f1"]
        == 1.0
    )
    assert output.is_file()


def test_calibrates_projection_support_feature_set(tmp_path: Path) -> None:
    fields = [field for field, _, _ in PROJECTION_RULE_GRID]
    path = tmp_path / "projection_labels.tsv"
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("model_id", "model_seqid", "label", *fields),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for index, (validation, label) in enumerate(
            ((False, 1), (False, 0), (True, 1), (True, 0))
        ):
            strong = bool(label)
            writer.writerow(
                {
                    "model_id": f"P{index}",
                    "model_seqid": _seqid_for_split(validation),
                    "label": label,
                    "model_identity": 0.99 if strong else 0.5,
                    "model_query_coverage": 1.0 if strong else 0.5,
                    "baseline_existing_cds_overlap_fraction": 0.0,
                    "support_model_count": 5 if strong else 1,
                    "support_query_count": 3 if strong else 1,
                    "support_source_count": 2 if strong else 1,
                    "support_min_identity": 0.95 if strong else 0.5,
                    "support_min_query_coverage": 0.95 if strong else 0.5,
                }
            )
    report = calibrate_synteny_tiers(
        labeled_tsv_path=path,
        output_json_path=tmp_path / "projection.json",
        label_column="label",
        feature_set="projection_support",
        high_precision_floor=0.5,
        high_precision_min_selected=1,
    )
    assert report["input"]["feature_set"] == "projection_support"
    assert report["split"]["column"] == "model_seqid"
    assert report["selected_policies"]["balanced_train_f1"]["validation"]["f1"] == 1.0
