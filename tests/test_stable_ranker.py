from __future__ import annotations

import copy

import pytest

from ploidypatch.copy_features import FEATURE_COLUMNS
from ploidypatch.copy_model import fit_copy_feature_contract
from ploidypatch.stable_ranker import (
    _weighted_median,
    fit_weighted_copy_feature_contract,
)


def _row(**updates: str) -> dict[str, str]:
    row = {field: "0" for field in FEATURE_COLUMNS}
    row.update(
        {
            "candidate_digest": "a" * 64,
            "seqid": "chr1",
            "strand": "+",
            "support_methods": "miniprot",
        }
    )
    row.update(updates)
    return row


def test_weighted_median_reproduces_even_and_odd_unweighted_median() -> None:
    assert _weighted_median([4.0, 1.0, 3.0], [1.0, 1.0, 1.0]) == 3.0
    assert _weighted_median([4.0, 1.0, 3.0, 2.0], [1.0] * 4) == 2.5


def test_equal_weights_reproduce_existing_full_contract() -> None:
    rows = [
        _row(span_bp="10", miniprot_identity="", support_methods="miniprot"),
        _row(span_bp="20", miniprot_identity="0.5", support_methods="gemoma"),
        _row(span_bp="40", miniprot_identity="0.75", support_methods="lifton"),
        _row(span_bp="80", miniprot_identity="1", support_methods="gemoma,lifton"),
    ]

    old = fit_copy_feature_contract(rows, feature_set="full")
    new = fit_weighted_copy_feature_contract(
        rows, sample_weight=[1.0] * len(rows), feature_set="full"
    )

    comparable = copy.deepcopy(new)
    comparable["missing_value_policy"] = old["missing_value_policy"]
    comparable["scaling_policy"] = old["scaling_policy"]
    comparable.pop("weight_policy")
    assert comparable == old


def test_species_balancing_changes_imputation_and_centering() -> None:
    rows = [
        _row(span_bp="1", miniprot_identity="0.1"),
        _row(span_bp="2", miniprot_identity="0.2"),
        _row(span_bp="100", miniprot_identity="0.9"),
    ]
    # First two rows form one species (total mass one); final row is another.
    contract = fit_weighted_copy_feature_contract(
        rows, sample_weight=[0.5, 0.5, 1.0]
    )

    identity = contract["numeric"]["miniprot_identity"]
    assert identity["impute_median"] == pytest.approx((0.2 + 0.9) / 2)
    assert identity["mean"] == pytest.approx((0.1 * 0.5 + 0.2 * 0.5 + 0.9) / 2)


@pytest.mark.parametrize(
    "weights",
    ([1.0], [1.0, 0.0], [1.0, float("nan")], [1.0, float("inf")]),
)
def test_invalid_weights_fail_closed(weights: list[float]) -> None:
    with pytest.raises(ValueError):
        fit_weighted_copy_feature_contract([_row(), _row()], sample_weight=weights)

