from __future__ import annotations

from ploidypatch.homeolog_ranker import TOPOLOGY_COMPONENT_FIELDS
from ploidypatch.support_ranker import (
    CORRECTION_FEATURE_NAMES,
    SUPPORT_PATTERNS,
    support_conditioned_correction_vector,
)


def _topology(available: int) -> dict[str, str]:
    row = {"candidate_digest": "candidate-a", "topology_available": str(available)}
    for index, field in enumerate(TOPOLOGY_COMPONENT_FIELDS, start=1):
        row[field] = str(index / 10) if available else ""
    return row


def test_support_conditioned_correction_is_zero_without_topology() -> None:
    values = support_conditioned_correction_vector(
        {"support_methods": "gemoma,miniprot"}, _topology(0)
    )

    assert len(values) == len(CORRECTION_FEATURE_NAMES)
    assert values == [0.0] * len(values)


def test_support_conditioned_correction_uses_only_matching_pattern() -> None:
    pattern = "gemoma,miniprot"
    values = support_conditioned_correction_vector(
        {"support_methods": pattern}, _topology(1)
    )

    component_count = len(TOPOLOGY_COMPONENT_FIELDS)
    assert values[:component_count] == [0.1, 0.2, 0.3, 0.4, 0.5]
    pattern_values = values[component_count:]
    expected = []
    for level in SUPPORT_PATTERNS:
        expected.extend((float(level == pattern), 0.3 * float(level == pattern)))
    assert pattern_values == expected
