from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_coffea_amendment_is_label_blind_and_bounded() -> None:
    text = (ROOT / "docs/COFFEA_BLIND_REPRODUCIBILITY_AMENDMENT_v1.1.md").read_text(
        encoding="utf-8"
    )
    normalized = " ".join(text.split())
    assert "exact complete decision row" in text
    assert "without inspecting labels" in text
    assert "may not contain a species gene allowlist or denylist" in normalized
    assert "will not trigger another species" in normalized


def test_reconciliation_entry_has_no_gene_specific_exceptions() -> None:
    text = (ROOT / "scripts/reconcile_coffea_blind_replicates_v1.1.py").read_text(
        encoding="utf-8"
    )
    assert "compare_decision_tables" in text
    assert "filter_candidate_gff_by_upstream_models" in text
    assert "label_access\": False" in text
    for unstable_id in (
        "Coeug002g041590.1",
        "Chr06g0171471.1",
        "Chr08g0227331.1",
        "Chr11g0277251.1_1",
    ):
        assert unstable_id not in text
