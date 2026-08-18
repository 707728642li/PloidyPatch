from __future__ import annotations

import json
from pathlib import Path

import pytest

from ploidypatch.species_applicability import evaluate_species_applicability


ROOT = Path(__file__).parents[1]
POLICY = ROOT / "config" / "species_applicability_policy_v0.6.tsv"


def metrics() -> dict[str, object]:
    return {
        "schema_version": "ploidypatch.species_applicability_metrics.v0.6",
        "species_id": "synthetic_quality_target",
        "controlled_holdout": True,
        "assembly": {
            "primary_seqid_count_matches_declared_karyotype": True,
            "primary_assembly_fraction": 0.99,
            "primary_non_N_fraction": 0.999,
            "assembly_BUSCO_complete_fraction": 0.98,
            "read_backed_QV": None,
        },
        "annotation": {
            "primary_protein_coding_gene_fraction": 0.99,
            "exact_unique_GFF_protein_mapping_fraction": 0.995,
            "valid_coding_hierarchy_fraction": 1.0,
            "valid_representative_translation_fraction": 0.98,
            "protein_BUSCO_complete_fraction": 0.96,
            "fuzzy_identifier_repairs": 0,
        },
        "backbone": {
            "independent_WGD_source_count": 2,
            "minimum_block_pairs": 20,
            "cross_primary_seqid_only": True,
            "primary_gene_midpoint_backbone_coverage": 0.72,
            "primary_chromosome_cell_coverage_fraction": 0.80,
            "minimum_cells_per_covered_chromosome_observed": 25,
            "unique_partner_chromosome_fraction_among_covered_genes": 0.79,
            "input_permutation_deterministic": True,
            "reverse_duplicate_invariant": True,
            "built_without_candidates": True,
            "built_without_truth_or_labels": True,
            "used_perturbed_annotation_if_controlled": True,
        },
        "source_manifests": {"assembly": "a" * 64, "annotation": "b" * 64},
    }


def write_metrics(tmp_path: Path, value: dict[str, object]) -> Path:
    path = tmp_path / "metrics.json"
    path.write_text(json.dumps(value) + "\n", encoding="utf-8", newline="")
    return path


def evaluate(tmp_path: Path, value: dict[str, object]):
    return evaluate_species_applicability(
        metrics_path=write_metrics(tmp_path, value),
        policy_path=POLICY,
        output_path=tmp_path / "report.json",
    )


def test_full_topology_evaluable_without_optional_qv(tmp_path: Path) -> None:
    report = evaluate(tmp_path, metrics())
    assert report["applicability_state"] == "full_topology_evaluable"
    assert report["input_quality_pass"] is True
    assert report["backbone_pass"] is True
    assert report["gates"]["assembly"]["optional_read_backed_QV"] is None
    assert report["claim_boundary"]["candidate_statistics_used"] is False


def test_backbone_failure_keeps_chain_preservation_scope(tmp_path: Path) -> None:
    value = metrics()
    value["backbone"]["unique_partner_chromosome_fraction_among_covered_genes"] = 0.69  # type: ignore[index]
    report = evaluate(tmp_path, value)
    assert report["applicability_state"] == "chain_preservation_only"
    assert report["input_quality_pass"] is True
    assert report["backbone_pass"] is False


def test_assembly_failure_is_not_a_negative_biological_result(tmp_path: Path) -> None:
    value = metrics()
    value["assembly"]["primary_assembly_fraction"] = 0.80  # type: ignore[index]
    report = evaluate(tmp_path, value)
    assert report["applicability_state"] == "not_evaluable_input_quality"
    assert report["claim_boundary"]["input_quality_failure_is_negative_biological_result"] is False


def test_available_low_qv_fails_but_missing_qv_does_not(tmp_path: Path) -> None:
    value = metrics()
    value["assembly"]["read_backed_QV"] = 29.9  # type: ignore[index]
    report = evaluate(tmp_path, value)
    assert report["applicability_state"] == "not_evaluable_input_quality"
    assert report["gates"]["assembly"]["optional_read_backed_QV"] is False


def test_forbids_post_candidate_or_truth_metrics(tmp_path: Path) -> None:
    for key in ("candidate_count", "labels", "truth_pairs", "ranking_AP"):
        value = metrics()
        value[key] = 1
        path = tmp_path / key
        path.mkdir()
        with pytest.raises(ValueError, match="Forbidden post-candidate metric"):
            evaluate(path, value)


def test_refuses_overwrite_and_symlinked_inputs(tmp_path: Path) -> None:
    value = metrics()
    evaluate(tmp_path, value)
    with pytest.raises(FileExistsError, match="overwrite"):
        evaluate_species_applicability(
            metrics_path=tmp_path / "metrics.json",
            policy_path=POLICY,
            output_path=tmp_path / "report.json",
        )

    target = tmp_path / "real.json"
    target.write_text(json.dumps(value), encoding="utf-8")
    link = tmp_path / "linked.json"
    try:
        link.symlink_to(target)
    except OSError:
        pytest.skip("Host does not permit symlink creation")
    with pytest.raises(ValueError, match="symlinked"):
        evaluate_species_applicability(
            metrics_path=link,
            policy_path=POLICY,
            output_path=tmp_path / "other.json",
        )
