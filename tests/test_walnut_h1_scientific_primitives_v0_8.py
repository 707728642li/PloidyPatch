from __future__ import annotations

import csv
import json
from pathlib import Path
import pytest

from ploidypatch.io import iter_fasta

from ploidypatch.walnut_h1 import (
    BOOTSTRAP_REPLICATES,
    BOOTSTRAP_SEED,
    EVALUATION_SCHEMA,
    PAIR_SCHEMA,
    RAW_SCHEMA,
    build_raw_predictions_manifest,
    evaluate_h1_scores,
    filter_self_pairs_by_closed_yn00_ks,
    infer_exact_two_outgroup_pair_consistent_truth,
    score_collateral_gate,
    seal_blind_pool_manifest,
    write_primary_candidate_proteins,
)


def test_core_h1_collateral_gate_is_lightweight_and_fail_closed() -> None:
    passing = {
        "quality_gate": {"grade": "pass"},
        "collateral_changes": {
            "baseline_transcript_structures_missing_from_candidate": 0
        },
    }
    assert score_collateral_gate(passing)
    assert not score_collateral_gate({})
    failing = json.loads(json.dumps(passing))
    failing["collateral_changes"][
        "baseline_transcript_structures_missing_from_candidate"
    ] = 1
    assert not score_collateral_gate(failing)


def write_tsv(path: Path, fields: tuple[str, ...], rows: list[tuple[object, ...]]) -> Path:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(fields)
        writer.writerows(rows)
    return path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_closed_yn00_filter_keeps_both_endpoints_and_abstains_missing(
    tmp_path: Path,
) -> None:
    self_pairs = write_tsv(
        tmp_path / "self.tsv",
        ("gene_id_a", "gene_id_b"),
        [("gene:A", "gene:B"), ("C", "D"), ("E", "F"), ("G", "H")],
    )
    ks = write_tsv(
        tmp_path / "ks.tsv",
        ("id1", "id2", "ks_YN00"),
        [("B", "A", "0.10"), ("C", "D", "0.75"), ("E", "F", "0.7501")],
    )
    output = tmp_path / "accepted.tsv"
    decisions = tmp_path / "decisions.tsv"
    manifest = filter_self_pairs_by_closed_yn00_ks(
        self_pairs_path=self_pairs,
        ks_path=ks,
        output_pairs_path=output,
        decisions_path=decisions,
    )

    assert [(row["gene_id_a"], row["gene_id_b"]) for row in read_tsv(output)] == [
        ("gene:A", "gene:B"),
        ("C", "D"),
    ]
    assert manifest["parameters"] == {
        "interval": "closed",
        "minimum": 0.1,
        "maximum": 0.75,
        "missing_or_nonfinite": "abstain",
        "thresholds_fitted_from_target": False,
        "pair_identifier_join": "exact_normalize_feature_id_then_unordered_pair",
    }
    reasons = {row["reason"] for row in read_tsv(decisions)}
    assert "missing_or_nonfinite_yn00_ks" in reasons
    assert "yn00_ks_outside_closed_0.10_0.75" in reasons


def test_candidate_proteins_are_exact_primary_only_header_aware_and_unique(
    tmp_path: Path,
) -> None:
    gff = tmp_path / "primary.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "chr1\tx\tgene\t1\t90\t.\t+\t.\tID=gene:G1\n"
        "chr1\tx\tmRNA\t1\t90\t.\t+\t.\tID=transcript:T1;Parent=gene:G1\n"
        "chr1\tx\tCDS\t1\t90\t.\t+\t0\tParent=transcript:T1;protein_id=P1\n"
        "chr1\tx\tgene\t101\t190\t.\t+\t.\tID=gene:G2\n"
        "chr1\tx\tmRNA\t101\t190\t.\t+\t.\tID=transcript:T2;Parent=gene:G2\n"
        "chr1\tx\tCDS\t101\t190\t.\t+\t0\tParent=transcript:T2\n",
        encoding="utf-8",
    )
    provider = tmp_path / "provider.fa"
    provider.write_text(
        ">P1 provider protein\nMAAAAA\n"
        ">ALT2 transcript:T2\nMCCCCCC\n"
        ">PLASTID1 transcript:plastid_t1\nMGGGGG\n",
        encoding="utf-8",
    )
    output = tmp_path / "primary.fa"
    whitelist = tmp_path / "whitelist.tsv"
    manifest_path = tmp_path / "manifest.json"
    manifest = write_primary_candidate_proteins(
        primary_gff_path=gff,
        provider_protein_path=provider,
        output_fasta_path=output,
        whitelist_tsv_path=whitelist,
        manifest_path=manifest_path,
    )
    assert [row["gene_id"] for row in read_tsv(whitelist)] == ["G1", "G2"]
    assert [record[0] for record in iter_fasta(output)] == ["P1", "ALT2"]
    assert manifest["counts"]["provider_records"] == 3
    assert manifest["counts"]["accepted_unique_representatives"] == 2
    assert manifest["counts"]["excluded_provider_records"] == 1

    duplicate = tmp_path / "duplicate.fa"
    duplicate.write_text(">P1\nMAAA\n>P1\nMCCC\n", encoding="utf-8")
    with pytest.raises(ValueError, match="Duplicate provider protein"):
        write_primary_candidate_proteins(
            primary_gff_path=gff,
            provider_protein_path=duplicate,
            output_fasta_path=tmp_path / "duplicate.out.fa",
            whitelist_tsv_path=tmp_path / "duplicate.tsv",
            manifest_path=tmp_path / "duplicate.json",
        )


def test_two_serial_wgd_alternate_pair_vetoes_event_even_in_one_group(
    tmp_path: Path,
) -> None:
    self_pairs = write_tsv(
        tmp_path / "self.tsv",
        ("gene_id_a", "gene_id_b", "yn00_ks"),
        [("A", "B", "0.2"), ("C", "D", "0.3")],
    )
    corylus = write_tsv(
        tmp_path / "corylus.tsv",
        ("gene_id_a", "gene_id_b", "status", "reason"),
        [
            ("A", "B", "accepted", "exact_1_to_2"),
            ("B", "A", "accepted", "second_evaluator_counterpart_same_target_pair"),
            ("A", "X", "rejected", "nonreciprocal_multiple_partners"),
            ("C", "D", "accepted", "exact_1_to_2"),
        ],
    )
    castanea = write_tsv(
        tmp_path / "castanea.tsv",
        ("gene_id_a", "gene_id_b", "status", "reason"),
        [("A", "B", "accepted", "exact_1_to_2"), ("C", "D", "accepted", "exact_1_to_2")],
    )
    output = tmp_path / "truth"
    manifest = infer_exact_two_outgroup_pair_consistent_truth(
        ks_filtered_self_pairs_path=self_pairs,
        evaluator_group_decisions={"corylus": corylus, "castanea": castanea},
        output_dir=output,
    )

    assert manifest["schema_version"] == PAIR_SCHEMA
    assert [(row["gene_id_a"], row["gene_id_b"]) for row in read_tsv(output / "pairs.tsv")] == [
        ("C", "D")
    ]
    rejected = next(row for row in read_tsv(output / "decisions.tsv") if row["gene_id_a"] == "A")
    assert rejected["reason"] == "missing_or_discordant_group:corylus"
    corylus_audit = manifest["counts"]["evaluator_groups"][0]
    assert corylus_audit["discordant_target_members"] == 1


def test_raw_prediction_manifest_and_both_pool_arms_are_custody_bound(
    tmp_path: Path,
) -> None:
    trees: dict[str, Path] = {}
    for method in ("miniprot", "gemoma", "lifton"):
        for bundle in ("candidate_mandshurica", "candidate_carya"):
            label = f"{method}__{bundle}"
            tree = tmp_path / "raw" / label
            tree.mkdir(parents=True)
            (tree / "prediction.gff3").write_text(
                f"##gff-version 3\n# {label}\n", encoding="utf-8"
            )
            trees[label] = tree
    raw_path = tmp_path / "raw_predictions.manifest.json"
    raw = build_raw_predictions_manifest(
        project_root=tmp_path,
        raw_prediction_trees=trees,
        input_hashes={"staged_inputs": "a" * 64, "blind_benchmark": "b" * 64},
        output=raw_path,
    )
    assert raw["schema_version"] == RAW_SCHEMA
    assert len(raw["raw_prediction_trees"]) == 6
    assert all(row["file_count"] == 1 for row in raw["raw_prediction_trees"].values())

    for arm in (
        "retain_distinct_phased_CDS_chains",
        "suppress_strongly_overlapping_alternative_chains",
    ):
        root = tmp_path / arm
        root.mkdir()
        candidate = root / "candidate.gff3"
        decisions = root / "decisions.tsv"
        candidate.write_text("##gff-version 3\n", encoding="utf-8")
        decisions.write_text("status\n", encoding="utf-8")
        import hashlib

        sha = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
        manifest_path = root / "candidate.gff3.manifest.json"
        schema = (
            "ploidypatch.method_candidate_pool.v2"
            if arm == "retain_distinct_phased_CDS_chains"
            else "ploidypatch.method_consensus.v1"
        )
        manifest_path.write_text(
            json.dumps(
                {
                    "schema_version": schema,
                    "inputs": {},
                    "outputs": {
                        "candidate_gff": {"sha256": sha(candidate)},
                        "decisions": {"sha256": sha(decisions)},
                    },
                }
            )
            + "\n",
            encoding="utf-8",
        )
        sealed = seal_blind_pool_manifest(
            manifest_path=manifest_path,
            candidate_gff_path=candidate,
            decisions_path=decisions,
            raw_predictions_manifest_path=raw_path,
            policy_arm=arm,
        )
        assert sealed["policy_arm"] == arm
        assert sealed["schema_version"] == schema
        assert sealed["truth_access"] is False
        assert sealed["ranker_access"] is False
        assert sealed["inputs"]["raw_predictions_manifest"]["file_name"] == (
            "raw_predictions.manifest.json"
        )


def score_payload(values: list[bool]) -> dict[str, object]:
    successes = sum(values)
    return {
        "schema_version": "ploidypatch.annotation_repair_score.v5",
        "quality_gate": {"grade": "pass"},
        "collateral_changes": {
            "baseline_transcript_structures_missing_from_candidate": 0
        },
        "event_recovery": {
            "events": len(values),
            "complete_cds_chain_recovery": successes,
            "complete_cds_chain_recall": successes / len(values),
        },
        "strict_cds_chain": {},
        "background_subtraction": {},
        "event_details": [
            {
                "event_id": f"event_{index}",
                "event_type": "annotation_copy_collapse",
                "complete_cds_chain_recovery": value,
            }
            for index, value in enumerate(values)
        ],
    }


def test_h1_evaluator_reports_only_paired_recall_and_zero_collateral(
    tmp_path: Path,
) -> None:
    retain = tmp_path / "retain.json"
    suppress = tmp_path / "suppress.json"
    retain.write_text(json.dumps(score_payload([True, True, True, False])) + "\n")
    suppress.write_text(json.dumps(score_payload([True, False, False, False])) + "\n")
    output = tmp_path / "evaluation.json"
    result = evaluate_h1_scores(
        retain_score_path=retain, suppress_score_path=suppress, output=output
    )

    assert result["schema_version"] == EVALUATION_SCHEMA
    assert result["paired_event_bootstrap"]["observed_delta"] == 0.5
    assert result["bootstrap_parameters"]["replicates"] == BOOTSTRAP_REPLICATES
    assert result["bootstrap_parameters"]["seed"] == BOOTSTRAP_SEED
    assert result["sentinels"] == {
        "all_arm_collateral_loss_zero": True,
        "automatic_approval": False,
    }
    text = output.read_text(encoding="utf-8").lower()
    assert result["protocol_profile"] == "core_H1_only_no_ranker"
    assert result["ranker_or_model_executed"] is False
    assert result["h2_or_topology_ranking_executed"] is False
    assert result["all_arm_collateral_loss"] == 0
    for forbidden in (
        "average_precision",
        "roc_auc",
        "topology_features",
        '"ranker":',
        '"scores":',
    ):
        assert forbidden not in text

    with pytest.raises(FileExistsError, match="overwrite"):
        evaluate_h1_scores(
            retain_score_path=retain, suppress_score_path=suppress, output=output
        )
