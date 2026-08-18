from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from ploidypatch.event_graph import infer_event_graph


CANDIDATES = (
    "event_id\tevent_type\tcandidate_id\tgene_ids\trelationship_state\t"
    "wgd_event\tsubgenome\n"
    "E1\tmissing_exon\tC1\tG1\t\t\tA\n"
    "E2\ttrue_absence\tC2\tG2\t\t\tC\n"
    "E3\tmissing_annotation\tC3\tG3\twgd_homeolog\t\tA\n"
    "E4\tmissing_annotation\tC4a\tG4\t\t\tA\n"
    "E4\tmissing_annotation\tC4b\tG4\t\t\tA\n"
)

EVIDENCE_HEADER = (
    "candidate_id\tevidence_id\tevidence_type\tdirection\tscope\tstrength\t"
    "reliability\tindependent_group\tsource\tdetails\n"
)
EVIDENCE = EVIDENCE_HEADER + (
    "C1\tP1\tprotein_projection\tsupport\tboth\t1\t1\thomology\tminiprot\t\n"
    "C1\tP2\tprotein_projection\tsupport\tboth\t0.8\t1\thomology\tgemoma\t\n"
    "C1\tR1\trna_junction\tsupport\tisoform\t1\t1\trna\tSTAR\t\n"
    "C1\tS1\tsplice_site\tsupport\tisoform\t1\t1\tsplice\tplant_model\t\n"
    "C2\tA1\twhole_genome_alignment_absence\tsupport\tcoding\t1\t1\twga\tCactus\t\n"
    "C2\tA2\tprotein_projection\tcontradict\tcoding\t1\t1\thomology\tminiprot\t\n"
    "C3\tW1\tprotein_projection\tsupport\tboth\t1\t1\thomology\tminiprot\t\n"
    "C3\tW2\tsynteny\tsupport\trelationship\t1\t1\tsynteny\tWGDI\t\n"
    "C4a\tX1\tprotein_projection\tsupport\tboth\t1\t1\thomology_a\tminiprot\t\n"
    "C4b\tX2\tprotein_projection\tsupport\tboth\t1\t1\thomology_b\tminiprot\t\n"
)


def _write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8", newline="")
    return path


def test_event_graph_retains_contradictions_and_applies_constraints(
    tmp_path: Path,
) -> None:
    candidates = _write(tmp_path / "candidates.tsv", CANDIDATES)
    evidence = _write(tmp_path / "evidence.tsv", EVIDENCE)
    output_json = tmp_path / "graph.json"
    decisions_tsv = tmp_path / "decisions.tsv"

    graph = infer_event_graph(
        candidates,
        evidence,
        output_json_path=output_json,
        decisions_tsv_path=decisions_tsv,
    )
    rows = {row["candidate_id"]: row for row in graph["decisions"]}

    assert graph["counts"]["events"] == 4
    assert graph["counts"]["candidates"] == 5
    assert rows["C1"]["decision"] == "accept_high_confidence"
    assert rows["C1"]["selected"] is True
    assert rows["C1"]["coding_model"]["logit"] == 2.0
    assert rows["C1"]["coding_model"]["contributing_edges"] == 2
    assert rows["C1"]["coding_model"]["independent_groups"] == 1
    assert rows["C1"]["event_coherence"]["available"] is False
    assert rows["C1"]["risk_output"] == {
        "tier": "abstain",
        "reason": "paired_before_after_evidence_missing",
        "automatic_approval": False,
    }

    assert rows["C2"]["decision"] == "abstain"
    assert any(
        item["constraint"] == "absence_vs_intact_locus_evidence"
        for item in rows["C2"]["constraint_violations"]
    )
    assert {edge["direction"] for edge in rows["C2"]["evidence_edges"]} == {
        "support",
        "contradict",
    }

    assert rows["C3"]["decision"] == "abstain"
    assert any(
        item["constraint"] == "wgd_homeolog_requires_named_event"
        for item in rows["C3"]["constraint_violations"]
    )

    assert rows["C4a"]["decision"] == "review"
    assert rows["C4b"]["decision"] == "review"
    assert rows["C4a"]["decision_reason"] == "competing_high_confidence_alternatives"
    assert rows["C4b"]["selected"] is False
    assert graph["counts"]["selected"] == 1

    persisted = json.loads(output_json.read_text(encoding="utf-8"))
    assert persisted["inputs"] == graph["inputs"]
    with decisions_tsv.open("r", encoding="utf-8", newline="") as handle:
        written = list(csv.DictReader(handle, delimiter="\t"))
    assert len(written) == 5
    assert {row["candidate_id"] for row in written} == set(rows)


def test_event_graph_rejects_unknown_evidence_and_refuses_overwrite(
    tmp_path: Path,
) -> None:
    candidates = _write(tmp_path / "candidates.tsv", CANDIDATES)
    bad_evidence = _write(
        tmp_path / "bad.tsv",
        EVIDENCE_HEADER
        + "C1\tbad\tnot_a_type\tsupport\tboth\t1\t1\tg\ts\t\n",
    )
    with pytest.raises(ValueError, match="Unknown evidence_type"):
        infer_event_graph(candidates, bad_evidence)

    evidence = _write(tmp_path / "evidence.tsv", EVIDENCE)
    output = tmp_path / "graph.json"
    infer_event_graph(candidates, evidence, output_json_path=output)
    with pytest.raises(FileExistsError, match="Refusing to overwrite"):
        infer_event_graph(candidates, evidence, output_json_path=output)


def test_multilocus_plant_repair_requires_disambiguating_synteny(
    tmp_path: Path,
) -> None:
    candidates = _write(
        tmp_path / "candidates.tsv",
        "event_id\tevent_type\tcandidate_id\tgene_ids\n"
        "E1\tmissing_annotation\tC1\tG1\n"
        "E2\tmissing_annotation\tC2\tG2\n",
    )
    evidence = _write(
        tmp_path / "evidence.tsv",
        EVIDENCE_HEADER
        + "C1\tC1P\tprotein_projection\tsupport\tboth\t1\t1\tp1\tminiprot\t\n"
        "C1\tC1R\trna_junction\tsupport\tboth\t1\t1\tr1\tRNA\t\n"
        "C1\tC1M\tmappability\tcontradict\tcoding\t1\t1\tm1\tminiprot\t\n"
        "C2\tC2P\tprotein_projection\tsupport\tboth\t1\t1\tp2\tminiprot\t\n"
        "C2\tC2R\trna_junction\tsupport\tboth\t1\t1\tr2\tRNA\t\n"
        "C2\tC2M\tmappability\tcontradict\tcoding\t1\t1\tm2\tminiprot\t\n"
        "C2\tC2S\tsynteny\tsupport\trelationship\t1\t1\ts2\tWGDI\t\n",
    )
    graph = infer_event_graph(candidates, evidence)
    rows = {row["candidate_id"]: row for row in graph["decisions"]}

    assert rows["C1"]["decision"] == "abstain"
    assert any(
        violation["constraint"] == "multilocus_repair_requires_synteny"
        for violation in rows["C1"]["constraint_violations"]
    )
    assert rows["C2"]["decision"] == "accept_high_confidence"


def test_event_coherence_is_before_after_and_dependency_adjusted(
    tmp_path: Path,
) -> None:
    candidates = _write(
        tmp_path / "candidates.tsv",
        "event_id\tevent_type\tcandidate_id\tgene_ids\trelationship_state\t"
        "wgd_event\n"
        "E1\tmissing_annotation\tC1\tG1\twgd_homeolog\tWGD_A\n",
    )
    evidence = _write(
        tmp_path / "evidence.tsv",
        EVIDENCE_HEADER.rstrip("\n")
        + "\tevidence_phase\talgorithm_family\tbiological_source_group\n"
        + "C1\tB1\tprotein_projection\tcontradict\tboth\t1\t1\tb1\tminiprot\t\t"
        "before_edit\tprojection\trefA\n"
        + "C1\tB2\tcds_chain_concordance\tcontradict\tboth\t0.9\t1\tb2\tGeMoMa\t\t"
        "before_edit\tgemoma\trefA\n"
        + "C1\tA1\tprotein_projection\tsupport\tboth\t1\t1\ta1\tminiprot\t\t"
        "after_edit\tprojection\trefA\n"
        + "C1\tA2\tcds_chain_concordance\tsupport\tboth\t0.9\t1\ta2\tminiprot\t\t"
        "after_edit\tprojection\trefB\n"
        + "C1\tB3\tsynteny\tcontradict\trelationship\t1\t1\tb3\tWGDI\t\t"
        "before_edit\twgdi\tgenome_pair\n"
        + "C1\tA3\tsynteny\tsupport\trelationship\t1\t1\ta3\tWGDI\t\t"
        "after_edit\twgdi\tgenome_pair\n",
    )

    graph = infer_event_graph(candidates, evidence)
    row = graph["decisions"][0]
    coherence = row["event_coherence"]

    assert coherence["available"] is True
    assert coherence["scopes"]["coding"]["gain"] == pytest.approx(4.0)
    assert coherence["scopes"]["coding"]["raw_edge_gain"] == pytest.approx(7.6)
    assert coherence["scopes"]["coding"]["dependency_adjustment"] == pytest.approx(
        3.6
    )
    assert coherence["scopes"]["coding"]["before_edit"][
        "dependency_components"
    ] == 1
    assert coherence["scopes"]["coding"]["after_edit"][
        "dependency_components"
    ] == 1
    assert coherence["scopes"]["relationship"]["gain"] == pytest.approx(3.0)
    assert coherence["overall_gain"] == pytest.approx(3.5)
    assert row["risk_output"] == {
        "tier": "review_ranked",
        "reason": "positive_dependency_adjusted_coherence_gain",
        "automatic_approval": False,
    }
    assert row["decision"] == "accept_high_confidence"
    assert graph["counts"]["risk_review_ranked"] == 1
