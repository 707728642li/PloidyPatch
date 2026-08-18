from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def _namespace(relative: str) -> dict[str, object]:
    value: dict[str, object] = {"__name__": "coffea_pipeline_test"}
    path = ROOT / relative
    exec(compile(path.read_text(encoding="utf-8"), str(path), "exec"), value)
    return value


def test_target_gene_cds_uses_exact_protein_transcript_relation(tmp_path: Path) -> None:
    gff = tmp_path / "supported.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "chr1\tx\tgene\t1\t30\t.\t+\t.\tID=gene:g1\n"
        "chr1\tx\tmRNA\t1\t30\t.\t+\t.\tID=transcript:t1;Parent=gene:g1\n"
        "chr1\tx\tCDS\t1\t30\t.\t+\t0\tID=cds:p1;Parent=transcript:t1;protein_id=p1\n",
        encoding="utf-8",
    )
    representatives = tmp_path / "representatives.tsv"
    representatives.write_text(
        "gene_id\tprotein_id\trelation_source\tprotein_length\n"
        "g1\tp1\tFASTA_first_token->GFF_CDS:protein_id\t10\n",
        encoding="utf-8",
    )
    cds = tmp_path / "cds.fa"
    cds.write_text(">t1\nATGATGATG\n", encoding="utf-8")
    output = tmp_path / "gene.cds.fa"
    _namespace("scripts/prepare_coffea_evaluator_wgdi_inputs_v1.0.py")[
        "build_target_gene_cds"
    ](
        supported_gff=gff,
        representatives=representatives,
        transcript_cds=cds,
        output=output,
    )
    assert output.read_text(encoding="utf-8") == ">g1\nATGATGATG\n"


def test_coffea_evaluator_pipeline_is_no_ranker_and_exact_three_stream() -> None:
    wgdi = (ROOT / "scripts/run_coffea_evaluator_wgdi_v1.0.sh").read_text(
        encoding="utf-8"
    )
    pair = (ROOT / "scripts/infer_coffea_external_pairs_v1.0.py").read_text(
        encoding="utf-8"
    )
    holdout = (ROOT / "scripts/build_coffea_structure_holdout_v1.0.py").read_text(
        encoding="utf-8"
    )
    assert all(name in wgdi for name in ("car_self", "car_vs_gja", "car_vs_opu"))
    assert "--max-target-seqs 20" in wgdi
    assert "min_block_pairs=20" in pair
    assert "attach_descriptive_yn00_ks" in pair
    assert '"yn00_ks_selection_use": False' in pair
    assert 'balance_group_field="homoeolog_group"' in holdout
    assert "MINIMUM_CHROMOSOMES = 17" in holdout
    assert "MINIMUM_HOMOEOLOG_GROUPS = 9" in holdout
    combined = (wgdi + pair + holdout).casefold()
    assert "score_candidates" not in combined
    assert "topology_features" not in combined
    assert "composite_model" not in combined


def test_coffea_protocol_freeze_precedes_all_pair_and_candidate_outputs() -> None:
    protocol = json.loads(
        Path(
            ROOT / "config/holdouts/coffea_et39_v1.0/contract.json"
        ).read_text(encoding="utf-8")
    )
    assert protocol["truth_blind"]["wgd_pairs_enumerated_before_protocol_freeze"] is False
    assert protocol["truth_blind"]["candidate_counts_computed_before_protocol_freeze"] is False
    assert protocol["scientific_parameters"]["h2_or_topology_ranking"] == "forbidden"
