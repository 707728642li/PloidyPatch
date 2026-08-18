from __future__ import annotations

import csv
import json
from pathlib import Path

from ploidypatch.consensus import _chain_digest
from ploidypatch.copy_features import (
    build_copy_candidate_features,
    label_copy_candidate_features,
)


def _write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8", newline="")
    return path


def test_build_and_label_copy_candidate_features(tmp_path: Path) -> None:
    chain = ("A01", "+", ((10, 30, "0"), (50, 80, "0")))
    digest = _chain_digest(chain)
    consensus = _write(
        tmp_path / "consensus.tsv",
        "consensus_digest\tseqid\tstart\tend\tstrand\tcds_segments\tcds_bp\t"
        "support_method_count\tsupport_methods\tupstream_models\tstatus\treason\t"
        "redundant_with_digest\n"
        f"{digest}\tA01\t10\t80\t+\t2\t52\t3\tgemoma,lifton,miniprot\t"
        "gemoma:GM1,lifton:LT1,miniprot:MP1\taccepted\tmethod_support_pass\t\n",
    )
    miniprot = _write(
        tmp_path / "miniprot.tsv",
        "model_id\tstatus\tidentity\tquery_coverage\trank\tscore\tpositive\t"
        "frameshifts\tstop_codons\n"
        "MP1\taccepted\t0.98\t1\t1\t900\t0.99\t0\t0\n",
    )
    gemoma = _write(
        tmp_path / "gemoma.tsv",
        "model_id\tstatus\tupstream_attributes\n"
        "GM1\taccepted\tID=GM1;pAA=0.91;iAA=0.87;score=90;bestScore=100;"
        "ce=2;rce=2;start=M;stop=*;nps=0;\n",
    )
    lifton = _write(
        tmp_path / "lifton.tsv",
        "model_id\tstatus\tupstream_attributes\n"
        "LT1\taccepted\tID=LT1;dna_identity=0.97;protein_identity=0.99;"
        "mutation=inframe_deletion;\n",
    )
    wgd = _write(
        tmp_path / "wgd.tsv",
        "consensus_digest\tsupport_block_count\tlongest_block_pairs\tstatus\n"
        f"{digest}\t2\t44\taccepted\n",
    )
    features = tmp_path / "features.tsv"
    manifest = build_copy_candidate_features(
        consensus_decisions_tsv_path=consensus,
        method_decision_inputs=(
            f"miniprot={miniprot}",
            f"gemoma={gemoma}",
            f"lifton={lifton}",
        ),
        wgd_selection_tsv_path=wgd,
        output_tsv_path=features,
    )

    assert manifest["truth_access"] is False
    assert manifest["counts"]["accepted_candidates"] == 1
    with features.open(encoding="utf-8", newline="") as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    assert row["candidate_digest"] == digest
    assert row["support_method_count"] == "3"
    assert row["miniprot_identity"] == "0.98"
    assert row["gemoma_score_ratio"] == "0.9"
    assert row["gemoma_start_complete"] == "1"
    assert row["lifton_inframe_indel"] == "1"
    assert row["wgd_existing_partner"] == "1"

    truth = {
        "events": [
            {
                "event_id": "E1",
                "event_type": "annotation_copy_collapse",
                "target": {
                    "gene_id": "G1",
                    "seqid": "A01",
                    "start": 10,
                    "end": 80,
                    "strand": "+",
                },
                "line_edits": [
                    {
                        "source_raw_line": (
                            "A01\tx\tmRNA\t10\t80\t.\t+\t.\tID=T1;Parent=G1\n"
                        )
                    },
                    {
                        "source_raw_line": (
                            "A01\tx\tCDS\t10\t30\t.\t+\t0\tParent=T1\n"
                        )
                    },
                    {
                        "source_raw_line": (
                            "A01\tx\tCDS\t50\t80\t.\t+\t0\tParent=T1\n"
                        )
                    },
                ],
            }
        ]
    }
    truth_path = tmp_path / "truth.json"
    truth_path.write_text(json.dumps(truth), encoding="utf-8")
    labels = tmp_path / "labels.tsv"
    label_manifest = label_copy_candidate_features(
        feature_tsv_path=features,
        hidden_truth_json_path=truth_path,
        output_tsv_path=labels,
    )

    assert label_manifest["evaluator_only"] is True
    assert label_manifest["counts"]["positive_exact_cds"] == 1
    with labels.open(encoding="utf-8", newline="") as handle:
        labeled = next(csv.DictReader(handle, delimiter="\t"))
    assert labeled["label_exact_cds"] == "1"
    assert labeled["truth_event_id"] == "E1"
    assert labeled["truth_locus_overlap_fraction_shorter"] == "1"
