from __future__ import annotations

import json
from pathlib import Path

from ploidypatch.baseline import _file_sha256
from ploidypatch.perturb import TRUTH_SCHEMA_VERSION, read_gff_document
from ploidypatch.structure_hypothesis import infer_structure_hypotheses
from ploidypatch.structure_hypothesis_score import score_structure_hypotheses


def write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8", newline="")
    return path


def test_score_structure_hypothesis_recovers_audit_candidate(
    tmp_path: Path,
) -> None:
    source_text = (
        "##gff-version 3\n"
        "chr1\ttest\tgene\t1\t90\t.\t+\t.\tID=G1\n"
        "chr1\ttest\tmRNA\t1\t90\t.\t+\t.\tID=T1;Parent=G1\n"
        "chr1\ttest\tCDS\t1\t20\t.\t+\t0\tParent=T1\n"
        "chr1\ttest\tCDS\t40\t50\t.\t+\t0\tParent=T1\n"
        "chr1\ttest\tCDS\t70\t90\t.\t+\t0\tParent=T1\n"
    )
    perturbed_text = source_text.replace(
        "chr1\ttest\tCDS\t40\t50\t.\t+\t0\tParent=T1\n", ""
    )
    candidate_text = (
        perturbed_text
        + "###\n"
        "chr1\tPloidyPatchCandidate\tgene\t1\t90\t.\t+\t.\t"
        "ID=CG;candidate_group=PPSC-truth;"
        "candidate_category=cds_extension_or_missing_segment;"
        "support_source_count=2;support_sources=a,b\n"
        "chr1\tPloidyPatchCandidate\tmRNA\t1\t90\t.\t+\t.\t"
        "ID=CT;Parent=CG;candidate_group=PPSC-truth;"
        "candidate_category=cds_extension_or_missing_segment;"
        "support_source_count=2;support_sources=a,b\n"
        "chr1\tPloidyPatchCandidate\tCDS\t1\t20\t.\t+\t0\tParent=CT\n"
        "chr1\tPloidyPatchCandidate\tCDS\t40\t50\t.\t+\t0\tParent=CT\n"
        "chr1\tPloidyPatchCandidate\tCDS\t70\t90\t.\t+\t0\tParent=CT\n"
    )
    source = write(tmp_path / "source.gff3", source_text)
    perturbed = write(tmp_path / "perturbed.gff3", perturbed_text)
    candidate = write(tmp_path / "candidate.gff3", candidate_text)
    hypotheses = tmp_path / "hypotheses.tsv"
    infer_structure_hypotheses(
        annotation_gff_path=perturbed,
        candidate_gff_path=candidate,
        output_tsv_path=hypotheses,
        candidate_topology_tsv_path=tmp_path / "topology.tsv",
    )
    truth = {
        "schema_version": TRUTH_SCHEMA_VERSION,
        "source": {
            "file_sha256": _file_sha256(source),
            "text_sha256": read_gff_document(source).text_sha256,
        },
        "perturbation": {
            "perturbed_text_sha256": read_gff_document(perturbed).text_sha256
        },
        "events": [
            {
                "event_id": "EVENT-1",
                "event_type": "annotation_missing_internal_exon",
                "target": {
                    "gene_ids": ["G1"],
                    "perturbed_gene_ids": ["G1"],
                    "transcript_ids": ["T1"],
                },
            }
        ],
    }
    truth_path = tmp_path / "truth.json"
    truth_path.write_text(json.dumps(truth), encoding="utf-8")

    report = score_structure_hypotheses(
        source_gff_path=source,
        perturbed_gff_path=perturbed,
        candidate_gff_path=candidate,
        hypotheses_tsv_path=hypotheses,
        truth_path=truth_path,
        include_event_details=True,
    )

    assert report["quality_gate"]["grade"] == "pass"
    assert report["raw_proposals"]["precision"] == 1.0
    assert report["event_recovery"]["recovered_events"] == 1
    assert report["event_recovery"]["cds_chain_recall"] == 1.0
    assert report["event_details"][0]["recovered"] is True
