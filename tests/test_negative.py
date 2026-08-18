from __future__ import annotations

import json
from pathlib import Path

from ploidypatch.io import iter_fasta, read_fai
from ploidypatch.negative import (
    MASK_TRUTH_SCHEMA_VERSION,
    audit_masked_gap_genome,
    create_masked_gap_control,
    score_masked_gap_abstention,
    summarize_masked_gap_selection,
)
from ploidypatch.perturb import generate_missing_gene_benchmark


GFF = (
    "##gff-version 3\n"
    "chr1\ttest\tgene\t1\t90\t.\t+\t.\tID=gene:G1\n"
    "chr1\ttest\tmRNA\t1\t90\t.\t+\t.\tID=transcript:T1;Parent=gene:G1\n"
    "chr1\ttest\texon\t1\t30\t.\t+\t.\tParent=transcript:T1\n"
    "chr1\ttest\tCDS\t1\t30\t.\t+\t0\tParent=transcript:T1\n"
    "chr1\ttest\tgene\t101\t190\t.\t-\t.\tID=gene:G2\n"
    "chr1\ttest\tmRNA\t101\t190\t.\t-\t.\tID=transcript:T2;Parent=gene:G2\n"
    "chr1\ttest\texon\t161\t190\t.\t-\t.\tParent=transcript:T2\n"
    "chr1\ttest\tCDS\t161\t190\t.\t-\t0\tParent=transcript:T2\n"
)


def _write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8", newline="")
    return path


def test_masked_gap_control_is_coordinate_stable_and_scores_abstention(
    tmp_path: Path,
) -> None:
    source = _write(tmp_path / "source.gff3", GFF)
    genome = _write(tmp_path / "genome.fa", ">chr1\n" + "A" * 300 + "\n")
    run = tmp_path / "run"
    generate_missing_gene_benchmark(source, run, count=1, seed=7)
    masked = tmp_path / "masked.fa"
    mask_truth = tmp_path / "mask_truth.json"

    manifest = create_masked_gap_control(
        source_genome_path=genome,
        hidden_truth_path=run / "hidden_truth.json",
        output_genome_path=masked,
        output_mask_truth_path=mask_truth,
    )

    event = json.loads(mask_truth.read_text(encoding="utf-8"))["mask"][
        "event_intervals"
    ][0]
    _, _, sequence = next(iter_fasta(masked))
    assert len(sequence) == 300
    assert set(sequence[event["start"] - 1 : event["end"]]) == {"N"}
    assert read_fai(str(masked) + ".fai") == {"chr1": 300}
    assert manifest["mask"]["events"] == 1
    audit = audit_masked_gap_genome(
        source_genome_path=genome,
        masked_genome_path=masked,
        mask_truth_path=mask_truth,
    )
    assert audit["quality_gate"]["grade"] == "pass"
    assert audit["mask"]["changed_bp"] == event["end"] - event["start"] + 1

    noop = score_masked_gap_abstention(
        perturbed_gff_path=run / "perturbed.gff3",
        candidate_gff_path=run / "perturbed.gff3",
        mask_truth_path=mask_truth,
    )
    oracle = score_masked_gap_abstention(
        perturbed_gff_path=run / "perturbed.gff3",
        candidate_gff_path=source,
        mask_truth_path=mask_truth,
    )
    assert noop["abstention"]["event_abstention_rate"] == 1.0
    assert oracle["abstention"]["event_false_repair_rate"] == 1.0
    assert (
        oracle["abstention"][
            "novel_structures_with_feature_claim_in_masked_loci"
        ]
        == 1
    )


def test_transcript_span_without_exonic_claim_is_not_a_false_repair(
    tmp_path: Path,
) -> None:
    perturbed = _write(tmp_path / "perturbed.gff3", GFF)
    candidate = _write(
        tmp_path / "candidate.gff3",
        GFF
        + "chr1\ttest\tgene\t31\t60\t.\t+\t.\tID=gene:bridge\n"
        + "chr1\ttest\tmRNA\t31\t60\t.\t+\t.\t"
        + "ID=transcript:bridge;Parent=gene:bridge\n"
        + "chr1\ttest\texon\t31\t40\t.\t+\t.\tParent=transcript:bridge\n"
        + "chr1\ttest\tCDS\t31\t40\t.\t+\t0\tParent=transcript:bridge\n"
        + "chr1\ttest\texon\t51\t60\t.\t+\t.\tParent=transcript:bridge\n"
        + "chr1\ttest\tCDS\t51\t60\t.\t+\t2\tParent=transcript:bridge\n",
    )
    mask_truth = tmp_path / "mask_truth.json"
    mask_truth.write_text(
        json.dumps(
            {
                "schema_version": MASK_TRUTH_SCHEMA_VERSION,
                "mask": {
                    "event_intervals": [
                        {
                            "event_id": "event-1",
                            "seqid": "chr1",
                            "start": 41,
                            "end": 50,
                        }
                    ]
                },
            }
        ),
        encoding="utf-8",
    )

    report = score_masked_gap_abstention(
        perturbed_gff_path=perturbed,
        candidate_gff_path=candidate,
        mask_truth_path=mask_truth,
        include_event_details=True,
    )

    assert report["abstention"]["events_with_false_repair"] == 0
    assert report["abstention"]["event_abstention_rate"] == 1.0
    assert (
        report["abstention"][
            "events_with_transcript_span_crossing_masked_locus"
        ]
        == 1
    )
    assert report["abstention"]["novel_structures_spanning_masked_loci"] == 1
    assert report["control_diagnostics"][
        "events_with_preexisting_feature_overlap"
    ] == 0
    assert report["event_details"][0]["candidate_structures_spanning"] == 1
    assert report["event_details"][0]["candidate_structures_claiming_features"] == 0


def test_masked_gap_control_excludes_background_feature_collisions(
    tmp_path: Path,
) -> None:
    source = _write(tmp_path / "source.gff3", GFF)
    genome = _write(tmp_path / "genome.fa", ">chr1\n" + "A" * 300 + "\n")
    run = tmp_path / "run"
    generate_missing_gene_benchmark(source, run, count=2, seed=7)
    background = _write(
        tmp_path / "background.gff3",
        (run / "perturbed.gff3").read_text(encoding="utf-8")
        + "chr1\ttest\tgene\t10\t20\t.\t+\t.\tID=gene:nested\n"
        + "chr1\ttest\tmRNA\t10\t20\t.\t+\t.\t"
        + "ID=transcript:nested;Parent=gene:nested\n"
        + "chr1\ttest\texon\t10\t20\t.\t+\t.\tParent=transcript:nested\n"
        + "chr1\ttest\tCDS\t10\t20\t.\t+\t0\tParent=transcript:nested\n",
    )

    manifest = create_masked_gap_control(
        source_genome_path=genome,
        hidden_truth_path=run / "hidden_truth.json",
        output_genome_path=tmp_path / "masked.fa",
        output_mask_truth_path=tmp_path / "mask_truth.json",
        background_gff_path=background,
    )
    truth = json.loads((tmp_path / "mask_truth.json").read_text(encoding="utf-8"))

    assert manifest["mask"]["requested_events"] == 2
    assert manifest["mask"]["events"] == 1
    assert manifest["mask"]["excluded_background_overlap_events"] == 1
    assert truth["selection"]["eligible_events"] == 1
    assert truth["selection"]["excluded_events"] == 1
    assert truth["source"]["background_gff"]["file_name"] == "background.gff3"

    hidden = json.loads((run / "hidden_truth.json").read_text(encoding="utf-8"))
    genes = [event["target"]["gene_id"] for event in hidden["events"]]
    strata = _write(
        tmp_path / "strata.tsv",
        "gene_id\tgroup\n" + "\n".join(f"{gene}\tgroup_{index}" for index, gene in enumerate(genes)) + "\n",
    )
    summary = summarize_masked_gap_selection(
        mask_truth_path=tmp_path / "mask_truth.json",
        hidden_truth_path=run / "hidden_truth.json",
        strata_tsv_path=strata,
        columns=("group",),
    )
    assert summary["totals"] == {
        "requested_events": 2,
        "eligible_events": 1,
        "excluded_events": 1,
    }
