from __future__ import annotations

import csv
import json
from pathlib import Path

from ploidypatch.cli import build_parser
from ploidypatch.event_graph import infer_event_graph
from ploidypatch.natural import (
    ALTERNATIVE_STRUCTURE_EVENT,
    AMBIGUOUS_OPPOSITE_STRAND_EVENT,
    annotate_natural_assembly_context,
    discover_natural_candidates,
    prepare_natural_graph_inputs,
    summarize_natural_validation,
    validate_natural_candidates_with_rna,
    validate_natural_candidates_with_secondary_groups,
)
from ploidypatch.perturb import (
    BOUNDARY_SHIFT_EVENT,
    MISSING_GENE_EVENT,
    MISSING_INTERNAL_EXON_EVENT,
    SPLIT_GENE_EVENT,
    _file_sha256,
)
from ploidypatch.rna import (
    JUNCTION_AGGREGATE_SCHEMA,
    JUNCTION_GROUP_AGGREGATE_SCHEMA,
)


def write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8", newline="")
    return path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def make_projection_inputs(tmp_path: Path) -> tuple[Path, Path, Path]:
    target = write(
        tmp_path / "target.gff3",
        "##gff-version 3\n"
        "chr1\ttest\tgene\t100\t250\t.\t+\t.\tID=G1\n"
        "chr1\ttest\tmRNA\t100\t250\t.\t+\t.\tID=T1;Parent=G1\n"
        "chr1\ttest\tCDS\t100\t130\t.\t+\t0\tParent=T1\n"
        "chr1\ttest\tCDS\t200\t250\t.\t+\t0\tParent=T1\n"
        "chr1\ttest\tgene\t400\t450\t.\t+\t.\tID=G2\n"
        "chr1\ttest\tmRNA\t400\t450\t.\t+\t.\tID=T2;Parent=G2\n"
        "chr1\ttest\tCDS\t400\t450\t.\t+\t0\tParent=T2\n"
        "chr1\ttest\tgene\t600\t650\t.\t+\t.\tID=G3\n"
        "chr1\ttest\tmRNA\t600\t650\t.\t+\t.\tID=T3;Parent=G3\n"
        "chr1\ttest\tCDS\t600\t650\t.\t+\t0\tParent=T3\n"
        "chr1\ttest\tgene\t700\t750\t.\t+\t.\tID=G4\n"
        "chr1\ttest\tmRNA\t700\t750\t.\t+\t.\tID=T4;Parent=G4\n"
        "chr1\ttest\tCDS\t700\t750\t.\t+\t0\tParent=T4\n"
        "chr1\ttest\tgene\t1200\t1250\t.\t-\t.\tID=G5\n"
        "chr1\ttest\tmRNA\t1200\t1250\t.\t-\t.\tID=T5;Parent=G5\n"
        "chr1\ttest\tCDS\t1200\t1250\t.\t-\t0\tParent=T5\n",
    )
    query_ids = [
        "hap1__exact",
        "hap1__internal",
        "hap2__internal",
        "hap1__boundary",
        "hap1__split",
        "hap1__missing",
        "hap2__missing",
        "hap1__opposite",
        "hap1__weak",
    ]
    protein_map = write(
        tmp_path / "proteins.tsv",
        "query_id\tsource\tsource_record_id\tlength_aa\tsource_header\n"
        + "".join(
            f"{query_id}\t{query_id.split('__')[0]}\t{query_id}\t100\t{query_id}\n"
            for query_id in query_ids
        ),
    )

    def model(
        model_id: str,
        query_id: str,
        start: int,
        end: int,
        segments: list[tuple[int, int]],
        identity: float = 0.95,
        rank: int = 1,
        score: int = 100,
    ) -> str:
        lines = [
            f"chr1\tminiprot\tmRNA\t{start}\t{end}\t{score}\t+\t.\t"
            f"ID={model_id};Rank={rank};Identity={identity};Positive={identity};"
            f"Target={query_id} 1 100\n"
        ]
        lines.extend(
            f"chr1\tminiprot\tCDS\t{left}\t{right}\t100\t+\t0\t"
            f"Parent={model_id}\n"
            for left, right in segments
        )
        return "".join(lines)

    miniprot = write(
        tmp_path / "miniprot.gff3",
        "##gff-version 3\n"
        + model("MP_exact", "hap1__exact", 100, 250, [(100, 130), (200, 250)])
        + model(
            "MP_internal_1",
            "hap1__internal",
            100,
            250,
            [(100, 130), (160, 170), (200, 250)],
        )
        + model(
            "MP_internal_2",
            "hap2__internal",
            100,
            250,
            [(100, 130), (160, 170), (200, 250)],
            0.94,
        )
        + model("MP_boundary", "hap1__boundary", 370, 450, [(370, 450)])
        + model(
            "MP_split", "hap1__split", 600, 750, [(600, 650), (700, 750)]
        )
        + model("MP_missing_1", "hap1__missing", 900, 950, [(900, 950)])
        + model(
            "MP_missing_2", "hap2__missing", 900, 950, [(900, 950)], 0.94
        )
        + model(
            "MP_missing_alt",
            "hap1__missing",
            1100,
            1150,
            [(1100, 1150)],
            0.94,
            2,
            96,
        )
        + model(
            "MP_opposite", "hap1__opposite", 1200, 1250, [(1200, 1250)]
        )
        + model("MP_weak", "hap1__weak", 1000, 1050, [(1000, 1050)], 0.5),
    )
    return target, protein_map, miniprot


def test_natural_discovery_classifies_and_clusters_haplotype_support(
    tmp_path: Path,
) -> None:
    target, protein_map, miniprot = make_projection_inputs(tmp_path)
    output = tmp_path / "candidates.tsv"
    manifest = discover_natural_candidates(
        target_gff_path=target,
        miniprot_gff_path=miniprot,
        protein_map_path=protein_map,
        output_tsv_path=output,
    )
    rows = {row["event_type"]: row for row in read_tsv(output)}

    assert set(rows) == {
        MISSING_GENE_EVENT,
        MISSING_INTERNAL_EXON_EVENT,
        BOUNDARY_SHIFT_EVENT,
        SPLIT_GENE_EVENT,
        AMBIGUOUS_OPPOSITE_STRAND_EVENT,
    }
    assert rows[MISSING_INTERNAL_EXON_EVENT]["support_source_count"] == "2"
    assert rows[MISSING_INTERNAL_EXON_EVENT]["discovery_tier"] == (
        "two_haplotype_sources"
    )
    assert rows[MISSING_INTERNAL_EXON_EVENT]["target_gene_ids"] == "G1"
    assert json.loads(rows[MISSING_INTERNAL_EXON_EVENT]["novel_junctions_json"]) == [
        ["chr1", 130, 160],
        ["chr1", 170, 200],
    ]
    assert rows[SPLIT_GENE_EVENT]["target_gene_ids"] == "G3,G4"
    assert rows[MISSING_GENE_EVENT]["support_source_count"] == "2"
    assert rows[MISSING_GENE_EVENT]["mapping_specificity_tier"] == (
        "one_or_more_queries_multilocus"
    )
    assert rows[MISSING_GENE_EVENT]["max_query_near_best_projection_count"] == "2"
    assert rows[BOUNDARY_SHIFT_EVENT]["novel_junction_count"] == "0"
    assert rows[AMBIGUOUS_OPPOSITE_STRAND_EVENT]["target_gene_ids"] == "G5"
    assert manifest["counts"]["candidate_loci"] == 5
    assert manifest["counts"]["filter_counts"] == {
        "concordant_target_chain": 1,
        "low_identity": 1,
        "non_primary_rank": 1,
    }
    assert ALTERNATIVE_STRUCTURE_EVENT not in manifest["counts"]["event_counts"]


def test_natural_rna_validation_is_held_out_and_absence_is_not_contradiction(
    tmp_path: Path,
) -> None:
    target, protein_map, miniprot = make_projection_inputs(tmp_path)
    candidates = tmp_path / "candidates.tsv"
    discover_natural_candidates(
        target_gff_path=target,
        miniprot_gff_path=miniprot,
        protein_map_path=protein_map,
        output_tsv_path=candidates,
    )
    junctions = write(
        tmp_path / "junctions.tsv",
        "seqid\tleft_exon_end\tright_exon_start\t"
        "primary_samples_ge_threshold\tprimary_read_count\t"
        "secondary_samples_ge_threshold\tsecondary_read_count\t"
        "total_samples_ge_threshold\ttotal_read_count\t"
        "primary_support\tsecondary_support\n"
        "chr1\t130\t160\t3\t12\t0\t0\t3\t12\ttrue\tfalse\n"
        "chr1\t170\t200\t2\t8\t0\t0\t2\t8\ttrue\tfalse\n",
    )
    junction_manifest = {
        "schema_version": JUNCTION_AGGREGATE_SCHEMA,
        "output": {"sha256": _file_sha256(junctions)},
    }
    Path(str(junctions) + ".manifest.json").write_text(
        json.dumps(junction_manifest), encoding="utf-8", newline=""
    )
    output = tmp_path / "validated.tsv"
    manifest = validate_natural_candidates_with_rna(
        candidate_tsv_path=candidates,
        junction_aggregate_tsv_path=junctions,
        output_tsv_path=output,
    )
    rows = {row["event_type"]: row for row in read_tsv(output)}

    assert rows[MISSING_INTERNAL_EXON_EVENT]["rna_validation_state"] == (
        "all_novel_junctions_supported"
    )
    assert rows[MISSING_INTERNAL_EXON_EVENT]["rna_validation_tier"] == (
        "two_haplotype_unique_projection_plus_primary_rna"
    )
    assert rows[SPLIT_GENE_EVENT]["rna_validation_state"] == (
        "no_qualifying_junction_observed"
    )
    assert rows[SPLIT_GENE_EVENT]["rna_negative_evidence_policy"] == (
        "absence_is_missing_not_contradiction"
    )
    assert rows[MISSING_GENE_EVENT]["rna_validation_state"] == (
        "not_assessable_no_novel_splice_junction"
    )
    assert manifest["policy"]["candidate_discovery_used_rna"] is False
    assert manifest["counts"]["validation_state_counts"] == {
        "all_novel_junctions_supported": 1,
        "no_qualifying_junction_observed": 1,
        "not_assessable_no_novel_splice_junction": 3,
    }

    summary = summarize_natural_validation(
        validation_tsv_path=output,
        output_json_path=tmp_path / "summary.json",
    )
    event_states = summary["joint_counts"][
        "event_type__x__rna_validation_state"
    ]
    assert {
        (row["event_type"], row["rna_validation_state"]): row["count"]
        for row in event_states
    }[(MISSING_INTERNAL_EXON_EVENT, "all_novel_junctions_supported")] == 1

    graph_candidates = tmp_path / "graph_candidates.tsv"
    graph_evidence = tmp_path / "graph_evidence.tsv"
    graph_manifest = prepare_natural_graph_inputs(
        validation_tsv_path=output,
        output_candidate_tsv_path=graph_candidates,
        output_evidence_tsv_path=graph_evidence,
    )
    graph = infer_event_graph(graph_candidates, graph_evidence)
    graph_rows = {
        row["candidate_id"]: row for row in graph["decisions"]
    }
    natural_rows = {row["event_type"]: row for row in read_tsv(output)}
    internal_id = natural_rows[MISSING_INTERNAL_EXON_EVENT]["candidate_id"]
    split_id = natural_rows[SPLIT_GENE_EVENT]["candidate_id"]
    missing_id = natural_rows[MISSING_GENE_EVENT]["candidate_id"]
    opposite_id = natural_rows[AMBIGUOUS_OPPOSITE_STRAND_EVENT]["candidate_id"]
    assert graph_rows[internal_id]["decision"] == "accept_high_confidence"
    assert graph_rows[split_id]["decision"] == "abstain"
    assert graph_rows[missing_id]["decision"] == "abstain"
    assert graph_rows[opposite_id]["decision"] == "review"
    assert graph_manifest["policy"][
        "haplotype_sources_share_one_independent_group"
    ] is True
    assert graph_manifest["policy"]["rna_absence_generates_contradiction"] is False

    sequence = list("A" * 1500)
    sequence[164] = "N"
    sequence[899] = "a"
    genome = write(tmp_path / "genome.fa", ">chr1\n" + "".join(sequence) + "\n")
    fai = write(tmp_path / "genome.fa.fai", "chr1\t1500\t6\t1500\t1501\n")
    assembly_validation = tmp_path / "validated.assembly.tsv"
    assembly_manifest = annotate_natural_assembly_context(
        validation_tsv_path=output,
        genome_fasta_path=genome,
        genome_fai_path=fai,
        output_tsv_path=assembly_validation,
        flank_bp=100,
    )
    assembly_rows = {
        row["event_type"]: row for row in read_tsv(assembly_validation)
    }
    assert assembly_rows[MISSING_INTERNAL_EXON_EVENT][
        "assembly_context_state"
    ] == "locus_ambiguous_sequence"
    assert assembly_rows[MISSING_GENE_EVENT]["locus_softmasked_bp"] == "1"
    assert assembly_manifest["counts"]["candidates"] == 5

    assembly_graph_candidates = tmp_path / "assembly_graph_candidates.tsv"
    assembly_graph_evidence = tmp_path / "assembly_graph_evidence.tsv"
    prepare_natural_graph_inputs(
        validation_tsv_path=assembly_validation,
        output_candidate_tsv_path=assembly_graph_candidates,
        output_evidence_tsv_path=assembly_graph_evidence,
    )
    assembly_graph = infer_event_graph(
        assembly_graph_candidates, assembly_graph_evidence
    )
    assembly_graph_rows = {
        row["candidate_id"]: row for row in assembly_graph["decisions"]
    }
    assert assembly_graph_rows[internal_id]["decision"] == "abstain"
    assert {
        violation["constraint"]
        for violation in assembly_graph_rows[internal_id]["constraint_violations"]
    } == {"repair_requires_resolved_assembly_context"}

    grouped_junctions = write(
        tmp_path / "grouped_junctions.tsv",
        "seqid\tleft_exon_end\tright_exon_start\tsecondary_read_count\t"
        "secondary_groups_ge_threshold\tsecondary_group_support\n"
        "chr1\t130\t160\t30\t3\ttrue\n"
        "chr1\t170\t200\t20\t2\ttrue\n",
    )
    grouped_manifest = {
        "schema_version": JUNCTION_GROUP_AGGREGATE_SCHEMA,
        "output": {"sha256": _file_sha256(grouped_junctions)},
    }
    Path(str(grouped_junctions) + ".manifest.json").write_text(
        json.dumps(grouped_manifest), encoding="utf-8", newline=""
    )
    secondary_validation = tmp_path / "validated.secondary.tsv"
    secondary_manifest = validate_natural_candidates_with_secondary_groups(
        primary_validation_tsv_path=output,
        grouped_junction_tsv_path=grouped_junctions,
        output_tsv_path=secondary_validation,
    )
    secondary_rows = {
        row["event_type"]: row for row in read_tsv(secondary_validation)
    }
    assert secondary_rows[MISSING_INTERNAL_EXON_EVENT][
        "secondary_group_validation_tier"
    ] == "primary_and_secondary_groups_all_junctions"
    assert secondary_rows[MISSING_INTERNAL_EXON_EVENT][
        "secondary_group_min_supporting_groups"
    ] == "2"
    assert secondary_manifest["policy"][
        "secondary_group_support_is_context_not_truth"
    ] is True
    assert secondary_manifest["policy"]["input_validation_schema"] == (
        "ploidypatch.natural_rna_validation.v2"
    )

    secondary_graph_candidates = tmp_path / "secondary_graph_candidates.tsv"
    secondary_graph_evidence = tmp_path / "secondary_graph_evidence.tsv"
    prepare_natural_graph_inputs(
        validation_tsv_path=secondary_validation,
        output_candidate_tsv_path=secondary_graph_candidates,
        output_evidence_tsv_path=secondary_graph_evidence,
    )
    secondary_evidence_rows = read_tsv(secondary_graph_evidence)
    secondary_edges = [
        row
        for row in secondary_evidence_rows
        if row["source"] == "57_nonprimary_BAM_group_aggregate"
    ]
    assert len(secondary_edges) == 1
    assert {row["direction"] for row in secondary_edges} == {"context"}
    assert secondary_edges[0]["candidate_id"] == internal_id
    secondary_graph = infer_event_graph(
        secondary_graph_candidates, secondary_graph_evidence
    )
    secondary_graph_rows = {
        row["candidate_id"]: row for row in secondary_graph["decisions"]
    }
    assert secondary_graph_rows[internal_id]["decision"] == (
        graph_rows[internal_id]["decision"]
    )

    assembly_secondary_validation = tmp_path / "validated.assembly.secondary.tsv"
    assembly_secondary_manifest = validate_natural_candidates_with_secondary_groups(
        primary_validation_tsv_path=assembly_validation,
        grouped_junction_tsv_path=grouped_junctions,
        output_tsv_path=assembly_secondary_validation,
    )
    assert assembly_secondary_manifest["policy"]["input_validation_schema"] == (
        "ploidypatch.natural_assembly_context.v1"
    )
    assembly_secondary_rows = {
        row["event_type"]: row for row in read_tsv(assembly_secondary_validation)
    }
    assert assembly_secondary_rows[MISSING_INTERNAL_EXON_EVENT][
        "assembly_context_state"
    ] == "locus_ambiguous_sequence"


def test_natural_cli_freezes_conservative_defaults() -> None:
    args = build_parser().parse_args(
        [
            "evidence",
            "discover-natural",
            "--target-gff",
            "target.gff3",
            "--miniprot-gff",
            "miniprot.gff3",
            "--protein-map",
            "proteins.tsv",
            "--output-tsv",
            "candidates.tsv",
        ]
    )
    assert args.min_identity == 0.8
    assert args.min_query_coverage == 0.8
    assert args.max_existing_cds_overlap == 0.1
    assert args.min_boundary_extension_bp == 30
    assert args.near_best_score_fraction == 0.95

    assembly_args = build_parser().parse_args(
        [
            "evidence",
            "annotate-natural-assembly",
            "--validation-tsv",
            "validated.tsv",
            "--genome",
            "genome.fa",
            "--fai",
            "genome.fa.fai",
            "--output-tsv",
            "context.tsv",
        ]
    )
    assert assembly_args.flank_bp == 5000
    assert assembly_args.max_ambiguous_fraction == 0.0
