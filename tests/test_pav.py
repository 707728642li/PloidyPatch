from __future__ import annotations

import csv
from pathlib import Path

import pytest

from ploidypatch.pav import (
    annotate_pav_with_wgdi,
    catalog_genes_in_target_deletions,
    collapse_high_confidence_pav_events,
    combine_chromosome_pafs,
    extract_paf_target_deletions,
    reconcile_pav_gene_catalogs,
    screen_pav_candidate_proteins,
    subset_catalog_proteins,
)


GFF = (
    "##gff-version 3\n"
    "t1\ttest\tgene\t5001\t7000\t.\t+\t.\tID=gene:G1\n"
    "t1\ttest\tmRNA\t5001\t7000\t.\t+\t.\tID=transcript:T1;Parent=gene:G1\n"
    "t1\ttest\texon\t5101\t5300\t.\t+\t.\tParent=transcript:T1\n"
    "t1\ttest\tCDS\t5101\t5300\t.\t+\t0\tParent=transcript:T1\n"
    "t1\ttest\tgene\t6801\t7800\t.\t+\t.\tID=gene:G2\n"
    "t1\ttest\tmRNA\t6801\t7800\t.\t+\t.\tID=transcript:T2;Parent=gene:G2\n"
    "t1\ttest\texon\t6801\t6900\t.\t+\t.\tParent=transcript:T2\n"
    "t1\ttest\tCDS\t6801\t6900\t.\t+\t0\tParent=transcript:T2\n"
)


def write_indexed_fasta(path: Path, records: dict[str, str]) -> Path:
    fai = Path(str(path) + ".fai")
    rows = []
    with path.open("wb") as handle:
        for seqid, sequence in records.items():
            handle.write(f">{seqid}\n".encode("ascii"))
            offset = handle.tell()
            for index in range(0, len(sequence), 100):
                handle.write(sequence[index : index + 100].encode("ascii") + b"\n")
            rows.append(f"{seqid}\t{len(sequence)}\t{offset}\t100\t101\n")
    fai.write_text("".join(rows), encoding="utf-8")
    return fai


def test_extracts_only_well_anchored_primary_target_deletions(
    tmp_path: Path,
) -> None:
    paf = tmp_path / "alignments.paf"
    paf.write_text(
        "q1\t10000\t0\t10000\t+\tt1\t12000\t0\t12000\t"
        "9800\t12000\t60\ttp:A:P\tcg:Z:5000M2000D5000M\n"
        "q2\t2000\t0\t2000\t+\tt2\t2300\t0\t2300\t"
        "1900\t2300\t60\ttp:A:P\tcg:Z:100M300D1900M\n"
        "q3\t10000\t0\t10000\t+\tt3\t12000\t0\t12000\t"
        "9800\t12000\t60\ttp:A:S\tcg:Z:5000M2000D5000M\n",
        encoding="utf-8",
    )
    output = tmp_path / "deletions.tsv"
    pairs = tmp_path / "pairs.tsv"
    pairs.write_text(
        "query_seqid\ttarget_seqid\nq1\tt1\nq2\tt2\nq3\tt3\n",
        encoding="utf-8",
    )

    manifest = extract_paf_target_deletions(
        paf_path=paf,
        output_tsv_path=output,
        min_deletion_bp=200,
        min_flanking_anchor_bp=1000,
        min_mapq=20,
        chromosome_pair_tsv_path=pairs,
    )

    with output.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 1
    assert rows[0]["target_seqid"] == "t1"
    assert rows[0]["target_start"] == "5001"
    assert rows[0]["target_end"] == "7000"
    assert rows[0]["query_breakpoint"] == "5000"
    assert rows[0]["left_anchor_bp"] == "5000"
    assert rows[0]["right_anchor_bp"] == "5000"
    assert manifest["counts"]["alignments_non_primary"] == 1
    assert manifest["counts"]["target_deletions_below_anchor"] == 1
    assert manifest["output"]["rows"] == 1
    assert manifest["input"]["chromosome_pairs"]["pairs"] == 3

    gff = tmp_path / "source.gff3"
    gff.write_text(GFF, encoding="utf-8")
    genes = tmp_path / "pav_genes.tsv"
    gene_manifest = catalog_genes_in_target_deletions(
        gff_path=gff,
        deletion_tsv_path=output,
        output_tsv_path=genes,
        min_gene_coverage=0.95,
        min_cds_coverage=1.0,
    )
    with genes.open("r", encoding="utf-8", newline="") as handle:
        gene_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [row["gene_id"] for row in gene_rows] == ["gene:G1"]
    assert gene_manifest["counts"]["coding_genes_overlapping_deletion"] == 2
    assert gene_manifest["counts"]["genes_below_gene_coverage"] == 1
    assert gene_manifest["counts"]["genes_accepted"] == 1


def test_combines_chromosome_pafs_in_declared_order(tmp_path: Path) -> None:
    paf_dir = tmp_path / "by_chr"
    paf_dir.mkdir()
    paf_dir.joinpath("q2.paf").write_text(
        "q2\t10\t0\t10\t+\tt2\t10\t0\t10\t10\t10\t60\ttp:A:P\tcg:Z:10M\n",
        encoding="utf-8",
    )
    paf_dir.joinpath("q1.paf").write_text(
        "q1\t10\t0\t10\t+\tt9\t10\t0\t10\t10\t10\t60\ttp:A:P\tcg:Z:10M\n",
        encoding="utf-8",
    )
    pairs = tmp_path / "pairs.tsv"
    pairs.write_text(
        "query_seqid\ttarget_seqid\nq2\tt2\nq1\tt1\n", encoding="utf-8"
    )
    output = tmp_path / "combined.paf"
    manifest = combine_chromosome_pafs(
        input_dir_path=paf_dir,
        chromosome_pair_tsv_path=pairs,
        output_paf_path=output,
    )
    lines = output.read_text(encoding="utf-8").splitlines()
    assert [line.split("\t", 1)[0] for line in lines] == ["q2", "q1"]
    assert manifest["counts"]["alignments_total"] == 2
    assert manifest["counts"]["alignments_expected_target"] == 1
    assert manifest["counts"]["alignments_other_target"] == 1
    assert [item["query_seqid"] for item in manifest["input"]["pafs"]] == [
        "q2",
        "q1",
    ]


def test_query_flanks_exclude_ambiguous_and_truncated_breakpoints(
    tmp_path: Path,
) -> None:
    genome = tmp_path / "query.fa"
    fai = write_indexed_fasta(
        genome,
        {
            "q_clean": "A" * 10000,
            "q_ambiguous": "A" * 4999 + "N" + "A" * 5000,
            "q_edge": "A" * 2000,
        },
    )
    paf = tmp_path / "alignments.paf"
    paf.write_text(
        "q_clean\t10000\t0\t10000\t+\tt1\t10200\t0\t10200\t"
        "9900\t10200\t60\ttp:A:P\tcg:Z:5000M200D5000M\n"
        "q_ambiguous\t10000\t0\t10000\t+\tt2\t10200\t0\t10200\t"
        "9900\t10200\t60\ttp:A:P\tcg:Z:5000M200D5000M\n"
        "q_edge\t2000\t0\t2000\t+\tt3\t2200\t0\t2200\t"
        "1950\t2200\t60\ttp:A:P\tcg:Z:500M200D1500M\n",
        encoding="utf-8",
    )
    pairs = tmp_path / "pairs.tsv"
    pairs.write_text(
        "query_seqid\ttarget_seqid\nq_clean\tt1\nq_ambiguous\tt2\nq_edge\tt3\n",
        encoding="utf-8",
    )
    output = tmp_path / "clean.tsv"
    manifest = extract_paf_target_deletions(
        paf_path=paf,
        output_tsv_path=output,
        min_deletion_bp=200,
        min_flanking_anchor_bp=100,
        chromosome_pair_tsv_path=pairs,
        query_genome_path=genome,
        query_fai_path=fai,
        query_flank_bp=1000,
    )
    with output.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [row["query_seqid"] for row in rows] == ["q_clean"]
    assert rows[0]["query_flank_status"] == "pass"
    assert rows[0]["query_flank_observed_bp"] == "2000"
    assert manifest["counts"]["target_deletions_query_flank_pass"] == 1
    assert manifest["counts"]["target_deletions_query_flank_ambiguous_sequence"] == 1
    assert (
        manifest["counts"]["target_deletions_query_flank_truncated_at_sequence_end"]
        == 1
    )
    assert manifest["counts"]["target_deletions_excluded_by_query_flank"] == 2
    assert manifest["input"]["query_fai"]["sequences"] == 3

    retained = tmp_path / "retained.tsv"
    extract_paf_target_deletions(
        paf_path=paf,
        output_tsv_path=retained,
        min_deletion_bp=200,
        min_flanking_anchor_bp=100,
        chromosome_pair_tsv_path=pairs,
        query_genome_path=genome,
        query_fai_path=fai,
        query_flank_bp=1000,
        require_clean_query_flanks=False,
    )
    with retained.open("r", encoding="utf-8", newline="") as handle:
        retained_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert {row["query_flank_status"] for row in retained_rows} == {
        "pass",
        "ambiguous_sequence",
        "truncated_at_sequence_end",
    }


def test_query_flank_inputs_must_be_paired(tmp_path: Path) -> None:
    paf = tmp_path / "empty.paf"
    paf.write_text("", encoding="utf-8")
    with pytest.raises(ValueError, match="supplied together"):
        extract_paf_target_deletions(
            paf_path=paf,
            output_tsv_path=tmp_path / "out.tsv",
            query_genome_path=tmp_path / "query.fa",
        )


def write_pav_catalog(path: Path, genes: list[str], deletion_offset: int) -> None:
    fieldnames = [
        "candidate_id",
        "gene_id",
        "seqid",
        "gene_start",
        "gene_end",
        "strand",
        "gene_span_bp",
        "transcript_structures",
        "cds_union_bp",
        "deletion_id",
        "deletion_start",
        "deletion_end",
        "deletion_bp",
        "gene_coverage",
        "cds_coverage",
        "query_seqid",
        "query_breakpoint",
        "alignment_strand",
        "mapq",
        "alignment_identity",
        "left_anchor_bp",
        "right_anchor_bp",
        "query_flank_ambiguous_fraction",
        "query_flank_status",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for index, gene_id in enumerate(genes, start=1):
            start = index * 1000
            writer.writerow(
                {
                    "candidate_id": f"candidate-{gene_id}",
                    "gene_id": gene_id,
                    "seqid": "chr1",
                    "gene_start": start,
                    "gene_end": start + 199,
                    "strand": "+",
                    "gene_span_bp": 200,
                    "transcript_structures": 1,
                    "cds_union_bp": 150,
                    "deletion_id": f"deletion-{gene_id}-{deletion_offset}",
                    "deletion_start": start - deletion_offset,
                    "deletion_end": start + 199 + deletion_offset,
                    "deletion_bp": 200 + 2 * deletion_offset,
                    "gene_coverage": "1.00000000",
                    "cds_coverage": "1.00000000",
                    "query_seqid": "q1",
                    "query_breakpoint": 5000 + index,
                    "alignment_strand": "+",
                    "mapq": 60,
                    "alignment_identity": "0.99000000",
                    "left_anchor_bp": 1000,
                    "right_anchor_bp": 1000,
                    "query_flank_ambiguous_fraction": "0.00000000",
                    "query_flank_status": "pass",
                }
            )


def test_reconciles_only_genes_supported_by_multiple_catalogs(tmp_path: Path) -> None:
    asm10 = tmp_path / "asm10.tsv"
    asm20 = tmp_path / "asm20.tsv"
    write_pav_catalog(asm10, ["G1", "G2"], 20)
    write_pav_catalog(asm20, ["G1", "G3"], 30)
    output = tmp_path / "reconciled.tsv"
    manifest = reconcile_pav_gene_catalogs(
        catalog_inputs=[f"asm10={asm10}", f"asm20={asm20}"],
        output_tsv_path=output,
        min_support=2,
    )
    with output.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [row["gene_id"] for row in rows] == ["G1"]
    assert rows[0]["support_sources"] == "asm10,asm20"
    assert rows[0]["asm10_deletion_start"] == "980"
    assert rows[0]["asm20_deletion_start"] == "970"
    assert rows[0]["asm10_query_flank_status"] == "pass"
    assert manifest["counts"]["genes_with_2_source_support"] == 1
    assert manifest["counts"]["genes_below_minimum_support"] == 2
    assert manifest["output"]["rows"] == 1


def test_reconciler_rejects_duplicate_labels(tmp_path: Path) -> None:
    catalog = tmp_path / "catalog.tsv"
    write_pav_catalog(catalog, ["G1"], 10)
    with pytest.raises(ValueError, match="Duplicate catalog label"):
        reconcile_pav_gene_catalogs(
            catalog_inputs=[f"same={catalog}", f"same={catalog}"],
            output_tsv_path=tmp_path / "out.tsv",
        )


def write_reconciled_catalog(path: Path, genes: list[str]) -> None:
    fieldnames = [
        "candidate_id",
        "gene_id",
        "seqid",
        "gene_start",
        "gene_end",
        "strand",
        "gene_span_bp",
        "transcript_structures",
        "cds_union_bp",
        "support_count",
        "support_sources",
        "asm10_query_seqid",
        "asm10_query_breakpoint",
        "asm20_query_seqid",
        "asm20_query_breakpoint",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for index, gene_id in enumerate(genes, start=1):
            writer.writerow(
                {
                    "candidate_id": f"PPRPAV-{gene_id}",
                    "gene_id": gene_id,
                    "seqid": "chr1",
                    "gene_start": index * 1000,
                    "gene_end": index * 1000 + 299,
                    "strand": "+",
                    "gene_span_bp": 300,
                    "transcript_structures": 1,
                    "cds_union_bp": 300,
                    "support_count": 2,
                    "support_sources": "asm10,asm20",
                    "asm10_query_seqid": "chr1",
                    "asm10_query_breakpoint": index * 1000,
                    "asm20_query_seqid": "chr1",
                    "asm20_query_breakpoint": index * 1000,
                }
            )


def test_subsets_candidate_proteins_in_catalog_order(tmp_path: Path) -> None:
    catalog = tmp_path / "catalog.tsv"
    write_reconciled_catalog(catalog, ["G1", "G2"])
    proteins = tmp_path / "representatives.fa"
    proteins.write_text(
        ">G2 source two\nMMMM*\n>irrelevant\nAAAA\n>G1 source one\nKKKK\n",
        encoding="utf-8",
    )
    output = tmp_path / "candidate.fa"
    manifest = subset_catalog_proteins(
        catalog_tsv_path=catalog,
        protein_fasta_path=proteins,
        output_fasta_path=output,
    )
    assert output.read_text(encoding="utf-8") == ">G1\nKKKK\n>G2\nMMMM\n"
    assert manifest["output"]["records"] == 2
    assert manifest["protein_lengths_aa"]["total"] == 8


def test_screens_full_partial_and_disrupted_miniprot_evidence(tmp_path: Path) -> None:
    catalog = tmp_path / "catalog.tsv"
    write_reconciled_catalog(catalog, ["G1", "G2", "G3", "G4", "G5"])
    protein_map = tmp_path / "proteins.map.tsv"
    protein_map.write_text(
        "query_id\tsource\tsource_record_id\tlength_aa\tsource_header\n"
        "src__G1\tsrc\tG1\t100\tG1\n"
        "src__G2\tsrc\tG2\t100\tG2\n"
        "src__G3\tsrc\tG3\t100\tG3\n"
        "src__G4\tsrc\tG4\t100\tG4\n"
        "src__G5\tsrc\tG5\t100\tG5\n",
        encoding="utf-8",
    )
    miniprot = tmp_path / "miniprot.gff3"
    miniprot.write_text(
        "##gff-version 3\n"
        "chr1\tminiprot\tmRNA\t1\t300\t100\t+\t.\tID=M1;Rank=1;Identity=0.9;Target=src__G1 1 100\n"
        "chr1\tminiprot\tCDS\t1\t300\t100\t+\t0\tParent=M1\n"
        "chr1\tminiprot\tmRNA\t1401\t1700\t90\t+\t.\tID=M2;Rank=1;Identity=0.8;Frameshift=1;Target=src__G2 1 100\n"
        "chr1\tminiprot\tCDS\t1401\t1700\t90\t+\t0\tParent=M2\n"
        "scaffold1\tminiprot\tmRNA\t1\t90\t80\t-\t.\tID=M3;Rank=1;Identity=0.8;Target=src__G3 1 30\n"
        "scaffold1\tminiprot\tCDS\t1\t90\t80\t-\t0\tParent=M3\n"
        "chr1\tminiprot\tmRNA\t801\t1100\t70\t+\t.\tID=M4;Rank=1;Identity=0.4;Target=src__G4 1 100\n"
        "chr1\tminiprot\tCDS\t801\t1100\t70\t+\t0\tParent=M4\n"
        "chr2\tminiprot\tmRNA\t1\t300\t100\t+\t.\tID=M5;Rank=1;Identity=0.9;Target=src__G5 1 100\n"
        "chr2\tminiprot\tCDS\t1\t300\t100\t+\t0\tParent=M5\n",
        encoding="utf-8",
    )
    pairs = tmp_path / "pairs.tsv"
    pairs.write_text(
        "query_seqid\ttarget_seqid\nchr1\tt1\nchr2\tt2\n", encoding="utf-8"
    )
    output = tmp_path / "screen.tsv"
    manifest = screen_pav_candidate_proteins(
        catalog_tsv_path=catalog,
        miniprot_gff_path=miniprot,
        protein_map_path=protein_map,
        chromosome_pair_tsv_path=pairs,
        output_tsv_path=output,
        min_identity=0.5,
        min_query_coverage=0.5,
        min_partial_query_coverage=0.2,
        max_local_distance_bp=1000,
    )
    with output.open("r", encoding="utf-8", newline="") as handle:
        rows = {row["gene_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
    assert rows["G1"]["screen_class"] == "local_intact_homolog"
    assert rows["G1"]["screen_decision"] == "exclude_local_sequence_present"
    assert rows["G2"]["screen_class"] == "local_disrupted_homolog"
    assert rows["G3"]["screen_class"] == "unlocalized_partial_homolog"
    assert rows["G3"]["best_seqid"] == "scaffold1"
    assert rows["G4"]["screen_class"] == "no_substantial_homolog"
    assert (
        rows["G4"]["screen_decision"] == "retain_sequence_absence_candidate"
    )
    assert rows["G5"]["screen_class"] == "off_locus_intact_homolog_only"
    assert rows["G5"]["screen_decision"] == "retain_sequence_absence_candidate"
    assert manifest["counts"]["screen_class_local_intact_homolog"] == 1
    assert manifest["counts"]["screen_class_local_disrupted_homolog"] == 1
    assert manifest["counts"]["screen_class_unlocalized_partial_homolog"] == 1
    assert manifest["counts"]["screen_class_no_substantial_homolog"] == 1
    assert manifest["counts"]["screen_class_off_locus_intact_homolog_only"] == 1


def test_protein_screen_abstains_on_breakpoint_disagreement(tmp_path: Path) -> None:
    catalog = tmp_path / "catalog.tsv"
    write_reconciled_catalog(catalog, ["G1"])
    with catalog.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = list(reader.fieldnames or [])
        rows = list(reader)
    rows[0]["asm20_query_seqid"] = "chr2"
    with catalog.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    protein_map = tmp_path / "proteins.map.tsv"
    protein_map.write_text(
        "query_id\tsource\tsource_record_id\tlength_aa\tsource_header\n"
        "src__G1\tsrc\tG1\t100\tG1\n",
        encoding="utf-8",
    )
    miniprot = tmp_path / "empty.gff3"
    miniprot.write_text("##gff-version 3\n", encoding="utf-8")
    pairs = tmp_path / "pairs.tsv"
    pairs.write_text(
        "query_seqid\ttarget_seqid\nchr1\tt1\nchr2\tt2\n", encoding="utf-8"
    )
    output = tmp_path / "screen.tsv"
    screen_pav_candidate_proteins(
        catalog_tsv_path=catalog,
        miniprot_gff_path=miniprot,
        protein_map_path=protein_map,
        chromosome_pair_tsv_path=pairs,
        output_tsv_path=output,
    )
    with output.open("r", encoding="utf-8", newline="") as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    assert row["screen_class"] == "breakpoint_disagreement"
    assert row["screen_decision"] == "abstain_breakpoint_disagreement"


def test_wgdi_annotation_only_promotes_expected_progenitor_support(
    tmp_path: Path,
) -> None:
    screen = tmp_path / "screen.tsv"
    screen.write_text(
        "gene_id\tscreen_class\tscreen_decision\n"
        "G1\toff_locus_intact_homolog_only\tretain_sequence_absence_candidate\n"
        "G2\toff_locus_intact_homolog_only\tretain_sequence_absence_candidate\n"
        "G3\tno_substantial_homolog\tretain_sequence_absence_candidate\n"
        "G4\tlocal_intact_homolog\texclude_local_sequence_present\n"
        "G5\tunlocalized_intact_homolog\tabstain_unlocalized_sequence_evidence\n",
        encoding="utf-8",
    )
    wgdi = tmp_path / "wgdi.tsv"
    wgdi.write_text(
        "gene_id\tsubgenome\texpected_progenitor\t"
        "expected_progenitor_supported\tcross_progenitor_supported\t"
        "synteny_stratum\n"
        "G1\tA\tbra_a\ttrue\ttrue\texpected_and_cross\n"
        "G2\tA\tbra_a\tfalse\ttrue\tcross_only\n"
        "G3\tC\tbol_c\tfalse\tfalse\tno_block\n"
        "G4\tC\tbol_c\ttrue\tfalse\texpected_only\n"
        "G5\tA\tbra_a\ttrue\ttrue\texpected_and_cross\n",
        encoding="utf-8",
    )
    output = tmp_path / "annotated.tsv"
    manifest = annotate_pav_with_wgdi(
        protein_screen_tsv_path=screen,
        wgdi_gene_tsv_path=wgdi,
        output_tsv_path=output,
    )
    with output.open("r", encoding="utf-8", newline="") as handle:
        rows = {row["gene_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
    assert rows["G1"]["plant_pav_decision"] == "retain_high_confidence_candidate"
    assert rows["G2"]["plant_pav_decision"] == "abstain_cross_only_synteny"
    assert rows["G3"]["plant_pav_decision"] == "abstain_no_synteny_support"
    assert rows["G4"]["plant_pav_evidence_class"] == "sequence_presence_contradiction"
    assert (
        rows["G5"]["plant_pav_decision"]
        == "abstain_unlocalized_sequence_evidence"
    )
    assert rows["G1"]["wgdi_synteny_stratum"] == "expected_and_cross"
    assert manifest["counts"]["plant_pav_decision_retain_high_confidence_candidate"] == 1


def test_collapses_gene_candidates_into_structural_events(tmp_path: Path) -> None:
    annotated = tmp_path / "annotated.tsv"
    fieldnames = [
        "gene_id",
        "gene_start",
        "seqid",
        "support_sources",
        "screen_class",
        "plant_pav_decision",
        "wgdi_subgenome",
        "wgdi_chromosome_label",
        "wgdi_synteny_stratum",
    ]
    for source in ("asm10", "asm20"):
        fieldnames.extend(
            f"{source}_{suffix}"
            for suffix in (
                "deletion_id",
                "deletion_start",
                "deletion_end",
                "deletion_bp",
                "query_seqid",
                "query_breakpoint",
            )
        )
    rows = []
    for gene_id, gene_start, event_number, decision in (
        ("G1", 150, 1, "retain_high_confidence_candidate"),
        ("G2", 250, 1, "retain_high_confidence_candidate"),
        ("G3", 1100, 2, "retain_high_confidence_candidate"),
        ("G4", 1200, 2, "abstain_no_synteny_support"),
    ):
        deletion_start, deletion_end = ((100, 500) if event_number == 1 else (1000, 1300))
        row = {
            "gene_id": gene_id,
            "gene_start": gene_start,
            "seqid": "t1",
            "support_sources": "asm10,asm20",
            "screen_class": "off_locus_intact_homolog_only",
            "plant_pav_decision": decision,
            "wgdi_subgenome": "A",
            "wgdi_chromosome_label": "A1",
            "wgdi_synteny_stratum": "expected_and_cross",
        }
        for source in ("asm10", "asm20"):
            row.update(
                {
                    f"{source}_deletion_id": f"{source}-D{event_number}",
                    f"{source}_deletion_start": deletion_start,
                    f"{source}_deletion_end": deletion_end,
                    f"{source}_deletion_bp": deletion_end - deletion_start + 1,
                    f"{source}_query_seqid": "q1",
                    f"{source}_query_breakpoint": event_number * 100,
                }
            )
        rows.append(row)
    with annotated.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    output = tmp_path / "events.tsv"
    manifest = collapse_high_confidence_pav_events(
        annotated_gene_tsv_path=annotated,
        output_tsv_path=output,
    )
    with output.open("r", encoding="utf-8", newline="") as handle:
        events = list(csv.DictReader(handle, delimiter="\t"))
    assert len(events) == 2
    assert [event["gene_count"] for event in events] == ["2", "1"]
    assert events[0]["gene_ids"] == "G1,G2"
    assert events[0]["source_count"] == "2"
    assert events[0]["target_boundary_spread_bp"] == "0"
    assert manifest["counts"]["events_total"] == 2
    assert manifest["counts"]["genes_total"] == 3
    assert manifest["counts"]["events_with_2_genes"] == 1
