from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import pytest

from ploidypatch.species_applicability import evaluate_species_applicability


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "build_species_applicability_metrics_v0.6.py"


def load_script():
    spec = importlib.util.spec_from_file_location("applicability_metrics_builder", SCRIPT)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def write_json(path: Path, value: dict) -> None:
    path.write_text(json.dumps(value) + "\n", encoding="utf-8")


def make_inputs(tmp_path: Path) -> dict[str, Path]:
    full_genome = tmp_path / "genome.fa"
    primary_genome = tmp_path / "primary.fa"
    full_genome.write_text(">chr1\n" + "A" * 1000 + "\n", encoding="utf-8")
    primary_genome.write_text(">chr1\n" + "A" * 1000 + "\n", encoding="utf-8")
    fai = tmp_path / "primary.fa.fai"
    fai.write_text("chr1\t1000\t6\t1000\t1001\n", encoding="utf-8")
    stats = tmp_path / "stats.tsv"
    stats.write_text(
        "file\tformat\ttype\tnum_seqs\tsum_len\tsum_n\n"
        f"{full_genome}\tFASTA\tDNA\t1\t1000\t0\n"
        f"{primary_genome}\tFASTA\tDNA\t1\t1000\t0\n",
        encoding="utf-8",
    )
    gff_text = (
        "##gff-version 3\n"
        "chr1\tx\tgene\t1\t900\t.\t+\t.\tID=g1\n"
        "chr1\tx\tmRNA\t1\t900\t.\t+\t.\tID=t1;Parent=g1\n"
        "chr1\tx\tCDS\t1\t300\t.\t+\t0\tParent=t1\n"
    )
    full_gff = tmp_path / "full.gff3"
    primary_gff = tmp_path / "primary.gff3"
    full_gff.write_text(gff_text, encoding="utf-8")
    primary_gff.write_text(gff_text, encoding="utf-8")
    protein = tmp_path / "protein.json"
    write_json(
        protein,
        {
            "schema_version": "ploidypatch.primary_protein_quality_subset.v0.6",
            "candidate_access": False,
            "truth_access": False,
            "mapping_policy": (
                "exact_transcript_relation_ID_after_repeated_recognized_namespace_stripping"
            ),
            "fractions": {"exact_unique_GFF_protein_mapping_fraction": 1.0},
        },
    )
    translation = tmp_path / "translation.json"
    write_json(
        translation,
        {
            "schema_version": "ploidypatch.primary_translation_quality.v0.6",
            "candidate_access": False,
            "truth_access": False,
            "comparison_policy": (
                "exact_amino_acid_identity_after_terminal_stop_removal_and_"
                "repeated_recognized_namespace_stripping"
            ),
            "fractions": {"valid_representative_translation_fraction": 1.0},
        },
    )
    protein_busco = tmp_path / "protein_busco.json"
    genome_busco = tmp_path / "genome_busco.json"
    write_json(
        protein_busco,
        {
            "parameters": {"mode": "proteins"},
            "versions": {"busco": "6.1.0"},
            "lineage_dataset": {"name": "test_odb12", "number_of_buscos": "10"},
            "results": {"Complete percentage": 98.0, "n_markers": 10},
        },
    )
    write_json(
        genome_busco,
        {
            "parameters": {
                "mode": "euk_genome_min",
                "gene_predictor": "miniprot",
                "use_miniprot": "True",
            },
            "versions": {"busco": "6.1.0"},
            "lineage_dataset": {"name": "test_odb12", "number_of_buscos": "10"},
            "results": {"Complete percentage": 98.0, "n_markers": 10},
        },
    )
    backbone = tmp_path / "backbone.json"
    write_json(
        backbone,
        {
            "schema_version": "ploidypatch.fixed_target_backbone.v0.6",
            "candidate_access": False,
            "truth_access": False,
            "applicable": True,
            "identifier_mapping": {
                "mapping_policy": (
                    "unique_exact_attribute_then_unique_normalized_feature_ID_fallback"
                ),
                "bidirectional_unique": True,
            },
            "policy": {"min_block_pairs": 20, "require_cross_seqid": True},
            "semantic_invariants": {
                "input_permutation_deterministic": True,
                "reverse_duplicate_statistics_audit_only": True,
                "ambiguous_endpoint_block_quarantine": True,
            },
            "applicability_metrics": {
                "primary_gene_midpoint_backbone_coverage": 0.8,
                "primary_chromosome_cell_coverage_fraction": 1.0,
                "minimum_cells_per_covered_chromosome_observed": 25,
                "unique_partner_chromosome_fraction_among_covered_genes": 0.9,
            },
        },
    )
    sources = tmp_path / "sources.tsv"
    sources.write_text(
        "source_id\tcitation\tdoi\turl\tevidence_scope\tindependent_evidence_class\tqualifies_for_applicability\n"
        "s1\tA\t10.1/a\thttps://doi.org/10.1/a\tx\tgenome\ttrue\n"
        "s2\tB\t10.1/b\thttps://doi.org/10.1/b\tx\tgene_tree\ttrue\n",
        encoding="utf-8",
    )
    lineage_root = tmp_path / "lineage"
    lineage_root.mkdir()
    lineage_file = lineage_root / "dataset.cfg"
    lineage_file.write_text("dataset=test\n", encoding="utf-8")
    import hashlib
    lineage_digest = hashlib.sha256(lineage_file.read_bytes()).hexdigest()
    lineage_sums = lineage_root / "SHA256SUMS"
    lineage_sums.write_text(lineage_digest + "  dataset.cfg\n", encoding="utf-8")
    implementation_archive = tmp_path / "source.tar.gz"
    implementation_archive.write_bytes(b"frozen source")
    implementation_digest = hashlib.sha256(implementation_archive.read_bytes()).hexdigest()
    implementation_sums = tmp_path / "implementation.SHA256SUMS"
    implementation_sums.write_text(
        implementation_digest + "  source.tar.gz\n", encoding="utf-8"
    )
    return locals()


def run_builder(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> Path:
    values = make_inputs(tmp_path)
    output = tmp_path / "metrics.json"
    argv = [
        str(SCRIPT),
        "--species-id", "synthetic_species",
        "--declared-primary-seqid-count", "1",
        "--full-genome", str(values["full_genome"]),
        "--primary-genome", str(values["primary_genome"]),
        "--primary-fai", str(values["fai"]),
        "--seqkit-stats", str(values["stats"]),
        "--full-gff", str(values["full_gff"]),
        "--primary-gff", str(values["primary_gff"]),
        "--protein-subset-manifest", str(values["protein"]),
        "--translation-quality", str(values["translation"]),
        "--protein-busco-summary", str(values["protein_busco"]),
        "--genome-busco-summary", str(values["genome_busco"]),
        "--busco-lineage-root", str(values["lineage_root"]),
        "--busco-lineage-sha256sums", str(values["lineage_sums"]),
        "--expected-busco-version", "6.1.0",
        "--expected-busco-lineage", "test_odb12",
        "--expected-busco-markers", "10",
        "--backbone-manifest", str(values["backbone"]),
        "--wgd-source-registry", str(values["sources"]),
        "--implementation-source-archive", str(values["implementation_archive"]),
        "--implementation-freeze-sha256sums", str(values["implementation_sums"]),
        "--output-json", str(output),
    ]
    monkeypatch.setattr(sys, "argv", argv)
    assert load_script().main() == 0
    return output


def test_builder_produces_metrics_accepted_by_three_state_evaluator(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    metrics = run_builder(monkeypatch, tmp_path)
    report = evaluate_species_applicability(
        metrics_path=metrics,
        policy_path=ROOT / "config" / "species_applicability_policy_v0.6.tsv",
        output_path=tmp_path / "report.json",
    )
    assert report["applicability_state"] == "full_topology_evaluable"
    assert report["metrics"]["annotation"]["valid_coding_hierarchy_fraction"] == 1.0
    assert report["metrics"]["backbone"]["independent_WGD_source_count"] == 2


def test_hierarchy_fraction_marks_out_of_bounds_cds_invalid(tmp_path: Path) -> None:
    gff = tmp_path / "bad.gff3"
    gff.write_text(
        "chr1\tx\tgene\t1\t100\t.\t+\t.\tID=g1\n"
        "chr1\tx\tmRNA\t1\t100\t.\t+\t.\tID=t1;Parent=g1\n"
        "chr1\tx\tCDS\t1\t101\t.\t+\t0\tParent=t1\n",
        encoding="utf-8",
    )
    report = load_script().parse_coding_hierarchy(gff, {"chr1": 100})
    assert report["coding_transcripts"] == 1
    assert report["valid_coding_hierarchy_fraction"] == 0.0
    assert report["invalid_reasons"] == {"invalid_CDS_geometry": 1}


def test_malformed_feature_line_fails_hierarchy_fraction_closed(tmp_path: Path) -> None:
    gff = tmp_path / "malformed.gff3"
    gff.write_text(
        "chr1\tx\tgene\t1\t100\t.\t+\t.\tID=g1\n"
        "chr1\tx\tmRNA\t1\t100\t.\t+\t.\tID=t1;Parent=g1\n"
        "chr1\tx\tCDS\t1\t99\t.\t+\t0\tParent=t1\n"
        "not-a-nine-column-GFF-row\n",
        encoding="utf-8",
    )
    report = load_script().parse_coding_hierarchy(gff, {"chr1": 100})
    assert report["valid_coding_transcripts"] == 1
    assert report["malformed_rows_or_attributes"] == 1
    assert report["valid_coding_hierarchy_fraction"] == 0.0


def test_wgd_registry_requires_independent_unique_sources(tmp_path: Path) -> None:
    registry = tmp_path / "sources.tsv"
    registry.write_text(
        "source_id\tcitation\tdoi\turl\tevidence_scope\tindependent_evidence_class\tqualifies_for_applicability\n"
        "s1\tA\t10.1/a\tu\tx\tgenome\ttrue\n"
        "s2\tB\t10.1/a\tu\tx\tgene_tree\ttrue\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="unique nonempty doi"):
        load_script().independent_wgd_sources(registry)
