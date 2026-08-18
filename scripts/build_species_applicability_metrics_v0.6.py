#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from ploidypatch.artifact_manifest import read_sha256sums, verify_sha256sums
from ploidypatch.gff import TRANSCRIPT_TYPES, parse_attributes
from ploidypatch.io import open_text


SCHEMA_VERSION = "ploidypatch.species_applicability_metrics.v0.6"
BACKBONE_SCHEMA = "ploidypatch.fixed_target_backbone.v0.6"
PROTEIN_SUBSET_SCHEMA = "ploidypatch.primary_protein_quality_subset.v0.6"
TRANSLATION_SCHEMA = "ploidypatch.primary_translation_quality.v0.6"


@dataclass(frozen=True)
class Feature:
    seqid: str
    start: int
    end: int
    strand: str
    parents: tuple[str, ...]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(8 * 1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def required_file(path: Path, label: str) -> Path:
    if not path.is_file() or path.is_symlink() or path.stat().st_size == 0:
        raise ValueError(f"Missing, empty or symlinked {label}: {path}")
    return path


def load_json(path: Path, label: str) -> dict[str, Any]:
    required_file(path, label)
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"{label} must be a JSON object")
    return value


def parse_seqkit_stats(
    path: Path, full_genome: Path, primary_genome: Path
) -> tuple[dict[str, str], dict[str, str]]:
    required_file(path, "seqkit assembly stats")
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
    required = {"file", "num_seqs", "sum_len", "sum_n"}
    if reader.fieldnames is None or not required <= set(reader.fieldnames):
        raise ValueError("seqkit stats lacks required file/num_seqs/sum_len/sum_n fields")

    def choose(input_path: Path, label: str) -> dict[str, str]:
        matches = [row for row in rows if Path(row["file"]).name == input_path.name]
        if len(matches) != 1:
            raise ValueError(f"seqkit stats must contain exactly one {label} row")
        return matches[0]

    return choose(full_genome, "full genome"), choose(primary_genome, "primary genome")


def parse_fai(path: Path) -> dict[str, int]:
    required_file(path, "primary genome FAI")
    lengths: dict[str, int] = {}
    with path.open(encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, start=1):
            fields = raw.rstrip("\r\n").split("\t")
            if len(fields) < 2:
                raise ValueError(f"Malformed FAI row {line_number}")
            seqid = fields[0]
            try:
                length = int(fields[1])
            except ValueError as exc:
                raise ValueError(f"Non-integer FAI length at row {line_number}") from exc
            if not seqid or length <= 0 or seqid in lengths:
                raise ValueError(f"Unsafe or duplicate FAI seqid at row {line_number}")
            lengths[seqid] = length
    if not lengths:
        raise ValueError("Primary FAI is empty")
    return lengths


def parse_coding_hierarchy(
    path: Path, sequence_lengths: dict[str, int] | None
) -> dict[str, Any]:
    required_file(path, "GFF3 hierarchy input")
    genes: dict[str, Feature] = {}
    transcripts: dict[str, Feature] = {}
    duplicate_genes: set[str] = set()
    duplicate_transcripts: set[str] = set()
    cds_by_parent: dict[str, list[Feature]] = defaultdict(list)
    malformed_rows = 0
    with open_text(path) as handle:
        for raw in handle:
            line = raw.rstrip("\r\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) != 9:
                malformed_rows += 1
                continue
            seqid, _, feature_type, start_text, end_text, _, strand, phase, attr = fields
            try:
                start = int(start_text)
                end = int(end_text)
            except ValueError:
                malformed_rows += 1
                continue
            attributes, malformed = parse_attributes(attr)
            malformed_rows += malformed
            parents = tuple(value for value in attributes.get("Parent", "").split(",") if value)
            feature = Feature(seqid, start, end, strand, parents)
            if feature_type == "gene":
                identifier = attributes.get("ID", "")
                if not identifier or identifier in genes:
                    duplicate_genes.add(identifier)
                else:
                    genes[identifier] = feature
            elif feature_type in TRANSCRIPT_TYPES:
                identifier = attributes.get("ID", "")
                if not identifier or identifier in transcripts:
                    duplicate_transcripts.add(identifier)
                else:
                    transcripts[identifier] = feature
            elif feature_type == "CDS":
                if len(parents) != 1 or phase not in {"0", "1", "2"}:
                    malformed_rows += 1
                for parent in parents:
                    cds_by_parent[parent].append(feature)

    coding_transcripts = set(cds_by_parent)
    valid_transcripts: set[str] = set()
    coding_genes: set[str] = set()
    invalid_reasons: Counter[str] = Counter()
    for transcript_id in coding_transcripts:
        transcript = transcripts.get(transcript_id)
        if transcript is None or transcript_id in duplicate_transcripts:
            invalid_reasons["missing_or_duplicate_transcript"] += 1
            continue
        if len(transcript.parents) != 1:
            invalid_reasons["nonunique_gene_parent"] += 1
            continue
        gene_id = transcript.parents[0]
        gene = genes.get(gene_id)
        if gene is None or gene_id in duplicate_genes:
            invalid_reasons["missing_or_duplicate_gene"] += 1
            continue
        if (
            transcript.start < 1
            or transcript.end < transcript.start
            or gene.start < 1
            or gene.end < gene.start
            or transcript.seqid != gene.seqid
            or transcript.start < gene.start
            or transcript.end > gene.end
            or (
                transcript.strand in {"+", "-"}
                and gene.strand in {"+", "-"}
                and transcript.strand != gene.strand
            )
        ):
            invalid_reasons["invalid_transcript_gene_geometry"] += 1
            continue
        seq_length = sequence_lengths.get(transcript.seqid) if sequence_lengths else None
        if sequence_lengths is not None and seq_length is None:
            invalid_reasons["unknown_primary_seqid"] += 1
            continue
        cds_valid = True
        for cds in cds_by_parent[transcript_id]:
            if (
                cds.start < transcript.start
                or cds.end > transcript.end
                or cds.end < cds.start
                or cds.seqid != transcript.seqid
                or (seq_length is not None and cds.end > seq_length)
                or (
                    cds.strand in {"+", "-"}
                    and transcript.strand in {"+", "-"}
                    and cds.strand != transcript.strand
                )
            ):
                cds_valid = False
                break
        if not cds_valid:
            invalid_reasons["invalid_CDS_geometry"] += 1
            continue
        valid_transcripts.add(transcript_id)
        coding_genes.add(gene_id)

    denominator = len(coding_transcripts)
    hierarchy_fraction = len(valid_transcripts) / denominator if denominator else 0.0
    if malformed_rows:
        hierarchy_fraction = 0.0
    return {
        "genes": len(genes),
        "coding_transcripts": denominator,
        "valid_coding_transcripts": len(valid_transcripts),
        "coding_genes": len(coding_genes),
        "valid_coding_hierarchy_fraction": hierarchy_fraction,
        "malformed_rows_or_attributes": malformed_rows,
        "invalid_reasons": dict(sorted(invalid_reasons.items())),
    }


def busco_fraction(
    path: Path,
    expected_mode: str,
    expected_version: str,
    expected_lineage: str,
    expected_markers: int,
) -> float:
    report = load_json(path, f"{expected_mode} BUSCO summary")
    parameters = report.get("parameters")
    results = report.get("results")
    versions = report.get("versions")
    lineage = report.get("lineage_dataset")
    if not all(isinstance(value, dict) for value in (parameters, results, versions, lineage)):
        raise ValueError("BUSCO summary lacks parameters/results/versions/lineage objects")
    mode = parameters.get("mode")
    expected_observed_mode = {
        "genome": "euk_genome_min",
        "proteins": "proteins",
    }.get(expected_mode)
    if mode != expected_observed_mode:
        raise ValueError(
            f"BUSCO mode differs: expected {expected_observed_mode}, observed {mode}"
        )
    if expected_mode == "genome" and (
        parameters.get("gene_predictor") != "miniprot"
        or parameters.get("use_miniprot") != "True"
    ):
        raise ValueError("Genome BUSCO did not use the frozen Miniprot route")
    try:
        marker_count = int(lineage.get("number_of_buscos"))
    except (TypeError, ValueError) as exc:
        raise ValueError("BUSCO lineage marker count is invalid") from exc
    if (
        versions.get("busco") != expected_version
        or lineage.get("name") != expected_lineage
        or marker_count != expected_markers
        or results.get("n_markers") != expected_markers
    ):
        raise ValueError("BUSCO version, lineage, or marker count differs")
    value = results.get("Complete percentage")
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError("BUSCO Complete percentage is absent or nonnumeric")
    fraction = float(value) / 100.0
    if not 0.0 <= fraction <= 1.0:
        raise ValueError("BUSCO complete fraction is outside [0, 1]")
    return fraction


def independent_wgd_sources(path: Path) -> int:
    required_file(path, "WGD source registry")
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
    expected = [
        "source_id",
        "citation",
        "doi",
        "url",
        "evidence_scope",
        "independent_evidence_class",
        "qualifies_for_applicability",
    ]
    if reader.fieldnames != expected:
        raise ValueError("WGD source registry header differs")
    qualifying = [row for row in rows if row["qualifies_for_applicability"] == "true"]
    if any(row["qualifies_for_applicability"] not in {"true", "false"} for row in rows):
        raise ValueError("WGD source registry has a nonboolean qualification")
    for key in ("source_id", "doi", "independent_evidence_class"):
        values = [row[key] for row in qualifying]
        if any(not value for value in values) or len(values) != len(set(values)):
            raise ValueError(f"Qualifying WGD sources require unique nonempty {key}")
    return len(qualifying)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Assemble truth-blind species-applicability metrics"
    )
    parser.add_argument("--species-id", required=True)
    parser.add_argument("--declared-primary-seqid-count", required=True, type=int)
    parser.add_argument("--full-genome", required=True, type=Path)
    parser.add_argument("--primary-genome", required=True, type=Path)
    parser.add_argument("--primary-fai", required=True, type=Path)
    parser.add_argument("--seqkit-stats", required=True, type=Path)
    parser.add_argument("--full-gff", required=True, type=Path)
    parser.add_argument("--primary-gff", required=True, type=Path)
    parser.add_argument("--protein-subset-manifest", required=True, type=Path)
    parser.add_argument("--translation-quality", required=True, type=Path)
    parser.add_argument("--protein-busco-summary", required=True, type=Path)
    parser.add_argument("--genome-busco-summary", required=True, type=Path)
    parser.add_argument("--busco-lineage-root", required=True, type=Path)
    parser.add_argument("--busco-lineage-sha256sums", required=True, type=Path)
    parser.add_argument("--expected-busco-version", required=True)
    parser.add_argument("--expected-busco-lineage", required=True)
    parser.add_argument("--expected-busco-markers", required=True, type=int)
    parser.add_argument("--backbone-manifest", required=True, type=Path)
    parser.add_argument("--wgd-source-registry", required=True, type=Path)
    parser.add_argument("--implementation-source-archive", required=True, type=Path)
    parser.add_argument("--implementation-freeze-sha256sums", required=True, type=Path)
    parser.add_argument("--read-backed-qv", type=float)
    parser.add_argument("--controlled-holdout", action="store_true")
    parser.add_argument("--output-json", required=True, type=Path)
    args = parser.parse_args()
    partial_output = Path(str(args.output_json) + ".partial")
    if args.output_json.exists() or partial_output.exists():
        raise FileExistsError("Refusing to overwrite species-applicability metrics")
    if (
        not args.species_id.strip()
        or args.declared_primary_seqid_count <= 0
        or not args.expected_busco_version.strip()
        or not args.expected_busco_lineage.strip()
        or args.expected_busco_markers <= 0
    ):
        raise ValueError("Species, karyotype, and expected BUSCO fields must be positive")
    inputs = [
        args.full_genome,
        args.primary_genome,
        args.primary_fai,
        args.seqkit_stats,
        args.full_gff,
        args.primary_gff,
        args.protein_subset_manifest,
        args.translation_quality,
        args.protein_busco_summary,
        args.genome_busco_summary,
        args.busco_lineage_sha256sums,
        args.backbone_manifest,
        args.wgd_source_registry,
        args.implementation_source_archive,
        args.implementation_freeze_sha256sums,
    ]
    for path in inputs:
        required_file(path, "applicability input")
    if not args.busco_lineage_root.is_dir() or args.busco_lineage_root.is_symlink():
        raise ValueError("BUSCO lineage root is missing or symlinked")
    if args.busco_lineage_sha256sums.parent.resolve() != args.busco_lineage_root.resolve():
        raise ValueError("BUSCO lineage SHA256SUMS must be inside the lineage root")
    verify_sha256sums(
        args.busco_lineage_root,
        manifest_path=args.busco_lineage_sha256sums,
        ignore_checksum_file=True,
    )
    implementation_entries = read_sha256sums(args.implementation_freeze_sha256sums)
    if (
        len(implementation_entries) != 1
        or implementation_entries[0].relative_path.as_posix() != args.implementation_source_archive.name
        or implementation_entries[0].sha256 != sha256(args.implementation_source_archive)
    ):
        raise ValueError("Implementation freeze manifest does not bind the source archive")

    full_stats, primary_stats = parse_seqkit_stats(
        args.seqkit_stats, args.full_genome, args.primary_genome
    )
    full_length = int(full_stats["sum_len"])
    primary_length = int(primary_stats["sum_len"])
    primary_n = int(primary_stats["sum_n"])
    primary_seqids = int(primary_stats["num_seqs"])
    if full_length <= 0 or primary_length <= 0 or primary_n < 0 or primary_n > primary_length:
        raise ValueError("Invalid seqkit length/N statistics")
    fai = parse_fai(args.primary_fai)
    if sum(fai.values()) != primary_length:
        raise ValueError("Primary FAI length sum differs from seqkit stats")

    full_hierarchy = parse_coding_hierarchy(args.full_gff, None)
    primary_hierarchy = parse_coding_hierarchy(args.primary_gff, fai)
    if full_hierarchy["coding_genes"] <= 0:
        raise ValueError("Full GFF has no valid protein-coding genes")
    primary_coding_gene_fraction = (
        primary_hierarchy["coding_genes"] / full_hierarchy["coding_genes"]
    )
    if not 0.0 <= primary_coding_gene_fraction <= 1.0:
        raise ValueError("Primary protein-coding genes exceed the full annotation universe")

    protein = load_json(args.protein_subset_manifest, "protein subset manifest")
    translation = load_json(args.translation_quality, "translation-quality report")
    backbone = load_json(args.backbone_manifest, "fixed-backbone manifest")
    if protein.get("schema_version") != PROTEIN_SUBSET_SCHEMA:
        raise ValueError("Protein-subset schema differs")
    if translation.get("schema_version") != TRANSLATION_SCHEMA:
        raise ValueError("Translation-quality schema differs")
    if backbone.get("schema_version") != BACKBONE_SCHEMA:
        raise ValueError("Fixed-backbone schema differs")
    if protein.get("candidate_access") is not False or protein.get("truth_access") is not False:
        raise ValueError("Protein quality report is not truth/candidate blind")
    if protein.get("mapping_policy") != (
        "exact_transcript_relation_ID_after_repeated_recognized_namespace_stripping"
    ):
        raise ValueError("Protein mapping policy is not the frozen exact-ID adapter")
    if translation.get("candidate_access") is not False or translation.get("truth_access") is not False:
        raise ValueError("Translation report is not truth/candidate blind")
    if translation.get("comparison_policy") != (
        "exact_amino_acid_identity_after_terminal_stop_removal_and_"
        "repeated_recognized_namespace_stripping"
    ):
        raise ValueError("Translation comparison policy differs")
    if backbone.get("candidate_access") is not False or backbone.get("truth_access") is not False:
        raise ValueError("Backbone is not truth/candidate blind")
    invariants = backbone.get("semantic_invariants")
    identifier_mapping = backbone.get("identifier_mapping")
    expected_invariants = {
        "input_permutation_deterministic": True,
        "reverse_duplicate_statistics_audit_only": True,
        "ambiguous_endpoint_block_quarantine": True,
    }
    if not isinstance(invariants, dict) or any(
        invariants.get(key) is not value for key, value in expected_invariants.items()
    ):
        raise ValueError("Fixed-backbone semantic invariants are absent or differ")
    if (
        not isinstance(identifier_mapping, dict)
        or identifier_mapping.get("mapping_policy")
        != "unique_exact_attribute_then_unique_normalized_feature_ID_fallback"
        or identifier_mapping.get("bidirectional_unique") is not True
    ):
        raise ValueError("Fixed-backbone identifier mapping is not exact and unique")
    backbone_metrics = backbone.get("applicability_metrics")
    policy = backbone.get("policy")
    if not isinstance(backbone_metrics, dict) or not isinstance(policy, dict):
        raise ValueError("Fixed-backbone metrics/policy are absent")
    if policy.get("require_cross_seqid") is not True:
        raise ValueError("Fixed backbone is not cross-primary-seqid only")
    minimum_block_pairs = policy.get("min_block_pairs")
    if isinstance(minimum_block_pairs, bool) or not isinstance(minimum_block_pairs, int):
        raise ValueError("Fixed-backbone minimum block pairs is invalid")

    mapping_fraction = protein.get("fractions", {}).get(
        "exact_unique_GFF_protein_mapping_fraction"
    )
    translation_fraction = translation.get("fractions", {}).get(
        "valid_representative_translation_fraction"
    )
    if not isinstance(mapping_fraction, (int, float)) or isinstance(mapping_fraction, bool):
        raise ValueError("Protein mapping fraction is invalid")
    if not isinstance(translation_fraction, (int, float)) or isinstance(
        translation_fraction, bool
    ):
        raise ValueError("Translation fraction is invalid")

    controlled_flag = bool(args.controlled_holdout)
    report = {
        "schema_version": SCHEMA_VERSION,
        "species_id": args.species_id,
        "controlled_holdout": controlled_flag,
        "assembly": {
            "primary_seqid_count_matches_declared_karyotype": (
                primary_seqids == args.declared_primary_seqid_count
                and len(fai) == args.declared_primary_seqid_count
            ),
            "primary_assembly_fraction": primary_length / full_length,
            "primary_non_N_fraction": (primary_length - primary_n) / primary_length,
            "assembly_BUSCO_complete_fraction": busco_fraction(
                args.genome_busco_summary,
                "genome",
                args.expected_busco_version,
                args.expected_busco_lineage,
                args.expected_busco_markers,
            ),
            "read_backed_QV": args.read_backed_qv,
        },
        "annotation": {
            "primary_protein_coding_gene_fraction": primary_coding_gene_fraction,
            "exact_unique_GFF_protein_mapping_fraction": float(mapping_fraction),
            "valid_coding_hierarchy_fraction": primary_hierarchy[
                "valid_coding_hierarchy_fraction"
            ],
            "valid_representative_translation_fraction": float(translation_fraction),
            "protein_BUSCO_complete_fraction": busco_fraction(
                args.protein_busco_summary,
                "proteins",
                args.expected_busco_version,
                args.expected_busco_lineage,
                args.expected_busco_markers,
            ),
            "fuzzy_identifier_repairs": 0,
        },
        "backbone": {
            "independent_WGD_source_count": independent_wgd_sources(
                args.wgd_source_registry
            ),
            "minimum_block_pairs": minimum_block_pairs,
            "cross_primary_seqid_only": True,
            "primary_gene_midpoint_backbone_coverage": backbone_metrics[
                "primary_gene_midpoint_backbone_coverage"
            ],
            "primary_chromosome_cell_coverage_fraction": backbone_metrics[
                "primary_chromosome_cell_coverage_fraction"
            ],
            "minimum_cells_per_covered_chromosome_observed": backbone_metrics[
                "minimum_cells_per_covered_chromosome_observed"
            ],
            "unique_partner_chromosome_fraction_among_covered_genes": backbone_metrics[
                "unique_partner_chromosome_fraction_among_covered_genes"
            ],
            "input_permutation_deterministic": True,
            "reverse_duplicate_invariant": True,
            "built_without_candidates": True,
            "built_without_truth_or_labels": True,
            "used_perturbed_annotation_if_controlled": (
                backbone.get("base_annotation_role") == "perturbed"
                if controlled_flag
                else False
            ),
        },
        "source_manifests": {
            path.name + ":" + str(index): sha256(path)
            for index, path in enumerate(inputs)
        },
        "audit": {
            "full_coding_hierarchy": full_hierarchy,
            "primary_coding_hierarchy": primary_hierarchy,
            "backbone_applicable": backbone.get("applicable"),
            "no_candidate_statistics_used": True,
            "no_truth_or_labels_used": True,
        },
    }
    for section_name in ("assembly", "annotation", "backbone"):
        for field, value in report[section_name].items():
            if isinstance(value, float) and not math.isfinite(value):
                raise ValueError(f"Non-finite applicability metric: {section_name}.{field}")
    serialized = json.dumps(report, indent=2, sort_keys=True) + "\n"
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    if args.output_json.parent.is_symlink():
        raise ValueError("Applicability output parent must not be a symlink")
    with partial_output.open("x", encoding="utf-8", newline="") as handle:
        handle.write(serialized)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(partial_output, args.output_json)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
