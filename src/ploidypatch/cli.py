from __future__ import annotations

import argparse
import json
from collections.abc import Sequence
from pathlib import Path

from . import __version__
from .audit import audit_bundle, write_json
from .baseline import (
    adapt_annotation_gff_baseline,
    adapt_miniprot_baseline,
    prepare_reference_proteins,
    summarize_projection_support,
)
from .bootstrap import (
    EVENT_BOOLEAN_METRICS,
    independent_event_bootstrap,
    paired_event_bootstrap,
)
from .catalog import write_missing_gene_candidate_catalog
from .catalog_summary import write_candidate_catalog_summary
from .calibration import calibrate_synteny_tiers
from .candidate_merge import merge_candidate_gffs
from .confusion_bootstrap import (
    SUPPORTED_CONFUSION_SECTIONS,
    independent_confusion_bootstrap,
)
from .consensus import select_method_consensus
from .conflict_guard import apply_conflict_winner_guard
from .copy_features import build_copy_candidate_features, label_copy_candidate_features
from .copy_model import (
    score_copy_candidate_features,
    select_scored_copy_candidates,
)
from .copy_model_eval import evaluate_copy_candidate_scores
from .copy_pair_sampling import sample_balanced_copy_pairs
from .copy_patch import (
    compile_copy_addition_patch_edits,
    compile_reviewed_copy_addition_patch_edits,
)
from .event_graph import infer_event_graph
from .fixed_backbone_projection import (
    FixedBackboneProjectionPolicy,
    project_candidates_to_fixed_backbone,
)
from .fixed_target_backbone import (
    FixedTargetBackbonePolicy,
    build_fixed_target_backbone,
)
from .homeolog_pairs import (
    infer_outgroup_duplicated_pairs,
    infer_wgdi_homeolog_pairs,
)
from .homeolog_ranker import (
    freeze_homeolog_review_rankings,
    score_homeolog_copy_candidates,
)
from .homeolog_topology import build_homeolog_topology_features
from .support_ranker import score_support_conditioned_candidates
from .isoseq_validation import (
    bootstrap_isoseq_review_yield,
    filter_candidate_query_paf,
    join_isoseq_review_rankings,
    prepare_b73_isoseq_transcripts,
    validate_isoseq_candidate_chains,
)
from .localization import score_synteny_localization, write_synteny_model_labels
from .normalize import (
    normalize_provider_gff3,
    prepare_ncbi_primary_bundle,
    prepare_primary_annotation_bundle,
    read_primary_seqid_table,
)
from .pair_consensus import intersect_copy_pair_evidence
from .negative import (
    audit_masked_gap_genome,
    create_masked_gap_control,
    score_masked_gap_abstention,
    summarize_masked_gap_selection,
)
from .natural import (
    annotate_natural_assembly_context,
    discover_natural_candidates,
    prepare_natural_graph_inputs,
    summarize_natural_validation,
    validate_natural_candidates_with_rna,
    validate_natural_candidates_with_secondary_groups,
)
from .natural_audit import audit_natural_candidates, export_natural_candidate_cds
from .pav import (
    annotate_pav_with_wgdi,
    catalog_genes_in_target_deletions,
    collapse_high_confidence_pav_events,
    combine_chromosome_pafs,
    extract_paf_target_deletions,
    reconcile_pav_gene_catalogs,
    screen_pav_candidate_proteins,
    subset_catalog_proteins,
)
from .patch import (
    apply_annotation_patch,
    create_annotation_patch,
    revert_annotation_patch,
)
from .projection_select import select_projection_support_models
from .rna import aggregate_junctions, extract_bam_junctions
from .sampling import sample_candidate_catalog
from .self_wgd_pairs import infer_self_wgdi_pairs
from .reference_anchored import aggregate_reference_anchored_projection
from .review_report import build_review_report
from .species_applicability import evaluate_species_applicability
from .perturb import (
    BOUNDARY_SHIFT_EVENT,
    COPY_COLLAPSE_EVENT,
    FUSED_GENE_EVENT,
    MISSING_GENE_EVENT,
    MISSING_INTERNAL_EXON_EVENT,
    SPLIT_GENE_EVENT,
    generate_missing_gene_benchmark,
    restore_gff_from_truth,
)
from .score import score_annotation_repair
from .structure_candidate import adapt_miniprot_structure_candidates
from .structure_bootstrap import paired_structure_hypothesis_bootstrap
from .structure_hypothesis import infer_structure_hypotheses
from .structure_hypothesis_score import score_structure_hypotheses
from .structure_patch import SUPPORTED_PATCH_EVENTS, compile_structure_patch_edits
from .structure_perturb import generate_structure_benchmark
from .synteny_io import prepare_wgdi_inputs
from .synteny_gap import infer_wgdi_synteny_gaps, select_synteny_gap_models
from .wgdi_summary import summarize_wgdi_gene_evidence
from .wgd_candidate_select import (
    propagate_wgd_selection_to_conflict_pool,
    select_wgd_supported_candidates,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="ploidypatch")
    parser.add_argument("--version", action="version", version=__version__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    audit_parser = subparsers.add_parser(
        "audit", help="Audit one assembly annotation bundle before inference"
    )
    audit_parser.add_argument("--gff", required=True, help="GFF3 or GFF3.gz path")
    audit_parser.add_argument("--protein", required=True, help="protein FASTA path")
    audit_parser.add_argument("--cds", required=True, help="CDS FASTA path")
    audit_parser.add_argument("--fai", help="optional FASTA index for coordinate checks")
    audit_parser.add_argument("--output", help="write JSON report to this path")
    audit_parser.add_argument(
        "--checksums", action="store_true", help="calculate SHA-256 for audited inputs"
    )

    report_parser = subparsers.add_parser(
        "report",
        help="Build a self-contained review report from a candidate-pool run",
    )
    report_parser.add_argument(
        "--candidate-gff", required=True,
        help="chain-preserving consensus GFF3 containing consensus_digest genes",
    )
    report_parser.add_argument(
        "--pool-decisions", required=True,
        help="candidate-pool decisions TSV bound by the pool manifest",
    )
    report_parser.add_argument(
        "--pool-manifest", required=True,
        help="ploidypatch.method_candidate_pool.v2 JSON manifest",
    )
    report_parser.add_argument(
        "--output-dir", required=True,
        help="new non-overwriting directory for HTML, JSON, TSV, and checksums",
    )
    report_parser.add_argument(
        "--review-decisions",
        help="optional explicit human review ledger TSV",
    )
    report_parser.add_argument(
        "--scores",
        help="optional exact-universe rank-score TSV",
    )
    report_parser.add_argument(
        "--copy-features",
        help="optional exact-universe copy-feature TSV for method overlap",
    )
    report_parser.add_argument(
        "--topology-features",
        help="optional exact-universe topology TSV for applicability context",
    )
    report_parser.add_argument(
        "--patch-edits",
        help="optional reviewed patch-edits JSON for accept/event consistency",
    )
    report_parser.add_argument(
        "--run-summary",
        help="optional run summary proving byte-identical reversion",
    )
    report_parser.add_argument(
        "--title", default="Annotation repair review",
        help="human-readable project or sample title",
    )
    report_parser.add_argument(
        "--max-embedded-candidates",
        type=int,
        default=5000,
        help=(
            "maximum ranked candidates embedded in the interactive HTML; "
            "report.json and candidates.tsv always retain the full universe"
        ),
    )
    report_parser.add_argument(
        "--fail-on-attention",
        action="store_true",
        help="exit with status 2 when the generated report requires attention",
    )

    benchmark_parser = subparsers.add_parser(
        "benchmark", help="Generate or validate controlled annotation perturbations"
    )
    benchmark_subparsers = benchmark_parser.add_subparsers(
        dest="benchmark_command", required=True
    )

    graph_parser = subparsers.add_parser(
        "graph", help="Build and score a plant event evidence graph"
    )
    graph_subparsers = graph_parser.add_subparsers(
        dest="graph_command", required=True
    )
    graph_infer_parser = graph_subparsers.add_parser(
        "infer", help="Apply transparent evidence scores and hard constraints"
    )
    graph_infer_parser.add_argument("--candidates", required=True)
    graph_infer_parser.add_argument("--evidence", required=True)
    graph_infer_parser.add_argument("--output-json", required=True)
    graph_infer_parser.add_argument("--decisions-tsv", required=True)

    patch_parser = subparsers.add_parser(
        "patch", help="Create, apply, or revert an immutable annotation patch"
    )
    patch_subparsers = patch_parser.add_subparsers(
        dest="patch_command", required=True
    )
    patch_create_parser = patch_subparsers.add_parser("create")
    patch_create_parser.add_argument("--source-gff", required=True)
    patch_create_parser.add_argument("--edits-json", required=True)
    patch_create_parser.add_argument("--output-patch", required=True)
    patch_compile_parser = patch_subparsers.add_parser(
        "compile-structure",
        help="Compile explicitly allowed exact-topology hypotheses to line edits",
    )
    patch_compile_parser.add_argument("--annotation-gff", required=True)
    patch_compile_parser.add_argument("--candidate-gff", required=True)
    patch_compile_parser.add_argument("--hypotheses-tsv", required=True)
    patch_compile_parser.add_argument("--output-edits-json", required=True)
    patch_compile_parser.add_argument(
        "--event-type",
        action="append",
        required=True,
        choices=sorted(SUPPORTED_PATCH_EVENTS),
    )
    patch_compile_parser.add_argument(
        "--min-support-group-count", type=int, default=2
    )
    patch_copy_compile_parser = patch_subparsers.add_parser(
        "compile-copy-additions",
        help="Compile a prefiltered candidate GFF into reversible EOF edits",
    )
    patch_copy_compile_parser.add_argument("--annotation-gff", required=True)
    patch_copy_compile_parser.add_argument("--candidate-gff", required=True)
    patch_copy_compile_parser.add_argument("--output-edits-json", required=True)
    patch_reviewed_copy_compile_parser = patch_subparsers.add_parser(
        "compile-reviewed-copy-additions",
        help="Compile only explicit human-accepted candidates from a frozen pool",
    )
    patch_reviewed_copy_compile_parser.add_argument("--annotation-gff", required=True)
    patch_reviewed_copy_compile_parser.add_argument("--candidate-gff", required=True)
    patch_reviewed_copy_compile_parser.add_argument("--pool-decisions", required=True)
    patch_reviewed_copy_compile_parser.add_argument("--pool-manifest", required=True)
    patch_reviewed_copy_compile_parser.add_argument("--review-decisions", required=True)
    patch_reviewed_copy_compile_parser.add_argument(
        "--output-edits-json", required=True
    )
    patch_apply_parser = patch_subparsers.add_parser("apply")
    patch_apply_parser.add_argument("--source-gff", required=True)
    patch_apply_parser.add_argument("--patch", required=True)
    patch_apply_parser.add_argument("--output-gff", required=True)
    patch_revert_parser = patch_subparsers.add_parser("revert")
    patch_revert_parser.add_argument("--patched-gff", required=True)
    patch_revert_parser.add_argument("--patch", required=True)
    patch_revert_parser.add_argument("--output-gff", required=True)
    perturb_parser = benchmark_subparsers.add_parser(
        "perturb", help="Create a deterministic hidden-truth benchmark"
    )
    perturb_parser.add_argument("--gff", required=True, help="source GFF3 or GFF3.gz")
    perturb_parser.add_argument(
        "--output-dir", required=True, help="new or empty benchmark output directory"
    )
    perturb_parser.add_argument(
        "--truth-dir",
        help="separate evaluator-only directory for hidden_truth.json",
    )
    perturb_parser.add_argument(
        "--event-type",
        choices=[
            MISSING_GENE_EVENT,
            MISSING_INTERNAL_EXON_EVENT,
            BOUNDARY_SHIFT_EVENT,
            SPLIT_GENE_EVENT,
            FUSED_GENE_EVENT,
            COPY_COLLAPSE_EVENT,
        ],
        default=MISSING_GENE_EVENT,
    )
    perturb_selection = perturb_parser.add_mutually_exclusive_group(required=True)
    perturb_selection.add_argument("--count", type=int)
    perturb_selection.add_argument(
        "--selection-tsv",
        help="evaluator-only catalog sample with an adjacent verified manifest",
    )
    perturb_parser.add_argument("--seed", type=int, default=1)
    perturb_parser.add_argument(
        "--pair-tsv",
        help=(
            "evaluator-only gene_id_a/gene_id_b table required for "
            "annotation_copy_collapse"
        ),
    )

    restore_parser = benchmark_subparsers.add_parser(
        "restore", help="Invert a perturbation using evaluator-only hidden truth"
    )
    restore_parser.add_argument("--perturbed-gff", required=True)
    restore_parser.add_argument("--truth", required=True)
    restore_parser.add_argument("--output-gff", required=True)

    score_parser = benchmark_subparsers.add_parser(
        "score", help="Score a candidate repaired GFF3 against evaluator-only truth"
    )
    score_parser.add_argument("--source-gff", required=True)
    score_parser.add_argument("--perturbed-gff", required=True)
    score_parser.add_argument("--candidate-gff", required=True)
    score_parser.add_argument("--truth", required=True)
    score_parser.add_argument(
        "--include-event-details",
        action="store_true",
        help="include target IDs; keep disabled for aggregate-only tuning",
    )
    hypothesis_score_parser = benchmark_subparsers.add_parser(
        "score-structure-hypotheses",
        help="Score typed audit hypotheses separately from applied repairs",
    )
    hypothesis_score_parser.add_argument("--source-gff", required=True)
    hypothesis_score_parser.add_argument("--perturbed-gff", required=True)
    hypothesis_score_parser.add_argument("--candidate-gff", required=True)
    hypothesis_score_parser.add_argument("--hypotheses-tsv", required=True)
    hypothesis_score_parser.add_argument("--truth", required=True)
    hypothesis_score_parser.add_argument("--control-hypotheses-tsv")
    hypothesis_score_parser.add_argument(
        "--include-event-details", action="store_true"
    )
    bootstrap_parser = benchmark_subparsers.add_parser(
        "bootstrap-events",
        help="Paired within-event-type bootstrap for complete recovery metrics",
    )
    bootstrap_parser.add_argument(
        "--score",
        action="append",
        required=True,
        help="repeatable LABEL=score.json input",
    )
    bootstrap_parser.add_argument("--output-json", required=True)
    bootstrap_parser.add_argument(
        "--metric",
        choices=sorted(EVENT_BOOLEAN_METRICS),
        default="complete_cds_chain_recovery",
    )
    bootstrap_parser.add_argument("--replicates", type=int, default=10_000)
    bootstrap_parser.add_argument("--seed", type=int, default=20260807)
    bootstrap_parser.add_argument("--alpha", type=float, default=0.05)
    independent_bootstrap_parser = benchmark_subparsers.add_parser(
        "bootstrap-independent-events",
        help="Independent event bootstrap across different benchmark sets",
    )
    independent_bootstrap_parser.add_argument(
        "--score",
        action="append",
        required=True,
        help="repeatable LABEL=score.json input",
    )
    independent_bootstrap_parser.add_argument("--output-json", required=True)
    independent_bootstrap_parser.add_argument(
        "--metric",
        choices=sorted(EVENT_BOOLEAN_METRICS),
        default="complete_cds_chain_recovery",
    )
    independent_bootstrap_parser.add_argument(
        "--replicates", type=int, default=10_000
    )
    independent_bootstrap_parser.add_argument("--seed", type=int, default=20260807)
    independent_bootstrap_parser.add_argument("--alpha", type=float, default=0.05)
    confusion_bootstrap_parser = benchmark_subparsers.add_parser(
        "bootstrap-confusion",
        help="Independent TP/FP/FN bootstrap for precision, recall, and F1",
    )
    confusion_bootstrap_parser.add_argument(
        "--score",
        action="append",
        required=True,
        help="repeatable LABEL=score.json input",
    )
    confusion_bootstrap_parser.add_argument("--output-json", required=True)
    confusion_bootstrap_parser.add_argument(
        "--section",
        choices=sorted(SUPPORTED_CONFUSION_SECTIONS),
        default="strict_cds_chain",
    )
    confusion_bootstrap_parser.add_argument("--replicates", type=int, default=10_000)
    confusion_bootstrap_parser.add_argument("--seed", type=int, default=20260807)
    confusion_bootstrap_parser.add_argument("--alpha", type=float, default=0.05)
    copy_label_parser = benchmark_subparsers.add_parser(
        "label-copy-features",
        help="Attach evaluator-only exact-CDS labels to blind copy features",
    )
    copy_label_parser.add_argument("--features", required=True)
    copy_label_parser.add_argument("--truth", required=True)
    copy_label_parser.add_argument("--output-tsv", required=True)
    copy_ranking_score_parser = benchmark_subparsers.add_parser(
        "score-copy-ranking",
        help="Evaluate frozen truth-blind copy scores against evaluator-only labels",
    )
    copy_ranking_score_parser.add_argument("--scores", required=True)
    copy_ranking_score_parser.add_argument("--labeled-features", required=True)
    copy_ranking_score_parser.add_argument("--output-json", required=True)
    copy_pair_sample_parser = benchmark_subparsers.add_parser(
        "sample-copy-pairs",
        help=(
            "Select evaluator copy-collapse pairs with deterministic chromosome "
            "and CDS-complexity balance"
        ),
    )
    copy_pair_sample_parser.add_argument("--source-gff", required=True)
    copy_pair_sample_parser.add_argument("--pairs", required=True)
    copy_pair_sample_parser.add_argument("--count", type=int, required=True)
    copy_pair_sample_parser.add_argument("--seed", type=int, required=True)
    copy_pair_sample_parser.add_argument("--output-pairs", required=True)
    copy_pair_sample_parser.add_argument("--decisions-tsv", required=True)
    structure_bootstrap_parser = benchmark_subparsers.add_parser(
        "bootstrap-structure-hypotheses",
        help="Paired event bootstrap for typed structure-hypothesis tiers",
    )
    structure_bootstrap_parser.add_argument(
        "--score",
        action="append",
        required=True,
        help="repeatable LABEL=score.json; a label may span event files",
    )
    structure_bootstrap_parser.add_argument("--output-json", required=True)
    structure_bootstrap_parser.add_argument(
        "--replicates", type=int, default=10_000
    )
    structure_bootstrap_parser.add_argument("--seed", type=int, default=20260807)
    structure_bootstrap_parser.add_argument("--alpha", type=float, default=0.05)
    score_parser.add_argument(
        "--control-candidate-gff",
        help=(
            "optional same-method candidate generated from the complete annotation; "
            "subtracts background proposals in controlled perturbation evaluation"
        ),
    )
    score_parser.add_argument(
        "--event-strata",
        help="evaluator-only TSV keyed by hidden gene_id",
    )
    score_parser.add_argument(
        "--stratum-column",
        action="append",
        default=[],
        help="repeatable marginal/joint recall column from --event-strata",
    )
    localization_score_parser = benchmark_subparsers.add_parser(
        "score-localization",
        help="Score blind synteny gaps and model spans against hidden loci",
    )
    localization_score_parser.add_argument(
        "--gaps", action="append", required=True, help="repeatable upstream gap TSV"
    )
    localization_score_parser.add_argument("--selection", required=True)
    localization_score_parser.add_argument("--truth", required=True)
    localization_score_parser.add_argument(
        "--include-event-details", action="store_true"
    )
    label_models_parser = benchmark_subparsers.add_parser(
        "label-synteny-models",
        help="Write evaluator-only per-model localization and CDS labels",
    )
    label_models_parser.add_argument("--source-gff", required=True)
    label_models_parser.add_argument("--candidate-gff", required=True)
    label_models_parser.add_argument("--selection", required=True)
    label_models_parser.add_argument("--baseline-decisions", required=True)
    label_models_parser.add_argument("--truth", required=True)
    label_models_parser.add_argument("--output-tsv", required=True)
    label_models_parser.add_argument("--control-candidate-gff")
    calibrate_parser = benchmark_subparsers.add_parser(
        "calibrate-synteny-tiers",
        help="Fit interpretable evidence tiers on chromosome-disjoint development data",
    )
    calibrate_parser.add_argument("--labels", required=True)
    calibrate_parser.add_argument("--output", required=True)
    calibrate_parser.add_argument(
        "--label-column", default="label_exact_cds_chain"
    )
    calibrate_parser.add_argument(
        "--eligible-column",
        help="optional binary column restricting rows eligible for calibration",
    )
    calibrate_parser.add_argument(
        "--feature-set",
        choices=("synteny", "projection_support"),
        default="synteny",
    )
    calibrate_parser.add_argument("--high-precision-floor", type=float, default=0.1)
    calibrate_parser.add_argument(
        "--high-precision-min-selected", type=int, default=20
    )

    catalog_parser = benchmark_subparsers.add_parser(
        "catalog", help="Create an evaluator-only catalog for stratified sampling"
    )
    catalog_parser.add_argument("--gff", required=True)
    catalog_parser.add_argument("--output-tsv", required=True)
    catalog_parser.add_argument(
        "--external-strata",
        help="optional TSV keyed by exact gene_id, for WGD/subgenome labels",
    )
    catalog_parser.add_argument(
        "--external-strata-prefix",
        default="",
        help="prefix added to every imported strata column",
    )
    catalog_parser.add_argument(
        "--primary-chromosomes-only",
        action="store_true",
        help="retain NCBI region features marked genome=chromosome",
    )
    catalog_summary_parser = benchmark_subparsers.add_parser(
        "summarize-catalog", help="Count auditable one-way and joint candidate strata"
    )
    catalog_summary_parser.add_argument("--catalog", required=True)
    catalog_summary_parser.add_argument("--output", required=True)
    catalog_summary_parser.add_argument(
        "--column",
        action="append",
        required=True,
        help="repeatable column to summarize; empty values are counted explicitly",
    )
    catalog_summary_parser.add_argument(
        "--cross",
        action="append",
        default=[],
        help="repeatable COLUMN_A,COLUMN_B joint count",
    )
    catalog_sample_parser = benchmark_subparsers.add_parser(
        "sample-catalog", help="Take a deterministic sample from declared strata"
    )
    catalog_sample_parser.add_argument("--catalog", required=True)
    catalog_sample_parser.add_argument("--plan", required=True)
    catalog_sample_parser.add_argument("--output-tsv", required=True)
    catalog_sample_parser.add_argument("--seed", type=int, required=True)
    catalog_sample_parser.add_argument(
        "--exclude-selection",
        action="append",
        default=[],
        help="repeatable prior selection TSV whose candidate IDs are excluded",
    )
    mask_parser = benchmark_subparsers.add_parser(
        "mask-genome", help="Create a paired assembly-gap negative control"
    )
    mask_parser.add_argument("--genome", required=True)
    mask_parser.add_argument("--truth", required=True)
    mask_parser.add_argument("--output-genome", required=True)
    mask_parser.add_argument("--output-mask-truth", required=True)
    mask_parser.add_argument(
        "--background-gff",
        help=(
            "exclude target spans overlapping retained exon/CDS features in this GFF3"
        ),
    )
    abstention_parser = benchmark_subparsers.add_parser(
        "score-abstention", help="Score false repairs at evaluator-only masked loci"
    )
    abstention_parser.add_argument("--perturbed-gff", required=True)
    abstention_parser.add_argument("--candidate-gff", required=True)
    abstention_parser.add_argument("--mask-truth", required=True)
    abstention_parser.add_argument(
        "--include-event-details",
        action="store_true",
        help="include evaluator-only event identifiers and coordinates",
    )
    mask_audit_parser = benchmark_subparsers.add_parser(
        "audit-mask", help="Verify a masked genome against evaluator-only mask truth"
    )
    mask_audit_parser.add_argument("--source-genome", required=True)
    mask_audit_parser.add_argument("--masked-genome", required=True)
    mask_audit_parser.add_argument("--mask-truth", required=True)
    mask_summary_parser = benchmark_subparsers.add_parser(
        "summarize-mask-selection",
        help="Summarize clean masked-control retention by evaluator strata",
    )
    mask_summary_parser.add_argument("--mask-truth", required=True)
    mask_summary_parser.add_argument("--hidden-truth", required=True)
    mask_summary_parser.add_argument("--strata-tsv", required=True)
    mask_summary_parser.add_argument("--column", action="append", required=True)

    evidence_parser = subparsers.add_parser(
        "evidence", help="Prepare auditable inputs for external evidence engines"
    )
    evidence_subparsers = evidence_parser.add_subparsers(
        dest="evidence_command", required=True
    )
    wgdi_parser = evidence_subparsers.add_parser(
        "prepare-wgdi", help="Normalize gene coordinates and representative proteins"
    )
    wgdi_parser.add_argument("--gff", required=True)
    wgdi_parser.add_argument("--protein", required=True)
    wgdi_parser.add_argument("--fai", required=True)
    wgdi_parser.add_argument("--output-dir", required=True)
    wgdi_parser.add_argument("--prefix", required=True)
    wgdi_parser.add_argument("--min-genes-per-seqid", type=int, default=1)
    wgdi_parser.add_argument(
        "--primary-chromosomes-only",
        action="store_true",
        help="retain NCBI region features marked genome=chromosome",
    )
    wgdi_summary_parser = evidence_subparsers.add_parser(
        "summarize-wgdi", help="Summarize WGDI blocks into per-gene plant strata"
    )
    wgdi_summary_parser.add_argument("--query-gff", required=True)
    wgdi_summary_parser.add_argument("--query-input-manifest", required=True)
    wgdi_summary_parser.add_argument(
        "--collinearity",
        action="append",
        required=True,
        help="repeatable SOURCE=PATH value",
    )
    wgdi_summary_parser.add_argument(
        "--expected-source",
        action="append",
        required=True,
        help="repeatable SUBGENOME=SOURCE mapping",
    )
    wgdi_summary_parser.add_argument("--output-gene-tsv", required=True)
    wgdi_summary_parser.add_argument("--output-block-tsv", required=True)
    homeolog_pair_parser = evidence_subparsers.add_parser(
        "infer-homeolog-pairs",
        help="Infer conservative named-WGD pairs from independent WGDI anchors",
    )
    homeolog_pair_parser.add_argument("--gene-evidence", required=True)
    homeolog_pair_parser.add_argument(
        "--collinearity",
        action="append",
        required=True,
        help="repeatable SOURCE=PATH value",
    )
    homeolog_pair_parser.add_argument("--source-group-map")
    homeolog_pair_parser.add_argument("--wgd-event", required=True)
    homeolog_pair_parser.add_argument(
        "--subgenome",
        action="append",
        required=True,
        help="repeat exactly twice in the desired A/B output order",
    )
    homeolog_pair_parser.add_argument(
        "--min-support-group-count", type=int, default=2
    )
    homeolog_pair_parser.add_argument(
        "--allow-multiple-partners",
        action="store_true",
        help="disable the default reciprocal-unique partner gate",
    )
    homeolog_pair_parser.add_argument("--output-pairs", required=True)
    homeolog_pair_parser.add_argument("--decisions-tsv", required=True)
    outgroup_duplicate_parser = evidence_subparsers.add_parser(
        "infer-outgroup-duplicated-pairs",
        help=(
            "Infer duplicated query descendants supported one-to-two by "
            "independent outgroups"
        ),
    )
    outgroup_duplicate_parser.add_argument("--query-wgdi-gff", required=True)
    outgroup_duplicate_parser.add_argument(
        "--source-gff",
        help="map normalized WGDI IDs uniquely to operable source gene feature IDs",
    )
    outgroup_duplicate_parser.add_argument(
        "--collinearity",
        action="append",
        required=True,
        help="repeatable SOURCE=PATH value",
    )
    outgroup_duplicate_parser.add_argument("--source-group-map")
    outgroup_duplicate_parser.add_argument("--wgd-event", required=True)
    outgroup_duplicate_parser.add_argument(
        "--min-support-group-count", type=int, default=2
    )
    outgroup_duplicate_parser.add_argument("--min-block-pairs", type=int, default=20)
    outgroup_duplicate_parser.add_argument(
        "--allow-same-seqid", action="store_true"
    )
    outgroup_duplicate_parser.add_argument(
        "--allow-multiple-partners", action="store_true"
    )
    outgroup_duplicate_parser.add_argument("--output-pairs", required=True)
    outgroup_duplicate_parser.add_argument("--decisions-tsv", required=True)
    self_wgd_pair_parser = evidence_subparsers.add_parser(
        "infer-self-wgd-pairs",
        help="Infer cross-chromosome reciprocal pairs from self-WGDI blocks",
    )
    self_wgd_pair_parser.add_argument("--query-wgdi-gff", required=True)
    self_wgd_pair_parser.add_argument("--collinearity", required=True)
    self_wgd_pair_parser.add_argument(
        "--source-gff",
        help="map WGDI IDs uniquely to operable source gene feature IDs",
    )
    self_wgd_pair_parser.add_argument("--wgd-event", required=True)
    self_wgd_pair_parser.add_argument("--min-block-pairs", type=int, default=20)
    self_wgd_pair_parser.add_argument(
        "--allow-same-seqid", action="store_true"
    )
    self_wgd_pair_parser.add_argument(
        "--allow-multiple-partners", action="store_true"
    )
    self_wgd_pair_parser.add_argument("--output-pairs", required=True)
    self_wgd_pair_parser.add_argument("--decisions-tsv", required=True)
    reference_anchored_parser = evidence_subparsers.add_parser(
        "aggregate-reference-anchored",
        help=(
            "Aggregate candidate-source provenance, accepted reference-WGD "
            "pairs, and source-homeolog-to-base evidence"
        ),
    )
    reference_anchored_parser.add_argument(
        "--candidate-source-provenance", required=True
    )
    reference_anchored_parser.add_argument(
        "--accepted-wgd-pairs",
        action="append",
        required=True,
        help="repeatable REFERENCE=PATH accepted directed reciprocal pair table",
    )
    reference_anchored_parser.add_argument(
        "--homeolog-base-evidence",
        action="append",
        required=True,
        help="repeatable REFERENCE=PATH arm/homeolog/base evidence table",
    )
    reference_anchored_parser.add_argument("--output-dir", required=True)
    reference_anchored_parser.add_argument(
        "--evidence-type", default="reference_anchored_projection"
    )
    pair_intersection_parser = evidence_subparsers.add_parser(
        "intersect-copy-pair-evidence",
        help="Intersect exact unordered duplicate pairs from independent evidence",
    )
    pair_intersection_parser.add_argument(
        "--pairs",
        action="append",
        required=True,
        help="repeatable LABEL=PATH pair table",
    )
    pair_intersection_parser.add_argument("--pair-set-label", required=True)
    pair_intersection_parser.add_argument(
        "--allow-multiple-partners", action="store_true"
    )
    pair_intersection_parser.add_argument("--output-pairs", required=True)
    pair_intersection_parser.add_argument("--decisions-tsv", required=True)
    synteny_gap_parser = evidence_subparsers.add_parser(
        "infer-synteny-gaps",
        help="Infer missing-gene hypotheses between adjacent WGDI anchors",
    )
    synteny_gap_parser.add_argument("--query-wgdi-gff", required=True)
    synteny_gap_parser.add_argument("--target-wgdi-gff", required=True)
    synteny_gap_parser.add_argument("--collinearity", required=True)
    synteny_gap_parser.add_argument("--source-label", required=True)
    synteny_gap_parser.add_argument("--output-tsv", required=True)
    synteny_gap_parser.add_argument(
        "--expected-chromosome-pairs",
        help="optional query/target/source TSV for subgenome-aware filtering",
    )
    synteny_gap_parser.add_argument(
        "--max-query-intervening-genes", type=int, default=0
    )
    synteny_gap_parser.add_argument(
        "--min-target-excess-genes", type=int, default=1
    )
    synteny_gap_parser.add_argument("--max-target-gap-genes", type=int, default=5)
    synteny_gap_parser.add_argument("--max-query-locus-bp", type=int, default=500000)
    synteny_select_parser = evidence_subparsers.add_parser(
        "select-synteny-gap-models",
        help="Retain accepted protein models that fall inside blind synteny gaps",
    )
    synteny_select_parser.add_argument(
        "--gaps",
        action="append",
        required=True,
        help="repeatable filtered synteny-gap TSV",
    )
    synteny_select_parser.add_argument("--baseline-decisions", required=True)
    synteny_select_parser.add_argument("--adapted-candidate-gff", required=True)
    synteny_select_parser.add_argument("--output-selection", required=True)
    synteny_select_parser.add_argument("--output-candidate-gff", required=True)
    projection_support_parser = evidence_subparsers.add_parser(
        "summarize-projection-support",
        help="Recover multi-query and multi-progenitor support after deduplication",
    )
    projection_support_parser.add_argument("--decisions", required=True)
    projection_support_parser.add_argument("--output-tsv", required=True)
    copy_feature_parser = evidence_subparsers.add_parser(
        "build-copy-features",
        help="Build truth-blind copy-candidate features across method families",
    )
    copy_feature_parser.add_argument("--consensus-decisions", required=True)
    copy_feature_parser.add_argument(
        "--method-decisions",
        action="append",
        required=True,
        help="repeatable METHOD=PATH decision table",
    )
    copy_feature_parser.add_argument("--wgd-selection")
    copy_feature_parser.add_argument("--output-tsv", required=True)
    topology_feature_parser = evidence_subparsers.add_parser(
        "build-homeolog-topology-features",
        help="Compare candidate CDS topology with an existing WGD partner",
    )
    topology_feature_parser.add_argument("--copy-features", required=True)
    topology_feature_parser.add_argument("--wgd-selection", required=True)
    topology_feature_parser.add_argument("--candidate-gff", required=True)
    topology_feature_parser.add_argument("--base-gff", required=True)
    topology_feature_parser.add_argument("--output-tsv", required=True)
    fixed_backbone_parser = evidence_subparsers.add_parser(
        "build-fixed-target-backbone",
        help="Build a candidate-free target self-synteny backbone for stable topology",
    )
    fixed_backbone_parser.add_argument("--base-gff", required=True)
    fixed_backbone_parser.add_argument("--query-wgdi-gff", required=True)
    fixed_backbone_parser.add_argument("--base-only-collinearity", required=True)
    fixed_backbone_parser.add_argument("--primary-seqid-table", required=True)
    fixed_backbone_parser.add_argument("--min-block-pairs", type=int, default=20)
    fixed_backbone_parser.add_argument("--output-dir", required=True)
    fixed_projection_parser = evidence_subparsers.add_parser(
        "project-fixed-backbone",
        help="Project each candidate independently onto a fixed target backbone",
    )
    fixed_projection_parser.add_argument("--backbone-dir", required=True)
    fixed_projection_parser.add_argument("--candidate-catalog", required=True)
    fixed_projection_parser.add_argument("--fixed-base-hits", required=True)
    fixed_projection_parser.add_argument("--fixed-base-hits-manifest", required=True)
    fixed_projection_parser.add_argument(
        "--minimum-query-coverage", type=float, default=0.5
    )
    fixed_projection_parser.add_argument("--maximum-evalue", type=float, default=1e-5)
    fixed_projection_parser.add_argument("--output-dir", required=True)
    applicability_parser = evidence_subparsers.add_parser(
        "evaluate-species-applicability",
        help="Apply frozen pre-candidate assembly, annotation and backbone gates",
    )
    applicability_parser.add_argument("--metrics", required=True)
    applicability_parser.add_argument("--policy", required=True)
    applicability_parser.add_argument("--output-json", required=True)
    homeolog_rank_parser = evidence_subparsers.add_parser(
        "score-homeolog-copy-candidates",
        help=(
            "Rank copy candidates with frozen plant homeolog topology and "
            "explicit automatic-repair abstention"
        ),
    )
    homeolog_rank_parser.add_argument("--copy-features", required=True)
    homeolog_rank_parser.add_argument("--topology-features", required=True)
    homeolog_rank_parser.add_argument("--model-json", required=True)
    homeolog_rank_parser.add_argument("--output-tsv", required=True)
    support_rank_parser = evidence_subparsers.add_parser(
        "score-support-conditioned-candidates",
        help="Apply the frozen v0.3 support-conditioned review ranker",
    )
    support_rank_parser.add_argument("--copy-features", required=True)
    support_rank_parser.add_argument("--topology-features", required=True)
    support_rank_parser.add_argument("--model-json", required=True)
    support_rank_parser.add_argument("--output-tsv", required=True)
    conflict_guard_parser = evidence_subparsers.add_parser(
        "apply-conflict-winner-guard",
        help="Apply the label-free v0.4 baseline-winner safety fallback",
    )
    conflict_guard_parser.add_argument("--v03-scores", required=True)
    conflict_guard_parser.add_argument("--pool-decisions", required=True)
    conflict_guard_parser.add_argument("--pool-manifest", required=True)
    conflict_guard_parser.add_argument("--output-tsv", required=True)
    homeolog_review_parser = evidence_subparsers.add_parser(
        "freeze-homeolog-review-rankings",
        help="Freeze deterministic truth-blind baseline and topology review queues",
    )
    homeolog_review_parser.add_argument("--scores", required=True)
    homeolog_review_parser.add_argument(
        "--review-budget",
        action="append",
        type=int,
        help="repeatable top-K review budget; defaults to 25,50,100,200",
    )
    homeolog_review_parser.add_argument("--output-tsv", required=True)
    isoseq_prepare_parser = evidence_subparsers.add_parser(
        "prepare-b73-isoseq",
        help="Select full-length transcripts using only pure-B73 read counts",
    )
    isoseq_prepare_parser.add_argument("--fasta", required=True)
    isoseq_prepare_parser.add_argument("--counts-csv", required=True)
    isoseq_prepare_parser.add_argument("--minimum-b73-reads", type=int, default=2)
    isoseq_prepare_parser.add_argument("--output-fasta", required=True)
    isoseq_prepare_parser.add_argument("--output-counts", required=True)
    isoseq_filter_parser = evidence_subparsers.add_parser(
        "filter-candidate-query-paf",
        help=(
            "Retain all alignments for reads whose strand-aware mapping can "
            "overlap a frozen candidate"
        ),
    )
    isoseq_filter_parser.add_argument("--candidate-gff", required=True)
    isoseq_filter_parser.add_argument(
        "--paf-input",
        action="append",
        required=True,
        help="repeatable ACCESSION=PATH PAF input",
    )
    isoseq_filter_parser.add_argument(
        "--alignment-strand-source",
        choices=("query_orientation", "minimap2_ts"),
        default="query_orientation",
    )
    isoseq_filter_parser.add_argument("--output-paf", required=True)
    isoseq_filter_parser.add_argument("--output-counts", required=True)
    isoseq_filter_parser.add_argument("--output-summary", required=True)
    isoseq_filter_parser.add_argument("--output-manifest", required=True)
    isoseq_validate_parser = evidence_subparsers.add_parser(
        "validate-isoseq-candidates",
        help="Validate frozen CDS chains with independently mapped full-length RNA",
    )
    isoseq_validate_parser.add_argument("--candidate-gff", required=True)
    isoseq_validate_parser.add_argument("--paf", required=True)
    isoseq_validate_parser.add_argument("--selected-counts", required=True)
    isoseq_validate_parser.add_argument("--genome-fasta", required=True)
    isoseq_validate_parser.add_argument("--output-evidence", required=True)
    isoseq_validate_parser.add_argument("--minimum-query-coverage", type=float, default=0.90)
    isoseq_validate_parser.add_argument("--minimum-identity", type=float, default=0.98)
    isoseq_validate_parser.add_argument("--minimum-mapq", type=int, default=20)
    isoseq_validate_parser.add_argument("--maximum-secondary-score-fraction", type=float, default=0.95)
    isoseq_validate_parser.add_argument("--minimum-candidate-cds-coverage", type=float, default=0.90)
    isoseq_validate_parser.add_argument("--flank-bp", type=int, default=5000)
    isoseq_validate_parser.add_argument(
        "--alignment-strand-source",
        choices=("query_orientation", "minimap2_ts"),
        default="query_orientation",
    )
    isoseq_join_parser = evidence_subparsers.add_parser(
        "join-isoseq-review-rankings",
        help="Attach frozen candidate ranks to full-length RNA evidence",
    )
    isoseq_join_parser.add_argument("--evidence", required=True)
    isoseq_join_parser.add_argument("--review-rankings", required=True)
    isoseq_join_parser.add_argument("--review-budget", action="append", type=int)
    isoseq_join_parser.add_argument("--comparator-estimator", default="baseline")
    isoseq_join_parser.add_argument("--primary-estimator", default="topology")
    isoseq_join_parser.add_argument("--output-tsv", required=True)
    isoseq_join_parser.add_argument("--output-summary", required=True)
    isoseq_bootstrap_parser = evidence_subparsers.add_parser(
        "bootstrap-isoseq-review-yield",
        help="Chromosome-aware top-K Iso-Seq enrichment and paired ranker bootstrap",
    )
    isoseq_bootstrap_parser.add_argument("--evidence", required=True)
    isoseq_bootstrap_parser.add_argument("--review-rankings", required=True)
    isoseq_bootstrap_parser.add_argument("--review-budget", action="append", type=int)
    isoseq_bootstrap_parser.add_argument("--replicates", type=int, default=20_000)
    isoseq_bootstrap_parser.add_argument("--seed", type=int, default=20260808)
    isoseq_bootstrap_parser.add_argument("--alpha", type=float, default=0.05)
    isoseq_bootstrap_parser.add_argument("--comparator-estimator", default="baseline")
    isoseq_bootstrap_parser.add_argument("--primary-estimator", default="topology")
    isoseq_bootstrap_parser.add_argument("--output-json", required=True)
    natural_cds_parser = evidence_subparsers.add_parser(
        "export-natural-candidate-cds",
        help="Export deterministic candidate CDS queries for genome self-mapping",
    )
    natural_cds_parser.add_argument("--candidate-gff", required=True)
    natural_cds_parser.add_argument("--genome-fasta", required=True)
    natural_cds_parser.add_argument("--output-fasta", required=True)
    natural_audit_parser = evidence_subparsers.add_parser(
        "audit-natural-candidates",
        help="Audit candidate ORFs, annotation collisions, copy number, and RNA support",
    )
    natural_audit_parser.add_argument("--candidate-gff", required=True)
    natural_audit_parser.add_argument("--base-gff", required=True)
    natural_audit_parser.add_argument("--genome-fasta", required=True)
    natural_audit_parser.add_argument("--review-rankings", required=True)
    natural_audit_parser.add_argument("--isoseq-evidence")
    natural_audit_parser.add_argument("--self-map-paf")
    natural_audit_parser.add_argument("--repeat-gff")
    natural_audit_parser.add_argument("--repeat-seqid-map")
    natural_audit_parser.add_argument("--repeat-flank-bp", type=int, default=2000)
    natural_audit_parser.add_argument(
        "--minimum-full-length-read-support", type=int, default=1
    )
    natural_audit_parser.add_argument("--review-budget", action="append", type=int)
    natural_audit_parser.add_argument(
        "--minimum-query-coverage", type=float, default=0.90
    )
    natural_audit_parser.add_argument("--minimum-identity", type=float, default=0.98)
    natural_audit_parser.add_argument(
        "--near-equal-score-fraction", type=float, default=0.95
    )
    natural_audit_parser.add_argument("--output-tsv", required=True)
    natural_audit_parser.add_argument("--output-summary", required=True)
    copy_score_parser = evidence_subparsers.add_parser(
        "score-copy-candidates",
        help="Score truth-blind copy candidates with a frozen portable model",
    )
    copy_score_parser.add_argument("--features", required=True)
    copy_score_parser.add_argument("--model-json", required=True)
    copy_score_parser.add_argument(
        "--mask-feature-group",
        action="append",
        choices=("wgd_context", "method_quality"),
        default=[],
        help=(
            "repeatable counterfactual evidence mask; applies the same frozen "
            "model and thresholds without refitting"
        ),
    )
    copy_score_parser.add_argument("--output-tsv", required=True)
    copy_model_select_parser = evidence_subparsers.add_parser(
        "select-scored-copy-candidates",
        help="Retain complete candidate hierarchies passing a frozen model policy",
    )
    copy_model_select_parser.add_argument("--base-gff", required=True)
    copy_model_select_parser.add_argument("--candidate-gff", required=True)
    copy_model_select_parser.add_argument("--scores", required=True)
    copy_model_select_parser.add_argument("--model-json", required=True)
    copy_model_select_parser.add_argument(
        "--policy", choices=("review", "high_confidence"), default="review"
    )
    copy_model_select_parser.add_argument("--output-gff", required=True)
    copy_model_select_parser.add_argument("--selection-tsv", required=True)
    projection_select_parser = evidence_subparsers.add_parser(
        "select-projection-support",
        help="Retain appended protein models supported by independent groups",
    )
    projection_select_parser.add_argument("--candidate-gff", required=True)
    projection_select_parser.add_argument("--projection-support", required=True)
    projection_select_parser.add_argument("--source-group-map")
    projection_select_parser.add_argument(
        "--min-support-group-count", type=int, default=2
    )
    projection_select_parser.add_argument("--output-gff", required=True)
    projection_select_parser.add_argument("--selection-tsv", required=True)
    wgd_candidate_select_parser = evidence_subparsers.add_parser(
        "select-wgd-supported-candidates",
        help="Retain candidates paired to an existing annotated WGD partner",
    )
    wgd_candidate_select_parser.add_argument("--base-gff", required=True)
    wgd_candidate_select_parser.add_argument("--candidate-gff", required=True)
    wgd_candidate_select_parser.add_argument("--pairs", required=True)
    wgd_candidate_select_parser.add_argument("--output-gff", required=True)
    wgd_candidate_select_parser.add_argument("--selection-tsv", required=True)
    wgd_propagate_parser = evidence_subparsers.add_parser(
        "propagate-wgd-conflict-partners",
        help="Propagate one unique locus-level WGD partner to retained chain alternatives",
    )
    wgd_propagate_parser.add_argument("--base-gff", required=True)
    wgd_propagate_parser.add_argument("--candidate-gff", required=True)
    wgd_propagate_parser.add_argument("--pool-decisions", required=True)
    wgd_propagate_parser.add_argument("--prior-wgd-selection", required=True)
    wgd_propagate_parser.add_argument("--output-selection", required=True)
    junction_extract_parser = evidence_subparsers.add_parser(
        "extract-bam-junctions",
        help="Extract auditable unstranded splice-junction counts from one BAM",
    )
    junction_extract_parser.add_argument("--bam", required=True)
    junction_extract_parser.add_argument("--sample", required=True)
    junction_extract_parser.add_argument("--samtools", required=True)
    junction_extract_parser.add_argument("--output-tsv", required=True)
    junction_extract_parser.add_argument("--threads", type=int, default=4)
    junction_extract_parser.add_argument("--min-mapq", type=int, default=20)
    junction_aggregate_parser = evidence_subparsers.add_parser(
        "aggregate-junctions",
        help="Aggregate primary and secondary RNA junction evidence tiers",
    )
    junction_aggregate_parser.add_argument(
        "--input-dir",
        action="append",
        required=True,
        help="repeatable directory containing *.junctions.tsv plus manifests",
    )
    junction_aggregate_parser.add_argument(
        "--primary-sample", action="append", required=True
    )
    junction_aggregate_parser.add_argument("--output-tsv", required=True)
    junction_aggregate_parser.add_argument(
        "--min-reads-per-sample", type=int, default=2
    )
    junction_aggregate_parser.add_argument(
        "--min-supporting-samples", type=int, default=2
    )
    junction_aggregate_parser.add_argument(
        "--sample-groups",
        help="optional explicit sample/filename_stem_group TSV",
    )
    junction_aggregate_parser.add_argument(
        "--min-samples-per-group", type=int, default=2
    )
    junction_aggregate_parser.add_argument(
        "--min-secondary-groups", type=int, default=2
    )
    natural_discovery_parser = evidence_subparsers.add_parser(
        "discover-natural",
        help="Discover RNA-blind plant structure hypotheses from projections",
    )
    natural_discovery_parser.add_argument("--target-gff", required=True)
    natural_discovery_parser.add_argument("--miniprot-gff", required=True)
    natural_discovery_parser.add_argument("--protein-map", required=True)
    natural_discovery_parser.add_argument("--output-tsv", required=True)
    natural_discovery_parser.add_argument("--min-identity", type=float, default=0.8)
    natural_discovery_parser.add_argument(
        "--min-query-coverage", type=float, default=0.8
    )
    natural_discovery_parser.add_argument(
        "--max-existing-cds-overlap", type=float, default=0.1
    )
    natural_discovery_parser.add_argument(
        "--min-boundary-extension-bp", type=int, default=30
    )
    natural_discovery_parser.add_argument(
        "--near-best-score-fraction", type=float, default=0.95
    )
    natural_rna_parser = evidence_subparsers.add_parser(
        "validate-natural-rna",
        help="Validate a frozen RNA-blind catalog against held-out junctions",
    )
    natural_rna_parser.add_argument("--candidates", required=True)
    natural_rna_parser.add_argument("--junctions", required=True)
    natural_rna_parser.add_argument("--output-tsv", required=True)
    natural_secondary_parser = evidence_subparsers.add_parser(
        "validate-natural-secondary-rna",
        help="Add conservative filename-stem-group recurrence to frozen candidates",
    )
    natural_secondary_parser.add_argument("--primary-validation", required=True)
    natural_secondary_parser.add_argument("--grouped-junctions", required=True)
    natural_secondary_parser.add_argument("--output-tsv", required=True)
    natural_assembly_parser = evidence_subparsers.add_parser(
        "annotate-natural-assembly",
        help="Add indexed assembly-edge, ambiguity, and soft-mask context",
    )
    natural_assembly_parser.add_argument("--validation-tsv", required=True)
    natural_assembly_parser.add_argument("--genome", required=True)
    natural_assembly_parser.add_argument("--fai", required=True)
    natural_assembly_parser.add_argument("--output-tsv", required=True)
    natural_assembly_parser.add_argument("--flank-bp", type=int, default=5000)
    natural_assembly_parser.add_argument(
        "--max-ambiguous-fraction", type=float, default=0.0
    )
    natural_summary_parser = evidence_subparsers.add_parser(
        "summarize-natural",
        help="Summarize natural candidate event-by-evidence cross counts",
    )
    natural_summary_parser.add_argument("--validation-tsv", required=True)
    natural_summary_parser.add_argument("--output-json", required=True)
    natural_graph_parser = evidence_subparsers.add_parser(
        "prepare-natural-graph",
        help="Adapt held-out natural evidence into correlation-aware graph inputs",
    )
    natural_graph_parser.add_argument("--validation-tsv", required=True)
    natural_graph_parser.add_argument("--output-candidates", required=True)
    natural_graph_parser.add_argument("--output-evidence", required=True)
    combine_paf_parser = evidence_subparsers.add_parser(
        "combine-chromosome-pafs",
        help="Validate and combine per-query PAFs in chromosome-pair order",
    )
    combine_paf_parser.add_argument("--input-dir", required=True)
    combine_paf_parser.add_argument("--chromosome-pairs", required=True)
    combine_paf_parser.add_argument("--output-paf", required=True)
    paf_deletion_parser = evidence_subparsers.add_parser(
        "extract-paf-deletions",
        help="Extract well-anchored target-only intervals from PAF cg CIGARs",
    )
    paf_deletion_parser.add_argument("--paf", required=True)
    paf_deletion_parser.add_argument("--output-tsv", required=True)
    paf_deletion_parser.add_argument("--min-deletion-bp", type=int, default=100)
    paf_deletion_parser.add_argument(
        "--min-flanking-anchor-bp", type=int, default=1000
    )
    paf_deletion_parser.add_argument("--min-mapq", type=int, default=20)
    paf_deletion_parser.add_argument(
        "--chromosome-pairs",
        help="optional query_seqid/target_seqid TSV for homologous-pair filtering",
    )
    paf_deletion_parser.add_argument(
        "--query-genome",
        help="uncompressed query FASTA for assembly-ambiguity checks at breakpoints",
    )
    paf_deletion_parser.add_argument(
        "--query-fai", help="samtools FAI for --query-genome"
    )
    paf_deletion_parser.add_argument("--query-flank-bp", type=int, default=5000)
    paf_deletion_parser.add_argument(
        "--max-query-flank-ambiguous-fraction", type=float, default=0.0
    )
    paf_deletion_parser.add_argument(
        "--retain-uncertain-query-flanks",
        action="store_true",
        help="report rather than exclude truncated or ambiguous query flanks",
    )
    paf_deletion_parser.add_argument(
        "--allow-secondary", action="store_true", help="retain non-primary PAF rows"
    )
    pav_gene_parser = evidence_subparsers.add_parser(
        "catalog-pav-genes",
        help="Catalog coding genes contained in one anchored target deletion",
    )
    pav_gene_parser.add_argument("--gff", required=True)
    pav_gene_parser.add_argument("--deletions", required=True)
    pav_gene_parser.add_argument("--output-tsv", required=True)
    pav_gene_parser.add_argument("--min-gene-coverage", type=float, default=0.95)
    pav_gene_parser.add_argument("--min-cds-coverage", type=float, default=1.0)
    pav_reconcile_parser = evidence_subparsers.add_parser(
        "reconcile-pav-catalogs",
        help="Retain gene-level PAV candidates supported by multiple runs",
    )
    pav_reconcile_parser.add_argument(
        "--catalog",
        action="append",
        required=True,
        help="repeatable LABEL=PATH catalog input",
    )
    pav_reconcile_parser.add_argument("--output-tsv", required=True)
    pav_reconcile_parser.add_argument("--min-support", type=int, default=2)
    pav_protein_subset_parser = evidence_subparsers.add_parser(
        "subset-pav-proteins",
        help="Extract gene-keyed representative proteins for a PAV catalog",
    )
    pav_protein_subset_parser.add_argument("--catalog", required=True)
    pav_protein_subset_parser.add_argument("--proteins", required=True)
    pav_protein_subset_parser.add_argument("--output-fasta", required=True)
    pav_protein_screen_parser = evidence_subparsers.add_parser(
        "screen-pav-proteins",
        help="Classify full-assembly protein evidence for PAV candidates",
    )
    pav_protein_screen_parser.add_argument("--catalog", required=True)
    pav_protein_screen_parser.add_argument("--miniprot-gff", required=True)
    pav_protein_screen_parser.add_argument("--protein-map", required=True)
    pav_protein_screen_parser.add_argument("--chromosome-pairs", required=True)
    pav_protein_screen_parser.add_argument("--output-tsv", required=True)
    pav_protein_screen_parser.add_argument("--min-identity", type=float, default=0.5)
    pav_protein_screen_parser.add_argument(
        "--min-query-coverage", type=float, default=0.5
    )
    pav_protein_screen_parser.add_argument(
        "--min-partial-query-coverage", type=float, default=0.2
    )
    pav_protein_screen_parser.add_argument(
        "--max-local-distance-bp", type=int, default=100000
    )
    pav_protein_screen_parser.add_argument(
        "--max-breakpoint-spread-bp", type=int, default=1000
    )
    pav_wgdi_parser = evidence_subparsers.add_parser(
        "annotate-pav-wgdi",
        help="Combine locus-aware PAV screening with plant WGD evidence",
    )
    pav_wgdi_parser.add_argument("--protein-screen", required=True)
    pav_wgdi_parser.add_argument("--wgdi-genes", required=True)
    pav_wgdi_parser.add_argument("--output-tsv", required=True)
    pav_event_parser = evidence_subparsers.add_parser(
        "collapse-pav-events",
        help="Collapse high-confidence PAV genes into structural events",
    )
    pav_event_parser.add_argument("--annotated-genes", required=True)
    pav_event_parser.add_argument("--output-tsv", required=True)

    normalize_parser = subparsers.add_parser(
        "normalize", help="Create auditable, internally consistent input bundles"
    )
    normalize_subparsers = normalize_parser.add_subparsers(
        dest="normalize_command", required=True
    )
    ncbi_primary_parser = normalize_subparsers.add_parser(
        "ncbi-primary", help="Subset an NCBI bundle to clean primary chromosomes"
    )
    ncbi_primary_parser.add_argument("--gff", required=True)
    ncbi_primary_parser.add_argument("--protein", required=True)
    ncbi_primary_parser.add_argument("--cds", required=True)
    ncbi_primary_parser.add_argument("--genome", required=True)
    ncbi_primary_parser.add_argument("--output-dir", required=True)
    primary_annotation_parser = normalize_subparsers.add_parser(
        "primary-annotation",
        help="Subset genome and GFF to primary chromosomes without filtering models",
    )
    primary_annotation_parser.add_argument("--gff", required=True)
    primary_annotation_parser.add_argument("--genome", required=True)
    primary_annotation_parser.add_argument("--output-dir", required=True)
    primary_annotation_parser.add_argument(
        "--primary-seqid-table",
        help=(
            "optional strict TSV with seqid and chromosome_label columns; "
            "use for provider GFF files without chromosome declaration records"
        ),
    )
    provider_gff_parser = normalize_subparsers.add_parser(
        "provider-gff3",
        help="Apply narrow, explicitly enabled compatibility fixes to provider GFF3",
    )
    provider_gff_parser.add_argument("--gff", required=True)
    provider_gff_parser.add_argument("--output-dir", required=True)
    provider_gff_parser.add_argument(
        "--repair-unescaped-note-semicolons",
        action="store_true",
        help="percent-escape bare semicolon continuations only inside Note values",
    )
    provider_gff_parser.add_argument(
        "--drop-invalid-intron-intervals",
        action="store_true",
        help="drop negative-length intron records without changing exon or CDS records",
    )
    provider_gff_parser.add_argument(
        "--strip-embedded-fasta",
        action="store_true",
        help="stop at ##FASTA when a separate genome FASTA is supplied downstream",
    )

    baseline_parser = subparsers.add_parser(
        "baseline", help="Prepare or adapt declared external baseline evidence"
    )
    baseline_subparsers = baseline_parser.add_subparsers(
        dest="baseline_command", required=True
    )
    reference_parser = baseline_subparsers.add_parser(
        "prepare-proteins", help="Prefix and combine external protein references"
    )
    reference_parser.add_argument(
        "--protein",
        action="append",
        required=True,
        help="repeatable SOURCE=PATH input",
    )
    reference_parser.add_argument("--output-fasta", required=True)
    reference_parser.add_argument("--output-map", required=True)
    merge_candidate_parser = baseline_subparsers.add_parser(
        "merge-candidate-gffs",
        help="Namespace and merge multiple references from one method family",
    )
    merge_candidate_parser.add_argument(
        "--candidate",
        action="append",
        required=True,
        help="repeatable REFERENCE_SOURCE=GFF3 input",
    )
    merge_candidate_parser.add_argument("--output-gff", required=True)
    merge_candidate_parser.add_argument("--provenance-tsv", required=True)
    miniprot_adapter_parser = baseline_subparsers.add_parser(
        "adapt-miniprot", help="Append conservative gap-only miniprot projections"
    )
    miniprot_adapter_parser.add_argument("--perturbed-gff", required=True)
    miniprot_adapter_parser.add_argument("--miniprot-gff", required=True)
    miniprot_adapter_parser.add_argument("--protein-map", required=True)
    miniprot_adapter_parser.add_argument("--output-gff", required=True)
    miniprot_adapter_parser.add_argument("--decisions-tsv", required=True)
    miniprot_adapter_parser.add_argument("--min-identity", type=float, default=0.5)
    miniprot_adapter_parser.add_argument(
        "--min-query-coverage", type=float, default=0.5
    )
    miniprot_adapter_parser.add_argument(
        "--max-existing-cds-overlap", type=float, default=0.2
    )
    miniprot_adapter_parser.add_argument(
        "--max-redundancy-overlap", type=float, default=0.5
    )
    miniprot_adapter_parser.add_argument(
        "--allow-disrupted",
        action="store_true",
        help="retain alignments with frameshifts or in-frame stops",
    )
    gff_adapter_parser = baseline_subparsers.add_parser(
        "adapt-gff",
        help="Append coding hypotheses from a general GFF3 in annotation gaps",
    )
    gff_adapter_parser.add_argument("--perturbed-gff", required=True)
    gff_adapter_parser.add_argument("--candidate-gff", required=True)
    gff_adapter_parser.add_argument("--source", required=True)
    gff_adapter_parser.add_argument("--output-gff", required=True)
    gff_adapter_parser.add_argument("--decisions-tsv", required=True)
    gff_adapter_parser.add_argument("--score-attribute")
    gff_adapter_parser.add_argument(
        "--max-existing-cds-overlap", type=float, default=0.2
    )
    gff_adapter_parser.add_argument(
        "--max-redundancy-overlap", type=float, default=0.5
    )
    gff_adapter_parser.add_argument(
        "--infer-missing-cds-phase",
        action="store_true",
        help=(
            "infer phases only when every CDS phase in a model is missing; "
            "assumes a complete CDS whose first phase is zero"
        ),
    )
    consensus_parser = baseline_subparsers.add_parser(
        "select-method-consensus",
        help="Select exact phased-CDS chains supported by method families",
    )
    consensus_parser.add_argument("--base-gff", required=True)
    consensus_parser.add_argument(
        "--candidate",
        action="append",
        required=True,
        help="repeatable METHOD=adapted_candidate.gff3 input",
    )
    consensus_parser.add_argument("--output-gff", required=True)
    consensus_parser.add_argument("--decisions-tsv", required=True)
    consensus_parser.add_argument("--min-method-support", type=int, default=2)
    consensus_parser.add_argument(
        "--max-redundancy-overlap", type=float, default=0.5
    )
    consensus_parser.add_argument(
        "--redundancy-policy",
        choices=("suppress_overlapping", "retain_distinct_chains"),
        default="suppress_overlapping",
        help=(
            "suppress strongly overlapping chains (v1 default), or retain "
            "every distinct phased CDS chain and label conflict sets"
        ),
    )
    structure_adapter_parser = baseline_subparsers.add_parser(
        "adapt-miniprot-structure",
        help="Append exact-chain projection disagreements as structure candidates",
    )
    structure_adapter_parser.add_argument("--annotation-gff", required=True)
    structure_adapter_parser.add_argument("--miniprot-gff", required=True)
    structure_adapter_parser.add_argument("--protein-map", required=True)
    structure_adapter_parser.add_argument(
        "--source-group-map",
        help="optional TSV mapping each protein source to an independent support group",
    )
    structure_adapter_parser.add_argument("--output-gff", required=True)
    structure_adapter_parser.add_argument("--decisions-tsv", required=True)
    structure_adapter_parser.add_argument("--min-identity", type=float, default=0.5)
    structure_adapter_parser.add_argument(
        "--min-query-coverage", type=float, default=0.5
    )
    structure_adapter_parser.add_argument(
        "--min-source-support", type=int, default=1
    )
    structure_adapter_parser.add_argument(
        "--min-gene-overlap-fraction", type=float, default=0.1
    )
    structure_adapter_parser.add_argument(
        "--allow-disrupted",
        action="store_true",
        help="retain alignments with frameshifts or in-frame stops",
    )
    structure_hypothesis_parser = baseline_subparsers.add_parser(
        "infer-structure-hypotheses",
        help="Infer conservative error types from exact CDS-chain topology",
    )
    structure_hypothesis_parser.add_argument("--annotation-gff", required=True)
    structure_hypothesis_parser.add_argument("--candidate-gff", required=True)
    structure_hypothesis_parser.add_argument("--output-tsv", required=True)
    structure_hypothesis_parser.add_argument(
        "--candidate-topology-tsv", required=True
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command == "audit":
        report = audit_bundle(
            gff_path=args.gff,
            protein_path=args.protein,
            cds_path=args.cds,
            fai_path=args.fai,
            checksums=args.checksums,
        )
        write_json(report, args.output)
        return 0 if report["quality_gate"]["grade"] != "fail" else 1
    if args.command == "report":
        report = build_review_report(
            candidate_gff_path=args.candidate_gff,
            pool_decisions_tsv_path=args.pool_decisions,
            pool_manifest_json_path=args.pool_manifest,
            output_dir=args.output_dir,
            review_decisions_tsv_path=args.review_decisions,
            scores_tsv_path=args.scores,
            copy_features_tsv_path=args.copy_features,
            topology_features_tsv_path=args.topology_features,
            patch_edits_json_path=args.patch_edits,
            run_summary_json_path=args.run_summary,
            title=args.title,
            max_embedded_candidates=args.max_embedded_candidates,
        )
        print(
            json.dumps(
                {
                    "output_dir": str(Path(args.output_dir).resolve()),
                    "state": report["summary"]["state"],
                    "counts": report["counts"],
                },
                indent=2,
                sort_keys=True,
            )
        )
        if args.fail_on_attention and report["summary"]["state"] == "attention_required":
            return 2
        return 0
    if args.command == "graph" and args.graph_command == "infer":
        graph = infer_event_graph(
            candidate_tsv_path=args.candidates,
            evidence_tsv_path=args.evidence,
            output_json_path=args.output_json,
            decisions_tsv_path=args.decisions_tsv,
        )
        print(json.dumps(graph["counts"], indent=2, sort_keys=True))
        return 0
    if args.command == "patch" and args.patch_command == "create":
        patch = create_annotation_patch(
            source_gff_path=args.source_gff,
            edits_json_path=args.edits_json,
            output_patch_path=args.output_patch,
        )
        print(json.dumps(patch, indent=2, sort_keys=True))
        return 0
    if args.command == "patch" and args.patch_command == "compile-copy-additions":
        report = compile_copy_addition_patch_edits(
            annotation_gff_path=args.annotation_gff,
            candidate_gff_path=args.candidate_gff,
            output_edits_json_path=args.output_edits_json,
        )
        print(json.dumps(report, indent=2, sort_keys=True))
        return 0
    if (
        args.command == "patch"
        and args.patch_command == "compile-reviewed-copy-additions"
    ):
        report = compile_reviewed_copy_addition_patch_edits(
            annotation_gff_path=args.annotation_gff,
            candidate_gff_path=args.candidate_gff,
            pool_decisions_tsv_path=args.pool_decisions,
            pool_manifest_json_path=args.pool_manifest,
            review_decisions_tsv_path=args.review_decisions,
            output_edits_json_path=args.output_edits_json,
        )
        print(json.dumps(report, indent=2, sort_keys=True))
        return 0
    if args.command == "patch" and args.patch_command == "compile-structure":
        edits = compile_structure_patch_edits(
            annotation_gff_path=args.annotation_gff,
            candidate_gff_path=args.candidate_gff,
            hypotheses_tsv_path=args.hypotheses_tsv,
            output_edits_json_path=args.output_edits_json,
            allowed_event_types=tuple(args.event_type),
            min_support_group_count=args.min_support_group_count,
        )
        print(json.dumps(edits["counts"], indent=2, sort_keys=True))
        return 0
    if args.command == "patch" and args.patch_command == "apply":
        report = apply_annotation_patch(
            source_gff_path=args.source_gff,
            patch_path=args.patch,
            output_gff_path=args.output_gff,
        )
        print(json.dumps(report, indent=2, sort_keys=True))
        return 0
    if args.command == "patch" and args.patch_command == "revert":
        report = revert_annotation_patch(
            patched_gff_path=args.patched_gff,
            patch_path=args.patch,
            output_gff_path=args.output_gff,
        )
        print(json.dumps(report, indent=2, sort_keys=True))
        return 0
    if args.command == "benchmark" and args.benchmark_command == "perturb":
        if args.event_type == MISSING_GENE_EVENT:
            if args.pair_tsv is not None:
                raise ValueError("--pair-tsv is only valid for annotation_copy_collapse")
            manifest = generate_missing_gene_benchmark(
                gff_path=args.gff,
                output_dir=args.output_dir,
                count=args.count,
                seed=args.seed,
                selection_tsv_path=args.selection_tsv,
                truth_output_dir=args.truth_dir,
            )
        else:
            if args.selection_tsv is not None:
                raise ValueError(
                    "--selection-tsv v1 is specific to annotation_missing_gene"
                )
            if args.count is None:
                raise ValueError("--count is required for structural perturbations")
            if (
                args.event_type == COPY_COLLAPSE_EVENT
                and args.pair_tsv is None
            ):
                raise ValueError("annotation_copy_collapse requires --pair-tsv")
            if (
                args.event_type != COPY_COLLAPSE_EVENT
                and args.pair_tsv is not None
            ):
                raise ValueError(
                    "--pair-tsv is only valid for annotation_copy_collapse"
                )
            manifest = generate_structure_benchmark(
                gff_path=args.gff,
                output_dir=args.output_dir,
                event_type=args.event_type,
                count=args.count,
                seed=args.seed,
                truth_output_dir=args.truth_dir,
                pair_tsv_path=args.pair_tsv,
            )
        print(json.dumps(manifest, indent=2, sort_keys=True))
        return 0
    if args.command == "benchmark" and args.benchmark_command == "restore":
        report = restore_gff_from_truth(
            perturbed_gff_path=args.perturbed_gff,
            truth_path=args.truth,
            output_gff_path=args.output_gff,
        )
        print(json.dumps(report, indent=2, sort_keys=True))
        return 0
    if args.command == "benchmark" and args.benchmark_command == "score":
        report = score_annotation_repair(
            source_gff_path=args.source_gff,
            perturbed_gff_path=args.perturbed_gff,
            candidate_gff_path=args.candidate_gff,
            truth_path=args.truth,
            include_event_details=args.include_event_details,
            event_strata_path=args.event_strata,
            stratum_columns=tuple(args.stratum_column),
            control_candidate_gff_path=args.control_candidate_gff,
        )
        print(json.dumps(report, indent=2, sort_keys=True))
        return 0
    if (
        args.command == "benchmark"
        and args.benchmark_command == "score-structure-hypotheses"
    ):
        report = score_structure_hypotheses(
            source_gff_path=args.source_gff,
            perturbed_gff_path=args.perturbed_gff,
            candidate_gff_path=args.candidate_gff,
            hypotheses_tsv_path=args.hypotheses_tsv,
            truth_path=args.truth,
            control_hypotheses_tsv_path=args.control_hypotheses_tsv,
            include_event_details=args.include_event_details,
        )
        print(json.dumps(report, indent=2, sort_keys=True))
        return 0
    if args.command == "benchmark" and args.benchmark_command == "bootstrap-events":
        score_inputs = []
        for value in args.score:
            label, separator, path = value.partition("=")
            if not separator or not label or not path:
                raise ValueError("--score must use LABEL=score.json")
            score_inputs.append((label, path))
        report = paired_event_bootstrap(
            score_inputs=score_inputs,
            output_json_path=args.output_json,
            metric=args.metric,
            replicates=args.replicates,
            seed=args.seed,
            alpha=args.alpha,
        )
        print(json.dumps(report, indent=2, sort_keys=True))
        return 0
    if (
        args.command == "benchmark"
        and args.benchmark_command == "bootstrap-independent-events"
    ):
        score_inputs = []
        for value in args.score:
            label, separator, path = value.partition("=")
            if not separator or not label or not path:
                raise ValueError("--score must use LABEL=score.json")
            score_inputs.append((label, path))
        report = independent_event_bootstrap(
            score_inputs=score_inputs,
            output_json_path=args.output_json,
            metric=args.metric,
            replicates=args.replicates,
            seed=args.seed,
            alpha=args.alpha,
        )
        print(json.dumps(report, indent=2, sort_keys=True))
        return 0
    if (
        args.command == "benchmark"
        and args.benchmark_command == "bootstrap-confusion"
    ):
        score_inputs = []
        for value in args.score:
            label, separator, path = value.partition("=")
            if not separator or not label or not path:
                raise ValueError("--score must use LABEL=score.json")
            score_inputs.append((label, path))
        report = independent_confusion_bootstrap(
            score_inputs=score_inputs,
            output_json_path=args.output_json,
            section=args.section,
            replicates=args.replicates,
            seed=args.seed,
            alpha=args.alpha,
        )
        print(json.dumps(report, indent=2, sort_keys=True))
        return 0
    if (
        args.command == "benchmark"
        and args.benchmark_command == "label-copy-features"
    ):
        manifest = label_copy_candidate_features(
            feature_tsv_path=args.features,
            hidden_truth_json_path=args.truth,
            output_tsv_path=args.output_tsv,
        )
        print(json.dumps(manifest, indent=2, sort_keys=True))
        return 0
    if (
        args.command == "benchmark"
        and args.benchmark_command == "score-copy-ranking"
    ):
        report = evaluate_copy_candidate_scores(
            scored_tsv_path=args.scores,
            labeled_feature_tsv_path=args.labeled_features,
            output_json_path=args.output_json,
        )
        print(json.dumps(report, indent=2, sort_keys=True))
        return 0
    if (
        args.command == "benchmark"
        and args.benchmark_command == "sample-copy-pairs"
    ):
        manifest = sample_balanced_copy_pairs(
            source_gff_path=args.source_gff,
            pair_tsv_path=args.pairs,
            output_pair_tsv_path=args.output_pairs,
            decisions_tsv_path=args.decisions_tsv,
            count=args.count,
            seed=args.seed,
        )
        print(json.dumps(manifest, indent=2, sort_keys=True))
        return 0
    if (
        args.command == "benchmark"
        and args.benchmark_command == "bootstrap-structure-hypotheses"
    ):
        score_inputs = []
        for value in args.score:
            label, separator, path = value.partition("=")
            if not separator or not label or not path:
                raise ValueError("--score must use LABEL=score.json")
            score_inputs.append((label, path))
        report = paired_structure_hypothesis_bootstrap(
            score_inputs=score_inputs,
            output_json_path=args.output_json,
            replicates=args.replicates,
            seed=args.seed,
            alpha=args.alpha,
        )
        print(json.dumps(report, indent=2, sort_keys=True))
        return 0
    if (
        args.command == "benchmark"
        and args.benchmark_command == "score-localization"
    ):
        report = score_synteny_localization(
            gap_tsv_paths=args.gaps,
            selection_tsv_path=args.selection,
            truth_path=args.truth,
            include_event_details=args.include_event_details,
        )
        print(json.dumps(report, indent=2, sort_keys=True))
        return 0
    if (
        args.command == "benchmark"
        and args.benchmark_command == "label-synteny-models"
    ):
        manifest = write_synteny_model_labels(
            source_gff_path=args.source_gff,
            candidate_gff_path=args.candidate_gff,
            selection_tsv_path=args.selection,
            baseline_decisions_tsv_path=args.baseline_decisions,
            truth_path=args.truth,
            output_tsv_path=args.output_tsv,
            control_candidate_gff_path=args.control_candidate_gff,
        )
        print(json.dumps(manifest, indent=2, sort_keys=True))
        return 0
    if (
        args.command == "benchmark"
        and args.benchmark_command == "calibrate-synteny-tiers"
    ):
        report = calibrate_synteny_tiers(
            labeled_tsv_path=args.labels,
            output_json_path=args.output,
            label_column=args.label_column,
            high_precision_floor=args.high_precision_floor,
            high_precision_min_selected=args.high_precision_min_selected,
            eligible_column=args.eligible_column,
            feature_set=args.feature_set,
        )
        print(json.dumps(report, indent=2, sort_keys=True))
        return 0
    if args.command == "benchmark" and args.benchmark_command == "catalog":
        manifest = write_missing_gene_candidate_catalog(
            gff_path=args.gff,
            output_tsv_path=args.output_tsv,
            external_strata_path=args.external_strata,
            external_strata_prefix=args.external_strata_prefix,
            primary_chromosomes_only=args.primary_chromosomes_only,
        )
        print(json.dumps(manifest, indent=2, sort_keys=True))
        return 0
    if args.command == "benchmark" and args.benchmark_command == "summarize-catalog":
        report = write_candidate_catalog_summary(
            catalog_path=args.catalog,
            output_path=args.output,
            columns=args.column,
            crossings=args.cross,
        )
        print(json.dumps(report, indent=2, sort_keys=True))
        return 0
    if args.command == "benchmark" and args.benchmark_command == "sample-catalog":
        manifest = sample_candidate_catalog(
            catalog_path=args.catalog,
            plan_path=args.plan,
            output_tsv_path=args.output_tsv,
            seed=args.seed,
            exclude_tsv_paths=tuple(args.exclude_selection),
        )
        print(json.dumps(manifest, indent=2, sort_keys=True))
        return 0
    if args.command == "benchmark" and args.benchmark_command == "mask-genome":
        manifest = create_masked_gap_control(
            source_genome_path=args.genome,
            hidden_truth_path=args.truth,
            output_genome_path=args.output_genome,
            output_mask_truth_path=args.output_mask_truth,
            background_gff_path=args.background_gff,
        )
        print(json.dumps(manifest, indent=2, sort_keys=True))
        return 0
    if args.command == "benchmark" and args.benchmark_command == "score-abstention":
        report = score_masked_gap_abstention(
            perturbed_gff_path=args.perturbed_gff,
            candidate_gff_path=args.candidate_gff,
            mask_truth_path=args.mask_truth,
            include_event_details=args.include_event_details,
        )
        print(json.dumps(report, indent=2, sort_keys=True))
        return 0
    if args.command == "benchmark" and args.benchmark_command == "audit-mask":
        report = audit_masked_gap_genome(
            source_genome_path=args.source_genome,
            masked_genome_path=args.masked_genome,
            mask_truth_path=args.mask_truth,
        )
        print(json.dumps(report, indent=2, sort_keys=True))
        return 0
    if (
        args.command == "benchmark"
        and args.benchmark_command == "summarize-mask-selection"
    ):
        report = summarize_masked_gap_selection(
            mask_truth_path=args.mask_truth,
            hidden_truth_path=args.hidden_truth,
            strata_tsv_path=args.strata_tsv,
            columns=tuple(args.column),
        )
        print(json.dumps(report, indent=2, sort_keys=True))
        return 0
    if args.command == "evidence" and args.evidence_command == "prepare-wgdi":
        manifest = prepare_wgdi_inputs(
            gff_path=args.gff,
            protein_path=args.protein,
            fai_path=args.fai,
            output_dir=args.output_dir,
            prefix=args.prefix,
            min_genes_per_seqid=args.min_genes_per_seqid,
            primary_chromosomes_only=args.primary_chromosomes_only,
        )
        print(json.dumps(manifest, indent=2, sort_keys=True))
        return 0
    if args.command == "evidence" and args.evidence_command == "summarize-wgdi":
        manifest = summarize_wgdi_gene_evidence(
            query_gff_path=args.query_gff,
            query_input_manifest_path=args.query_input_manifest,
            collinearity_inputs=args.collinearity,
            expected_source_inputs=args.expected_source,
            output_gene_tsv_path=args.output_gene_tsv,
            output_block_tsv_path=args.output_block_tsv,
        )
        print(json.dumps(manifest, indent=2, sort_keys=True))
        return 0
    if (
        args.command == "evidence"
        and args.evidence_command == "infer-homeolog-pairs"
    ):
        if len(args.subgenome) != 2:
            raise ValueError("--subgenome must be repeated exactly twice")
        manifest = infer_wgdi_homeolog_pairs(
            gene_evidence_tsv_path=args.gene_evidence,
            collinearity_inputs=args.collinearity,
            output_pair_tsv_path=args.output_pairs,
            decisions_tsv_path=args.decisions_tsv,
            wgd_event=args.wgd_event,
            subgenomes=tuple(args.subgenome),
            source_group_map_path=args.source_group_map,
            min_support_group_count=args.min_support_group_count,
            require_reciprocal_unique=not args.allow_multiple_partners,
        )
        print(json.dumps(manifest, indent=2, sort_keys=True))
        return 0
    if (
        args.command == "evidence"
        and args.evidence_command == "infer-outgroup-duplicated-pairs"
    ):
        manifest = infer_outgroup_duplicated_pairs(
            query_wgdi_gff_path=args.query_wgdi_gff,
            source_gff_path=args.source_gff,
            collinearity_inputs=args.collinearity,
            output_pair_tsv_path=args.output_pairs,
            decisions_tsv_path=args.decisions_tsv,
            wgd_event=args.wgd_event,
            source_group_map_path=args.source_group_map,
            min_support_group_count=args.min_support_group_count,
            min_block_pairs=args.min_block_pairs,
            require_cross_seqid=not args.allow_same_seqid,
            require_reciprocal_unique=not args.allow_multiple_partners,
        )
        print(json.dumps(manifest, indent=2, sort_keys=True))
        return 0
    if (
        args.command == "evidence"
        and args.evidence_command == "infer-self-wgd-pairs"
    ):
        manifest = infer_self_wgdi_pairs(
            query_wgdi_gff_path=args.query_wgdi_gff,
            collinearity_path=args.collinearity,
            source_gff_path=args.source_gff,
            output_pair_tsv_path=args.output_pairs,
            decisions_tsv_path=args.decisions_tsv,
            wgd_event=args.wgd_event,
            min_block_pairs=args.min_block_pairs,
            require_different_seqids=not args.allow_same_seqid,
            require_reciprocal_unique=not args.allow_multiple_partners,
        )
        print(json.dumps(manifest, indent=2, sort_keys=True))
        return 0
    if (
        args.command == "evidence"
        and args.evidence_command == "aggregate-reference-anchored"
    ):
        manifest = aggregate_reference_anchored_projection(
            candidate_source_provenance_path=args.candidate_source_provenance,
            accepted_wgd_pair_inputs=args.accepted_wgd_pairs,
            homeolog_base_evidence_inputs=args.homeolog_base_evidence,
            output_dir_path=args.output_dir,
            evidence_type=args.evidence_type,
        )
        print(json.dumps(manifest, indent=2, sort_keys=True))
        return 0
    if (
        args.command == "evidence"
        and args.evidence_command == "intersect-copy-pair-evidence"
    ):
        manifest = intersect_copy_pair_evidence(
            pair_inputs=args.pairs,
            output_pair_tsv_path=args.output_pairs,
            decisions_tsv_path=args.decisions_tsv,
            pair_set_label=args.pair_set_label,
            require_reciprocal_unique=not args.allow_multiple_partners,
        )
        print(json.dumps(manifest, indent=2, sort_keys=True))
        return 0
    if args.command == "evidence" and args.evidence_command == "infer-synteny-gaps":
        manifest = infer_wgdi_synteny_gaps(
            query_wgdi_gff_path=args.query_wgdi_gff,
            target_wgdi_gff_path=args.target_wgdi_gff,
            collinearity_path=args.collinearity,
            source_label=args.source_label,
            output_tsv_path=args.output_tsv,
            expected_chromosome_pair_tsv_path=args.expected_chromosome_pairs,
            max_query_intervening_genes=args.max_query_intervening_genes,
            min_target_excess_genes=args.min_target_excess_genes,
            max_target_gap_genes=args.max_target_gap_genes,
            max_query_locus_bp=args.max_query_locus_bp,
        )
        print(json.dumps(manifest, indent=2, sort_keys=True))
        return 0
    if args.command == "evidence" and args.evidence_command == "summarize-natural":
        report = summarize_natural_validation(
            validation_tsv_path=args.validation_tsv,
            output_json_path=args.output_json,
        )
        print(json.dumps(report, indent=2, sort_keys=True))
        return 0
    if (
        args.command == "evidence"
        and args.evidence_command == "prepare-natural-graph"
    ):
        manifest = prepare_natural_graph_inputs(
            validation_tsv_path=args.validation_tsv,
            output_candidate_tsv_path=args.output_candidates,
            output_evidence_tsv_path=args.output_evidence,
        )
        print(json.dumps(manifest, indent=2, sort_keys=True))
        return 0
    if (
        args.command == "evidence"
        and args.evidence_command == "select-synteny-gap-models"
    ):
        manifest = select_synteny_gap_models(
            gap_tsv_paths=args.gaps,
            baseline_decisions_tsv_path=args.baseline_decisions,
            adapted_candidate_gff_path=args.adapted_candidate_gff,
            output_selection_tsv_path=args.output_selection,
            output_candidate_gff_path=args.output_candidate_gff,
        )
        print(json.dumps(manifest, indent=2, sort_keys=True))
        return 0
    if (
        args.command == "evidence"
        and args.evidence_command == "summarize-projection-support"
    ):
        manifest = summarize_projection_support(
            decisions_tsv_path=args.decisions,
            output_tsv_path=args.output_tsv,
        )
        print(json.dumps(manifest, indent=2, sort_keys=True))
        return 0
    if (
        args.command == "evidence"
        and args.evidence_command == "build-copy-features"
    ):
        manifest = build_copy_candidate_features(
            consensus_decisions_tsv_path=args.consensus_decisions,
            method_decision_inputs=args.method_decisions,
            wgd_selection_tsv_path=args.wgd_selection,
            output_tsv_path=args.output_tsv,
        )
        print(json.dumps(manifest, indent=2, sort_keys=True))
        return 0
    if (
        args.command == "evidence"
        and args.evidence_command == "build-homeolog-topology-features"
    ):
        manifest = build_homeolog_topology_features(
            copy_feature_tsv_path=args.copy_features,
            wgd_selection_tsv_path=args.wgd_selection,
            candidate_gff_path=args.candidate_gff,
            base_gff_path=args.base_gff,
            output_tsv_path=args.output_tsv,
        )
        print(json.dumps(manifest, indent=2, sort_keys=True))
        return 0
    if (
        args.command == "evidence"
        and args.evidence_command == "build-fixed-target-backbone"
    ):
        primary_seqids, _ = read_primary_seqid_table(args.primary_seqid_table)
        manifest = build_fixed_target_backbone(
            base_gff_path=args.base_gff,
            query_wgdi_gff_path=args.query_wgdi_gff,
            base_only_collinearity_path=args.base_only_collinearity,
            policy=FixedTargetBackbonePolicy(
                primary_seqids=tuple(sorted(primary_seqids)),
                min_block_pairs=args.min_block_pairs,
            ),
            output_dir_path=args.output_dir,
        )
        print(json.dumps(manifest, indent=2, sort_keys=True))
        return 0
    if (
        args.command == "evidence"
        and args.evidence_command == "project-fixed-backbone"
    ):
        manifest = project_candidates_to_fixed_backbone(
            backbone_dir_path=args.backbone_dir,
            candidate_catalog_tsv_path=args.candidate_catalog,
            fixed_base_hits_tsv_path=args.fixed_base_hits,
            fixed_base_hits_manifest_path=args.fixed_base_hits_manifest,
            policy=FixedBackboneProjectionPolicy(
                minimum_query_coverage=args.minimum_query_coverage,
                maximum_evalue=args.maximum_evalue,
            ),
            output_dir_path=args.output_dir,
        )
        print(json.dumps(manifest, indent=2, sort_keys=True))
        return 0
    if (
        args.command == "evidence"
        and args.evidence_command == "evaluate-species-applicability"
    ):
        report = evaluate_species_applicability(
            metrics_path=args.metrics,
            policy_path=args.policy,
            output_path=args.output_json,
        )
        print(json.dumps(report, indent=2, sort_keys=True))
        return 0
    if (
        args.command == "evidence"
        and args.evidence_command == "score-homeolog-copy-candidates"
    ):
        manifest = score_homeolog_copy_candidates(
            copy_feature_tsv_path=args.copy_features,
            topology_tsv_path=args.topology_features,
            model_json_path=args.model_json,
            output_tsv_path=args.output_tsv,
        )
        print(json.dumps(manifest, indent=2, sort_keys=True))
        return 0
    if (
        args.command == "evidence"
        and args.evidence_command == "freeze-homeolog-review-rankings"
    ):
        manifest = freeze_homeolog_review_rankings(
            score_tsv_path=args.scores,
            output_tsv_path=args.output_tsv,
            review_budgets=(
                args.review_budget
                if args.review_budget is not None
                else (25, 50, 100, 200)
            ),
        )
        print(json.dumps(manifest, indent=2, sort_keys=True))
        return 0
    if (
        args.command == "evidence"
        and args.evidence_command == "score-support-conditioned-candidates"
    ):
        manifest = score_support_conditioned_candidates(
            copy_feature_tsv_path=args.copy_features,
            topology_tsv_path=args.topology_features,
            model_json_path=args.model_json,
            output_tsv_path=args.output_tsv,
        )
        print(json.dumps(manifest, indent=2, sort_keys=True))
        return 0
    if (
        args.command == "evidence"
        and args.evidence_command == "apply-conflict-winner-guard"
    ):
        manifest = apply_conflict_winner_guard(
            v03_score_tsv_path=args.v03_scores,
            pool_decisions_tsv_path=args.pool_decisions,
            pool_manifest_json_path=args.pool_manifest,
            output_tsv_path=args.output_tsv,
        )
        print(json.dumps(manifest, indent=2, sort_keys=True))
        return 0
    if args.command == "evidence" and args.evidence_command == "prepare-b73-isoseq":
        manifest = prepare_b73_isoseq_transcripts(
            fasta_path=args.fasta,
            count_csv_path=args.counts_csv,
            output_fasta_path=args.output_fasta,
            output_count_tsv_path=args.output_counts,
            minimum_b73_full_length_reads=args.minimum_b73_reads,
        )
        print(json.dumps(manifest, indent=2, sort_keys=True))
        return 0
    if (
        args.command == "evidence"
        and args.evidence_command == "filter-candidate-query-paf"
    ):
        paf_inputs = []
        for specification in args.paf_input:
            accession, separator, path = specification.partition("=")
            if not separator or not accession or not path:
                raise ValueError(
                    "--paf-input must use the exact ACCESSION=PATH form"
                )
            paf_inputs.append((accession, path))
        manifest = filter_candidate_query_paf(
            candidate_gff_path=args.candidate_gff,
            paf_inputs=paf_inputs,
            output_paf_path=args.output_paf,
            output_count_tsv_path=args.output_counts,
            output_summary_tsv_path=args.output_summary,
            output_manifest_json_path=args.output_manifest,
            alignment_strand_source=args.alignment_strand_source,
        )
        print(json.dumps(manifest, indent=2, sort_keys=True))
        return 0
    if (
        args.command == "evidence"
        and args.evidence_command == "validate-isoseq-candidates"
    ):
        manifest = validate_isoseq_candidate_chains(
            candidate_gff_path=args.candidate_gff,
            paf_path=args.paf,
            selected_count_tsv_path=args.selected_counts,
            genome_fasta_path=args.genome_fasta,
            output_evidence_tsv_path=args.output_evidence,
            minimum_query_coverage=args.minimum_query_coverage,
            minimum_identity=args.minimum_identity,
            minimum_mapq=args.minimum_mapq,
            maximum_secondary_score_fraction=args.maximum_secondary_score_fraction,
            minimum_candidate_cds_coverage=args.minimum_candidate_cds_coverage,
            flank_bp=args.flank_bp,
            alignment_strand_source=args.alignment_strand_source,
        )
        print(json.dumps(manifest, indent=2, sort_keys=True))
        return 0
    if (
        args.command == "evidence"
        and args.evidence_command == "join-isoseq-review-rankings"
    ):
        summary = join_isoseq_review_rankings(
            evidence_tsv_path=args.evidence,
            review_rankings_tsv_path=args.review_rankings,
            output_tsv_path=args.output_tsv,
            output_summary_json_path=args.output_summary,
            review_budgets=args.review_budget or (25, 50, 100, 200),
            comparator_estimator=args.comparator_estimator,
            primary_estimator=args.primary_estimator,
        )
        print(json.dumps(summary, indent=2, sort_keys=True))
        return 0
    if (
        args.command == "evidence"
        and args.evidence_command == "bootstrap-isoseq-review-yield"
    ):
        report = bootstrap_isoseq_review_yield(
            evidence_tsv_path=args.evidence,
            review_rankings_tsv_path=args.review_rankings,
            output_json_path=args.output_json,
            review_budgets=args.review_budget or (25, 50, 100, 200),
            replicates=args.replicates,
            seed=args.seed,
            alpha=args.alpha,
            comparator_estimator=args.comparator_estimator,
            primary_estimator=args.primary_estimator,
        )
        print(json.dumps(report, indent=2, sort_keys=True))
        return 0
    if (
        args.command == "evidence"
        and args.evidence_command == "export-natural-candidate-cds"
    ):
        manifest = export_natural_candidate_cds(
            candidate_gff_path=args.candidate_gff,
            genome_fasta_path=args.genome_fasta,
            output_fasta_path=args.output_fasta,
        )
        print(json.dumps(manifest, indent=2, sort_keys=True))
        return 0
    if (
        args.command == "evidence"
        and args.evidence_command == "audit-natural-candidates"
    ):
        report = audit_natural_candidates(
            candidate_gff_path=args.candidate_gff,
            base_gff_path=args.base_gff,
            genome_fasta_path=args.genome_fasta,
            review_rankings_tsv_path=args.review_rankings,
            isoseq_evidence_tsv_path=args.isoseq_evidence,
            self_map_paf_path=args.self_map_paf,
            repeat_gff_path=args.repeat_gff,
            repeat_seqid_map_tsv_path=args.repeat_seqid_map,
            repeat_flank_bp=args.repeat_flank_bp,
            minimum_full_length_read_support=args.minimum_full_length_read_support,
            output_tsv_path=args.output_tsv,
            output_summary_json_path=args.output_summary,
            review_budgets=args.review_budget or (25, 50, 100, 200),
            minimum_query_coverage=args.minimum_query_coverage,
            minimum_identity=args.minimum_identity,
            near_equal_score_fraction=args.near_equal_score_fraction,
        )
        print(json.dumps(report, indent=2, sort_keys=True))
        return 0
    if (
        args.command == "evidence"
        and args.evidence_command == "score-copy-candidates"
    ):
        manifest = score_copy_candidate_features(
            feature_tsv_path=args.features,
            model_json_path=args.model_json,
            output_tsv_path=args.output_tsv,
            mask_feature_groups=args.mask_feature_group,
        )
        print(json.dumps(manifest, indent=2, sort_keys=True))
        return 0
    if (
        args.command == "evidence"
        and args.evidence_command == "select-scored-copy-candidates"
    ):
        manifest = select_scored_copy_candidates(
            base_gff_path=args.base_gff,
            candidate_gff_path=args.candidate_gff,
            scored_tsv_path=args.scores,
            model_json_path=args.model_json,
            output_gff_path=args.output_gff,
            selection_tsv_path=args.selection_tsv,
            policy=args.policy,
        )
        print(json.dumps(manifest, indent=2, sort_keys=True))
        return 0
    if (
        args.command == "evidence"
        and args.evidence_command == "select-projection-support"
    ):
        manifest = select_projection_support_models(
            candidate_gff_path=args.candidate_gff,
            projection_support_tsv_path=args.projection_support,
            output_gff_path=args.output_gff,
            selection_tsv_path=args.selection_tsv,
            min_support_group_count=args.min_support_group_count,
            source_group_map_path=args.source_group_map,
        )
        print(json.dumps(manifest, indent=2, sort_keys=True))
        return 0
    if (
        args.command == "evidence"
        and args.evidence_command == "select-wgd-supported-candidates"
    ):
        manifest = select_wgd_supported_candidates(
            base_gff_path=args.base_gff,
            candidate_gff_path=args.candidate_gff,
            pair_tsv_path=args.pairs,
            output_gff_path=args.output_gff,
            selection_tsv_path=args.selection_tsv,
        )
        print(json.dumps(manifest, indent=2, sort_keys=True))
        return 0
    if (
        args.command == "evidence"
        and args.evidence_command == "propagate-wgd-conflict-partners"
    ):
        manifest = propagate_wgd_selection_to_conflict_pool(
            base_gff_path=args.base_gff,
            candidate_gff_path=args.candidate_gff,
            pool_decisions_tsv_path=args.pool_decisions,
            prior_wgd_selection_tsv_path=args.prior_wgd_selection,
            output_selection_tsv_path=args.output_selection,
        )
        print(json.dumps(manifest, indent=2, sort_keys=True))
        return 0
    if (
        args.command == "evidence"
        and args.evidence_command == "extract-bam-junctions"
    ):
        manifest = extract_bam_junctions(
            bam_path=args.bam,
            output_tsv_path=args.output_tsv,
            sample=args.sample,
            samtools_path=args.samtools,
            threads=args.threads,
            min_mapq=args.min_mapq,
        )
        print(json.dumps(manifest, indent=2, sort_keys=True))
        return 0
    if (
        args.command == "evidence"
        and args.evidence_command == "aggregate-junctions"
    ):
        input_paths = sorted(
            path
            for input_dir in args.input_dir
            for path in Path(input_dir).glob("*.junctions.tsv")
        )
        manifest = aggregate_junctions(
            junction_tsv_paths=input_paths,
            output_tsv_path=args.output_tsv,
            primary_samples=frozenset(args.primary_sample),
            min_reads_per_sample=args.min_reads_per_sample,
            min_supporting_samples=args.min_supporting_samples,
            sample_group_tsv_path=args.sample_groups,
            min_samples_per_group=args.min_samples_per_group,
            min_secondary_groups=args.min_secondary_groups,
        )
        print(json.dumps(manifest, indent=2, sort_keys=True))
        return 0
    if args.command == "evidence" and args.evidence_command == "discover-natural":
        manifest = discover_natural_candidates(
            target_gff_path=args.target_gff,
            miniprot_gff_path=args.miniprot_gff,
            protein_map_path=args.protein_map,
            output_tsv_path=args.output_tsv,
            min_identity=args.min_identity,
            min_query_coverage=args.min_query_coverage,
            max_existing_cds_overlap=args.max_existing_cds_overlap,
            min_boundary_extension_bp=args.min_boundary_extension_bp,
            near_best_score_fraction=args.near_best_score_fraction,
        )
        print(json.dumps(manifest, indent=2, sort_keys=True))
        return 0
    if (
        args.command == "evidence"
        and args.evidence_command == "validate-natural-rna"
    ):
        manifest = validate_natural_candidates_with_rna(
            candidate_tsv_path=args.candidates,
            junction_aggregate_tsv_path=args.junctions,
            output_tsv_path=args.output_tsv,
        )
        print(json.dumps(manifest, indent=2, sort_keys=True))
        return 0
    if (
        args.command == "evidence"
        and args.evidence_command == "annotate-natural-assembly"
    ):
        manifest = annotate_natural_assembly_context(
            validation_tsv_path=args.validation_tsv,
            genome_fasta_path=args.genome,
            genome_fai_path=args.fai,
            output_tsv_path=args.output_tsv,
            flank_bp=args.flank_bp,
            max_ambiguous_fraction=args.max_ambiguous_fraction,
        )
        print(json.dumps(manifest, indent=2, sort_keys=True))
        return 0
    if (
        args.command == "evidence"
        and args.evidence_command == "validate-natural-secondary-rna"
    ):
        manifest = validate_natural_candidates_with_secondary_groups(
            primary_validation_tsv_path=args.primary_validation,
            grouped_junction_tsv_path=args.grouped_junctions,
            output_tsv_path=args.output_tsv,
        )
        print(json.dumps(manifest, indent=2, sort_keys=True))
        return 0
    if (
        args.command == "evidence"
        and args.evidence_command == "extract-paf-deletions"
    ):
        manifest = extract_paf_target_deletions(
            paf_path=args.paf,
            output_tsv_path=args.output_tsv,
            min_deletion_bp=args.min_deletion_bp,
            min_flanking_anchor_bp=args.min_flanking_anchor_bp,
            min_mapq=args.min_mapq,
            require_primary=not args.allow_secondary,
            chromosome_pair_tsv_path=args.chromosome_pairs,
            query_genome_path=args.query_genome,
            query_fai_path=args.query_fai,
            query_flank_bp=args.query_flank_bp,
            max_query_flank_ambiguous_fraction=(
                args.max_query_flank_ambiguous_fraction
            ),
            require_clean_query_flanks=not args.retain_uncertain_query_flanks,
        )
        print(json.dumps(manifest, indent=2, sort_keys=True))
        return 0
    if (
        args.command == "evidence"
        and args.evidence_command == "combine-chromosome-pafs"
    ):
        manifest = combine_chromosome_pafs(
            input_dir_path=args.input_dir,
            chromosome_pair_tsv_path=args.chromosome_pairs,
            output_paf_path=args.output_paf,
        )
        print(json.dumps(manifest, indent=2, sort_keys=True))
        return 0
    if args.command == "evidence" and args.evidence_command == "catalog-pav-genes":
        manifest = catalog_genes_in_target_deletions(
            gff_path=args.gff,
            deletion_tsv_path=args.deletions,
            output_tsv_path=args.output_tsv,
            min_gene_coverage=args.min_gene_coverage,
            min_cds_coverage=args.min_cds_coverage,
        )
        print(json.dumps(manifest, indent=2, sort_keys=True))
        return 0
    if (
        args.command == "evidence"
        and args.evidence_command == "reconcile-pav-catalogs"
    ):
        manifest = reconcile_pav_gene_catalogs(
            catalog_inputs=args.catalog,
            output_tsv_path=args.output_tsv,
            min_support=args.min_support,
        )
        print(json.dumps(manifest, indent=2, sort_keys=True))
        return 0
    if args.command == "evidence" and args.evidence_command == "subset-pav-proteins":
        manifest = subset_catalog_proteins(
            catalog_tsv_path=args.catalog,
            protein_fasta_path=args.proteins,
            output_fasta_path=args.output_fasta,
        )
        print(json.dumps(manifest, indent=2, sort_keys=True))
        return 0
    if args.command == "evidence" and args.evidence_command == "screen-pav-proteins":
        manifest = screen_pav_candidate_proteins(
            catalog_tsv_path=args.catalog,
            miniprot_gff_path=args.miniprot_gff,
            protein_map_path=args.protein_map,
            chromosome_pair_tsv_path=args.chromosome_pairs,
            output_tsv_path=args.output_tsv,
            min_identity=args.min_identity,
            min_query_coverage=args.min_query_coverage,
            min_partial_query_coverage=args.min_partial_query_coverage,
            max_local_distance_bp=args.max_local_distance_bp,
            max_breakpoint_spread_bp=args.max_breakpoint_spread_bp,
        )
        print(json.dumps(manifest, indent=2, sort_keys=True))
        return 0
    if args.command == "evidence" and args.evidence_command == "annotate-pav-wgdi":
        manifest = annotate_pav_with_wgdi(
            protein_screen_tsv_path=args.protein_screen,
            wgdi_gene_tsv_path=args.wgdi_genes,
            output_tsv_path=args.output_tsv,
        )
        print(json.dumps(manifest, indent=2, sort_keys=True))
        return 0
    if args.command == "evidence" and args.evidence_command == "collapse-pav-events":
        manifest = collapse_high_confidence_pav_events(
            annotated_gene_tsv_path=args.annotated_genes,
            output_tsv_path=args.output_tsv,
        )
        print(json.dumps(manifest, indent=2, sort_keys=True))
        return 0
    if args.command == "normalize" and args.normalize_command == "ncbi-primary":
        manifest = prepare_ncbi_primary_bundle(
            gff_path=args.gff,
            protein_path=args.protein,
            cds_path=args.cds,
            genome_path=args.genome,
            output_dir=args.output_dir,
        )
        print(json.dumps(manifest, indent=2, sort_keys=True))
        return 0
    if args.command == "normalize" and args.normalize_command == "provider-gff3":
        manifest = normalize_provider_gff3(
            gff_path=args.gff,
            output_dir=args.output_dir,
            repair_unescaped_note_semicolons=(
                args.repair_unescaped_note_semicolons
            ),
            drop_invalid_intron_intervals=args.drop_invalid_intron_intervals,
            strip_embedded_fasta=args.strip_embedded_fasta,
        )
        print(json.dumps(manifest, indent=2, sort_keys=True))
        return 0
    if args.command == "normalize" and args.normalize_command == "primary-annotation":
        manifest = prepare_primary_annotation_bundle(
            gff_path=args.gff,
            genome_path=args.genome,
            output_dir=args.output_dir,
            primary_seqid_table_path=args.primary_seqid_table,
        )
        print(json.dumps(manifest, indent=2, sort_keys=True))
        return 0
    if args.command == "baseline" and args.baseline_command == "prepare-proteins":
        manifest = prepare_reference_proteins(
            protein_inputs=args.protein,
            output_fasta_path=args.output_fasta,
            output_map_path=args.output_map,
        )
        print(json.dumps(manifest, indent=2, sort_keys=True))
        return 0
    if args.command == "baseline" and args.baseline_command == "merge-candidate-gffs":
        manifest = merge_candidate_gffs(
            candidate_inputs=args.candidate,
            output_gff_path=args.output_gff,
            provenance_tsv_path=args.provenance_tsv,
        )
        print(json.dumps(manifest, indent=2, sort_keys=True))
        return 0
    if args.command == "baseline" and args.baseline_command == "adapt-miniprot":
        manifest = adapt_miniprot_baseline(
            perturbed_gff_path=args.perturbed_gff,
            miniprot_gff_path=args.miniprot_gff,
            protein_map_path=args.protein_map,
            output_gff_path=args.output_gff,
            decisions_tsv_path=args.decisions_tsv,
            min_identity=args.min_identity,
            min_query_coverage=args.min_query_coverage,
            max_existing_cds_overlap=args.max_existing_cds_overlap,
            max_redundancy_overlap=args.max_redundancy_overlap,
            require_intact=not args.allow_disrupted,
        )
        print(json.dumps(manifest, indent=2, sort_keys=True))
        return 0
    if args.command == "baseline" and args.baseline_command == "adapt-gff":
        manifest = adapt_annotation_gff_baseline(
            perturbed_gff_path=args.perturbed_gff,
            candidate_gff_path=args.candidate_gff,
            source=args.source,
            output_gff_path=args.output_gff,
            decisions_tsv_path=args.decisions_tsv,
            score_attribute=args.score_attribute,
            max_existing_cds_overlap=args.max_existing_cds_overlap,
            max_redundancy_overlap=args.max_redundancy_overlap,
            infer_missing_cds_phase=args.infer_missing_cds_phase,
        )
        print(json.dumps(manifest, indent=2, sort_keys=True))
        return 0
    if (
        args.command == "baseline"
        and args.baseline_command == "select-method-consensus"
    ):
        candidate_inputs = []
        for value in args.candidate:
            label, separator, path = value.partition("=")
            if not separator or not label or not path:
                raise ValueError("--candidate must use METHOD=GFF3")
            candidate_inputs.append((label, path))
        manifest = select_method_consensus(
            base_gff_path=args.base_gff,
            candidate_inputs=candidate_inputs,
            output_gff_path=args.output_gff,
            decisions_tsv_path=args.decisions_tsv,
            min_method_support=args.min_method_support,
            max_redundancy_overlap=args.max_redundancy_overlap,
            redundancy_policy=args.redundancy_policy,
        )
        print(json.dumps(manifest, indent=2, sort_keys=True))
        return 0
    if (
        args.command == "baseline"
        and args.baseline_command == "adapt-miniprot-structure"
    ):
        manifest = adapt_miniprot_structure_candidates(
            annotation_gff_path=args.annotation_gff,
            miniprot_gff_path=args.miniprot_gff,
            protein_map_path=args.protein_map,
            output_gff_path=args.output_gff,
            decisions_tsv_path=args.decisions_tsv,
            min_identity=args.min_identity,
            min_query_coverage=args.min_query_coverage,
            min_source_support=args.min_source_support,
            min_gene_overlap_fraction=args.min_gene_overlap_fraction,
            require_intact=not args.allow_disrupted,
            source_group_map_path=args.source_group_map,
        )
        print(json.dumps(manifest, indent=2, sort_keys=True))
        return 0
    if (
        args.command == "baseline"
        and args.baseline_command == "infer-structure-hypotheses"
    ):
        manifest = infer_structure_hypotheses(
            annotation_gff_path=args.annotation_gff,
            candidate_gff_path=args.candidate_gff,
            output_tsv_path=args.output_tsv,
            candidate_topology_tsv_path=args.candidate_topology_tsv,
        )
        print(json.dumps(manifest, indent=2, sort_keys=True))
        return 0
    raise AssertionError(f"Unhandled command: {args.command}")


if __name__ == "__main__":
    raise SystemExit(main())
