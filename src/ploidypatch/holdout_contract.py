"""Strict declarative contracts for untouched external holdouts.

The contract is intentionally narrower than a general workflow configuration.
It describes one target, two candidate-only references and two evaluator-only
references while freezing the scientific rules that must not become tuning
knobs in a confirmatory replication.
"""
from __future__ import annotations

from dataclasses import dataclass, fields
import json
import math
from pathlib import Path, PurePosixPath
import re
from typing import Any, Iterable, Mapping


SCHEMA_VERSION = "ploidypatch.external_holdout_contract.v0.5"
MODEL_VERSION = "PloidyPatch_ranker_v0.4"
CORE_H1_MODEL_VERSION = "none_core_H1_only_v0.8"
KNOWN_SUBGENOME_CORE_H1_MODEL_VERSION = "none_core_H1_known_subgenome_v1.0"
CORE_H1_MODEL_VERSIONS = frozenset(
    {CORE_H1_MODEL_VERSION, KNOWN_SUBGENOME_CORE_H1_MODEL_VERSION}
)
ALLOWED_MODEL_VERSIONS = frozenset({MODEL_VERSION, *CORE_H1_MODEL_VERSIONS})

TARGET_ROLE = "target"
CANDIDATE_ROLE = "candidate_reference"
EVALUATOR_ROLE = "evaluator_reference"
REFERENCE_ROLES = frozenset({TARGET_ROLE, CANDIDATE_ROLE, EVALUATOR_ROLE})
ARTIFACT_NAMES = ("genome", "gff3", "protein")
ROLE_COUNTS = {TARGET_ROLE: 1, CANDIDATE_ROLE: 2, EVALUATOR_ROLE: 2}
ALLOWED_TEST_ROLES = frozenset(
    {
        "untouched_confirmatory_external_species",
        "target_level_predeclared_untouched_secondary_replication",
        "untouched_scale_stress_external_species",
    }
)

_ID_PATTERN = re.compile(r"[a-z0-9][a-z0-9._-]{2,127}")
_SPECIES_PATTERN = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,127}")
_RELEASE_PATTERN = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.+-]{0,159}")
_BUNDLE_PATTERN = re.compile(r"[a-z0-9][a-z0-9_-]{0,63}")
_PREFIX_PATTERN = re.compile(r"[a-z][a-z0-9]{1,7}")
_SHA256_PATTERN = re.compile(r"[0-9a-f]{64}")


@dataclass(frozen=True)
class TarMemberSource:
    """An exact regular member extracted from a tar-gzip source."""

    format: str
    member_name: PurePosixPath
    member_bytes: int
    member_sha256: str


@dataclass(frozen=True)
class ArtifactSource:
    """One immutable public input artifact."""

    source_relative_path: PurePosixPath
    bytes: int
    sha256: str
    container: TarMemberSource | None = None

    @property
    def staged_bytes(self) -> int:
        return self.container.member_bytes if self.container else self.bytes

    @property
    def staged_sha256(self) -> str:
        return self.container.member_sha256 if self.container else self.sha256

    @property
    def staged_filename(self) -> str:
        return (
            self.container.member_name.name
            if self.container
            else self.source_relative_path.name
        )


@dataclass(frozen=True)
class ReferenceContract:
    """One species and its three consistently role-bound artifacts."""

    role: str
    species_id: str
    release: str
    bundle_id: str
    wgdi_prefix: str
    primary_seqid_table: PurePosixPath
    genome: ArtifactSource
    gff3: ArtifactSource
    protein: ArtifactSource

    def artifact_items(self) -> tuple[tuple[str, ArtifactSource], ...]:
        return (
            ("genome", self.genome),
            ("gff3", self.gff3),
            ("protein", self.protein),
        )


@dataclass(frozen=True)
class HoldoutSeeds:
    truth_sampler: int
    h1_bootstrap: int
    h2_bootstrap: int
    guard_v03_bootstrap: int


@dataclass(frozen=True)
class CoreH1Seeds:
    truth_sampler: int
    h1_bootstrap: int


@dataclass(frozen=True)
class ScientificParameters:
    candidate_method_families: tuple[str, ...]
    multiple_references_per_method_vote: str
    candidate_topology_identity: str
    primary_candidate_policy: str
    legacy_candidate_comparator: str
    adapter_min_identity: float
    adapter_min_query_coverage: float
    adapter_require_intact: bool
    adapter_max_existing_cds_overlap: float
    adapter_max_redundancy_overlap: float
    truth_pair_self_wgdi_min_block_pairs: int
    truth_pair_self_wgdi_require_cross_seqid: bool
    truth_pair_self_wgdi_require_reciprocal_unique: bool
    truth_pair_outgroup_min_block_pairs: int
    truth_pair_outgroup_counterpart_multiplicity: str
    truth_pair_outgroup_min_support_groups: int
    truth_pair_outgroup_require_cross_seqid: bool
    truth_pair_outgroup_require_reciprocal_unique: bool
    truth_pair_final_rule: str
    truth_event_count: int
    minimum_formal_event_count: int
    minimum_events_per_complexity_bin: int
    truth_removal_policy: str
    truth_sampler_balance: str
    bootstrap_replicates: int
    minimum_chromosome_bootstrap_valid_replicates: int
    minimum_topology_coverage_among_positive_candidates: float
    minimum_v03_AP_gain_retained_fraction: float
    review_fraction_budgets: tuple[float, ...]
    review_absolute_budgets: tuple[int, ...]
    automatic_copy_addition_approval: bool
    score_interpretation: str


@dataclass(frozen=True)
class CoreH1ScientificParameters:
    """Exact no-ranker H1 profile for chain-preservation replication."""

    protocol_profile: str
    candidate_method_families: tuple[str, ...]
    multiple_references_per_method_vote: str
    candidate_topology_identity: str
    primary_candidate_policy: str
    legacy_candidate_comparator: str
    adapter_min_identity: float
    adapter_min_query_coverage: float
    adapter_require_intact: bool
    adapter_max_existing_cds_overlap: float
    adapter_max_redundancy_overlap: float
    truth_pair_self_wgdi_min_block_pairs: int
    truth_pair_self_wgdi_require_cross_seqid: bool
    truth_pair_self_wgdi_require_reciprocal_unique: bool
    truth_pair_yn00_ks_minimum: float
    truth_pair_yn00_ks_maximum: float
    truth_pair_missing_or_out_of_range_ks: str
    truth_pair_outgroup_min_block_pairs: int
    truth_pair_outgroup_counterpart_multiplicity: str
    truth_pair_outgroup_min_support_groups: int
    truth_pair_outgroup_require_cross_seqid: bool
    truth_pair_outgroup_require_reciprocal_unique: bool
    truth_pair_final_rule: str
    truth_event_count: int
    minimum_formal_event_count: int
    minimum_events_per_complexity_bin: int
    complexity_bins: tuple[str, ...]
    truth_removal_policy: str
    truth_sampler_balance: str
    primary_hypothesis: str
    h1_metric: str
    h1_success_gate: str
    bootstrap_replicates: int
    bootstrap_unit: str
    bootstrap_interval: str
    h2_or_topology_ranking: str
    retired_ranker: str
    no_op_and_oracle_sentinels: str
    restoration_sentinel: str
    genome_sha256_sentinel: str
    all_arm_collateral_loss_maximum: int
    automatic_copy_addition_approval: bool


@dataclass(frozen=True)
class KnownSubgenomeCoreH1ScientificParameters:
    """No-ranker H1 profile for a predeclared allopolyploid subgenome event."""

    protocol_profile: str
    candidate_method_families: tuple[str, ...]
    multiple_references_per_method_vote: str
    candidate_topology_identity: str
    primary_candidate_policy: str
    legacy_candidate_comparator: str
    adapter_min_identity: float
    adapter_min_query_coverage: float
    adapter_require_intact: bool
    adapter_max_existing_cds_overlap: float
    adapter_max_redundancy_overlap: float
    truth_pair_self_wgdi_min_block_pairs: int
    truth_pair_self_wgdi_require_cross_seqid: bool
    truth_pair_self_wgdi_require_reciprocal_unique: bool
    truth_pair_event_discriminator: str
    truth_pair_target_subgenome_pairing: str
    truth_pair_yn00_ks_policy: str
    truth_pair_outgroup_min_block_pairs: int
    truth_pair_outgroup_counterpart_multiplicity: str
    truth_pair_outgroup_min_support_groups: int
    truth_pair_outgroup_require_cross_seqid: bool
    truth_pair_outgroup_require_reciprocal_unique: bool
    truth_pair_final_rule: str
    truth_event_count: int
    minimum_formal_event_count: int
    minimum_events_per_complexity_bin: int
    complexity_bins: tuple[str, ...]
    truth_removal_policy: str
    truth_sampler_balance: str
    primary_hypothesis: str
    h1_metric: str
    h1_success_gate: str
    bootstrap_replicates: int
    bootstrap_unit: str
    bootstrap_interval: str
    h2_or_topology_ranking: str
    retired_ranker: str
    no_op_and_oracle_sentinels: str
    restoration_sentinel: str
    genome_sha256_sentinel: str
    all_arm_collateral_loss_maximum: int
    automatic_copy_addition_approval: bool


@dataclass(frozen=True)
class TargetResolvedParameters:
    """Target-size-dependent gates resolved and frozen before enumeration."""

    primary_chromosome_count: int
    minimum_target_chromosomes_fraction: float
    minimum_target_chromosomes: int


FIXED_SCIENTIFIC_PARAMETERS = ScientificParameters(
    candidate_method_families=(
        "miniprot_0.18-r281",
        "GeMoMa_1.9",
        "LiftOn_1.0.11",
    ),
    multiple_references_per_method_vote="one_method_family_vote",
    candidate_topology_identity="seqid+strand+exact_phased_CDS_chain",
    primary_candidate_policy="retain_distinct_phased_CDS_chains",
    legacy_candidate_comparator="suppress_strongly_overlapping_alternative_chains",
    adapter_min_identity=0.5,
    adapter_min_query_coverage=0.5,
    adapter_require_intact=True,
    adapter_max_existing_cds_overlap=0.2,
    adapter_max_redundancy_overlap=0.5,
    truth_pair_self_wgdi_min_block_pairs=20,
    truth_pair_self_wgdi_require_cross_seqid=True,
    truth_pair_self_wgdi_require_reciprocal_unique=True,
    truth_pair_outgroup_min_block_pairs=20,
    truth_pair_outgroup_counterpart_multiplicity="exactly_two_target_genes",
    truth_pair_outgroup_min_support_groups=2,
    truth_pair_outgroup_require_cross_seqid=True,
    truth_pair_outgroup_require_reciprocal_unique=True,
    truth_pair_final_rule=(
        "exact_unordered_pair_intersection_of_self_wgdi_and_two_outgroup_support"
    ),
    truth_event_count=800,
    minimum_formal_event_count=500,
    minimum_events_per_complexity_bin=20,
    truth_removal_policy=(
        "remove_exactly_one_deterministically_selected_partner_per_pair"
    ),
    truth_sampler_balance=(
        "global_chromosome_round_robin_then_four_CDS_complexity_bins_then_sha256_rank"
    ),
    bootstrap_replicates=20_000,
    minimum_chromosome_bootstrap_valid_replicates=19_000,
    minimum_topology_coverage_among_positive_candidates=0.70,
    minimum_v03_AP_gain_retained_fraction=0.90,
    review_fraction_budgets=(0.005, 0.01, 0.02),
    review_absolute_budgets=(100, 250, 500),
    automatic_copy_addition_approval=False,
    score_interpretation="review_priority_only",
)


FIXED_CORE_H1_SCIENTIFIC_PARAMETERS = CoreH1ScientificParameters(
    protocol_profile="core_H1_only_no_ranker",
    candidate_method_families=(
        "miniprot_0.18-r281",
        "GeMoMa_1.9",
        "LiftOn_1.0.11",
    ),
    multiple_references_per_method_vote="one_method_family_vote",
    candidate_topology_identity="seqid+strand+exact_phased_CDS_chain",
    primary_candidate_policy="retain_distinct_phased_CDS_chains",
    legacy_candidate_comparator="suppress_strongly_overlapping_alternative_chains",
    adapter_min_identity=0.5,
    adapter_min_query_coverage=0.5,
    adapter_require_intact=True,
    adapter_max_existing_cds_overlap=0.2,
    adapter_max_redundancy_overlap=0.5,
    truth_pair_self_wgdi_min_block_pairs=20,
    truth_pair_self_wgdi_require_cross_seqid=True,
    truth_pair_self_wgdi_require_reciprocal_unique=True,
    truth_pair_yn00_ks_minimum=0.10,
    truth_pair_yn00_ks_maximum=0.75,
    truth_pair_missing_or_out_of_range_ks="abstain",
    truth_pair_outgroup_min_block_pairs=20,
    truth_pair_outgroup_counterpart_multiplicity="exactly_two_target_genes",
    truth_pair_outgroup_min_support_groups=2,
    truth_pair_outgroup_require_cross_seqid=True,
    truth_pair_outgroup_require_reciprocal_unique=True,
    truth_pair_final_rule=(
        "exact_unordered_pair_intersection_of_ks_filtered_self_wgdi_and_"
        "two_outgroup_support"
    ),
    truth_event_count=800,
    minimum_formal_event_count=500,
    minimum_events_per_complexity_bin=20,
    complexity_bins=("one", "two_to_three", "four_to_six", "seven_plus"),
    truth_removal_policy=(
        "remove_exactly_one_deterministically_selected_partner_per_pair"
    ),
    truth_sampler_balance=(
        "global_chromosome_round_robin_then_four_CDS_complexity_bins_then_sha256_rank"
    ),
    primary_hypothesis="H1_retain_distinct_vs_suppress_overlap_only",
    h1_metric=(
        "event_exact_phased_CDS_recall_retain_distinct_minus_suppress_overlap"
    ),
    h1_success_gate="delta_gt_0_and_paired_bootstrap_CI_lower_gt_0",
    bootstrap_replicates=20_000,
    bootstrap_unit="paired_event",
    bootstrap_interval="two_sided_95pct_percentile",
    h2_or_topology_ranking="forbidden",
    retired_ranker="ploidypatch.stable_reference_ranker.v0.9",
    no_op_and_oracle_sentinels="required",
    restoration_sentinel="byte_identical_required",
    genome_sha256_sentinel="identical_required",
    all_arm_collateral_loss_maximum=0,
    automatic_copy_addition_approval=False,
)


FIXED_KNOWN_SUBGENOME_CORE_H1_SCIENTIFIC_PARAMETERS = (
    KnownSubgenomeCoreH1ScientificParameters(
        protocol_profile="core_H1_known_subgenome_no_ranker",
        candidate_method_families=(
            "miniprot_0.18-r281",
            "GeMoMa_1.9",
            "LiftOn_1.0.11",
        ),
        multiple_references_per_method_vote="one_method_family_vote",
        candidate_topology_identity="seqid+strand+exact_phased_CDS_chain",
        primary_candidate_policy="retain_distinct_phased_CDS_chains",
        legacy_candidate_comparator=(
            "suppress_strongly_overlapping_alternative_chains"
        ),
        adapter_min_identity=0.5,
        adapter_min_query_coverage=0.5,
        adapter_require_intact=True,
        adapter_max_existing_cds_overlap=0.2,
        adapter_max_redundancy_overlap=0.5,
        truth_pair_self_wgdi_min_block_pairs=20,
        truth_pair_self_wgdi_require_cross_seqid=True,
        truth_pair_self_wgdi_require_reciprocal_unique=True,
        truth_pair_event_discriminator=(
            "predeclared_homoeolog_group_and_subgenome_labels"
        ),
        truth_pair_target_subgenome_pairing=(
            "same_group_exactly_one_member_from_each_subgenome"
        ),
        truth_pair_yn00_ks_policy="descriptive_only_not_used_for_selection",
        truth_pair_outgroup_min_block_pairs=20,
        truth_pair_outgroup_counterpart_multiplicity=(
            "exactly_two_target_genes"
        ),
        truth_pair_outgroup_min_support_groups=2,
        truth_pair_outgroup_require_cross_seqid=True,
        truth_pair_outgroup_require_reciprocal_unique=True,
        truth_pair_final_rule=(
            "exact_unordered_pair_intersection_of_predeclared_subgenome_"
            "self_wgdi_and_two_outgroup_support"
        ),
        truth_event_count=800,
        minimum_formal_event_count=500,
        minimum_events_per_complexity_bin=20,
        complexity_bins=("one", "two_to_three", "four_to_six", "seven_plus"),
        truth_removal_policy=(
            "remove_exactly_one_deterministically_selected_partner_per_pair"
        ),
        truth_sampler_balance=(
            "global_chromosome_round_robin_then_four_CDS_complexity_bins_"
            "then_sha256_rank"
        ),
        primary_hypothesis="H1_retain_distinct_vs_suppress_overlap_only",
        h1_metric=(
            "event_exact_phased_CDS_recall_retain_distinct_minus_suppress_overlap"
        ),
        h1_success_gate="delta_gt_0_and_paired_bootstrap_CI_lower_gt_0",
        bootstrap_replicates=20_000,
        bootstrap_unit="paired_event",
        bootstrap_interval="two_sided_95pct_percentile",
        h2_or_topology_ranking="forbidden",
        retired_ranker="ploidypatch.stable_reference_ranker.v0.9",
        no_op_and_oracle_sentinels="required",
        restoration_sentinel="byte_identical_required",
        genome_sha256_sentinel="identical_required",
        all_arm_collateral_loss_maximum=0,
        automatic_copy_addition_approval=False,
    )
)


TRUTH_BLIND_DECLARATIONS: Mapping[str, bool] = {
    "selection_by_pair_yield_candidate_count_or_model_performance": False,
    "reference_role_change_after_freeze": False,
    "wgd_pairs_enumerated_before_protocol_freeze": False,
    "candidate_counts_computed_before_protocol_freeze": False,
    "truth_labels_accessed_before_protocol_freeze": False,
    "candidate_truth_access": False,
    "candidate_evaluator_reference_access": False,
    "blind_truth_mount": False,
    "blind_complete_target_annotation_mount": False,
    "blind_evaluator_reference_mount": False,
    "blind_nas_data_mount": False,
    "blind_network_access": False,
    "complete_control_generated_after_blind_raw_prediction_freeze": True,
    "automatic_copy_addition_approval": False,
}


@dataclass(frozen=True)
class HoldoutContract:
    schema_version: str
    holdout_id: str
    policy_id: str
    test_role: str
    model_version: str
    references: tuple[ReferenceContract, ...]
    seeds: HoldoutSeeds | CoreH1Seeds
    target_resolved_parameters: TargetResolvedParameters
    scientific_parameters: (
        ScientificParameters
        | CoreH1ScientificParameters
        | KnownSubgenomeCoreH1ScientificParameters
    )
    truth_blind: Mapping[str, bool]

    def references_for_role(self, role: str) -> tuple[ReferenceContract, ...]:
        if role not in REFERENCE_ROLES:
            raise ValueError(f"Unknown holdout reference role: {role}")
        return tuple(reference for reference in self.references if reference.role == role)

    @property
    def target(self) -> ReferenceContract:
        return self.references_for_role(TARGET_ROLE)[0]


def _reject_duplicate_keys(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"Duplicate JSON object key: {key}")
        result[key] = value
    return result


def _object(value: Any, context: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise ValueError(f"{context} must be a JSON object")
    return value


def _exact_keys(value: Mapping[str, Any], expected: Iterable[str], context: str) -> None:
    expected_set = set(expected)
    actual = set(value)
    if actual != expected_set:
        missing = sorted(expected_set - actual)
        extra = sorted(actual - expected_set)
        raise ValueError(
            f"{context} fields differ; missing={missing}, extra={extra}"
        )


def _string(value: Any, context: str, pattern: re.Pattern[str] | None = None) -> str:
    if not isinstance(value, str) or not value:
        raise ValueError(f"{context} must be a non-empty string")
    if pattern is not None and pattern.fullmatch(value) is None:
        raise ValueError(f"Unsafe or malformed {context}: {value!r}")
    return value


def _integer(value: Any, context: str, *, minimum: int = 0) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise ValueError(f"{context} must be an integer >= {minimum}")
    return value


def _number(value: Any, context: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{context} must be numeric")
    return float(value)


def safe_relative_path(value: Any, context: str) -> PurePosixPath:
    """Return a canonical, traversal-free POSIX relative path."""

    raw = _string(value, context)
    if (
        "\\" in raw
        or "\x00" in raw
        or raw.startswith("/")
        or raw.endswith("/")
        or "//" in raw
        or ":" in raw
    ):
        raise ValueError(f"Unsafe {context}: {raw!r}")
    parts = raw.split("/")
    if any(part in {"", ".", ".."} for part in parts):
        raise ValueError(f"Unsafe {context}: {raw!r}")
    path = PurePosixPath(*parts)
    if path.is_absolute() or path.as_posix() != raw:
        raise ValueError(f"Non-canonical {context}: {raw!r}")
    return path


def _load_artifact(
    value: Any, context: str, *, allow_container: bool = False
) -> ArtifactSource:
    item = _object(value, context)
    base_keys = {"source_relative_path", "bytes", "sha256"}
    if set(item) not in (base_keys, base_keys | {"container"}):
        expected = base_keys | ({"container"} if "container" in item else set())
        _exact_keys(item, expected, context)
    digest = _string(item["sha256"], f"{context}.sha256", _SHA256_PATTERN)
    source_path = safe_relative_path(
        item["source_relative_path"], f"{context}.source_relative_path"
    )
    container: TarMemberSource | None = None
    if "container" in item:
        if not allow_container:
            raise ValueError("Only a genome artifact may use a tar container")
        container_item = _object(item["container"], f"{context}.container")
        _exact_keys(
            container_item,
            ("format", "member_name", "member_bytes", "member_sha256"),
            f"{context}.container",
        )
        container_format = _string(
            container_item["format"], f"{context}.container.format"
        )
        if container_format != "tar.gz" or not source_path.name.endswith(".tar.gz"):
            raise ValueError(f"{context}.container must bind a .tar.gz source")
        member_name = safe_relative_path(
            container_item["member_name"], f"{context}.container.member_name"
        )
        if not member_name.name.lower().endswith((".fa", ".fasta", ".fna")):
            raise ValueError(f"{context}.container member must be a FASTA file")
        container = TarMemberSource(
            format=container_format,
            member_name=member_name,
            member_bytes=_integer(
                container_item["member_bytes"],
                f"{context}.container.member_bytes",
                minimum=1,
            ),
            member_sha256=_string(
                container_item["member_sha256"],
                f"{context}.container.member_sha256",
                _SHA256_PATTERN,
            ),
        )
    return ArtifactSource(
        source_relative_path=source_path,
        bytes=_integer(item["bytes"], f"{context}.bytes", minimum=1),
        sha256=digest,
        container=container,
    )


def _load_reference(value: Any, index: int) -> ReferenceContract:
    context = f"references[{index}]"
    item = _object(value, context)
    _exact_keys(
        item,
        (
            "role",
            "species_id",
            "release",
            "bundle_id",
            "wgdi_prefix",
            "primary_seqid_table",
            "artifacts",
        ),
        context,
    )
    role = _string(item["role"], f"{context}.role")
    if role not in REFERENCE_ROLES:
        raise ValueError(f"Unknown {context}.role: {role}")
    primary = safe_relative_path(
        item["primary_seqid_table"], f"{context}.primary_seqid_table"
    )
    if primary.suffix != ".tsv":
        raise ValueError(f"{context}.primary_seqid_table must end in .tsv")
    artifacts = _object(item["artifacts"], f"{context}.artifacts")
    _exact_keys(artifacts, ARTIFACT_NAMES, f"{context}.artifacts")
    return ReferenceContract(
        role=role,
        species_id=_string(
            item["species_id"], f"{context}.species_id", _SPECIES_PATTERN
        ),
        release=_string(item["release"], f"{context}.release", _RELEASE_PATTERN),
        bundle_id=_string(
            item["bundle_id"], f"{context}.bundle_id", _BUNDLE_PATTERN
        ),
        wgdi_prefix=_string(
            item["wgdi_prefix"], f"{context}.wgdi_prefix", _PREFIX_PATTERN
        ),
        primary_seqid_table=primary,
        genome=_load_artifact(
            artifacts["genome"], f"{context}.artifacts.genome", allow_container=True
        ),
        gff3=_load_artifact(artifacts["gff3"], f"{context}.artifacts.gff3"),
        protein=_load_artifact(
            artifacts["protein"], f"{context}.artifacts.protein"
        ),
    )


def _load_seeds(
    value: Any, model_version: str
) -> HoldoutSeeds | CoreH1Seeds:
    item = _object(value, "seeds")
    seed_type = CoreH1Seeds if model_version in CORE_H1_MODEL_VERSIONS else HoldoutSeeds
    names = tuple(field.name for field in fields(seed_type))
    _exact_keys(item, names, "seeds")
    values = {
        name: _integer(item[name], f"seeds.{name}", minimum=1) for name in names
    }
    if any(seed > (2**63 - 1) for seed in values.values()):
        raise ValueError("Holdout seeds must fit a signed 64-bit integer")
    if len(set(values.values())) != len(values):
        raise ValueError("Every holdout randomization stream requires a unique seed")
    return seed_type(**values)


def _load_target_resolved_parameters(value: Any) -> TargetResolvedParameters:
    item = _object(value, "target_resolved_parameters")
    names = tuple(field.name for field in fields(TargetResolvedParameters))
    _exact_keys(item, names, "target_resolved_parameters")
    count = _integer(
        item["primary_chromosome_count"],
        "target_resolved_parameters.primary_chromosome_count",
        minimum=1,
    )
    fraction = _number(
        item["minimum_target_chromosomes_fraction"],
        "target_resolved_parameters.minimum_target_chromosomes_fraction",
    )
    minimum = _integer(
        item["minimum_target_chromosomes"],
        "target_resolved_parameters.minimum_target_chromosomes",
        minimum=1,
    )
    if fraction != 0.75:
        raise ValueError(
            "minimum_target_chromosomes_fraction differs from the frozen 0.75 rule"
        )
    expected = math.ceil(fraction * count)
    if minimum != expected:
        raise ValueError(
            "minimum_target_chromosomes must equal "
            f"ceil({fraction} * {count}) = {expected}, observed={minimum}"
        )
    return TargetResolvedParameters(
        primary_chromosome_count=count,
        minimum_target_chromosomes_fraction=fraction,
        minimum_target_chromosomes=minimum,
    )


def _load_scientific_parameters(value: Any) -> ScientificParameters:
    item = _object(value, "scientific_parameters")
    names = tuple(field.name for field in fields(ScientificParameters))
    _exact_keys(item, names, "scientific_parameters")
    sequence_strings = item["candidate_method_families"]
    if not isinstance(sequence_strings, list) or not all(
        isinstance(entry, str) and entry for entry in sequence_strings
    ):
        raise ValueError("candidate_method_families must be a non-empty string list")
    fractions = item["review_fraction_budgets"]
    absolutes = item["review_absolute_budgets"]
    if not isinstance(fractions, list) or not isinstance(absolutes, list):
        raise ValueError("Review budgets must be JSON arrays")
    parameters = ScientificParameters(
        candidate_method_families=tuple(sequence_strings),
        multiple_references_per_method_vote=_string(
            item["multiple_references_per_method_vote"],
            "scientific_parameters.multiple_references_per_method_vote",
        ),
        candidate_topology_identity=_string(
            item["candidate_topology_identity"],
            "scientific_parameters.candidate_topology_identity",
        ),
        primary_candidate_policy=_string(
            item["primary_candidate_policy"],
            "scientific_parameters.primary_candidate_policy",
        ),
        legacy_candidate_comparator=_string(
            item["legacy_candidate_comparator"],
            "scientific_parameters.legacy_candidate_comparator",
        ),
        adapter_min_identity=_number(
            item["adapter_min_identity"], "scientific_parameters.adapter_min_identity"
        ),
        adapter_min_query_coverage=_number(
            item["adapter_min_query_coverage"],
            "scientific_parameters.adapter_min_query_coverage",
        ),
        adapter_require_intact=item["adapter_require_intact"],
        adapter_max_existing_cds_overlap=_number(
            item["adapter_max_existing_cds_overlap"],
            "scientific_parameters.adapter_max_existing_cds_overlap",
        ),
        adapter_max_redundancy_overlap=_number(
            item["adapter_max_redundancy_overlap"],
            "scientific_parameters.adapter_max_redundancy_overlap",
        ),
        truth_pair_self_wgdi_min_block_pairs=_integer(
            item["truth_pair_self_wgdi_min_block_pairs"],
            "scientific_parameters.truth_pair_self_wgdi_min_block_pairs",
        ),
        truth_pair_self_wgdi_require_cross_seqid=item[
            "truth_pair_self_wgdi_require_cross_seqid"
        ],
        truth_pair_self_wgdi_require_reciprocal_unique=item[
            "truth_pair_self_wgdi_require_reciprocal_unique"
        ],
        truth_pair_outgroup_min_block_pairs=_integer(
            item["truth_pair_outgroup_min_block_pairs"],
            "scientific_parameters.truth_pair_outgroup_min_block_pairs",
        ),
        truth_pair_outgroup_counterpart_multiplicity=_string(
            item["truth_pair_outgroup_counterpart_multiplicity"],
            "scientific_parameters.truth_pair_outgroup_counterpart_multiplicity",
        ),
        truth_pair_outgroup_min_support_groups=_integer(
            item["truth_pair_outgroup_min_support_groups"],
            "scientific_parameters.truth_pair_outgroup_min_support_groups",
        ),
        truth_pair_outgroup_require_cross_seqid=item[
            "truth_pair_outgroup_require_cross_seqid"
        ],
        truth_pair_outgroup_require_reciprocal_unique=item[
            "truth_pair_outgroup_require_reciprocal_unique"
        ],
        truth_pair_final_rule=_string(
            item["truth_pair_final_rule"],
            "scientific_parameters.truth_pair_final_rule",
        ),
        truth_event_count=_integer(
            item["truth_event_count"], "scientific_parameters.truth_event_count"
        ),
        minimum_formal_event_count=_integer(
            item["minimum_formal_event_count"],
            "scientific_parameters.minimum_formal_event_count",
        ),
        minimum_events_per_complexity_bin=_integer(
            item["minimum_events_per_complexity_bin"],
            "scientific_parameters.minimum_events_per_complexity_bin",
        ),
        truth_removal_policy=_string(
            item["truth_removal_policy"],
            "scientific_parameters.truth_removal_policy",
        ),
        truth_sampler_balance=_string(
            item["truth_sampler_balance"],
            "scientific_parameters.truth_sampler_balance",
        ),
        bootstrap_replicates=_integer(
            item["bootstrap_replicates"],
            "scientific_parameters.bootstrap_replicates",
        ),
        minimum_chromosome_bootstrap_valid_replicates=_integer(
            item["minimum_chromosome_bootstrap_valid_replicates"],
            "scientific_parameters.minimum_chromosome_bootstrap_valid_replicates",
        ),
        minimum_topology_coverage_among_positive_candidates=_number(
            item["minimum_topology_coverage_among_positive_candidates"],
            "scientific_parameters.minimum_topology_coverage_among_positive_candidates",
        ),
        minimum_v03_AP_gain_retained_fraction=_number(
            item["minimum_v03_AP_gain_retained_fraction"],
            "scientific_parameters.minimum_v03_AP_gain_retained_fraction",
        ),
        review_fraction_budgets=tuple(
            _number(value, "scientific_parameters.review_fraction_budgets[]")
            for value in fractions
        ),
        review_absolute_budgets=tuple(
            _integer(
                value,
                "scientific_parameters.review_absolute_budgets[]",
                minimum=1,
            )
            for value in absolutes
        ),
        automatic_copy_addition_approval=item[
            "automatic_copy_addition_approval"
        ],
        score_interpretation=_string(
            item["score_interpretation"],
            "scientific_parameters.score_interpretation",
        ),
    )
    boolean_names = (
        "adapter_require_intact",
        "truth_pair_self_wgdi_require_cross_seqid",
        "truth_pair_self_wgdi_require_reciprocal_unique",
        "truth_pair_outgroup_require_cross_seqid",
        "truth_pair_outgroup_require_reciprocal_unique",
        "automatic_copy_addition_approval",
    )
    if any(type(getattr(parameters, name)) is not bool for name in boolean_names):
        raise ValueError("Scientific boolean parameters must be JSON booleans")
    if parameters != FIXED_SCIENTIFIC_PARAMETERS:
        changed = [
            field.name
            for field in fields(ScientificParameters)
            if getattr(parameters, field.name)
            != getattr(FIXED_SCIENTIFIC_PARAMETERS, field.name)
        ]
        raise ValueError(
            "Scientific holdout parameters differ from frozen v0.5 rules: "
            + ", ".join(changed)
        )
    return parameters


def _load_core_h1_scientific_parameters(
    value: Any,
) -> CoreH1ScientificParameters:
    item = _object(value, "scientific_parameters")
    names = tuple(field.name for field in fields(CoreH1ScientificParameters))
    _exact_keys(item, names, "scientific_parameters")
    values: dict[str, Any] = {}
    for definition in fields(CoreH1ScientificParameters):
        name = definition.name
        expected = getattr(FIXED_CORE_H1_SCIENTIFIC_PARAMETERS, name)
        raw = item[name]
        context = f"scientific_parameters.{name}"
        if isinstance(expected, tuple):
            if not isinstance(raw, list) or not raw or not all(
                isinstance(entry, str) and entry for entry in raw
            ):
                raise ValueError(f"{context} must be a non-empty string list")
            values[name] = tuple(raw)
        elif isinstance(expected, bool):
            if type(raw) is not bool:
                raise ValueError(f"{context} must be a JSON boolean")
            values[name] = raw
        elif isinstance(expected, int):
            values[name] = _integer(raw, context)
        elif isinstance(expected, float):
            values[name] = _number(raw, context)
        else:
            values[name] = _string(raw, context)
    parameters = CoreH1ScientificParameters(**values)
    if parameters != FIXED_CORE_H1_SCIENTIFIC_PARAMETERS:
        changed = [
            definition.name
            for definition in fields(CoreH1ScientificParameters)
            if getattr(parameters, definition.name)
            != getattr(FIXED_CORE_H1_SCIENTIFIC_PARAMETERS, definition.name)
        ]
        raise ValueError(
            "Core H1-only parameters differ from frozen rules: "
            + ", ".join(changed)
        )
    return parameters


def _load_known_subgenome_core_h1_scientific_parameters(
    value: Any,
) -> KnownSubgenomeCoreH1ScientificParameters:
    item = _object(value, "scientific_parameters")
    parameter_type = KnownSubgenomeCoreH1ScientificParameters
    frozen = FIXED_KNOWN_SUBGENOME_CORE_H1_SCIENTIFIC_PARAMETERS
    names = tuple(definition.name for definition in fields(parameter_type))
    _exact_keys(item, names, "scientific_parameters")
    values: dict[str, Any] = {}
    for definition in fields(parameter_type):
        name = definition.name
        expected = getattr(frozen, name)
        raw = item[name]
        context = f"scientific_parameters.{name}"
        if isinstance(expected, tuple):
            if not isinstance(raw, list) or not raw or not all(
                isinstance(entry, str) and entry for entry in raw
            ):
                raise ValueError(f"{context} must be a non-empty string list")
            values[name] = tuple(raw)
        elif isinstance(expected, bool):
            if type(raw) is not bool:
                raise ValueError(f"{context} must be a JSON boolean")
            values[name] = raw
        elif isinstance(expected, int):
            values[name] = _integer(raw, context)
        elif isinstance(expected, float):
            values[name] = _number(raw, context)
        else:
            values[name] = _string(raw, context)
    parameters = parameter_type(**values)
    if parameters != frozen:
        changed = [
            definition.name
            for definition in fields(parameter_type)
            if getattr(parameters, definition.name) != getattr(frozen, definition.name)
        ]
        raise ValueError(
            "Known-subgenome core H1-only parameters differ from frozen rules: "
            + ", ".join(changed)
        )
    return parameters


def _load_truth_blind(value: Any) -> Mapping[str, bool]:
    item = _object(value, "truth_blind")
    _exact_keys(item, TRUTH_BLIND_DECLARATIONS, "truth_blind")
    if any(type(value) is not bool for value in item.values()):
        raise ValueError("truth_blind declarations must be JSON booleans")
    changed = [
        key
        for key, expected in TRUTH_BLIND_DECLARATIONS.items()
        if item[key] is not expected
    ]
    if changed:
        raise ValueError(
            "Truth-blind declarations differ from the fail-closed contract: "
            + ", ".join(changed)
        )
    return dict(item)


def staged_relative_path(
    reference: ReferenceContract, artifact_name: str
) -> PurePosixPath:
    """Return the only permitted staged role path for one source artifact."""

    if artifact_name not in ARTIFACT_NAMES:
        raise ValueError(f"Unknown holdout artifact: {artifact_name}")
    source = dict(reference.artifact_items())[artifact_name]
    filename = source.staged_filename
    species = reference.species_id
    if reference.role == TARGET_ROLE and artifact_name == "genome":
        root = PurePosixPath("shared_target")
    elif reference.role == TARGET_ROLE:
        root = PurePosixPath("evaluator_only/target_complete")
    elif reference.role == CANDIDATE_ROLE:
        root = PurePosixPath("candidate_only")
    elif reference.role == EVALUATOR_ROLE:
        root = PurePosixPath("evaluator_only/truth_references")
    else:  # pragma: no cover - ReferenceContract is constructed by the loader.
        raise ValueError(f"Unknown reference role: {reference.role}")
    return root / species / filename


def validate_holdout_contract(contract: HoldoutContract) -> None:
    """Validate cross-reference invariants not expressible per JSON object."""

    counts = {
        role: sum(reference.role == role for reference in contract.references)
        for role in REFERENCE_ROLES
    }
    if counts != ROLE_COUNTS:
        raise ValueError(
            "A holdout requires exactly one target, two candidate references and "
            f"two evaluator references; observed={counts}"
        )
    for label, values in (
        ("species_id", [reference.species_id for reference in contract.references]),
        ("bundle_id", [reference.bundle_id for reference in contract.references]),
        ("wgdi_prefix", [reference.wgdi_prefix for reference in contract.references]),
        (
            "primary_seqid_table",
            [reference.primary_seqid_table.as_posix() for reference in contract.references],
        ),
    ):
        folded = [str(value).casefold() for value in values]
        if len(folded) != len(set(folded)):
            raise ValueError(f"Holdout {label} values must be unique across roles")
    species = {reference.species_id.casefold() for reference in contract.references}
    prefixes = {reference.wgdi_prefix.casefold() for reference in contract.references}
    if species & prefixes:
        raise ValueError("Species IDs and WGDI prefixes must use disjoint identifiers")

    source_paths: list[str] = []
    staged_paths: list[str] = []
    for reference in contract.references:
        for artifact_name, artifact in reference.artifact_items():
            source_paths.append(artifact.source_relative_path.as_posix().casefold())
            staged_paths.append(
                staged_relative_path(reference, artifact_name).as_posix().casefold()
            )
    if len(source_paths) != len(set(source_paths)):
        raise ValueError("Every holdout artifact requires a distinct source path")
    if len(staged_paths) != len(set(staged_paths)):
        raise ValueError("Two holdout artifacts map to the same staged path")


def load_holdout_contract(path: str | Path) -> HoldoutContract:
    """Load a contract JSON with duplicate-key and unknown-field rejection."""

    contract_path = Path(path)
    if (
        not contract_path.is_file()
        or contract_path.is_symlink()
        or contract_path.stat().st_size == 0
    ):
        raise ValueError(f"Missing, empty or symlinked holdout contract: {contract_path}")
    try:
        value = json.loads(
            contract_path.read_text(encoding="utf-8"),
            object_pairs_hook=_reject_duplicate_keys,
        )
    except json.JSONDecodeError as error:
        raise ValueError(f"Malformed holdout contract JSON: {contract_path}") from error
    item = _object(value, "contract")
    _exact_keys(
        item,
        (
            "schema_version",
            "holdout_id",
            "policy_id",
            "test_role",
            "model_version",
            "references",
            "seeds",
            "target_resolved_parameters",
            "scientific_parameters",
            "truth_blind",
        ),
        "contract",
    )
    if item["schema_version"] != SCHEMA_VERSION:
        raise ValueError(f"Unsupported holdout contract schema: {item['schema_version']!r}")
    model_version = _string(item["model_version"], "model_version")
    if model_version not in ALLOWED_MODEL_VERSIONS:
        raise ValueError("Unsupported holdout model/profile version")
    test_role = _string(item["test_role"], "test_role")
    if test_role not in ALLOWED_TEST_ROLES:
        raise ValueError(f"Unsupported untouched holdout test role: {test_role}")
    raw_references = item["references"]
    if not isinstance(raw_references, list):
        raise ValueError("references must be a JSON array")
    references = tuple(
        _load_reference(reference, index)
        for index, reference in enumerate(raw_references)
    )
    contract = HoldoutContract(
        schema_version=SCHEMA_VERSION,
        holdout_id=_string(item["holdout_id"], "holdout_id", _ID_PATTERN),
        policy_id=_string(item["policy_id"], "policy_id", _ID_PATTERN),
        test_role=test_role,
        model_version=model_version,
        references=references,
        seeds=_load_seeds(item["seeds"], model_version),
        target_resolved_parameters=_load_target_resolved_parameters(
            item["target_resolved_parameters"]
        ),
        scientific_parameters=(
            _load_core_h1_scientific_parameters(item["scientific_parameters"])
            if model_version == CORE_H1_MODEL_VERSION
            else _load_known_subgenome_core_h1_scientific_parameters(
                item["scientific_parameters"]
            )
            if model_version == KNOWN_SUBGENOME_CORE_H1_MODEL_VERSION
            else _load_scientific_parameters(item["scientific_parameters"])
        ),
        truth_blind=_load_truth_blind(item["truth_blind"]),
    )
    validate_holdout_contract(contract)
    return contract
