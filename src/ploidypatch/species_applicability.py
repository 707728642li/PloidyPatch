from __future__ import annotations

import csv
import hashlib
import json
import os
from pathlib import Path
from typing import Any, Mapping


SCHEMA_VERSION = "ploidypatch.species_applicability_report.v0.6"
METRICS_SCHEMA = "ploidypatch.species_applicability_metrics.v0.6"
POLICY_SCHEMA = "ploidypatch.species_applicability_policy.v0.6"
FORBIDDEN_KEY_PREFIXES = (
    "candidate_",
    "label_",
    "truth_",
    "ranking_",
    "performance_",
    "endpoint_",
)
FORBIDDEN_KEYS = {
    "candidates",
    "events",
    "labels",
    "positives",
    "truth_pairs",
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(8 * 1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def _load_json(path: Path) -> dict[str, Any]:
    if not path.is_file() or path.is_symlink() or path.stat().st_size == 0:
        raise ValueError(f"Missing, empty or symlinked applicability metrics: {path}")
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError("Applicability metrics must be a JSON object")
    return value


def _read_policy(path: Path) -> dict[str, str]:
    if not path.is_file() or path.is_symlink() or path.stat().st_size == 0:
        raise ValueError(f"Missing, empty or symlinked applicability policy: {path}")
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ["field", "value"]:
            raise ValueError("Applicability policy requires exact field/value columns")
        rows = list(reader)
    values = {row["field"]: row["value"] for row in rows}
    if len(values) != len(rows) or "" in values or any(value == "" for value in values.values()):
        raise ValueError("Applicability policy contains empty or duplicate fields")
    if values.get("schema_version") != POLICY_SCHEMA:
        raise ValueError("Applicability policy schema differs")
    return values


def _reject_forbidden_keys(value: Any, path: str = "metrics") -> None:
    if isinstance(value, Mapping):
        for raw_key, child in value.items():
            if not isinstance(raw_key, str):
                raise ValueError(f"Non-string metrics key at {path}")
            key = raw_key.lower()
            if key in FORBIDDEN_KEYS or key.startswith(FORBIDDEN_KEY_PREFIXES):
                raise ValueError(f"Forbidden post-candidate metric at {path}.{raw_key}")
            _reject_forbidden_keys(child, f"{path}.{raw_key}")
    elif isinstance(value, list):
        for index, child in enumerate(value):
            _reject_forbidden_keys(child, f"{path}[{index}]")


def _section(metrics: Mapping[str, Any], name: str, fields: set[str]) -> dict[str, Any]:
    value = metrics.get(name)
    if not isinstance(value, dict) or set(value) != fields:
        raise ValueError(f"Applicability {name} section must contain {sorted(fields)}")
    return value


def _fraction(value: Any, name: str) -> float:
    if isinstance(value, bool):
        raise ValueError(f"Applicability fraction is boolean: {name}")
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Applicability fraction is not numeric: {name}") from exc
    if not 0.0 <= number <= 1.0:
        raise ValueError(f"Applicability fraction is outside [0, 1]: {name}")
    return number


def _integer(value: Any, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError(f"Applicability count is not a non-negative integer: {name}")
    return value


def _boolean(value: Any, name: str) -> bool:
    if not isinstance(value, bool):
        raise ValueError(f"Applicability flag is not boolean: {name}")
    return value


def _minimum(policy: Mapping[str, str], key: str) -> float:
    try:
        value = float(policy[key])
    except (KeyError, ValueError) as exc:
        raise ValueError(f"Missing or invalid applicability threshold: {key}") from exc
    return value


def evaluate_species_applicability(
    *, metrics_path: str | Path, policy_path: str | Path, output_path: str | Path
) -> dict[str, Any]:
    metrics_file = Path(metrics_path)
    policy_file = Path(policy_path)
    output = Path(output_path)
    partial = Path(str(output) + ".partial")
    if output.exists() or partial.exists():
        raise FileExistsError("Refusing to overwrite species-applicability report")

    metrics = _load_json(metrics_file)
    policy = _read_policy(policy_file)
    _reject_forbidden_keys(metrics)
    if metrics.get("schema_version") != METRICS_SCHEMA:
        raise ValueError("Applicability metrics schema differs")
    species_id = metrics.get("species_id")
    if not isinstance(species_id, str) or not species_id.strip():
        raise ValueError("Applicability metrics require a non-empty species_id")
    controlled = _boolean(metrics.get("controlled_holdout"), "controlled_holdout")

    assembly = _section(
        metrics,
        "assembly",
        {
            "primary_seqid_count_matches_declared_karyotype",
            "primary_assembly_fraction",
            "primary_non_N_fraction",
            "assembly_BUSCO_complete_fraction",
            "read_backed_QV",
        },
    )
    annotation = _section(
        metrics,
        "annotation",
        {
            "primary_protein_coding_gene_fraction",
            "exact_unique_GFF_protein_mapping_fraction",
            "valid_coding_hierarchy_fraction",
            "valid_representative_translation_fraction",
            "protein_BUSCO_complete_fraction",
            "fuzzy_identifier_repairs",
        },
    )
    backbone = _section(
        metrics,
        "backbone",
        {
            "independent_WGD_source_count",
            "minimum_block_pairs",
            "cross_primary_seqid_only",
            "primary_gene_midpoint_backbone_coverage",
            "primary_chromosome_cell_coverage_fraction",
            "minimum_cells_per_covered_chromosome_observed",
            "unique_partner_chromosome_fraction_among_covered_genes",
            "input_permutation_deterministic",
            "reverse_duplicate_invariant",
            "built_without_candidates",
            "built_without_truth_or_labels",
            "used_perturbed_annotation_if_controlled",
        },
    )
    sources = metrics.get("source_manifests")
    if not isinstance(sources, dict) or not sources or any(
        not isinstance(key, str)
        or not key
        or not isinstance(value, str)
        or len(value) != 64
        for key, value in sources.items()
    ):
        raise ValueError("Applicability metrics require source manifest SHA-256 values")

    qv = assembly["read_backed_QV"]
    if qv is not None and (isinstance(qv, bool) or not isinstance(qv, (int, float))):
        raise ValueError("read_backed_QV must be numeric or null")
    assembly_gates: dict[str, bool | None] = {
        "primary_seqid_count_matches_declared_karyotype": _boolean(
            assembly["primary_seqid_count_matches_declared_karyotype"],
            "primary_seqid_count_matches_declared_karyotype",
        ),
        "primary_assembly_fraction": _fraction(
            assembly["primary_assembly_fraction"], "primary_assembly_fraction"
        )
        >= _minimum(policy, "minimum_primary_assembly_fraction"),
        "primary_non_N_fraction": _fraction(
            assembly["primary_non_N_fraction"], "primary_non_N_fraction"
        )
        >= _minimum(policy, "minimum_primary_non_N_fraction"),
        "assembly_BUSCO_complete_fraction": _fraction(
            assembly["assembly_BUSCO_complete_fraction"],
            "assembly_BUSCO_complete_fraction",
        )
        >= _minimum(policy, "minimum_assembly_BUSCO_complete_fraction"),
        "optional_read_backed_QV": (
            None
            if qv is None
            else float(qv) >= _minimum(policy, "minimum_optional_read_backed_QV")
        ),
    }
    annotation_gates = {
        "primary_protein_coding_gene_fraction": _fraction(
            annotation["primary_protein_coding_gene_fraction"],
            "primary_protein_coding_gene_fraction",
        )
        >= _minimum(policy, "minimum_primary_protein_coding_gene_fraction"),
        "exact_unique_GFF_protein_mapping_fraction": _fraction(
            annotation["exact_unique_GFF_protein_mapping_fraction"],
            "exact_unique_GFF_protein_mapping_fraction",
        )
        >= _minimum(policy, "minimum_exact_unique_GFF_protein_mapping_fraction"),
        "valid_coding_hierarchy_fraction": _fraction(
            annotation["valid_coding_hierarchy_fraction"],
            "valid_coding_hierarchy_fraction",
        )
        >= _minimum(policy, "minimum_valid_coding_hierarchy_fraction"),
        "valid_representative_translation_fraction": _fraction(
            annotation["valid_representative_translation_fraction"],
            "valid_representative_translation_fraction",
        )
        >= _minimum(policy, "minimum_valid_representative_translation_fraction"),
        "protein_BUSCO_complete_fraction": _fraction(
            annotation["protein_BUSCO_complete_fraction"],
            "protein_BUSCO_complete_fraction",
        )
        >= _minimum(policy, "minimum_protein_BUSCO_complete_fraction"),
        "no_fuzzy_identifier_repairs": _integer(
            annotation["fuzzy_identifier_repairs"], "fuzzy_identifier_repairs"
        )
        == 0,
    }
    backbone_gates = {
        "independent_WGD_sources": _integer(
            backbone["independent_WGD_source_count"],
            "independent_WGD_source_count",
        )
        >= int(_minimum(policy, "minimum_independent_WGD_sources")),
        "minimum_WGDI_block_pairs": _integer(
            backbone["minimum_block_pairs"], "minimum_block_pairs"
        )
        >= int(_minimum(policy, "minimum_WGDI_block_pairs")),
        "cross_primary_seqid_only": _boolean(
            backbone["cross_primary_seqid_only"], "cross_primary_seqid_only"
        ),
        "primary_gene_midpoint_backbone_coverage": _fraction(
            backbone["primary_gene_midpoint_backbone_coverage"],
            "primary_gene_midpoint_backbone_coverage",
        )
        >= _minimum(policy, "minimum_primary_gene_midpoint_backbone_coverage"),
        "primary_chromosome_cell_coverage_fraction": _fraction(
            backbone["primary_chromosome_cell_coverage_fraction"],
            "primary_chromosome_cell_coverage_fraction",
        )
        >= _minimum(policy, "minimum_primary_chromosome_cell_coverage_fraction"),
        "minimum_cells_per_covered_chromosome": _integer(
            backbone["minimum_cells_per_covered_chromosome_observed"],
            "minimum_cells_per_covered_chromosome_observed",
        )
        >= int(_minimum(policy, "minimum_directed_cells_per_covered_chromosome")),
        "unique_partner_chromosome_fraction": _fraction(
            backbone["unique_partner_chromosome_fraction_among_covered_genes"],
            "unique_partner_chromosome_fraction_among_covered_genes",
        )
        >= _minimum(
            policy,
            "minimum_unique_partner_chromosome_fraction_among_covered_genes",
        ),
        "input_permutation_deterministic": _boolean(
            backbone["input_permutation_deterministic"],
            "input_permutation_deterministic",
        ),
        "reverse_duplicate_invariant": _boolean(
            backbone["reverse_duplicate_invariant"], "reverse_duplicate_invariant"
        ),
        "built_without_candidates": _boolean(
            backbone["built_without_candidates"], "built_without_candidates"
        ),
        "built_without_truth_or_labels": _boolean(
            backbone["built_without_truth_or_labels"],
            "built_without_truth_or_labels",
        ),
        "controlled_annotation_is_perturbed": (
            not controlled
            or _boolean(
                backbone["used_perturbed_annotation_if_controlled"],
                "used_perturbed_annotation_if_controlled",
            )
        ),
    }

    input_quality_pass = all(value is not False for value in assembly_gates.values()) and all(
        annotation_gates.values()
    )
    backbone_pass = all(backbone_gates.values())
    if not input_quality_pass:
        state = "not_evaluable_input_quality"
    elif not backbone_pass:
        state = "chain_preservation_only"
    else:
        state = "full_topology_evaluable"

    report: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "policy_id": policy.get("policy_id"),
        "species_id": species_id,
        "controlled_holdout": controlled,
        "applicability_state": state,
        "input_quality_pass": input_quality_pass,
        "backbone_pass": backbone_pass,
        "gates": {
            "assembly": assembly_gates,
            "annotation": annotation_gates,
            "backbone": backbone_gates,
        },
        "metrics": {
            "assembly": assembly,
            "annotation": annotation,
            "backbone": backbone,
        },
        "source_manifests": dict(sorted(sources.items())),
        "inputs": {
            "metrics_sha256": _sha256(metrics_file),
            "policy_sha256": _sha256(policy_file),
        },
        "claim_boundary": {
            "candidate_statistics_used": False,
            "truth_or_labels_used": False,
            "backbone_failure_is_topology_success": False,
            "input_quality_failure_is_negative_biological_result": False,
            "automatic_annotation_approval": False,
        },
    }
    partial.parent.mkdir(parents=True, exist_ok=True)
    partial.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="",
    )
    os.replace(partial, output)
    return report

