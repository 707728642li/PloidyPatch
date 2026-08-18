from __future__ import annotations

import csv
import hashlib
import json
import os
from collections import Counter, defaultdict
from pathlib import Path
import shutil
import tempfile
from typing import Any, Iterable, Mapping

from .artifact_manifest import sha256_file, verify_sha256sums, write_sha256sums
from .audit import PROTEIN_ALPHABET
from .bootstrap import paired_event_bootstrap
from .io import iter_fasta, normalize_feature_id, parse_fasta_header_fields
from .perturb import read_gff_document
from .synteny_io import _normalized_gene_records


HOLDOUT_ID = "walnut_walnut2_v0.8"
POLICY_ID = "ploidypatch_walnut_external_core_h1_v0.8"
RETAIN_POOL_SCHEMA = "ploidypatch.method_candidate_pool.v2"
SUPPRESS_POOL_SCHEMA = "ploidypatch.method_consensus.v1"
POOL_SCHEMAS = {
    "retain_distinct_phased_CDS_chains": RETAIN_POOL_SCHEMA,
    "suppress_strongly_overlapping_alternative_chains": SUPPRESS_POOL_SCHEMA,
}
RAW_SCHEMA = "ploidypatch.walnut_h1_raw_predictions.v0.8"
PAIR_SCHEMA = "ploidypatch.walnut_exact_pair_truth.v0.8"
EVALUATION_SCHEMA = "ploidypatch.walnut_external_core_h1_evaluation.v0.8"
STATUS_SCHEMA = "ploidypatch.walnut_h1_reveal_status.v0.8"
FORMAL_STATUSES = frozenset({"ready", "not_evaluable", "invalid"})
KS_MINIMUM = 0.10
KS_MAXIMUM = 0.75
MIN_BLOCK_PAIRS = 20
BOOTSTRAP_REPLICATES = 20_000
BOOTSTRAP_SEED = 20260822
TARGET_EVENT = "Juglandaceae_ancestral_juglandoid_WGD"


def score_collateral_gate(score: Mapping[str, Any]) -> bool:
    """Apply the core-H1 zero-collateral quality gate without ranker imports."""

    return (
        score.get("quality_gate", {}).get("grade") == "pass"
        and score.get("collateral_changes", {}).get(
            "baseline_transcript_structures_missing_from_candidate"
        )
        == 0
    )


def write_primary_candidate_proteins(
    *,
    primary_gff_path: str | Path,
    provider_protein_path: str | Path,
    output_fasta_path: str | Path,
    whitelist_tsv_path: str | Path,
    manifest_path: str | Path,
) -> dict[str, Any]:
    """Write one exact, header-aware provider protein per primary coding gene."""

    output = Path(output_fasta_path)
    whitelist = Path(whitelist_tsv_path)
    manifest_output = Path(manifest_path)
    if any(path.exists() or path.is_symlink() for path in (output, whitelist, manifest_output)):
        raise FileExistsError("Refusing to overwrite Walnut primary protein artifacts")
    document = read_gff_document(primary_gff_path)
    genes, transcript_to_gene, hierarchy = _normalized_gene_records(document.records)
    coding_genes: set[str] = set()
    cds_protein_ids: set[str] = set()
    protein_to_genes: dict[str, set[str]] = defaultdict(set)
    unmapped_cds = 0
    for record in document.records:
        if record.feature_type != "CDS":
            continue
        protein_id = record.attributes.get("protein_id")
        if protein_id:
            cds_protein_ids.add(normalize_feature_id(protein_id))
        parents = {normalize_feature_id(parent) for parent in record.parents}
        mapped = {
            transcript_to_gene[parent]
            for parent in parents
            if parent in transcript_to_gene
        } | {parent for parent in parents if parent in genes}
        if len(mapped) != 1:
            unmapped_cds += 1
        else:
            gene_id = next(iter(mapped))
            coding_genes.add(gene_id)
            if protein_id:
                protein_to_genes[normalize_feature_id(protein_id)].add(gene_id)
    if not coding_genes or unmapped_cds:
        raise ValueError(
            f"Primary CDS hierarchy is not uniquely gene-resolved: genes={len(coding_genes)}, unmapped={unmapped_cds}"
        )

    ambiguous_protein_ids = {
        protein_id for protein_id, mapped in protein_to_genes.items() if len(mapped) != 1
    }
    if ambiguous_protein_ids:
        raise ValueError("Primary CDS protein_id maps to multiple genes")
    protein_to_gene = {
        protein_id: next(iter(mapped)) for protein_id, mapped in protein_to_genes.items()
    }
    provider_ids: set[str] = set()
    provider_records = 0
    excluded = 0
    relation_sources: Counter[str] = Counter()
    representatives: dict[str, tuple[str, str, str]] = {}
    for protein_id, header, raw_sequence in iter_fasta(provider_protein_path):
        provider_records += 1
        if protein_id in provider_ids:
            raise ValueError(f"Duplicate provider protein identifier: {protein_id}")
        provider_ids.add(protein_id)
        fields = parse_fasta_header_fields(header)
        candidates: dict[str, str] = {}
        explicit_gene = fields.get("gene") or fields.get("gene_id")
        if explicit_gene:
            gene_id = normalize_feature_id(explicit_gene)
            if gene_id in genes:
                candidates[gene_id] = "header:gene"
        explicit_transcript = fields.get("transcript") or fields.get("transcript_id")
        if explicit_transcript:
            transcript_id = normalize_feature_id(explicit_transcript)
            if transcript_id in transcript_to_gene:
                candidates[transcript_to_gene[transcript_id]] = "header:transcript->gff:Parent"
        if protein_id in protein_to_gene:
            candidates[protein_to_gene[protein_id]] = "FASTA_first_token->GFF_CDS:protein_id"
        if protein_id in transcript_to_gene:
            candidates[transcript_to_gene[protein_id]] = "FASTA_first_token->GFF_transcript"
        if len(candidates) > 1:
            raise ValueError(
                f"Provider protein has discordant exact primary mappings: {protein_id}"
            )
        if not candidates:
            excluded += 1
            continue
        gene_id, relation_source = next(iter(candidates.items()))
        if gene_id not in coding_genes:
            excluded += 1
            continue
        sequence = raw_sequence.upper().rstrip("*")
        if not sequence or set(sequence) - PROTEIN_ALPHABET:
            raise ValueError(f"Invalid primary provider protein sequence: {protein_id}")
        relation_sources[relation_source] += 1
        current = representatives.get(gene_id)
        candidate = (protein_id, relation_source, sequence)
        if current is None or (-len(sequence), protein_id) < (-len(current[2]), current[0]):
            representatives[gene_id] = candidate
    selected = {gene: representatives[gene] for gene in sorted(coding_genes) if gene in representatives}
    missing = sorted(coding_genes - set(selected))
    if missing:
        raise ValueError(
            "Primary coding genes lack a unique header-aware provider protein: "
            + ",".join(missing[:10])
        )
    selected_ids = [item[0] for item in selected.values()]
    if len(selected_ids) != len(set(selected_ids)):
        raise ValueError("One provider protein maps to multiple primary coding genes")

    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x", encoding="utf-8", newline="") as handle:
        for gene_id in sorted(selected):
            item = selected[gene_id]
            handle.write(f">{item[0]} gene:{gene_id}\n")
            for offset in range(0, len(item[2]), 60):
                handle.write(item[2][offset:offset + 60] + "\n")
    _write_tsv(
        whitelist,
        ("gene_id", "protein_id", "relation_source"),
        (
            {
                "gene_id": gene_id,
                "protein_id": selected[gene_id][0],
                "relation_source": selected[gene_id][1],
            }
            for gene_id in sorted(selected)
        ),
    )
    observed = [protein_id for protein_id, _, _ in iter_fasta(output)]
    if observed != selected_ids or len(observed) != len(set(observed)):
        raise ValueError("Walnut primary protein output has missing, duplicate, or extra records")
    manifest = {
        "schema_version": "ploidypatch.walnut_primary_candidate_proteins.v0.8",
        "holdout_id": HOLDOUT_ID,
        "policy_id": POLICY_ID,
        "mapping": "primary_GFF_CDS_gene_plus_header_aware_exact_transcript_relation",
        "hierarchy": hierarchy,
        "counts": {
            "provider_records": provider_records,
            "primary_coding_genes": len(coding_genes),
            "primary_CDS_protein_ids": len(cds_protein_ids),
            "accepted_unique_representatives": len(selected),
            "excluded_provider_records": provider_records - len(selected),
            "unmapped_or_nonprimary_provider_records": excluded,
            "mapped_alternative_isoforms_not_selected": (
                provider_records - excluded - len(selected)
            ),
            "relation_source_counts": dict(sorted(relation_sources.items())),
        },
        "inputs": {
            "primary_gff_sha256": sha256_file(primary_gff_path),
            "provider_protein_sha256": sha256_file(provider_protein_path),
        },
        "outputs": {
            "protein_fasta_sha256": sha256_file(output),
            "whitelist_tsv_sha256": sha256_file(whitelist),
        },
        "truth_access": False,
        "ranker_access": False,
    }
    manifest_output.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return manifest


def load_json_object(path: str | Path) -> dict[str, Any]:
    source = Path(path)
    if not source.is_file() or source.is_symlink() or source.stat().st_size == 0:
        raise ValueError(f"Missing, empty or symlinked JSON: {source}")

    def reject_duplicates(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        value: dict[str, Any] = {}
        for key, item in pairs:
            if key in value:
                raise ValueError(f"Duplicate JSON key in {source}: {key}")
            value[key] = item
        return value

    result = json.loads(
        source.read_text(encoding="utf-8"), object_pairs_hook=reject_duplicates
    )
    if not isinstance(result, dict):
        raise ValueError(f"Expected JSON object: {source}")
    return result


def read_tsv(path: str | Path, required: Iterable[str]) -> list[dict[str, str]]:
    source = Path(path)
    if not source.is_file() or source.is_symlink() or source.stat().st_size == 0:
        raise ValueError(f"Missing, empty or symlinked TSV: {source}")
    with source.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        missing = set(required) - set(fields)
        if missing:
            raise ValueError(f"TSV lacks fields {sorted(missing)}: {source}")
        return list(reader)


def canonical_pair(left: str, right: str) -> tuple[str, str]:
    if not left or not right or left == right:
        raise ValueError(f"Invalid unordered gene pair: {left!r}, {right!r}")
    return tuple(sorted((left, right)))


def _write_tsv(
    path: Path, fields: tuple[str, ...], rows: Iterable[Mapping[str, Any]]
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def _pair_map(
    rows: Iterable[dict[str, str]], left_field: str, right_field: str
) -> dict[tuple[str, str], dict[str, str]]:
    result: dict[tuple[str, str], dict[str, str]] = {}
    for row in rows:
        pair = canonical_pair(row[left_field], row[right_field])
        if pair in result:
            raise ValueError(f"Duplicate unordered pair: {pair}")
        result[pair] = row
    return result


def _normalized_pair_map(
    rows: Iterable[dict[str, str]], left_field: str, right_field: str
) -> dict[tuple[str, str], dict[str, str]]:
    result: dict[tuple[str, str], dict[str, str]] = {}
    for row in rows:
        pair = canonical_pair(
            normalize_feature_id(row[left_field]),
            normalize_feature_id(row[right_field]),
        )
        if pair in result:
            raise ValueError(f"Duplicate normalized unordered pair: {pair}")
        result[pair] = row
    return result


def filter_self_pairs_by_closed_yn00_ks(
    *,
    self_pairs_path: str | Path,
    ks_path: str | Path,
    output_pairs_path: str | Path,
    decisions_path: str | Path,
    minimum: float = KS_MINIMUM,
    maximum: float = KS_MAXIMUM,
) -> dict[str, Any]:
    """Apply the frozen closed YN00 interval to reciprocal target self-pairs."""

    if (minimum, maximum) != (KS_MINIMUM, KS_MAXIMUM):
        raise ValueError("Walnut v0.8 requires the exact closed YN00 Ks 0.10-0.75")
    pair_output = Path(output_pairs_path)
    decision_output = Path(decisions_path)
    manifest_output = Path(str(pair_output) + ".manifest.json")
    if any(path.exists() or path.is_symlink() for path in (
        pair_output, decision_output, manifest_output
    )):
        raise FileExistsError("Refusing to overwrite Walnut Ks pair outputs")

    self_rows = read_tsv(self_pairs_path, ("gene_id_a", "gene_id_b"))
    self_by_pair = _normalized_pair_map(self_rows, "gene_id_a", "gene_id_b")
    ks_rows = read_tsv(ks_path, ("id1", "id2", "ks_YN00"))
    ks_by_pair = _normalized_pair_map(ks_rows, "id1", "id2")

    fields = ("gene_id_a", "gene_id_b", "yn00_ks", "status", "reason")
    decisions: list[dict[str, str]] = []
    accepted: list[dict[str, str]] = []
    for pair, source in sorted(self_by_pair.items()):
        ks_row = ks_by_pair.get(pair)
        raw = "" if ks_row is None else ks_row.get("ks_YN00", "").strip()
        try:
            value = float(raw)
        except ValueError:
            value = float("nan")
        if ks_row is None or value != value:
            status, reason = "rejected", "missing_or_nonfinite_yn00_ks"
        elif not minimum <= value <= maximum:
            status, reason = "rejected", "yn00_ks_outside_closed_0.10_0.75"
        else:
            status, reason = "accepted", "closed_yn00_ks_0.10_0.75"
        source_pair = canonical_pair(source["gene_id_a"], source["gene_id_b"])
        row = {
            "gene_id_a": source_pair[0],
            "gene_id_b": source_pair[1],
            "yn00_ks": raw,
            "status": status,
            "reason": reason,
        }
        decisions.append(row)
        if status == "accepted":
            accepted.append(row)

    _write_tsv(pair_output, fields[:3], accepted)
    _write_tsv(decision_output, fields, decisions)
    reasons = Counter(row["reason"] for row in decisions)
    manifest = {
        "schema_version": "ploidypatch.walnut_closed_yn00_pair_filter.v0.8",
        "holdout_id": HOLDOUT_ID,
        "policy_id": POLICY_ID,
        "target_event": TARGET_EVENT,
        "parameters": {
            "interval": "closed",
            "minimum": minimum,
            "maximum": maximum,
            "missing_or_nonfinite": "abstain",
            "thresholds_fitted_from_target": False,
            "pair_identifier_join": "exact_normalize_feature_id_then_unordered_pair",
        },
        "counts": {
            "self_pairs": len(self_by_pair),
            "accepted_pairs": len(accepted),
            "reasons": dict(sorted(reasons.items())),
        },
        "inputs": {
            "self_pairs_sha256": sha256_file(self_pairs_path),
            "yn00_ks_sha256": sha256_file(ks_path),
        },
        "outputs": {
            "pairs_sha256": sha256_file(pair_output),
            "decisions_sha256": sha256_file(decision_output),
        },
        "truth_labels_accessed": False,
    }
    manifest_output.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return manifest


def infer_exact_two_outgroup_pair_consistent_truth(
    *,
    ks_filtered_self_pairs_path: str | Path,
    evaluator_group_decisions: Mapping[str, str | Path],
    output_dir: str | Path,
) -> dict[str, Any]:
    """Intersect target Ks pairs with two exact pair-consistent outgroups.

    Every row in a group decision table is a qualifying exact-1:2 proposal
    after the frozen block and multiplicity filters.  An alternate proposal
    incident to either target member vetoes that pair even when the alternate
    proposal was already rejected by a reciprocal-uniqueness check.
    """

    if set(evaluator_group_decisions) != {"corylus", "castanea"}:
        raise ValueError("Walnut truth requires exactly Corylus and Castanea")
    output = Path(output_dir)
    if output.exists() or output.is_symlink():
        raise FileExistsError(f"Refusing to overwrite Walnut pair truth: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    working = Path(tempfile.mkdtemp(prefix=f".{output.name}.working.", dir=output.parent))
    try:
        self_rows = read_tsv(
            ks_filtered_self_pairs_path, ("gene_id_a", "gene_id_b", "yn00_ks")
        )
        self_pairs = _pair_map(self_rows, "gene_id_a", "gene_id_b")
        group_pairs: dict[str, set[tuple[str, str]]] = {}
        group_audits: list[dict[str, Any]] = []
        for label in ("corylus", "castanea"):
            path = Path(evaluator_group_decisions[label])
            rows = read_tsv(path, ("gene_id_a", "gene_id_b", "status", "reason"))
            # Several distinct evaluator counterparts may independently yield
            # the same target pair.  Collapse those counterparts to one target
            # proposal before the cross-pair consistency audit; duplicates are
            # supporting evidence, not a malformed target-pair table.
            proposals: set[tuple[str, str]] = set()
            accepted_proposals: set[tuple[str, str]] = set()
            for row in rows:
                pair = canonical_pair(row["gene_id_a"], row["gene_id_b"])
                proposals.add(pair)
                if row["status"] == "accepted":
                    accepted_proposals.add(pair)
            partners: dict[str, set[str]] = defaultdict(set)
            for left, right in proposals:
                partners[left].add(right)
                partners[right].add(left)
            discordant = {gene for gene, values in partners.items() if len(values) > 1}
            accepted = {
                pair
                for pair in accepted_proposals
                if pair[0] not in discordant
                and pair[1] not in discordant
            }
            group_pairs[label] = accepted
            group_audits.append(
                {
                    "group": label,
                    "qualifying_exact_1_to_2_pairs": len(proposals),
                    "pair_consistent_pairs": len(accepted),
                    "discordant_target_members": len(discordant),
                    "input_sha256": sha256_file(path),
                }
            )

        final_pairs = set(self_pairs) & group_pairs["corylus"] & group_pairs["castanea"]
        decision_fields = (
            "gene_id_a", "gene_id_b", "yn00_ks", "status", "reason"
        )
        decision_rows: list[dict[str, str]] = []
        for pair, self_row in sorted(self_pairs.items()):
            missing = [label for label in ("corylus", "castanea") if pair not in group_pairs[label]]
            status = "accepted" if not missing else "rejected"
            reason = (
                "exact_self_corylus_castanea_intersection"
                if not missing
                else "missing_or_discordant_group:" + ",".join(missing)
            )
            decision_rows.append(
                {
                    "gene_id_a": pair[0],
                    "gene_id_b": pair[1],
                    "yn00_ks": self_row["yn00_ks"],
                    "status": status,
                    "reason": reason,
                }
            )
        accepted_rows = [row for row in decision_rows if row["status"] == "accepted"]
        _write_tsv(working / "pairs.tsv", decision_fields[:3], accepted_rows)
        _write_tsv(working / "decisions.tsv", decision_fields, decision_rows)
        manifest = {
            "schema_version": PAIR_SCHEMA,
            "holdout_id": HOLDOUT_ID,
            "policy_id": POLICY_ID,
            "target_event": TARGET_EVENT,
            "formal_rule": (
                "exact_unordered_pair_intersection_of_closed_ks_target_self_wgdi_"
                "and_pair_consistent_corylus_and_castanea"
            ),
            "parameters": {
                "minimum_block_pairs": MIN_BLOCK_PAIRS,
                "cross_primary_chromosome": True,
                "reciprocal_unique": True,
                "counterpart_target_multiplicity": "exactly_two",
                "ks_interval": "closed_0.10_to_0.75",
                "best_pair_rescue": False,
            },
            "counts": {
                "ks_filtered_self_pairs": len(self_pairs),
                "accepted_truth_pairs": len(final_pairs),
                "evaluator_groups": group_audits,
            },
            "inputs": {
                "ks_filtered_self_pairs_sha256": sha256_file(
                    ks_filtered_self_pairs_path
                )
            },
            "outputs": {
                "pairs_sha256": sha256_file(working / "pairs.tsv"),
                "decisions_sha256": sha256_file(working / "decisions.tsv"),
            },
            "candidate_references_used_for_truth": False,
            "truth_labels_accessed": False,
        }
        (working / "manifest.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        write_sha256sums(working)
        verify_sha256sums(working, ignore_checksum_file=True)
        os.replace(working, output)
        return manifest
    except BaseException:
        shutil.rmtree(working, ignore_errors=True)
        raise


def seal_blind_pool_manifest(
    *,
    manifest_path: str | Path,
    candidate_gff_path: str | Path,
    decisions_path: str | Path,
    raw_predictions_manifest_path: str | Path,
    policy_arm: str,
) -> dict[str, Any]:
    if policy_arm not in POOL_SCHEMAS:
        raise ValueError(f"Unsupported Walnut H1 arm: {policy_arm}")
    path = Path(manifest_path)
    manifest = load_json_object(path)
    if manifest.get("schema_version") != POOL_SCHEMAS[policy_arm]:
        raise ValueError("Walnut H1 pool manifest schema differs")
    raw = load_json_object(raw_predictions_manifest_path)
    if raw.get("schema_version") != RAW_SCHEMA:
        raise ValueError("Walnut H1 raw-prediction manifest schema differs")
    if raw.get("truth_access") is not False or raw.get("ranker_access") is not False:
        raise ValueError("Raw projections violate the Walnut blind firewall")
    inputs = manifest.get("inputs")
    outputs = manifest.get("outputs")
    if not isinstance(inputs, dict) or not isinstance(outputs, dict):
        raise ValueError("Walnut H1 pool manifest lacks inputs/outputs objects")
    candidate_output = outputs.get("candidate_gff")
    decision_output = outputs.get("decisions")
    if (
        not isinstance(candidate_output, dict)
        or candidate_output.get("sha256") != sha256_file(candidate_gff_path)
        or not isinstance(decision_output, dict)
        or decision_output.get("sha256") != sha256_file(decisions_path)
    ):
        raise ValueError("Walnut H1 pool output hashes differ before sealing")
    inputs["raw_predictions_manifest"] = {
        "file_name": "raw_predictions.manifest.json",
        "sha256": sha256_file(raw_predictions_manifest_path),
    }
    manifest.update(
        {
            "holdout_id": HOLDOUT_ID,
            "policy_id": POLICY_ID,
            "policy_arm": policy_arm,
            "truth_access": False,
            "ranker_access": False,
            "automatic_approval": False,
            "within_method_reference_vote_count": 1,
            "candidate_references_per_method": 2,
        }
    )
    temporary = path.with_name(path.name + ".working")
    if temporary.exists() or temporary.is_symlink():
        raise FileExistsError(f"Refusing stale pool-manifest working file: {temporary}")
    temporary.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    os.replace(temporary, path)
    return manifest


def canonical_tree_record(root: str | Path) -> dict[str, Any]:
    tree = Path(root)
    if not tree.is_dir() or tree.is_symlink():
        raise ValueError(f"Raw prediction tree is missing or symlinked: {tree}")
    files = sorted(
        path for path in tree.rglob("*") if path.is_file() and not path.is_symlink()
    )
    if not files:
        raise ValueError(f"Raw prediction tree is empty: {tree}")
    digest = hashlib.sha256()
    total = 0
    for path in files:
        relative = path.relative_to(tree).as_posix()
        file_digest = sha256_file(path)
        digest.update(f"{file_digest}  {relative}\n".encode("utf-8"))
        total += path.stat().st_size
    return {
        "file_count": len(files),
        "bytes": total,
        "tree_sha256": digest.hexdigest(),
    }


def build_raw_predictions_manifest(
    *,
    project_root: str | Path,
    raw_prediction_trees: Mapping[str, str | Path],
    input_hashes: Mapping[str, str],
    output: str | Path,
) -> dict[str, Any]:
    required = {
        f"{method}__{bundle}"
        for method in ("miniprot", "gemoma", "lifton")
        for bundle in ("candidate_mandshurica", "candidate_carya")
    }
    if set(raw_prediction_trees) != required:
        raise ValueError("Walnut raw predictions require three methods by two references")
    destination = Path(output)
    if destination.exists() or destination.is_symlink():
        raise FileExistsError(f"Refusing to overwrite raw prediction manifest: {destination}")
    root = Path(project_root).resolve()
    records: dict[str, Any] = {}
    for label, raw_path in sorted(raw_prediction_trees.items()):
        tree = Path(raw_path).resolve()
        try:
            relative = tree.relative_to(root).as_posix()
        except ValueError as error:
            raise ValueError(f"Raw prediction tree escapes project root: {tree}") from error
        records[label] = {"relative_path": relative, **canonical_tree_record(tree)}
    if not input_hashes or any(
        len(value) != 64 or any(character not in "0123456789abcdef" for character in value)
        for value in input_hashes.values()
    ):
        raise ValueError("Walnut raw prediction input hashes are incomplete or malformed")
    payload = {
        "schema_version": RAW_SCHEMA,
        "holdout_id": HOLDOUT_ID,
        "policy_id": POLICY_ID,
        "truth_access": False,
        "ranker_access": False,
        "method_families": ["miniprot", "gemoma", "lifton"],
        "candidate_references": ["candidate_mandshurica", "candidate_carya"],
        "within_method_reference_vote_count": 1,
        "tree_hash_algorithm": "sha256_of_sorted_sha256_two_space_relative_path_newline",
        "input_hashes": dict(sorted(input_hashes.items())),
        "raw_prediction_trees": records,
    }
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return payload


def evaluate_h1_scores(
    *, retain_score_path: str | Path, suppress_score_path: str | Path, output: str | Path
) -> dict[str, Any]:
    destination = Path(output)
    bootstrap_path = destination.with_name("paired_event_bootstrap.json")
    if destination.exists() or destination.is_symlink():
        raise FileExistsError(f"Refusing to overwrite Walnut H1 evaluation: {destination}")
    if bootstrap_path.exists() or bootstrap_path.is_symlink():
        raise FileExistsError(
            f"Refusing to overwrite Walnut H1 bootstrap: {bootstrap_path}"
        )
    destination.parent.mkdir(parents=True, exist_ok=True)
    retain = load_json_object(retain_score_path)
    suppress = load_json_object(suppress_score_path)
    for label, score in (("retain_distinct", retain), ("suppress_overlap", suppress)):
        if not score_collateral_gate(score):
            raise ValueError(f"Walnut {label} collateral-loss or quality sentinel failed")
    bootstrap = paired_event_bootstrap(
        score_inputs=(
            ("retain_distinct", retain_score_path),
            ("suppress_overlap", suppress_score_path),
        ),
        output_json_path=bootstrap_path,
        metric="complete_cds_chain_recovery",
        replicates=BOOTSTRAP_REPLICATES,
        seed=BOOTSTRAP_SEED,
    )
    difference = bootstrap["paired_differences"][0]
    if (
        difference["left"] != "retain_distinct"
        or difference["right"] != "suppress_overlap"
    ):
        raise ValueError("Walnut H1 bootstrap contrast direction differs")
    def arm_summary(score: dict[str, Any]) -> dict[str, Any]:
        recovery = score.get("event_recovery", {})
        collateral = score.get("collateral_changes", {})
        return {
            "events": recovery.get("events"),
            "exact_phased_cds_chain_recovered": recovery.get(
                "complete_cds_chain_recovery"
            ),
            "exact_phased_cds_chain_recall": recovery.get(
                "complete_cds_chain_recall"
            ),
            "baseline_transcript_structure_loss": collateral.get(
                "baseline_transcript_structures_missing_from_candidate"
            ),
        }

    result = {
        "schema_version": EVALUATION_SCHEMA,
        "holdout_id": HOLDOUT_ID,
        "policy_id": POLICY_ID,
        "status": "ready",
        "reason_codes": [],
        "hypothesis": "H1_retain_distinct_vs_suppress_overlap_only",
        "metric": "event_exact_phased_CDS_recall_retain_distinct_minus_suppress_overlap",
        "arms": {
            "retain_distinct": arm_summary(retain),
            "suppress_overlap": arm_summary(suppress),
        },
        "paired_event_bootstrap": difference,
        "bootstrap_parameters": bootstrap["parameters"],
        "bootstrap_replicates": BOOTSTRAP_REPLICATES,
        "bootstrap_unit": "paired_event",
        "success": (
            difference["observed_delta"] > 0 and difference["ci_lower"] > 0
        ),
        "sentinels": {
            "all_arm_collateral_loss_zero": True,
            "automatic_approval": False,
        },
        "protocol_profile": "core_H1_only_no_ranker",
        "ranker_or_model_executed": False,
        "h2_or_topology_ranking_executed": False,
        "all_arm_collateral_loss": 0,
        "inputs": {
            "retain_score_sha256": sha256_file(retain_score_path),
            "suppress_score_sha256": sha256_file(suppress_score_path),
        },
    }
    temporary = destination.with_name(f".{destination.name}.working.{os.getpid()}")
    if temporary.exists() or temporary.is_symlink():
        raise FileExistsError(f"Refusing stale evaluation working file: {temporary}")
    temporary.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    os.replace(temporary, destination)
    return result
