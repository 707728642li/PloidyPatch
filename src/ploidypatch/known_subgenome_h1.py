"""Fail-closed truth primitives for known allopolyploid subgenome events.

These routines are deliberately evaluator-only.  They select target self-WGDI
pairs using a predeclared chromosome-level homoeolog map, then require exact
pair-consistent support from two event-outside evaluator groups.  No candidate
reference, score, label, fitted Ks peak, or best-pair rescue is accepted.
"""
from __future__ import annotations

from collections import Counter, defaultdict
import csv
import json
import os
from pathlib import Path
import shutil
import tempfile
from typing import Any, Iterable, Mapping

from .artifact_manifest import sha256_file, verify_sha256sums, write_sha256sums


HOLDOUT_ID = "coffea_et39_v1.0"
POLICY_ID = "ploidypatch_coffea_external_core_h1_v1.0"
TARGET_EVENT = "recent_Arabica_C_E_subgenome_homoeolog_pairing"
SELF_FILTER_SCHEMA = "ploidypatch.known_subgenome_self_pair_filter.v1.0"
PAIR_SCHEMA = "ploidypatch.known_subgenome_exact_pair_truth.v1.0"
EVALUATOR_GROUPS = ("gardenia", "ophiorrhiza")
MIN_BLOCK_PAIRS = 20


def attach_descriptive_yn00_ks(
    *,
    self_pairs_path: str | Path,
    ks_path: str | Path,
    output_path: str | Path,
) -> dict[str, Any]:
    """Attach exact WGDI-pair YN00 Ks without filtering or imputing values."""

    destination = Path(output_path)
    manifest_path = Path(str(destination) + ".manifest.json")
    if any(path.exists() or path.is_symlink() for path in (destination, manifest_path)):
        raise FileExistsError("Refusing to overwrite descriptive YN00 artifacts")
    rows = _read_tsv(
        self_pairs_path,
        ("wgdi_gene_id_a", "wgdi_gene_id_b", "gene_id_a", "gene_id_b"),
    )
    with Path(ks_path).open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not {"id1", "id2", "ks_YN00"} <= set(reader.fieldnames or ()):
            raise ValueError("YN00 table lacks id1, id2 and ks_YN00")
        ks_rows = list(reader)
    values: dict[tuple[str, str], str] = {}
    for row in ks_rows:
        pair = _canonical_pair(row["id1"], row["id2"])
        if pair in values:
            raise ValueError(f"Duplicate YN00 pair: {pair}")
        values[pair] = row["ks_YN00"]
    if not rows:
        raise ValueError("Self-pair table has no accepted pairs")
    fields = tuple(rows[0])
    if "yn00_ks" in fields:
        raise ValueError("Self-pair table already contains yn00_ks")
    output_rows: list[dict[str, str]] = []
    matched = 0
    for row in rows:
        pair = _canonical_pair(row["wgdi_gene_id_a"], row["wgdi_gene_id_b"])
        value = values.get(pair, "")
        matched += bool(value)
        output_rows.append({**row, "yn00_ks": value})
    destination.parent.mkdir(parents=True, exist_ok=True)
    _write_tsv(destination, (*fields, "yn00_ks"), output_rows)
    manifest = {
        "schema_version": "ploidypatch.descriptive_yn00_attachment.v1.0",
        "holdout_id": HOLDOUT_ID,
        "policy_id": POLICY_ID,
        "selection_use": False,
        "missing_value_policy": "report_empty_never_filter_or_impute",
        "counts": {
            "self_pairs": len(rows),
            "yn00_rows": len(ks_rows),
            "matched_self_pairs": matched,
            "missing_self_pairs": len(rows) - matched,
        },
        "inputs": {
            "self_pairs_sha256": sha256_file(self_pairs_path),
            "yn00_sha256": sha256_file(ks_path),
        },
        "output_sha256": sha256_file(destination),
    }
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return manifest


def _read_tsv(path: str | Path, required: Iterable[str]) -> list[dict[str, str]]:
    source = Path(path)
    if not source.is_file() or source.is_symlink() or source.stat().st_size == 0:
        raise ValueError(f"Missing, empty or symlinked TSV: {source}")
    with source.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        if len(fields) != len(set(fields)):
            raise ValueError(f"Duplicate TSV field name: {source}")
        missing = set(required) - set(fields)
        if missing:
            raise ValueError(f"TSV lacks fields {sorted(missing)}: {source}")
        rows = list(reader)
    return rows


def _write_tsv(
    path: Path, fields: tuple[str, ...], rows: Iterable[Mapping[str, Any]]
) -> None:
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


def _canonical_pair(left: str, right: str) -> tuple[str, str]:
    if not left or not right or left == right:
        raise ValueError(f"Invalid exact unordered target pair: {left!r}, {right!r}")
    return tuple(sorted((left, right)))


def _atomic_working_directory(output: Path) -> Path:
    if output.exists() or output.is_symlink():
        raise FileExistsError(f"Refusing to overwrite known-subgenome output: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    return Path(tempfile.mkdtemp(prefix=f".{output.name}.working.", dir=output.parent))


def _load_homoeolog_map(
    path: str | Path,
) -> tuple[dict[str, tuple[str, str]], dict[str, tuple[str, str]]]:
    rows = _read_tsv(
        path,
        (
            "homoeolog_group",
            "c_subgenome_seqid",
            "e_subgenome_seqid",
            "c_label",
            "e_label",
        ),
    )
    if not rows:
        raise ValueError("Known-subgenome table has no homoeolog groups")
    by_seqid: dict[str, tuple[str, str]] = {}
    by_group: dict[str, tuple[str, str]] = {}
    for row in rows:
        group = row["homoeolog_group"]
        c_seqid = row["c_subgenome_seqid"]
        e_seqid = row["e_subgenome_seqid"]
        if not group or not c_seqid or not e_seqid or c_seqid == e_seqid:
            raise ValueError("Malformed known-subgenome homoeolog group")
        if group in by_group or c_seqid in by_seqid or e_seqid in by_seqid:
            raise ValueError("Duplicate homoeolog group or chromosome assignment")
        if row["c_label"] != f"{group}c" or row["e_label"] != f"{group}e":
            raise ValueError("Homoeolog labels differ from frozen group+subgenome labels")
        by_group[group] = (c_seqid, e_seqid)
        by_seqid[c_seqid] = (group, "c")
        by_seqid[e_seqid] = (group, "e")
    return by_seqid, by_group


def filter_self_pairs_by_known_subgenome(
    *,
    self_pairs_path: str | Path,
    homoeolog_groups_path: str | Path,
    output_dir: str | Path,
) -> dict[str, Any]:
    """Keep only exact same-group C/E pairs from accepted target self-WGDI.

    The upstream self-pair table is expected to have already applied the
    frozen cross-chromosome, minimum-20-block and reciprocal-unique rules.
    YN00 Ks may be present and is copied for audit, but is never consulted.
    """

    output = Path(output_dir)
    working = _atomic_working_directory(output)
    try:
        by_seqid, by_group = _load_homoeolog_map(homoeolog_groups_path)
        source_rows = _read_tsv(
            self_pairs_path,
            (
                "pair_id",
                "gene_id_a",
                "gene_id_b",
                "seqid_a",
                "seqid_b",
                "support_block_count",
                "longest_block_pairs",
            ),
        )
        fields = (
            "pair_id",
            "gene_id_a",
            "gene_id_b",
            "seqid_a",
            "seqid_b",
            "homoeolog_group",
            "subgenome_a",
            "subgenome_b",
            "support_block_count",
            "longest_block_pairs",
            "yn00_ks",
        )
        decision_fields = (*fields, "status", "reason")
        accepted: list[dict[str, str]] = []
        decisions: list[dict[str, str]] = []
        seen_pairs: set[tuple[str, str]] = set()
        seen_pair_ids: set[str] = set()
        for source in source_rows:
            pair_id = source["pair_id"]
            if not pair_id or pair_id in seen_pair_ids:
                raise ValueError(f"Missing or duplicate target self pair_id: {pair_id!r}")
            seen_pair_ids.add(pair_id)
            members = sorted(
                (
                    (source["gene_id_a"], source["seqid_a"]),
                    (source["gene_id_b"], source["seqid_b"]),
                )
            )
            pair = _canonical_pair(members[0][0], members[1][0])
            if pair in seen_pairs:
                raise ValueError(f"Duplicate exact unordered target self pair: {pair}")
            seen_pairs.add(pair)
            annotations = [by_seqid.get(members[0][1]), by_seqid.get(members[1][1])]
            if any(annotation is None for annotation in annotations):
                group, subgenomes = "", ("", "")
                status, reason = "rejected", "nonprimary_or_unregistered_seqid"
            else:
                first = annotations[0]
                second = annotations[1]
                assert first is not None and second is not None
                group = first[0] if first[0] == second[0] else ""
                subgenomes = (first[1], second[1])
                if first[0] != second[0]:
                    status, reason = "rejected", "different_predeclared_homoeolog_group"
                elif set(subgenomes) != {"c", "e"}:
                    status, reason = "rejected", "same_subgenome_or_duplicate_chromosome"
                elif members[0][1] == members[1][1]:
                    status, reason = "rejected", "same_primary_chromosome"
                else:
                    status, reason = "accepted", "predeclared_same_group_exact_C_E_pair"
            row = {
                "pair_id": pair_id,
                "gene_id_a": members[0][0],
                "gene_id_b": members[1][0],
                "seqid_a": members[0][1],
                "seqid_b": members[1][1],
                "homoeolog_group": group,
                "subgenome_a": subgenomes[0],
                "subgenome_b": subgenomes[1],
                "support_block_count": source["support_block_count"],
                "longest_block_pairs": source["longest_block_pairs"],
                "yn00_ks": source.get("yn00_ks", ""),
                "status": status,
                "reason": reason,
            }
            decisions.append(row)
            if status == "accepted":
                accepted.append(row)

        accepted.sort(key=lambda row: (row["gene_id_a"], row["gene_id_b"]))
        decisions.sort(key=lambda row: (row["gene_id_a"], row["gene_id_b"]))
        _write_tsv(working / "pairs.tsv", fields, accepted)
        _write_tsv(working / "decisions.tsv", decision_fields, decisions)
        reasons = Counter(row["reason"] for row in decisions)
        manifest = {
            "schema_version": SELF_FILTER_SCHEMA,
            "holdout_id": HOLDOUT_ID,
            "policy_id": POLICY_ID,
            "target_event": TARGET_EVENT,
            "formal_rule": "predeclared_same_group_exactly_one_C_and_one_E",
            "parameters": {
                "minimum_block_pairs": MIN_BLOCK_PAIRS,
                "cross_primary_chromosome": True,
                "reciprocal_unique_upstream_required": True,
                "yn00_ks_policy": "descriptive_only_not_used_for_selection",
                "candidate_references_used": False,
                "best_pair_rescue": False,
            },
            "counts": {
                "input_self_pairs": len(source_rows),
                "accepted_known_subgenome_pairs": len(accepted),
                "homoeolog_groups": len(by_group),
                "reasons": dict(sorted(reasons.items())),
            },
            "inputs": {
                "self_pairs_sha256": sha256_file(self_pairs_path),
                "homoeolog_groups_sha256": sha256_file(homoeolog_groups_path),
            },
            "outputs": {
                "pairs_sha256": sha256_file(working / "pairs.tsv"),
                "decisions_sha256": sha256_file(working / "decisions.tsv"),
            },
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


def infer_exact_two_outgroup_pair_consistent_truth(
    *,
    known_subgenome_self_pairs_path: str | Path,
    evaluator_group_decisions: Mapping[str, str | Path],
    output_dir: str | Path,
) -> dict[str, Any]:
    """Require exact target-pair support from Gardenia and Ophiorrhiza.

    All qualifying proposals, including rejected alternate pairs, participate
    in the pair-consistency veto.  Multiple counterparts supporting the same
    unordered target pair are collapsed as concordant evidence.
    """

    if set(evaluator_group_decisions) != set(EVALUATOR_GROUPS):
        raise ValueError("Known-subgenome truth requires Gardenia and Ophiorrhiza")
    output = Path(output_dir)
    working = _atomic_working_directory(output)
    try:
        self_rows = _read_tsv(
            known_subgenome_self_pairs_path,
            ("pair_id", "gene_id_a", "gene_id_b", "homoeolog_group"),
        )
        self_pairs: dict[tuple[str, str], dict[str, str]] = {}
        for row in self_rows:
            pair = _canonical_pair(row["gene_id_a"], row["gene_id_b"])
            if pair in self_pairs:
                raise ValueError(f"Duplicate known-subgenome self pair: {pair}")
            self_pairs[pair] = row

        group_pairs: dict[str, set[tuple[str, str]]] = {}
        group_audits: list[dict[str, Any]] = []
        for label in EVALUATOR_GROUPS:
            path = Path(evaluator_group_decisions[label])
            rows = _read_tsv(path, ("gene_id_a", "gene_id_b", "status", "reason"))
            proposals: set[tuple[str, str]] = set()
            accepted_proposals: set[tuple[str, str]] = set()
            for row in rows:
                if row["status"] not in {"accepted", "rejected"}:
                    raise ValueError(f"Unexpected evaluator proposal status: {row['status']}")
                pair = _canonical_pair(row["gene_id_a"], row["gene_id_b"])
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
                if pair[0] not in discordant and pair[1] not in discordant
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

        decision_fields = (
            "pair_id",
            "gene_id_a",
            "gene_id_b",
            "homoeolog_group",
            "status",
            "reason",
        )
        decisions: list[dict[str, str]] = []
        for pair, source in sorted(self_pairs.items()):
            missing = [label for label in EVALUATOR_GROUPS if pair not in group_pairs[label]]
            status = "accepted" if not missing else "rejected"
            reason = (
                "exact_self_gardenia_ophiorrhiza_intersection"
                if not missing
                else "missing_or_discordant_group:" + ",".join(missing)
            )
            decisions.append(
                {
                    "pair_id": source["pair_id"],
                    "gene_id_a": pair[0],
                    "gene_id_b": pair[1],
                    "homoeolog_group": source["homoeolog_group"],
                    "status": status,
                    "reason": reason,
                }
            )
        accepted = [row for row in decisions if row["status"] == "accepted"]
        _write_tsv(working / "pairs.tsv", decision_fields[:4], accepted)
        _write_tsv(working / "decisions.tsv", decision_fields, decisions)
        manifest = {
            "schema_version": PAIR_SCHEMA,
            "holdout_id": HOLDOUT_ID,
            "policy_id": POLICY_ID,
            "target_event": TARGET_EVENT,
            "formal_rule": (
                "exact_unordered_pair_intersection_of_predeclared_subgenome_"
                "target_self_wgdi_and_pair_consistent_gardenia_and_ophiorrhiza"
            ),
            "parameters": {
                "minimum_block_pairs": MIN_BLOCK_PAIRS,
                "cross_primary_chromosome": True,
                "reciprocal_unique": True,
                "counterpart_target_multiplicity": "exactly_two",
                "yn00_ks_policy": "descriptive_only_not_used_for_selection",
                "best_pair_rescue": False,
            },
            "counts": {
                "known_subgenome_self_pairs": len(self_pairs),
                "accepted_truth_pairs": len(accepted),
                "evaluator_groups": group_audits,
            },
            "inputs": {
                "known_subgenome_self_pairs_sha256": sha256_file(
                    known_subgenome_self_pairs_path
                ),
                "evaluator_group_decisions": {
                    label: sha256_file(evaluator_group_decisions[label])
                    for label in EVALUATOR_GROUPS
                },
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
