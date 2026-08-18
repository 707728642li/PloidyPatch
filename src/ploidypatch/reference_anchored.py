from __future__ import annotations

import csv
import hashlib
import json
import os
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

from . import __version__


REFERENCE_ANCHORED_SCHEMA = "ploidypatch.reference_anchored_projection.v1"
CANDIDATE_SOURCE_FIELDS = (
    "candidate_digest",
    "candidate_gene_id",
    "arm_id",
    "reference_id",
    "source_gene_id",
)
WGD_PAIR_FIELDS = (
    "source_gene_id",
    "partner_source_gene_id",
    "support_block_count",
    "longest_block_pairs",
)
HOMELOG_BASE_FIELDS = (
    "arm_id",
    "source_homeolog_gene_id",
    "evidence_status",
    "target_base_gene_id",
    "evidence_reason",
)
ARM_DECISION_FIELDS = (
    "candidate_digest",
    "candidate_gene_id",
    "arm_id",
    "reference_id",
    "source_gene_id",
    "source_homeolog_gene_id",
    "arm_status",
    "target_base_gene_id",
    "arm_reason",
    "support_block_count",
    "longest_block_pairs",
)
SELECTION_FIELDS = (
    "gene_id",
    "consensus_digest",
    "pair_id",
    "partner_gene_id",
    "support_block_count",
    "longest_block_pairs",
    "status",
    "reason",
    "evidence_type",
    "source_base_gene_id",
    "anchor_cell_count",
    "compatible_partner_count",
    "qualifying_base_hit_count",
    "backbone_id",
)
SAFE_LABEL = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")


def _file_sha256(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        while block := handle.read(8 * 1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def _read_exact_tsv(path: Path, expected: tuple[str, ...]) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != expected:
            raise ValueError(
                f"TSV columns differ for {path}: expected {expected}, "
                f"observed {tuple(reader.fieldnames or ())}"
            )
        return list(reader)


def _parse_named_paths(values: list[str], option: str) -> dict[str, Path]:
    output: dict[str, Path] = {}
    for value in values:
        reference, separator, path_text = value.partition("=")
        if not separator or not SAFE_LABEL.fullmatch(reference) or not path_text:
            raise ValueError(f"invalid {option} value: {value!r}")
        if reference in output:
            raise ValueError(f"duplicate reference in {option}: {reference}")
        output[reference] = Path(path_text)
    if not output:
        raise ValueError(f"at least one {option} value is required")
    return output


def _positive_integer(value: str, field: str) -> int:
    try:
        parsed = int(value)
    except ValueError as exc:
        raise ValueError(f"invalid integer in {field}: {value!r}") from exc
    if parsed < 1:
        raise ValueError(f"{field} must be positive")
    return parsed


def aggregate_reference_anchored_projection(
    *,
    candidate_source_provenance_path: str | Path,
    accepted_wgd_pair_inputs: list[str],
    homeolog_base_evidence_inputs: list[str],
    output_dir_path: str | Path,
    evidence_type: str = "reference_anchored_projection",
) -> dict[str, Any]:
    """Aggregate candidate-independent reference anchors without species rules.

    WGD event selection and source-homeolog projection are upstream adapter
    responsibilities. This function only joins their frozen outputs and applies
    the generic fail-closed consensus rule: missing arms do not veto, at least
    one evidence arm is required, and any conflict or distinct target partners
    causes abstention.
    """

    if not SAFE_LABEL.fullmatch(evidence_type):
        raise ValueError("evidence_type must be a safe label")
    candidate_path = Path(candidate_source_provenance_path)
    pair_paths = _parse_named_paths(accepted_wgd_pair_inputs, "--accepted-wgd-pairs")
    evidence_paths = _parse_named_paths(
        homeolog_base_evidence_inputs, "--homeolog-base-evidence"
    )
    if set(pair_paths) != set(evidence_paths):
        raise ValueError("WGD-pair and homeolog-base reference sets differ")
    for path in (candidate_path, *pair_paths.values(), *evidence_paths.values()):
        if not path.is_file() or path.is_symlink():
            raise ValueError(f"input must be a regular file: {path}")

    output = Path(output_dir_path)
    working = Path(f"{output}.working")
    collisions = [path for path in (output, working) if path.exists()]
    if collisions:
        raise FileExistsError(
            "refusing to overwrite reference-anchored output: "
            + ", ".join(str(path) for path in collisions)
        )

    candidate_rows = _read_exact_tsv(candidate_path, CANDIDATE_SOURCE_FIELDS)
    candidates: dict[str, dict[str, Any]] = {}
    seen_candidate_arms: set[tuple[str, str, str, str]] = set()
    global_source_arm_owner: dict[tuple[str, str, str], str] = {}
    for line_number, row in enumerate(candidate_rows, start=2):
        if any(not row[field] for field in CANDIDATE_SOURCE_FIELDS):
            raise ValueError(f"empty candidate provenance field on line {line_number}")
        reference = row["reference_id"]
        if reference not in pair_paths:
            raise ValueError(f"unknown provenance reference: {reference}")
        arm_key = (
            row["candidate_digest"],
            row["arm_id"],
            reference,
            row["source_gene_id"],
        )
        if arm_key in seen_candidate_arms:
            raise ValueError(f"duplicate candidate source arm: {arm_key}")
        seen_candidate_arms.add(arm_key)
        global_arm_key = (reference, row["arm_id"], row["source_gene_id"])
        previous_owner = global_source_arm_owner.get(global_arm_key)
        if previous_owner is not None and previous_owner != row["candidate_digest"]:
            raise ValueError(
                "reference/arm/source provenance is reused by different "
                f"candidates: {global_arm_key} -> {previous_owner}, "
                f"{row['candidate_digest']}"
            )
        global_source_arm_owner[global_arm_key] = row["candidate_digest"]
        candidate = candidates.setdefault(
            row["candidate_digest"],
            {"candidate_gene_id": row["candidate_gene_id"], "arms": []},
        )
        if candidate["candidate_gene_id"] != row["candidate_gene_id"]:
            raise ValueError(
                f"candidate gene differs within digest: {row['candidate_digest']}"
            )
        candidate["arms"].append(row)
    if not candidates:
        raise ValueError("candidate provenance is empty")

    pair_maps: dict[str, dict[str, tuple[str, int, int]]] = {}
    for reference, path in pair_paths.items():
        rows = _read_exact_tsv(path, WGD_PAIR_FIELDS)
        mapping: dict[str, tuple[str, int, int]] = {}
        for line_number, row in enumerate(rows, start=2):
            source = row["source_gene_id"]
            partner = row["partner_source_gene_id"]
            if not source or not partner or source == partner or source in mapping:
                raise ValueError(
                    f"invalid/duplicate accepted WGD pair on {path}:{line_number}"
                )
            mapping[source] = (
                partner,
                _positive_integer(row["support_block_count"], "support_block_count"),
                _positive_integer(row["longest_block_pairs"], "longest_block_pairs"),
            )
        for source, (partner, block_count, longest) in mapping.items():
            reverse = mapping.get(partner)
            if reverse is None or reverse != (source, block_count, longest):
                raise ValueError(
                    f"accepted WGD pairs are not exact reciprocal records: "
                    f"{reference}:{source}:{partner}"
                )
        pair_maps[reference] = mapping

    evidence_maps: dict[
        str, dict[tuple[str, str], tuple[str, str | None, str]]
    ] = {}
    for reference, path in evidence_paths.items():
        rows = _read_exact_tsv(path, HOMELOG_BASE_FIELDS)
        mapping: dict[tuple[str, str], tuple[str, str | None, str]] = {}
        for line_number, row in enumerate(rows, start=2):
            key = (row["arm_id"], row["source_homeolog_gene_id"])
            status = row["evidence_status"]
            target = row["target_base_gene_id"] or None
            reason = row["evidence_reason"]
            if (
                not key[0]
                or not key[1]
                or key in mapping
                or status not in {"evidence", "missing", "conflict"}
                or not reason
                or (status == "evidence") != (target is not None)
            ):
                raise ValueError(
                    f"invalid/duplicate homeolog-base evidence on {path}:{line_number}"
                )
            mapping[key] = (status, target, reason)
        evidence_maps[reference] = mapping

    policy = {
        "minimum_evidence_arms": 1,
        "missing_evidence_arms_veto": False,
        "conflicting_evidence_arm_veto": True,
        "all_evidence_arms_must_agree_on_target_base_gene": True,
        "wgd_event_selection_is_upstream_adapter_responsibility": True,
        "homeolog_base_projection_is_upstream_adapter_responsibility": True,
        "support_summary": "minimum_across_unique_evidence_events",
    }
    backbone_payload = {
        "schema_version": REFERENCE_ANCHORED_SCHEMA,
        "policy": policy,
        "references": {
            reference: {
                "accepted_wgd_pairs_sha256": _file_sha256(pair_paths[reference]),
                "homeolog_base_evidence_sha256": _file_sha256(
                    evidence_paths[reference]
                ),
            }
            for reference in sorted(pair_paths)
        },
    }
    backbone_digest = hashlib.sha256(
        json.dumps(
            backbone_payload, sort_keys=True, separators=(",", ":")
        ).encode()
    ).hexdigest()
    backbone_id = f"PPRA-BB-{backbone_digest[:24]}"

    working.mkdir(parents=True)
    try:
        arm_decisions: list[dict[str, str]] = []
        selections: list[dict[str, str]] = []
        reasons: Counter[str] = Counter()
        accepted_count = 0
        for digest in sorted(candidates):
            candidate = candidates[digest]
            evidence_arms: list[tuple[str, int, int, tuple[str, str, str]]] = []
            conflict = False
            for row in sorted(
                candidate["arms"],
                key=lambda item: (
                    item["reference_id"], item["arm_id"], item["source_gene_id"]
                ),
            ):
                reference = row["reference_id"]
                pair = pair_maps[reference].get(row["source_gene_id"])
                if pair is None:
                    homeolog = ""
                    status, target, reason = (
                        "missing",
                        None,
                        "source_gene_has_no_accepted_wgd_partner",
                    )
                    block_count = longest = ""
                else:
                    homeolog, block_value, longest_value = pair
                    block_count, longest = str(block_value), str(longest_value)
                    status, target, reason = evidence_maps[reference].get(
                        (row["arm_id"], homeolog),
                        (
                            "missing",
                            None,
                            "source_homeolog_has_no_base_evidence_record",
                        ),
                    )
                    if status == "evidence":
                        assert target is not None
                        evidence_arms.append(
                            (
                                target,
                                block_value,
                                longest_value,
                                (reference, row["source_gene_id"], homeolog),
                            )
                        )
                    elif status == "conflict":
                        conflict = True
                arm_decisions.append(
                    {
                        "candidate_digest": digest,
                        "candidate_gene_id": candidate["candidate_gene_id"],
                        "arm_id": row["arm_id"],
                        "reference_id": reference,
                        "source_gene_id": row["source_gene_id"],
                        "source_homeolog_gene_id": homeolog,
                        "arm_status": status,
                        "target_base_gene_id": target or "",
                        "arm_reason": reason,
                        "support_block_count": block_count,
                        "longest_block_pairs": longest,
                    }
                )

            targets = {item[0] for item in evidence_arms}
            if conflict:
                final_reason = "abstain_conflicting_evidence_arm"
            elif not evidence_arms:
                final_reason = "abstain_no_evidence_arm"
            elif len(targets) != 1:
                final_reason = "abstain_discordant_target_base_partners"
            else:
                final_reason = "accepted_reference_anchored_unique_consensus"
            reasons[final_reason] += 1
            accepted = final_reason.startswith("accepted_")
            partner = next(iter(targets)) if accepted else ""
            if accepted:
                accepted_count += 1
                unique_events = {
                    item[3]: (item[1], item[2]) for item in evidence_arms
                }
                block_count = str(min(item[0] for item in unique_events.values()))
                longest = str(min(item[1] for item in unique_events.values()))
                pair_id = "PPRA-" + hashlib.sha256(
                    f"{backbone_id}\t{digest}\t{partner}".encode()
                ).hexdigest()[:20]
            else:
                block_count = longest = pair_id = ""
            selections.append(
                {
                    "gene_id": candidate["candidate_gene_id"],
                    "consensus_digest": digest,
                    "pair_id": pair_id,
                    "partner_gene_id": partner,
                    "support_block_count": block_count,
                    "longest_block_pairs": longest,
                    "status": "accepted" if accepted else "rejected",
                    "reason": final_reason,
                    "evidence_type": evidence_type,
                    "source_base_gene_id": "",
                    "anchor_cell_count": str(len(evidence_arms)),
                    "compatible_partner_count": str(len(targets)),
                    "qualifying_base_hit_count": str(len(evidence_arms)),
                    "backbone_id": backbone_id,
                }
            )

        def write_tsv(
            path: Path, fields: tuple[str, ...], rows: list[dict[str, str]]
        ) -> None:
            with path.open("x", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
                )
                writer.writeheader()
                writer.writerows(rows)

        write_tsv(working / "arm_decisions.tsv", ARM_DECISION_FIELDS, arm_decisions)
        write_tsv(working / "selection.tsv", SELECTION_FIELDS, selections)
        manifest: dict[str, Any] = {
            "schema_version": REFERENCE_ANCHORED_SCHEMA,
            "generator": {"name": "PloidyPatch", "version": __version__},
            "truth_access": False,
            "labels_used": False,
            "candidate_or_target_union_graph_used": False,
            "backbone_id": backbone_id,
            "backbone_definition_sha256": backbone_digest,
            "policy": policy,
            "counts": {
                "candidates": len(candidates),
                "candidate_source_arms": len(candidate_rows),
                "accepted": accepted_count,
                "rejected": len(candidates) - accepted_count,
                "reasons": dict(sorted(reasons.items())),
            },
            "inputs": {
                "candidate_source_provenance": {
                    "path": candidate_path.name,
                    "bytes": candidate_path.stat().st_size,
                    "sha256": _file_sha256(candidate_path),
                },
                "accepted_wgd_pairs": {
                    reference: {
                        "path": path.name,
                        "bytes": path.stat().st_size,
                        "sha256": _file_sha256(path),
                    }
                    for reference, path in sorted(pair_paths.items())
                },
                "homeolog_base_evidence": {
                    reference: {
                        "path": path.name,
                        "bytes": path.stat().st_size,
                        "sha256": _file_sha256(path),
                    }
                    for reference, path in sorted(evidence_paths.items())
                },
            },
            "outputs": {},
        }
        for name in ("arm_decisions.tsv", "selection.tsv"):
            path = working / name
            manifest["outputs"][name] = {
                "bytes": path.stat().st_size,
                "sha256": _file_sha256(path),
            }
        (working / "manifest.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        with (working / "SHA256SUMS").open(
            "x", encoding="utf-8", newline=""
        ) as handle:
            for name in ("arm_decisions.tsv", "manifest.json", "selection.tsv"):
                handle.write(f"{_file_sha256(working / name)}  {name}\n")
        os.replace(working, output)
        return manifest
    except BaseException:
        # The non-overwriting .working tree is intentionally retained for audit.
        raise
