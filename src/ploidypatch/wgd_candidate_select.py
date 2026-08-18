from __future__ import annotations

import csv
import io
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

from . import __version__
from .baseline import _file_sha256
from .consensus import _consume_candidate_prefix
from .gff import parse_attributes
from .io import open_text


WGD_CANDIDATE_SELECTION_SCHEMA = "ploidypatch.wgd_candidate_selection.v1"
WGD_CONFLICT_PROPAGATION_SCHEMA = "ploidypatch.wgd_conflict_propagation.v1"


def _base_gene_ids(path: Path) -> set[str]:
    gene_ids: set[str] = set()
    with open_text(path) as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            stripped = raw_line.rstrip("\r\n")
            if not stripped or stripped.startswith("#"):
                continue
            fields = stripped.split("\t")
            if len(fields) != 9:
                raise ValueError(f"Malformed base GFF line {line_number}")
            if fields[2] != "gene":
                continue
            attributes, malformed = parse_attributes(fields[8])
            gene_id = attributes.get("ID", "")
            if malformed or not gene_id or gene_id in gene_ids:
                raise ValueError(f"Invalid or duplicate base gene ID at line {line_number}")
            gene_ids.add(gene_id)
    return gene_ids


def _candidate_blocks(
    candidate_path: Path, base_path: Path
) -> tuple[dict[str, list[str]], dict[str, str]]:
    blocks: dict[str, list[str]] = {}
    digests: dict[str, str] = {}
    current_gene = ""
    current_transcripts: set[str] = set()
    current_has_cds = False

    def finish_current() -> None:
        nonlocal current_gene, current_transcripts, current_has_cds
        if not current_gene:
            return
        if not current_transcripts or not current_has_cds:
            raise ValueError(f"Incomplete consensus hierarchy: {current_gene}")
        current_gene = ""
        current_transcripts = set()
        current_has_cds = False

    with candidate_path.open("rb") as raw_handle:
        _consume_candidate_prefix(base_path, raw_handle, candidate_path)
        text_handle = io.TextIOWrapper(
            raw_handle, encoding="utf-8", errors="strict", newline=None
        )
        try:
            for suffix_line, raw_line in enumerate(text_handle, start=1):
                stripped = raw_line.rstrip("\r\n")
                if not stripped or stripped.startswith("#"):
                    continue
                fields = stripped.split("\t")
                if len(fields) != 9:
                    raise ValueError(
                        f"Malformed consensus GFF line {suffix_line}: {candidate_path}"
                    )
                if fields[1] != "PloidyPatchConsensus":
                    raise ValueError(
                        f"Unexpected appended source at line {suffix_line}: {candidate_path}"
                    )
                attributes, malformed = parse_attributes(fields[8])
                if malformed:
                    raise ValueError(
                        f"Malformed consensus attributes at line {suffix_line}"
                    )
                feature_type = fields[2]
                if feature_type == "gene":
                    finish_current()
                    current_gene = attributes.get("ID", "")
                    digest = attributes.get("consensus_digest", "")
                    if (
                        not current_gene
                        or current_gene in blocks
                        or not digest
                        or len(digest) != 64
                    ):
                        raise ValueError(
                            f"Invalid consensus gene at line {suffix_line}"
                        )
                    blocks[current_gene] = [raw_line]
                    digests[current_gene] = digest
                elif feature_type in {"mRNA", "transcript"}:
                    if not current_gene:
                        raise ValueError("Consensus transcript precedes its gene")
                    parents = {
                        value
                        for value in attributes.get("Parent", "").split(",")
                        if value
                    }
                    transcript_id = attributes.get("ID", "")
                    if parents != {current_gene} or not transcript_id:
                        raise ValueError(f"Invalid consensus transcript hierarchy")
                    current_transcripts.add(transcript_id)
                    blocks[current_gene].append(raw_line)
                elif feature_type in {"exon", "CDS"}:
                    if not current_gene:
                        raise ValueError("Consensus child precedes its gene")
                    parents = {
                        value
                        for value in attributes.get("Parent", "").split(",")
                        if value
                    }
                    if len(parents) != 1 or not parents <= current_transcripts:
                        raise ValueError(f"Invalid consensus child hierarchy")
                    if feature_type == "CDS":
                        current_has_cds = True
                    blocks[current_gene].append(raw_line)
                else:
                    raise ValueError(
                        f"Unsupported consensus feature type: {feature_type}"
                    )
            finish_current()
        finally:
            text_handle.detach()
    if not blocks:
        raise ValueError("Consensus candidate GFF contains no appended models")
    return blocks, digests


def select_wgd_supported_candidates(
    *,
    base_gff_path: str | Path,
    candidate_gff_path: str | Path,
    pair_tsv_path: str | Path,
    output_gff_path: str | Path,
    selection_tsv_path: str | Path,
) -> dict[str, Any]:
    """Retain candidates paired uniquely to an existing annotated WGD partner."""

    base_path = Path(base_gff_path)
    candidate_path = Path(candidate_gff_path)
    pair_path = Path(pair_tsv_path)
    output_path = Path(output_gff_path)
    selection_path = Path(selection_tsv_path)
    manifest_path = Path(str(output_path) + ".manifest.json")
    for required in (base_path, candidate_path, pair_path):
        if not required.is_file() or required.stat().st_size == 0:
            raise ValueError(f"Missing or empty WGD candidate input: {required}")
        if required.suffix == ".gz":
            raise ValueError("WGD candidate selection requires uncompressed inputs")
    collisions = [
        path for path in (output_path, selection_path, manifest_path) if path.exists()
    ]
    if collisions:
        raise FileExistsError(
            "Refusing to overwrite WGD candidate artifact(s): "
            + ", ".join(str(path) for path in collisions)
        )

    base_genes = _base_gene_ids(base_path)
    candidate_blocks, candidate_digests = _candidate_blocks(candidate_path, base_path)
    candidate_genes = set(candidate_blocks)
    if base_genes & candidate_genes:
        raise ValueError("Consensus candidate gene IDs collide with the base annotation")
    known_genes = base_genes | candidate_genes

    candidate_base_pairs: dict[str, list[dict[str, str]]] = defaultdict(list)
    candidate_candidate_pairs: dict[str, list[dict[str, str]]] = defaultdict(list)
    pair_rows = 0
    pair_classes: Counter[str] = Counter()
    with pair_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required_fields = {
            "pair_id",
            "gene_id_a",
            "gene_id_b",
            "support_block_count",
            "longest_block_pairs",
        }
        if reader.fieldnames is None or not required_fields <= set(reader.fieldnames):
            raise ValueError("WGD pair table lacks required columns")
        seen_pair_ids: set[str] = set()
        for line_number, row in enumerate(reader, start=2):
            pair_rows += 1
            pair_id = row["pair_id"]
            gene_a = row["gene_id_a"]
            gene_b = row["gene_id_b"]
            if not pair_id or pair_id in seen_pair_ids or gene_a == gene_b:
                raise ValueError(f"Invalid or duplicate WGD pair at line {line_number}")
            seen_pair_ids.add(pair_id)
            unknown = {gene_a, gene_b} - known_genes
            if unknown:
                raise ValueError(
                    f"WGD pair references genes absent from candidate GFF at line {line_number}"
                )
            endpoints = {gene_a, gene_b}
            candidate_endpoints = endpoints & candidate_genes
            if len(candidate_endpoints) == 1:
                candidate_gene = next(iter(candidate_endpoints))
                partner = next(iter(endpoints - candidate_endpoints))
                if partner not in base_genes:
                    raise ValueError("Candidate WGD partner is not a base gene")
                candidate_base_pairs[candidate_gene].append({**row, "partner": partner})
                pair_classes["candidate_base"] += 1
            elif len(candidate_endpoints) == 2:
                for candidate_gene in candidate_endpoints:
                    partner = next(iter(endpoints - {candidate_gene}))
                    candidate_candidate_pairs[candidate_gene].append(
                        {**row, "partner": partner}
                    )
                pair_classes["candidate_candidate"] += 1
            else:
                pair_classes["base_base"] += 1

    selection_rows: list[dict[str, str]] = []
    selected_genes: set[str] = set()
    reason_counts: Counter[str] = Counter()
    for gene_id in sorted(candidate_genes):
        base_support = candidate_base_pairs.get(gene_id, [])
        candidate_support = candidate_candidate_pairs.get(gene_id, [])
        if len(base_support) == 1 and not candidate_support:
            support = base_support[0]
            status = "accepted"
            reason = "reciprocal_wgd_partner_is_existing_gene"
            selected_genes.add(gene_id)
        elif len(base_support) > 1:
            raise ValueError(f"Candidate has multiple accepted base WGD partners: {gene_id}")
        elif candidate_support:
            support = candidate_support[0]
            status = "rejected"
            reason = "candidate_partner_only"
        else:
            support = {}
            status = "rejected"
            reason = "no_accepted_wgd_partner"
        reason_counts[reason] += 1
        selection_rows.append(
            {
                "gene_id": gene_id,
                "consensus_digest": candidate_digests[gene_id],
                "pair_id": support.get("pair_id", ""),
                "partner_gene_id": support.get("partner", ""),
                "support_block_count": support.get("support_block_count", ""),
                "longest_block_pairs": support.get("longest_block_pairs", ""),
                "status": status,
                "reason": reason,
            }
        )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with base_path.open("rb") as source_handle, output_path.open("xb") as output_handle:
        last_byte = b""
        while chunk := source_handle.read(8 * 1024 * 1024):
            output_handle.write(chunk)
            last_byte = chunk[-1:]
        if last_byte not in {b"\n", b"\r"}:
            output_handle.write(b"\n")
        output_handle.write(b"###\n")
        for gene_id in sorted(selected_genes, key=lambda value: candidate_digests[value]):
            output_handle.write("".join(candidate_blocks[gene_id]).encode("utf-8"))

    selection_path.parent.mkdir(parents=True, exist_ok=True)
    selection_fields = (
        "gene_id",
        "consensus_digest",
        "pair_id",
        "partner_gene_id",
        "support_block_count",
        "longest_block_pairs",
        "status",
        "reason",
    )
    with selection_path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=selection_fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(selection_rows)

    manifest: dict[str, Any] = {
        "schema_version": WGD_CANDIDATE_SELECTION_SCHEMA,
        "generator": {"name": "PloidyPatch", "version": __version__},
        "inputs": {
            "base_gff": {
                "file_name": base_path.name,
                "sha256": _file_sha256(base_path),
                "genes": len(base_genes),
            },
            "candidate_gff": {
                "file_name": candidate_path.name,
                "sha256": _file_sha256(candidate_path),
                "models": len(candidate_genes),
            },
            "wgd_pairs": {
                "file_name": pair_path.name,
                "sha256": _file_sha256(pair_path),
                "rows": pair_rows,
            },
        },
        "policy": {
            "accepted_relation": "one_candidate_to_one_existing_gene",
            "candidate_candidate_support": "rejected_as_circular",
            "pair_table_requirement": "upstream_reciprocal_unique",
            "claim_boundary": (
                "blind candidate WGD context, not gene-tree-confirmed homeology"
            ),
        },
        "counts": {
            "candidate_models": len(candidate_genes),
            "selected_models": len(selected_genes),
            "rejected_models": len(candidate_genes) - len(selected_genes),
            "pair_classes": dict(sorted(pair_classes.items())),
            "decision_counts": dict(sorted(reason_counts.items())),
        },
        "outputs": {
            "candidate_gff": {
                "file_name": output_path.name,
                "sha256": _file_sha256(output_path),
            },
            "selection": {
                "file_name": selection_path.name,
                "sha256": _file_sha256(selection_path),
                "rows": len(selection_rows),
            },
        },
    }
    with manifest_path.open("x", encoding="utf-8", newline="") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return manifest


def propagate_wgd_selection_to_conflict_pool(
    *,
    base_gff_path: str | Path,
    candidate_gff_path: str | Path,
    pool_decisions_tsv_path: str | Path,
    prior_wgd_selection_tsv_path: str | Path,
    output_selection_tsv_path: str | Path,
) -> dict[str, Any]:
    """Propagate a unique locus-level WGD partner to retained chain alternatives."""

    base_path = Path(base_gff_path)
    candidate_path = Path(candidate_gff_path)
    decisions_path = Path(pool_decisions_tsv_path)
    prior_path = Path(prior_wgd_selection_tsv_path)
    output_path = Path(output_selection_tsv_path)
    manifest_path = Path(str(output_path) + ".manifest.json")
    for required in (base_path, candidate_path, decisions_path, prior_path):
        if not required.is_file() or required.stat().st_size == 0:
            raise ValueError(f"Missing or empty WGD propagation input: {required}")
    collisions = [path for path in (output_path, manifest_path) if path.exists()]
    if collisions:
        raise FileExistsError(
            "Refusing to overwrite WGD propagation artifact(s): "
            + ", ".join(str(path) for path in collisions)
        )

    _, digest_by_gene = _candidate_blocks(candidate_path, base_path)
    gene_by_digest = {digest: gene for gene, digest in digest_by_gene.items()}
    if len(gene_by_digest) != len(digest_by_gene):
        raise ValueError("Candidate pool contains duplicate digests")
    with decisions_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required_fields = {
            "consensus_digest", "status", "conflict_set_digest",
            "conflict_member_count",
        }
        if reader.fieldnames is None or not required_fields <= set(reader.fieldnames):
            raise ValueError("Candidate-pool decisions lack conflict fields")
        pool_rows = list(reader)
    accepted_pool = [row for row in pool_rows if row["status"] == "accepted"]
    if {row["consensus_digest"] for row in accepted_pool} != set(gene_by_digest):
        raise ValueError("Accepted pool decisions disagree with candidate GFF")
    pool_by_digest = {row["consensus_digest"]: row for row in accepted_pool}
    if len(pool_by_digest) != len(accepted_pool):
        raise ValueError("Candidate-pool decisions contain duplicate digests")

    with prior_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required_fields = {
            "gene_id", "consensus_digest", "pair_id", "partner_gene_id",
            "support_block_count", "longest_block_pairs", "status", "reason",
        }
        if reader.fieldnames is None or not required_fields <= set(reader.fieldnames):
            raise ValueError("Prior WGD selection lacks required fields")
        prior_rows = list(reader)
    prior_by_digest = {row["consensus_digest"]: row for row in prior_rows}
    if len(prior_by_digest) != len(prior_rows):
        raise ValueError("Prior WGD selection contains duplicate digests")
    if set(prior_by_digest) - set(pool_by_digest):
        raise ValueError("Prior WGD selection is not a subset of the candidate pool")

    accepted_prior_by_conflict: dict[str, list[dict[str, str]]] = defaultdict(list)
    for digest, prior in prior_by_digest.items():
        conflict_set = pool_by_digest[digest]["conflict_set_digest"]
        if prior["status"] == "accepted" and conflict_set:
            accepted_prior_by_conflict[conflict_set].append(prior)

    output_rows: list[dict[str, str]] = []
    reasons: Counter[str] = Counter()
    inherited = 0
    for digest in sorted(pool_by_digest):
        pool = pool_by_digest[digest]
        exact = prior_by_digest.get(digest)
        source_digest = digest if exact is not None else ""
        evidence_origin = "prior_exact_candidate"
        if exact is not None:
            status = exact["status"]
            reason = exact["reason"]
            pair_id = exact["pair_id"]
            partner = exact["partner_gene_id"]
            support_blocks = exact["support_block_count"]
            longest_block = exact["longest_block_pairs"]
        else:
            conflict_set = pool["conflict_set_digest"]
            supports = accepted_prior_by_conflict.get(conflict_set, []) if conflict_set else []
            partners = {row["partner_gene_id"] for row in supports if row["partner_gene_id"]}
            if len(partners) == 1:
                support = max(
                    supports,
                    key=lambda row: (
                        int(row["support_block_count"] or 0),
                        int(row["longest_block_pairs"] or 0),
                        row["consensus_digest"],
                    ),
                )
                status = "accepted"
                reason = "conflict_set_unique_existing_wgd_partner"
                pair_id = support["pair_id"]
                partner = support["partner_gene_id"]
                support_blocks = support["support_block_count"]
                longest_block = support["longest_block_pairs"]
                source_digest = support["consensus_digest"]
                evidence_origin = "inherited_from_conflict_set"
                inherited += 1
            elif len(partners) > 1:
                status = "rejected"
                reason = "conflict_set_multiple_existing_wgd_partners"
                pair_id = partner = support_blocks = longest_block = ""
                evidence_origin = "ambiguous_conflict_set"
            else:
                status = "rejected"
                reason = "no_inherited_existing_wgd_partner"
                pair_id = partner = support_blocks = longest_block = ""
                evidence_origin = "no_prior_conflict_support"
        reasons[reason] += 1
        output_rows.append(
            {
                "gene_id": gene_by_digest[digest],
                "consensus_digest": digest,
                "pair_id": pair_id,
                "partner_gene_id": partner,
                "support_block_count": support_blocks,
                "longest_block_pairs": longest_block,
                "status": status,
                "reason": reason,
                "evidence_origin": evidence_origin,
                "source_candidate_digest": source_digest,
                "conflict_set_digest": pool["conflict_set_digest"],
                "conflict_member_count": pool["conflict_member_count"],
            }
        )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fields = tuple(output_rows[0]) if output_rows else ()
    with output_path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        if output_rows:
            writer.writeheader()
            writer.writerows(output_rows)
    manifest: dict[str, Any] = {
        "schema_version": WGD_CONFLICT_PROPAGATION_SCHEMA,
        "generator": {"name": "PloidyPatch", "version": __version__},
        "inputs": {
            "base_gff_sha256": _file_sha256(base_path),
            "candidate_gff_sha256": _file_sha256(candidate_path),
            "pool_decisions_sha256": _file_sha256(decisions_path),
            "prior_wgd_selection_sha256": _file_sha256(prior_path),
        },
        "policy": {
            "exact_candidate": "preserve_prior_decision",
            "new_alternative": "inherit_only_one_unique_existing_partner_within_conflict_set",
            "multiple_partners": "reject_as_ambiguous",
            "candidate_candidate_support": "not_propagated",
            "automatic_approval": False,
        },
        "counts": {
            "pool_candidates": len(output_rows),
            "prior_candidates": len(prior_rows),
            "accepted": sum(row["status"] == "accepted" for row in output_rows),
            "inherited_accepted": inherited,
            "decision_reasons": dict(sorted(reasons.items())),
        },
        "output": {
            "file_name": output_path.name,
            "sha256": _file_sha256(output_path),
            "rows": len(output_rows),
        },
    }
    with manifest_path.open("x", encoding="utf-8", newline="") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return manifest
