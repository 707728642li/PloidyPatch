from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import pytest

from ploidypatch.fixed_backbone_projection import (
    CATALOG_FIELDS,
    HIT_FIELDS,
    FixedBackboneProjectionPolicy,
    project_candidates_to_fixed_backbone,
)
from ploidypatch.fixed_target_backbone import (
    FixedTargetBackbonePolicy,
    build_fixed_target_backbone,
)
from ploidypatch.cli import build_parser


GENES = (
    ("chr1", "A1", 100, 199, 1),
    ("chr1", "A2", 400, 499, 2),
    ("chr1", "A3", 700, 799, 3),
    ("chr2", "B1", 1100, 1199, 1),
    ("chr2", "BX", 1250, 1349, 2),
    ("chr2", "BY", 1370, 1469, 3),
    ("chr2", "B2", 1500, 1599, 4),
    ("chr2", "B3", 1800, 1899, 5),
)


def write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8", newline="")
    return path


def make_backbone(tmp_path: Path) -> Path:
    gff_lines = ["##gff-version 3\n"]
    query_lines: list[str] = []
    for seqid, gene, start, end, order in GENES:
        transcript = f"{gene}.t1"
        gff_lines.extend(
            (
                f"{seqid}\ttest\tgene\t{start}\t{end}\t.\t+\t.\tID={gene}\n",
                f"{seqid}\ttest\tmRNA\t{start}\t{end}\t.\t+\t.\tID={transcript};Parent={gene}\n",
                f"{seqid}\ttest\tCDS\t{start}\t{end}\t.\t+\t0\tID={gene}.cds;Parent={transcript}\n",
            )
        )
        query_lines.append(f"{seqid}\t{gene}\t{start}\t{end}\t+\t{order}\n")
    pairs = (("A1", "B1"), ("A2", "B2"), ("A3", "B3"))
    collinearity = (
        "# Alignment 1: score=100 pvalue=0.01 N=3 chr1&chr2 plus\n"
        + "".join(f"{left} 1 {right} 1 1\n" for left, right in pairs)
    )
    output = tmp_path / "backbone"
    build_fixed_target_backbone(
        base_gff_path=write(tmp_path / "base.gff3", "".join(gff_lines)),
        query_wgdi_gff_path=write(tmp_path / "query.gff", "".join(query_lines)),
        base_only_collinearity_path=write(tmp_path / "self.collinearity", collinearity),
        policy=FixedTargetBackbonePolicy(
            primary_seqids=("chr1", "chr2"), min_block_pairs=2
        ),
        output_dir_path=output,
    )
    return output


def write_tsv(path: Path, fields: tuple[str, ...], rows: list[dict[str, object]]) -> Path:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return path


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_hits_manifest(
    path: Path, backbone: Path, catalog: Path, hits: Path, rows: int
) -> Path:
    backbone_id = json.loads((backbone / "manifest.json").read_text(encoding="utf-8"))[
        "backbone_id"
    ]
    value = {
        "schema_version": "ploidypatch.fixed_base_protein_hits.v0.6",
        "truth_access": False,
        "candidate_candidate_hits": False,
        "backbone_id": backbone_id,
        "database_scope": "fixed_backbone_representative_proteins_only",
        "exhaustive_target_search": True,
        "candidate_catalog_sha256": sha256(catalog),
        "base_gene_count": len(GENES),
        "query_protein_fasta_sha256": "a" * 64,
        "base_protein_fasta_sha256": "b" * 64,
        "parameters": {
            "maximum_target_sequences": 0,
            "maximum_evalue": 1e-5,
            "more_sensitive": True,
            "outfmt_fields": list(HIT_FIELDS),
        },
        "output": {"hits_sha256": sha256(hits), "rows": rows},
    }
    path.write_text(json.dumps(value) + "\n", encoding="utf-8", newline="")
    return path


def candidate(digest_character: str, gene: str, protein: str, start: int, end: int) -> dict[str, object]:
    return {
        "candidate_digest": digest_character * 64,
        "candidate_gene_id": gene,
        "query_protein_id": protein,
        "seqid": "chr1",
        "start": start,
        "end": end,
        "strand": "+",
    }


def hit(protein: str, subject: str) -> dict[str, object]:
    return {
        "query_protein_id": protein,
        "subject_gene_id": subject,
        "pident": 80,
        "alignment_length": 90,
        "query_length": 100,
        "evalue": "1e-20",
        "bitscore": 200,
    }


def run_projection(
    tmp_path: Path,
    backbone: Path,
    candidates: list[dict[str, object]],
    hits: list[dict[str, object]],
    name: str,
) -> Path:
    catalog = write_tsv(tmp_path / f"{name}.catalog.tsv", CATALOG_FIELDS, candidates)
    hit_path = write_tsv(tmp_path / f"{name}.hits.tsv", HIT_FIELDS, hits)
    hit_manifest = write_hits_manifest(
        tmp_path / f"{name}.hits.manifest.json",
        backbone,
        catalog,
        hit_path,
        len(hits),
    )
    output = tmp_path / name
    project_candidates_to_fixed_backbone(
        backbone_dir_path=backbone,
        candidate_catalog_tsv_path=catalog,
        fixed_base_hits_tsv_path=hit_path,
        fixed_base_hits_manifest_path=hit_manifest,
        policy=FixedBackboneProjectionPolicy(),
        output_dir_path=output,
    )
    return output


def read_selection(path: Path) -> list[dict[str, str]]:
    with (path / "selection.tsv").open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_gap_projection_requires_one_hit_inside_strict_fixed_bracket(tmp_path: Path) -> None:
    backbone = make_backbone(tmp_path)
    output = run_projection(
        tmp_path,
        backbone,
        [candidate("a", "C1", "P1", 250, 349)],
        [hit("P1", "BX")],
        "projection",
    )
    row = read_selection(output)[0]
    assert row["status"] == "accepted"
    assert row["partner_gene_id"] == "BX"
    assert row["reason"] == "accepted_strict_gap_unique_fixed_partner"
    assert row["evidence_type"] == "strict_gap_anchor_projection"


def test_candidate_result_is_invariant_to_unrelated_candidate_universe(tmp_path: Path) -> None:
    backbone = make_backbone(tmp_path)
    first = run_projection(
        tmp_path,
        backbone,
        [candidate("a", "C1", "P1", 250, 349)],
        [hit("P1", "BX")],
        "one",
    )
    second = run_projection(
        tmp_path,
        backbone,
        [
            candidate("a", "C1", "P1", 250, 349),
            candidate("b", "C2", "P2", 550, 649),
        ],
        [hit("P1", "BX"), hit("P2", "BY")],
        "two",
    )
    first_row = read_selection(first)[0]
    second_row = next(
        row for row in read_selection(second) if row["consensus_digest"] == "a" * 64
    )
    assert first_row == second_row


def test_multiple_compatible_gap_partners_abstain_without_winner_selection(tmp_path: Path) -> None:
    backbone = make_backbone(tmp_path)
    output = run_projection(
        tmp_path,
        backbone,
        [candidate("a", "C1", "P1", 250, 349)],
        [hit("P1", "BX"), hit("P1", "BY")],
        "projection",
    )
    row = read_selection(output)[0]
    assert row["status"] == "rejected"
    assert row["partner_gene_id"] == ""
    assert row["reason"] == "multiple_compatible_fixed_partners"
    assert row["compatible_partner_count"] == "2"


def test_direct_overlap_uses_only_unique_fixed_edge_and_never_gap_fallback(tmp_path: Path) -> None:
    backbone = make_backbone(tmp_path)
    accepted = run_projection(
        tmp_path,
        backbone,
        [candidate("a", "C1", "P1", 120, 180)],
        [hit("P1", "B1")],
        "accepted",
    )
    row = read_selection(accepted)[0]
    assert row["status"] == "accepted"
    assert row["partner_gene_id"] == "B1"
    assert row["evidence_type"] == "direct_overlap_fixed_edge"

    rejected = run_projection(
        tmp_path,
        backbone,
        [candidate("b", "C2", "P2", 120, 180)],
        [hit("P2", "BX")],
        "rejected",
    )
    row = read_selection(rejected)[0]
    assert row["status"] == "rejected"
    assert row["reason"] == "direct_fixed_partner_lacks_qualifying_base_hit"


def test_unknown_base_hit_fails_closed_without_partial_output(tmp_path: Path) -> None:
    backbone = make_backbone(tmp_path)
    catalog = write_tsv(
        tmp_path / "catalog.tsv", CATALOG_FIELDS, [candidate("a", "C1", "P1", 250, 349)]
    )
    hits = write_tsv(tmp_path / "hits.tsv", HIT_FIELDS, [hit("P1", "UNKNOWN")])
    hit_manifest = write_hits_manifest(
        tmp_path / "hits.manifest.json", backbone, catalog, hits, 1
    )
    output = tmp_path / "projection"
    with pytest.raises(ValueError, match="unknown candidate/base ID"):
        project_candidates_to_fixed_backbone(
            backbone_dir_path=backbone,
            candidate_catalog_tsv_path=catalog,
            fixed_base_hits_tsv_path=hits,
            fixed_base_hits_manifest_path=hit_manifest,
            policy=FixedBackboneProjectionPolicy(),
            output_dir_path=output,
        )
    assert not output.exists()
    assert not Path(f"{output}.working").exists()


def test_cli_exposes_fixed_backbone_and_projection_contracts() -> None:
    parser = build_parser()
    backbone = parser.parse_args(
        [
            "evidence", "build-fixed-target-backbone", "--base-gff", "b.gff3",
            "--query-wgdi-gff", "q.gff", "--base-only-collinearity", "c.tsv",
            "--primary-seqid-table", "p.tsv", "--output-dir", "out",
        ]
    )
    assert backbone.evidence_command == "build-fixed-target-backbone"
    projection = parser.parse_args(
        [
            "evidence", "project-fixed-backbone", "--backbone-dir", "b",
            "--candidate-catalog", "c.tsv", "--fixed-base-hits", "h.tsv",
            "--fixed-base-hits-manifest", "h.json", "--output-dir", "out",
        ]
    )
    assert projection.evidence_command == "project-fixed-backbone"
