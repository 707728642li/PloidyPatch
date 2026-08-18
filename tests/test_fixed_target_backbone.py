from __future__ import annotations

import csv
import hashlib
import inspect
import json
from pathlib import Path

import pytest

from ploidypatch.fixed_target_backbone import (
    FixedTargetBackbonePolicy,
    build_fixed_target_backbone,
)


GENES = (
    ("chr1", "A1", 10, 99, "+", 1),
    ("chr1", "A2", 200, 299, "+", 2),
    ("chr1", "A3", 400, 499, "+", 3),
    ("chr2", "B1", 1010, 1099, "+", 1),
    ("chr2", "B2", 1200, 1299, "+", 2),
    ("chr2", "B3", 1400, 1499, "+", 3),
)


def write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8", newline="")
    return path


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def base_gff_text(*, reverse: bool = False) -> str:
    groups: list[str] = []
    for seqid, gene_id, start, end, strand, _ in GENES:
        transcript_id = f"{gene_id}.t1"
        groups.append(
            f"{seqid}\ttest\tgene\t{start}\t{end}\t.\t{strand}\t.\tID={gene_id}\n"
            f"{seqid}\ttest\tmRNA\t{start}\t{end}\t.\t{strand}\t.\t"
            f"ID={transcript_id};Parent={gene_id}\n"
            f"{seqid}\ttest\tCDS\t{start + 9}\t{end - 9}\t.\t{strand}\t0\t"
            f"ID={gene_id}.cds;Parent={transcript_id}\n"
        )
    if reverse:
        groups.reverse()
    return "##gff-version 3\n" + "".join(groups)


def prefixed_base_gff_text() -> str:
    groups: list[str] = []
    for seqid, gene_id, start, end, strand, _ in GENES:
        feature_id = f"gene:{gene_id}"
        transcript_id = f"{feature_id}.t1"
        groups.append(
            f"{seqid}\ttest\tgene\t{start}\t{end}\t.\t{strand}\t.\t"
            f"ID={feature_id};gene_id={gene_id}\n"
            f"{seqid}\ttest\tmRNA\t{start}\t{end}\t.\t{strand}\t.\t"
            f"ID={transcript_id};Parent={feature_id}\n"
            f"{seqid}\ttest\tCDS\t{start + 9}\t{end - 9}\t.\t{strand}\t0\t"
            f"ID={feature_id}.cds;Parent={transcript_id}\n"
        )
    return "##gff-version 3\n" + "".join(groups)


def query_gff_text(*, reverse: bool = False) -> str:
    rows = [
        f"{seqid}\t{gene_id}\t{start}\t{end}\t{strand}\t{order}\n"
        for seqid, gene_id, start, end, strand, order in GENES
    ]
    if reverse:
        rows.reverse()
    return "".join(rows)


def block(
    block_id: str,
    pairs: list[tuple[str, str]],
    *,
    query_seqid: str = "chr1",
    target_seqid: str = "chr2",
    orientation: str = "plus",
) -> str:
    lines = [
        f"# Alignment {block_id}: score=100 pvalue=0.01 N={len(pairs)} "
        f"{query_seqid}&{target_seqid} {orientation}\n"
    ]
    lines.extend(f"{left} 1 {right} 1 1\n" for left, right in pairs)
    return "".join(lines)


def inputs(
    tmp_path: Path,
    collinearity: str,
    *,
    reverse_base: bool = False,
    reverse_query: bool = False,
) -> tuple[Path, Path, Path]:
    return (
        write(tmp_path / "base.gff3", base_gff_text(reverse=reverse_base)),
        write(tmp_path / "query.gff", query_gff_text(reverse=reverse_query)),
        write(tmp_path / "self.tsv", collinearity),
    )


def policy() -> FixedTargetBackbonePolicy:
    return FixedTargetBackbonePolicy(
        primary_seqids=("chr2", "chr1"),
        min_block_pairs=2,
    )


def build(
    tmp_path: Path,
    collinearity: str,
    *,
    output_name: str = "backbone",
    reverse_base: bool = False,
    reverse_query: bool = False,
):
    base, query, self_wgdi = inputs(
        tmp_path,
        collinearity,
        reverse_base=reverse_base,
        reverse_query=reverse_query,
    )
    output = tmp_path / output_name
    manifest = build_fixed_target_backbone(
        base_gff_path=base,
        query_wgdi_gff_path=query,
        base_only_collinearity_path=self_wgdi,
        policy=policy(),
        output_dir_path=output,
    )
    return output, manifest, (base, query, self_wgdi)


def test_builder_has_no_candidate_input_and_binds_base_only_inputs(tmp_path: Path) -> None:
    assert "candidate" not in inspect.signature(build_fixed_target_backbone).parameters
    collinearity = block("1", [("A1", "B1"), ("A2", "B2"), ("A3", "B3")])
    output, manifest, source_paths = build(tmp_path, collinearity)

    assert manifest["truth_access"] is False
    assert manifest["candidate_access"] is False
    assert manifest["semantic_invariants"] == {
        "input_permutation_deterministic": True,
        "reverse_duplicate_statistics_audit_only": True,
        "ambiguous_endpoint_block_quarantine": True,
    }
    assert set(manifest["inputs"]) == {
        "base_gff",
        "query_wgdi_gff",
        "base_only_collinearity",
    }
    assert all("candidate" not in key for key in manifest["inputs"])
    for key, path in zip(manifest["inputs"], source_paths, strict=True):
        assert manifest["inputs"][key]["sha256"] == sha256(path)
    for key in ("genes", "edges", "cells"):
        artifact = output / manifest["outputs"][key]["file_name"]
        assert sha256(artifact) == manifest["outputs"][key]["sha256"]
    assert json.loads((output / "manifest.json").read_text(encoding="utf-8")) == manifest
    assert (output / "SHA256SUMS").is_file()


def test_plus_backbone_emits_reciprocal_edges_and_bidirectional_cells(
    tmp_path: Path,
) -> None:
    collinearity = block("1", [("A1", "B1"), ("A2", "B2"), ("A3", "B3")])
    output, manifest, _ = build(tmp_path, collinearity)
    genes = read_tsv(output / "genes.tsv")
    edges = read_tsv(output / "edges.tsv")
    cells = read_tsv(output / "cells.tsv")

    assert manifest["counts"]["genes"] == 6
    assert manifest["counts"]["accepted_edges"] == 3
    assert manifest["counts"]["directed_cells"] == 4
    assert {row["status"] for row in edges} == {"accepted"}
    assert {row["orientation"] for row in cells} == {"plus"}
    assert {(row["source_seqid"], row["partner_seqid"]) for row in cells} == {
        ("chr1", "chr2"),
        ("chr2", "chr1"),
    }
    a1 = next(row for row in genes if row["gene_id"] == "A1")
    assert a1["midpoint2"] == str(10 + 99)
    assert a1["coding_start"] == "19"
    assert a1["coding_end"] == "90"


def test_minus_backbone_preserves_local_orientation_in_both_directions(
    tmp_path: Path,
) -> None:
    collinearity = block(
        "1",
        [("A1", "B3"), ("A2", "B2"), ("A3", "B1")],
        orientation="minus",
    )
    output, manifest, _ = build(tmp_path, collinearity)
    cells = read_tsv(output / "cells.tsv")

    assert manifest["counts"]["directed_cells"] == 4
    assert {row["orientation"] for row in cells} == {"minus"}
    for row in cells:
        assert int(row["source_left_order"]) < int(row["source_right_order"])
        assert int(row["partner_at_source_left_order"]) > int(
            row["partner_at_source_right_order"]
        )


def test_reverse_duplicate_and_input_permutation_leave_semantic_outputs_identical(
    tmp_path: Path,
) -> None:
    forward_pairs = [("A1", "B1"), ("A2", "B2"), ("A3", "B3")]
    reverse_pairs = [(right, left) for left, right in reversed(forward_pairs)]
    duplicate = block("1", forward_pairs) + block(
        "2",
        reverse_pairs,
        query_seqid="chr2",
        target_seqid="chr1",
    )
    first_root = tmp_path / "first"
    second_root = tmp_path / "second"
    first_root.mkdir()
    second_root.mkdir()
    first, first_manifest, _ = build(first_root, duplicate)
    second, second_manifest, _ = build(
        second_root,
        block("9", list(reversed(forward_pairs))),
        reverse_base=True,
        reverse_query=True,
    )

    assert first_manifest["backbone_id"] == second_manifest["backbone_id"]
    assert first_manifest["counts"]["reverse_or_exact_duplicate_blocks_collapsed"] == 1
    for name in ("genes.tsv", "edges.tsv", "cells.tsv"):
        assert (first / name).read_bytes() == (second / name).read_bytes()


@pytest.mark.parametrize(
    ("query_text", "collinearity", "message"),
    [
        (
            query_gff_text().replace("chr1\tA2\t200\t299\t+\t2", "chr1\tA2\t200\t299\t+\t1"),
            block("1", [("A1", "B1"), ("A2", "B2")]),
            "Invalid or duplicate WGDI query gene",
        ),
        (
            query_gff_text(),
            "# not a WGDI header\nA1 1 B1 1 1\n",
            "Malformed WGDI block header",
        ),
    ],
)
def test_malformed_or_duplicate_inputs_fail_without_partial_output(
    tmp_path: Path,
    query_text: str,
    collinearity: str,
    message: str,
) -> None:
    base = write(tmp_path / "base.gff3", base_gff_text())
    query = write(tmp_path / "query.gff", query_text)
    self_wgdi = write(tmp_path / "self.tsv", collinearity)
    output = tmp_path / "backbone"
    with pytest.raises(ValueError, match=message):
        build_fixed_target_backbone(
            base_gff_path=base,
            query_wgdi_gff_path=query,
            base_only_collinearity_path=self_wgdi,
            policy=policy(),
            output_dir_path=output,
        )
    assert not output.exists()
    assert not Path(f"{output}.working").exists()


@pytest.mark.parametrize("alias_key", ["gene_id", "locus_tag", "Name"])
def test_exact_unique_alias_maps_wgdi_ids_to_feature_ids(
    tmp_path: Path, alias_key: str
) -> None:
    base_text = prefixed_base_gff_text().replace(";gene_id=", f";{alias_key}=")
    base = write(tmp_path / "base.gff3", base_text)
    query = write(tmp_path / "query.gff", query_gff_text())
    self_wgdi = write(
        tmp_path / "self.tsv",
        block("1", [("A1", "B1"), ("A2", "B2"), ("A3", "B3")]),
    )
    output = tmp_path / "backbone"
    manifest = build_fixed_target_backbone(
        base_gff_path=base,
        query_wgdi_gff_path=query,
        base_only_collinearity_path=self_wgdi,
        policy=policy(),
        output_dir_path=output,
    )
    genes = read_tsv(output / "genes.tsv")
    assert {row["gene_id"] for row in genes} == {
        "gene:A1", "gene:A2", "gene:A3", "gene:B1", "gene:B2", "gene:B3"
    }
    assert {row["wgdi_gene_id"] for row in genes} == {
        "A1", "A2", "A3", "B1", "B2", "B3"
    }
    assert manifest["identifier_mapping"]["mapping_mode_counts"] == {
        "exact_attribute": 6
    }


def test_exact_alias_must_be_unique_in_both_directions(tmp_path: Path) -> None:
    base_text = prefixed_base_gff_text().replace(
        "ID=gene:A1;gene_id=A1", "ID=gene:A1;gene_id=A1;Name=A2"
    ).replace("ID=gene:A2;gene_id=A2", "ID=gene:A2;unused=A2")
    base = write(tmp_path / "base.gff3", base_text)
    query = write(tmp_path / "query.gff", query_gff_text())
    self_wgdi = write(
        tmp_path / "self.tsv",
        block("1", [("A1", "B1"), ("A2", "B2"), ("A3", "B3")]),
    )
    with pytest.raises(ValueError, match="Multiple WGDI gene IDs map to one"):
        build_fixed_target_backbone(
            base_gff_path=base,
            query_wgdi_gff_path=query,
            base_only_collinearity_path=self_wgdi,
            policy=policy(),
            output_dir_path=tmp_path / "backbone",
        )


def test_one_wgdi_alias_matching_multiple_features_is_rejected(tmp_path: Path) -> None:
    extra = (
        "chr1\ttest\tgene\t600\t699\t.\t+\t.\tID=extra;Name=A1\n"
        "chr1\ttest\tmRNA\t600\t699\t.\t+\t.\tID=extra.t1;Parent=extra\n"
        "chr1\ttest\tCDS\t609\t690\t.\t+\t0\tParent=extra.t1\n"
    )
    base = write(tmp_path / "base.gff3", prefixed_base_gff_text() + extra)
    query = write(tmp_path / "query.gff", query_gff_text())
    self_wgdi = write(
        tmp_path / "self.tsv",
        block("1", [("A1", "B1"), ("A2", "B2"), ("A3", "B3")]),
    )
    with pytest.raises(ValueError, match="do not map uniquely"):
        build_fixed_target_backbone(
            base_gff_path=base,
            query_wgdi_gff_path=query,
            base_only_collinearity_path=self_wgdi,
            policy=policy(),
            output_dir_path=tmp_path / "backbone",
        )


def test_ambiguous_anchor_rows_are_removed_before_block_qualification(tmp_path: Path) -> None:
    collinearity = block(
        "1", [("A1", "B1"), ("A1", "B2"), ("A2", "B3"), ("A3", "B2")]
    )
    base, query, self_wgdi = inputs(tmp_path, collinearity)
    output = tmp_path / "backbone"
    manifest = build_fixed_target_backbone(
        base_gff_path=base,
        query_wgdi_gff_path=query,
        base_only_collinearity_path=self_wgdi,
        policy=policy(),
        output_dir_path=output,
    )
    assert manifest["applicable"] is False
    assert manifest["counts"]["ambiguous_endpoint_blocks_quarantined"] == 1
    assert read_tsv(output / "edges.tsv") == []
    assert read_tsv(output / "cells.tsv") == []


def test_ambiguous_endpoint_block_is_quarantined_without_altering_clean_blocks(
    tmp_path: Path,
) -> None:
    clean_pairs = [("A1", "B1"), ("A2", "B2"), ("A3", "B3")]
    ambiguous = block("1", [("A1", "B1"), ("A1", "B2")])
    first_root = tmp_path / "with_quarantine"
    second_root = tmp_path / "clean_only"
    first_root.mkdir()
    second_root.mkdir()
    first, first_manifest, _ = build(first_root, ambiguous + block("2", clean_pairs))
    second, second_manifest, _ = build(second_root, block("9", clean_pairs))

    assert first_manifest["applicable"] is True
    assert first_manifest["counts"]["ambiguous_endpoint_blocks_quarantined"] == 1
    assert first_manifest["backbone_id"] == second_manifest["backbone_id"]
    for name in ("genes.tsv", "edges.tsv", "cells.tsv"):
        assert (first / name).read_bytes() == (second / name).read_bytes()


def test_duplicate_geometry_with_conflicting_statistics_is_audit_only(
    tmp_path: Path,
) -> None:
    pairs = [("A1", "B1"), ("A2", "B2"), ("A3", "B3")]
    conflicting = block("1", pairs) + block(
        "2", [(right, left) for left, right in reversed(pairs)],
        query_seqid="chr2", target_seqid="chr1",
    ).replace("score=100", "score=99").replace("pvalue=0.01", "pvalue=0.02")
    first_root = tmp_path / "first"
    second_root = tmp_path / "second"
    first_root.mkdir()
    second_root.mkdir()
    first, first_manifest, _ = build(first_root, conflicting)
    second, second_manifest, _ = build(second_root, block("9", pairs))
    assert first_manifest["backbone_id"] == second_manifest["backbone_id"]
    assert first_manifest["counts"]["duplicate_geometry_statistic_discordances"] == 1
    statistic_audit = first_manifest["block_statistic_audit"]
    assert statistic_audit["semantic_use"] is False
    assert statistic_audit["duplicate_geometry_discordant_group_count"] == 1
    assert statistic_audit["score_min"] == "99"
    assert statistic_audit["score_max"] == "100"
    assert statistic_audit["pvalue_min"] == "0.01"
    assert statistic_audit["pvalue_max"] == "0.02"
    for name in ("genes.tsv", "edges.tsv", "cells.tsv"):
        assert (first / name).read_bytes() == (second / name).read_bytes()
    assert "score" not in read_tsv(first / "edges.tsv")[0]
    assert "score" not in read_tsv(first / "cells.tsv")[0]


def test_normalized_feature_id_fallback_is_unique_and_audited(tmp_path: Path) -> None:
    base = write(
        tmp_path / "base.gff3",
        prefixed_base_gff_text().replace(";gene_id=A", ";unused=A").replace(
            ";gene_id=B", ";unused=B"
        ),
    )
    query = write(tmp_path / "query.gff", query_gff_text())
    self_wgdi = write(
        tmp_path / "self.tsv",
        block("1", [("A1", "B1"), ("A2", "B2"), ("A3", "B3")]),
    )
    manifest = build_fixed_target_backbone(
        base_gff_path=base,
        query_wgdi_gff_path=query,
        base_only_collinearity_path=self_wgdi,
        policy=policy(),
        output_dir_path=tmp_path / "backbone",
    )
    assert manifest["identifier_mapping"]["mapping_mode_counts"] == {
        "normalized_feature_ID_fallback": 6
    }


def test_manifest_emits_pre_candidate_applicability_metrics(tmp_path: Path) -> None:
    output, manifest, _ = build(
        tmp_path, block("1", [("A1", "B1"), ("A2", "B2"), ("A3", "B3")])
    )
    metrics = manifest["applicability_metrics"]
    assert metrics["primary_gene_midpoint_backbone_coverage"] == 1.0
    assert metrics["primary_chromosome_cell_coverage_fraction"] == 1.0
    assert metrics["minimum_cells_per_covered_chromosome_observed"] == 2
    assert metrics["unique_partner_chromosome_fraction_among_covered_genes"] == 1.0


def test_refuses_symlink_input(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    base, query, self_wgdi = inputs(
        tmp_path,
        block("1", [("A1", "B1"), ("A2", "B2")]),
    )
    original_is_symlink = Path.is_symlink
    monkeypatch.setattr(
        Path,
        "is_symlink",
        lambda path: path == base.absolute() or original_is_symlink(path),
    )
    with pytest.raises(ValueError, match="symbolic link"):
        build_fixed_target_backbone(
            base_gff_path=base,
            query_wgdi_gff_path=query,
            base_only_collinearity_path=self_wgdi,
            policy=policy(),
            output_dir_path=tmp_path / "backbone",
        )


def test_refuses_overwrite_and_leaves_no_working_directory(tmp_path: Path) -> None:
    collinearity = block("1", [("A1", "B1"), ("A2", "B2"), ("A3", "B3")])
    output, _, source_paths = build(tmp_path, collinearity)
    with pytest.raises(FileExistsError, match="Refusing to overwrite"):
        build_fixed_target_backbone(
            base_gff_path=source_paths[0],
            query_wgdi_gff_path=source_paths[1],
            base_only_collinearity_path=source_paths[2],
            policy=policy(),
            output_dir_path=output,
        )
    assert not Path(f"{output}.working").exists()
