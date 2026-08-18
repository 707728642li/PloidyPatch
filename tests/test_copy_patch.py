from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from ploidypatch.copy_patch import (
    compile_copy_addition_patch_edits,
    compile_reviewed_copy_addition_patch_edits,
)
from ploidypatch.cli import main
from ploidypatch.patch import (
    apply_annotation_patch,
    create_annotation_patch,
    revert_annotation_patch,
)


def write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8", newline="")
    return path


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_copy_addition_patch_matches_candidate_and_reverts(tmp_path: Path) -> None:
    source_text = (
        "##gff-version 3\n"
        "chr1\ttest\tgene\t1\t30\t.\t+\t.\tID=G1\n"
        "chr1\ttest\tmRNA\t1\t30\t.\t+\t.\tID=T1;Parent=G1\n"
        "chr1\ttest\tCDS\t1\t30\t.\t+\t0\tParent=T1\n"
    )
    digest = "d" * 64
    appended = (
        f"chr2\tPloidyPatchConsensus\tgene\t101\t130\t.\t+\t.\tID=PPCONS_gene_d;consensus_digest={digest}\n"
        "chr2\tPloidyPatchConsensus\tmRNA\t101\t130\t.\t+\t.\tID=PPCONS_tx_d;Parent=PPCONS_gene_d\n"
        "chr2\tPloidyPatchConsensus\texon\t101\t130\t.\t+\t.\tID=PPCONS_exon_d;Parent=PPCONS_tx_d\n"
        "chr2\tPloidyPatchConsensus\tCDS\t101\t130\t.\t+\t0\tParent=PPCONS_tx_d\n"
    )
    source = write(tmp_path / "source.gff3", source_text)
    candidate = write(tmp_path / "candidate.gff3", source_text + "###\n" + appended)
    edits = tmp_path / "edits.json"

    report = compile_copy_addition_patch_edits(
        annotation_gff_path=source,
        candidate_gff_path=candidate,
        output_edits_json_path=edits,
    )
    assert report["counts"]["copy_addition_events"] == 1
    assert report["policy"]["automatic_approval"] is False

    patch = tmp_path / "patch.json"
    create_annotation_patch(source, edits, patch)
    patched = tmp_path / "patched.gff3"
    apply_annotation_patch(source, patch, patched)
    assert patched.read_bytes() == candidate.read_bytes()

    reverted = tmp_path / "reverted.gff3"
    revert_annotation_patch(patched, patch, reverted)
    assert reverted.read_bytes() == source.read_bytes()


def test_copy_addition_patch_handles_missing_terminal_newline(tmp_path: Path) -> None:
    source_text = "##gff-version 3"
    digest = "e" * 64
    source = write(tmp_path / "source.gff3", source_text)
    candidate = write(
        tmp_path / "candidate.gff3",
        source_text
        + "\n###\n"
        + f"chr2\tPloidyPatchConsensus\tgene\t1\t3\t.\t+\t.\tID=PPCONS_gene_e;consensus_digest={digest}\n"
        + "chr2\tPloidyPatchConsensus\tmRNA\t1\t3\t.\t+\t.\tID=PPCONS_tx_e;Parent=PPCONS_gene_e\n"
        + "chr2\tPloidyPatchConsensus\tCDS\t1\t3\t.\t+\t0\tParent=PPCONS_tx_e\n",
    )
    edits = tmp_path / "edits.json"
    compile_copy_addition_patch_edits(
        annotation_gff_path=source,
        candidate_gff_path=candidate,
        output_edits_json_path=edits,
    )
    patch = tmp_path / "patch.json"
    create_annotation_patch(source, edits, patch)
    patched = tmp_path / "patched.gff3"
    apply_annotation_patch(source, patch, patched)
    assert patched.read_bytes() == candidate.read_bytes()


def reviewed_fixture(tmp_path: Path) -> tuple[Path, Path, Path, Path]:
    source_text = (
        "##gff-version 3\n"
        "chr1\ttest\tgene\t1\t30\t.\t+\t.\tID=G1\n"
        "chr1\ttest\tmRNA\t1\t30\t.\t+\t.\tID=T1;Parent=G1\n"
        "chr1\ttest\tCDS\t1\t30\t.\t+\t0\tParent=T1\n"
    )
    source = write(tmp_path / "source.gff3", source_text)
    digests = {letter: letter * 64 for letter in "abc"}
    blocks = []
    for index, letter in enumerate("abc", start=1):
        digest = digests[letter]
        start = 100 * index
        blocks.append(
            f"chr2\tPloidyPatchConsensus\tgene\t{start}\t{start + 29}\t.\t+\t.\t"
            f"ID=PPCONS_gene_{letter};consensus_digest={digest}\n"
            f"chr2\tPloidyPatchConsensus\tmRNA\t{start}\t{start + 29}\t.\t+\t.\t"
            f"ID=PPCONS_tx_{letter};Parent=PPCONS_gene_{letter}\n"
            f"chr2\tPloidyPatchConsensus\tCDS\t{start}\t{start + 29}\t.\t+\t0\t"
            f"Parent=PPCONS_tx_{letter}\n"
        )
    candidate = write(
        tmp_path / "candidate.gff3", source_text + "###\n" + "".join(blocks)
    )
    decisions = write(
        tmp_path / "decisions.tsv",
        "consensus_digest\tstatus\tconflict_set_digest\tconflict_member_count\n"
        f"{digests['a']}\taccepted\tconflict1\t2\n"
        f"{digests['b']}\taccepted\tconflict1\t2\n"
        f"{digests['c']}\taccepted\t\t1\n",
    )
    manifest = tmp_path / "candidate.gff3.manifest.json"
    manifest.write_text(
        json.dumps(
            {
                "schema_version": "ploidypatch.method_candidate_pool.v2",
                "parameters": {
                    "min_method_support": 1,
                    "max_redundancy_overlap": 0.5,
                    "redundancy_policy": "retain_distinct_chains",
                    "conflict_action": "retain_all_for_ranking_and_review",
                },
                "counts": {
                    "accepted_models": 3,
                    "conflict_sets": 1,
                    "conflicted_chains": 2,
                },
                "outputs": {
                    "candidate_gff": {"sha256": sha256(candidate)},
                    "decisions": {"sha256": sha256(decisions), "rows": 3},
                },
            }
        ),
        encoding="utf-8",
    )
    return source, candidate, decisions, manifest


def test_reviewed_copy_patch_requires_explicit_human_acceptance_and_reverts(
    tmp_path: Path,
) -> None:
    source, candidate, decisions, manifest = reviewed_fixture(tmp_path)
    reviews = write(
        tmp_path / "reviews.tsv",
        "candidate_digest\treview_decision\treviewer\treviewed_at_utc\trationale\n"
        f"{'a' * 64}\taccept_addition\treviewer-1\t2026-08-08T12:00:00Z\tfull-chain evidence\n"
        f"{'b' * 64}\tdefer\treviewer-1\t2026-08-08T12:01:00Z\tconflicting chain\n"
        f"{'c' * 64}\taccept_addition\treviewer-2\t2026-08-08T12:02:00+00:00\tindependent support\n",
    )
    edits = tmp_path / "reviewed.edits.json"
    report = compile_reviewed_copy_addition_patch_edits(
        annotation_gff_path=source,
        candidate_gff_path=candidate,
        pool_decisions_tsv_path=decisions,
        pool_manifest_json_path=manifest,
        review_decisions_tsv_path=reviews,
        output_edits_json_path=edits,
    )
    assert report["counts"]["copy_addition_events"] == 2
    assert report["counts"]["review_decision_counts"] == {
        "accept_addition": 2,
        "defer": 1,
    }
    assert report["policy"]["automatic_approval"] is False
    assert report["policy"]["candidate_selection"] == (
        "explicit_human_accept_addition_only"
    )
    assert all(event["review"]["reviewer"] for event in report["events"])

    patch = tmp_path / "patch.json"
    create_annotation_patch(source, edits, patch)
    patched = tmp_path / "patched.gff3"
    apply_annotation_patch(source, patch, patched)
    text = patched.read_text(encoding="utf-8")
    assert "PPCONS_gene_a" in text and "PPCONS_gene_c" in text
    assert "PPCONS_gene_b" not in text
    reverted = tmp_path / "reverted.gff3"
    revert_annotation_patch(patched, patch, reverted)
    assert reverted.read_bytes() == source.read_bytes()


def test_reviewed_copy_patch_is_exposed_through_public_cli(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    source, candidate, decisions, manifest = reviewed_fixture(tmp_path)
    reviews = write(
        tmp_path / "reviews.tsv",
        "candidate_digest\treview_decision\treviewer\treviewed_at_utc\trationale\n"
        f"{'a' * 64}\taccept_addition\treviewer-1\t2026-08-08T12:00:00Z\tfull-chain evidence\n"
        f"{'b' * 64}\treject\treviewer-1\t2026-08-08T12:01:00Z\tconflicting chain\n",
    )
    edits = tmp_path / "reviewed.edits.json"
    assert (
        main(
            [
                "patch",
                "compile-reviewed-copy-additions",
                "--annotation-gff",
                str(source),
                "--candidate-gff",
                str(candidate),
                "--pool-decisions",
                str(decisions),
                "--pool-manifest",
                str(manifest),
                "--review-decisions",
                str(reviews),
                "--output-edits-json",
                str(edits),
            ]
        )
        == 0
    )
    report = json.loads(capsys.readouterr().out)
    assert report["counts"]["copy_addition_events"] == 1
    assert report["policy"]["automatic_approval"] is False
    assert edits.is_file()


def test_reviewed_copy_patch_rejects_two_conflicting_acceptances(
    tmp_path: Path,
) -> None:
    source, candidate, decisions, manifest = reviewed_fixture(tmp_path)
    reviews = write(
        tmp_path / "reviews.tsv",
        "candidate_digest\treview_decision\treviewer\treviewed_at_utc\trationale\n"
        f"{'a' * 64}\taccept_addition\tr1\t2026-08-08T12:00:00Z\tchain a\n"
        f"{'b' * 64}\taccept_addition\tr2\t2026-08-08T12:01:00Z\tchain b\n",
    )
    with pytest.raises(ValueError, match="mutually exclusive"):
        compile_reviewed_copy_addition_patch_edits(
            annotation_gff_path=source,
            candidate_gff_path=candidate,
            pool_decisions_tsv_path=decisions,
            pool_manifest_json_path=manifest,
            review_decisions_tsv_path=reviews,
            output_edits_json_path=tmp_path / "edits.json",
        )


def test_reviewed_copy_patch_rejects_changed_pool_manifest(tmp_path: Path) -> None:
    source, candidate, decisions, manifest = reviewed_fixture(tmp_path)
    value = json.loads(manifest.read_text(encoding="utf-8"))
    value["outputs"]["candidate_gff"]["sha256"] = "0" * 64
    manifest.write_text(json.dumps(value), encoding="utf-8")
    reviews = write(
        tmp_path / "reviews.tsv",
        "candidate_digest\treview_decision\treviewer\treviewed_at_utc\trationale\n"
        f"{'a' * 64}\taccept_addition\tr1\t2026-08-08T12:00:00Z\tchain a\n",
    )
    with pytest.raises(ValueError, match="schema or checksum"):
        compile_reviewed_copy_addition_patch_edits(
            annotation_gff_path=source,
            candidate_gff_path=candidate,
            pool_decisions_tsv_path=decisions,
            pool_manifest_json_path=manifest,
            review_decisions_tsv_path=reviews,
            output_edits_json_path=tmp_path / "edits.json",
        )


def test_reviewed_copy_patch_requires_utc_review_timestamp(tmp_path: Path) -> None:
    source, candidate, decisions, manifest = reviewed_fixture(tmp_path)
    reviews = write(
        tmp_path / "reviews.tsv",
        "candidate_digest\treview_decision\treviewer\treviewed_at_utc\trationale\n"
        f"{'a' * 64}\taccept_addition\tr1\t2026-08-08T20:00:00+08:00\tchain a\n",
    )
    with pytest.raises(ValueError, match="not UTC"):
        compile_reviewed_copy_addition_patch_edits(
            annotation_gff_path=source,
            candidate_gff_path=candidate,
            pool_decisions_tsv_path=decisions,
            pool_manifest_json_path=manifest,
            review_decisions_tsv_path=reviews,
            output_edits_json_path=tmp_path / "edits.json",
        )


def test_copy_patch_compilers_refuse_stale_working_files(tmp_path: Path) -> None:
    source, candidate, decisions, manifest = reviewed_fixture(tmp_path)
    output = tmp_path / "edits.json"
    Path(f"{output}.working").write_text("interrupted\n", encoding="utf-8")
    with pytest.raises(FileExistsError, match="Refusing to overwrite"):
        compile_copy_addition_patch_edits(
            annotation_gff_path=source,
            candidate_gff_path=candidate,
            output_edits_json_path=output,
        )
    reviews = write(
        tmp_path / "reviews.tsv",
        "candidate_digest\treview_decision\treviewer\treviewed_at_utc\trationale\n"
        f"{'a' * 64}\taccept_addition\tr1\t2026-08-08T12:00:00Z\tchain a\n",
    )
    with pytest.raises(FileExistsError, match="Refusing to overwrite"):
        compile_reviewed_copy_addition_patch_edits(
            annotation_gff_path=source,
            candidate_gff_path=candidate,
            pool_decisions_tsv_path=decisions,
            pool_manifest_json_path=manifest,
            review_decisions_tsv_path=reviews,
            output_edits_json_path=output,
        )
