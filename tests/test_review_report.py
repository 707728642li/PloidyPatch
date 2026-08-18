from __future__ import annotations

import csv
import copy
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

import pytest

from ploidypatch.artifact_manifest import verify_sha256sums
from ploidypatch.copy_patch import compile_reviewed_copy_addition_patch_edits
from ploidypatch.review_report import _render_html, build_review_report


ROOT = Path(__file__).resolve().parents[1]
EXAMPLE = ROOT / "examples/minimal_reviewed_patch"
DIGESTS = ("a" * 64, "b" * 64, "c" * 64)


def _write_table(path: Path, fields: tuple[str, ...], rows: list[dict[str, object]]) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _optional_inputs(
    tmp_path: Path, *, automatic_digest: str | None = None
) -> tuple[Path, Path, Path, Path]:
    scores = tmp_path / "scores.tsv"
    _write_table(
        scores,
        (
            "candidate_digest",
            "v04_primary_rank_score",
            "v04_primary_rank_percentile",
            "v04_automatic_approval",
        ),
        [
            {
                "candidate_digest": digest,
                "v04_primary_rank_score": score,
                "v04_primary_rank_percentile": percentile,
                "v04_automatic_approval": int(digest == automatic_digest),
            }
            for digest, score, percentile in zip(
                DIGESTS, (0.91, 0.62, 0.79), (1.0, 1 / 3, 2 / 3), strict=True
            )
        ],
    )
    features = tmp_path / "copy_features.tsv"
    _write_table(
        features,
        ("candidate_digest", "has_miniprot", "has_gemoma", "has_lifton"),
        [
            {"candidate_digest": DIGESTS[0], "has_miniprot": 1, "has_gemoma": 1, "has_lifton": 1},
            {"candidate_digest": DIGESTS[1], "has_miniprot": 1, "has_gemoma": 0, "has_lifton": 1},
            {"candidate_digest": DIGESTS[2], "has_miniprot": 0, "has_gemoma": 1, "has_lifton": 1},
        ],
    )
    topology = tmp_path / "topology.tsv"
    _write_table(
        topology,
        (
            "candidate_digest",
            "topology_available",
            "topology_reason",
            "partner_gene_id",
            "topology_coherence_score",
        ),
        [
            {
                "candidate_digest": DIGESTS[0],
                "topology_available": 1,
                "topology_reason": "accepted_existing_wgd_partner",
                "partner_gene_id": "BASE_homeolog_a",
                "topology_coherence_score": 0.88,
            },
            {
                "candidate_digest": DIGESTS[1],
                "topology_available": 0,
                "topology_reason": "nonreciprocal_multiple_partners",
                "partner_gene_id": "",
                "topology_coherence_score": "",
            },
            {
                "candidate_digest": DIGESTS[2],
                "topology_available": 1,
                "topology_reason": "accepted_existing_wgd_partner",
                "partner_gene_id": "BASE_homeolog_c",
                "topology_coherence_score": 0.73,
            },
        ],
    )
    patch = tmp_path / "reviewed.edits.json"
    compile_reviewed_copy_addition_patch_edits(
        annotation_gff_path=EXAMPLE / "source.gff3",
        candidate_gff_path=EXAMPLE / "candidate.gff3",
        pool_decisions_tsv_path=EXAMPLE / "decisions.tsv",
        pool_manifest_json_path=EXAMPLE / "candidate.gff3.manifest.json",
        review_decisions_tsv_path=EXAMPLE / "review_decisions.tsv",
        output_edits_json_path=patch,
    )
    return scores, features, topology, patch


def test_report_builds_a_ready_self_contained_review_package(tmp_path: Path) -> None:
    scores, features, topology, patch = _optional_inputs(tmp_path)
    output = tmp_path / "report"
    report = build_review_report(
        candidate_gff_path=EXAMPLE / "candidate.gff3",
        pool_decisions_tsv_path=EXAMPLE / "decisions.tsv",
        pool_manifest_json_path=EXAMPLE / "candidate.gff3.manifest.json",
        review_decisions_tsv_path=EXAMPLE / "review_decisions.tsv",
        scores_tsv_path=scores,
        copy_features_tsv_path=features,
        topology_features_tsv_path=topology,
        patch_edits_json_path=patch,
        output_dir=output,
        title="Minimal reviewed patch",
    )
    assert report["truth_access"] is False
    assert report["summary"]["state"] == "ready_to_apply"
    assert report["counts"]["candidates"] == 3
    assert report["counts"]["reviewed"] == 3
    assert report["counts"]["accepted_review"] == 2
    assert report["counts"]["patch_events"] == 2
    assert report["counts"]["automatic_approval"] == 0
    assert report["counts"]["automatic_approval_checked"] is True
    assert report["counts"]["topology_available"] == 2
    assert report["method_overlap"]["all_three"] == 1
    assert report["distributions"]["by_seqid"] == [
        {"seqid": "chr2", "count": 2},
        {"seqid": "chr3", "count": 1},
    ]
    assert sum(
        item["count"] for item in report["distributions"]["score_histogram"]["bins"]
    ) == 3
    assert {path.name for path in output.iterdir()} == {
        "index.html",
        "report.json",
        "candidates.tsv",
        "review_ledger_template.tsv",
        "SHA256SUMS",
    }
    verify_sha256sums(output, ignore_checksum_file=True)
    rendered = (output / "index.html").read_text(encoding="utf-8")
    assert "Candidate-to-patch funnel" in rendered
    assert "Projection-method intersections" in rendered
    assert "Genomic candidate distribution" in rendered
    assert "Review-priority distribution" in rendered
    assert "Candidate explorer" in rendered
    assert "Integrity and provenance" in rendered
    assert "fetch(" not in rendered
    assert "cdn." not in rendered
    assert "review_rank_not_probability" in rendered
    assert "Content-Security-Policy" in rendered
    assert "review_ledger_template.tsv" in rendered
    candidates = report["candidates"]
    assert candidates[0]["score"] == pytest.approx(0.91)
    candidate_a = next(row for row in candidates if row["candidate_digest"] == DIGESTS[0])
    assert candidate_a["review"]["decision"] == "accept_addition"
    assert candidate_a["transcripts"][0]["cds"] == [[101, 190, "0"]]
    round_trip_edits = tmp_path / "ledger_round_trip.edits.json"
    compiled = compile_reviewed_copy_addition_patch_edits(
        annotation_gff_path=EXAMPLE / "source.gff3",
        candidate_gff_path=EXAMPLE / "candidate.gff3",
        pool_decisions_tsv_path=EXAMPLE / "decisions.tsv",
        pool_manifest_json_path=EXAMPLE / "candidate.gff3.manifest.json",
        review_decisions_tsv_path=output / "review_ledger_template.tsv",
        output_edits_json_path=round_trip_edits,
    )
    assert compiled["counts"]["copy_addition_events"] == 2


def test_report_without_optional_evidence_explains_absence(tmp_path: Path) -> None:
    output = tmp_path / "report"
    report = build_review_report(
        candidate_gff_path=EXAMPLE / "candidate.gff3",
        pool_decisions_tsv_path=EXAMPLE / "decisions.tsv",
        pool_manifest_json_path=EXAMPLE / "candidate.gff3.manifest.json",
        output_dir=output,
    )
    assert report["summary"]["state"] == "review_in_progress"
    assert report["counts"]["reviewed"] == 0
    assert report["counts"]["topology_available"] is None
    assert report["counts"]["automatic_approval_checked"] is False
    assert report["method_overlap"]["available"] is False
    assert any(
        note["title"] == "Evolutionary context" and "not supplied" in note["text"]
        for note in report["interpretation"]
    )


def test_report_rejects_evaluator_fields_and_universe_drift(tmp_path: Path) -> None:
    contaminated = tmp_path / "scores.tsv"
    _write_table(
        contaminated,
        ("candidate_digest", "label_exact_cds", "v04_primary_rank_score"),
        [
            {"candidate_digest": digest, "label_exact_cds": 1, "v04_primary_rank_score": 0.5}
            for digest in DIGESTS
        ],
    )
    with pytest.raises(ValueError, match="evaluator-only"):
        build_review_report(
            candidate_gff_path=EXAMPLE / "candidate.gff3",
            pool_decisions_tsv_path=EXAMPLE / "decisions.tsv",
            pool_manifest_json_path=EXAMPLE / "candidate.gff3.manifest.json",
            scores_tsv_path=contaminated,
            output_dir=tmp_path / "report",
        )
    assert not (tmp_path / "report").exists()
    harmless = tmp_path / "harmless.tsv"
    _write_table(
        harmless,
        ("candidate_digest", "chr_label"),
        [
            {"candidate_digest": digest, "chr_label": chromosome}
            for digest, chromosome in zip(DIGESTS, ("chr2", "chr2", "chr3"), strict=True)
        ],
    )
    report = build_review_report(
        candidate_gff_path=EXAMPLE / "candidate.gff3",
        pool_decisions_tsv_path=EXAMPLE / "decisions.tsv",
        pool_manifest_json_path=EXAMPLE / "candidate.gff3.manifest.json",
        scores_tsv_path=harmless,
        output_dir=tmp_path / "harmless-report",
    )
    assert report["truth_access"] is False


def test_report_cli_emits_compact_machine_summary(tmp_path: Path) -> None:
    output = tmp_path / "report"
    environment = os.environ.copy()
    environment["PYTHONPATH"] = str(ROOT / "src")
    completed = subprocess.run(
        [
            sys.executable,
            "-m",
            "ploidypatch.cli",
            "report",
            "--candidate-gff",
            str(EXAMPLE / "candidate.gff3"),
            "--pool-decisions",
            str(EXAMPLE / "decisions.tsv"),
            "--pool-manifest",
            str(EXAMPLE / "candidate.gff3.manifest.json"),
            "--review-decisions",
            str(EXAMPLE / "review_decisions.tsv"),
            "--output-dir",
            str(output),
            "--title",
            "CLI report smoke",
        ],
        cwd=ROOT,
        env=environment,
        text=True,
        capture_output=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr
    summary = json.loads(completed.stdout)
    assert summary["state"] == "review_complete_patch_pending"
    assert summary["counts"]["candidates"] == 3
    assert (output / "index.html").is_file()


def test_report_cli_can_fail_on_attention(tmp_path: Path) -> None:
    rows: list[dict[str, object]] = []
    with (EXAMPLE / "review_decisions.tsv").open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["candidate_digest"] == DIGESTS[1]:
                row["review_decision"] = "accept_addition"
            rows.append(row)
    review = tmp_path / "conflicting_review.tsv"
    _write_table(review, REVIEW_FIELDS_FOR_TEST, rows)
    environment = os.environ.copy()
    environment["PYTHONPATH"] = str(ROOT / "src")
    completed = subprocess.run(
        [
            sys.executable,
            "-m",
            "ploidypatch.cli",
            "report",
            "--candidate-gff",
            str(EXAMPLE / "candidate.gff3"),
            "--pool-decisions",
            str(EXAMPLE / "decisions.tsv"),
            "--pool-manifest",
            str(EXAMPLE / "candidate.gff3.manifest.json"),
            "--review-decisions",
            str(review),
            "--output-dir",
            str(tmp_path / "attention-report"),
            "--fail-on-attention",
        ],
        cwd=ROOT,
        env=environment,
        text=True,
        capture_output=True,
        check=False,
    )
    assert completed.returncode == 2, completed.stderr
    assert json.loads(completed.stdout)["state"] == "attention_required"


def test_report_validated_state_requires_matching_reversible_run(tmp_path: Path) -> None:
    scores, features, topology, patch = _optional_inputs(tmp_path)
    run_summary = tmp_path / "run_summary.json"
    run_summary.write_text(
        json.dumps(
            {
                "accepted_additions": 2,
                "automatic_approval": False,
                "byte_identical_reversion": True,
                "source_sha256": "1" * 64,
                "reverted_sha256": "1" * 64,
                "selected_candidate_digests": [DIGESTS[0], DIGESTS[2]],
            }
        ),
        encoding="utf-8",
    )
    report = build_review_report(
        candidate_gff_path=EXAMPLE / "candidate.gff3",
        pool_decisions_tsv_path=EXAMPLE / "decisions.tsv",
        pool_manifest_json_path=EXAMPLE / "candidate.gff3.manifest.json",
        review_decisions_tsv_path=EXAMPLE / "review_decisions.tsv",
        scores_tsv_path=scores,
        copy_features_tsv_path=features,
        topology_features_tsv_path=topology,
        patch_edits_json_path=patch,
        run_summary_json_path=run_summary,
        output_dir=tmp_path / "report",
    )
    assert report["summary"]["state"] == "validated_reversible_run"
    assert report["summary"]["detail"].startswith("Archive this checksum-bound")


def test_report_surfaces_automatic_approval_and_patch_mismatch(tmp_path: Path) -> None:
    scores, _, _, patch = _optional_inputs(tmp_path, automatic_digest=DIGESTS[0])
    patch_data = json.loads(patch.read_text(encoding="utf-8"))
    patch_data["events"] = patch_data["events"][:1]
    mismatched_patch = tmp_path / "mismatched.edits.json"
    mismatched_patch.write_text(json.dumps(patch_data), encoding="utf-8")
    report = build_review_report(
        candidate_gff_path=EXAMPLE / "candidate.gff3",
        pool_decisions_tsv_path=EXAMPLE / "decisions.tsv",
        pool_manifest_json_path=EXAMPLE / "candidate.gff3.manifest.json",
        review_decisions_tsv_path=EXAMPLE / "review_decisions.tsv",
        scores_tsv_path=scores,
        patch_edits_json_path=mismatched_patch,
        output_dir=tmp_path / "report",
    )
    assert report["summary"]["state"] == "attention_required"
    assert report["counts"]["automatic_approval"] == 1
    patch_note = next(note for note in report["interpretation"] if note["title"] == "Patch verification")
    assert "does not match" in patch_note["text"]


def test_report_surfaces_conflicting_human_acceptance(tmp_path: Path) -> None:
    rows: list[dict[str, object]] = []
    with (EXAMPLE / "review_decisions.tsv").open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["candidate_digest"] == DIGESTS[1]:
                row["review_decision"] = "accept_addition"
                row["rationale"] = "Deliberate test collision"
            rows.append(row)
    review = tmp_path / "conflicting_review.tsv"
    _write_table(review, REVIEW_FIELDS_FOR_TEST, rows)
    report = build_review_report(
        candidate_gff_path=EXAMPLE / "candidate.gff3",
        pool_decisions_tsv_path=EXAMPLE / "decisions.tsv",
        pool_manifest_json_path=EXAMPLE / "candidate.gff3.manifest.json",
        review_decisions_tsv_path=review,
        output_dir=tmp_path / "report",
    )
    assert report["summary"]["state"] == "attention_required"
    conflict_note = next(note for note in report["interpretation"] if note["title"] == "Conflict safety")
    assert conflict_note["level"] == "critical"


def test_report_escapes_untrusted_display_text_and_has_exact_checksum_universe(
    tmp_path: Path,
) -> None:
    rows: list[dict[str, object]] = []
    injection = "</script><img src=x onerror=alert(1)>"
    with (EXAMPLE / "review_decisions.tsv").open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            row["rationale"] = injection
            rows.append(row)
    review = tmp_path / "safe_review.tsv"
    _write_table(review, REVIEW_FIELDS_FOR_TEST, rows)
    output = tmp_path / "report"
    build_review_report(
        candidate_gff_path=EXAMPLE / "candidate.gff3",
        pool_decisions_tsv_path=EXAMPLE / "decisions.tsv",
        pool_manifest_json_path=EXAMPLE / "candidate.gff3.manifest.json",
        review_decisions_tsv_path=review,
        output_dir=output,
        title="Safe <report> & review __PLOIDYPATCH_VERSION__ __REPORT_PAYLOAD__",
    )
    rendered = (output / "index.html").read_text(encoding="utf-8")
    assert injection not in rendered
    assert (
        "Safe &lt;report&gt; &amp; review __PLOIDYPATCH_VERSION__ __REPORT_PAYLOAD__"
        in rendered
    )
    checksum_lines = (output / "SHA256SUMS").read_text(encoding="utf-8").splitlines()
    assert {line.split("  ", 1)[1] for line in checksum_lines} == {
        "candidates.tsv",
        "index.html",
        "report.json",
        "review_ledger_template.tsv",
    }


def test_report_is_timestamp_reproducible_and_caps_html_payload(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("SOURCE_DATE_EPOCH", "0")
    output = tmp_path / "report"
    report = build_review_report(
        candidate_gff_path=EXAMPLE / "candidate.gff3",
        pool_decisions_tsv_path=EXAMPLE / "decisions.tsv",
        pool_manifest_json_path=EXAMPLE / "candidate.gff3.manifest.json",
        output_dir=output,
        max_embedded_candidates=1,
    )
    assert report["run"]["generated_at_utc"] == "1970-01-01T00:00:00Z"
    assert report["embedding"] == {
        "maximum_requested": 1,
        "embedded": 2,
        "total": 3,
        "truncated": True,
        "conflict_closure_extra": 1,
        "order": "descending_review_priority_with_conflict_set_closure",
        "full_universe_files": [
            "report.json",
            "candidates.tsv",
            "review_ledger_template.tsv",
        ],
    }
    html = (output / "index.html").read_text(encoding="utf-8")
    assert DIGESTS[0] in html and DIGESTS[1] in html
    assert DIGESTS[2] not in html
    ledger_rows = list(
        csv.DictReader(
            (output / "review_ledger_template.tsv").open(encoding="utf-8", newline=""),
            delimiter="\t",
        )
    )
    assert len(ledger_rows) == 3
    assert tuple(ledger_rows[0]) == REVIEW_FIELDS_FOR_TEST


def test_report_conflict_workbench_contains_exact_chain_differences(tmp_path: Path) -> None:
    report = build_review_report(
        candidate_gff_path=EXAMPLE / "candidate.gff3",
        pool_decisions_tsv_path=EXAMPLE / "decisions.tsv",
        pool_manifest_json_path=EXAMPLE / "candidate.gff3.manifest.json",
        output_dir=tmp_path / "report",
    )
    assert len(report["conflict_sets"]) == 1
    conflict = report["conflict_sets"][0]
    assert conflict["member_digests"] == [DIGESTS[0], DIGESTS[1]]
    difference = conflict["pairwise_differences"][0]
    assert difference["same_cds_count"] is True
    assert difference["summed_boundary_displacement_bp"] == 100
    assert difference["phase_mismatch_segments"] == 0


def test_html_renderer_caps_a_twenty_thousand_candidate_universe(tmp_path: Path) -> None:
    report = build_review_report(
        candidate_gff_path=EXAMPLE / "candidate.gff3",
        pool_decisions_tsv_path=EXAMPLE / "decisions.tsv",
        pool_manifest_json_path=EXAMPLE / "candidate.gff3.manifest.json",
        output_dir=tmp_path / "seed_report",
    )
    prototype = report["candidates"][2]
    expanded: list[dict[str, object]] = []
    for index in range(20_000):
        row = copy.deepcopy(prototype)
        digest = hashlib.sha256(f"candidate-{index}".encode()).hexdigest()
        row["candidate_digest"] = digest
        row["gene_id"] = f"PPCONS_gene_{index}"
        row["rank_position"] = index + 1
        expanded.append(row)
    report["candidates"] = expanded
    report["conflict_sets"] = []
    report["triage"] = {
        "policy": "inspection_order_only_never_acceptance_advice",
        "counts": {"unreviewed": 20_000, "deferred": 0, "unresolved_conflict_sets": 0},
        "top_unreviewed": [],
        "deferred": [],
        "unresolved_conflict_sets": [],
    }
    report["counts"]["candidates"] = 20_000
    report["embedding"] = {
        "maximum_requested": 5000,
        "embedded": 5000,
        "total": 20_000,
        "truncated": True,
        "conflict_closure_extra": 0,
        "order": "descending_review_priority_with_conflict_set_closure",
        "full_universe_files": [
            "report.json",
            "candidates.tsv",
            "review_ledger_template.tsv",
        ],
    }
    embedded = {row["candidate_digest"] for row in expanded[:5000]}
    html = _render_html(report, embedded_digests=embedded)
    assert html.count('"candidate_digest"') == 5000
    assert len(html.encode("utf-8")) < 12_000_000


REVIEW_FIELDS_FOR_TEST = (
    "candidate_digest",
    "review_decision",
    "reviewer",
    "reviewed_at_utc",
    "rationale",
)
