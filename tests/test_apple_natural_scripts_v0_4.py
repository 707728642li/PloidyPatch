from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def _text(name: str) -> str:
    return (ROOT / "scripts" / name).read_text(encoding="utf-8")


def test_natural_discovery_uses_complete_control_without_rna() -> None:
    self_wgd = _text("run_apple_natural_self_wgd_v0.4.sh")
    score = _text("score_apple_natural_candidates_v0.4.sh")
    assert "consensus/primary_union/complete_control" in self_wgd
    assert "primary_chromosomes.gff3" in self_wgd
    assert "RNA_access\\tfalse" in self_wgd
    assert "candidate_and_rank_freeze_precedes_validation_access\\ttrue" in self_wgd
    assert "methods/miniprot/complete_control/decisions.tsv" in score
    assert "methods/gemoma/complete_control/decisions.tsv" in score
    assert "methods/lifton/complete_control/decisions.tsv" in score
    assert "apply-conflict-winner-guard" in score
    assert "v04_primary_rank_score" in score
    assert "automatic_approval\\tfalse" in score
    for relative in (
        "scripts/build_wgdi_source_alias_gff.py",
        "src/ploidypatch/self_wgd_pairs.py",
        "src/ploidypatch/wgd_candidate_select.py",
        "src/ploidypatch/synteny_io.py",
        "src/ploidypatch/copy_features.py",
        "src/ploidypatch/homeolog_topology.py",
        "src/ploidypatch/support_ranker.py",
        "src/ploidypatch/conflict_guard.py",
    ):
        assert (ROOT / relative).is_file(), relative


def test_ont_download_is_barriered_and_checksum_locked() -> None:
    script = _text("download_apple_golden_delicious_ont_v0.1.sh")
    assert "candidate_and_rank_freeze_precedes_validation_access" in script
    assert "RNA_access" in script
    assert "CRA021523" in script
    assert script.count("https://download.cncb.ac.cn/gsa4/CRA021523/CRR") == 7
    for accession in range(1429911, 1429918):
        assert f"CRR{accession}" in script
    assert "md5sum" in script
    assert "official_md5" in script
    assert "file contract differs from official GSA MD5" in script
    assert "gzip -t" in script
    assert "sha256sum" in script
    assert "download_script_sha256" in script
    assert "axel -n 8" in script


def test_apple_ont_resume_is_sequential_and_preserves_contract() -> None:
    script = _text("resume_apple_golden_delicious_ont_v0.1.sh")
    assert "parallel_files\\t1" in script
    assert "source_changed\\tfalse" in script
    assert "checksums_changed\\tfalse" in script
    assert "candidate_coordinates_or_ranks_modified\\tfalse" in script
    assert "${output}.st" in script
    assert "attempt -lt 20" in script
    assert "transport_attempts.tsv" in script
    assert "resume_invocations.tsv" in script
    assert "md5sum" in script
    assert "gzip -t" in script


def test_apple_final_accession_parallel_resume_is_audited() -> None:
    script = _text("parallel_resume_apple_ont_final_accession_v0.1.sh")
    assert "accession=CRR1429917" in script
    assert "parallel_files_total\\t2" in script
    assert "source_changed\\tfalse" in script
    assert "checksums_changed\\tfalse" in script
    assert "candidate_coordinates_or_ranks_modified\\tfalse" in script
    assert "attempt -lt 20" in script
    assert "md5sum" in script
    assert "gzip -t" in script


def test_csv_limit_failure_is_preserved_before_parser_only_patch() -> None:
    freeze = _text("freeze_apple_ont_csv_limit_failed_attempt_v0.1.sh")
    resume = _text("resume_apple_ont_post_alignment_patch1_v0.1.sh")
    assert "field larger than field limit (131072)" in freeze
    assert "labels_accessed\\ttrue" in freeze
    assert "scientific_thresholds_changed\\tfalse" in freeze
    assert "ont_raw_v0.1_failed_csv_field_limit" in freeze
    assert "ont_raw_v0.1_patch1" in resume
    assert "labels_seen_before_patch\\ttrue" in resume
    assert "alignment_recomputed\\tfalse" in resume
    assert "strict_chain_evidence_recomputed\\tfalse" in resume
    assert "--minimum-full-length-read-support 2" in resume
    assert "--comparator-estimator baseline --primary-estimator v04_guard" in resume
    assert "--replicates 20000 --seed 20261004" in resume


def test_discovery_and_download_are_nonoverwriting_atomic_freezes() -> None:
    for name in (
        "run_apple_natural_self_wgd_v0.4.sh",
        "score_apple_natural_candidates_v0.4.sh",
        "download_apple_golden_delicious_ont_v0.1.sh",
        "resume_apple_golden_delicious_ont_v0.1.sh",
        "download_apple_gddh13_te_annotation_v0.1.sh",
        "run_apple_golden_delicious_ont_validation_v0.1.sh",
        "resume_apple_ont_transcript_strand_patch2_v0.1.sh",
    ):
        script = _text(name)
        assert ".working" in script
        assert "refusing to overwrite" in script
        assert "chmod -R a-w" in script
        assert "mv \"$working_root\" \"$result_root\"" in script


def test_te_annotation_is_version_and_byte_locked() -> None:
    script = _text("download_apple_gddh13_te_annotation_v0.1.sh")
    assert "GDDH13_1-1_TE.gff3.bz2" in script
    assert "expected_bytes=34251728" in script
    assert "bzip2 -t" in script
    assert "source_etag" in script
    assert "candidate_and_rank_freeze_precedes_validation_access" in script


def test_raw_ont_validation_uses_exact_chain_and_biological_audits() -> None:
    script = _text("run_apple_golden_delicious_ont_validation_v0.1.sh")
    assert "parallel_runs\\t7" in script
    assert "threads_per_run\\t16" in script
    assert "shared_minimap2_index\\ttrue" in script
    assert "-x splice -k14 -d \"$genome_index\"" in script
    assert "-c --cs=long \"$genome_index\"" in script
    assert "--secondary=yes -N 10" in script
    assert "--minimum-query-coverage 0.85" in script
    assert "--minimum-identity 0.90" in script
    assert "--minimum-mapq 20" in script
    assert "--maximum-secondary-score-fraction 0.95" in script
    assert "--minimum-candidate-cds-coverage 0.90" in script
    assert "--minimum-full-length-read-support 2" in script
    assert "candidate_query_prefilter\\tlossless_all_alignments" in script
    assert "alignment_strand_source\\tminimap2_ts" in script
    assert "reference_transcript_strand=paf_query_target_strand*" in script
    assert "filter-candidate-query-paf" in script
    assert script.count("--alignment-strand-source minimap2_ts") == 2
    assert "candidate_query_universe.paf" in script
    assert "query_filter_manifest.json" in script
    assert "--comparator-estimator baseline --primary-estimator v04_guard" in script
    assert "--repeat-gff" in script
    assert "candidate_cds_to_genome.paf" in script
    assert "--replicates 20000 --seed 20261004" in script
    assert "automatic_annotation_patch\\tfalse" in script


def test_transcript_strand_patch2_is_nonoverwriting_and_reuses_frozen_paf() -> None:
    script = _text("resume_apple_ont_transcript_strand_patch2_v0.1.sh")
    assert "raw_reads_remapped\\tfalse" in script
    assert "frozen_raw_PAF_reused\\ttrue" in script
    assert "alignment_strand_source\\tminimap2_ts" in script
    assert script.count("--alignment-strand-source minimap2_ts") == 2
    assert "patch1_status\\tdescriptive_invalid_transcript_strand_interpretation" in script
    assert "sha256sum -c SHA256SUMS" in script
    assert "patch1_vs_patch2.json" in script
