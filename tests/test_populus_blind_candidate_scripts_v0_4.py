from pathlib import Path


ROOT = Path(__file__).parents[1]
SCRIPTS = {
    name: (ROOT / "scripts" / name).read_text(encoding="utf-8")
    for name in (
        "run_populus_miniprot_upstream_v0.4.sh",
        "build_populus_method_trio_candidate_pools_v0.4.sh",
        "run_populus_blind_union_self_wgd_v0.4.sh",
        "score_populus_candidates_blind_v0.4.sh",
    )
}


def test_all_candidate_stages_enforce_canonical_blind_custody() -> None:
    for name, script in SCRIPTS.items():
        assert "data/derived/external_inputs/populus_v0.4" in script, name
        assert "PLOIDYPATCH_BLIND_RUNNER" in script, name
        assert '"$input_root/evaluator_only"' in script, name
        assert "/nas_data/NFS" not in script, name
        assert "populus_external_v0.4_execution" in script, name
        assert "verify_implementation" in script, name
        assert "blind_runner.mountinfo" in script, name


def test_miniprot_uses_only_target_genome_and_salix_candidate_roles() -> None:
    script = SCRIPTS["run_populus_miniprot_upstream_v0.4.sh"]
    assert "shared_target/Populus_trichocarpa/Ptrichocarpa_533_v4.0.fa.gz" in script
    assert "candidate_only/Salix_purpurea" in script
    assert "candidate_only/Salix_suchowensis" in script
    assert "Ptrichocarpa_533_v4.1.gene.gff3" not in script
    assert "Ptrichocarpa_533_v4.1.protein.fa" not in script
    assert "candidate_reference" in script
    assert "--protein \"salix_purpurea=" in script
    assert "--protein \"salix_suchowensis=" in script
    assert "within_method_reference_vote_count\\t1" in script
    assert "from ploidypatch.io import iter_fasta" in script
    assert "range(0, len(encoded), 60)" in script
    assert 'export PYTHONPATH="$code_root/src${PYTHONPATH:+:$PYTHONPATH}"' in script
    assert "synthesize_missing_transcript_exons.py" in script
    assert "coordinate_or_cds_changes" in script
    assert "unresolved_transcripts" in script


def test_method_trio_has_one_vote_per_family_and_both_pool_policies() -> None:
    script = SCRIPTS["build_populus_method_trio_candidate_pools_v0.4.sh"]
    for method in ("miniprot", "gemoma", "lifton"):
        assert f'--candidate "{method}=' in script
    assert script.count('--candidate "miniprot=') == 1
    assert script.count('--candidate "gemoma=') == 1
    assert script.count('--candidate "lifton=') == 1
    assert "primary_union) support=1; redundancy=retain_distinct_chains" in script
    assert "legacy_union) support=1; redundancy=suppress_overlapping" in script
    assert "--min-identity 0.5" in script
    assert "--min-query-coverage 0.5" in script
    assert "--max-existing-cds-overlap 0.2" in script
    assert "--max-redundancy-overlap 0.5" in script


def test_blind_self_wgd_recomputes_only_from_primary_union() -> None:
    script = SCRIPTS["run_populus_blind_union_self_wgd_v0.4.sh"]
    assert "consensus/primary_union/blind/candidate.gff3" in script
    assert "salicoid_wgd_blind_recomputed" in script
    assert "--min-block-pairs 20" in script
    assert "--evalue 1e-5" in script
    assert "--max-target-seqs 20" in script
    assert "repeat_number = 20" in script
    assert "candidate_candidate_pair_policy\\treject_as_circular" in script


def test_v03_production_score_precedes_v04_production_guard() -> None:
    script = SCRIPTS["score_populus_candidates_blind_v0.4.sh"]
    scorer = script.index("score-support-conditioned-candidates")
    guard = script.index("apply-conflict-winner-guard")
    assert scorer < guard
    assert "results/models/ploidypatch_ranker_v0.4" in script
    assert "model_v0.3.json" in script
    assert '--pool-manifest "$pool_manifest"' in script
    assert '--output-tsv "$working_root/scores/v04.tsv"' in script
    assert "baseline_mapping_sha256" in script
    assert "v04_guard_mapping_sha256" in script
    assert "winner_mismatch_count" in script
    assert "v04_automatic_approval" in script
