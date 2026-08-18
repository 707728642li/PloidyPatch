from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).parents[1]
SCRIPTS = {
    name: (ROOT / "scripts" / name).read_text(encoding="utf-8")
    for name in (
        "prepare_populus_external_normalized_inputs_v0.4.sh",
        "prepare_populus_evaluator_wgdi_inputs_v0.4.sh",
        "run_populus_evaluator_wgdi_v0.4.sh",
        "infer_populus_external_pairs_v0.4.sh",
        "run_populus_copy_collapse_benchmark_v0.4.sh",
        "build_populus_complete_control_reveal_inputs_v0.4.sh",
    )
}


def test_canonical_normalized_and_wgdi_paths_are_consistent() -> None:
    normalized = "data/derived/external_inputs/populus_v0.4"
    wgdi = "data/derived/external_evaluator/populus_v0.4_wgdi_inputs"
    assert normalized in SCRIPTS["prepare_populus_external_normalized_inputs_v0.4.sh"]
    for name in (
        "prepare_populus_evaluator_wgdi_inputs_v0.4.sh",
        "infer_populus_external_pairs_v0.4.sh",
        "run_populus_copy_collapse_benchmark_v0.4.sh",
        "build_populus_complete_control_reveal_inputs_v0.4.sh",
    ):
        assert normalized in SCRIPTS[name]
    for name in (
        "prepare_populus_evaluator_wgdi_inputs_v0.4.sh",
        "run_populus_evaluator_wgdi_v0.4.sh",
        "infer_populus_external_pairs_v0.4.sh",
    ):
        assert wgdi in SCRIPTS[name]
        assert "external_evaluator/populus/v0.4_wgdi_inputs" not in SCRIPTS[name]


def test_normalization_preserves_role_separated_staging_as_plain_hashed_files() -> None:
    script = SCRIPTS["prepare_populus_external_normalized_inputs_v0.4.sh"]
    assert "shutil.copyfile(incoming, outgoing)" in script
    assert "incoming.is_symlink()" in script and "outgoing.is_symlink()" in script
    assert 'for name in ("role_manifest.tsv", "role_contract.json")' in script
    assert "staging_SHA256SUMS" in script
    assert "BLIND_SHA256SUMS" in script
    assert "EVALUATOR_SHA256SUMS" in script
    for role_root in ("shared_target", "candidate_only", "evaluator_only"):
        assert role_root in script


def test_evaluator_stages_do_not_checksum_candidate_reference_subtree() -> None:
    for name in (
        "prepare_populus_evaluator_wgdi_inputs_v0.4.sh",
        "infer_populus_external_pairs_v0.4.sh",
        "run_populus_copy_collapse_benchmark_v0.4.sh",
        "build_populus_complete_control_reveal_inputs_v0.4.sh",
    ):
        assert "EVALUATOR_SHA256SUMS" in SCRIPTS[name]


def test_evaluator_scripts_lock_self_to_execution_manifest() -> None:
    for name, script in SCRIPTS.items():
        assert f"self_relative=scripts/{name}" in script
        assert "code_root=$execution_root/source" in script
        assert 'export PYTHONPATH="$code_root/src${PYTHONPATH:+:$PYTHONPATH}"' in script
        assert "environment_bindings=$execution_root/environment_bindings.tsv" in script
        assert "execution_root/implementation_manifest.tsv" in script
        assert "stat -Lc %s" in script
        assert "sha256sum \"$code_root/$relative\"" in script
        assert ".working" in script
        assert "refusing to overwrite" in script


def test_wgdi_and_pair_contract_freezes_populus_parameters() -> None:
    wgdi = SCRIPTS["run_populus_evaluator_wgdi_v0.4.sh"]
    pairs = SCRIPTS["infer_populus_external_pairs_v0.4.sh"]
    for prefix in ("ptr", "mes", "rco"):
        assert prefix in wgdi
    assert "--evalue 1e-5 --max-target-seqs 20 --more-sensitive" in wgdi
    assert "score = 100" in wgdi and "repeat_number = 20" in wgdi
    assert "--min-support-group-count 2 --min-block-pairs 20" in pairs
    assert "intersect-copy-pair-evidence" in pairs
    assert "candidate_reference_access\\tfalse" in pairs


def test_benchmark_preserves_evaluability_and_nonleaking_blind_manifest() -> None:
    script = SCRIPTS["run_populus_copy_collapse_benchmark_v0.4.sh"]
    assert "seed=20260930" in script and "maximum_count=800" in script
    assert "not_evaluable_without_rule_relaxation" in script
    assert "invalid_run" in script
    assert "blind_noop_exact_recovery" in script
    assert "complete_oracle_exact_recovery" in script
    assert "blind/manifest.json\" \"$working_root/evaluator/perturbation_manifest.json" in script
    safe_start = script.index('"schema_version": "ploidypatch.blind_benchmark_input.v0.4"')
    safe_end = script.index("with output.open", safe_start)
    safe_manifest = script[safe_start:safe_end]
    assert "hidden_truth" not in safe_manifest
    assert "complete_target_annotation_access\": False" in safe_manifest


def test_reveal_builder_reads_only_isolated_blind_project_and_requires_barrier() -> None:
    script = SCRIPTS["build_populus_complete_control_reveal_inputs_v0.4.sh"]
    assert "blind_project=$blind_run_root/project" in script
    assert "upstream=$blind_project/results/baselines/populus_v0.4" in script
    assert "method_root=$blind_project/results/copy_collapse/external/populus_v0.4_method_trio" in script
    assert "ranking_root=$blind_project/results/copy_collapse/external/populus_v0.4_blind_rankings" in script
    assert "PLOIDYPATCH_BLIND_RUN_ROOT" in script
    assert "PLOIDYPATCH_REVEAL_AUTHORIZATION" in script
    assert "PLOIDYPATCH_EXECUTION_FREEZE_OVERRIDE" in script
    assert "ploidypatch.populus_reveal_authorization.v0.4" in script
    assert 'authorization.get("truth_opened") is not False' in script
    barrier = script.index('"$python_bin" - "$authorization"')
    complete_access = script.index("# Custody has now passed")
    assert barrier < complete_access


def test_reveal_builder_emits_canonical_evaluator_contract() -> None:
    script = SCRIPTS["build_populus_complete_control_reveal_inputs_v0.4.sh"]
    assert "results/evaluator/populus/v0.4/reveal_inputs" in script
    assert '"schema_version": "ploidypatch.populus_reveal_inputs.v0.4"' in script
    assert '"generated_after_blind_freeze": True' in script
    assert '"evaluator_only": True' in script
    for key in (
        "blind_run_SHA256SUMS_sha256",
        "custody_manifest_sha256",
        "blind_scores_sha256",
        "blind_score_manifest_sha256",
        "pool_decisions_sha256",
        "pool_manifest_sha256",
    ):
        assert f'"{key}"' in script
    for relative in (
        "labels/candidate_labels.tsv",
        "scores/consensus/primary_union.json",
        "scores/consensus/legacy_union.json",
        "scores/methods/miniprot.json",
        "scores/methods/gemoma.json",
        "scores/methods/lifton.json",
        "scores/consensus/support2.json",
        "scores/consensus/support3.json",
    ):
        assert relative in script
    assert '"formal_evaluator_environment": "ploidypatch-model"' in script
    assert '"status": status' in script


def test_reveal_builder_has_fail_closed_three_state_dispatch() -> None:
    script = SCRIPTS["build_populus_complete_control_reveal_inputs_v0.4.sh"]
    assert "if [[ $formal_outcome == invalid_run ]]" in script
    assert "write_manifest invalid_run" in script
    assert "elif [[ $formal_outcome == not_evaluable_without_rule_relaxation ]]" in script
    assert "write_manifest not_evaluable_without_rule_relaxation" in script
    assert 'label_status=$("$python_bin"' in script
    assert 'print("ready_for_evaluation" if' in script
    assert 'write_manifest "$label_status"' in script
    disk = script.index('du -sb "$working_root" > "$working_root/disk_bytes.txt"')
    checksum = script.index("xargs -0 sha256sum > SHA256SUMS", disk)
    assert disk < checksum
