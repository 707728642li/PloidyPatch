from __future__ import annotations

import csv
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys

import pytest

from ploidypatch.artifact_manifest import (
    read_sha256sums,
    verify_sha256sums,
    write_sha256sums,
)


ROOT = Path(__file__).parents[1]
NAMES = (
    "prepare_actinidia_external_normalized_inputs_v0.5.sh",
    "prepare_actinidia_evaluator_wgdi_inputs_v0.5.sh",
    "run_actinidia_evaluator_wgdi_v0.5.sh",
    "infer_actinidia_external_pairs_v0.5.sh",
    "run_actinidia_copy_collapse_benchmark_v0.5.sh",
    "build_actinidia_complete_control_reveal_inputs_v0.5.sh",
)
SCRIPTS = {
    name: (ROOT / "scripts" / name).read_text(encoding="utf-8") for name in NAMES
}


def test_actinidia_evaluator_paths_and_frozen_implementation_are_consistent() -> None:
    normalized = "data/derived/external_inputs/actinidia/v0.5"
    wgdi_inputs = "data/derived/external_evaluator/actinidia_v0.5_wgdi_inputs"
    for name in NAMES:
        script = SCRIPTS[name]
        assert "data/derived/external_inputs/actinidia_v0.5" not in script
        assert f"self_relative=scripts/{name}" in script
        assert "code_root=$execution_root/source" in script
        assert 'export PYTHONPATH="$code_root/src${PYTHONPATH:+:$PYTHONPATH}"' in script
        assert "execution_root/implementation_manifest.tsv" in script
        assert "stat -Lc %s" in script
        assert 'sha256sum "$code_root/$relative"' in script
        assert ".working" in script
        assert "refusing to overwrite" in script
    for name in (
        "prepare_actinidia_external_normalized_inputs_v0.5.sh",
        "prepare_actinidia_evaluator_wgdi_inputs_v0.5.sh",
        "infer_actinidia_external_pairs_v0.5.sh",
        "run_actinidia_copy_collapse_benchmark_v0.5.sh",
        "build_actinidia_complete_control_reveal_inputs_v0.5.sh",
    ):
        assert normalized in SCRIPTS[name]
    for name in (
        "prepare_actinidia_evaluator_wgdi_inputs_v0.5.sh",
        "run_actinidia_evaluator_wgdi_v0.5.sh",
        "infer_actinidia_external_pairs_v0.5.sh",
    ):
        assert wgdi_inputs in SCRIPTS[name]


def test_normalization_verifies_all_15_role_artifacts_and_keeps_candidate_firewall() -> None:
    script = SCRIPTS["prepare_actinidia_external_normalized_inputs_v0.5.sh"]
    assert "ploidypatch.external_holdout_input_stage.v0.5" in script
    assert 'protocol_root / "role_manifest.tsv"' in script
    assert "len(frozen) != 15" in script
    assert 'for field in ("release", "bundle_id", "wgdi_prefix", "bytes", "sha256", "staged_relative_path")' in script
    assert "shutil.copyfile(incoming, outgoing)" in script
    assert "incoming.is_symlink()" in script and "outgoing.is_symlink()" in script
    assert "BLIND_SHA256SUMS" in script and "EVALUATOR_SHA256SUMS" in script
    for role_root in ("shared_target", "candidate_only", "evaluator_only"):
        assert role_root in script
    assert "candidate_eriantha" not in script
    assert "candidate_rufa" not in script
    assert "normalized/evaluator_rhododendron" in script
    assert "normalized/evaluator_diospyros" in script


def test_red5_peptides_are_exactly_whitelisted_by_chromosome_gff() -> None:
    normalize = SCRIPTS["prepare_actinidia_external_normalized_inputs_v0.5.sh"]
    prepare = SCRIPTS["prepare_actinidia_evaluator_wgdi_inputs_v0.5.sh"]
    for token in (
        "chromosome_protein_id_whitelist.tsv",
        "exact_GFF_CDS_protein_id_to_FASTA_first_token",
        'attrs.get("protein_id")',
        "len(protein_ids) != 33021",
        '"all_provider_proteins": 33115',
        '"excluded_nonchromosome_proteins": 94',
        "seen != protein_ids",
    ):
        assert token in normalize
    assert "if [[ $name == target_red5 ]]" in prepare
    assert "chromosome_whitelist.pep.fa" in prepare
    assert "gffread" in prepare  # evaluator proteins remain genome/GFF derived
    for prefix, genes, chromosomes in (
        ("red5", 32950, 29),
        ("rhs", 30419, 13),
        ("dol", 29268, 15),
    ):
        assert f"{prefix}) expected_genes={genes}; expected_chromosomes={chromosomes}" in prepare


def test_wgdi_contract_uses_only_red5_and_the_two_frozen_evaluators() -> None:
    script = SCRIPTS["run_actinidia_evaluator_wgdi_v0.5.sh"]
    assert "red5 rhs dol" in script
    assert "red5_self" in script
    assert "red5_vs_rhs" in script and "red5_vs_dol" in script
    assert "Actinidia_chinensis_Red5_PS1_1.69.0_EnsemblPlants45" in script
    assert "Actinidia_specific_Ad_alpha" in script
    assert "--evalue 1e-5 --max-target-seqs 20 --more-sensitive" in script
    assert "score = 100" in script and "repeat_number = 20" in script
    assert "position = order" in script and "process = 64" in script
    assert "candidate_reference_access\\tfalse" in script


def test_pair_wrapper_enforces_each_group_before_exact_final_intersection() -> None:
    script = SCRIPTS["infer_actinidia_external_pairs_v0.5.sh"]
    assert "red5_self.tsv" in script
    assert "red5_vs_rhs.tsv" in script and "red5_vs_dol.tsv" in script
    assert script.count("infer-outgroup-duplicated-pairs") == 2
    assert script.count("--min-support-group-count 1 --min-block-pairs 20") == 2
    assert "evaluator_groups/rhododendron/decisions.tsv" in script
    assert "evaluator_groups/diospyros/decisions.tsv" in script
    assert "pair_consistency_audit.tsv" in script
    assert "discordant_target_members" in script
    assert "alternate exact-1:2 pair in either group" in script
    assert "intersect-copy-pair-evidence" in script
    assert "Actinidia_specific_Ad_alpha_self_and_two_evaluator_groups" in script


def _pair_audit_code() -> str:
    script = SCRIPTS["infer_actinidia_external_pairs_v0.5.sh"]
    marker = '"$working_root/pair_consistency_audit.json" <<\'PY\'\n'
    start = script.index(marker) + len(marker)
    end = script.index("\nPY\n", start)
    return script[start:end]


def _write_decisions(path: Path, rows: list[tuple[str, str, str, str]]) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("gene_id_a", "gene_id_b", "status", "reason"))
        writer.writerows(rows)


def _run_pair_audit(tmp_path: Path, *, unsafe_accept: bool = False) -> dict[str, object]:
    rh = tmp_path / "rhododendron.tsv"
    dol = tmp_path / "diospyros.tsv"
    accepted = tmp_path / "two_groups.tsv"
    audit_tsv = tmp_path / "audit.tsv"
    audit_json = tmp_path / "audit.json"
    _write_decisions(
        rh,
        [
            ("A", "B", "accepted" if unsafe_accept else "rejected", "nonreciprocal_multiple_partners"),
            ("A", "C", "rejected", "nonreciprocal_multiple_partners"),
        ],
    )
    _write_decisions(
        dol,
        [("A", "B", "accepted", "independent_one_to_two_outgroup_support")],
    )
    with accepted.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("gene_id_a", "gene_id_b"))
    old_argv = sys.argv
    sys.argv = [
        "pair_consistency_audit",
        str(rh),
        str(dol),
        str(accepted),
        str(audit_tsv),
        str(audit_json),
    ]
    try:
        exec(compile(_pair_audit_code(), "pair_consistency_audit", "exec"), {})
    finally:
        sys.argv = old_argv
    return json.loads(audit_json.read_text(encoding="utf-8"))


def test_two_serial_wgd_toy_rejects_one_group_alternative_pair(tmp_path: Path) -> None:
    report = _run_pair_audit(tmp_path)
    assert report["accepted_exact_two_group_intersection_pairs"] == 0
    by_group = {row["evaluator_group"]: row for row in report["evaluator_groups"]}
    assert by_group["rhododendron"]["qualifying_exact_1_to_2_pairs"] == 2
    assert by_group["rhododendron"]["discordant_target_members"] == 1
    assert by_group["rhododendron"]["accepted_pair_consistent_pairs"] == 0


def test_pair_audit_fails_closed_if_core_erroneously_accepts_incident_pair(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="accepted a pair incident"):
        _run_pair_audit(tmp_path, unsafe_accept=True)


def test_benchmark_binds_red5_gates_without_truth_leakage() -> None:
    script = SCRIPTS["run_actinidia_copy_collapse_benchmark_v0.5.sh"]
    assert "seed=20261010" in script and "maximum_count=800" in script
    assert "red5_ps1_1.69.0" in script
    assert "minimum_target_chromosomes) == 22" in script
    assert "Actinidia_specific_Ad_alpha" in script
    assert "not_evaluable_without_rule_relaxation" in script and "invalid_run" in script
    assert "blind_noop_exact_recovery" in script
    assert "complete_oracle_exact_recovery" in script
    safe_start = script.index('"schema_version": "ploidypatch.blind_benchmark_input.v0.5"')
    safe_end = script.index("with output.open", safe_start)
    safe_manifest = script[safe_start:safe_end]
    assert "hidden_truth" not in safe_manifest
    assert '"complete_target_annotation_access": False' in safe_manifest


def test_reveal_builder_uses_generic_freeze_environment_and_three_states() -> None:
    script = SCRIPTS["build_actinidia_complete_control_reveal_inputs_v0.5.sh"]
    for name in (
        "PLOIDYPATCH_HOLDOUT_CONTRACT",
        "PLOIDYPATCH_PROTOCOL_FREEZE",
        "PLOIDYPATCH_EXECUTION_FREEZE",
        "PLOIDYPATCH_COMPOSITE_MODEL_FREEZE",
        "PLOIDYPATCH_BLIND_RUN_ROOT",
        "PLOIDYPATCH_CUSTODY_MANIFEST",
        "PLOIDYPATCH_REVEAL_AUTHORIZATION",
        "PLOIDYPATCH_EVALUATOR_ONLY_ROOT",
        "PLOIDYPATCH_REVEAL_INPUTS_OUTPUT",
    ):
        assert name in script
    assert "ploidypatch.external_holdout_blind_custody.v0.5" in script
    assert "ploidypatch.external_holdout_reveal_authorization.v0.5" in script
    assert '"holdout_id": "actinidia_red5_v0.5"' in script
    assert '"policy_id": "ploidypatch_actinidia_external_validation_v0.5"' in script
    assert '"formal_status": status' in script
    assert '"evaluability": {' in script and '"relative_path": "evaluability.json"' in script
    assert "reveal_input_manifest.json" in script and "status.json" in script
    assert "ready_for_evaluation" in script
    assert "not_evaluable_without_rule_relaxation" in script
    assert "invalid_run" in script
    assert "from ploidypatch.artifact_manifest import verify_sha256sums, write_sha256sums" in script
    assert "write_sha256sums(root)" in script
    assert "verify_sha256sums(root, ignore_checksum_file=True)" in script


def test_reveal_builder_checksum_writer_is_strict_and_canonical(tmp_path: Path) -> None:
    root = tmp_path / "reveal_inputs"
    (root / "labels").mkdir(parents=True)
    (root / "scores" / "consensus").mkdir(parents=True)
    (root / "status.json").write_text(
        '{"formal_status":"ready_for_evaluation"}\n', encoding="utf-8"
    )
    (root / "labels" / "candidate_labels.tsv").write_text(
        "candidate_digest\tlabel_exact_cds\na\t1\n", encoding="utf-8"
    )
    (root / "scores" / "consensus" / "primary_union.json").write_text(
        "{}\n", encoding="utf-8"
    )

    write_sha256sums(root)
    observed = verify_sha256sums(root, ignore_checksum_file=True)
    paths = [
        entry.relative_path.as_posix()
        for entry in read_sha256sums(root / "SHA256SUMS")
    ]

    assert set(observed) == {
        "labels/candidate_labels.tsv",
        "scores/consensus/primary_union.json",
        "status.json",
    }
    assert paths == sorted(paths)
    assert all(not path.startswith("./") and "/./" not in path for path in paths)


def test_no_forbidden_populus_residue_in_actinidia_scripts() -> None:
    forbidden = re.compile(
        r"Populus|populus|ptr_v4|\bptr\b|Salix|Manihot|Ricinus|"
        r"20260930|2026100[1-3]|Actinidia_trichocarpa|ptr_self|ptr_vs|two_outgroup"
    )
    for name, script in SCRIPTS.items():
        assert forbidden.search(script) is None, name


@pytest.mark.skipif(shutil.which("bash") is None, reason="bash is unavailable")
def test_all_actinidia_evaluator_shell_scripts_parse() -> None:
    for name in NAMES:
        subprocess.run(
            ["bash", "-n"],
            input=SCRIPTS[name].replace("\r\n", "\n").encode("utf-8"),
            check=True,
        )
