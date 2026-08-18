from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import subprocess
import sys

from ploidypatch.walnut_h1_framework import BLIND_OUTPUTS, PIPELINE_ENTRIES


ROOT = Path(__file__).parents[1]
PREPARE_SPEC = importlib.util.spec_from_file_location(
    "prepare_walnut_blind_candidate_inputs_v0_8",
    ROOT / "scripts/prepare_walnut_blind_candidate_inputs_v0.8.py",
)
assert PREPARE_SPEC is not None and PREPARE_SPEC.loader is not None
PREPARE = importlib.util.module_from_spec(PREPARE_SPEC)
PREPARE_SPEC.loader.exec_module(PREPARE)


def test_blind_primary_seqids_are_table_ordered_and_json_serializable(
    tmp_path: Path,
) -> None:
    table = tmp_path / "primary.tsv"
    table.write_text(
        "seqid\tchromosome_label\nchr2\t2\nchr1\t1\n",
        encoding="utf-8",
        newline="",
    )
    fasta = tmp_path / "genome.fa"
    fasta.write_text(">chr1\nACGT\n>chr2\nTGCA\n", encoding="utf-8", newline="")
    output = tmp_path / "primary.fa"
    seqids = PREPARE.write_primary_genome(fasta, table, output)
    assert seqids == ["chr2", "chr1"]
    assert json.loads(json.dumps({"target_primary_seqids": seqids})) == {
        "target_primary_seqids": ["chr2", "chr1"]
    }
    assert output.read_text(encoding="utf-8") == ">chr1\nACGT\n>chr2\nTGCA\n"


def test_exact_species_entries_and_blind_output_paths_are_implemented() -> None:
    assert PIPELINE_ENTRIES == {
        "blind_pipeline": "scripts/run_walnut_blind_pipeline_v0.8.sh",
        "reveal_input_builder": "scripts/build_walnut_complete_control_reveal_inputs_v0.8.py",
        "evaluator": "scripts/evaluate_walnut_external_h1_v0.8.py",
    }
    assert all((ROOT / relative).is_file() for relative in PIPELINE_ENTRIES.values())
    assert set(BLIND_OUTPUTS) == {
        "raw_predictions_manifest",
        "retain_pool",
        "retain_decisions",
        "retain_manifest",
        "suppress_pool",
        "suppress_decisions",
        "suppress_manifest",
        "command_log",
    }
    assert BLIND_OUTPUTS["retain_pool"].endswith(
        "walnut_v0.8_h1/retain_distinct/blind/candidate.gff3"
    )
    assert BLIND_OUTPUTS["suppress_pool"].endswith(
        "walnut_v0.8_h1/suppress_overlap/blind/candidate.gff3"
    )


def test_blind_entry_is_h1_only_and_has_no_scoring_stage() -> None:
    text = (ROOT / PIPELINE_ENTRIES["blind_pipeline"]).read_text(encoding="utf-8")
    for required in (
        "prepare_walnut_blind_candidate_inputs_v0.8.py",
        "run_walnut_candidate_methods_v0.8.sh",
        "build_walnut_h1_candidate_pools_v0.8.py",
        "PLOIDYPATCH_NETWORK_ACCESS",
    ):
        assert required in text
    for forbidden in (
        "score_walnut",
        "topology_features",
        "candidate_labels",
        "run_walnut_evaluator_wgdi_v0.8.sh",
    ):
        assert forbidden not in text


def test_miniprot_runner_does_not_expand_a_new_local_in_the_same_command() -> None:
    text = (ROOT / "scripts/run_walnut_candidate_methods_v0.8.sh").read_text(
        encoding="utf-8"
    )
    assert "local bundle=$1 output=" not in text
    assert "local bundle=$1\n" in text
    assert "local output=$baseline_root/miniprot/$bundle\n" in text
    assert "local working=${output}.working\n" in text


def test_candidate_runner_validates_each_publishers_native_manifest() -> None:
    text = (ROOT / "scripts/run_walnut_candidate_methods_v0.8.sh").read_text(
        encoding="utf-8"
    )
    assert "if [[ $method == miniprot ]]" in text
    assert "sha256sum -c SHA256SUMS" in text
    assert "verify_published_method_output_v0.8.py" in text
    assert "final_annotation.gff" in text
    assert "lifton.gff3" in text

    pool_builder = (ROOT / "scripts/build_walnut_h1_candidate_pools_v0.8.py").read_text(
        encoding="utf-8"
    )
    assert "verify_published_method_output(tree, PUBLISHED_ROLES[method])" in pool_builder
    assert 'if method == "miniprot"' in pool_builder
    assert 'RESULT_RELATIVE = Path("results/copy_collapse/external/walnut_v0.8_h1")' in pool_builder
    assert "extract_candidate_suffix" in pool_builder
    assert "assemble_method_candidate_gff" in pool_builder
    assert "collapse_duplicate_exact_chains" in pool_builder


def test_walnut_python_entries_compile_and_shell_entries_parse() -> None:
    python_entries = [
        "scripts/prepare_walnut_external_normalized_inputs_v0.8.py",
        "scripts/prepare_walnut_blind_candidate_inputs_v0.8.py",
        "scripts/build_walnut_structure_holdout_v0.8.py",
        "scripts/prepare_walnut_evaluator_wgdi_inputs_v0.8.py",
        "scripts/prepare_walnut_target_ks_inputs_v0.8.py",
        "scripts/infer_walnut_external_pairs_v0.8.py",
        "scripts/build_walnut_h1_candidate_pools_v0.8.py",
        "scripts/build_walnut_complete_control_reveal_inputs_v0.8.py",
        "scripts/evaluate_walnut_external_h1_v0.8.py",
        "scripts/verify_published_method_output_v0.8.py",
    ]
    subprocess.run(
        [sys.executable, "-m", "py_compile", *python_entries], cwd=ROOT, check=True
    )
    if sys.platform != "win32":
        subprocess.run(
            [
                "bash",
                "-n",
                "scripts/run_walnut_blind_pipeline_v0.8.sh",
                "scripts/run_walnut_candidate_methods_v0.8.sh",
                "scripts/run_walnut_evaluator_wgdi_v0.8.sh",
            ],
            cwd=ROOT,
            check=True,
        )


def test_complete_control_builder_and_evaluator_are_strict_h1_only() -> None:
    builder = (ROOT / PIPELINE_ENTRIES["reveal_input_builder"]).read_text(encoding="utf-8")
    evaluator = (ROOT / PIPELINE_ENTRIES["evaluator"]).read_text(encoding="utf-8")
    for text in (builder, evaluator):
        assert "PLOIDYPATCH_COMPOSITE_MODEL_FREEZE" not in text
        assert "PLOIDYPATCH_RANKER_FREEZE" not in text
    assert "control_candidate_gff_path" in builder
    assert "score_collateral_gate" in builder
    assert "evaluate_h1_scores" in evaluator
    assert "average_precision" not in evaluator


def test_pool_custody_preserves_each_arms_native_schema() -> None:
    core = (ROOT / "src/ploidypatch/walnut_h1.py").read_text(encoding="utf-8")
    custody = (ROOT / "scripts/finalize_walnut_blind_custody_v0.8.py").read_text(
        encoding="utf-8"
    )
    assert 'RETAIN_POOL_SCHEMA = "ploidypatch.method_candidate_pool.v2"' in core
    assert 'SUPPRESS_POOL_SCHEMA = "ploidypatch.method_consensus.v1"' in core
    assert "from ploidypatch.walnut_h1 import POOL_SCHEMAS" in custody
    assert "manifest.get(\"schema_version\") != POOL_SCHEMAS[arm]" in custody


def test_evaluator_truth_tree_is_never_claimed_as_blind_mounted() -> None:
    finalizer = (
        ROOT / "scripts/finalize_walnut_evaluator_only_inputs_v0.8.py"
    ).read_text(encoding="utf-8")
    assert "evaluator_only_never_mounted_to_blind_runner" in finalizer
    assert '"truth_reveal_requires_sealed_blind_custody": True' in finalizer
    assert "evaluator_only_after_blind_custody" not in finalizer


def test_evaluator_wgdi_records_global_parallel_binary_provenance() -> None:
    runner = (
        ROOT / "scripts/run_walnut_evaluator_wgdi_v0.8.sh"
    ).read_text(encoding="utf-8")
    assert 'tool_manifest.tsv' in runner
    assert 'sha256sum "$parallel"' in runner
    assert '"$parallel" --version' in runner


def test_evaluator_outgroups_use_primary_header_aware_protein_mapping() -> None:
    preparer = (
        ROOT / "scripts/prepare_walnut_evaluator_wgdi_inputs_v0.8.py"
    ).read_text(encoding="utf-8")
    assert "write_primary_candidate_proteins" in preparer
    assert 'if label != "target"' in preparer
    assert "primary.protein_whitelist.tsv" in preparer
    assert "primary_protein_manifest_sha256" in preparer
