#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "usage: $0 ISOLATED_PROJECT_ROOT" >&2
    exit 2
fi
project_root=$(realpath "$1")
code_root=$project_root/code
python_bin=$project_root/envs/ploidypatch-dev/bin/python
input_root=${PLOIDYPATCH_STAGED_INPUT_ROOT:?PLOIDYPATCH_STAGED_INPUT_ROOT is required}
blind_benchmark_root=${PLOIDYPATCH_BLIND_BENCHMARK_ROOT:?PLOIDYPATCH_BLIND_BENCHMARK_ROOT is required}
protocol_root=${PLOIDYPATCH_PROTOCOL_FREEZE:?PLOIDYPATCH_PROTOCOL_FREEZE is required}
execution_root=${PLOIDYPATCH_EXECUTION_FREEZE:?PLOIDYPATCH_EXECUTION_FREEZE is required}
contract_path=${PLOIDYPATCH_HOLDOUT_CONTRACT:?PLOIDYPATCH_HOLDOUT_CONTRACT is required}
composite_root=${PLOIDYPATCH_COMPOSITE_MODEL_FREEZE:?PLOIDYPATCH_COMPOSITE_MODEL_FREEZE is required}
blind_output_root=${PLOIDYPATCH_BLIND_OUTPUT_ROOT:?PLOIDYPATCH_BLIND_OUTPUT_ROOT is required}
context_verifier=$code_root/scripts/verify_external_holdout_blind_context_v0.5.py
command_log=$project_root/pipeline_commands.tsv
command_log_working=${command_log}.working
baseline_root=$project_root/results/baselines/actinidia_v0.5
normalized=$baseline_root/miniprot/normalized
target=$normalized/target/primary_chromosomes.genome.fa

[[ ${PLOIDYPATCH_BLIND_RUNNER:-} == 1 ]] || {
    echo "Actinidia pipeline requires the frozen blind runner" >&2
    exit 1
}
[[ ${PLOIDYPATCH_NETWORK_ACCESS:-} == none ]] || {
    echo "Actinidia pipeline requires a network-disabled namespace" >&2
    exit 1
}

# A second check inside the namespace makes accidental mount broadening fail
# before an upstream program starts.
for forbidden in \
    "$input_root/evaluator_only" "$blind_benchmark_root/evaluator" \
    "$blind_benchmark_root/truth" "$blind_benchmark_root/complete" \
    "$project_root/truth" "$project_root/labels" /nas_data; do
    [[ ! -e $forbidden ]] || {
        echo "forbidden evaluator/NAS path visible in blind namespace: $forbidden" >&2
        exit 1
    }
done
if [[ -r /proc/self/mountinfo ]] && \
    grep -Eq '/nas_data|/evaluator_only|/target_complete|/truth_references|/complete($|/)|/labels($|/)' \
        /proc/self/mountinfo; then
    echo "forbidden evaluator, complete, truth, label, or NAS mount detected" >&2
    exit 1
fi
for required in "$python_bin" "$context_verifier" "$contract_path" \
    "$input_root/role_manifest.tsv" "$input_root/role_contract.json" \
    "$blind_benchmark_root/perturbed.gff3" "$blind_benchmark_root/blind_manifest.json" \
    "$blind_benchmark_root/SHA256SUMS" "$protocol_root/SHA256SUMS" \
    "$execution_root/SHA256SUMS" "$composite_root/SHA256SUMS"; do
    [[ -s $required ]] || { echo "missing frozen blind pipeline input: $required" >&2; exit 1; }
done
(cd "$protocol_root" && sha256sum -c SHA256SUMS >/dev/null)
(cd "$execution_root" && sha256sum -c SHA256SUMS >/dev/null)
(cd "$composite_root" && sha256sum -c SHA256SUMS >/dev/null)
PYTHONPATH="$code_root/src" "$python_bin" "$context_verifier" \
    --input-root "$input_root" --contract "$contract_path" \
    --protocol-freeze "$protocol_root" --execution-freeze "$execution_root" \
    --blind-benchmark-root "$blind_benchmark_root" \
    --expected-holdout-id actinidia_red5_v0.5 --expected-primary-chromosomes 29 \
    --model-freeze "$composite_root" \
    --output-json "$blind_output_root/pipeline_blind_context.json"

[[ ! -e $command_log && ! -e $command_log_working ]] || {
    echo "refusing to overwrite frozen pipeline command log" >&2
    exit 1
}
implementation_sha() {
    local relative=$1 expected observed
    expected=$(awk -F '\t' -v path="$relative" \
        'NR > 1 && $1 == path {if (seen++) exit 3; print $3}' \
        "$execution_root/implementation_manifest.tsv")
    observed=$(sha256sum "$code_root/$relative" | awk '{print $1}')
    [[ $expected =~ ^[0-9a-f]{64}$ && $observed == "$expected" ]] || {
        echo "pipeline implementation differs from execution freeze: $relative" >&2
        exit 1
    }
    printf '%s' "$observed"
}
security_env="PLOIDYPATCH_BLIND_RUNNER=1;PLOIDYPATCH_NETWORK_ACCESS=none;PLOIDYPATCH_STAGED_INPUT_ROOT=$input_root;PLOIDYPATCH_BLIND_BENCHMARK_ROOT=$blind_benchmark_root;PLOIDYPATCH_PROTOCOL_FREEZE=$protocol_root;PLOIDYPATCH_EXECUTION_FREEZE=$execution_root;PLOIDYPATCH_COMPOSITE_MODEL_FREEZE=$composite_root"
for relative in \
    scripts/verify_external_holdout_blind_context_v0.5.py \
    scripts/run_actinidia_blind_pipeline_v0.5.sh \
    scripts/run_actinidia_miniprot_upstream_v0.5.sh \
    scripts/run_gemoma_homology.sh scripts/run_lifton_transfer.sh \
    scripts/build_actinidia_method_trio_candidate_pools_v0.5.sh \
    scripts/run_actinidia_blind_union_self_wgd_v0.5.sh \
    scripts/score_actinidia_candidates_blind_v0.5.sh; do
    implementation_sha "$relative" >/dev/null
done
{
    printf 'stage_order\tinvocation_scope\tfrozen_relative_script\tscript_sha256\tsafety_environment\n'
    while IFS=$'\t' read -r stage scope relative; do
        printf '%s\t%s\t%s\t%s\t%s\n' "$stage" "$scope" "$relative" \
            "$(implementation_sha "$relative")" "$security_env"
    done <<'COMMANDS'
0	context_check	scripts/verify_external_holdout_blind_context_v0.5.py
1	pipeline	scripts/run_actinidia_blind_pipeline_v0.5.sh
10	candidate_upstream	scripts/run_actinidia_miniprot_upstream_v0.5.sh
20	candidate_eriantha	scripts/run_gemoma_homology.sh
20	candidate_rufa	scripts/run_gemoma_homology.sh
20	candidate_eriantha	scripts/run_lifton_transfer.sh
20	candidate_rufa	scripts/run_lifton_transfer.sh
30	method_family_pool	scripts/build_actinidia_method_trio_candidate_pools_v0.5.sh
40	blind_self_wgd	scripts/run_actinidia_blind_union_self_wgd_v0.5.sh
50	frozen_v04_score	scripts/score_actinidia_candidates_blind_v0.5.sh
COMMANDS
} > "$command_log_working"
if grep -Eqi '/nas_data|evaluator|truth|labels|target_complete' "$command_log_working"; then
    rm -f "$command_log_working"
    echo "pipeline command log contains a forbidden role token" >&2
    exit 1
fi
mv "$command_log_working" "$command_log"

bash "$code_root/scripts/run_actinidia_miniprot_upstream_v0.5.sh" "$project_root"
[[ -s $target && -s ${target}.fai ]] || {
    echo "candidate-safe target FASTA/index missing after normalization" >&2
    exit 1
}
"$python_bin" - "$target" "$blind_benchmark_root/blind_manifest.json" <<'PY'
import hashlib
import json
import sys
from pathlib import Path

target, manifest_path = map(Path, sys.argv[1:])
digest = hashlib.sha256()
with target.open("rb") as handle:
    while block := handle.read(8 * 1024 * 1024):
        digest.update(block)
manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
if manifest.get("target_genome", {}).get("sha256") != digest.hexdigest():
    raise SystemExit("normalized Red5 target differs from sealed blind benchmark")
PY

pids=()
labels=()
for reference in actinidia_eriantha actinidia_rufa; do
    reference_fasta=$normalized/$reference/primary_chromosomes.genome.fa
    reference_gff=$normalized/$reference/primary_chromosomes.gff3
    lifton_reference_gff=$normalized/$reference/primary_chromosomes.lifton.gff3
    for input in "$reference_fasta" "$reference_gff" "$lifton_reference_gff"; do
        [[ -s $input ]] || { echo "missing normalized candidate reference: $input" >&2; exit 1; }
    done
    (
        bash "$code_root/scripts/run_gemoma_homology.sh" \
            "$project_root/envs/ploidypatch-gemoma" \
            "$baseline_root/gemoma/$reference" \
            "$target" "$reference_fasta" "$reference_gff" "$reference" 32
    ) &
    pids+=("$!")
    labels+=("gemoma:$reference")
    (
        bash "$code_root/scripts/run_lifton_transfer.sh" \
            "$project_root/envs/ploidypatch-lifton" \
            "$baseline_root/lifton/$reference" \
            "$target" "$reference_fasta" "$lifton_reference_gff" "$reference" 16
    ) &
    pids+=("$!")
    labels+=("lifton:$reference")
done
failed=0
for index in "${!pids[@]}"; do
    if ! wait "${pids[$index]}"; then
        printf 'blind upstream failed: %s\n' "${labels[$index]}" >&2
        failed=1
    fi
done
[[ $failed -eq 0 ]] || exit 1

bash "$code_root/scripts/build_actinidia_method_trio_candidate_pools_v0.5.sh" "$project_root"
bash "$code_root/scripts/run_actinidia_blind_union_self_wgd_v0.5.sh" "$project_root"
bash "$code_root/scripts/score_actinidia_candidates_blind_v0.5.sh" "$project_root"

printf 'Actinidia v0.5 blind pipeline completed inside isolated namespace\n'
