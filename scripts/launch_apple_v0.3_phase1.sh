#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then echo "usage: $0 PROJECT_ROOT" >&2; exit 2; fi
project_root=$(realpath "$1")
code_root=$project_root/code
source_root=$project_root/data/derived/external_inputs/apple_v0.3
protocol_root=$project_root/results/protocol_freezes/apple_external_v0.3
state_root=$project_root/logs/apple_v0.3/phase1
target=$source_root/target_apple/primary_chromosomes.genome.fa
gemoma_env=$project_root/envs/ploidypatch-gemoma
lifton_env=$project_root/envs/ploidypatch-lifton

for required in "$protocol_root/SHA256SUMS" "$target" "${target}.fai" \
    "$code_root/scripts/prepare_apple_evaluator_wgdi_inputs_v0.3.sh" \
    "$code_root/scripts/run_apple_miniprot_upstream_v0.3.sh" \
    "$code_root/scripts/run_gemoma_homology.sh" \
    "$code_root/scripts/run_lifton_transfer.sh"; do
    [[ -s $required ]] || { echo "missing apple phase-1 prerequisite: $required" >&2; exit 1; }
done
(cd "$protocol_root" && sha256sum -c SHA256SUMS >/dev/null)
[[ ! -e $state_root ]] || { echo "refusing to reuse apple phase-1 state" >&2; exit 1; }
mkdir -p "$state_root"
{
    printf 'field\tvalue\n'
    printf 'code_commit\t%s\n' "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'truth_pair_access\tfalse\nexternal_label_access\tfalse\n'
    printf 'parallel_cpu_budget\t128\n'
    printf 'protocol_freeze_sha256sums_sha256\t%s\n' \
        "$(sha256sum "$protocol_root/SHA256SUMS" | awk '{print $1}')"
} > "$state_root/run_contract.tsv"

start_job() {
    local name=$1; shift
    nohup setsid env PLOIDYPATCH_CODE_COMMIT="${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}" \
        "$@" > "$state_root/$name.log" 2>&1 < /dev/null &
    printf '%s\n' "$!" > "$state_root/$name.pid"
}

start_job evaluator_inputs bash "$code_root/scripts/prepare_apple_evaluator_wgdi_inputs_v0.3.sh" "$project_root"
start_job miniprot bash "$code_root/scripts/run_apple_miniprot_upstream_v0.3.sh" "$project_root"
for reference in pear peach; do
    case $reference in
        pear) bundle=candidate_pear; reference_id=pyrus_bartlettdh_v2 ;;
        peach) bundle=candidate_peach; reference_id=prunus_ncbiv2 ;;
    esac
    ref_root=$source_root/$bundle
    start_job "gemoma_$reference" bash "$code_root/scripts/run_gemoma_homology.sh" \
        "$gemoma_env" "$project_root/results/baselines/apple_v0.3/gemoma/$reference" \
        "$target" "$ref_root/primary_chromosomes.genome.fa" \
        "$ref_root/primary_chromosomes.gff3" "$reference_id" 32
    start_job "lifton_$reference" bash "$code_root/scripts/run_lifton_transfer.sh" \
        "$lifton_env" "$project_root/results/baselines/apple_v0.3/lifton/$reference" \
        "$target" "$ref_root/primary_chromosomes.genome.fa" \
        "$ref_root/primary_chromosomes.gff3" "$reference_id" 16
done

sleep 3
printf 'job\tpid\tstate\n' > "$state_root/startup.tsv"
for pid_file in "$state_root"/*.pid; do
    name=$(basename "$pid_file" .pid); pid=$(cat "$pid_file")
    if kill -0 "$pid" 2>/dev/null; then state=running; else state=exited_during_startup; fi
    printf '%s\t%s\t%s\n' "$name" "$pid" "$state" >> "$state_root/startup.tsv"
    [[ $state == running ]] || exit 1
done
cat "$state_root/startup.tsv"
