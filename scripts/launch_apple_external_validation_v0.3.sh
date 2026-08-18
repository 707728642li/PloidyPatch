#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then echo "usage: $0 PROJECT_ROOT" >&2; exit 2; fi
project_root=$(realpath "$1")
code_root=$project_root/code
export PYTHONPATH=$code_root/src${PYTHONPATH:+:$PYTHONPATH}
execution_root=$project_root/results/protocol_freezes/apple_external_v0.3_execution
state_root=$project_root/logs/apple_v0.3/external_validation
state=$state_root/stages.tsv
max_start_load=${PLOIDYPATCH_MAX_START_LOAD:-160}

[[ -s $execution_root/SHA256SUMS ]] || {
    echo "missing apple execution freeze" >&2; exit 1;
}
(cd "$execution_root" && sha256sum -c SHA256SUMS >/dev/null)
frozen_commit=$(awk -F '\t' '$1 == "code_commit" {print $2}' \
    "$execution_root/run_contract.tsv")
[[ $frozen_commit =~ ^[0-9a-f]{40}$ ]] || {
    echo "execution freeze lacks a concrete code commit" >&2; exit 1;
}
export PLOIDYPATCH_CODE_COMMIT=$frozen_commit
mkdir -p "$state_root"
exec 9> "$state_root/orchestrator.lock"
flock -n 9 || { echo "apple external orchestrator is already running" >&2; exit 1; }
if [[ ! -e $state ]]; then
    printf 'timestamp\tstage\tstatus\tdetail\n' > "$state"
fi
record() {
    printf '%s\t%s\t%s\t%s\n' "$(date -Is)" "$1" "$2" "$3" >> "$state"
}
wait_for_capacity() {
    local stage=$1 load
    while true; do
        load=$(awk '{print $1}' /proc/loadavg)
        if awk -v load="$load" -v maximum="$max_start_load" \
            'BEGIN {exit !(load < maximum)}'; then
            record "$stage" capacity "load1=$load threshold=$max_start_load"
            return
        fi
        record "$stage" waiting_capacity "load1=$load threshold=$max_start_load"
        sleep 30
    done
}
wait_for_upstream() {
    local label=$1 output=$2 pid_file=$3 pid
    while [[ ! -s $output ]]; do
        if [[ -s $pid_file ]]; then
            pid=$(cat "$pid_file")
            if ! ps -p "$pid" >/dev/null 2>&1; then
                record "$label" failed "upstream exited without $output"
                return 1
            fi
        fi
        record "$label" waiting_upstream "$output"
        sleep 30
    done
    record "$label" ready "$output"
}
run_stage() {
    local label=$1 output=$2
    shift 2
    if [[ -s $output/SHA256SUMS ]]; then
        (cd "$output" && sha256sum -c SHA256SUMS >/dev/null)
        record "$label" reused "$output"
        return
    fi
    [[ ! -e $output && ! -e ${output}.working ]] || {
        record "$label" failed "incomplete or unhashed output already exists"
        return 1
    }
    record "$label" started "$*"
    "$@"
    [[ -s $output/SHA256SUMS ]] || {
        record "$label" failed "stage returned without SHA256SUMS"
        return 1
    }
    (cd "$output" && sha256sum -c SHA256SUMS >/dev/null)
    record "$label" completed "$output"
}

wgdi=$project_root/results/evaluator/apple_v0.3/wgdi
truth_pairs=$project_root/results/evaluator/apple_v0.3/truth_pairs
benchmark=$project_root/benchmark/structure/copy_collapse_v0.3/mdx_gddh13/annotation_copy_collapse_seed20260831
method_pool=$project_root/results/copy_collapse/external/apple_v0.3_method_trio
blind_wgd=$project_root/results/copy_collapse/external/apple_v0.3_blind_self_wgd
rankings=$project_root/results/copy_collapse/external/apple_v0.3_blind_rankings
reveal=$project_root/results/copy_collapse/external/apple_v0.3_reveal

record orchestrator started "frozen_commit=$frozen_commit max_start_load=$max_start_load"
wait_for_capacity evaluator_wgdi
run_stage evaluator_wgdi "$wgdi" \
    "$code_root/scripts/run_apple_evaluator_wgdi_v0.3.sh" "$project_root"
run_stage truth_pairs "$truth_pairs" \
    "$code_root/scripts/infer_apple_external_pairs_v0.3.sh" "$project_root"
run_stage hidden_benchmark "$benchmark" \
    "$code_root/scripts/run_apple_copy_collapse_benchmark_v0.3.sh" "$project_root"

phase1=$project_root/logs/apple_v0.3/phase1
wait_for_upstream gemoma_pear \
    "$project_root/results/baselines/apple_v0.3/gemoma/pear/upstream/final_annotation.gff" \
    "$phase1/gemoma_pear.pid"
wait_for_upstream gemoma_peach \
    "$project_root/results/baselines/apple_v0.3/gemoma/peach/upstream/final_annotation.gff" \
    "$phase1/gemoma_peach.pid"
wait_for_upstream lifton_pear \
    "$project_root/results/baselines/apple_v0.3/lifton/pear/upstream/lifton.gff3" \
    "$phase1/lifton_pear.pid"
wait_for_upstream lifton_peach \
    "$project_root/results/baselines/apple_v0.3/lifton/peach/upstream/lifton.gff3" \
    "$phase1/lifton_peach.pid"
run_stage method_pool "$method_pool" \
    "$code_root/scripts/build_apple_method_trio_candidate_pools_v0.3.sh" "$project_root"
wait_for_capacity blind_self_wgdi
run_stage blind_self_wgdi "$blind_wgd" \
    "$code_root/scripts/run_apple_blind_union_self_wgd_v0.3.sh" "$project_root"
run_stage blind_rankings "$rankings" \
    "$code_root/scripts/score_apple_candidates_blind_v0.3.sh" "$project_root"
run_stage reveal "$reveal" \
    "$code_root/scripts/run_apple_external_reveal_v0.3.sh" "$project_root"
record orchestrator completed "$reveal"
printf 'apple external v0.3 pipeline complete: %s\n' "$reveal"
