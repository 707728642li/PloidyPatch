#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then echo "usage: $0 PROJECT_ROOT" >&2; exit 2; fi
project_root=$(realpath "$1")
code_root=$project_root/code
method_root=$project_root/results/copy_collapse/holdout/cotton_method_trio_v0.1
method_watcher_pid_file=$project_root/logs/copy_collapse/cotton_candidate_pipeline_v0.1.pid
state_root=$project_root/logs/copy_collapse/cotton_post_methods_v0.1
if [[ -e $state_root ]]; then echo "refusing to reuse cotton post-method state" >&2; exit 1; fi
mkdir -p "$state_root"
printf 'started\t%s\ncode_commit\t%s\n' "$(date --iso-8601=seconds)" \
    "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}" > "$state_root/status.tsv"

while [[ ! -s $method_root/SHA256SUMS ]]; do
    if [[ -s $method_watcher_pid_file ]]; then
        method_watcher_pid=$(cat "$method_watcher_pid_file")
        if ! kill -0 "$method_watcher_pid" 2>/dev/null; then
            printf 'failed\t%s\nreason\tmethod_watcher_stopped\n' \
                "$(date --iso-8601=seconds)" >> "$state_root/status.tsv"
            exit 1
        fi
    fi
    sleep 30
done
(cd "$method_root" && sha256sum -c SHA256SUMS >/dev/null)
printf 'method_trio_frozen\t%s\n' "$(date --iso-8601=seconds)" >> "$state_root/status.tsv"

pids=(); labels=()
for mode in blind complete_control; do
    env PLOIDYPATCH_CODE_COMMIT="${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}" \
        bash "$code_root/scripts/run_cotton_union_self_wgd_one_v0.1.sh" \
        "$project_root" "$mode" > "$state_root/self_wgd.$mode.log" 2>&1 &
    pids+=("$!"); labels+=("$mode")
done
failed=0
for i in "${!pids[@]}"; do
    if ! wait "${pids[$i]}"; then
        printf 'self_wgd_%s_failed\t%s\n' "${labels[$i]}" \
            "$(date --iso-8601=seconds)" >> "$state_root/status.tsv"
        failed=1
    else
        printf 'self_wgd_%s_completed\t%s\n' "${labels[$i]}" \
            "$(date --iso-8601=seconds)" >> "$state_root/status.tsv"
    fi
done
[[ $failed -eq 0 ]] || exit 1

env PLOIDYPATCH_CODE_COMMIT="${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}" \
    bash "$code_root/scripts/evaluate_cotton_glycine_copy_model_transfer_v0.1.sh" \
    "$project_root" > "$state_root/model_transfer.log" 2>&1
printf 'model_transfer_completed\t%s\n' "$(date --iso-8601=seconds)" >> "$state_root/status.tsv"
env PLOIDYPATCH_CODE_COMMIT="${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}" \
    bash "$code_root/scripts/run_cotton_full_statistics_v0.1.sh" \
    "$project_root" > "$state_root/statistics.log" 2>&1
printf 'statistics_completed\t%s\ncompleted\t%s\n' \
    "$(date --iso-8601=seconds)" "$(date --iso-8601=seconds)" >> "$state_root/status.tsv"
