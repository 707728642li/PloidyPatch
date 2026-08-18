#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then echo "usage: $0 PROJECT_ROOT" >&2; exit 2; fi
project_root=$(realpath "$1")
code_root=$project_root/code
prepared=$project_root/data/derived/candidate_inputs/cotton_syngap_v0.1
prepare_pid_file=$project_root/logs/baseline/cotton_holdout_syngap_prepare_v0.1/pid
launch_state=$project_root/logs/baseline/cotton_holdout_syngap_v0.1
upstream=$project_root/results/baselines/cotton_holdout_v0.1/syngap
state=$project_root/logs/copy_collapse/cotton_syngap_pipeline_v0.1
if [[ -e $state ]]; then echo "refusing to reuse cotton SynGAP pipeline state" >&2; exit 1; fi
mkdir -p "$state"
printf 'started\t%s\ncode_commit\t%s\n' "$(date --iso-8601=seconds)" \
    "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}" > "$state/status.tsv"
while [[ ! -s $prepared/SHA256SUMS ]]; do
    pid=$(cat "$prepare_pid_file")
    kill -0 "$pid" 2>/dev/null || { echo "cotton SynGAP preparation stopped" >&2; exit 1; }
    sleep 30
done
printf 'inputs_frozen\t%s\n' "$(date --iso-8601=seconds)" >> "$state/status.tsv"
env PLOIDYPATCH_CODE_COMMIT="${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}" \
    bash "$code_root/scripts/launch_cotton_syngap_v0.1.sh" "$project_root" \
    > "$state/launch.log" 2>&1
printf 'upstreams_launched\t%s\n' "$(date --iso-8601=seconds)" >> "$state/status.tsv"
for mode in blind complete_control; do
    for ref in gar_a gra_d; do
        manifest=$upstream/$ref/$mode/output_manifest.tsv
        pid_file=$launch_state/$ref.$mode.pid
        while [[ ! -s $manifest ]]; do
            pid=$(cat "$pid_file")
            kill -0 "$pid" 2>/dev/null || { echo "cotton SynGAP arm stopped: $ref $mode" >&2; exit 1; }
            sleep 30
        done
        printf '%s_%s_complete\t%s\n' "$ref" "$mode" \
            "$(date --iso-8601=seconds)" >> "$state/status.tsv"
    done
done
env PLOIDYPATCH_CODE_COMMIT="${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}" \
    bash "$code_root/scripts/evaluate_cotton_syngap_v0.1.sh" "$project_root" \
    > "$state/evaluation.log" 2>&1
printf 'evaluation_complete\t%s\n' "$(date --iso-8601=seconds)" >> "$state/status.tsv"
env PLOIDYPATCH_CODE_COMMIT="${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}" \
    bash "$code_root/scripts/run_cotton_extended_statistics_when_ready_v0.1.sh" "$project_root" \
    > "$state/extended_statistics.log" 2>&1
printf 'extended_statistics_complete\t%s\ncompleted\t%s\n' \
    "$(date --iso-8601=seconds)" "$(date --iso-8601=seconds)" >> "$state/status.tsv"
