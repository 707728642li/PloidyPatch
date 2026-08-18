#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "usage: $0 PROJECT_ROOT" >&2
    exit 2
fi

project_root=$(realpath "$1")
runner=$project_root/code/scripts/run_glycine_candidate_self_wgd_one_v0.1.sh
result_root=$project_root/results/copy_collapse/wgd_reanchor_glycine_v0.1/self_wgdi
state_root=$project_root/logs/copy_collapse/wgd_reanchor_glycine_v0.1

if [[ ! -s $runner ]]; then
    echo "run script is absent: $runner" >&2
    exit 2
fi
if [[ -e $state_root || -e $result_root/blind || \
      -e $result_root/complete_control ]]; then
    echo "refusing to reuse candidate self-WGD launcher state or results" >&2
    exit 2
fi
mkdir -p "$state_root"

for mode in blind complete_control; do
    nohup setsid env PLOIDYPATCH_CODE_COMMIT="${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}" \
        bash "$runner" "$project_root" "$mode" \
        > "$state_root/$mode.launcher.log" 2>&1 &
    pid=$!
    printf '%s\n' "$pid" > "$state_root/$mode.pid"
done

sleep 2
for mode in blind complete_control; do
    pid=$(cat "$state_root/$mode.pid")
    if ! kill -0 "$pid"; then
        echo "candidate self-WGD $mode launcher exited during startup" >&2
        exit 1
    fi
    printf '%s_pid\t%s\n' "$mode" "$pid"
done
