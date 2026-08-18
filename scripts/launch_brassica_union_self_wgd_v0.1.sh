#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then echo "usage: $0 PROJECT_ROOT" >&2; exit 2; fi
project_root=$(realpath "$1")
runner=$project_root/code/scripts/run_brassica_union_self_wgd_one_v0.1.sh
result_root=$project_root/results/copy_collapse/model_development/brassica_union_self_wgd_v0.1
state_root=$project_root/logs/copy_collapse/brassica_union_self_wgd_v0.1
if [[ -e $state_root || -e $result_root/blind || -e $result_root/complete_control ]]; then
    echo "refusing to reuse Brassica self-WGD state" >&2; exit 1
fi
mkdir -p "$state_root"
for mode in blind complete_control; do
    nohup setsid env PLOIDYPATCH_CODE_COMMIT="${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}" \
        bash "$runner" "$project_root" "$mode" \
        > "$state_root/$mode.launcher.log" 2>&1 &
    printf '%s\n' "$!" > "$state_root/$mode.pid"
done
sleep 3
for mode in blind complete_control; do
    pid=$(cat "$state_root/$mode.pid")
    if ! kill -0 "$pid" 2>/dev/null; then
        echo "Brassica self-WGD $mode exited during startup" >&2; exit 1
    fi
    printf '%s_pid\t%s\n' "$mode" "$pid"
done
