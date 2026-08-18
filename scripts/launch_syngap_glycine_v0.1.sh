#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "usage: $0 PROJECT_ROOT" >&2
    exit 2
fi

project_root=$(realpath "$1")
runner=$project_root/code/scripts/run_syngap_glycine_v0.1.sh
output_root=$project_root/results/baselines/syngap_v1.2.5/glycine_v0.1/genblastg
state_root=$project_root/logs/baseline/syngap_glycine_v0.1_genblastg

if [[ ! -s $runner ]]; then
    echo "run script is absent: $runner" >&2
    exit 2
fi
if [[ -e $state_root ]]; then
    echo "refusing to reuse launcher state: $state_root" >&2
    exit 2
fi
if [[ -e $output_root/blind || -e $output_root/complete_control ]]; then
    echo "refusing to reuse SynGAP result roots" >&2
    exit 2
fi
mkdir -p "$state_root"

nohup setsid bash "$runner" "$project_root" blind \
    > "$state_root/blind.launcher.log" 2>&1 &
blind_pid=$!
printf '%s\n' "$blind_pid" > "$state_root/blind.pid"

nohup setsid bash "$runner" "$project_root" complete_control \
    > "$state_root/complete_control.launcher.log" 2>&1 &
control_pid=$!
printf '%s\n' "$control_pid" > "$state_root/complete_control.pid"

sleep 2
if ! kill -0 "$blind_pid" || ! kill -0 "$control_pid"; then
    echo "one or more SynGAP launchers exited during startup" >&2
    exit 1
fi
printf 'blind_pid\t%s\ncomplete_control_pid\t%s\n' \
    "$blind_pid" "$control_pid"
