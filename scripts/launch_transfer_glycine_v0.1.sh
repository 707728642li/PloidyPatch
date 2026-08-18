#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "usage: $0 PROJECT_ROOT" >&2
    exit 2
fi

project_root=$(realpath "$1")
runner=$project_root/code/scripts/run_transfer_glycine_v0.1.sh
state_root=$project_root/logs/baseline/transfer_glycine_v0.1
gemoma_root=$project_root/results/baselines/gemoma_v1.9/glycine_v0.1/raw
lifton_root=$project_root/results/baselines/lifton_v1.0.11/glycine_v0.1/raw

if [[ ! -s $runner ]]; then
    echo "run script is absent: $runner" >&2
    exit 2
fi
if [[ -e $state_root ]]; then
    echo "refusing to reuse launcher state: $state_root" >&2
    exit 2
fi
for run_root in "$gemoma_root" "${gemoma_root}.working" \
                "$lifton_root" "${lifton_root}.working"; do
    if [[ -e $run_root ]]; then
        echo "refusing to reuse transfer result root: $run_root" >&2
        exit 2
    fi
done
mkdir -p "$state_root"

nohup setsid bash "$runner" "$project_root" gemoma \
    > "$state_root/gemoma.launcher.log" 2>&1 &
gemoma_pid=$!
printf '%s\n' "$gemoma_pid" > "$state_root/gemoma.pid"

nohup setsid bash "$runner" "$project_root" lifton \
    > "$state_root/lifton.launcher.log" 2>&1 &
lifton_pid=$!
printf '%s\n' "$lifton_pid" > "$state_root/lifton.pid"

sleep 2
if ! kill -0 "$gemoma_pid" || ! kill -0 "$lifton_pid"; then
    echo "one or more transfer launchers exited during startup" >&2
    exit 1
fi
printf 'gemoma_pid\t%s\nlifton_pid\t%s\n' "$gemoma_pid" "$lifton_pid"
