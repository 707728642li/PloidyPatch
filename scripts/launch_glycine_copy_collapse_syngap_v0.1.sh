#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "usage: $0 PROJECT_ROOT" >&2
    exit 2
fi

project_root=$(realpath "$1")
runner=$project_root/code/scripts/run_glycine_copy_collapse_syngap_v0.1.sh
output_root=$project_root/results/copy_collapse/syngap_glycine_v0.1/upstream/blind
state_root=$project_root/logs/copy_collapse/syngap_glycine_v0.1

if [[ ! -s $runner ]]; then
    echo "run script is absent: $runner" >&2
    exit 2
fi
if [[ -e $state_root ]]; then
    echo "refusing to reuse launcher state: $state_root" >&2
    exit 2
fi
if [[ -e $output_root ]]; then
    echo "refusing to reuse SynGAP copy-collapse result root" >&2
    exit 2
fi
mkdir -p "$state_root"

nohup setsid bash "$runner" "$project_root" \
    > "$state_root/launcher.log" 2>&1 &
pid=$!
printf '%s\n' "$pid" > "$state_root/pid"

sleep 2
if ! kill -0 "$pid"; then
    echo "SynGAP copy-collapse launcher exited during startup" >&2
    exit 1
fi
printf 'pid\t%s\n' "$pid"
