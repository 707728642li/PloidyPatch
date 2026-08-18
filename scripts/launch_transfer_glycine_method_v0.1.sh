#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 3 ]]; then
    echo "usage: $0 PROJECT_ROOT gemoma|lifton ATTEMPT_LABEL" >&2
    exit 2
fi

project_root=$(realpath "$1")
method=$2
attempt=$3
if [[ $method != gemoma && $method != lifton ]]; then
    echo "method must be gemoma or lifton" >&2
    exit 2
fi
if [[ ! $attempt =~ ^[A-Za-z0-9][A-Za-z0-9_.-]*$ ]]; then
    echo "attempt label must be a safe identifier" >&2
    exit 2
fi

runner=$project_root/code/scripts/run_transfer_glycine_v0.1.sh
state_root=$project_root/logs/baseline/transfer_glycine_v0.1
case $method in
    gemoma)
        run_root=$project_root/results/baselines/gemoma_v1.9/glycine_v0.1/raw
        ;;
    lifton)
        run_root=$project_root/results/baselines/lifton_v1.0.11/glycine_v0.1/raw
        ;;
esac
log=$state_root/${method}.${attempt}.launcher.log
pid_file=$state_root/${method}.${attempt}.pid

if [[ ! -s $runner ]]; then
    echo "run script is absent: $runner" >&2
    exit 2
fi
if [[ -e $run_root || -e ${run_root}.working ]]; then
    echo "refusing to reuse a completed or working result: $run_root" >&2
    exit 2
fi
if [[ -e $log || -e $pid_file ]]; then
    echo "refusing to overwrite attempt state: $method $attempt" >&2
    exit 2
fi
mkdir -p "$state_root"

nohup setsid bash "$runner" "$project_root" "$method" > "$log" 2>&1 &
pid=$!
printf '%s\n' "$pid" > "$pid_file"
sleep 2
if ! kill -0 "$pid"; then
    echo "$method exited during startup; inspect $log" >&2
    exit 1
fi
printf '%s_pid\t%s\nattempt\t%s\n' "$method" "$pid" "$attempt"
