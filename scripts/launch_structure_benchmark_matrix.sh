#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "usage: $0 PROJECT_ROOT ATTEMPT_LABEL" >&2
    exit 2
fi

project_root=$(realpath "$1")
attempt=$2
if [[ ! $attempt =~ ^[A-Za-z0-9][A-Za-z0-9_.-]*$ ]]; then
    echo "attempt label must be a safe identifier" >&2
    exit 2
fi
runner=$project_root/code/scripts/run_structure_benchmark_matrix.sh
result_root=$project_root/benchmark/structure/v0.1
state_root=$project_root/logs/benchmark/structure_matrix_v0.1
log=$state_root/${attempt}.launcher.log
pid_file=$state_root/${attempt}.pid
if [[ ! -s $runner ]]; then
    echo "structure matrix runner is absent: $runner" >&2
    exit 2
fi
if [[ -e $result_root ]]; then
    echo "refusing to reuse structure matrix root: $result_root" >&2
    exit 2
fi
if [[ -e $log || -e $pid_file ]]; then
    echo "refusing to overwrite structure matrix attempt: $attempt" >&2
    exit 2
fi
mkdir -p "$state_root"
nohup setsid env \
    PLOIDYPATCH_CODE_COMMIT="${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}" \
    bash "$runner" "$project_root" > "$log" 2>&1 &
pid=$!
printf '%s\n' "$pid" > "$pid_file"
sleep 2
if ! kill -0 "$pid"; then
    echo "structure matrix exited during startup; inspect $log" >&2
    exit 1
fi
printf 'structure_matrix_pid\t%s\nattempt\t%s\n' "$pid" "$attempt"
