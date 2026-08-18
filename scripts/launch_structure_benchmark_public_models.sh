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
config=$project_root/code/config/structure_benchmark_public_models_v0.1.tsv
source_root=$project_root/data/derived/structure_sources_v0.1
result_root=$project_root/benchmark/structure/public_models_v0.1
state_root=$project_root/logs/benchmark/structure_public_models_v0.1
log=$state_root/${attempt}.launcher.log
pid_file=$state_root/${attempt}.pid
for required in "$runner" "$config" "$source_root/ath_tair10/source.gff3" \
                "$source_root/osa_irgsp10/source.gff3"; do
    if [[ ! -s $required ]]; then
        echo "public-model structure prerequisite is absent: $required" >&2
        exit 2
    fi
done
if [[ -e $result_root ]]; then
    echo "refusing to reuse public-model structure root: $result_root" >&2
    exit 2
fi
if [[ -e $log || -e $pid_file ]]; then
    echo "refusing to overwrite public-model attempt: $attempt" >&2
    exit 2
fi
mkdir -p "$state_root"
nohup setsid env \
    PLOIDYPATCH_CODE_COMMIT="${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}" \
    PLOIDYPATCH_STRUCTURE_NAMESPACE=public_models_v0.1 \
    PLOIDYPATCH_STRUCTURE_CONFIG="$config" \
    bash "$runner" "$project_root" > "$log" 2>&1 &
pid=$!
printf '%s\n' "$pid" > "$pid_file"
sleep 2
if ! kill -0 "$pid"; then
    echo "public-model structure matrix exited during startup; inspect $log" >&2
    exit 1
fi
printf 'structure_public_models_pid\t%s\nattempt\t%s\n' "$pid" "$attempt"
