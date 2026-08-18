#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "usage: $0 PROJECT_ROOT" >&2
    exit 2
fi

project_root=$(realpath "$1")
config=${PLOIDYPATCH_STRUCTURE_CONFIG:-$project_root/code/config/structure_benchmark_v0.1.tsv}
runner=$project_root/code/scripts/run_structure_benchmark_one.sh
parallel_bin=${PLOIDYPATCH_PARALLEL_BIN:-$(command -v parallel || true)}
if [[ -z $parallel_bin ]] && command -v conda >/dev/null 2>&1; then
    conda_base=$(conda info --base)
    if [[ -x $conda_base/bin/parallel ]]; then
        parallel_bin=$conda_base/bin/parallel
    fi
fi
for required in "$config" "$runner"; do
    if [[ ! -s $required ]]; then
        echo "missing structure matrix prerequisite: $required" >&2
        exit 1
    fi
done
if [[ -z $parallel_bin ]]; then
    echo "GNU parallel is required for the frozen structure matrix" >&2
    exit 1
fi

tail -n +2 "$config" | "$parallel_bin" \
    --colsep '\t' --jobs 4 --delay 0.5 --halt soon,fail=1 \
    env PLOIDYPATCH_CODE_COMMIT="${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}" \
    PLOIDYPATCH_STRUCTURE_NAMESPACE="${PLOIDYPATCH_STRUCTURE_NAMESPACE:-v0.1}" \
    bash "$runner" "$project_root" {1} {2} {3} {4} {5}
printf 'structure benchmark matrix completed\n'
