#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 7 ]]; then
    echo "usage: $0 ENV_PREFIX INDEX QUERY_DIR OUTPUT_DIR LOG_DIR PRESET JOBS" >&2
    exit 2
fi

export PP_ENV_PREFIX=$1
export PP_INDEX=$2
export PP_QUERY_DIR=$3
export PP_OUTPUT_DIR=$4
export PP_LOG_DIR=$5
export PP_PRESET=$6
export PP_THREADS_PER_JOB=6
jobs=$7
script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)

mkdir -p "$PP_OUTPUT_DIR" "$PP_LOG_DIR"

find "$PP_QUERY_DIR" -maxdepth 1 -type f -name '*.fa' -print0 \
    | sort -z \
    | conda run -n base --no-capture-output parallel \
        -0 --jobs "$jobs" --delay 0.3 --halt soon,fail=1 \
        "$script_dir/run_one_chromosome_alignment.sh" \
        "$PP_ENV_PREFIX" "$PP_INDEX" {} "$PP_OUTPUT_DIR" "$PP_LOG_DIR" \
        "$PP_PRESET" "$PP_THREADS_PER_JOB"
