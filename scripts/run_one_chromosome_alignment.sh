#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 7 ]]; then
    echo "usage: $0 ENV_PREFIX INDEX QUERY OUTPUT_DIR LOG_DIR PRESET THREADS" >&2
    exit 2
fi

env_prefix=$1
index=$2
query=$3
output_dir=$4
log_dir=$5
preset=$6
threads=$7
stem=$(basename "$query" .fa)
output="$output_dir/$stem.paf"
log="$log_dir/$stem.log"

if [[ -e "$output" || -e "$log" ]]; then
    echo "refusing to overwrite alignment artifact for $stem" >&2
    exit 1
fi

/usr/bin/time -v conda run -p "$env_prefix" --no-capture-output \
    minimap2 -x "$preset" -t"$threads" --secondary=no -c --cs=short \
    "$index" "$query" > "$output" 2> "$log"

