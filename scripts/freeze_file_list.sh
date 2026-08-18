#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 3 ]]; then
    echo "usage: $0 ROOT FILE_LIST OUTPUT_SHA256SUMS" >&2
    exit 2
fi

root=$(cd -- "$1" && pwd)
list=$2
output=$3

if [[ -e "$output" ]]; then
    echo "refusing to overwrite checksum artifact: $output" >&2
    exit 1
fi

mkdir -p "$(dirname -- "$output")"
while IFS= read -r relative_path || [[ -n "$relative_path" ]]; do
    relative_path=${relative_path%$'\r'}
    [[ -z "$relative_path" || "$relative_path" == \#* ]] && continue
    if [[ "$relative_path" = /* || "$relative_path" == *".."* ]]; then
        echo "unsafe relative path in freeze list: $relative_path" >&2
        exit 1
    fi
    if [[ ! -f "$root/$relative_path" ]]; then
        echo "missing freeze input: $root/$relative_path" >&2
        exit 1
    fi
    (cd -- "$root" && sha256sum -- "$relative_path")
done < "$list" > "$output"
