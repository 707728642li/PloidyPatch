#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 STAGED_DATA_ROOT OUTPUT_ROOT" >&2
    exit 2
fi

data_root=$(realpath "$1")
output_root=$(realpath -m "$2")
run_id=$(date -u +%Y%m%dT%H%M%SZ)
run_dir="$output_root/$run_id"
manifest="$run_dir/audit_manifest.tsv"

mkdir -p "$run_dir"
printf 'dataset_id\treport_path\treport_sha256\tgrade\texit_code\n' > "$manifest"

resolve_one() {
    local search_dir=$1
    local pattern=$2
    local -a matches=()
    mapfile -t matches < <(find "$search_dir" -maxdepth 1 -type f -name "$pattern" | sort)
    if [[ ${#matches[@]} -ne 1 ]]; then
        echo "Expected one match for $search_dir/$pattern; found ${#matches[@]}" >&2
        return 2
    fi
    printf '%s\n' "${matches[0]}"
}

for dataset_dir in "$data_root"/*; do
    [[ -d "$dataset_dir" ]] || continue
    dataset_id=$(basename "$dataset_dir")
    gff=$(resolve_one "$dataset_dir/annotation" '*.gff3.gz')
    protein=$(resolve_one "$dataset_dir/annotation" '*.pep.all.fa.gz')
    cds=$(resolve_one "$dataset_dir/annotation" '*.cds.all.fa.gz')
    fai=$(resolve_one "$dataset_dir/genome" '*.fai')
    report="$run_dir/$dataset_id.json"

    exit_code=0
    ploidypatch audit \
        --gff "$gff" \
        --protein "$protein" \
        --cds "$cds" \
        --fai "$fai" \
        --checksums \
        --output "$report" || exit_code=$?

    if [[ ! -s "$report" ]]; then
        echo "Audit did not create a report for $dataset_id" >&2
        exit 1
    fi
    grade=$(python -c 'import json,sys; print(json.load(open(sys.argv[1]))["quality_gate"]["grade"])' "$report")
    report_sha=$(sha256sum "$report" | cut -d ' ' -f 1)
    printf '%s\t%s\t%s\t%s\t%s\n' \
        "$dataset_id" "$report" "$report_sha" "$grade" "$exit_code" >> "$manifest"
done

printf 'Audit complete. Manifest: %s\n' "$manifest"
