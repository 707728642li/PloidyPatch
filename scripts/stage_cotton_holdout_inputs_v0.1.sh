#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "usage: $0 PUBLIC_ROOT PROJECT_ROOT" >&2
    exit 2
fi

public_root=$(realpath -m "$1")
project_root=$(realpath -m "$2")
case "$project_root" in
    /data/codexli/projects/*) ;;
    *)
        echo "refusing non-standard project root: $project_root" >&2
        exit 2
        ;;
esac

source_root=$public_root/PGCP/data
destination_root=$project_root/data/public/pgcp_cotton_holdout_v0.1
manifest_dir=$project_root/logs/staging
code_commit=${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}
script_sha=$(sha256sum "$0" | awk '{print $1}')
run_id=$(date -u +%Y%m%dT%H%M%SZ)
manifest=$manifest_dir/pgcp_cotton_holdout_v0.1_${run_id}.tsv
manifest_partial=${manifest}.partial.$$
copy_partial=""

cleanup() {
    if [[ -n $copy_partial ]]; then
        rm -f -- "$copy_partial"
    fi
    rm -f -- "$manifest_partial"
}
trap cleanup EXIT

datasets=(
    gossypium_hirsutum
    gossypium_arboreum
    gossypium_raimondii
    gossypium_barbadense
)

mkdir -p "$destination_root" "$manifest_dir"
printf 'dataset_id\tcode_commit\tstage_script_sha256\trole\tsource_path\tdestination_path\tbytes\tsha256\tgzip_test\tstatus\n' \
    > "$manifest_partial"

for dataset_id in "${datasets[@]}"; do
    source_dir=$source_root/$dataset_id
    destination_dir=$destination_root/$dataset_id
    mkdir -p "$destination_dir"
    for role in genome gff3; do
        case "$role" in
            genome) suffix=genomic.fa.gz ;;
            gff3) suffix=genomic.gff.gz ;;
        esac
        source_path=$source_dir/$dataset_id.$suffix
        destination_path=$destination_dir/$dataset_id.$suffix
        if [[ ! -s $source_path ]]; then
            echo "missing or empty public cotton input: $source_path" >&2
            exit 1
        fi
        source_sha=$(sha256sum "$source_path" | awk '{print $1}')
        status=copied_verified
        if [[ -e $destination_path ]]; then
            destination_sha=$(sha256sum "$destination_path" | awk '{print $1}')
            if [[ $source_sha != "$destination_sha" ]]; then
                echo "existing staged input differs from source: $destination_path" >&2
                exit 1
            fi
            status=already_verified
        else
            copy_partial=${destination_path}.partial.$$
            cp -p "$source_path" "$copy_partial"
            destination_sha=$(sha256sum "$copy_partial" | awk '{print $1}')
            if [[ $source_sha != "$destination_sha" ]]; then
                echo "checksum mismatch after staging: $copy_partial" >&2
                exit 1
            fi
            mv -- "$copy_partial" "$destination_path"
            copy_partial=""
        fi
        if ! gzip -t "$destination_path"; then
            echo "gzip integrity test failed: $destination_path" >&2
            exit 1
        fi
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tpass\t%s\n' \
            "$dataset_id" "$code_commit" "$script_sha" "$role" \
            "$source_path" "$destination_path" \
            "$(stat -Lc %s "$destination_path")" "$source_sha" "$status" \
            >> "$manifest_partial"
    done
done

mv -- "$manifest_partial" "$manifest"
trap - EXIT
printf 'cotton holdout inputs staged: %s\n' "$manifest"
