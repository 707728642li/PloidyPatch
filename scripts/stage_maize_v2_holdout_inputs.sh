#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then echo "usage: $0 PROJECT_ROOT" >&2; exit 2; fi
project_root=$(realpath -m "$1")
case "$project_root" in
    /data/codexli/projects/*) ;;
    *) echo "refusing non-standard project root: $project_root" >&2; exit 2 ;;
esac

code_root=$project_root/code
source_contract=$code_root/config/maize_v2_source_files.tsv
destination_root=$project_root/data/public/maize_v2_holdout
manifest_root=$project_root/logs/staging
code_commit=${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}
script_sha=$(sha256sum "$0" | awk '{print $1}')
run_id=$(date -u +%Y%m%dT%H%M%SZ)
manifest=$manifest_root/maize_v2_holdout_${run_id}.tsv
manifest_partial=${manifest}.partial.$$
copy_partial=""

[[ -s $source_contract ]] || { echo "missing maize source contract" >&2; exit 1; }
mkdir -p "$destination_root" "$manifest_root"
cleanup() {
    [[ -z $copy_partial ]] || rm -f -- "$copy_partial"
    rm -f -- "$manifest_partial"
}
trap cleanup EXIT
printf 'species\trole\tcode_commit\tstage_script_sha256\tsource_path\tdestination_path\tbytes\tsha256\tgzip_test\tstatus\n' > "$manifest_partial"

while IFS=$'\t' read -r species role expected_bytes expected_sha source_path destination_name; do
    [[ -n $species && -n $role && -n $source_path && -n $destination_name ]] || {
        echo "malformed maize source contract row" >&2; exit 1;
    }
    case "$source_path" in
        /nas_data/NFS/Public_genome_data/*) ;;
        *) echo "source outside declared public tree: $source_path" >&2; exit 1 ;;
    esac
    [[ -s $source_path ]] || { echo "missing public input: $source_path" >&2; exit 1; }
    observed_bytes=$(stat -Lc %s "$source_path")
    [[ $observed_bytes == "$expected_bytes" ]] || {
        echo "source byte count changed: $source_path" >&2; exit 1;
    }
    observed_sha=$(sha256sum "$source_path" | awk '{print $1}')
    [[ $observed_sha == "$expected_sha" ]] || {
        echo "source checksum changed: $source_path" >&2; exit 1;
    }
    gzip -t "$source_path"
    destination_dir=$destination_root/$species
    destination_path=$destination_dir/$destination_name
    mkdir -p "$destination_dir"
    status=copied_verified
    if [[ -e $destination_path ]]; then
        destination_sha=$(sha256sum "$destination_path" | awk '{print $1}')
        [[ $destination_sha == "$expected_sha" ]] || {
            echo "existing staged file differs: $destination_path" >&2; exit 1;
        }
        status=already_verified
    else
        copy_partial=${destination_path}.partial.$$
        cp -p "$source_path" "$copy_partial"
        destination_sha=$(sha256sum "$copy_partial" | awk '{print $1}')
        [[ $destination_sha == "$expected_sha" ]] || {
            echo "checksum mismatch after copy: $destination_path" >&2; exit 1;
        }
        mv -- "$copy_partial" "$destination_path"
        copy_partial=""
    fi
    gzip -t "$destination_path"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tpass\t%s\n' \
        "$species" "$role" "$code_commit" "$script_sha" "$source_path" \
        "$destination_path" "$expected_bytes" "$expected_sha" "$status" \
        >> "$manifest_partial"
done < <(tail -n +2 "$source_contract")

rows=$(tail -n +2 "$manifest_partial" | wc -l)
[[ $rows -eq 9 ]] || { echo "expected nine staged inputs; observed $rows" >&2; exit 1; }
mv -- "$manifest_partial" "$manifest"
trap - EXIT
printf 'maize v0.2 holdout inputs staged: %s\n' "$manifest"
