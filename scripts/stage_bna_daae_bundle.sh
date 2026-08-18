#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 PUBLIC_ROOT PROJECT_ROOT" >&2
    exit 2
fi

public_root=$(realpath -m "$1")
project_root=$(realpath -m "$2")
case "$project_root" in
    /data/codexli/projects/*) ;;
    *)
        echo "Refusing non-standard project root: $project_root" >&2
        exit 2
        ;;
esac

source_root="$public_root/genome_database/zhangshulin_db/clean_data"
destination_root="$project_root/data/public/refseq_bna_daae_gcf0203794851"
manifest_dir="$project_root/logs/staging"
run_id=$(date -u +%Y%m%dT%H%M%SZ)
manifest="$manifest_dir/refseq_bna_daae_gcf0203794851_${run_id}.tsv"
manifest_partial="$manifest.partial.$$"
copy_partial=""

cleanup() {
    [[ -z "$copy_partial" ]] || rm -f -- "$copy_partial"
    rm -f -- "$manifest_partial"
}
trap cleanup EXIT

declare -a roles=(genome fai gff3 protein cds)
declare -a source_paths=(
    "$source_root/fa_dir/Brassica_napus.fa"
    "$source_root/fa_dir/Brassica_napus.fa.fai"
    "$source_root/gff_dir/Brassica_napus.gff3"
    "$source_root/cds_prot_dir/Brassica_napus.prot.fa"
    "$source_root/cds_prot_dir/Brassica_napus.cds.fa"
)
declare -a relative_paths=(
    "genome/Brassica_napus.fa"
    "genome/Brassica_napus.fa.fai"
    "annotation/Brassica_napus.gff3"
    "annotation/Brassica_napus.prot.fa"
    "annotation/Brassica_napus.cds.fa"
)

for source_path in "${source_paths[@]}"; do
    if [[ ! -f "$source_path" ]]; then
        echo "Required public file is missing: $source_path" >&2
        exit 1
    fi
done

mkdir -p "$destination_root/genome" "$destination_root/annotation" "$manifest_dir"
printf 'dataset_id\tassembly_accession\trole\tsource_path\tdestination_path\tbytes\tsha256\tstatus\n' > "$manifest_partial"

for index in "${!source_paths[@]}"; do
    role=${roles[$index]}
    source_path=${source_paths[$index]}
    destination_path="$destination_root/${relative_paths[$index]}"
    source_sha=$(sha256sum "$source_path" | awk '{print $1}')
    status=copied_verified
    if [[ -e "$destination_path" ]]; then
        destination_sha=$(sha256sum "$destination_path" | awk '{print $1}')
        if [[ "$source_sha" != "$destination_sha" ]]; then
            echo "Existing destination differs from source: $destination_path" >&2
            exit 1
        fi
        status=already_verified
    else
        copy_partial="$destination_path.partial.$$"
        cp -p "$source_path" "$copy_partial"
        destination_sha=$(sha256sum "$copy_partial" | awk '{print $1}')
        if [[ "$source_sha" != "$destination_sha" ]]; then
            echo "Checksum mismatch after copy: $copy_partial" >&2
            exit 1
        fi
        mv -- "$copy_partial" "$destination_path"
        copy_partial=""
    fi
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        bna_daae GCF_020379485.1 "$role" "$source_path" "$destination_path" \
        "$(stat -c %s "$source_path")" "$source_sha" "$status" >> "$manifest_partial"
done

mv -- "$manifest_partial" "$manifest"
printf 'Staging complete. Manifest: %s\n' "$manifest"
