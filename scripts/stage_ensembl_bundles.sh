#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 3 ]]; then
    echo "Usage: $0 PUBLIC_ROOT PROJECT_ROOT DATASET_SLUG [DATASET_SLUG ...]" >&2
    exit 2
fi

public_root=$1
project_root=$2
shift 2

resolved_project_root=$(realpath -m "$project_root")
case "$resolved_project_root" in
    /data/codexli/projects/*) ;;
    *)
        echo "Refusing non-standard project root: $resolved_project_root" >&2
        exit 2
        ;;
esac

genome_root="$public_root/genome_database"
gff_root="$public_root/Ensemble/gff3"
staging_root="$resolved_project_root/data/public/ensembl_plants_62"
manifest_dir="$resolved_project_root/logs/staging"
run_id=$(date -u +%Y%m%dT%H%M%SZ)
manifest="$manifest_dir/ensembl_plants_62_${run_id}.tsv"

mkdir -p "$staging_root" "$manifest_dir"
printf 'dataset_id\trole\tsource_path\tdestination_path\tbytes\tsha256\tstatus\n' > "$manifest"

stage_unique() {
    local dataset_id=$1
    local role=$2
    local search_dir=$3
    local pattern=$4
    local destination_dir=$5
    local -a matches=()

    mapfile -t matches < <(find "$search_dir" -maxdepth 1 -type f -name "$pattern" | sort)
    if [[ ${#matches[@]} -ne 1 ]]; then
        echo "Expected one $role file for $dataset_id; found ${#matches[@]}" >&2
        exit 1
    fi

    local source_path=${matches[0]}
    local destination_path="$destination_dir/$(basename "$source_path")"
    mkdir -p "$destination_dir"

    if [[ -e "$destination_path" ]]; then
        local source_sha destination_sha
        source_sha=$(sha256sum "$source_path" | awk '{print $1}')
        destination_sha=$(sha256sum "$destination_path" | awk '{print $1}')
        if [[ "$source_sha" != "$destination_sha" ]]; then
            echo "Existing destination differs from source: $destination_path" >&2
            exit 1
        fi
        printf '%s\t%s\t%s\t%s\t%s\t%s\talready_verified\n' \
            "$dataset_id" "$role" "$source_path" "$destination_path" \
            "$(stat -c %s "$source_path")" "$source_sha" >> "$manifest"
        return
    fi

    cp -p "$source_path" "$destination_path"
    local source_sha destination_sha
    source_sha=$(sha256sum "$source_path" | awk '{print $1}')
    destination_sha=$(sha256sum "$destination_path" | awk '{print $1}')
    if [[ "$source_sha" != "$destination_sha" ]]; then
        echo "Checksum mismatch after copy: $destination_path" >&2
        exit 1
    fi
    printf '%s\t%s\t%s\t%s\t%s\t%s\tcopied_verified\n' \
        "$dataset_id" "$role" "$source_path" "$destination_path" \
        "$(stat -c %s "$source_path")" "$source_sha" >> "$manifest"
}

for dataset_id in "$@"; do
    dataset_destination="$staging_root/$dataset_id"
    stage_unique "$dataset_id" genome \
        "$genome_root/$dataset_id/dna_index" '*.dna.toplevel.fa.gz' \
        "$dataset_destination/genome"
    stage_unique "$dataset_id" fai \
        "$genome_root/$dataset_id/dna_index" '*.dna.toplevel.fa.gz.fai' \
        "$dataset_destination/genome"
    stage_unique "$dataset_id" gzi \
        "$genome_root/$dataset_id/dna_index" '*.dna.toplevel.fa.gz.gzi' \
        "$dataset_destination/genome"
    stage_unique "$dataset_id" gff3 \
        "$gff_root/$dataset_id" '*.62.gff3.gz' \
        "$dataset_destination/annotation"
    stage_unique "$dataset_id" cds \
        "$genome_root/$dataset_id/cds" '*.cds.all.fa.gz' \
        "$dataset_destination/annotation"
    stage_unique "$dataset_id" protein \
        "$genome_root/$dataset_id/pep" '*.pep.all.fa.gz' \
        "$dataset_destination/annotation"
done

printf 'Staging complete. Manifest: %s\n' "$manifest"
