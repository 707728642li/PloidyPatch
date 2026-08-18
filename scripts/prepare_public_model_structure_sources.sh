#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "usage: $0 PROJECT_ROOT" >&2
    exit 2
fi

project_root=$(realpath "$1")
public_root=$project_root/data/public/ensembl_plants_62
result_root=$project_root/data/derived/structure_sources_v0.1
working_root=${result_root}.working
if [[ -e $result_root || -e $working_root ]]; then
    echo "refusing to overwrite public-model structure sources: $result_root" >&2
    exit 1
fi

declare -A sources=(
    [ath_tair10]="$public_root/arabidopsis_thaliana/annotation/Arabidopsis_thaliana.TAIR10.62.gff3.gz"
    [osa_irgsp10]="$public_root/oryza_sativa/annotation/Oryza_sativa.IRGSP-1.0.62.gff3.gz"
)
for source in "${sources[@]}"; do
    if [[ ! -s $source ]]; then
        echo "missing staged public-model GFF3: $source" >&2
        exit 1
    fi
    gzip -t "$source"
done

mkdir -p "$working_root"
{
    printf 'field\tvalue\n'
    printf 'code_commit\t%s\n' "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'source_collection\tEnsembl_Plants_62\n'
    printf 'transformation\tgzip_decompression_only\n'
    printf 'source_data_modified\tfalse\n'
} > "$working_root/run_contract.tsv"
printf 'dataset_id\trole\tbytes\tsha256\tpath\n' > "$working_root/file_manifest.tsv"
for dataset_id in ath_tair10 osa_irgsp10; do
    source=${sources[$dataset_id]}
    dataset_root=$working_root/$dataset_id
    output=$dataset_root/source.gff3
    mkdir -p "$dataset_root"
    gzip -cd -- "$source" > "$output"
    if [[ ! -s $output ]]; then
        echo "decompressed public-model GFF3 is empty: $dataset_id" >&2
        exit 1
    fi
    printf '%s\tcompressed_source\t%s\t%s\t%s\n' \
        "$dataset_id" "$(stat -Lc %s "$source")" \
        "$(sha256sum "$source" | awk '{print $1}')" "$source" \
        >> "$working_root/file_manifest.tsv"
    printf '%s\tdecompressed_source\t%s\t%s\t%s\n' \
        "$dataset_id" "$(stat -Lc %s "$output")" \
        "$(sha256sum "$output" | awk '{print $1}')" \
        "$result_root/$dataset_id/source.gff3" \
        >> "$working_root/file_manifest.tsv"
done
(
    cd "$working_root"
    find . -type f \( -name '*.gff3' -o -name '*.tsv' \) -print0 \
        | sort -z | xargs -0 sha256sum > SHA256SUMS
)
du -sb "$working_root" > "$working_root/disk_bytes.txt"
mv "$working_root" "$result_root"
printf 'public-model structure sources prepared: %s\n' "$result_root"
