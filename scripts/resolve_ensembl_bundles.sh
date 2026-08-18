#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 2 ]]; then
    echo "Usage: $0 PUBLIC_ROOT DATASET_SLUG [DATASET_SLUG ...]" >&2
    exit 2
fi

public_root=$1
shift

genome_root="$public_root/genome_database"
gff_root="$public_root/Ensemble/gff3"

printf 'dataset_id\trole\tsource_path\tbytes\tstatus\n'

emit_match() {
    local dataset_id=$1
    local role=$2
    local search_dir=$3
    local pattern=$4
    local -a matches=()

    if [[ -d "$search_dir" ]]; then
        mapfile -t matches < <(find "$search_dir" -maxdepth 1 -type f -name "$pattern" | sort)
    fi
    if [[ ${#matches[@]} -eq 1 ]]; then
        printf '%s\t%s\t%s\t%s\tok\n' \
            "$dataset_id" "$role" "${matches[0]}" "$(stat -c %s "${matches[0]}")"
    elif [[ ${#matches[@]} -eq 0 ]]; then
        printf '%s\t%s\tNA\t0\tmissing\n' "$dataset_id" "$role"
    else
        printf '%s\t%s\tNA\t0\tambiguous_%d_matches\n' \
            "$dataset_id" "$role" "${#matches[@]}"
    fi
}

for dataset_id in "$@"; do
    emit_match "$dataset_id" genome \
        "$genome_root/$dataset_id/dna_index" '*.dna.toplevel.fa.gz'
    emit_match "$dataset_id" fai \
        "$genome_root/$dataset_id/dna_index" '*.dna.toplevel.fa.gz.fai'
    emit_match "$dataset_id" gzi \
        "$genome_root/$dataset_id/dna_index" '*.dna.toplevel.fa.gz.gzi'
    emit_match "$dataset_id" gff3 \
        "$gff_root/$dataset_id" '*.62.gff3.gz'
    emit_match "$dataset_id" cds \
        "$genome_root/$dataset_id/cds" '*.cds.all.fa.gz'
    emit_match "$dataset_id" protein \
        "$genome_root/$dataset_id/pep" '*.pep.all.fa.gz'
done
