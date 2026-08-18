#!/usr/bin/env bash
set -euo pipefail

project_root=${1:?"usage: prepare_brassica_wgdi_inputs.sh PROJECT_ROOT"}
data_root="$project_root/data/public/ensembl_plants_62"
output_root="$project_root/data/derived/wgdi_inputs/v0.1"
log_root="$project_root/logs/wgdi_inputs/v0.1"
dev_env="$project_root/envs/ploidypatch-dev"
parallel_bin=$(command -v parallel || true)
if [[ -z "$parallel_bin" ]]; then
    parallel_bin="$(conda info --base)/bin/parallel"
fi
if [[ ! -x "$parallel_bin" ]]; then
    echo "GNU parallel is required but was not found" >&2
    exit 1
fi

mkdir -p "$output_root" "$log_root"

prepare_one() {
    shopt -s nullglob
    local dataset=$1
    local prefix=$2
    local dataset_root="$data_root/$dataset"
    local output_dir="$output_root/$prefix"
    local gff_files=("$dataset_root"/annotation/*.gff3.gz)
    local protein_files=("$dataset_root"/annotation/*.pep.all.fa.gz)
    local fai_files=("$dataset_root"/genome/*.fa.gz.fai)

    if [[ ${#gff_files[@]} -ne 1 || ${#protein_files[@]} -ne 1 || ${#fai_files[@]} -ne 1 ]]; then
        echo "Expected exactly one GFF3, protein FASTA, and FAI for $dataset" >&2
        return 1
    fi
    if [[ -e "$output_dir" ]]; then
        echo "Refusing to reuse existing output directory: $output_dir" >&2
        return 1
    fi

    /usr/bin/time -v \
        conda run -p "$dev_env" --no-capture-output \
        ploidypatch evidence prepare-wgdi \
        --gff "${gff_files[0]}" \
        --protein "${protein_files[0]}" \
        --fai "${fai_files[0]}" \
        --output-dir "$output_dir" \
        --prefix "$prefix" \
        --min-genes-per-seqid 100 \
        > "$log_root/$prefix.stdout.json" \
        2> "$log_root/$prefix.time.txt"
}

export -f prepare_one
export data_root output_root log_root dev_env

printf '%s\t%s\n' \
    brassica_rapa bra_a \
    brassica_oleracea bol_c \
    brassica_napus bna_ac \
    | "$parallel_bin" --jobs 3 --delay 0.5 --colsep '\t' prepare_one {1} {2}
