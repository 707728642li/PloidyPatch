#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "usage: $0 PROJECT_ROOT BLIND_BENCHMARK_ROOT" >&2
    exit 2
fi

project_root=$1
benchmark_root=$2
query_dir="$benchmark_root/blind/evidence/wgdi_inputs/bna_daae_blind"
target_root="$project_root/data/derived/wgdi_inputs/v0.1"
database_root="$project_root/results/evidence/wgdi/brassica_v0.1/db"
output_root="$benchmark_root/blind/evidence/wgdi"
log_root="$benchmark_root/logs/blind_wgdi_v0.1"
synteny_env="$project_root/envs/ploidypatch-synteny"
pair_table="$project_root/code/config/wgdi_brassica_blind_pairs_v0.1.tsv"
parallel_bin="$(conda info --base)/bin/parallel"

if [[ -e "$output_root" ]]; then
    echo "refusing to reuse blind WGDI output directory: $output_root" >&2
    exit 1
fi
test -s "$query_dir/bna_daae_blind.wgdi.pep.fa"
test -s "$query_dir/bna_daae_blind.wgdi.gff"
test -s "$query_dir/bna_daae_blind.wgdi.lens"
test -x "$parallel_bin"
mkdir -p "$output_root/blast" "$output_root/config" \
    "$output_root/collinearity" "$log_root"

run_similarity() {
    local pair_id=$1
    local target_prefix=$2
    local output="$output_root/blast/$pair_id.tsv"
    local time_log="$log_root/$pair_id.diamond.time.txt"
    test -s "$database_root/$target_prefix.dmnd"
    /usr/bin/time -v \
        conda run -p "$synteny_env" --no-capture-output \
        diamond blastp \
        --query "$query_dir/bna_daae_blind.wgdi.pep.fa" \
        --db "$database_root/$target_prefix" \
        --out "$output" \
        --outfmt 6 \
        --evalue 1e-5 \
        --max-target-seqs 20 \
        --more-sensitive \
        --threads 48 \
        > "$log_root/$pair_id.diamond.stdout.log" \
        2> "$time_log"
}

write_config() {
    local pair_id=$1
    local target_prefix=$2
    local multiple=$3
    local config="$output_root/config/$pair_id.collinearity.conf"
    {
        echo '[collinearity]'
        echo "gff1 = $query_dir/bna_daae_blind.wgdi.gff"
        echo "gff2 = $target_root/$target_prefix/$target_prefix.wgdi.gff"
        echo "lens1 = $query_dir/bna_daae_blind.wgdi.lens"
        echo "lens2 = $target_root/$target_prefix/$target_prefix.wgdi.lens"
        echo "blast = $output_root/blast/$pair_id.tsv"
        echo 'blast_reverse = false'
        echo 'comparison = genomes'
        echo "multiple = $multiple"
        echo 'process = 32'
        echo 'evalue = 1e-5'
        echo 'score = 100'
        echo 'grading = 50,40,25'
        echo 'mg = 40,40'
        echo 'pvalue = 0.2'
        echo 'repeat_number = 20'
        echo 'position = order'
        echo "savefile = $output_root/collinearity/$pair_id.collinearity.tsv"
    } > "$config"
}

run_collinearity() {
    local pair_id=$1
    /usr/bin/time -v \
        conda run -p "$synteny_env" --no-capture-output \
        wgdi -icl "$output_root/config/$pair_id.collinearity.conf" \
        > "$log_root/$pair_id.wgdi.stdout.log" \
        2> "$log_root/$pair_id.wgdi.time.txt"
}

export -f run_similarity write_config run_collinearity
export query_dir target_root database_root output_root log_root synteny_env

tail -n +2 "$pair_table" \
    | "$parallel_bin" --jobs 2 --delay 1 --colsep '\t' \
        run_similarity {1} {2}
tail -n +2 "$pair_table" \
    | "$parallel_bin" --jobs 2 --colsep '\t' write_config {1} {2} {3}
tail -n +2 "$pair_table" \
    | "$parallel_bin" --jobs 2 --delay 1 --colsep '\t' run_collinearity {1}

{
    echo -e 'artifact\tbytes\tsha256'
    find "$output_root" -type f -print0 \
        | sort -z \
        | while IFS= read -r -d '' file; do
            bytes=$(stat -c '%s' "$file")
            checksum=$(sha256sum "$file" | awk '{print $1}')
            echo -e "${file#"$output_root/"}\t$bytes\t$checksum"
        done
} > "$output_root/run_manifest.tsv"
