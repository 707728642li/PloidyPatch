#!/usr/bin/env bash
set -euo pipefail

project_root=${1:?"usage: run_brassica_wgdi.sh PROJECT_ROOT [fresh|resume]"}
mode=${2:-fresh}
input_root="$project_root/data/derived/wgdi_inputs/v0.1"
output_root="$project_root/results/evidence/wgdi/brassica_v0.1"
log_root="$project_root/logs/wgdi/brassica_v0.1"
synteny_env="$project_root/envs/ploidypatch-synteny"
pair_table="$project_root/code/config/wgdi_brassica_pairs_v0.1.tsv"
parallel_bin=$(command -v parallel || true)
if [[ -z "$parallel_bin" ]]; then
    parallel_bin="$(conda info --base)/bin/parallel"
fi
if [[ ! -x "$parallel_bin" ]]; then
    echo "GNU parallel is required but was not found" >&2
    exit 1
fi
if [[ "$mode" != fresh && "$mode" != resume ]]; then
    echo "Mode must be fresh or resume" >&2
    exit 1
fi
if [[ "$mode" == fresh && -e "$output_root" ]]; then
    echo "Refusing to reuse existing WGDI output directory: $output_root" >&2
    exit 1
fi
if [[ "$mode" == resume && ! -d "$output_root" ]]; then
    echo "Cannot resume missing WGDI output directory: $output_root" >&2
    exit 1
fi

mkdir -p "$output_root/db" "$output_root/blast" "$output_root/collinearity" \
    "$output_root/config" "$log_root"

is_complete() {
    local artifact=$1
    local time_log=$2
    [[ -s "$artifact" ]] && [[ -s "$time_log" ]] \
        && grep -q 'Exit status: 0' "$time_log"
}

make_database() {
    local prefix=$1
    local protein="$input_root/$prefix/$prefix.wgdi.pep.fa"
    local database="$output_root/db/$prefix.dmnd"
    local time_log="$log_root/$prefix.makedb.time.txt"
    test -s "$protein"
    if is_complete "$database" "$time_log"; then
        echo "SKIP complete database: $prefix" >&2
        return 0
    fi
    /usr/bin/time -v \
        conda run -p "$synteny_env" --no-capture-output \
        diamond makedb \
        --in "$protein" \
        --db "$output_root/db/$prefix" \
        > "$log_root/$prefix.makedb.stdout.log" \
        2> "$time_log"
}

export -f is_complete make_database
export input_root output_root log_root synteny_env
printf '%s\n' bra_a bol_c \
    | "$parallel_bin" --jobs 2 --delay 0.5 make_database {}

run_similarity() {
    local pair_id=$1
    local query_prefix=$2
    local target_prefix=$3
    local query="$input_root/$query_prefix/$query_prefix.wgdi.pep.fa"
    local database="$output_root/db/$target_prefix"
    local output="$output_root/blast/$pair_id.tsv"
    local time_log="$log_root/$pair_id.diamond.time.txt"
    test -s "$query"
    test -s "$database.dmnd"
    if is_complete "$output" "$time_log"; then
        echo "SKIP complete similarity: $pair_id" >&2
        return 0
    fi
    /usr/bin/time -v \
        conda run -p "$synteny_env" --no-capture-output \
        diamond blastp \
        --query "$query" \
        --db "$database" \
        --out "$output" \
        --outfmt 6 \
        --evalue 1e-5 \
        --max-target-seqs 20 \
        --more-sensitive \
        --threads 48 \
        > "$log_root/$pair_id.diamond.stdout.log" \
        2> "$time_log"
}

export -f run_similarity
tail -n +2 "$pair_table" \
    | "$parallel_bin" --jobs 2 --delay 1 --colsep '\t' \
        run_similarity {1} {2} {3}

write_config() {
    local pair_id=$1
    local query_prefix=$2
    local target_prefix=$3
    local multiple=$4
    local config="$output_root/config/$pair_id.collinearity.conf"
    {
        echo '[collinearity]'
        echo "gff1 = $input_root/$query_prefix/$query_prefix.wgdi.gff"
        echo "gff2 = $input_root/$target_prefix/$target_prefix.wgdi.gff"
        echo "lens1 = $input_root/$query_prefix/$query_prefix.wgdi.lens"
        echo "lens2 = $input_root/$target_prefix/$target_prefix.wgdi.lens"
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
    local config="$output_root/config/$pair_id.collinearity.conf"
    local output="$output_root/collinearity/$pair_id.collinearity.tsv"
    local time_log="$log_root/$pair_id.wgdi.time.txt"
    test -s "$config"
    if is_complete "$output" "$time_log"; then
        echo "SKIP complete collinearity: $pair_id" >&2
        return 0
    fi
    /usr/bin/time -v \
        conda run -p "$synteny_env" --no-capture-output \
        wgdi -icl "$config" \
        > "$log_root/$pair_id.wgdi.stdout.log" \
        2> "$time_log"
}

export -f write_config run_collinearity
export pair_table
tail -n +2 "$pair_table" \
    | "$parallel_bin" --jobs 2 --colsep '\t' write_config {1} {2} {3} {4}
tail -n +2 "$pair_table" \
    | "$parallel_bin" --jobs 2 --delay 1 --colsep '\t' run_collinearity {1}

{
    echo -e 'artifact\tbytes\tsha256'
    find "$output_root/blast" "$output_root/collinearity" "$output_root/config" \
        -type f -print0 \
        | sort -z \
        | while IFS= read -r -d '' file; do
            bytes=$(stat -c '%s' "$file")
            checksum=$(sha256sum "$file" | awk '{print $1}')
            echo -e "${file#"$output_root/"}\t$bytes\t$checksum"
        done
} > "$output_root/run_manifest.tsv"

conda run -p "$synteny_env" --no-capture-output wgdi --version \
    > "$output_root/wgdi.version.txt"
conda run -p "$synteny_env" --no-capture-output diamond version \
    > "$output_root/diamond.version.txt"
