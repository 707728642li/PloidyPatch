#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then echo "usage: $0 PROJECT_ROOT" >&2; exit 2; fi
project_root=$(realpath "$1")
input_root=$project_root/data/derived/holdout_evaluator/cotton_wgdi_inputs_v0.1
synteny_env=$project_root/envs/ploidypatch-synteny
diamond_bin=$synteny_env/bin/diamond
wgdi_bin=$synteny_env/bin/wgdi
parallel_bin=/data/codexli/software/conda/miniforge3/bin/parallel
result_root=$project_root/results/evaluator/cotton_holdout_v0.1/wgdi
working_root=${result_root}.working

for required in "$input_root/SHA256SUMS" "$diamond_bin" "$wgdi_bin" "$parallel_bin"; do
    if [[ ! -s $required ]]; then echo "missing cotton WGDI prerequisite: $required" >&2; exit 1; fi
done
(cd "$input_root" && sha256sum -c SHA256SUMS >/dev/null)
if [[ -e $result_root || -e $working_root ]]; then echo "refusing to overwrite cotton WGDI" >&2; exit 1; fi
mkdir -p "$working_root"/{db,blast,config,collinearity,logs}
{
    printf 'field\tvalue\ncode_commit\t%s\n' "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'access\tevaluator_only\nquery\tGossypium_hirsutum_AD\n'
    printf 'references\tGossypium_arboreum_A,Gossypium_raimondii_D\n'
    printf 'diamond_threads_per_job\t64\nwgdi_processes_per_job\t64\n'
} > "$working_root/run_contract.tsv"

make_db() {
    local prefix=$1
    /usr/bin/time -v -o "$working_root/logs/$prefix.makedb.time.txt" \
        "$diamond_bin" makedb --in "$input_root/$prefix/$prefix.wgdi.pep.fa" \
            --db "$working_root/db/$prefix" \
            > "$working_root/logs/$prefix.makedb.stdout.log" \
            2> "$working_root/logs/$prefix.makedb.stderr.log"
}
export -f make_db
export working_root input_root diamond_bin
printf '%s\n' gar_a gra_d | "$parallel_bin" --jobs 2 make_db {}

run_diamond() {
    local pair=$1 ref=$2
    /usr/bin/time -v -o "$working_root/logs/$pair.diamond.time.txt" \
        "$diamond_bin" blastp --query "$input_root/ghi_ad/ghi_ad.wgdi.pep.fa" \
            --db "$working_root/db/$ref" --out "$working_root/blast/$pair.tsv" \
            --outfmt 6 --evalue 1e-5 --max-target-seqs 20 --more-sensitive --threads 64 \
            > "$working_root/logs/$pair.diamond.stdout.log" \
            2> "$working_root/logs/$pair.diamond.stderr.log"
}
export -f run_diamond
printf '%s\t%s\n' ghi_ad_vs_gar_a gar_a ghi_ad_vs_gra_d gra_d \
    | "$parallel_bin" --jobs 2 --delay 1 --colsep '\t' run_diamond {1} {2}

for entry in ghi_ad_vs_gar_a:gar_a ghi_ad_vs_gra_d:gra_d; do
    pair=${entry%%:*}; ref=${entry#*:}; config=$working_root/config/$pair.conf
    {
        printf '[collinearity]\n'
        printf 'gff1 = %s\ngff2 = %s\n' "$input_root/ghi_ad/ghi_ad.wgdi.gff" "$input_root/$ref/$ref.wgdi.gff"
        printf 'lens1 = %s\nlens2 = %s\n' "$input_root/ghi_ad/ghi_ad.wgdi.lens" "$input_root/$ref/$ref.wgdi.lens"
        printf 'blast = %s\n' "$working_root/blast/$pair.tsv"
        printf 'blast_reverse = false\ncomparison = genomes\nmultiple = 2\nprocess = 64\n'
        printf 'evalue = 1e-5\nscore = 100\ngrading = 50,40,25\nmg = 40,40\n'
        printf 'pvalue = 0.2\nrepeat_number = 20\nposition = order\n'
        printf 'savefile = %s\n' "$working_root/collinearity/$pair.tsv"
    } > "$config"
done
run_wgdi() {
    local pair=$1
    /usr/bin/time -v -o "$working_root/logs/$pair.wgdi.time.txt" \
        "$wgdi_bin" -icl "$working_root/config/$pair.conf" \
        > "$working_root/logs/$pair.wgdi.stdout.log" \
        2> "$working_root/logs/$pair.wgdi.stderr.log"
}
export -f run_wgdi
export wgdi_bin
printf '%s\n' ghi_ad_vs_gar_a ghi_ad_vs_gra_d \
    | "$parallel_bin" --jobs 2 --delay 1 run_wgdi {}

for pair in ghi_ad_vs_gar_a ghi_ad_vs_gra_d; do
    if [[ ! -s $working_root/collinearity/$pair.tsv ]]; then echo "missing WGDI result: $pair" >&2; exit 1; fi
done
(
    cd "$working_root"
    find . -type f \( -name '*.tsv' -o -name '*.conf' -o -name '*.json' \) \
        -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS
    sha256sum -c SHA256SUMS >/dev/null
)
du -sb "$working_root" > "$working_root/disk_bytes.txt"
mv "$working_root" "$result_root"
printf 'cotton evaluator WGDI frozen: %s\n' "$result_root"
