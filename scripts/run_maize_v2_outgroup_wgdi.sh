#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then echo "usage: $0 PROJECT_ROOT" >&2; exit 2; fi
project_root=$(realpath "$1")
input_root=$project_root/data/derived/holdout_evaluator/maize_v2_wgdi_inputs
synteny_env=$project_root/envs/ploidypatch-synteny
diamond_bin=$synteny_env/bin/diamond
wgdi_bin=$synteny_env/bin/wgdi
parallel_bin=/data/codexli/software/conda/miniforge3/bin/parallel
policy=$project_root/code/config/maize_v2_zero_retuning_policy.tsv
result_root=$project_root/results/evaluator/maize_v2/wgdi_outgroups
working_root=${result_root}.working

for required in "$input_root/SHA256SUMS" "$diamond_bin" "$wgdi_bin" \
                "$parallel_bin" "$policy"; do
    [[ -s $required ]] || { echo "missing maize WGDI prerequisite: $required" >&2; exit 1; }
done
(cd "$input_root" && sha256sum -c SHA256SUMS >/dev/null)
if [[ -e $result_root || -e $working_root ]]; then
    echo "refusing to overwrite maize outgroup WGDI" >&2; exit 1
fi
mkdir -p "$working_root"/{db,blast,config,collinearity,logs}
{
    printf 'field\tvalue\n'
    printf 'code_commit\t%s\n' "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'access\tevaluator_only\n'
    printf 'query\tZea_mays_B73_NAM5\n'
    printf 'references\tSorghum_bicolor_NCBIv3,Setaria_italica_v2\n'
    printf 'expected_query_copy_multiplicity\t2\n'
    printf 'diamond_threads_per_job\t64\nwgdi_processes_per_job\t64\n'
    printf 'policy_sha256\t%s\n' "$(sha256sum "$policy" | awk '{print $1}')"
} > "$working_root/run_contract.tsv"

make_db() {
    local short=$1
    /usr/bin/time -v -o "$working_root/logs/$short.makedb.time.txt" \
        "$diamond_bin" makedb --in "$input_root/$short/$short.wgdi.pep.fa" \
        --db "$working_root/db/$short" \
        > "$working_root/logs/$short.makedb.stdout.log" \
        2> "$working_root/logs/$short.makedb.stderr.log"
}
export -f make_db
export working_root input_root diamond_bin
printf '%s\n' sbi sit | "$parallel_bin" --jobs 2 make_db {}

run_diamond() {
    local pair=$1 ref=$2
    /usr/bin/time -v -o "$working_root/logs/$pair.diamond.time.txt" \
        "$diamond_bin" blastp --query "$input_root/zma/zma.wgdi.pep.fa" \
        --db "$working_root/db/$ref" --out "$working_root/blast/$pair.tsv" \
        --outfmt 6 --evalue 1e-5 --max-target-seqs 20 --more-sensitive \
        --threads 64 > "$working_root/logs/$pair.diamond.stdout.log" \
        2> "$working_root/logs/$pair.diamond.stderr.log"
}
export -f run_diamond
printf '%s\t%s\n' zma_vs_sbi sbi zma_vs_sit sit \
    | "$parallel_bin" --jobs 2 --delay 1 --colsep '\t' run_diamond {1} {2}

for entry in zma_vs_sbi:sbi zma_vs_sit:sit; do
    pair=${entry%%:*}; ref=${entry#*:}; config=$working_root/config/$pair.conf
    {
        printf '[collinearity]\n'
        printf 'gff1 = %s\ngff2 = %s\n' "$input_root/zma/zma.wgdi.gff" \
            "$input_root/$ref/$ref.wgdi.gff"
        printf 'lens1 = %s\nlens2 = %s\n' "$input_root/zma/zma.wgdi.lens" \
            "$input_root/$ref/$ref.wgdi.lens"
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
printf '%s\n' zma_vs_sbi zma_vs_sit \
    | "$parallel_bin" --jobs 2 --delay 1 run_wgdi {}

for pair in zma_vs_sbi zma_vs_sit; do
    [[ -s $working_root/collinearity/$pair.tsv ]] || {
        echo "missing maize WGDI result: $pair" >&2; exit 1;
    }
done
(
    cd "$working_root"
    find . -type f \( -name '*.tsv' -o -name '*.conf' -o -name '*.json' \) \
        -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS
    sha256sum -c SHA256SUMS >/dev/null
)
du -sb "$working_root" > "$working_root/disk_bytes.txt"
mv "$working_root" "$result_root"
printf 'maize outgroup WGDI frozen: %s\n' "$result_root"
