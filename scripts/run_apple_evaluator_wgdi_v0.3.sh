#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then echo "usage: $0 PROJECT_ROOT" >&2; exit 2; fi
project_root=$(realpath "$1")
code_root=$project_root/code
input_root=$project_root/data/derived/external_evaluator/apple_v0.3_wgdi_inputs
protocol_root=$project_root/results/protocol_freezes/apple_external_v0.3
execution_root=$project_root/results/protocol_freezes/apple_external_v0.3_execution
synteny_env=$project_root/envs/ploidypatch-synteny
diamond_bin=$synteny_env/bin/diamond
wgdi_bin=$synteny_env/bin/wgdi
parallel_bin=/data/codexli/software/conda/miniforge3/bin/parallel
result_root=$project_root/results/evaluator/apple_v0.3/wgdi
working_root=${result_root}.working

for required in "$input_root/SHA256SUMS" "$protocol_root/SHA256SUMS" \
                "$execution_root/SHA256SUMS" \
                "$diamond_bin" "$wgdi_bin" "$parallel_bin"; do
    [[ -s $required ]] || { echo "missing apple WGDI prerequisite: $required" >&2; exit 1; }
done
(cd "$input_root" && sha256sum -c SHA256SUMS >/dev/null)
(cd "$protocol_root" && sha256sum -c SHA256SUMS >/dev/null)
(cd "$execution_root" && sha256sum -c SHA256SUMS >/dev/null)
expected_self=$(awk -F '\t' '$1 == "scripts/run_apple_evaluator_wgdi_v0.3.sh" {print $3}' \
    "$execution_root/implementation_manifest.tsv")
[[ -n $expected_self && $(sha256sum "$code_root/scripts/run_apple_evaluator_wgdi_v0.3.sh" | awk '{print $1}') == "$expected_self" ]] || {
    echo "apple evaluator WGDI script differs from execution freeze" >&2; exit 1;
}
[[ ! -e $result_root && ! -e $working_root ]] || {
    echo "refusing to overwrite apple evaluator WGDI" >&2; exit 1;
}
mkdir -p "$working_root"/{db,blast,config,collinearity,logs}
{
    printf 'field\tvalue\n'
    printf 'code_commit\t%s\n' "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'access\tevaluator_only\nquery\tMalus_x_domestica_GDDH13_v1.1\n'
    printf 'self_comparison\ttrue\nreferences\tRosa_chinensis,Fragaria_vesca\n'
    printf 'candidate_reference_access\tfalse\nexternal_label_access\tfalse\n'
    printf 'diamond_threads_per_job\t64\nwgdi_processes_per_job\t64\n'
    printf 'protocol_freeze_sha256sums_sha256\t%s\n' \
        "$(sha256sum "$protocol_root/SHA256SUMS" | awk '{print $1}')"
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
printf '%s\n' mdx rch fve | "$parallel_bin" --jobs 3 --delay 1 make_db {}

run_diamond() {
    local pair=$1 ref=$2
    /usr/bin/time -v -o "$working_root/logs/$pair.diamond.time.txt" \
        "$diamond_bin" blastp --query "$input_root/mdx/mdx.wgdi.pep.fa" \
        --db "$working_root/db/$ref" --out "$working_root/blast/$pair.tsv" \
        --outfmt 6 --evalue 1e-5 --max-target-seqs 20 --more-sensitive \
        --threads 64 > "$working_root/logs/$pair.diamond.stdout.log" \
        2> "$working_root/logs/$pair.diamond.stderr.log"
}
export -f run_diamond
printf '%s\t%s\n' mdx_self mdx mdx_vs_rch rch mdx_vs_fve fve \
    | "$parallel_bin" --jobs 2 --delay 1 --colsep '\t' run_diamond {1} {2}

for entry in mdx_self:mdx mdx_vs_rch:rch mdx_vs_fve:fve; do
    pair=${entry%%:*}; ref=${entry#*:}; config=$working_root/config/$pair.conf
    {
        printf '[collinearity]\n'
        printf 'gff1 = %s\ngff2 = %s\n' "$input_root/mdx/mdx.wgdi.gff" \
            "$input_root/$ref/$ref.wgdi.gff"
        printf 'lens1 = %s\nlens2 = %s\n' "$input_root/mdx/mdx.wgdi.lens" \
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
printf '%s\n' mdx_self mdx_vs_rch mdx_vs_fve \
    | "$parallel_bin" --jobs 2 --delay 1 run_wgdi {}

for pair in mdx_self mdx_vs_rch mdx_vs_fve; do
    [[ -s $working_root/collinearity/$pair.tsv ]] || {
        echo "missing apple WGDI result: $pair" >&2; exit 1;
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
printf 'apple evaluator WGDI frozen: %s\n' "$result_root"
