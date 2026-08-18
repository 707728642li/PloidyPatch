#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then echo "usage: $0 PROJECT_ROOT" >&2; exit 2; fi
project_root=$(realpath "$1")
input_root=$project_root/data/derived/external_evaluator/populus_v0.4_wgdi_inputs
protocol_root=$project_root/results/protocol_freezes/populus_external_v0.4
execution_root=$project_root/results/protocol_freezes/populus_external_v0.4_execution
code_root=$execution_root/source
export PYTHONPATH="$code_root/src${PYTHONPATH:+:$PYTHONPATH}"
environment_bindings=$execution_root/environment_bindings.tsv
[[ -s $environment_bindings ]] || { echo "missing frozen environment bindings" >&2; exit 1; }
synteny_env=$(awk -F '\t' '$1 == "ploidypatch-synteny" {print $2}' "$environment_bindings")
[[ $synteny_env == /* ]] || { echo "invalid frozen synteny binding" >&2; exit 1; }
diamond_bin=$synteny_env/bin/diamond
wgdi_bin=$synteny_env/bin/wgdi
parallel_bin=/data/codexli/software/conda/miniforge3/bin/parallel
result_root=$project_root/results/evaluator/populus/v0.4/wgdi
working_root=${result_root}.working
self_relative=scripts/run_populus_evaluator_wgdi_v0.4.sh

verify_implementation() {
    local relative=$1 manifest=$execution_root/implementation_manifest.tsv
    local rows=()
    mapfile -t rows < <(awk -F '\t' -v path="$relative" '$1 == path {print $2 "\t" $3}' "$manifest")
    [[ ${#rows[@]} -eq 1 ]] || { echo "execution freeze has no unique row for $relative" >&2; return 1; }
    local expected_bytes expected_sha
    IFS=$'\t' read -r expected_bytes expected_sha <<< "${rows[0]}"
    [[ $expected_bytes =~ ^[0-9]+$ && $expected_sha =~ ^[0-9a-f]{64}$ ]] || {
        echo "malformed execution implementation row for $relative" >&2; return 1;
    }
    [[ $(stat -Lc %s "$code_root/$relative") == "$expected_bytes" \
        && $(sha256sum "$code_root/$relative" | awk '{print $1}') == "$expected_sha" ]] || {
        echo "implementation differs from execution freeze: $relative" >&2; return 1;
    }
}

for required in "$input_root/SHA256SUMS" "$protocol_root/SHA256SUMS" \
                "$execution_root/SHA256SUMS" "$execution_root/implementation_manifest.tsv" \
                "$diamond_bin" "$wgdi_bin" "$parallel_bin"; do
    [[ -s $required ]] || { echo "missing Populus WGDI prerequisite: $required" >&2; exit 1; }
done
(cd "$input_root" && sha256sum -c SHA256SUMS >/dev/null)
(cd "$protocol_root" && sha256sum -c SHA256SUMS >/dev/null)
(cd "$execution_root" && sha256sum -c SHA256SUMS >/dev/null)
verify_implementation "$self_relative"
[[ ! -e $result_root && ! -e $working_root ]] || {
    echo "refusing to overwrite Populus evaluator WGDI" >&2; exit 1;
}
mkdir -p "$working_root"/{db,blast,config,collinearity,logs}
record_invalid() {
    local status=$?
    if [[ -d $working_root ]]; then
        printf 'field\tvalue\nformal_status\tinvalid_run\nstage\tevaluator_wgdi\nexit_status\t%s\n' \
            "$status" > "$working_root/invalid_run.tsv" || true
    fi
    exit "$status"
}
trap record_invalid ERR
{
    printf 'field\tvalue\n'
    printf 'code_commit\t%s\n' "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'access\tevaluator_only_after_protocol_and_execution_freeze\n'
    printf 'query\tPopulus_trichocarpa_JGI_v4.0_annotation_v4.1\n'
    printf 'self_comparison\ttrue\nreferences\tManihot_esculenta_v6,Ricinus_communis_Wild_castor\n'
    printf 'candidate_reference_access\tfalse\ntruth_label_generation\tfalse\n'
    printf 'truth_pair_table_enumeration\tfalse\n'
    printf 'diamond_threads_per_job\t64\nwgdi_processes_per_job\t64\n'
    printf 'protocol_freeze_sha256sums_sha256\t%s\n' \
        "$(sha256sum "$protocol_root/SHA256SUMS" | awk '{print $1}')"
    printf 'execution_freeze_sha256sums_sha256\t%s\n' \
        "$(sha256sum "$execution_root/SHA256SUMS" | awk '{print $1}')"
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
printf '%s\n' ptr mes rco | "$parallel_bin" --jobs 3 --delay 1 make_db {}

run_diamond() {
    local pair=$1 ref=$2
    /usr/bin/time -v -o "$working_root/logs/$pair.diamond.time.txt" \
        "$diamond_bin" blastp --query "$input_root/ptr/ptr.wgdi.pep.fa" \
        --db "$working_root/db/$ref" --out "$working_root/blast/$pair.tsv" \
        --outfmt 6 --evalue 1e-5 --max-target-seqs 20 --more-sensitive \
        --threads 64 > "$working_root/logs/$pair.diamond.stdout.log" \
        2> "$working_root/logs/$pair.diamond.stderr.log"
}
export -f run_diamond
printf '%s\t%s\n' ptr_self ptr ptr_vs_mes mes ptr_vs_rco rco \
    | "$parallel_bin" --jobs 2 --delay 1 --colsep '\t' run_diamond {1} {2}

for entry in ptr_self:ptr ptr_vs_mes:mes ptr_vs_rco:rco; do
    pair=${entry%%:*}; ref=${entry#*:}; config=$working_root/config/$pair.conf
    {
        printf '[collinearity]\n'
        printf 'gff1 = %s\ngff2 = %s\n' "$input_root/ptr/ptr.wgdi.gff" \
            "$input_root/$ref/$ref.wgdi.gff"
        printf 'lens1 = %s\nlens2 = %s\n' "$input_root/ptr/ptr.wgdi.lens" \
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
printf '%s\n' ptr_self ptr_vs_mes ptr_vs_rco \
    | "$parallel_bin" --jobs 2 --delay 1 run_wgdi {}

printf 'comparison\tbytes\tlines\n' > "$working_root/collinearity_counts.tsv"
for pair in ptr_self ptr_vs_mes ptr_vs_rco; do
    output=$working_root/collinearity/$pair.tsv
    [[ -s $output ]] || { echo "missing Populus WGDI result: $pair" >&2; exit 1; }
    printf '%s\t%s\t%s\n' "$pair" "$(stat -Lc %s "$output")" "$(wc -l < "$output")" \
        >> "$working_root/collinearity_counts.tsv"
done
printf 'field\tvalue\nformal_status\tvalid_pre_pair_evidence\nstage\tevaluator_wgdi\n' \
    > "$working_root/stage_status.tsv"
du -sb "$working_root" > "$working_root/disk_bytes.txt"
(
    cd "$working_root"
    find . -type f \( -name '*.tsv' -o -name '*.conf' -o -name '*.json' \
        -o -name disk_bytes.txt \) \
        ! -name SHA256SUMS -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS
    sha256sum -c SHA256SUMS >/dev/null
)
trap - ERR
mv "$working_root" "$result_root"
printf 'Populus evaluator WGDI frozen: %s\n' "$result_root"
