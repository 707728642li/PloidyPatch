#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "usage: $0 PROJECT_ROOT" >&2
    exit 2
fi

project_root=$(realpath "$1")
code_root=$project_root/code
python_bin=$project_root/envs/ploidypatch-dev/bin/python
synteny_env=$project_root/envs/ploidypatch-synteny
diamond_bin=$synteny_env/bin/diamond
wgdi_bin=$synteny_env/bin/wgdi
benchmark_root=$project_root/benchmark/structure/copy_collapse_v0.1/bna_daae/annotation_copy_collapse_seed20260814
bundle_root=$project_root/data/derived/normalized_bundles/v0.1/bna_daae_primary
source_gff=$bundle_root/primary_chromosomes.gff3
target_protein=$bundle_root/primary_chromosomes.protein.fa
target_fai=$bundle_root/primary_chromosomes.genome.fa.fai
perturbed_gff=$benchmark_root/blind/perturbed.gff3
truth=$benchmark_root/evaluator/truth/hidden_truth.json
target_wgdi_root=$project_root/data/derived/wgdi_inputs/v0.1
database_root=$project_root/results/evidence/wgdi/brassica_v0.1/db
baseline_root=$project_root/results/copy_collapse/miniprot_brassica_v0.1
control_candidate=$project_root/results/baselines/miniprot_brassica_v0.1/complete_control_synteny_v0.2/selected/candidate.gff3
chromosome_pairs=$code_root/config/bna_daae_progenitor_chromosome_pairs_v0.1.tsv
result_root=$project_root/results/copy_collapse/wgdi_gap_brassica_v0.1
working_root=${result_root}.working

for required in "$python_bin" "$diamond_bin" "$wgdi_bin" "$source_gff" \
                "$target_protein" "$target_fai" "$perturbed_gff" "$truth" \
                "$database_root/bra_a.dmnd" "$database_root/bol_c.dmnd" \
                "$target_wgdi_root/bra_a/bra_a.wgdi.gff" \
                "$target_wgdi_root/bol_c/bol_c.wgdi.gff" \
                "$baseline_root/blind/candidate.gff3" \
                "$baseline_root/blind/decisions.tsv" "$control_candidate" \
                "$chromosome_pairs"; do
    if [[ ! -s $required ]]; then
        echo "missing or empty copy-collapse WGDI input: $required" >&2
        exit 1
    fi
done
if [[ -e $result_root || -e $working_root ]]; then
    echo "refusing to overwrite copy-collapse WGDI result: $result_root" >&2
    exit 1
fi

mkdir -p "$working_root/query" "$working_root/blast" \
    "$working_root/config" "$working_root/collinearity" \
    "$working_root/gaps" "$working_root/selected" "$working_root/logs"
{
    printf 'field\tvalue\n'
    printf 'code_commit\t%s\n' "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'candidate_truth_access\tfalse\n'
    printf 'evaluator_truth_access\ttrue\n'
    printf 'pair_tsv_candidate_access\tfalse\n'
    printf 'query_annotation\tblind_perturbed\n'
    printf 'query_protein_selection\tannotation_present_representatives_only\n'
    printf 'expected_chromosome_policy\tassembly_defined_A_or_C_to_declared_progenitor\n'
    printf 'max_query_intervening_genes\t0\n'
    printf 'min_target_excess_genes\t1\n'
    printf 'max_target_gap_genes\t5\n'
    printf 'max_query_locus_bp\t500000\n'
    printf 'synteny_gap_module_sha256\t%s\n' \
        "$(sha256sum "$code_root/src/ploidypatch/synteny_gap.py" | awk '{print $1}')"
} > "$working_root/run_contract.tsv"

cd "$code_root"
/usr/bin/time -v -o "$working_root/query/resource.time.txt" \
    "$python_bin" -m ploidypatch.cli evidence prepare-wgdi \
        --gff "$perturbed_gff" \
        --protein "$target_protein" \
        --fai "$target_fai" \
        --output-dir "$working_root/query" \
        --prefix bna_daae_blind \
        --min-genes-per-seqid 100 \
        --primary-chromosomes-only \
        > "$working_root/query/stdout.json" \
        2> "$working_root/query/stderr.log"

run_similarity() {
    local source=$1
    /usr/bin/time -v -o "$working_root/logs/${source}.diamond.time.txt" \
        "$diamond_bin" blastp \
            --query "$working_root/query/bna_daae_blind.wgdi.pep.fa" \
            --db "$database_root/$source" \
            --out "$working_root/blast/bna_daae_blind_vs_${source}.tsv" \
            --outfmt 6 --evalue 1e-5 --max-target-seqs 20 \
            --more-sensitive --threads 48 \
            > "$working_root/logs/${source}.diamond.stdout.log" \
            2> "$working_root/logs/${source}.diamond.stderr.log"
}
run_similarity bra_a & pid_bra=$!
run_similarity bol_c & pid_bol=$!
status=0
wait "$pid_bra" || status=1
wait "$pid_bol" || status=1
if [[ $status -ne 0 ]]; then
    echo "copy-collapse WGDI similarity search failed" >&2
    exit 1
fi

write_and_run_wgdi() {
    local source=$1
    local config=$working_root/config/bna_daae_blind_vs_${source}.conf
    {
        printf '[collinearity]\n'
        printf 'gff1 = %s\n' "$working_root/query/bna_daae_blind.wgdi.gff"
        printf 'gff2 = %s\n' "$target_wgdi_root/$source/$source.wgdi.gff"
        printf 'lens1 = %s\n' "$working_root/query/bna_daae_blind.wgdi.lens"
        printf 'lens2 = %s\n' "$target_wgdi_root/$source/$source.wgdi.lens"
        printf 'blast = %s\n' "$working_root/blast/bna_daae_blind_vs_${source}.tsv"
        printf 'blast_reverse = false\ncomparison = genomes\nmultiple = 2\n'
        printf 'process = 32\nevalue = 1e-5\nscore = 100\n'
        printf 'grading = 50,40,25\nmg = 40,40\npvalue = 0.2\n'
        printf 'repeat_number = 20\nposition = order\n'
        printf 'savefile = %s\n' "$working_root/collinearity/bna_daae_blind_vs_${source}.tsv"
    } > "$config"
    /usr/bin/time -v -o "$working_root/logs/${source}.wgdi.time.txt" \
        "$wgdi_bin" -icl "$config" \
        > "$working_root/logs/${source}.wgdi.stdout.log" \
        2> "$working_root/logs/${source}.wgdi.stderr.log"
}
write_and_run_wgdi bra_a & pid_bra=$!
write_and_run_wgdi bol_c & pid_bol=$!
status=0
wait "$pid_bra" || status=1
wait "$pid_bol" || status=1
if [[ $status -ne 0 ]]; then
    echo "copy-collapse WGDI collinearity failed" >&2
    exit 1
fi

infer_gap() {
    local source=$1
    /usr/bin/time -v -o "$working_root/gaps/${source}.resource.time.txt" \
        "$python_bin" -m ploidypatch.cli evidence infer-synteny-gaps \
            --query-wgdi-gff "$working_root/query/bna_daae_blind.wgdi.gff" \
            --target-wgdi-gff "$target_wgdi_root/$source/$source.wgdi.gff" \
            --collinearity "$working_root/collinearity/bna_daae_blind_vs_${source}.tsv" \
            --source-label "$source" \
            --expected-chromosome-pairs "$chromosome_pairs" \
            --max-query-intervening-genes 0 \
            --min-target-excess-genes 1 \
            --max-target-gap-genes 5 \
            --max-query-locus-bp 500000 \
            --output-tsv "$working_root/gaps/${source}.tsv" \
            > "$working_root/gaps/${source}.stdout.json" \
            2> "$working_root/gaps/${source}.stderr.log"
}
infer_gap bra_a & pid_bra=$!
infer_gap bol_c & pid_bol=$!
status=0
wait "$pid_bra" || status=1
wait "$pid_bol" || status=1
if [[ $status -ne 0 ]]; then
    echo "copy-collapse WGDI gap inference failed" >&2
    exit 1
fi

/usr/bin/time -v -o "$working_root/selected/resource.time.txt" \
    "$python_bin" -m ploidypatch.cli evidence select-synteny-gap-models \
        --gaps "$working_root/gaps/bra_a.tsv" \
        --gaps "$working_root/gaps/bol_c.tsv" \
        --baseline-decisions "$baseline_root/blind/decisions.tsv" \
        --adapted-candidate-gff "$baseline_root/blind/candidate.gff3" \
        --output-selection "$working_root/selected/selection.tsv" \
        --output-candidate-gff "$working_root/selected/candidate.gff3" \
        > "$working_root/selected/stdout.json" \
        2> "$working_root/selected/stderr.log"

{
    printf 'role\tbytes\tsha256\tpath\n'
    path=$working_root/selected/candidate.gff3
    printf 'blind_candidate\t%s\t%s\t%s\n' "$(stat -Lc %s "$path")" \
        "$(sha256sum "$path" | awk '{print $1}')" "$path"
    printf 'complete_control_candidate\t%s\t%s\t%s\n' \
        "$(stat -Lc %s "$control_candidate")" \
        "$(sha256sum "$control_candidate" | awk '{print $1}')" \
        "$control_candidate"
} > "$working_root/candidate_freeze.tsv"

/usr/bin/time -v -o "$working_root/score.resource.time.txt" \
    "$python_bin" -m ploidypatch.cli benchmark score \
        --source-gff "$source_gff" \
        --perturbed-gff "$perturbed_gff" \
        --candidate-gff "$working_root/selected/candidate.gff3" \
        --control-candidate-gff "$control_candidate" \
        --truth "$truth" \
        --include-event-details \
        > "$working_root/score.json" \
        2> "$working_root/score.stderr.log"
if ! grep -q '"grade": "pass"' "$working_root/score.json"; then
    echo "copy-collapse WGDI-gap score quality gate failed" >&2
    exit 1
fi

for output in "$working_root/query/bna_daae_blind.wgdi.gff" \
              "$working_root/collinearity/bna_daae_blind_vs_bra_a.tsv" \
              "$working_root/collinearity/bna_daae_blind_vs_bol_c.tsv" \
              "$working_root/gaps/bra_a.tsv" "$working_root/gaps/bol_c.tsv" \
              "$working_root/selected/selection.tsv" \
              "$working_root/selected/candidate.gff3" \
              "$working_root/candidate_freeze.tsv" "$working_root/score.json"; do
    if [[ ! -s $output ]]; then
        echo "missing or empty copy-collapse WGDI output: $output" >&2
        exit 1
    fi
done
(
    cd "$working_root"
    find . -type f \( -name '*.gff3' -o -name '*.json' -o -name '*.tsv' \) \
        -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS
)
du -sb "$working_root" > "$working_root/disk_bytes.txt"
mv "$working_root" "$result_root"
printf 'copy-collapse WGDI-gap tier frozen: %s\n' "$result_root"
