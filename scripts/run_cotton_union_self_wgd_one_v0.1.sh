#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "usage: $0 PROJECT_ROOT blind|complete_control" >&2
    exit 2
fi

project_root=$(realpath "$1")
mode=$2
code_root=$project_root/code
python_bin=$project_root/envs/ploidypatch-dev/bin/python
synteny_env=$project_root/envs/ploidypatch-synteny
syngap_env=$project_root/envs/ploidypatch-syngap
diamond_bin=$synteny_env/bin/diamond
wgdi_bin=$synteny_env/bin/wgdi
bundle_root=$project_root/data/derived/holdout_inputs/cotton_v0.1/hirsutum
genome=$bundle_root/primary_chromosomes.genome.fa
fai=$bundle_root/primary_chromosomes.genome.fa.fai
benchmark_root=$project_root/benchmark/structure/copy_collapse_v0.1/ghi_ad/annotation_copy_collapse_seed20260817
method_root=$project_root/results/copy_collapse/holdout/cotton_method_trio_v0.1
case $mode in
    blind) base_gff=$benchmark_root/blind/perturbed.gff3 ;;
    complete_control) base_gff=$bundle_root/primary_chromosomes.gff3 ;;
    *) echo "mode must be blind or complete_control" >&2; exit 2 ;;
esac
candidate_gff=$method_root/consensus/union/$mode/candidate.gff3
result_root=$project_root/results/copy_collapse/holdout/cotton_union_self_wgd_v0.1/$mode
working_root=${result_root}.working
prefix=ghi_union_candidate_$mode

if [[ ! -s $method_root/SHA256SUMS ]]; then
    echo "cotton method-trio result is not frozen" >&2
    exit 1
fi
(cd "$method_root" && sha256sum -c SHA256SUMS >/dev/null)
for required in "$python_bin" "$diamond_bin" "$wgdi_bin" \
                "$syngap_env/bin/gffread" "$genome" "$fai" \
                "$base_gff" "$candidate_gff"; do
    if [[ ! -s $required ]]; then
        echo "missing or empty cotton candidate self-WGD input: $required" >&2
        exit 1
    fi
done
if [[ -e $result_root || -e $working_root ]]; then
    echo "refusing to overwrite cotton candidate self-WGD result" >&2
    exit 1
fi
mkdir -p "$working_root/input" "$working_root/db" "$working_root/blast" \
    "$working_root/config" "$working_root/collinearity" \
    "$working_root/pairs" "$working_root/selected" "$working_root/logs"
{
    printf 'field\tvalue\n'
    printf 'code_commit\t%s\n' "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'mode\t%s\n' "$mode"
    printf 'split\texternal_zero_retuning_holdout\n'
    printf 'formal_holdout_claim_allowed\ttrue\n'
    printf 'candidate_stage_truth_access\tfalse\n'
    printf 'preperturbation_pair_access\tfalse\n'
    printf 'candidate_source\tmethod_consensus_support1_union\n'
    printf 'protein_source\tgffread_from_mode_specific_candidate_GFF_and_target_genome\n'
    printf 'wgd_event\tgossypium_hirsutum_allotetraploidization\n'
    printf 'min_block_pairs\t20\n'
    printf 'require_different_seqids\ttrue\n'
    printf 'require_reciprocal_unique\ttrue\n'
    printf 'candidate_candidate_pair_policy\treject_as_circular\n'
    printf 'diamond_threads\t64\nwgdi_processes\t64\n'
} > "$working_root/run_contract.tsv"
{
    printf 'role\tbytes\tsha256\tpath\n'
    for entry in "base_gff:$base_gff" "candidate_gff:$candidate_gff" \
                 "genome:$genome" "fai:$fai"; do
        role=${entry%%:*}; path=${entry#*:}
        printf '%s\t%s\t%s\t%s\n' "$role" "$(stat -Lc %s "$path")" \
            "$(sha256sum "$path" | awk '{print $1}')" "$path"
    done
} > "$working_root/candidate_input_manifest.tsv"

/usr/bin/time -v -o "$working_root/logs/gffread_pep.time.txt" \
    conda run -p "$syngap_env" --no-capture-output \
    gffread "$candidate_gff" -g "$genome" \
        -y "$working_root/input/$prefix.all.pep.fa" -S \
        > "$working_root/logs/gffread_pep.stdout.log" \
        2> "$working_root/logs/gffread_pep.stderr.log"
cd "$code_root"
/usr/bin/time -v -o "$working_root/logs/prepare.time.txt" \
    "$python_bin" -m ploidypatch.cli evidence prepare-wgdi \
        --gff "$candidate_gff" \
        --protein "$working_root/input/$prefix.all.pep.fa" \
        --fai "$fai" --output-dir "$working_root/input" \
        --prefix "$prefix" --min-genes-per-seqid 100 \
        > "$working_root/logs/prepare.stdout.json" \
        2> "$working_root/logs/prepare.stderr.log"
/usr/bin/time -v -o "$working_root/logs/diamond_makedb.time.txt" \
    "$diamond_bin" makedb --in "$working_root/input/$prefix.wgdi.pep.fa" \
        --db "$working_root/db/$prefix" \
        > "$working_root/logs/diamond_makedb.stdout.log" \
        2> "$working_root/logs/diamond_makedb.stderr.log"
/usr/bin/time -v -o "$working_root/logs/diamond_self.time.txt" \
    "$diamond_bin" blastp --query "$working_root/input/$prefix.wgdi.pep.fa" \
        --db "$working_root/db/$prefix" \
        --out "$working_root/blast/${prefix}_self.tsv" \
        --outfmt 6 --evalue 1e-5 --max-target-seqs 20 \
        --more-sensitive --threads 64 \
        > "$working_root/logs/diamond_self.stdout.log" \
        2> "$working_root/logs/diamond_self.stderr.log"
config=$working_root/config/${prefix}_self.conf
{
    printf '[collinearity]\n'
    printf 'gff1 = %s\n' "$working_root/input/$prefix.wgdi.gff"
    printf 'gff2 = %s\n' "$working_root/input/$prefix.wgdi.gff"
    printf 'lens1 = %s\n' "$working_root/input/$prefix.wgdi.lens"
    printf 'lens2 = %s\n' "$working_root/input/$prefix.wgdi.lens"
    printf 'blast = %s\n' "$working_root/blast/${prefix}_self.tsv"
    printf 'blast_reverse = false\ncomparison = genomes\nmultiple = 2\n'
    printf 'process = 64\nevalue = 1e-5\nscore = 100\n'
    printf 'grading = 50,40,25\nmg = 40,40\npvalue = 0.2\n'
    printf 'repeat_number = 20\nposition = order\n'
    printf 'savefile = %s\n' "$working_root/collinearity/${prefix}_self.tsv"
} > "$config"
/usr/bin/time -v -o "$working_root/logs/wgdi_self.time.txt" \
    "$wgdi_bin" -icl "$config" \
        > "$working_root/logs/wgdi_self.stdout.log" \
        2> "$working_root/logs/wgdi_self.stderr.log"
/usr/bin/time -v -o "$working_root/logs/pair_inference.time.txt" \
    "$python_bin" -m ploidypatch.cli evidence infer-self-wgd-pairs \
        --query-wgdi-gff "$working_root/input/$prefix.wgdi.gff" \
        --collinearity "$working_root/collinearity/${prefix}_self.tsv" \
        --source-gff "$candidate_gff" \
        --wgd-event gossypium_hirsutum_allotetraploidization \
        --min-block-pairs 20 \
        --output-pairs "$working_root/pairs/${prefix}.self_wgd_pairs.tsv" \
        --decisions-tsv "$working_root/pairs/decisions.tsv" \
        > "$working_root/pairs/stdout.json" \
        2> "$working_root/pairs/stderr.log"
/usr/bin/time -v -o "$working_root/logs/candidate_selection.time.txt" \
    "$python_bin" -m ploidypatch.cli evidence select-wgd-supported-candidates \
        --base-gff "$base_gff" --candidate-gff "$candidate_gff" \
        --pairs "$working_root/pairs/${prefix}.self_wgd_pairs.tsv" \
        --output-gff "$working_root/selected/candidate.gff3" \
        --selection-tsv "$working_root/selected/selection.tsv" \
        > "$working_root/selected/stdout.json" \
        2> "$working_root/selected/stderr.log"
for output in "$working_root/input/$prefix.wgdi.gff" \
              "$working_root/input/$prefix.wgdi.pep.fa" \
              "$working_root/blast/${prefix}_self.tsv" \
              "$working_root/collinearity/${prefix}_self.tsv" \
              "$working_root/pairs/${prefix}.self_wgd_pairs.tsv" \
              "$working_root/selected/candidate.gff3" \
              "$working_root/selected/selection.tsv"; do
    if [[ ! -s $output ]]; then
        echo "missing cotton self-WGD output: $output" >&2
        exit 1
    fi
done
(
    cd "$working_root"
    find . -type f \( -name '*.json' -o -name '*.tsv' -o -name '*.gff' \
        -o -name '*.gff3' -o -name '*.fa' \) \
        -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS
    sha256sum -c SHA256SUMS >/dev/null
)
du -sb "$working_root" > "$working_root/disk_bytes.txt"
mv "$working_root" "$result_root"
printf 'cotton candidate self-WGD frozen: %s\n' "$result_root"
