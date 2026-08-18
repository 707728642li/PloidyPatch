#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then echo "usage: $0 PROJECT_ROOT" >&2; exit 2; fi
project_root=$(realpath "$1")
code_root=$project_root/code
python_bin=$project_root/envs/ploidypatch-dev/bin/python
synteny_env=$project_root/envs/ploidypatch-synteny
syngap_env=$project_root/envs/ploidypatch-syngap
diamond_bin=$synteny_env/bin/diamond
wgdi_bin=$synteny_env/bin/wgdi
bundle=$project_root/data/derived/external_inputs/apple_v0.3/target_apple
genome=$bundle/primary_chromosomes.genome.fa
fai=$genome.fai
base_gff=$bundle/primary_chromosomes.gff3
method_root=$project_root/results/copy_collapse/external/apple_v0.3_method_trio
candidate_root=$method_root/consensus/primary_union/complete_control
candidate_gff=$candidate_root/candidate.gff3
composite_root=$project_root/results/models/ploidypatch_ranker_v0.4
result_root=$project_root/results/natural/apple_gddh13_v0.4/discovery/self_wgd
working_root=${result_root}.working
prefix=mdx_gddh13_natural_union_v04

verify_tree() { (cd "$1" && sha256sum -c SHA256SUMS >/dev/null); }
for required in "$python_bin" "$diamond_bin" "$wgdi_bin" \
    "$syngap_env/bin/gffread" "$genome" "$fai" "$base_gff" \
    "$candidate_gff" "$candidate_root/decisions.tsv" \
    "$candidate_root/candidate.gff3.manifest.json" "$method_root/SHA256SUMS" \
    "$composite_root/SHA256SUMS" \
    "$code_root/scripts/build_wgdi_source_alias_gff.py"; do
    [[ -s $required ]] || { echo "missing apple natural self-WGD input: $required" >&2; exit 1; }
done
verify_tree "$method_root"
verify_tree "$composite_root"
[[ ! -e $result_root && ! -e $working_root ]] || {
    echo "refusing to overwrite apple natural self-WGD" >&2; exit 1;
}
mkdir -p "$working_root"/{input,db,blast,config,collinearity,pairs,selected,logs,freeze}
{
    printf 'field\tvalue\ncode_commit\t%s\n' \
        "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'split\tnatural_current_annotation_discovery_v0.4\n'
    printf 'target\tMalus_domestica_GDDH13_v1.1\n'
    printf 'candidate_source\tchecksum_frozen_complete_annotation_method_union\n'
    printf 'validation_sequence_access\tfalse\nRNA_access\tfalse\n'
    printf 'candidate_and_rank_freeze_precedes_validation_access\ttrue\n'
    printf 'wgd_event\tmalus_lineage_wgd_natural_candidate_recomputed\n'
    printf 'min_block_pairs\t20\nrequire_different_seqids\ttrue\n'
    printf 'require_reciprocal_unique\ttrue\n'
    printf 'candidate_candidate_pair_policy\treject_as_circular\n'
    printf 'diamond_threads\t64\nwgdi_processes\t64\n'
    printf 'automatic_approval\tfalse\n'
    printf 'method_trio_sha256sums_sha256\t%s\n' \
        "$(sha256sum "$method_root/SHA256SUMS" | awk '{print $1}')"
    printf 'composite_model_sha256sums_sha256\t%s\n' \
        "$(sha256sum "$composite_root/SHA256SUMS" | awk '{print $1}')"
} > "$working_root/run_contract.tsv"
{
    printf 'role\tbytes\tsha256\tpath\n'
    for entry in "base_gff:$base_gff" "candidate_gff:$candidate_gff" \
        "pool_decisions:$candidate_root/decisions.tsv" \
        "pool_manifest:$candidate_root/candidate.gff3.manifest.json" \
        "genome:$genome" "fai:$fai" \
        "method_trio_freeze:$method_root/SHA256SUMS" \
        "composite_freeze:$composite_root/SHA256SUMS"; do
        role=${entry%%:*}; path=${entry#*:}
        printf '%s\t%s\t%s\t%s\n' "$role" "$(stat -Lc %s "$path")" \
            "$(sha256sum "$path" | awk '{print $1}')" "$path"
    done
} > "$working_root/input_manifest.tsv"
{
    printf 'path\tsha256\n'
    for relative in scripts/run_apple_natural_self_wgd_v0.4.sh \
        scripts/build_wgdi_source_alias_gff.py \
        src/ploidypatch/self_wgd_pairs.py \
        src/ploidypatch/wgd_candidate_select.py \
        src/ploidypatch/synteny_io.py src/ploidypatch/cli.py; do
        printf '%s\t%s\n' "$relative" \
            "$(sha256sum "$code_root/$relative" | awk '{print $1}')"
    done
} > "$working_root/freeze/code_manifest.tsv"

/usr/bin/time -v -o "$working_root/logs/gffread_pep.time.txt" \
    conda run -p "$syngap_env" --no-capture-output \
    gffread "$candidate_gff" -g "$genome" \
    -y "$working_root/input/$prefix.all.pep.fa" -S \
    > "$working_root/logs/gffread_pep.stdout.log" \
    2> "$working_root/logs/gffread_pep.stderr.log"
cd "$code_root"
"$python_bin" -m ploidypatch.cli evidence prepare-wgdi \
    --gff "$candidate_gff" --protein "$working_root/input/$prefix.all.pep.fa" \
    --fai "$fai" --output-dir "$working_root/input" \
    --prefix "$prefix" --min-genes-per-seqid 100 \
    > "$working_root/logs/prepare.stdout.json" \
    2> "$working_root/logs/prepare.stderr.log"
source_alias_gff=$working_root/input/$prefix.source_alias.gff3
"$python_bin" scripts/build_wgdi_source_alias_gff.py \
    --source-gff "$candidate_gff" \
    --representatives "$working_root/input/$prefix.representatives.tsv" \
    --output-gff "$source_alias_gff" \
    > "$working_root/logs/source_alias.stdout.json" \
    2> "$working_root/logs/source_alias.stderr.log"
"$diamond_bin" makedb --in "$working_root/input/$prefix.wgdi.pep.fa" \
    --db "$working_root/db/$prefix" \
    > "$working_root/logs/diamond_makedb.stdout.log" \
    2> "$working_root/logs/diamond_makedb.stderr.log"
/usr/bin/time -v -o "$working_root/logs/diamond_self.time.txt" \
    "$diamond_bin" blastp --query "$working_root/input/$prefix.wgdi.pep.fa" \
    --db "$working_root/db/$prefix" \
    --out "$working_root/blast/${prefix}_self.tsv" --outfmt 6 --evalue 1e-5 \
    --max-target-seqs 20 --more-sensitive --threads 64 \
    > "$working_root/logs/diamond_self.stdout.log" \
    2> "$working_root/logs/diamond_self.stderr.log"
config=$working_root/config/${prefix}_self.conf
{
    printf '[collinearity]\n'
    printf 'gff1 = %s\ngff2 = %s\n' "$working_root/input/$prefix.wgdi.gff" \
        "$working_root/input/$prefix.wgdi.gff"
    printf 'lens1 = %s\nlens2 = %s\n' "$working_root/input/$prefix.wgdi.lens" \
        "$working_root/input/$prefix.wgdi.lens"
    printf 'blast = %s\n' "$working_root/blast/${prefix}_self.tsv"
    printf 'blast_reverse = false\ncomparison = genomes\nmultiple = 2\nprocess = 64\n'
    printf 'evalue = 1e-5\nscore = 100\ngrading = 50,40,25\nmg = 40,40\n'
    printf 'pvalue = 0.2\nrepeat_number = 20\nposition = order\n'
    printf 'savefile = %s\n' "$working_root/collinearity/${prefix}_self.tsv"
} > "$config"
/usr/bin/time -v -o "$working_root/logs/wgdi_self.time.txt" \
    "$wgdi_bin" -icl "$config" > "$working_root/logs/wgdi_self.stdout.log" \
    2> "$working_root/logs/wgdi_self.stderr.log"
"$python_bin" -m ploidypatch.cli evidence infer-self-wgd-pairs \
    --query-wgdi-gff "$working_root/input/$prefix.wgdi.gff" \
    --collinearity "$working_root/collinearity/${prefix}_self.tsv" \
    --source-gff "$source_alias_gff" \
    --wgd-event malus_lineage_wgd_natural_candidate_recomputed \
    --min-block-pairs 20 \
    --output-pairs "$working_root/pairs/${prefix}.self_wgd_pairs.tsv" \
    --decisions-tsv "$working_root/pairs/decisions.tsv" \
    > "$working_root/pairs/stdout.json" 2> "$working_root/pairs/stderr.log"
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
    [[ -s $output ]] || { echo "missing apple natural self-WGD output: $output" >&2; exit 1; }
done
(
    cd "$working_root"
    find . -type f ! -name SHA256SUMS -print0 | sort -z | xargs -0 sha256sum \
        > SHA256SUMS
    sha256sum -c SHA256SUMS >/dev/null
)
du -sb "$working_root" > "$working_root/disk_bytes.txt"
chmod -R a-w "$working_root"
mv "$working_root" "$result_root"
printf 'apple natural self-WGD frozen: %s\n' "$result_root"
