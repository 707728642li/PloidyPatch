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
input_root=${PLOIDYPATCH_STAGED_INPUT_ROOT:?PLOIDYPATCH_STAGED_INPUT_ROOT is required}
miniprot_root=$project_root/results/baselines/actinidia_v0.5/miniprot
bundle=$miniprot_root/normalized/target
genome=$bundle/primary_chromosomes.genome.fa
fai=$bundle/primary_chromosomes.genome.fa.fai
blind_benchmark_root=${PLOIDYPATCH_BLIND_BENCHMARK_ROOT:?PLOIDYPATCH_BLIND_BENCHMARK_ROOT is required}
base_gff=$blind_benchmark_root/perturbed.gff3
method_root=$project_root/results/copy_collapse/external/actinidia_v0.5_method_trio
candidate_gff=$method_root/consensus/primary_union/blind/candidate.gff3
candidate_manifest=$method_root/consensus/primary_union/blind/candidate.gff3.manifest.json
protocol_root=${PLOIDYPATCH_PROTOCOL_FREEZE:?PLOIDYPATCH_PROTOCOL_FREEZE is required}
execution_root=${PLOIDYPATCH_EXECUTION_FREEZE:?PLOIDYPATCH_EXECUTION_FREEZE is required}
contract_path=${PLOIDYPATCH_HOLDOUT_CONTRACT:?PLOIDYPATCH_HOLDOUT_CONTRACT is required}
context_verifier=$code_root/scripts/verify_external_holdout_blind_context_v0.5.py
result_root=$project_root/results/copy_collapse/external/actinidia_v0.5_blind_self_wgd
working_root=${result_root}.working
prefix=red5_union_candidate_blind

[[ ${PLOIDYPATCH_BLIND_RUNNER:-} == 1 ]] || {
    echo "Actinidia self-WGD must run inside the frozen blind runner" >&2; exit 1;
}
[[ ${PLOIDYPATCH_NETWORK_ACCESS:-} == none ]] || {
    echo "Actinidia self-WGD requires a network-disabled namespace" >&2; exit 1;
}
for forbidden in /nas_data "$input_root/evaluator_only" "$blind_benchmark_root/evaluator" \
    "$blind_benchmark_root/truth" "$blind_benchmark_root/complete"; do
    [[ ! -e $forbidden ]] || { echo "forbidden blind-runner path is visible: $forbidden" >&2; exit 1; }
done
if [[ -r /proc/self/mountinfo ]] && grep -Eq '/nas_data|/evaluator_only|/target_complete|/truth_references' /proc/self/mountinfo; then
    echo "forbidden evaluator or NAS mount detected in blind namespace" >&2; exit 1
fi

verify_tree() { (cd "$1" && sha256sum -c SHA256SUMS >/dev/null); }
verify_implementation() {
    local relative=$1 expected observed
    expected=$(awk -F '\t' -v path="$relative" 'NR > 1 && $1 == path {print $3}' \
        "$execution_root/implementation_manifest.tsv")
    observed=$(sha256sum "$code_root/$relative" | awk '{print $1}')
    [[ $expected =~ ^[0-9a-f]{64}$ && $observed == "$expected" ]] || {
        echo "implementation differs from execution freeze: $relative" >&2; exit 1;
    }
}

for required in "$python_bin" "$diamond_bin" "$wgdi_bin" \
    "$blind_benchmark_root/blind_manifest.json" "$blind_benchmark_root/SHA256SUMS" \
    "$syngap_env/bin/gffread" "$genome" "$fai" "$base_gff" \
    "$candidate_gff" "$candidate_manifest" "$method_root/SHA256SUMS" \
    "$miniprot_root/SHA256SUMS" "$protocol_root/SHA256SUMS" \
    "$execution_root/SHA256SUMS" "$contract_path" "$input_root/role_manifest.tsv" \
    "$input_root/role_contract.json" "$context_verifier" \
    "$code_root/scripts/build_wgdi_source_alias_gff.py"; do
    [[ -s $required ]] || { echo "missing Actinidia blind self-WGD input: $required" >&2; exit 1; }
done
verify_tree "$protocol_root"
verify_tree "$execution_root"
verify_tree "$miniprot_root"
verify_tree "$method_root"
verify_implementation scripts/run_actinidia_blind_union_self_wgd_v0.5.sh
verify_implementation scripts/verify_external_holdout_blind_context_v0.5.py
verify_implementation scripts/build_wgdi_source_alias_gff.py
verify_implementation src/ploidypatch/self_wgd_pairs.py
verify_implementation src/ploidypatch/wgd_candidate_select.py
[[ ! -e $result_root && ! -e $working_root ]] || {
    echo "refusing to overwrite Actinidia blind self-WGD" >&2; exit 1;
}
mkdir -p "$working_root"/{input,db,blast,config,collinearity,pairs,selected,logs,freeze}
cp /proc/self/mountinfo "$working_root/freeze/blind_runner.mountinfo"
PYTHONPATH="$code_root/src" "$python_bin" "$context_verifier" \
    --input-root "$input_root" --contract "$contract_path" \
    --protocol-freeze "$protocol_root" --execution-freeze "$execution_root" \
    --blind-benchmark-root "$blind_benchmark_root" \
    --expected-holdout-id actinidia_red5_v0.5 --expected-primary-chromosomes 29 \
    --output-json "$working_root/freeze/blind_context.json"
{
    printf 'field\tvalue\ncode_commit\t%s\n' "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'mode\tblind\nsplit\ttarget_level_predeclared_untouched_secondary_replication\n'
    printf 'candidate_stage_truth_access\tfalse\npreperturbation_pair_access\tfalse\n'
    printf 'complete_target_annotation_access\tfalse\nevaluator_reference_access\tfalse\n'
    printf 'candidate_source\tchain_preserving_method_union\n'
    printf 'wgd_event\tactinidia_specific_ad_alpha_blind_recomputed\nmin_block_pairs\t20\n'
    printf 'require_different_seqids\ttrue\nrequire_reciprocal_unique\ttrue\n'
    printf 'candidate_candidate_pair_policy\treject_as_circular\n'
    printf 'diamond_threads\t64\nwgdi_processes\t64\nautomatic_approval\tfalse\n'
    printf 'protocol_freeze_sha256sums_sha256\t%s\n' \
        "$(sha256sum "$protocol_root/SHA256SUMS" | awk '{print $1}')"
    printf 'execution_freeze_sha256sums_sha256\t%s\n' \
        "$(sha256sum "$execution_root/SHA256SUMS" | awk '{print $1}')"
    printf 'identifier_alias_adapter_sha256\t%s\n' \
        "$(sha256sum "$code_root/scripts/build_wgdi_source_alias_gff.py" | awk '{print $1}')"
} > "$working_root/run_contract.tsv"
{
    printf 'role\tbytes\tsha256\tpath\n'
    for entry in "base_gff:$base_gff" "candidate_gff:$candidate_gff" \
        "candidate_pool_manifest:$candidate_manifest" "genome:$genome" "fai:$fai"; do
        role=${entry%%:*}; path=${entry#*:}
        printf '%s\t%s\t%s\t%s\n' "$role" "$(stat -Lc %s "$path")" \
            "$(sha256sum "$path" | awk '{print $1}')" "$path"
    done
} > "$working_root/input_manifest.tsv"

/usr/bin/time -v -o "$working_root/logs/gffread_pep.time.txt" \
    conda run -p "$syngap_env" --no-capture-output \
    gffread "$candidate_gff" -g "$genome" -y "$working_root/input/$prefix.all.pep.fa" -S \
    > "$working_root/logs/gffread_pep.stdout.log" \
    2> "$working_root/logs/gffread_pep.stderr.log"
cd "$code_root"
/usr/bin/time -v -o "$working_root/logs/prepare.time.txt" \
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
/usr/bin/time -v -o "$working_root/logs/diamond_makedb.time.txt" \
    "$diamond_bin" makedb --in "$working_root/input/$prefix.wgdi.pep.fa" \
    --db "$working_root/db/$prefix" \
    > "$working_root/logs/diamond_makedb.stdout.log" \
    2> "$working_root/logs/diamond_makedb.stderr.log"
/usr/bin/time -v -o "$working_root/logs/diamond_self.time.txt" \
    "$diamond_bin" blastp --query "$working_root/input/$prefix.wgdi.pep.fa" \
    --db "$working_root/db/$prefix" --out "$working_root/blast/${prefix}_self.tsv" \
    --outfmt 6 --evalue 1e-5 --max-target-seqs 20 --more-sensitive --threads 64 \
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
/usr/bin/time -v -o "$working_root/logs/pair_inference.time.txt" \
    "$python_bin" -m ploidypatch.cli evidence infer-self-wgd-pairs \
    --query-wgdi-gff "$working_root/input/$prefix.wgdi.gff" \
    --collinearity "$working_root/collinearity/${prefix}_self.tsv" \
    --source-gff "$source_alias_gff" \
    --wgd-event actinidia_specific_ad_alpha_blind_recomputed --min-block-pairs 20 \
    --output-pairs "$working_root/pairs/${prefix}.self_wgd_pairs.tsv" \
    --decisions-tsv "$working_root/pairs/decisions.tsv" \
    > "$working_root/pairs/stdout.json" 2> "$working_root/pairs/stderr.log"
/usr/bin/time -v -o "$working_root/logs/candidate_selection.time.txt" \
    "$python_bin" -m ploidypatch.cli evidence select-wgd-supported-candidates \
    --base-gff "$base_gff" --candidate-gff "$candidate_gff" \
    --pairs "$working_root/pairs/${prefix}.self_wgd_pairs.tsv" \
    --output-gff "$working_root/selected/candidate.gff3" \
    --selection-tsv "$working_root/selected/selection.tsv" \
    > "$working_root/selected/stdout.json" 2> "$working_root/selected/stderr.log"
for output in "$working_root/input/$prefix.wgdi.gff" \
    "$working_root/input/$prefix.wgdi.pep.fa" "$working_root/blast/${prefix}_self.tsv" \
    "$working_root/collinearity/${prefix}_self.tsv" \
    "$working_root/pairs/${prefix}.self_wgd_pairs.tsv" \
    "$working_root/selected/candidate.gff3" "$working_root/selected/selection.tsv"; do
    [[ -s $output ]] || { echo "missing Actinidia blind self-WGD output: $output" >&2; exit 1; }
done
(
    cd "$working_root"
    find . -type f \( -name '*.json' -o -name '*.tsv' -o -name '*.gff' \
        -o -name '*.gff3' -o -name '*.fa' -o -name '*.conf' \
        -o -name '*.mountinfo' \) -print0 \
        | sort -z | xargs -0 sha256sum > SHA256SUMS
    sha256sum -c SHA256SUMS >/dev/null
)
du -sb "$working_root" > "$working_root/disk_bytes.txt"
mv "$working_root" "$result_root"
printf 'Actinidia blind candidate self-WGD frozen: %s\n' "$result_root"
