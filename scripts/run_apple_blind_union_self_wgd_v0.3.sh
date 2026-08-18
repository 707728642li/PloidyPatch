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
fai=$bundle/primary_chromosomes.genome.fa.fai
benchmark=$project_root/benchmark/structure/copy_collapse_v0.3/mdx_gddh13/annotation_copy_collapse_seed20260831
base_gff=$benchmark/blind/perturbed.gff3
method_root=$project_root/results/copy_collapse/external/apple_v0.3_method_trio
candidate_gff=$method_root/consensus/primary_union/blind/candidate.gff3
protocol_root=$project_root/results/protocol_freezes/apple_external_v0.3
execution_root=$project_root/results/protocol_freezes/apple_external_v0.3_execution
result_root=$project_root/results/copy_collapse/external/apple_v0.3_blind_self_wgd
working_root=${result_root}.working
prefix=mdx_union_candidate_blind

for required in "$python_bin" "$diamond_bin" "$wgdi_bin" \
    "$syngap_env/bin/gffread" "$genome" "$fai" "$base_gff" \
    "$candidate_gff" "$method_root/SHA256SUMS" "$protocol_root/SHA256SUMS" \
    "$execution_root/SHA256SUMS" \
    "$code_root/scripts/build_wgdi_source_alias_gff.py"; do
    [[ -s $required ]] || { echo "missing apple blind self-WGD input: $required" >&2; exit 1; }
done
(cd "$method_root" && sha256sum -c SHA256SUMS >/dev/null)
(cd "$protocol_root" && sha256sum -c SHA256SUMS >/dev/null)
(cd "$execution_root" && sha256sum -c SHA256SUMS >/dev/null)
expected_self=$(awk -F '\t' '$1 == "scripts/run_apple_blind_union_self_wgd_v0.3.sh" {print $3}' \
    "$execution_root/implementation_manifest.tsv")
[[ -n $expected_self && $(sha256sum "$code_root/scripts/run_apple_blind_union_self_wgd_v0.3.sh" | awk '{print $1}') == "$expected_self" ]] || {
    echo "apple blind self-WGDI script differs from execution freeze" >&2; exit 1;
}
for relative in src/ploidypatch/self_wgd_pairs.py src/ploidypatch/wgd_candidate_select.py; do
    expected=$(awk -F '\t' -v path="$relative" '$1 == path {print $2}' "$protocol_root/code_manifest.tsv")
    observed=$(sha256sum "$code_root/$relative" | awk '{print $1}')
    [[ -n $expected && $observed == "$expected" ]] || {
        echo "post-freeze module change detected: $relative" >&2; exit 1;
    }
done
[[ ! -e $result_root && ! -e $working_root ]] || {
    echo "refusing to overwrite apple blind self-WGD" >&2; exit 1;
}
mkdir -p "$working_root"/{input,db,blast,config,collinearity,pairs,selected,logs}
{
    printf 'field\tvalue\ncode_commit\t%s\n' \
        "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'mode\tblind\nsplit\tuntouched_external_v0.3\n'
    printf 'candidate_stage_truth_access\tfalse\npreperturbation_pair_access\tfalse\n'
    printf 'candidate_source\tchain_preserving_method_union\n'
    printf 'wgd_event\tmalus_lineage_wgd_blind_recomputed\nmin_block_pairs\t20\n'
    printf 'require_different_seqids\ttrue\nrequire_reciprocal_unique\ttrue\n'
    printf 'candidate_candidate_pair_policy\treject_as_circular\n'
    printf 'diamond_threads\t64\nwgdi_processes\t64\nautomatic_approval\tfalse\n'
    printf 'protocol_freeze_sha256sums_sha256\t%s\n' \
        "$(sha256sum "$protocol_root/SHA256SUMS" | awk '{print $1}')"
    printf 'identifier_alias_adapter_sha256\t%s\n' \
        "$(sha256sum "$code_root/scripts/build_wgdi_source_alias_gff.py" | awk '{print $1}')"
} > "$working_root/run_contract.tsv"
{
    printf 'role\tbytes\tsha256\tpath\n'
    for entry in "base_gff:$base_gff" "candidate_gff:$candidate_gff" \
                 "genome:$genome" "fai:$fai"; do
        role=${entry%%:*}; path=${entry#*:}
        printf '%s\t%s\t%s\t%s\n' "$role" "$(stat -Lc %s "$path")" \
            "$(sha256sum "$path" | awk '{print $1}')" "$path"
    done
} > "$working_root/input_manifest.tsv"

/usr/bin/time -v -o "$working_root/logs/gffread_pep.time.txt" \
    conda run -p "$syngap_env" --no-capture-output \
    gffread "$candidate_gff" -g "$genome" \
    -y "$working_root/input/$prefix.all.pep.fa" -S \
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
/usr/bin/time -v -o "$working_root/logs/pair_inference.time.txt" \
    "$python_bin" -m ploidypatch.cli evidence infer-self-wgd-pairs \
    --query-wgdi-gff "$working_root/input/$prefix.wgdi.gff" \
    --collinearity "$working_root/collinearity/${prefix}_self.tsv" \
    --source-gff "$source_alias_gff" \
    --wgd-event malus_lineage_wgd_blind_recomputed --min-block-pairs 20 \
    --output-pairs "$working_root/pairs/${prefix}.self_wgd_pairs.tsv" \
    --decisions-tsv "$working_root/pairs/decisions.tsv" \
    > "$working_root/pairs/stdout.json" 2> "$working_root/pairs/stderr.log"
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
    [[ -s $output ]] || { echo "missing apple blind self-WGD output: $output" >&2; exit 1; }
done
(
    cd "$working_root"
    find . -type f \( -name '*.json' -o -name '*.tsv' -o -name '*.gff' \
        -o -name '*.gff3' -o -name '*.fa' \) -print0 \
        | sort -z | xargs -0 sha256sum > SHA256SUMS
    sha256sum -c SHA256SUMS >/dev/null
)
du -sb "$working_root" > "$working_root/disk_bytes.txt"
mv "$working_root" "$result_root"
printf 'apple blind candidate self-WGD frozen: %s\n' "$result_root"
