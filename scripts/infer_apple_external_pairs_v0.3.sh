#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then echo "usage: $0 PROJECT_ROOT" >&2; exit 2; fi
project_root=$(realpath "$1")
code_root=$project_root/code
python_bin=$project_root/envs/ploidypatch-dev/bin/python
input_root=$project_root/data/derived/external_evaluator/apple_v0.3_wgdi_inputs
normalized_root=$project_root/data/derived/external_inputs/apple_v0.3
wgdi_root=$project_root/results/evaluator/apple_v0.3/wgdi
protocol_root=$project_root/results/protocol_freezes/apple_external_v0.3
execution_root=$project_root/results/protocol_freezes/apple_external_v0.3_execution
result_root=$project_root/results/evaluator/apple_v0.3/truth_pairs
working_root=${result_root}.working
query_gff=$input_root/mdx/mdx.wgdi.gff
representatives=$input_root/mdx/mdx.representatives.tsv
source_gff=$normalized_root/target_apple/primary_chromosomes.gff3
self_collinearity=$wgdi_root/collinearity/mdx_self.tsv
rose_collinearity=$wgdi_root/collinearity/mdx_vs_rch.tsv
strawberry_collinearity=$wgdi_root/collinearity/mdx_vs_fve.tsv

for required in "$python_bin" "$input_root/SHA256SUMS" "$wgdi_root/SHA256SUMS" \
    "$protocol_root/SHA256SUMS" "$execution_root/SHA256SUMS" "$query_gff" \
    "$representatives" "$source_gff" \
    "$self_collinearity" "$rose_collinearity" "$strawberry_collinearity" \
    "$code_root/src/ploidypatch/pair_consensus.py" \
    "$code_root/scripts/build_wgdi_source_alias_gff.py"; do
    [[ -s $required ]] || { echo "missing apple pair prerequisite: $required" >&2; exit 1; }
done
(cd "$input_root" && sha256sum -c SHA256SUMS >/dev/null)
(cd "$wgdi_root" && sha256sum -c SHA256SUMS >/dev/null)
(cd "$protocol_root" && sha256sum -c SHA256SUMS >/dev/null)
(cd "$execution_root" && sha256sum -c SHA256SUMS >/dev/null)
expected_self=$(awk -F '\t' '$1 == "scripts/infer_apple_external_pairs_v0.3.sh" {print $3}' \
    "$execution_root/implementation_manifest.tsv")
[[ -n $expected_self && $(sha256sum "$code_root/scripts/infer_apple_external_pairs_v0.3.sh" | awk '{print $1}') == "$expected_self" ]] || {
    echo "apple pair script differs from execution freeze" >&2; exit 1;
}
[[ ! -e $result_root && ! -e $working_root ]] || {
    echo "refusing to overwrite apple truth pairs" >&2; exit 1;
}
mkdir -p "$working_root"/{mapping,self_wgdi,two_outgroups,intersection}
{
    printf 'field\tvalue\n'
    printf 'code_commit\t%s\n' "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'access\tevaluator_only\nhidden_event_generation\tfalse\n'
    printf 'candidate_reference_access\tfalse\nexternal_label_access\tfalse\n'
    printf 'final_rule\texact_unordered_pair_intersection_of_self_wgdi_and_two_outgroup_support\n'
    printf 'min_block_pairs\t20\noutgroup_min_support_groups\t2\n'
    printf 'require_cross_seqid\ttrue\nrequire_reciprocal_unique\ttrue\n'
    printf 'protocol_freeze_sha256sums_sha256\t%s\n' \
        "$(sha256sum "$protocol_root/SHA256SUMS" | awk '{print $1}')"
    printf 'pair_intersection_module_sha256\t%s\n' \
        "$(sha256sum "$code_root/src/ploidypatch/pair_consensus.py" | awk '{print $1}')"
    printf 'identifier_alias_adapter_sha256\t%s\n' \
        "$(sha256sum "$code_root/scripts/build_wgdi_source_alias_gff.py" | awk '{print $1}')"
} > "$working_root/run_contract.tsv"

cd "$code_root"
source_alias_gff=$working_root/mapping/mdx.source_alias.gff3
"$python_bin" scripts/build_wgdi_source_alias_gff.py \
    --source-gff "$source_gff" --representatives "$representatives" \
    --output-gff "$source_alias_gff" \
    > "$working_root/mapping/stdout.json" \
    2> "$working_root/mapping/stderr.log"
"$python_bin" -m ploidypatch.cli evidence infer-self-wgd-pairs \
    --query-wgdi-gff "$query_gff" --collinearity "$self_collinearity" \
    --source-gff "$source_alias_gff" --wgd-event malus_lineage_wgd_self \
    --min-block-pairs 20 \
    --output-pairs "$working_root/self_wgdi/pairs.tsv" \
    --decisions-tsv "$working_root/self_wgdi/decisions.tsv" \
    > "$working_root/self_wgdi/stdout.json" \
    2> "$working_root/self_wgdi/stderr.log"
"$python_bin" -m ploidypatch.cli evidence infer-outgroup-duplicated-pairs \
    --query-wgdi-gff "$query_gff" --source-gff "$source_alias_gff" \
    --collinearity "rose=$rose_collinearity" \
    --collinearity "strawberry=$strawberry_collinearity" \
    --wgd-event malus_lineage_wgd_two_outgroups \
    --min-support-group-count 2 --min-block-pairs 20 \
    --output-pairs "$working_root/two_outgroups/pairs.tsv" \
    --decisions-tsv "$working_root/two_outgroups/decisions.tsv" \
    > "$working_root/two_outgroups/stdout.json" \
    2> "$working_root/two_outgroups/stderr.log"
"$python_bin" -m ploidypatch.cli evidence intersect-copy-pair-evidence \
    --pairs "self_wgdi=$working_root/self_wgdi/pairs.tsv" \
    --pairs "two_outgroups=$working_root/two_outgroups/pairs.tsv" \
    --pair-set-label malus_wgd_self_and_two_outgroups \
    --output-pairs "$working_root/intersection/pairs.tsv" \
    --decisions-tsv "$working_root/intersection/decisions.tsv" \
    > "$working_root/intersection/stdout.json" \
    2> "$working_root/intersection/stderr.log"

for path in self_wgdi/pairs.tsv two_outgroups/pairs.tsv intersection/pairs.tsv; do
    [[ -s $working_root/$path ]] || { echo "missing apple pair output: $path" >&2; exit 1; }
done
{
    printf 'pair_set\taccepted_pairs\n'
    for label in self_wgdi two_outgroups intersection; do
        printf '%s\t%s\n' "$label" "$(( $(wc -l < "$working_root/$label/pairs.tsv") - 1 ))"
    done
} > "$working_root/pair_counts.tsv"
(
    cd "$working_root"
    find . -type f \( -name '*.tsv' -o -name '*.json' -o -name '*.gff3' \) -print0 \
        | sort -z | xargs -0 sha256sum > SHA256SUMS
    sha256sum -c SHA256SUMS >/dev/null
)
mv "$working_root" "$result_root"
printf 'apple evaluator truth pairs frozen: %s\n' "$result_root"
