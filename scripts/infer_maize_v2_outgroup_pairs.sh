#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then echo "usage: $0 PROJECT_ROOT" >&2; exit 2; fi
project_root=$(realpath "$1")
code_root=$project_root/code
python_bin=$project_root/envs/ploidypatch-dev/bin/python
input_root=$project_root/data/derived/holdout_evaluator/maize_v2_wgdi_inputs
wgdi_root=$project_root/results/evaluator/maize_v2/wgdi_outgroups
policy=$code_root/config/maize_v2_zero_retuning_policy.tsv
result_root=$project_root/results/evaluator/maize_v2/outgroup_duplicated_pairs
working_root=${result_root}.working
query_gff=$input_root/zma/zma.wgdi.gff
source_gff=$project_root/data/derived/holdout_inputs/maize_v2/zea_mays/primary_chromosomes.gff3
sorghum=$wgdi_root/collinearity/zma_vs_sbi.tsv
setaria=$wgdi_root/collinearity/zma_vs_sit.tsv

for required in "$python_bin" "$input_root/SHA256SUMS" "$wgdi_root/SHA256SUMS" \
                "$query_gff" "$source_gff" "$sorghum" "$setaria" "$policy"; do
    [[ -s $required ]] || { echo "missing maize pair prerequisite: $required" >&2; exit 1; }
done
(cd "$input_root" && sha256sum -c SHA256SUMS >/dev/null)
(cd "$wgdi_root" && sha256sum -c SHA256SUMS >/dev/null)
if [[ -e $result_root || -e $working_root ]]; then
    echo "refusing to overwrite maize outgroup pairs" >&2; exit 1
fi
mkdir -p "$working_root"
{
    printf 'field\tvalue\n'
    printf 'code_commit\t%s\n' "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'access\tevaluator_only\n'
    printf 'wgd_event\tzea_mays_lineage_tetraploidy\n'
    printf 'support_unit\tindependent_outgroup_species\n'
    printf 'min_support_group_count\t2\n'
    printf 'min_block_pairs\t20\n'
    printf 'pair_policy\texactly_two_cross_chromosome_maize_genes_per_outgroup_counterpart_then_reciprocal_unique\n'
    printf 'policy_sha256\t%s\n' "$(sha256sum "$policy" | awk '{print $1}')"
} > "$working_root/run_contract.tsv"

cd "$code_root"
"$python_bin" -m ploidypatch.cli evidence infer-outgroup-duplicated-pairs \
    --query-wgdi-gff "$query_gff" \
    --source-gff "$source_gff" \
    --collinearity "sorghum=$sorghum" \
    --collinearity "setaria=$setaria" \
    --wgd-event zea_mays_lineage_tetraploidy \
    --min-support-group-count 2 --min-block-pairs 20 \
    --output-pairs "$working_root/zma.outgroup_duplicated_pairs.tsv" \
    --decisions-tsv "$working_root/decisions.tsv" \
    > "$working_root/stdout.json" 2> "$working_root/stderr.log"

for output in zma.outgroup_duplicated_pairs.tsv \
              zma.outgroup_duplicated_pairs.tsv.manifest.json decisions.tsv; do
    [[ -s $working_root/$output ]] || { echo "missing maize pair artifact: $output" >&2; exit 1; }
done
(
    cd "$working_root"
    find . -type f \( -name '*.tsv' -o -name '*.json' \) -print0 \
        | sort -z | xargs -0 sha256sum > SHA256SUMS
    sha256sum -c SHA256SUMS >/dev/null
)
du -sb "$working_root" > "$working_root/disk_bytes.txt"
mv "$working_root" "$result_root"
printf 'maize evaluator duplicated pairs frozen: %s\n' "$result_root"
