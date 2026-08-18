#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then echo "usage: $0 PROJECT_ROOT" >&2; exit 2; fi
project_root=$(realpath "$1")
code_root=$project_root/code
python_bin=$project_root/envs/ploidypatch-dev/bin/python
input_root=$project_root/data/derived/holdout_evaluator/cotton_wgdi_inputs_v0.1
wgdi_root=$project_root/results/evaluator/cotton_holdout_v0.1/wgdi
result_root=$project_root/results/evaluator/cotton_holdout_v0.1/homeolog_pairs
working_root=${result_root}.working
query_gff=$input_root/ghi_ad/ghi_ad.wgdi.gff
query_manifest=$project_root/data/derived/holdout_inputs/cotton_v0.1/hirsutum/manifest.json
gar=$wgdi_root/collinearity/ghi_ad_vs_gar_a.tsv
gra=$wgdi_root/collinearity/ghi_ad_vs_gra_d.tsv

for required in "$python_bin" "$wgdi_root/SHA256SUMS" "$query_gff" "$query_manifest" "$gar" "$gra"; do
    if [[ ! -s $required ]]; then echo "missing cotton pair prerequisite: $required" >&2; exit 1; fi
done
(cd "$wgdi_root" && sha256sum -c SHA256SUMS >/dev/null)
if [[ -e $result_root || -e $working_root ]]; then echo "refusing to overwrite cotton pairs" >&2; exit 1; fi
mkdir -p "$working_root"
{
    printf 'field\tvalue\ncode_commit\t%s\n' "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'access\tevaluator_only\nwgd_event\tgossypium_hirsutum_allopolyploid_AD\n'
    printf 'support_unit\tindependent_diploid_reference_species\nmin_support_group_count\t2\n'
    printf 'pair_policy\texactly_one_A_and_one_D_per_counterpart_then_reciprocal_unique\n'
} > "$working_root/run_contract.tsv"
cd "$code_root"
"$python_bin" -m ploidypatch.cli evidence summarize-wgdi \
    --query-gff "$query_gff" --query-input-manifest "$query_manifest" \
    --collinearity "gar_a=$gar" --collinearity "gra_d=$gra" \
    --expected-source A=gar_a --expected-source D=gra_d \
    --output-gene-tsv "$working_root/gene_evidence.tsv" \
    --output-block-tsv "$working_root/blocks.tsv" \
    > "$working_root/summary.stdout.json" 2> "$working_root/summary.stderr.log"
"$python_bin" -m ploidypatch.cli evidence infer-homeolog-pairs \
    --gene-evidence "$working_root/gene_evidence.tsv" \
    --collinearity "gar_a=$gar" --collinearity "gra_d=$gra" \
    --wgd-event gossypium_hirsutum_allopolyploid_AD --subgenome A --subgenome D \
    --min-support-group-count 2 \
    --output-pairs "$working_root/ghi_ad.AD.homeolog_pairs.tsv" \
    --decisions-tsv "$working_root/decisions.tsv" \
    > "$working_root/pairs.stdout.json" 2> "$working_root/pairs.stderr.log"
for output in gene_evidence.tsv blocks.tsv ghi_ad.AD.homeolog_pairs.tsv \
              ghi_ad.AD.homeolog_pairs.tsv.manifest.json decisions.tsv; do
    if [[ ! -s $working_root/$output ]]; then echo "missing cotton pair result: $output" >&2; exit 1; fi
done
(
    cd "$working_root"
    find . -type f \( -name '*.tsv' -o -name '*.json' \) -print0 | sort -z \
        | xargs -0 sha256sum > SHA256SUMS
    sha256sum -c SHA256SUMS >/dev/null
)
du -sb "$working_root" > "$working_root/disk_bytes.txt"
mv "$working_root" "$result_root"
printf 'cotton evaluator homeolog pairs frozen: %s\n' "$result_root"
