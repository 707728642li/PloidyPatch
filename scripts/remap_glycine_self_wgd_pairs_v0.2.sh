#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "usage: $0 PROJECT_ROOT" >&2
    exit 2
fi

project_root=$(realpath "$1")
code_root=$project_root/code
python_bin=$project_root/envs/ploidypatch-dev/bin/python
source_root=$project_root/results/evidence/wgdi/glycine_self_v0.1
query_wgdi_gff=$source_root/input/gma_v21.wgdi.gff
collinearity=$source_root/collinearity/gma_v21_self.tsv
source_gff=$project_root/data/derived/syngap_glycine_v0.1/gma_complete/primary_chromosomes.gff3
result_root=$project_root/results/evidence/wgdi/glycine_self_v0.2
working_root=${result_root}.working

for required in "$python_bin" "$query_wgdi_gff" "$collinearity" "$source_gff"; do
    if [[ ! -s $required ]]; then
        echo "missing or empty Glycine self-WGD remap input: $required" >&2
        exit 1
    fi
done
if [[ -e $result_root || -e $working_root ]]; then
    echo "refusing to overwrite Glycine self-WGD remap: $result_root" >&2
    exit 1
fi

mkdir -p "$working_root/pairs"
{
    printf 'field\tvalue\n'
    printf 'code_commit\t%s\n' "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'source_evidence_namespace\tglycine_self_v0.1\n'
    printf 'pair_schema\tploidypatch.self_wgd_pairs.v2\n'
    printf 'identifier_policy\tunique_source_gff_ID_gene_id_locus_tag_mapping\n'
    printf 'wgd_event\tglycine_max_paleopolyploid_WGD\n'
    printf 'min_block_pairs\t20\n'
    printf 'require_different_seqids\ttrue\n'
    printf 'require_reciprocal_unique\ttrue\n'
    printf 'benchmark_access\tevaluator_only\n'
    printf 'self_wgd_module_sha256\t%s\n' \
        "$(sha256sum "$code_root/src/ploidypatch/self_wgd_pairs.py" | awk '{print $1}')"
} > "$working_root/run_contract.tsv"

cd "$code_root"
/usr/bin/time -v -o "$working_root/resource.time.txt" \
    "$python_bin" -m ploidypatch.cli evidence infer-self-wgd-pairs \
        --query-wgdi-gff "$query_wgdi_gff" \
        --collinearity "$collinearity" \
        --source-gff "$source_gff" \
        --wgd-event glycine_max_paleopolyploid_WGD \
        --min-block-pairs 20 \
        --output-pairs "$working_root/pairs/gma_v21.self_wgd_pairs.tsv" \
        --decisions-tsv "$working_root/pairs/decisions.tsv" \
        > "$working_root/pairs/stdout.json" \
        2> "$working_root/pairs/stderr.log"

for output in "$working_root/pairs/gma_v21.self_wgd_pairs.tsv" \
              "$working_root/pairs/gma_v21.self_wgd_pairs.tsv.manifest.json" \
              "$working_root/pairs/decisions.tsv" "$working_root/pairs/stdout.json"; do
    if [[ ! -s $output ]]; then
        echo "missing or empty Glycine self-WGD remap output: $output" >&2
        exit 1
    fi
done
(
    cd "$working_root"
    find . -type f \( -name '*.json' -o -name '*.tsv' \) \
        -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS
)
du -sb "$working_root" > "$working_root/disk_bytes.txt"
mv "$working_root" "$result_root"
printf 'Glycine self-WGD feature-ID pairs frozen: %s\n' "$result_root"
