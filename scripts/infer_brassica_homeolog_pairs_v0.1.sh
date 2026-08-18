#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "usage: $0 PROJECT_ROOT" >&2
    exit 2
fi

project_root=$(realpath "$1")
code_root=$project_root/code
python_bin=$project_root/envs/ploidypatch-dev/bin/python
evidence_root=$project_root/results/evidence/wgdi/brassica_v0.1
gene_evidence=$evidence_root/summary/bna_daae.gene_evidence.tsv
bra_collinearity=$evidence_root/collinearity/bna_daae_vs_bra_a.collinearity.tsv
bol_collinearity=$evidence_root/collinearity/bna_daae_vs_bol_c.collinearity.tsv
result_root=$project_root/results/evidence/homeolog_pairs/brassica_v0.1
working_root=${result_root}.working

for required in "$python_bin" "$gene_evidence" "$bra_collinearity" \
                "$bol_collinearity"; do
    if [[ ! -s $required ]]; then
        echo "missing or empty homeolog-pair input: $required" >&2
        exit 1
    fi
done
if [[ -e $result_root || -e $working_root ]]; then
    echo "refusing to overwrite homeolog-pair result: $result_root" >&2
    exit 1
fi

mkdir -p "$working_root"
{
    printf 'field\tvalue\n'
    printf 'code_commit\t%s\n' "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'wgd_event\tbrassica_napus_allopolyploid_AC\n'
    printf 'subgenomes\tA,C\n'
    printf 'min_support_group_count\t2\n'
    printf 'support_unit\tindependent_reference_species\n'
    printf 'within_source_policy\texactly_one_gene_per_subgenome\n'
    printf 'across_source_policy\treciprocal_unique_partner\n'
    printf 'benchmark_access\tevaluator_only\n'
    printf 'homeolog_module_sha256\t%s\n' \
        "$(sha256sum "$code_root/src/ploidypatch/homeolog_pairs.py" | awk '{print $1}')"
} > "$working_root/run_contract.tsv"
{
    printf 'role\tbytes\tsha256\tpath\n'
    for entry in \
        "gene_evidence:$gene_evidence" \
        "bra_a_collinearity:$bra_collinearity" \
        "bol_c_collinearity:$bol_collinearity"; do
        role=${entry%%:*}
        path=${entry#*:}
        printf '%s\t%s\t%s\t%s\n' "$role" "$(stat -Lc %s "$path")" \
            "$(sha256sum "$path" | awk '{print $1}')" "$path"
    done
} > "$working_root/input_manifest.tsv"

cd "$code_root"
/usr/bin/time -v -o "$working_root/resource.time.txt" \
    "$python_bin" -m ploidypatch.cli evidence infer-homeolog-pairs \
        --gene-evidence "$gene_evidence" \
        --collinearity "bra_a=$bra_collinearity" \
        --collinearity "bol_c=$bol_collinearity" \
        --wgd-event brassica_napus_allopolyploid_AC \
        --subgenome A \
        --subgenome C \
        --min-support-group-count 2 \
        --output-pairs "$working_root/bna_daae.AC.homeolog_pairs.tsv" \
        --decisions-tsv "$working_root/decisions.tsv" \
        > "$working_root/stdout.json" \
        2> "$working_root/stderr.log"

for output in "$working_root/bna_daae.AC.homeolog_pairs.tsv" \
              "$working_root/bna_daae.AC.homeolog_pairs.tsv.manifest.json" \
              "$working_root/decisions.tsv" "$working_root/stdout.json"; do
    if [[ ! -s $output ]]; then
        echo "missing or empty homeolog-pair output: $output" >&2
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
printf 'homeolog pairs frozen: %s\n' "$result_root"
