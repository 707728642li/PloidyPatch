#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "usage: $0 PROJECT_ROOT" >&2
    exit 2
fi

project_root=$(realpath "$1")
code_root=$project_root/code
python_bin=$project_root/envs/ploidypatch-dev/bin/python
source_gff=$project_root/data/derived/syngap_glycine_v0.1/gma_complete/primary_chromosomes.gff3
pair_root=$project_root/results/evidence/wgdi/glycine_self_v0.2/pairs
pair_tsv=$pair_root/gma_v21.self_wgd_pairs.tsv
pair_manifest=${pair_tsv}.manifest.json
count=800
seed=20260815
result_root=$project_root/benchmark/structure/copy_collapse_v0.1/gma_v21/annotation_copy_collapse_seed${seed}
working_root=${result_root}.working

for required in "$python_bin" "$source_gff" "$pair_tsv" "$pair_manifest"; do
    if [[ ! -s $required ]]; then
        echo "missing or empty Glycine copy-collapse prerequisite: $required" >&2
        exit 1
    fi
done
if [[ -e $result_root || -e $working_root ]]; then
    echo "refusing to overwrite Glycine copy-collapse benchmark: $result_root" >&2
    exit 1
fi

mkdir -p "$working_root/blind" "$working_root/evaluator/truth" \
    "$working_root/evaluator/validation"
{
    printf 'field\tvalue\n'
    printf 'code_commit\t%s\n' "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'dataset_id\tgma_v21\n'
    printf 'event_type\tannotation_copy_collapse\n'
    printf 'count\t%s\n' "$count"
    printf 'seed\t%s\n' "$seed"
    printf 'split\theldout\n'
    printf 'pair_access\tevaluator_only\n'
    printf 'pair_policy\tself_wgdi_cross_seqid_reciprocal_unique_block20_feature_id_v2\n'
    printf 'pair_claim_boundary\tcandidate_named_wgd_syntenic_duplicate_not_gene_tree_homeology\n'
    printf 'source_gff_sha256\t%s\n' "$(sha256sum "$source_gff" | awk '{print $1}')"
    printf 'pair_tsv_sha256\t%s\n' "$(sha256sum "$pair_tsv" | awk '{print $1}')"
    printf 'pair_manifest_sha256\t%s\n' "$(sha256sum "$pair_manifest" | awk '{print $1}')"
    printf 'structure_module_sha256\t%s\n' \
        "$(sha256sum "$code_root/src/ploidypatch/structure_perturb.py" | awk '{print $1}')"
    printf 'score_module_sha256\t%s\n' \
        "$(sha256sum "$code_root/src/ploidypatch/score.py" | awk '{print $1}')"
} > "$working_root/run_contract.tsv"

cd "$code_root"
/usr/bin/time -v -o "$working_root/perturb.resource.time.txt" \
    "$python_bin" -m ploidypatch.cli benchmark perturb \
        --gff "$source_gff" \
        --output-dir "$working_root/blind" \
        --truth-dir "$working_root/evaluator/truth" \
        --event-type annotation_copy_collapse \
        --pair-tsv "$pair_tsv" \
        --count "$count" \
        --seed "$seed" \
        > "$working_root/perturb.stdout.json" \
        2> "$working_root/perturb.stderr.log"

perturbed=$working_root/blind/perturbed.gff3
truth=$working_root/evaluator/truth/hidden_truth.json
restored=$working_root/evaluator/restored.gff3
/usr/bin/time -v -o "$working_root/restore.resource.time.txt" \
    "$python_bin" -m ploidypatch.cli benchmark restore \
        --perturbed-gff "$perturbed" \
        --truth "$truth" \
        --output-gff "$restored" \
        > "$working_root/evaluator/restoration_report.json" \
        2> "$working_root/evaluator/restoration.stderr.log"
if ! cmp -s "$source_gff" "$restored"; then
    echo "Glycine copy-collapse restoration is not byte-identical to source" >&2
    exit 1
fi

for mode in noop oracle; do
    if [[ $mode == noop ]]; then
        candidate=$perturbed
    else
        candidate=$source_gff
    fi
    /usr/bin/time -v -o "$working_root/evaluator/validation/${mode}.resource.time.txt" \
        "$python_bin" -m ploidypatch.cli benchmark score \
            --source-gff "$source_gff" \
            --perturbed-gff "$perturbed" \
            --candidate-gff "$candidate" \
            --truth "$truth" \
            --include-event-details \
            > "$working_root/evaluator/validation/score_${mode}.json" \
            2> "$working_root/evaluator/validation/score_${mode}.stderr.log"
done

if ! grep -q '"complete_cds_chain_recovery": 0' \
        "$working_root/evaluator/validation/score_noop.json"; then
    echo "Glycine copy-collapse no-op CDS sentinel failed" >&2
    exit 1
fi
if ! grep -q '"complete_cds_chain_recovery": 800' \
        "$working_root/evaluator/validation/score_oracle.json"; then
    echo "Glycine copy-collapse oracle CDS sentinel failed" >&2
    exit 1
fi
for score in "$working_root/evaluator/validation/score_noop.json" \
             "$working_root/evaluator/validation/score_oracle.json"; do
    if ! grep -q '"grade": "pass"' "$score"; then
        echo "Glycine copy-collapse evaluator quality gate failed: $score" >&2
        exit 1
    fi
done

for output in "$perturbed" "$working_root/blind/manifest.json" "$truth" \
              "$restored" "$working_root/evaluator/restoration_report.json" \
              "$working_root/evaluator/validation/score_noop.json" \
              "$working_root/evaluator/validation/score_oracle.json"; do
    if [[ ! -s $output ]]; then
        echo "missing or empty Glycine copy-collapse output: $output" >&2
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
printf 'Glycine copy-collapse holdout frozen: %s\n' "$result_root"
