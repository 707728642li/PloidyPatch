#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then echo "usage: $0 PROJECT_ROOT" >&2; exit 2; fi
project_root=$(realpath "$1")
code_root=$project_root/code
feature_root=$project_root/results/copy_collapse/model_development/homeolog_topology_v0.2
result_root=$project_root/results/copy_collapse/model_development/homeolog_topology_nested_evaluation_v0.2
working_root=${result_root}.working
model_env=$project_root/envs/ploidypatch-model

[[ -s $feature_root/SHA256SUMS ]] || { echo "unfrozen topology features" >&2; exit 1; }
(cd "$feature_root" && sha256sum -c SHA256SUMS >/dev/null)
if [[ -e $result_root || -e $working_root ]]; then
    echo "refusing to overwrite topology evaluation" >&2
    exit 1
fi
mkdir -p "$working_root"

gly_labels=$project_root/results/copy_collapse/model_development/glycine_feature_matrix_v0.1/evaluator/labeled_features.tsv
gly_scores=$project_root/results/copy_collapse/model_development/glycine_copy_model_evaluation_v0.1/blind/scores.tsv
bra_labels=$project_root/results/copy_collapse/model_development/brassica_glycine_model_transfer_v0.1/evaluator/labeled_features.tsv
bra_scores=$project_root/results/copy_collapse/model_development/brassica_glycine_model_transfer_v0.1/blind/scores.tsv
cot_labels=$project_root/results/copy_collapse/holdout/cotton_glycine_model_transfer_v0.1/evaluator/labeled_features.tsv
cot_scores=$project_root/results/copy_collapse/holdout/cotton_glycine_model_transfer_v0.1/blind/scores.tsv

{
    printf 'field\tvalue\n'
    printf 'code_commit\t%s\n' "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'development_species\tGlycine_max,Brassica_napus\n'
    printf 'diagnostic_species\tGossypium_hirsutum\n'
    printf 'cotton_used_for_fitting\tfalse\n'
    printf 'validation\tfive_fold_stratified_chromosome_grouped_five_seeds\n'
    printf 'bootstrap\t2000_chromosome_group_replicates\n'
    printf 'ablation\tbaseline,normalized_WGD_context,topology,combined\n'
    printf 'stacking_on_in_sample_model_scores\tfalse\n'
    printf 'feature_contract_fit_scope\ttraining_fold_or_training_species_only\n'
    printf 'automatic_approval\tfalse\n'
} > "$working_root/run_contract.tsv"
{
    printf 'species\trole\tbytes\tsha256\tpath\n'
    for entry in \
        "glycine:topology:$feature_root/glycine/features.tsv" \
        "glycine:labels:$gly_labels" "glycine:scores:$gly_scores" \
        "brassica:topology:$feature_root/brassica/features.tsv" \
        "brassica:labels:$bra_labels" "brassica:scores:$bra_scores" \
        "cotton:topology:$feature_root/cotton/features.tsv" \
        "cotton:labels:$cot_labels" "cotton:scores:$cot_scores"; do
        species=${entry%%:*}; rest=${entry#*:}; role=${rest%%:*}; path=${rest#*:}
        [[ -s $path ]] || { echo "missing topology evaluation input: $path" >&2; exit 1; }
        printf '%s\t%s\t%s\t%s\t%s\n' "$species" "$role" \
            "$(stat -Lc %s "$path")" "$(sha256sum "$path" | awk '{print $1}')" \
            "$path"
    done
} > "$working_root/input_manifest.tsv"

cd "$code_root"
/usr/bin/time -v -o "$working_root/resource.time.txt" \
    env PYTHONPATH="$code_root/src" \
    conda run -p "$model_env" --no-capture-output python \
    scripts/evaluate_homeolog_topology_v0.2.py \
    --dataset "glycine=development,$feature_root/glycine/features.tsv,$gly_labels,$gly_scores" \
    --dataset "brassica=development,$feature_root/brassica/features.tsv,$bra_labels,$bra_scores" \
    --dataset "cotton=post_holdout_diagnostic,$feature_root/cotton/features.tsv,$cot_labels,$cot_scores" \
    --output-json "$working_root/evaluation.json" \
    > "$working_root/stdout.json" 2> "$working_root/stderr.log"
[[ -s $working_root/evaluation.json ]] || { echo "missing topology evaluation" >&2; exit 1; }
(
    cd "$working_root"
    find . -type f \( -name '*.json' -o -name '*.tsv' \) -print0 \
        | sort -z | xargs -0 sha256sum > SHA256SUMS
    sha256sum -c SHA256SUMS >/dev/null
)
du -sb "$working_root" > "$working_root/disk_bytes.txt"
mv "$working_root" "$result_root"
printf 'homeolog topology evaluation frozen: %s\n' "$result_root"
