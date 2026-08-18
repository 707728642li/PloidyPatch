#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then echo "usage: $0 PROJECT_ROOT" >&2; exit 2; fi
project_root=$(realpath "$1")
code_root=$project_root/code
python_bin=$project_root/envs/ploidypatch-dev/bin/python
method_root=$project_root/results/copy_collapse/holdout/cotton_method_trio_v0.1
transfer_root=$project_root/results/copy_collapse/holdout/cotton_glycine_model_transfer_v0.1
result_root=$project_root/results/copy_collapse/statistics/cotton_zero_retuning_comparison_v0.1
working_root=${result_root}.working

for frozen in "$method_root" "$transfer_root"; do
    if [[ ! -s $frozen/SHA256SUMS ]]; then echo "unfrozen cotton input: $frozen" >&2; exit 1; fi
    (cd "$frozen" && sha256sum -c SHA256SUMS >/dev/null)
done
declare -A scores=(
    [miniprot]="$method_root/methods/miniprot/score.json"
    [gemoma]="$method_root/methods/gemoma/score.json"
    [lifton]="$method_root/methods/lifton/score.json"
    [union]="$method_root/consensus/union/score.json"
    [consensus2]="$method_root/consensus/support2/score.json"
    [consensus3]="$method_root/consensus/support3/score.json"
    [glycine_rank_model]="$transfer_root/evaluator/paired_score.json"
    [model_mask_wgd]="$transfer_root/evaluator/paired_score.mask_wgd_context.json"
    [model_mask_method_quality]="$transfer_root/evaluator/paired_score.mask_method_quality.json"
)
for required in "$python_bin" "${scores[@]}"; do
    if [[ ! -s $required ]]; then echo "missing cotton statistical input: $required" >&2; exit 1; fi
done
if [[ -e $result_root || -e $working_root ]]; then
    echo "refusing to overwrite cotton statistics" >&2; exit 1
fi
mkdir -p "$working_root"
{
    printf 'field\tvalue\n'
    printf 'code_commit\t%s\n' "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'benchmark\tcotton_annotation_copy_collapse_v0.1\n'
    printf 'split\texternal_zero_retuning_holdout\n'
    printf 'formal_holdout_claim_allowed\ttrue\n'
    printf 'candidate_freeze_precedes_truth_scoring\ttrue\n'
    printf 'replicates\t20000\nseed\t20260807\nalpha\t0.05\n'
    printf 'design\tpaired_events_stratified_by_event_type\n'
    printf 'primary_metric\tcomplete_cds_chain_recovery\n'
    printf 'glycine_model_parameter_retuning\tfalse\n'
} > "$working_root/run_contract.tsv"
{
    printf 'label\tbytes\tsha256\tpath\n'
    for label in miniprot gemoma lifton union consensus2 consensus3 \
                 glycine_rank_model model_mask_wgd model_mask_method_quality; do
        path=${scores[$label]}
        printf '%s\t%s\t%s\t%s\n' "$label" "$(stat -Lc %s "$path")" \
            "$(sha256sum "$path" | awk '{print $1}')" "$path"
    done
} > "$working_root/input_manifest.tsv"

cd "$code_root"
score_args=()
for label in miniprot gemoma lifton union consensus2 consensus3 \
             glycine_rank_model model_mask_wgd model_mask_method_quality; do
    score_args+=(--score "$label=${scores[$label]}")
done
confusion_score_args=()
# A masked-feature counterfactual can legitimately select zero candidates at
# the frozen threshold, leaving precision and F1 undefined. Keep those arms in
# the paired event analysis, but exclude them from the confusion bootstrap.
for label in miniprot gemoma lifton union consensus2 consensus3 \
             glycine_rank_model; do
    confusion_score_args+=(--score "$label=${scores[$label]}")
done
"$python_bin" -m ploidypatch.cli benchmark bootstrap-events \
    "${score_args[@]}" --output-json "$working_root/event_recovery.json" \
    --metric complete_cds_chain_recovery \
    --replicates 20000 --seed 20260807 --alpha 0.05 \
    > "$working_root/event_recovery.stdout.json"
"$python_bin" -m ploidypatch.cli benchmark bootstrap-confusion \
    "${confusion_score_args[@]}" --output-json "$working_root/strict_cds_confusion.json" \
    --section strict_cds_chain \
    --replicates 20000 --seed 20260807 --alpha 0.05 \
    > "$working_root/strict_cds_confusion.stdout.json"
for output in "$working_root/event_recovery.json" \
              "$working_root/strict_cds_confusion.json"; do
    if [[ ! -s $output ]] || ! grep -q '"schema_version"' "$output"; then
        echo "cotton statistical output failed validation: $output" >&2
        exit 1
    fi
done
(
    cd "$working_root"
    find . -type f \( -name '*.json' -o -name '*.tsv' \) \
        -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS
    sha256sum -c SHA256SUMS >/dev/null
)
du -sb "$working_root" > "$working_root/disk_bytes.txt"
mv "$working_root" "$result_root"
printf 'cotton zero-retuning statistics frozen: %s\n' "$result_root"
