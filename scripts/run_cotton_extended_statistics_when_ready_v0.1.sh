#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then echo "usage: $0 PROJECT_ROOT" >&2; exit 2; fi
project_root=$(realpath "$1")
python_bin=$project_root/envs/ploidypatch-dev/bin/python
method_root=$project_root/results/copy_collapse/holdout/cotton_method_trio_v0.1
model_root=$project_root/results/copy_collapse/holdout/cotton_glycine_model_transfer_v0.1
syngap_root=$project_root/results/copy_collapse/holdout/cotton_syngap_v0.1
post_pid_file=$project_root/logs/copy_collapse/cotton_post_methods_v0.1.pid
result_root=$project_root/results/copy_collapse/statistics/cotton_zero_retuning_extended_comparison_v0.1
working_root=${result_root}.working

while [[ ! -s $model_root/SHA256SUMS ]]; do
    pid=$(cat "$post_pid_file")
    kill -0 "$pid" 2>/dev/null || { echo "cotton post-method pipeline stopped" >&2; exit 1; }
    sleep 30
done
for frozen in "$method_root" "$model_root" "$syngap_root"; do
    [[ -s $frozen/SHA256SUMS ]] || { echo "unfrozen extended input: $frozen" >&2; exit 1; }
    (cd "$frozen" && sha256sum -c SHA256SUMS >/dev/null)
done
declare -A scores=(
    [miniprot]="$method_root/methods/miniprot/score.json"
    [gemoma]="$method_root/methods/gemoma/score.json"
    [lifton]="$method_root/methods/lifton/score.json"
    [union]="$method_root/consensus/union/score.json"
    [consensus2]="$method_root/consensus/support2/score.json"
    [consensus3]="$method_root/consensus/support3/score.json"
    [glycine_rank_model]="$model_root/evaluator/paired_score.json"
    [model_mask_wgd]="$model_root/evaluator/paired_score.mask_wgd_context.json"
    [model_mask_method_quality]="$model_root/evaluator/paired_score.mask_method_quality.json"
    [syngap]="$syngap_root/evaluator/score.json"
)
if [[ -e $result_root || -e $working_root ]]; then echo "extended statistics exists" >&2; exit 1; fi
mkdir -p "$working_root"
{
    printf 'field\tvalue\ncode_commit\t%s\n' "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'split\texternal_zero_retuning_holdout\nformal_holdout_claim_allowed\ttrue\n'
    printf 'replicates\t20000\nseed\t20260807\nalpha\t0.05\n'
    printf 'design\tpaired_events_stratified_by_event_type\n'
    printf 'syngap_role\texternal_comparator_not_primary_policy\n'
} > "$working_root/run_contract.tsv"
score_args=()
{
    printf 'label\tbytes\tsha256\tpath\n'
    for label in miniprot gemoma lifton union consensus2 consensus3 \
        glycine_rank_model model_mask_wgd model_mask_method_quality syngap; do
        path=${scores[$label]}; [[ -s $path ]] || { echo "missing score: $path" >&2; exit 1; }
        score_args+=(--score "$label=$path")
        printf '%s\t%s\t%s\t%s\n' "$label" "$(stat -Lc %s "$path")" \
            "$(sha256sum "$path" | awk '{print $1}')" "$path"
    done
} > "$working_root/input_manifest.tsv"
confusion_score_args=()
# Masked-feature counterfactuals may select zero candidates at the frozen
# threshold. Their precision/F1 are undefined, so use them only in the paired
# event analysis and keep the confusion bootstrap restricted to scored arms.
for label in miniprot gemoma lifton union consensus2 consensus3 \
             glycine_rank_model syngap; do
    confusion_score_args+=(--score "$label=${scores[$label]}")
done
cd "$project_root/code"
"$python_bin" -m ploidypatch.cli benchmark bootstrap-events \
    "${score_args[@]}" --output-json "$working_root/event_recovery.json" \
    --metric complete_cds_chain_recovery --replicates 20000 --seed 20260807 --alpha 0.05 \
    > "$working_root/event_recovery.stdout.json"
"$python_bin" -m ploidypatch.cli benchmark bootstrap-confusion \
    "${confusion_score_args[@]}" --output-json "$working_root/strict_cds_confusion.json" \
    --section strict_cds_chain --replicates 20000 --seed 20260807 --alpha 0.05 \
    > "$working_root/strict_cds_confusion.stdout.json"
(
    cd "$working_root"
    find . -type f \( -name '*.json' -o -name '*.tsv' \) -print0 \
        | sort -z | xargs -0 sha256sum > SHA256SUMS
    sha256sum -c SHA256SUMS >/dev/null
)
du -sb "$working_root" > "$working_root/disk_bytes.txt"
mv "$working_root" "$result_root"
printf 'cotton extended statistics frozen: %s\n' "$result_root"
