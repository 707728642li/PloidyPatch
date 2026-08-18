#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then echo "usage: $0 PROJECT_ROOT" >&2; exit 2; fi
project_root=$(realpath "$1")
code_root=$project_root/code
python_bin=$project_root/envs/ploidypatch-dev/bin/python
benchmark_root=$project_root/benchmark/structure/copy_collapse_v0.1/ghi_ad/annotation_copy_collapse_seed20260817
blind_gff=$benchmark_root/blind/perturbed.gff3
truth=$benchmark_root/evaluator/truth/hidden_truth.json
source_gff=$project_root/data/derived/holdout_inputs/cotton_v0.1/hirsutum/primary_chromosomes.gff3
method_root=$project_root/results/copy_collapse/holdout/cotton_method_trio_v0.1
self_wgd_root=$project_root/results/copy_collapse/holdout/cotton_union_self_wgd_v0.1
model_root=$project_root/results/copy_collapse/model_development/glycine_copy_model_v0.1
model_json=$model_root/model.json
policy=$code_root/config/copy_collapse_zero_retuning_policy_v0.1.tsv
result_root=$project_root/results/copy_collapse/holdout/cotton_glycine_model_transfer_v0.1
working_root=${result_root}.working

for frozen in "$method_root" "$model_root" "$self_wgd_root/blind" \
              "$self_wgd_root/complete_control"; do
    if [[ ! -s $frozen/SHA256SUMS ]]; then echo "unfrozen input: $frozen" >&2; exit 1; fi
    (cd "$frozen" && sha256sum -c SHA256SUMS >/dev/null)
done
for required in "$python_bin" "$blind_gff" "$source_gff" "$model_json" "$policy"; do
    if [[ ! -s $required ]]; then echo "missing cotton transfer input: $required" >&2; exit 1; fi
done
if [[ -e $result_root || -e $working_root ]]; then
    echo "refusing to overwrite cotton model transfer" >&2; exit 1
fi
mkdir -p "$working_root/blind" "$working_root/complete_control" \
    "$working_root/evaluator" \
    "$working_root/counterfactual/wgd_context/blind" \
    "$working_root/counterfactual/wgd_context/complete_control" \
    "$working_root/counterfactual/method_quality/blind" \
    "$working_root/counterfactual/method_quality/complete_control"
{
    printf 'field\tvalue\n'
    printf 'code_commit\t%s\n' "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'source_model_species\tGlycine_max\n'
    printf 'target_species\tGossypium_hirsutum\n'
    printf 'parameter_retuning\tfalse\n'
    printf 'threshold_retuning\tfalse\n'
    printf 'model_sha256\t%s\n' "$(sha256sum "$model_json" | awk '{print $1}')"
    printf 'policy_sha256\t%s\n' "$(sha256sum "$policy" | awk '{print $1}')"
    printf 'split\texternal_zero_retuning_holdout\n'
    printf 'formal_external_holdout_claim_allowed\ttrue\n'
    printf 'candidate_truth_access\tfalse\n'
    printf 'evaluator_truth_access\ttrue\n'
    printf 'paired_complete_annotation_control\ttrue\n'
    printf 'automatic_approval\tfalse\n'
    printf 'selection_policy\treview_only\n'
    printf 'counterfactual_masks\twgd_context,method_quality\n'
    printf 'counterfactual_refitting\tfalse\n'
} > "$working_root/run_contract.tsv"
{
    printf 'role\tbytes\tsha256\tpath\n'
    for entry in "model:$model_json" "policy:$policy" \
                 "blind_gff:$blind_gff" "source_gff:$source_gff"; do
        role=${entry%%:*}; path=${entry#*:}
        printf '%s\t%s\t%s\t%s\n' "$role" "$(stat -Lc %s "$path")" \
            "$(sha256sum "$path" | awk '{print $1}')" "$path"
    done
} > "$working_root/input_manifest.tsv"

cd "$code_root"
for mode in blind complete_control; do
    "$python_bin" -m ploidypatch.cli evidence build-copy-features \
        --consensus-decisions "$method_root/consensus/union/$mode/decisions.tsv" \
        --method-decisions "miniprot=$method_root/methods/miniprot/$mode/decisions.tsv" \
        --method-decisions "gemoma=$method_root/methods/gemoma/$mode/decisions.tsv" \
        --method-decisions "lifton=$method_root/methods/lifton/$mode/decisions.tsv" \
        --wgd-selection "$self_wgd_root/$mode/selected/selection.tsv" \
        --output-tsv "$working_root/$mode/features.tsv" \
        > "$working_root/$mode/features.stdout.json" \
        2> "$working_root/$mode/features.stderr.log"
    "$python_bin" -m ploidypatch.cli evidence score-copy-candidates \
        --features "$working_root/$mode/features.tsv" --model-json "$model_json" \
        --output-tsv "$working_root/$mode/scores.tsv" \
        > "$working_root/$mode/scores.stdout.json" \
        2> "$working_root/$mode/scores.stderr.log"
    if [[ $mode == blind ]]; then base=$blind_gff; else base=$source_gff; fi
    "$python_bin" -m ploidypatch.cli evidence select-scored-copy-candidates \
        --base-gff "$base" \
        --candidate-gff "$method_root/consensus/union/$mode/candidate.gff3" \
        --scores "$working_root/$mode/scores.tsv" --model-json "$model_json" \
        --policy review --output-gff "$working_root/$mode/candidate.gff3" \
        --selection-tsv "$working_root/$mode/selection.tsv" \
        > "$working_root/$mode/selection.stdout.json" \
        2> "$working_root/$mode/selection.stderr.log"
    for mask in wgd_context method_quality; do
        mask_root=$working_root/counterfactual/$mask/$mode
        "$python_bin" -m ploidypatch.cli evidence score-copy-candidates \
            --features "$working_root/$mode/features.tsv" --model-json "$model_json" \
            --mask-feature-group "$mask" --output-tsv "$mask_root/scores.tsv" \
            > "$mask_root/scores.stdout.json" \
            2> "$mask_root/scores.stderr.log"
        "$python_bin" -m ploidypatch.cli evidence select-scored-copy-candidates \
            --base-gff "$base" \
            --candidate-gff "$method_root/consensus/union/$mode/candidate.gff3" \
            --scores "$mask_root/scores.tsv" --model-json "$model_json" \
            --policy review --output-gff "$mask_root/candidate.gff3" \
            --selection-tsv "$mask_root/selection.tsv" \
            > "$mask_root/selection.stdout.json" \
            2> "$mask_root/selection.stderr.log"
    done
done

# Candidate generation is cryptographically frozen before evaluator-only truth access.
{
    printf 'role\tbytes\tsha256\tpath\n'
    find "$working_root/blind" "$working_root/complete_control" \
        "$working_root/counterfactual" -type f \
        \( -name 'features.tsv' -o -name 'scores.tsv' -o -name 'candidate.gff3' \
           -o -name 'selection.tsv' \) -print0 | sort -z \
        | while IFS= read -r -d '' path; do
            role=${path#"$working_root/"}; role=${role//\//_}
            printf '%s\t%s\t%s\t%s\n' "$role" "$(stat -Lc %s "$path")" \
                "$(sha256sum "$path" | awk '{print $1}')" "$path"
        done
} > "$working_root/candidate_freeze.tsv"
date --iso-8601=seconds > "$working_root/candidate_frozen_at.txt"

if [[ ! -s $truth ]]; then
    echo "missing evaluator-only cotton truth: $truth" >&2
    exit 1
fi
{
    printf 'role\tbytes\tsha256\tpath\n'
    printf 'truth\t%s\t%s\t%s\n' "$(stat -Lc %s "$truth")" \
        "$(sha256sum "$truth" | awk '{print $1}')" "$truth"
} > "$working_root/evaluator_input_manifest.tsv"
"$python_bin" -m ploidypatch.cli benchmark label-copy-features \
    --features "$working_root/blind/features.tsv" --truth "$truth" \
    --output-tsv "$working_root/evaluator/labeled_features.tsv" \
    > "$working_root/evaluator/labels.stdout.json" \
    2> "$working_root/evaluator/labels.stderr.log"
"$python_bin" -m ploidypatch.cli benchmark score-copy-ranking \
    --scores "$working_root/blind/scores.tsv" \
    --labeled-features "$working_root/evaluator/labeled_features.tsv" \
    --output-json "$working_root/evaluator/ranking.json" \
    > "$working_root/evaluator/ranking.stdout.json" \
    2> "$working_root/evaluator/ranking.stderr.log"
for mask in wgd_context method_quality; do
    "$python_bin" -m ploidypatch.cli benchmark score-copy-ranking \
        --scores "$working_root/counterfactual/$mask/blind/scores.tsv" \
        --labeled-features "$working_root/evaluator/labeled_features.tsv" \
        --output-json "$working_root/evaluator/ranking.mask_$mask.json" \
        > "$working_root/evaluator/ranking.mask_$mask.stdout.json" \
        2> "$working_root/evaluator/ranking.mask_$mask.stderr.log"
done
/usr/bin/time -v -o "$working_root/evaluator/paired_score.resource.time.txt" \
    "$python_bin" -m ploidypatch.cli benchmark score \
        --source-gff "$source_gff" --perturbed-gff "$blind_gff" \
        --candidate-gff "$working_root/blind/candidate.gff3" \
        --control-candidate-gff "$working_root/complete_control/candidate.gff3" \
        --truth "$truth" --include-event-details \
        > "$working_root/evaluator/paired_score.json" \
        2> "$working_root/evaluator/paired_score.stderr.log"
if ! grep -q '"grade": "pass"' "$working_root/evaluator/paired_score.json"; then
    echo "cotton transfer paired score gate failed" >&2; exit 1
fi
for mask in wgd_context method_quality; do
    /usr/bin/time -v -o "$working_root/evaluator/paired_score.mask_$mask.resource.time.txt" \
        "$python_bin" -m ploidypatch.cli benchmark score \
            --source-gff "$source_gff" --perturbed-gff "$blind_gff" \
            --candidate-gff "$working_root/counterfactual/$mask/blind/candidate.gff3" \
            --control-candidate-gff "$working_root/counterfactual/$mask/complete_control/candidate.gff3" \
            --truth "$truth" --include-event-details \
            > "$working_root/evaluator/paired_score.mask_$mask.json" \
            2> "$working_root/evaluator/paired_score.mask_$mask.stderr.log"
    if ! grep -q '"grade": "pass"' "$working_root/evaluator/paired_score.mask_$mask.json"; then
        echo "cotton transfer $mask counterfactual score gate failed" >&2; exit 1
    fi
done
{
    printf 'role\tbytes\tsha256\tpath\n'
    find "$working_root/blind" "$working_root/complete_control" \
        "$working_root/counterfactual" "$working_root/evaluator" \
        -type f \( -name '*.tsv' -o -name '*.json' -o -name '*.gff3' \) \
        -print0 | sort -z | while IFS= read -r -d '' path; do
            role=${path#"$working_root/"}; role=${role//\//_}
            printf '%s\t%s\t%s\t%s\n' "$role" "$(stat -Lc %s "$path")" \
                "$(sha256sum "$path" | awk '{print $1}')" "$path"
        done
} > "$working_root/output_freeze.tsv"
(
    cd "$working_root"
    find . -type f \( -name '*.gff3' -o -name '*.json' -o -name '*.tsv' \) \
        -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS
    sha256sum -c SHA256SUMS >/dev/null
)
du -sb "$working_root" > "$working_root/disk_bytes.txt"
mv "$working_root" "$result_root"
printf 'cotton zero-retuning Glycine model transfer frozen: %s\n' "$result_root"
