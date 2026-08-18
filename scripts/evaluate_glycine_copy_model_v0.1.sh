#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "usage: $0 PROJECT_ROOT" >&2
    exit 2
fi

project_root=$(realpath "$1")
code_root=$project_root/code
python_bin=$project_root/envs/ploidypatch-dev/bin/python
benchmark_root=$project_root/benchmark/structure/copy_collapse_v0.1/gma_v21/annotation_copy_collapse_seed20260815
blind_gff=$benchmark_root/blind/perturbed.gff3
source_gff=$project_root/data/derived/syngap_glycine_v0.1/gma_complete/primary_chromosomes.gff3
truth=$benchmark_root/evaluator/truth/hidden_truth.json
feature_root=$project_root/results/copy_collapse/model_development/glycine_feature_matrix_v0.1
union_root=$project_root/results/copy_collapse/consensus_union_glycine_dev_v0.1
model_root=$project_root/results/copy_collapse/model_development/glycine_copy_model_v0.1
model_json=$model_root/model.json
result_root=$project_root/results/copy_collapse/model_development/glycine_copy_model_evaluation_v0.1
working_root=${result_root}.working

for required in \
    "$python_bin" "$blind_gff" "$source_gff" "$truth" "$model_json" \
    "$model_root/SHA256SUMS" \
    "$feature_root/blind/features.tsv" \
    "$feature_root/blind/features.tsv.manifest.json" \
    "$feature_root/complete_control/features.tsv" \
    "$feature_root/complete_control/features.tsv.manifest.json" \
    "$union_root/blind/candidate.gff3" \
    "$union_root/complete_control/candidate.gff3"; do
    if [[ ! -s $required ]]; then
        echo "missing or empty Glycine copy-model evaluation input: $required" >&2
        exit 1
    fi
done
(cd "$model_root" && sha256sum -c SHA256SUMS >/dev/null)
if [[ -e $result_root || -e $working_root ]]; then
    echo "refusing to overwrite Glycine copy-model evaluation" >&2
    exit 1
fi
mkdir -p "$working_root/blind" "$working_root/complete_control" "$working_root/evaluator"

{
    printf 'field\tvalue\n'
    printf 'code_commit\t%s\n' "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'split\tpost_holdout_model_development\n'
    printf 'formal_holdout_claim_allowed\tfalse\n'
    printf 'reason\tsoybean_truth_was_used_for_grouped_model_development\n'
    printf 'primary_unbiased_development_estimate\tnested_grouped_oof_in_training_report\n'
    printf 'paired_evaluation_role\tfull_fit_engineering_comparison_only\n'
    printf 'candidate_truth_access\tfalse\n'
    printf 'evaluator_truth_access\ttrue\n'
    printf 'paired_complete_annotation_control\ttrue\n'
    printf 'selection_policy\treview\n'
    printf 'automatic_approval\tfalse\n'
    printf 'copy_model_module_sha256\t%s\n' \
        "$(sha256sum "$code_root/src/ploidypatch/copy_model.py" | awk '{print $1}')"
    printf 'score_module_sha256\t%s\n' \
        "$(sha256sum "$code_root/src/ploidypatch/score.py" | awk '{print $1}')"
} > "$working_root/run_contract.tsv"
{
    printf 'role\tbytes\tsha256\tpath\n'
    for entry in \
        "model:$model_json" \
        "training_report:$model_root/training_report.json" \
        "blind_features:$feature_root/blind/features.tsv" \
        "control_features:$feature_root/complete_control/features.tsv" \
        "blind_union_candidate:$union_root/blind/candidate.gff3" \
        "control_union_candidate:$union_root/complete_control/candidate.gff3" \
        "blind_gff:$blind_gff" \
        "source_gff:$source_gff"; do
        role=${entry%%:*}
        path=${entry#*:}
        printf '%s\t%s\t%s\t%s\n' "$role" "$(stat -Lc %s "$path")" \
            "$(sha256sum "$path" | awk '{print $1}')" "$path"
    done
} > "$working_root/candidate_input_manifest.tsv"

cd "$code_root"
for mode in blind complete_control; do
    if [[ $mode == blind ]]; then
        base=$blind_gff
    else
        base=$source_gff
    fi
    /usr/bin/time -v -o "$working_root/$mode/score_features.resource.time.txt" \
        "$python_bin" -m ploidypatch.cli evidence score-copy-candidates \
            --features "$feature_root/$mode/features.tsv" \
            --model-json "$model_json" \
            --output-tsv "$working_root/$mode/scores.tsv" \
            > "$working_root/$mode/score_features.stdout.json" \
            2> "$working_root/$mode/score_features.stderr.log"
    /usr/bin/time -v -o "$working_root/$mode/select.resource.time.txt" \
        "$python_bin" -m ploidypatch.cli evidence select-scored-copy-candidates \
            --base-gff "$base" \
            --candidate-gff "$union_root/$mode/candidate.gff3" \
            --scores "$working_root/$mode/scores.tsv" \
            --model-json "$model_json" \
            --policy review \
            --output-gff "$working_root/$mode/candidate.gff3" \
            --selection-tsv "$working_root/$mode/selection.tsv" \
            > "$working_root/$mode/select.stdout.json" \
            2> "$working_root/$mode/select.stderr.log"
done

{
    printf 'role\tbytes\tsha256\tpath\n'
    for mode in blind complete_control; do
        for name in scores.tsv selection.tsv candidate.gff3; do
            path=$working_root/$mode/$name
            printf '%s_%s\t%s\t%s\t%s\n' "$mode" "$name" \
                "$(stat -Lc %s "$path")" \
                "$(sha256sum "$path" | awk '{print $1}')" "$path"
        done
    done
} > "$working_root/candidate_freeze.tsv"
{
    printf 'role\tbytes\tsha256\tpath\n'
    printf 'hidden_truth\t%s\t%s\t%s\n' "$(stat -Lc %s "$truth")" \
        "$(sha256sum "$truth" | awk '{print $1}')" "$truth"
} > "$working_root/evaluator/input_manifest.tsv"

/usr/bin/time -v -o "$working_root/evaluator/score.resource.time.txt" \
    "$python_bin" -m ploidypatch.cli benchmark score \
        --source-gff "$source_gff" \
        --perturbed-gff "$blind_gff" \
        --candidate-gff "$working_root/blind/candidate.gff3" \
        --control-candidate-gff "$working_root/complete_control/candidate.gff3" \
        --truth "$truth" \
        --include-event-details \
        > "$working_root/evaluator/score.json" \
        2> "$working_root/evaluator/score.stderr.log"
if ! grep -q '"grade": "pass"' "$working_root/evaluator/score.json"; then
    echo "Glycine copy-model evaluator quality gate failed" >&2
    exit 1
fi
(
    cd "$working_root"
    find . -type f \( -name '*.gff3' -o -name '*.json' -o -name '*.tsv' \) \
        -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS
    sha256sum -c SHA256SUMS >/dev/null
)
du -sb "$working_root" > "$working_root/disk_bytes.txt"
mv "$working_root" "$result_root"
printf 'Glycine copy-model engineering evaluation frozen: %s\n' "$result_root"
