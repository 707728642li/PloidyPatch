#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then echo "usage: $0 PROJECT_ROOT" >&2; exit 2; fi
project_root=$(realpath "$1")
code_root=$project_root/code
python_bin=$project_root/envs/ploidypatch-dev/bin/python
source_root=$project_root/data/derived/holdout_inputs/maize_v2
source_gff=$source_root/zea_mays/primary_chromosomes.gff3
pair_root=$project_root/results/evaluator/maize_v2/outgroup_duplicated_pairs
pair_tsv=$pair_root/zma.outgroup_duplicated_pairs.tsv
policy=$code_root/config/maize_v2_zero_retuning_policy.tsv
model=$project_root/results/copy_collapse/model_development/homeolog_ranker_v0.2/model.json
seed=20260829
maximum_count=800
result_root=$project_root/benchmark/structure/copy_collapse_v0.2/zma_maize1/annotation_copy_collapse_seed${seed}
working_root=${result_root}.working

for required in "$python_bin" "$source_root/SHA256SUMS" "$source_gff" \
                "$pair_root/SHA256SUMS" "$pair_tsv" "$policy" "$model"; do
    [[ -s $required ]] || { echo "missing maize holdout prerequisite: $required" >&2; exit 1; }
done
(cd "$source_root" && sha256sum -c SHA256SUMS >/dev/null)
(cd "$pair_root" && sha256sum -c SHA256SUMS >/dev/null)
expected_model_sha=$(awk -F '\t' '$1 == "model_sha256" {print $2}' "$policy")
observed_model_sha=$(sha256sum "$model" | awk '{print $1}')
[[ -n $expected_model_sha && $observed_model_sha == "$expected_model_sha" ]] || {
    echo "frozen model hash disagrees with maize policy" >&2; exit 1;
}
if [[ -e $result_root || -e $working_root ]]; then
    echo "refusing to overwrite maize v0.2 holdout" >&2; exit 1
fi
mkdir -p "$working_root/blind" "$working_root/evaluator/pair_selection" \
    "$working_root/evaluator/truth" "$working_root/evaluator/validation"
{
    printf 'field\tvalue\n'
    printf 'code_commit\t%s\n' "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'dataset_id\tzma_maize1\nevent_type\tannotation_copy_collapse\n'
    printf 'requested_count\t%s\nseed\t%s\nsplit\texternal_zero_retuning_holdout_v0.2\n' \
        "$maximum_count" "$seed"
    printf 'pair_access\tevaluator_only\npolicy_frozen_before_truth\ttrue\n'
    printf 'candidate_model_frozen_before_truth\ttrue\nautomatic_approval\tfalse\n'
    printf 'policy_sha256\t%s\n' "$(sha256sum "$policy" | awk '{print $1}')"
    printf 'model_sha256\t%s\n' "$observed_model_sha"
    printf 'pair_tsv_sha256\t%s\n' "$(sha256sum "$pair_tsv" | awk '{print $1}')"
    printf 'pair_sampler_module_sha256\t%s\n' \
        "$(sha256sum "$code_root/src/ploidypatch/copy_pair_sampling.py" | awk '{print $1}')"
    printf 'structure_module_sha256\t%s\n' \
        "$(sha256sum "$code_root/src/ploidypatch/structure_perturb.py" | awk '{print $1}')"
    printf 'score_module_sha256\t%s\n' \
        "$(sha256sum "$code_root/src/ploidypatch/score.py" | awk '{print $1}')"
} > "$working_root/run_contract.tsv"

selected=$working_root/evaluator/pair_selection/selected_pairs.tsv
selection_decisions=$working_root/evaluator/pair_selection/decisions.tsv
cd "$code_root"
"$python_bin" -m ploidypatch.cli benchmark sample-copy-pairs \
    --source-gff "$source_gff" --pairs "$pair_tsv" \
    --count "$maximum_count" --seed "$seed" \
    --output-pairs "$selected" --decisions-tsv "$selection_decisions" \
    > "$working_root/evaluator/pair_selection/stdout.json" \
    2> "$working_root/evaluator/pair_selection/stderr.log"
selected_count=$(( $(wc -l < "$selected") - 1 ))
[[ $selected_count -gt 0 ]] || { echo "maize selected pair set is empty" >&2; exit 1; }
printf 'selected_count\t%s\n' "$selected_count" >> "$working_root/run_contract.tsv"

"$python_bin" -m ploidypatch.cli benchmark perturb \
    --gff "$source_gff" --output-dir "$working_root/blind" \
    --truth-dir "$working_root/evaluator/truth" \
    --event-type annotation_copy_collapse --pair-tsv "$selected" \
    --count "$selected_count" --seed "$seed" \
    > "$working_root/perturb.stdout.json" 2> "$working_root/perturb.stderr.log"
perturbed=$working_root/blind/perturbed.gff3
truth=$working_root/evaluator/truth/hidden_truth.json
restored=$working_root/evaluator/restored.gff3
"$python_bin" scripts/audit_copy_pair_selection_truth.py \
    --selected-pairs "$selected" --truth "$truth" \
    --output-json "$working_root/evaluator/validation/pair_truth_audit.json"
"$python_bin" -m ploidypatch.cli benchmark restore \
    --perturbed-gff "$perturbed" --truth "$truth" --output-gff "$restored" \
    > "$working_root/evaluator/restoration_report.json" \
    2> "$working_root/evaluator/restoration.stderr.log"
cmp -s "$source_gff" "$restored" || {
    echo "maize restoration is not byte-identical" >&2; exit 1;
}
for mode in noop oracle; do
    candidate=$perturbed; [[ $mode == oracle ]] && candidate=$source_gff
    "$python_bin" -m ploidypatch.cli benchmark score \
        --source-gff "$source_gff" --perturbed-gff "$perturbed" \
        --candidate-gff "$candidate" --truth "$truth" --include-event-details \
        > "$working_root/evaluator/validation/score_$mode.json" \
        2> "$working_root/evaluator/validation/score_$mode.stderr.log"
done
grep -q '"complete_cds_chain_recovery": 0' \
    "$working_root/evaluator/validation/score_noop.json"
grep -q "\"complete_cds_chain_recovery\": $selected_count" \
    "$working_root/evaluator/validation/score_oracle.json"
grep -q '"grade": "pass"' "$working_root/evaluator/validation/score_noop.json"
grep -q '"grade": "pass"' "$working_root/evaluator/validation/score_oracle.json"
grep -q '"grade": "pass"' \
    "$working_root/evaluator/validation/pair_truth_audit.json"
(
    cd "$working_root"
    find . -type f \( -name '*.gff3' -o -name '*.json' -o -name '*.tsv' \) \
        -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS
    sha256sum -c SHA256SUMS >/dev/null
)
du -sb "$working_root" > "$working_root/disk_bytes.txt"
mv "$working_root" "$result_root"
printf 'maize v0.2 zero-retuning holdout frozen: %s\n' "$result_root"
