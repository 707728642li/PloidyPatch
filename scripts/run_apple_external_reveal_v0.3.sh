#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then echo "usage: $0 PROJECT_ROOT" >&2; exit 2; fi
project_root=$(realpath "$1")
code_root=$project_root/code
dev_python=$project_root/envs/ploidypatch-dev/bin/python
model_python=$project_root/envs/ploidypatch-model/bin/python
source_gff=$project_root/data/derived/external_inputs/apple_v0.3/target_apple/primary_chromosomes.gff3
benchmark=$project_root/benchmark/structure/copy_collapse_v0.3/mdx_gddh13/annotation_copy_collapse_seed20260831
blind_gff=$benchmark/blind/perturbed.gff3
truth=$benchmark/evaluator/truth/hidden_truth.json
evaluability=$benchmark/evaluator/pair_selection/evaluability.json
method_root=$project_root/results/copy_collapse/external/apple_v0.3_method_trio
ranking_root=$project_root/results/copy_collapse/external/apple_v0.3_blind_rankings
protocol_root=$project_root/results/protocol_freezes/apple_external_v0.3
execution_root=$project_root/results/protocol_freezes/apple_external_v0.3_execution
result_root=$project_root/results/copy_collapse/external/apple_v0.3_reveal
working_root=${result_root}.working

for required in "$dev_python" "$model_python" "$source_gff" "$blind_gff" \
    "$truth" "$evaluability" "$benchmark/SHA256SUMS" \
    "$method_root/SHA256SUMS" "$ranking_root/SHA256SUMS" \
    "$protocol_root/SHA256SUMS" "$execution_root/SHA256SUMS" \
    "$ranking_root/scores/v03.tsv" \
    "$ranking_root/scores/v03.tsv.manifest.json" \
    "$ranking_root/freeze/blind_score_freeze.tsv" \
    "$ranking_root/features/copy_features.tsv" \
    "$method_root/consensus/primary_union/blind/decisions.tsv"; do
    [[ -s $required ]] || { echo "missing apple reveal input: $required" >&2; exit 1; }
done
for root in "$benchmark" "$method_root" "$ranking_root" \
            "$protocol_root" "$execution_root"; do
    (cd "$root" && sha256sum -c SHA256SUMS >/dev/null)
done
while IFS=$'\t' read -r relative bytes expected; do
    [[ $relative == path ]] && continue
    path=$code_root/$relative
    [[ -s $path && $(stat -Lc %s "$path") == "$bytes" \
        && $(sha256sum "$path" | awk '{print $1}') == "$expected" ]] || {
        echo "post-execution-freeze implementation change: $relative" >&2; exit 1;
    }
done < "$execution_root/implementation_manifest.tsv"
[[ ! -e $result_root && ! -e $working_root ]] || {
    echo "refusing to overwrite apple v0.3 reveal" >&2; exit 1;
}
mkdir -p "$working_root"/{labels,scores/methods,scores/consensus,evaluation}
{
    printf 'field\tvalue\n'
    printf 'code_commit\t%s\n' "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'access\tevaluator_only_after_blind_score_freeze\n'
    printf 'blind_score_sha256\t%s\n' \
        "$(sha256sum "$ranking_root/scores/v03.tsv" | awk '{print $1}')"
    printf 'blind_score_freeze_sha256\t%s\n' \
        "$(sha256sum "$ranking_root/freeze/blind_score_freeze.tsv" | awk '{print $1}')"
    printf 'truth_opened_after_score_freeze\ttrue\nmodel_refit\tfalse\n'
    printf 'threshold_tuning\tfalse\nautomatic_approval\tfalse\n'
    printf 'failure_retained_without_apple_retuning\ttrue\n'
} > "$working_root/run_contract.tsv"

cd "$code_root"
"$dev_python" -m ploidypatch.cli benchmark label-copy-features \
    --features "$ranking_root/features/copy_features.tsv" --truth "$truth" \
    --output-tsv "$working_root/labels/labeled_features.tsv" \
    > "$working_root/labels/stdout.json" 2> "$working_root/labels/stderr.log"

pids=(); labels=()
for method in miniprot gemoma lifton; do
    (
        "$dev_python" -m ploidypatch.cli benchmark score \
            --source-gff "$source_gff" --perturbed-gff "$blind_gff" \
            --candidate-gff "$method_root/methods/$method/blind/candidate.gff3" \
            --control-candidate-gff "$method_root/methods/$method/complete_control/candidate.gff3" \
            --truth "$truth" --include-event-details \
            > "$working_root/scores/methods/$method.json" \
            2> "$working_root/scores/methods/$method.stderr.log"
    ) & pids+=("$!"); labels+=("method:$method")
done
for pool in primary_union legacy_union support2 support3; do
    (
        "$dev_python" -m ploidypatch.cli benchmark score \
            --source-gff "$source_gff" --perturbed-gff "$blind_gff" \
            --candidate-gff "$method_root/consensus/$pool/blind/candidate.gff3" \
            --control-candidate-gff "$method_root/consensus/$pool/complete_control/candidate.gff3" \
            --truth "$truth" --include-event-details \
            > "$working_root/scores/consensus/$pool.json" \
            2> "$working_root/scores/consensus/$pool.stderr.log"
    ) & pids+=("$!"); labels+=("consensus:$pool")
done
failed=0
for index in "${!pids[@]}"; do
    if ! wait "${pids[$index]}"; then
        echo "failed apple evaluator score: ${labels[$index]}" >&2; failed=1
    fi
done
[[ $failed -eq 0 ]] || exit 1
for score in "$working_root"/scores/{methods,consensus}/*.json; do
    grep -q '"grade": "pass"' "$score" || {
        echo "apple score quality gate failed: $score" >&2; exit 1;
    }
done

secondary_args=(
    --secondary-score "miniprot=$working_root/scores/methods/miniprot.json"
    --secondary-score "gemoma=$working_root/scores/methods/gemoma.json"
    --secondary-score "lifton=$working_root/scores/methods/lifton.json"
    --secondary-score "support2=$working_root/scores/consensus/support2.json"
    --secondary-score "support3=$working_root/scores/consensus/support3.json"
)
v02_args=()
if [[ -s $ranking_root/scores/v02.tsv ]]; then
    v02_args=(--v02-scores "$ranking_root/scores/v02.tsv")
fi
"$model_python" "$code_root/scripts/evaluate_apple_external_v0.3.py" \
    --scores "$ranking_root/scores/v03.tsv" \
    --labels "$working_root/labels/labeled_features.tsv" \
    --pool-decisions "$method_root/consensus/primary_union/blind/decisions.tsv" \
    --primary-pool-score "$working_root/scores/consensus/primary_union.json" \
    --legacy-pool-score "$working_root/scores/consensus/legacy_union.json" \
    --evaluability "$evaluability" --policy "$protocol_root/policy.tsv" \
    --execution-policy "$execution_root/policy_supplement.tsv" \
    --protocol "$protocol_root/protocol.md" \
    "${secondary_args[@]}" "${v02_args[@]}" \
    --output-dir "$working_root/evaluation/final"
(
    cd "$working_root"
    find . -type f \( -name '*.json' -o -name '*.tsv' \) -print0 \
        | sort -z | xargs -0 sha256sum > SHA256SUMS
    sha256sum -c SHA256SUMS >/dev/null
)
mv "$working_root" "$result_root"
printf 'apple v0.3 external result revealed and frozen: %s\n' "$result_root"
