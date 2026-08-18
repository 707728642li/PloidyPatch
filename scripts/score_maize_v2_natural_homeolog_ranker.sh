#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then echo "usage: $0 PROJECT_ROOT" >&2; exit 2; fi
project_root=$(realpath "$1")
code_root=$project_root/code
python_bin=$project_root/envs/ploidypatch-dev/bin/python
bundle=$project_root/data/derived/holdout_inputs/maize_v2/zea_mays
base_gff=$bundle/primary_chromosomes.gff3
method_root=$project_root/results/natural/maize_v2/discovery/method_trio
self_wgd=$project_root/results/natural/maize_v2/discovery/union_self_wgd
selection=$self_wgd/selected/selection.tsv
candidate_gff=$method_root/consensus/union/natural/candidate.gff3
model_root=$project_root/results/copy_collapse/model_development/homeolog_ranker_v0.2
model=$model_root/model.json
policy=$code_root/config/maize_v2_zero_retuning_policy.tsv
protocol=$code_root/docs/MAIZE_NATURAL_VALIDATION_PROTOCOL_v0.1.md
result_root=$project_root/results/natural/maize_v2/discovery/homeolog_ranker
working_root=${result_root}.working

for frozen in "$method_root" "$self_wgd" "$model_root"; do
    [[ -s $frozen/SHA256SUMS ]] || { echo "unfrozen maize natural ranker input: $frozen" >&2; exit 1; }
    (cd "$frozen" && sha256sum -c SHA256SUMS >/dev/null)
done
for required in "$python_bin" "$base_gff" "$selection" "$candidate_gff" \
                "$model" "$policy" "$protocol"; do
    [[ -s $required ]] || { echo "missing maize natural ranker input: $required" >&2; exit 1; }
done
expected_model_sha=$(awk -F '\t' '$1 == "model_sha256" {print $2}' "$policy")
observed_model_sha=$(sha256sum "$model" | awk '{print $1}')
[[ -n $expected_model_sha && $observed_model_sha == "$expected_model_sha" ]] || {
    echo "maize natural ranker model hash disagrees with frozen policy" >&2; exit 1;
}
if [[ -e $result_root || -e $working_root ]]; then
    echo "refusing to overwrite maize natural ranker" >&2; exit 1
fi
mkdir -p "$working_root/natural"
{
    printf 'field\tvalue\ncode_commit\t%s\n' \
        "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'split\tnatural_current_annotation_v0.1\n'
    printf 'truth_access\tfalse\nvalidation_evidence_access\tfalse\n'
    printf 'candidate_freeze_precedes_validation_access\ttrue\n'
    printf 'model_sha256\t%s\n' "$observed_model_sha"
    printf 'primary_estimator\ttopology\ncalibrated_probability\tfalse\n'
    printf 'portable_threshold\tnone\nautomatic_approval\tfalse\n'
    printf 'review_budgets\t25,50,100,200\n'
    printf 'protocol_sha256\t%s\n' "$(sha256sum "$protocol" | awk '{print $1}')"
} > "$working_root/run_contract.tsv"
{
    printf 'role\tbytes\tsha256\tpath\n'
    for entry in "base_gff:$base_gff" "candidate_gff:$candidate_gff" \
                 "wgd_selection:$selection" "model:$model"; do
        role=${entry%%:*}; path=${entry#*:}
        printf '%s\t%s\t%s\t%s\n' "$role" "$(stat -Lc %s "$path")" \
            "$(sha256sum "$path" | awk '{print $1}')" "$path"
    done
} > "$working_root/input_manifest.tsv"

cd "$code_root"
"$python_bin" -m ploidypatch.cli evidence build-copy-features \
    --consensus-decisions "$method_root/consensus/union/natural/decisions.tsv" \
    --method-decisions "miniprot=$method_root/methods/miniprot/natural/decisions.tsv" \
    --method-decisions "gemoma=$method_root/methods/gemoma/natural/decisions.tsv" \
    --method-decisions "lifton=$method_root/methods/lifton/natural/decisions.tsv" \
    --wgd-selection "$selection" \
    --output-tsv "$working_root/natural/copy_features.tsv" \
    > "$working_root/natural/copy_features.stdout.json" \
    2> "$working_root/natural/copy_features.stderr.log"
"$python_bin" -m ploidypatch.cli evidence build-homeolog-topology-features \
    --copy-features "$working_root/natural/copy_features.tsv" \
    --wgd-selection "$selection" --candidate-gff "$candidate_gff" \
    --base-gff "$base_gff" \
    --output-tsv "$working_root/natural/topology_features.tsv" \
    > "$working_root/natural/topology.stdout.json" \
    2> "$working_root/natural/topology.stderr.log"
"$python_bin" -m ploidypatch.cli evidence score-homeolog-copy-candidates \
    --copy-features "$working_root/natural/copy_features.tsv" \
    --topology-features "$working_root/natural/topology_features.tsv" \
    --model-json "$model" --output-tsv "$working_root/natural/scores.tsv" \
    > "$working_root/natural/scores.stdout.json" \
    2> "$working_root/natural/scores.stderr.log"
"$python_bin" -m ploidypatch.cli evidence freeze-homeolog-review-rankings \
    --scores "$working_root/natural/scores.tsv" \
    --review-budget 25 --review-budget 50 --review-budget 100 --review-budget 200 \
    --output-tsv "$working_root/natural/review_rankings.tsv" \
    > "$working_root/natural/review_rankings.stdout.json" \
    2> "$working_root/natural/review_rankings.stderr.log"
for output in copy_features.tsv copy_features.tsv.manifest.json \
              topology_features.tsv topology_features.tsv.manifest.json \
              scores.tsv scores.tsv.manifest.json \
              review_rankings.tsv review_rankings.tsv.manifest.json; do
    [[ -s $working_root/natural/$output ]] || { echo "missing maize natural rank output: $output" >&2; exit 1; }
done
{
    printf 'role\tbytes\tsha256\tpath\n'
    for name in copy_features.tsv topology_features.tsv scores.tsv review_rankings.tsv; do
        path=$working_root/natural/$name
        printf '%s\t%s\t%s\t%s\n' "$name" "$(stat -Lc %s "$path")" \
            "$(sha256sum "$path" | awk '{print $1}')" "$path"
    done
} > "$working_root/candidate_score_freeze.tsv"
(
    cd "$working_root"
    find . -type f \( -name '*.json' -o -name '*.tsv' \) -print0 \
        | sort -z | xargs -0 sha256sum > SHA256SUMS
    sha256sum -c SHA256SUMS >/dev/null
)
du -sb "$working_root" > "$working_root/disk_bytes.txt"
mv "$working_root" "$result_root"
printf 'maize natural review rankings frozen: %s\n' "$result_root"

