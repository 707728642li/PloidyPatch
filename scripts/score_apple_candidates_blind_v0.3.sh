#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then echo "usage: $0 PROJECT_ROOT" >&2; exit 2; fi
project_root=$(realpath "$1")
code_root=$project_root/code
python_bin=$project_root/envs/ploidypatch-dev/bin/python
benchmark=$project_root/benchmark/structure/copy_collapse_v0.3/mdx_gddh13/annotation_copy_collapse_seed20260831
base_gff=$benchmark/blind/perturbed.gff3
method_root=$project_root/results/copy_collapse/external/apple_v0.3_method_trio
candidate_gff=$method_root/consensus/primary_union/blind/candidate.gff3
pool_decisions=$method_root/consensus/primary_union/blind/decisions.tsv
self_wgd_root=$project_root/results/copy_collapse/external/apple_v0.3_blind_self_wgd
prior_wgd=$self_wgd_root/selected/selection.tsv
protocol_root=$project_root/results/protocol_freezes/apple_external_v0.3
execution_root=$project_root/results/protocol_freezes/apple_external_v0.3_execution
model_v03=$project_root/results/models/support_conditioned_ranker_v0.3/model.json
model_v02=$project_root/results/copy_collapse/model_development/homeolog_ranker_v0.2/model.json
result_root=$project_root/results/copy_collapse/external/apple_v0.3_blind_rankings
working_root=${result_root}.working

for required in "$python_bin" "$benchmark/SHA256SUMS" "$method_root/SHA256SUMS" \
    "$self_wgd_root/SHA256SUMS" "$protocol_root/SHA256SUMS" \
    "$execution_root/SHA256SUMS" "$base_gff" \
    "$candidate_gff" "$pool_decisions" "$prior_wgd" "$model_v03" \
    "$method_root/methods/miniprot/blind/decisions.tsv" \
    "$method_root/methods/gemoma/blind/decisions.tsv" \
    "$method_root/methods/lifton/blind/decisions.tsv"; do
    [[ -s $required ]] || { echo "missing apple blind ranking input: $required" >&2; exit 1; }
done
for root in "$benchmark" "$method_root" "$self_wgd_root" "$protocol_root" \
            "$execution_root"; do
    (cd "$root" && sha256sum -c SHA256SUMS >/dev/null)
done
expected_self=$(awk -F '\t' '$1 == "scripts/score_apple_candidates_blind_v0.3.sh" {print $3}' \
    "$execution_root/implementation_manifest.tsv")
[[ -n $expected_self && $(sha256sum "$code_root/scripts/score_apple_candidates_blind_v0.3.sh" | awk '{print $1}') == "$expected_self" ]] || {
    echo "apple blind scoring script differs from execution freeze" >&2; exit 1;
}
expected_model=$(awk -F '\t' '$1 == "model_sha256" {print $2}' "$protocol_root/policy.tsv")
[[ $(sha256sum "$model_v03" | awk '{print $1}') == "$expected_model" ]] || {
    echo "apple v0.3 model hash disagrees with frozen policy" >&2; exit 1;
}
for relative in src/ploidypatch/copy_features.py \
    src/ploidypatch/homeolog_topology.py src/ploidypatch/wgd_candidate_select.py \
    src/ploidypatch/support_ranker.py; do
    expected=$(awk -F '\t' -v path="$relative" '$1 == path {print $2}' "$protocol_root/code_manifest.tsv")
    observed=$(sha256sum "$code_root/$relative" | awk '{print $1}')
    [[ -n $expected && $observed == "$expected" ]] || {
        echo "post-freeze module change detected: $relative" >&2; exit 1;
    }
done
[[ ! -e $result_root && ! -e $working_root ]] || {
    echo "refusing to overwrite apple blind rankings" >&2; exit 1;
}
mkdir -p "$working_root"/{features,scores,freeze}
{
    printf 'field\tvalue\ncode_commit\t%s\n' \
        "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'split\tuntouched_external_v0.3\ntruth_access\tfalse\n'
    printf 'hidden_pair_access\tfalse\nexternal_label_access\tfalse\n'
    printf 'candidate_policy\tretain_distinct_phased_CDS_chains\n'
    printf 'wgd_source\tblind_candidate_recomputation_only\n'
    printf 'primary_model_sha256\t%s\n' "$expected_model"
    if [[ -s $model_v02 ]]; then
        printf 'secondary_v02_model_sha256\t%s\n' "$(sha256sum "$model_v02" | awk '{print $1}')"
    else
        printf 'secondary_v02_model_sha256\tunavailable\n'
    fi
    printf 'automatic_approval\tfalse\n'
    printf 'protocol_freeze_sha256sums_sha256\t%s\n' \
        "$(sha256sum "$protocol_root/SHA256SUMS" | awk '{print $1}')"
} > "$working_root/run_contract.tsv"

cd "$code_root"
"$python_bin" -m ploidypatch.cli evidence propagate-wgd-conflict-partners \
    --base-gff "$base_gff" --candidate-gff "$candidate_gff" \
    --pool-decisions "$pool_decisions" --prior-wgd-selection "$prior_wgd" \
    --output-selection "$working_root/features/wgd_selection.tsv" \
    > "$working_root/features/wgd_selection.stdout.json" \
    2> "$working_root/features/wgd_selection.stderr.log"
"$python_bin" -m ploidypatch.cli evidence build-copy-features \
    --consensus-decisions "$pool_decisions" \
    --method-decisions "miniprot=$method_root/methods/miniprot/blind/decisions.tsv" \
    --method-decisions "gemoma=$method_root/methods/gemoma/blind/decisions.tsv" \
    --method-decisions "lifton=$method_root/methods/lifton/blind/decisions.tsv" \
    --wgd-selection "$working_root/features/wgd_selection.tsv" \
    --output-tsv "$working_root/features/copy_features.tsv" \
    > "$working_root/features/copy_features.stdout.json" \
    2> "$working_root/features/copy_features.stderr.log"
"$python_bin" -m ploidypatch.cli evidence build-homeolog-topology-features \
    --copy-features "$working_root/features/copy_features.tsv" \
    --wgd-selection "$working_root/features/wgd_selection.tsv" \
    --candidate-gff "$candidate_gff" --base-gff "$base_gff" \
    --output-tsv "$working_root/features/topology_features.tsv" \
    > "$working_root/features/topology_features.stdout.json" \
    2> "$working_root/features/topology_features.stderr.log"
{
    printf 'feature\tbytes\tsha256\n'
    for path in "$working_root/features/wgd_selection.tsv" \
        "$working_root/features/copy_features.tsv" \
        "$working_root/features/topology_features.tsv"; do
        printf '%s\t%s\t%s\n' "$(basename "$path")" "$(stat -Lc %s "$path")" \
            "$(sha256sum "$path" | awk '{print $1}')"
    done
} > "$working_root/freeze/blind_feature_freeze.tsv"

"$python_bin" -m ploidypatch.cli evidence score-support-conditioned-candidates \
    --copy-features "$working_root/features/copy_features.tsv" \
    --topology-features "$working_root/features/topology_features.tsv" \
    --model-json "$model_v03" --output-tsv "$working_root/scores/v03.tsv" \
    > "$working_root/scores/v03.stdout.json" \
    2> "$working_root/scores/v03.stderr.log"
if [[ -s $model_v02 ]]; then
    "$python_bin" -m ploidypatch.cli evidence score-homeolog-copy-candidates \
        --copy-features "$working_root/features/copy_features.tsv" \
        --topology-features "$working_root/features/topology_features.tsv" \
        --model-json "$model_v02" --output-tsv "$working_root/scores/v02.tsv" \
        > "$working_root/scores/v02.stdout.json" \
        2> "$working_root/scores/v02.stderr.log"
fi
for output in features/wgd_selection.tsv features/copy_features.tsv \
    features/topology_features.tsv scores/v03.tsv freeze/blind_feature_freeze.tsv; do
    [[ -s $working_root/$output ]] || { echo "missing apple blind score output: $output" >&2; exit 1; }
done
{
    printf 'score\tbytes\tsha256\n'
    for path in "$working_root"/scores/*.tsv; do
        printf '%s\t%s\t%s\n' "$(basename "$path")" "$(stat -Lc %s "$path")" \
            "$(sha256sum "$path" | awk '{print $1}')"
    done
} > "$working_root/freeze/blind_score_freeze.tsv"
(
    cd "$working_root"
    find . -type f \( -name '*.json' -o -name '*.tsv' \) -print0 \
        | sort -z | xargs -0 sha256sum > SHA256SUMS
    sha256sum -c SHA256SUMS >/dev/null
)
mv "$working_root" "$result_root"
printf 'apple truth-blind v0.3 rankings frozen: %s\n' "$result_root"
