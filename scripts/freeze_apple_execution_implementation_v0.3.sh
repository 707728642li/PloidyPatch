#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then echo "usage: $0 PROJECT_ROOT" >&2; exit 2; fi
project_root=$(realpath "$1")
code_root=$project_root/code
protocol_root=$project_root/results/protocol_freezes/apple_external_v0.3
result_root=$project_root/results/protocol_freezes/apple_external_v0.3_execution
working_root=${result_root}.working

files=(
    src/ploidypatch/pair_consensus.py
    scripts/freeze_apple_execution_implementation_v0.3.sh
    scripts/build_wgdi_source_alias_gff.py
    scripts/run_apple_evaluator_wgdi_v0.3.sh
    scripts/infer_apple_external_pairs_v0.3.sh
    scripts/run_apple_copy_collapse_benchmark_v0.3.sh
    scripts/build_apple_method_trio_candidate_pools_v0.3.sh
    scripts/run_apple_blind_union_self_wgd_v0.3.sh
    scripts/score_apple_candidates_blind_v0.3.sh
    scripts/evaluate_apple_external_v0.3.py
    scripts/run_apple_external_reveal_v0.3.sh
    tests/test_pair_consensus.py
    tests/test_build_wgdi_source_alias_gff.py
    tests/test_evaluate_apple_external_v0_3.py
    docs/APPLE_EXTERNAL_VALIDATION_EXECUTION_AMENDMENT_v0.3.1.md
    docs/APPLE_EXTERNAL_VALIDATION_EXECUTION_AMENDMENT_v0.3.2.md
)
for required in "$protocol_root/SHA256SUMS" "${files[@]/#/$code_root/}"; do
    [[ -s $required ]] || { echo "missing execution-freeze input: $required" >&2; exit 1; }
done
(cd "$protocol_root" && sha256sum -c SHA256SUMS >/dev/null)
truth_root=$project_root/results/evaluator/apple_v0.3/truth_pairs
benchmark_root=$project_root/benchmark/structure/copy_collapse_v0.3/mdx_gddh13/annotation_copy_collapse_seed20260831
method_root=$project_root/results/copy_collapse/external/apple_v0.3_method_trio
blind_root=$project_root/results/copy_collapse/external/apple_v0.3_blind_self_wgd
ranking_root=$project_root/results/copy_collapse/external/apple_v0.3_blind_rankings
reveal_root=$project_root/results/copy_collapse/external/apple_v0.3_reveal
preexisting=("$truth_root" "$benchmark_root" "$method_root")
preexisting_count=0
for root in "${preexisting[@]}"; do
    [[ -e $root ]] && preexisting_count=$((preexisting_count + 1))
done
if [[ $preexisting_count -eq 0 ]]; then
    amendment_stage=false
elif [[ $preexisting_count -eq ${#preexisting[@]} ]]; then
    amendment_stage=true
    for root in "${preexisting[@]}"; do
        [[ -s $root/SHA256SUMS ]] || {
            echo "pre-amendment artifact lacks SHA256SUMS: $root" >&2; exit 1;
        }
        (cd "$root" && sha256sum -c SHA256SUMS >/dev/null)
    done
else
    echo "partial pre-amendment apple artifact set" >&2
    exit 1
fi
for forbidden in "$blind_root" "$ranking_root" "$reveal_root"; do
    [[ ! -e $forbidden ]] || {
        echo "execution freeze is too late; artifact already exists: $forbidden" >&2
        exit 1
    }
done
[[ ! -e $result_root && ! -e $working_root ]] || {
    echo "refusing to overwrite apple execution freeze" >&2; exit 1;
}
mkdir -p "$working_root"
{
    printf 'field\tvalue\n'
    printf 'code_commit\t%s\n' "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'created_before_pair_enumeration\t%s\n' \
        "$([[ $amendment_stage == false ]] && echo true || echo false)"
    printf 'created_before_hidden_event_sampling\t%s\n' \
        "$([[ $amendment_stage == false ]] && echo true || echo false)"
    printf 'created_before_blind_topology_features\ttrue\n'
    printf 'created_before_candidate_labels\ttrue\n'
    printf 'created_before_candidate_scores\ttrue\n'
    printf 'pre_label_identifier_amendment\t%s\n' "$amendment_stage"
    printf 'amendment_trigger\tidentifier_mapping_failure_before_blind_pair_output\n'
    printf 'protocol_freeze_sha256sums_sha256\t%s\n' \
        "$(sha256sum "$protocol_root/SHA256SUMS" | awk '{print $1}')"
} > "$working_root/run_contract.tsv"
{
    printf 'field\tvalue\n'
    printf 'execution_policy_id\tploidypatch_apple_external_execution_v0.3.2\n'
    printf 'execution_amendment\tunique_exact_gene_attribute_alias_only\n'
    printf 'H1_bootstrap_seed\t20260902\n'
    printf 'H2_bootstrap_seed\t20260901\n'
    printf 'bootstrap_replicates\t20000\n'
    printf 'H2_fixed_sequence\tonly_confirmatory_if_H1_passes\n'
    printf 'failure_retention\ttrue\n'
    printf 'automatic_copy_addition_approval\tfalse\n'
} > "$working_root/policy_supplement.tsv"
{
    printf 'artifact\tbytes\tsha256sums_sha256\tpath\n'
    if [[ $amendment_stage == true ]]; then
        for root in "${preexisting[@]}"; do
            printf '%s\t%s\t%s\t%s\n' "$(basename "$root")" \
                "$(du -sb "$root" | awk '{print $1}')" \
                "$(sha256sum "$root/SHA256SUMS" | awk '{print $1}')" "$root"
        done
    fi
} > "$working_root/frozen_preexisting_artifacts.tsv"
{
    printf 'path\tbytes\tsha256\n'
    for relative in "${files[@]}"; do
        path=$code_root/$relative
        printf '%s\t%s\t%s\n' "$relative" "$(stat -Lc %s "$path")" \
            "$(sha256sum "$path" | awk '{print $1}')"
    done
} > "$working_root/implementation_manifest.tsv"
(
    cd "$working_root"
    sha256sum frozen_preexisting_artifacts.tsv implementation_manifest.tsv \
        policy_supplement.tsv run_contract.tsv > SHA256SUMS
    sha256sum -c SHA256SUMS >/dev/null
)
mv "$working_root" "$result_root"
printf 'apple execution implementation frozen: %s\n' "$result_root"
