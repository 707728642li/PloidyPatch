#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "usage: $0 PROJECT_ROOT" >&2
    exit 2
fi

project_root=$(realpath "$1")
code_root=$project_root/code
python_bin=$project_root/envs/ploidypatch-dev/bin/python
union_root=$project_root/results/copy_collapse/consensus_union_glycine_dev_v0.1
wgd_root=$project_root/results/copy_collapse/wgd_reanchor_union_glycine_dev_v0.1/self_wgdi
miniprot_root=$project_root/results/copy_collapse/miniprot_glycine_v0.1
gemoma_root=$project_root/results/copy_collapse/gemoma_glycine_v0.1
lifton_root=$project_root/results/copy_collapse/lifton_glycine_v0.1
truth=$project_root/benchmark/structure/copy_collapse_v0.1/gma_v21/annotation_copy_collapse_seed20260815/evaluator/truth/hidden_truth.json
result_root=$project_root/results/copy_collapse/model_development/glycine_feature_matrix_v0.1
working_root=${result_root}.working

for required in "$python_bin" "$truth"; do
    if [[ ! -s $required ]]; then
        echo "missing or empty copy-feature prerequisite: $required" >&2
        exit 1
    fi
done
for mode in blind complete_control; do
    for required in \
        "$union_root/$mode/decisions.tsv" \
        "$wgd_root/$mode/selected/selection.tsv" \
        "$miniprot_root/$mode/decisions.tsv" \
        "$gemoma_root/$mode/decisions.tsv" \
        "$lifton_root/$mode/decisions.tsv"; do
        if [[ ! -s $required ]]; then
            echo "missing or empty copy-feature input: $required" >&2
            exit 1
        fi
    done
done
if [[ -e $result_root || -e $working_root ]]; then
    echo "refusing to overwrite copy-feature development result: $result_root" >&2
    exit 1
fi
mkdir -p "$working_root/blind" "$working_root/complete_control" \
    "$working_root/evaluator"

{
    printf 'field\tvalue\n'
    printf 'code_commit\t%s\n' "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'split\tpost_holdout_exploratory_development\n'
    printf 'formal_holdout_claim_allowed\tfalse\n'
    printf 'candidate_feature_truth_access\tfalse\n'
    printf 'label_stage_truth_access\tevaluator_only\n'
    printf 'candidate_universe\texact_method_family_union_after_redundancy_filter\n'
    printf 'wgd_context\tmode_specific_blind_recomputed_self_WGD\n'
    printf 'copy_feature_module_sha256\t%s\n' \
        "$(sha256sum "$code_root/src/ploidypatch/copy_features.py" | awk '{print $1}')"
} > "$working_root/run_contract.tsv"

{
    printf 'mode\trole\tbytes\tsha256\tpath\n'
    for mode in blind complete_control; do
        for entry in \
            "consensus_decisions:$union_root/$mode/decisions.tsv" \
            "wgd_selection:$wgd_root/$mode/selected/selection.tsv" \
            "miniprot_decisions:$miniprot_root/$mode/decisions.tsv" \
            "gemoma_decisions:$gemoma_root/$mode/decisions.tsv" \
            "lifton_decisions:$lifton_root/$mode/decisions.tsv"; do
            role=${entry%%:*}
            path=${entry#*:}
            printf '%s\t%s\t%s\t%s\t%s\n' "$mode" "$role" \
                "$(stat -Lc %s "$path")" "$(sha256sum "$path" | awk '{print $1}')" \
                "$path"
        done
    done
} > "$working_root/candidate_input_manifest.tsv"

cd "$code_root"
for mode in blind complete_control; do
    /usr/bin/time -v -o "$working_root/$mode/features.resource.time.txt" \
        "$python_bin" -m ploidypatch.cli evidence build-copy-features \
            --consensus-decisions "$union_root/$mode/decisions.tsv" \
            --method-decisions "miniprot=$miniprot_root/$mode/decisions.tsv" \
            --method-decisions "gemoma=$gemoma_root/$mode/decisions.tsv" \
            --method-decisions "lifton=$lifton_root/$mode/decisions.tsv" \
            --wgd-selection "$wgd_root/$mode/selected/selection.tsv" \
            --output-tsv "$working_root/$mode/features.tsv" \
            > "$working_root/$mode/features.stdout.json" \
            2> "$working_root/$mode/features.stderr.log"
done

{
    printf 'mode\tbytes\tsha256\tpath\n'
    for mode in blind complete_control; do
        path=$working_root/$mode/features.tsv
        printf '%s\t%s\t%s\t%s\n' "$mode" "$(stat -Lc %s "$path")" \
            "$(sha256sum "$path" | awk '{print $1}')" "$path"
    done
} > "$working_root/candidate_feature_freeze.tsv"

{
    printf 'role\tbytes\tsha256\tpath\n'
    printf 'hidden_truth\t%s\t%s\t%s\n' "$(stat -Lc %s "$truth")" \
        "$(sha256sum "$truth" | awk '{print $1}')" "$truth"
} > "$working_root/evaluator/input_manifest.tsv"
/usr/bin/time -v -o "$working_root/evaluator/labels.resource.time.txt" \
    "$python_bin" -m ploidypatch.cli benchmark label-copy-features \
        --features "$working_root/blind/features.tsv" \
        --truth "$truth" \
        --output-tsv "$working_root/evaluator/labeled_features.tsv" \
        > "$working_root/evaluator/labels.stdout.json" \
        2> "$working_root/evaluator/labels.stderr.log"

if ! grep -q '"accepted_candidates": 16638' \
        "$working_root/blind/features.stdout.json"; then
    echo "blind copy-feature candidate-count sentinel failed" >&2
    exit 1
fi
if ! grep -q '"positive_exact_cds": 619' \
        "$working_root/evaluator/labels.stdout.json"; then
    echo "copy-feature exact-CDS label sentinel failed" >&2
    exit 1
fi
(
    cd "$working_root"
    find . -type f \( -name '*.json' -o -name '*.tsv' \) \
        -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS
)
du -sb "$working_root" > "$working_root/disk_bytes.txt"
mv "$working_root" "$result_root"
printf 'Glycine copy-feature development matrix frozen: %s\n' "$result_root"
