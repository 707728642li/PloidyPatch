#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then echo "usage: $0 PROJECT_ROOT" >&2; exit 2; fi
project_root=$(realpath "$1")
code_root=$project_root/code
python_bin=$project_root/envs/ploidypatch-dev/bin/python
prepared=$project_root/data/derived/candidate_inputs/cotton_syngap_v0.1
upstream_root=$project_root/results/baselines/cotton_holdout_v0.1/syngap
benchmark=$project_root/benchmark/structure/copy_collapse_v0.1/ghi_ad/annotation_copy_collapse_seed20260817
blind_gff=$benchmark/blind/perturbed.gff3
source_gff=$project_root/data/derived/holdout_inputs/cotton_v0.1/hirsutum/primary_chromosomes.gff3
truth=$benchmark/evaluator/truth/hidden_truth.json
result_root=$project_root/results/copy_collapse/holdout/cotton_syngap_v0.1
working_root=${result_root}.working

resolve_full_gff() {
    local run_root=$1 expected_input=$2
    local manifest=$run_root/output_manifest.tsv contract=$run_root/run_contract.tsv
    [[ -s $manifest && -s $contract ]] || return 1
    local configured_input
    configured_input=$(awk -F '\t' '$1 == "sp1_gff" { print $2 }' "$contract")
    [[ -s $configured_input ]] || return 1
    [[ $(sha256sum "$configured_input" | awk '{print $1}') == \
       $(sha256sum "$expected_input" | awk '{print $1}') ]] || return 1
    local path expected_sha
    path=$(awk -F '\t' '$1 == "Ghi.SynGAP.gff3" { print $4 }' "$manifest")
    expected_sha=$(awk -F '\t' '$1 == "Ghi.SynGAP.gff3" { print $3 }' "$manifest")
    [[ -s $path && $(sha256sum "$path" | awk '{print $1}') == "$expected_sha" ]] || return 1
    printf '%s\n' "$path"
}

declare -A target_input=(
    [blind]="$prepared/target/blind/annotation.compat.gff3"
    [complete_control]="$prepared/target/complete_control/annotation.compat.gff3"
)
declare -A base_gff=([blind]="$blind_gff" [complete_control]="$source_gff")
declare -A upstream
for mode in blind complete_control; do
    for ref in gar_a gra_d; do
        upstream[$ref.$mode]=$(resolve_full_gff \
            "$upstream_root/$ref/$mode" "${target_input[$mode]}") || {
            echo "unvalidated cotton SynGAP arm: $ref $mode" >&2; exit 1;
        }
    done
done
for required in "$python_bin" "$blind_gff" "$source_gff" "${upstream[@]}"; do
    [[ -s $required ]] || { echo "missing cotton SynGAP evaluation input: $required" >&2; exit 1; }
done
if [[ -e $result_root || -e $working_root ]]; then
    echo "refusing to overwrite cotton SynGAP evaluation" >&2; exit 1
fi
mkdir -p "$working_root/methods/gar_a/blind" \
    "$working_root/methods/gar_a/complete_control" \
    "$working_root/methods/gra_d/blind" \
    "$working_root/methods/gra_d/complete_control" \
    "$working_root/consensus/blind" "$working_root/consensus/complete_control" \
    "$working_root/evaluator"
{
    printf 'field\tvalue\n'
    printf 'code_commit\t%s\n' "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'method\tSynGAP-1.2.5-genblastg\n'
    printf 'references\tGossypium_arboreum,Gossypium_raimondii\n'
    printf 'within_method_reference_policy\texact_union_one_method_family\n'
    printf 'split\texternal_zero_retuning_holdout\n'
    printf 'candidate_truth_access\tfalse\nevaluator_truth_access\ttrue\n'
    printf 'primary_policy_membership\tfalse_external_comparator_only\n'
    printf 'max_existing_cds_overlap\t0.2\nmax_redundancy_overlap\t0.5\n'
} > "$working_root/run_contract.tsv"

cd "$code_root"
for mode in blind complete_control; do
    for ref in gar_a gra_d; do
        "$python_bin" -m ploidypatch.cli baseline adapt-gff \
            --perturbed-gff "${base_gff[$mode]}" \
            --candidate-gff "${upstream[$ref.$mode]}" --source "syngap_$ref" \
            --max-existing-cds-overlap 0.2 --max-redundancy-overlap 0.5 \
            --infer-missing-cds-phase \
            --output-gff "$working_root/methods/$ref/$mode/candidate.gff3" \
            --decisions-tsv "$working_root/methods/$ref/$mode/decisions.tsv" \
            > "$working_root/methods/$ref/$mode/stdout.json" \
            2> "$working_root/methods/$ref/$mode/stderr.log"
    done
    "$python_bin" -m ploidypatch.cli baseline select-method-consensus \
        --base-gff "${base_gff[$mode]}" \
        --candidate "gar_a=$working_root/methods/gar_a/$mode/candidate.gff3" \
        --candidate "gra_d=$working_root/methods/gra_d/$mode/candidate.gff3" \
        --min-method-support 1 --max-redundancy-overlap 0.5 \
        --output-gff "$working_root/consensus/$mode/candidate.gff3" \
        --decisions-tsv "$working_root/consensus/$mode/decisions.tsv" \
        > "$working_root/consensus/$mode/stdout.json" \
        2> "$working_root/consensus/$mode/stderr.log"
done
{
    printf 'role\tbytes\tsha256\tpath\n'
    find "$working_root/methods" "$working_root/consensus" -type f \
        \( -name candidate.gff3 -o -name decisions.tsv \) -print0 | sort -z \
        | while IFS= read -r -d '' path; do
            role=${path#"$working_root/"}; role=${role//\//_}
            printf '%s\t%s\t%s\t%s\n' "$role" "$(stat -Lc %s "$path")" \
                "$(sha256sum "$path" | awk '{print $1}')" "$path"
        done
} > "$working_root/candidate_freeze.tsv"
date --iso-8601=seconds > "$working_root/candidate_frozen_at.txt"

[[ -s $truth ]] || { echo "missing evaluator-only cotton truth" >&2; exit 1; }
printf 'role\tbytes\tsha256\tpath\nhidden_truth\t%s\t%s\t%s\n' \
    "$(stat -Lc %s "$truth")" "$(sha256sum "$truth" | awk '{print $1}')" "$truth" \
    > "$working_root/evaluator_input_manifest.tsv"
/usr/bin/time -v -o "$working_root/evaluator/score.resource.time.txt" \
    "$python_bin" -m ploidypatch.cli benchmark score \
        --source-gff "$source_gff" --perturbed-gff "$blind_gff" \
        --candidate-gff "$working_root/consensus/blind/candidate.gff3" \
        --control-candidate-gff "$working_root/consensus/complete_control/candidate.gff3" \
        --truth "$truth" --include-event-details \
        > "$working_root/evaluator/score.json" \
        2> "$working_root/evaluator/score.stderr.log"
grep -q '"grade": "pass"' "$working_root/evaluator/score.json" || {
    echo "cotton SynGAP score gate failed" >&2; exit 1;
}
(
    cd "$working_root"
    find . -type f \( -name '*.gff3' -o -name '*.json' -o -name '*.tsv' \) \
        -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS
    sha256sum -c SHA256SUMS >/dev/null
)
du -sb "$working_root" > "$working_root/disk_bytes.txt"
mv "$working_root" "$result_root"
printf 'cotton SynGAP external comparator frozen: %s\n' "$result_root"
