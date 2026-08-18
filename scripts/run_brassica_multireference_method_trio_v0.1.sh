#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "usage: $0 PROJECT_ROOT" >&2
    exit 2
fi

project_root=$(realpath "$1")
code_root=$project_root/code
python_bin=$project_root/envs/ploidypatch-dev/bin/python
benchmark_root=$project_root/benchmark/structure/copy_collapse_v0.1/bna_daae/annotation_copy_collapse_seed20260814
blind_gff=$benchmark_root/blind/perturbed.gff3
truth=$benchmark_root/evaluator/truth/hidden_truth.json
source_gff=$project_root/data/derived/normalized_bundles/v0.1/bna_daae_primary/primary_chromosomes.gff3
upstream_root=$project_root/results/baselines/multireference_brassica_v0.1
miniprot_root=$project_root/results/copy_collapse/miniprot_brassica_v0.1
result_root=$project_root/results/copy_collapse/model_development/brassica_method_trio_v0.1
working_root=${result_root}.working
resume_stage=${PLOIDYPATCH_RESUME_STAGE:-fresh}
case $resume_stage in
    fresh|consensus) ;;
    *) echo "PLOIDYPATCH_RESUME_STAGE must be fresh or consensus" >&2; exit 2 ;;
esac

declare -A raw_gff=(
    [gemoma_brapa]="$upstream_root/gemoma/brapa/upstream/final_annotation.gff"
    [gemoma_bol]="$upstream_root/gemoma/bol/upstream/final_annotation.gff"
    [lifton_brapa]="$upstream_root/lifton/brapa/upstream/lifton.gff3"
    [lifton_bol]="$upstream_root/lifton/bol/upstream/lifton.gff3"
)
for required in "$python_bin" "$blind_gff" "$truth" "$source_gff" \
                "$miniprot_root/SHA256SUMS" \
                "$miniprot_root/blind/candidate.gff3" \
                "$miniprot_root/complete_control/candidate.gff3" \
                "${raw_gff[gemoma_brapa]}" "${raw_gff[gemoma_bol]}" \
                "${raw_gff[lifton_brapa]}" "${raw_gff[lifton_bol]}"; do
    if [[ ! -s $required ]]; then
        echo "missing or empty Brassica method-trio input: $required" >&2
        exit 1
    fi
done
(cd "$miniprot_root" && sha256sum -c SHA256SUMS >/dev/null)
if [[ $resume_stage == fresh ]]; then
if [[ -e $result_root || -e $working_root ]]; then
    echo "refusing to overwrite Brassica method-trio result" >&2
    exit 1
fi
mkdir -p "$working_root/merged" "$working_root/methods/gemoma/blind" \
    "$working_root/methods/gemoma/complete_control" \
    "$working_root/methods/lifton/blind" \
    "$working_root/methods/lifton/complete_control" \
    "$working_root/consensus/union/blind" \
    "$working_root/consensus/union/complete_control" \
    "$working_root/consensus/support2/blind" \
    "$working_root/consensus/support2/complete_control" \
    "$working_root/consensus/support3/blind" \
    "$working_root/consensus/support3/complete_control" \
    "$working_root/evaluator"

{
    printf 'field\tvalue\n'
    printf 'code_commit\t%s\n' "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'split\tpost_holdout_secondary_development\n'
    printf 'formal_holdout_claim_allowed\tfalse\n'
    printf 'candidate_truth_access\tfalse\n'
    printf 'evaluator_truth_access\ttrue\n'
    printf 'method_families\tminiprot,gemoma,lifton\n'
    printf 'reference_lineages\tBrassica_rapa_A,Brassica_oleracea_C\n'
    printf 'within_method_reference_vote_count\t1\n'
    printf 'paired_complete_annotation_control\ttrue\n'
    printf 'candidate_merge_module_sha256\t%s\n' \
        "$(sha256sum "$code_root/src/ploidypatch/candidate_merge.py" | awk '{print $1}')"
    printf 'consensus_module_sha256\t%s\n' \
        "$(sha256sum "$code_root/src/ploidypatch/consensus.py" | awk '{print $1}')"
} > "$working_root/run_contract.tsv"
{
    printf 'role\tbytes\tsha256\tpath\n'
    for key in gemoma_brapa gemoma_bol lifton_brapa lifton_bol; do
        path=${raw_gff[$key]}
        printf '%s\t%s\t%s\t%s\n' "$key" "$(stat -Lc %s "$path")" \
            "$(sha256sum "$path" | awk '{print $1}')" "$path"
    done
    for entry in \
        "miniprot_blind:$miniprot_root/blind/candidate.gff3" \
        "miniprot_control:$miniprot_root/complete_control/candidate.gff3" \
        "blind_gff:$blind_gff" "source_gff:$source_gff"; do
        role=${entry%%:*}; path=${entry#*:}
        printf '%s\t%s\t%s\t%s\n' "$role" "$(stat -Lc %s "$path")" \
            "$(sha256sum "$path" | awk '{print $1}')" "$path"
    done
} > "$working_root/candidate_input_manifest.tsv"

cd "$code_root"
for method in gemoma lifton; do
    "$python_bin" -m ploidypatch.cli baseline merge-candidate-gffs \
        --candidate "brapa=${raw_gff[${method}_brapa]}" \
        --candidate "bol=${raw_gff[${method}_bol]}" \
        --output-gff "$working_root/merged/$method.gff3" \
        --provenance-tsv "$working_root/merged/$method.provenance.tsv" \
        > "$working_root/merged/$method.stdout.json" \
        2> "$working_root/merged/$method.stderr.log"
done

pids=(); labels=()
for method in gemoma lifton; do
    for mode in blind complete_control; do
        if [[ $mode == blind ]]; then annotation=$blind_gff; else annotation=$source_gff; fi
        (
            /usr/bin/time -v -o "$working_root/methods/$method/$mode/resource.time.txt" \
                "$python_bin" -m ploidypatch.cli baseline adapt-gff \
                    --perturbed-gff "$annotation" \
                    --candidate-gff "$working_root/merged/$method.gff3" \
                    --source "$method" \
                    --output-gff "$working_root/methods/$method/$mode/candidate.gff3" \
                    --decisions-tsv "$working_root/methods/$method/$mode/decisions.tsv" \
                    --max-existing-cds-overlap 0.2 \
                    --max-redundancy-overlap 0.5 \
                    > "$working_root/methods/$method/$mode/stdout.json" \
                    2> "$working_root/methods/$method/$mode/stderr.log"
        ) &
        pids+=("$!"); labels+=("adapt:$method:$mode")
    done
done
failed=0
for index in "${!pids[@]}"; do
    if ! wait "${pids[$index]}"; then echo "failed ${labels[$index]}" >&2; failed=1; fi
done
[[ $failed -eq 0 ]] || exit 1

pids=(); labels=()
for method in gemoma lifton; do
    (
        /usr/bin/time -v -o "$working_root/methods/$method/score.resource.time.txt" \
            "$python_bin" -m ploidypatch.cli benchmark score \
                --source-gff "$source_gff" --perturbed-gff "$blind_gff" \
                --candidate-gff "$working_root/methods/$method/blind/candidate.gff3" \
                --control-candidate-gff "$working_root/methods/$method/complete_control/candidate.gff3" \
                --truth "$truth" --include-event-details \
                > "$working_root/methods/$method/score.json" \
                2> "$working_root/methods/$method/score.stderr.log"
    ) &
    pids+=("$!"); labels+=("score:$method")
done
failed=0
for index in "${!pids[@]}"; do
    if ! wait "${pids[$index]}"; then echo "failed ${labels[$index]}" >&2; failed=1; fi
done
[[ $failed -eq 0 ]] || exit 1

else
    if [[ -e $result_root || ! -d $working_root ]]; then
        echo "consensus resume requires only the existing method-trio working root" >&2
        exit 1
    fi
    for required in \
        "$working_root/run_contract.tsv" \
        "$working_root/candidate_input_manifest.tsv" \
        "$working_root/merged/gemoma.gff3" \
        "$working_root/merged/lifton.gff3" \
        "$working_root/methods/gemoma/blind/candidate.gff3" \
        "$working_root/methods/gemoma/complete_control/candidate.gff3" \
        "$working_root/methods/gemoma/score.json" \
        "$working_root/methods/lifton/blind/candidate.gff3" \
        "$working_root/methods/lifton/complete_control/candidate.gff3" \
        "$working_root/methods/lifton/score.json" \
        "$working_root/consensus/support2/blind/candidate.gff3" \
        "$working_root/consensus/support2/complete_control/candidate.gff3" \
        "$working_root/consensus/support3/blind/candidate.gff3" \
        "$working_root/consensus/support3/complete_control/candidate.gff3"; do
        if [[ ! -s $required ]]; then
            echo "missing consensus-resume checkpoint: $required" >&2
            exit 1
        fi
    done
    mkdir -p "$working_root/consensus/union/blind" \
        "$working_root/consensus/union/complete_control"
    {
        printf 'field\tvalue\n'
        printf 'resume_code_commit\t%s\n' \
            "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
        printf 'resume_stage\tconsensus\n'
        printf 'resume_reason\treplace_quadratic_redundancy_scan_with_byte_equivalent_interval_index\n'
        printf 'consensus_module_sha256\t%s\n' \
            "$(sha256sum "$code_root/src/ploidypatch/consensus.py" | awk '{print $1}')"
        printf 'completed_support_tiers_reused\tsupport2,support3\n'
        printf 'support2_byte_equivalence_preflight\tpass\n'
    } > "$working_root/resume_contract.tsv"
    cd "$code_root"
fi

pids=(); labels=()
if [[ $resume_stage == fresh ]]; then
    consensus_tiers=(union support2 support3)
else
    consensus_tiers=(union)
fi
for tier in "${consensus_tiers[@]}"; do
    case $tier in union) support=1 ;; support2) support=2 ;; support3) support=3 ;; esac
    for mode in blind complete_control; do
        if [[ $mode == blind ]]; then base=$blind_gff; else base=$source_gff; fi
        (
            /usr/bin/time -v -o "$working_root/consensus/$tier/$mode/resource.time.txt" \
                "$python_bin" -m ploidypatch.cli baseline select-method-consensus \
                    --base-gff "$base" \
                    --candidate "miniprot=$miniprot_root/$mode/candidate.gff3" \
                    --candidate "gemoma=$working_root/methods/gemoma/$mode/candidate.gff3" \
                    --candidate "lifton=$working_root/methods/lifton/$mode/candidate.gff3" \
                    --output-gff "$working_root/consensus/$tier/$mode/candidate.gff3" \
                    --decisions-tsv "$working_root/consensus/$tier/$mode/decisions.tsv" \
                    --min-method-support "$support" --max-redundancy-overlap 0.5 \
                    > "$working_root/consensus/$tier/$mode/stdout.json" \
                    2> "$working_root/consensus/$tier/$mode/stderr.log"
        ) &
        pids+=("$!"); labels+=("consensus:$tier:$mode")
    done
done
failed=0
for index in "${!pids[@]}"; do
    if ! wait "${pids[$index]}"; then echo "failed ${labels[$index]}" >&2; failed=1; fi
done
[[ $failed -eq 0 ]] || exit 1

pids=(); labels=()
for tier in union support2 support3; do
    (
        /usr/bin/time -v -o "$working_root/consensus/$tier/score.resource.time.txt" \
            "$python_bin" -m ploidypatch.cli benchmark score \
                --source-gff "$source_gff" --perturbed-gff "$blind_gff" \
                --candidate-gff "$working_root/consensus/$tier/blind/candidate.gff3" \
                --control-candidate-gff "$working_root/consensus/$tier/complete_control/candidate.gff3" \
                --truth "$truth" --include-event-details \
                > "$working_root/consensus/$tier/score.json" \
                2> "$working_root/consensus/$tier/score.stderr.log"
    ) &
    pids+=("$!"); labels+=("consensus-score:$tier")
done
failed=0
for index in "${!pids[@]}"; do
    if ! wait "${pids[$index]}"; then echo "failed ${labels[$index]}" >&2; failed=1; fi
done
[[ $failed -eq 0 ]] || exit 1

for score in "$working_root"/methods/*/score.json "$working_root"/consensus/*/score.json; do
    if ! grep -q '"grade": "pass"' "$score"; then
        echo "Brassica method-trio evaluator quality gate failed: $score" >&2
        exit 1
    fi
done
{
    printf 'role\tbytes\tsha256\tpath\n'
    find "$working_root/methods" "$working_root/consensus" -type f \
        \( -name 'candidate.gff3' -o -name 'decisions.tsv' -o -name 'score.json' \) \
        -print0 | sort -z | while IFS= read -r -d '' path; do
            role=${path#"$working_root/"}; role=${role//\//_}
            printf '%s\t%s\t%s\t%s\n' "$role" "$(stat -Lc %s "$path")" \
                "$(sha256sum "$path" | awk '{print $1}')" "$path"
        done
} > "$working_root/output_freeze.tsv"
{
    printf 'role\tbytes\tsha256\tpath\n'
    printf 'hidden_truth\t%s\t%s\t%s\n' "$(stat -Lc %s "$truth")" \
        "$(sha256sum "$truth" | awk '{print $1}')" "$truth"
} > "$working_root/evaluator/input_manifest.tsv"
(
    cd "$working_root"
    find . -type f \( -name '*.gff3' -o -name '*.json' -o -name '*.tsv' \) \
        -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS
    sha256sum -c SHA256SUMS >/dev/null
)
du -sb "$working_root" > "$working_root/disk_bytes.txt"
mv "$working_root" "$result_root"
printf 'Brassica multireference method trio frozen: %s\n' "$result_root"
