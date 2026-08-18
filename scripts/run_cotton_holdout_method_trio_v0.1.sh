#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then echo "usage: $0 PROJECT_ROOT" >&2; exit 2; fi
project_root=$(realpath "$1")
code_root=$project_root/code
python_bin=$project_root/envs/ploidypatch-dev/bin/python
benchmark=$project_root/benchmark/structure/copy_collapse_v0.1/ghi_ad/annotation_copy_collapse_seed20260817
blind_gff=$benchmark/blind/perturbed.gff3
truth=$benchmark/evaluator/truth/hidden_truth.json
source_gff=$project_root/data/derived/holdout_inputs/cotton_v0.1/hirsutum/primary_chromosomes.gff3
upstream=$project_root/results/baselines/cotton_holdout_v0.1
policy=$code_root/config/copy_collapse_zero_retuning_policy_v0.1.tsv
result_root=$project_root/results/copy_collapse/holdout/cotton_method_trio_v0.1
working_root=${result_root}.working
declare -A raw=(
    [gemoma_gar]="$upstream/gemoma/gar_a/upstream/final_annotation.gff"
    [gemoma_gra]="$upstream/gemoma/gra_d/upstream/final_annotation.gff"
    [lifton_gar]="$upstream/lifton/gar_a/upstream/lifton.gff3"
    [lifton_gra]="$upstream/lifton/gra_d/upstream/lifton.gff3"
)
miniprot_gff=$upstream/miniprot/raw/miniprot.gff3
protein_map=$upstream/miniprot/reference/cotton_diploids.map.tsv
for required in "$python_bin" "$blind_gff" "$truth" "$source_gff" "$policy" \
                "$miniprot_gff" "$protein_map" "${raw[@]}"; do
    if [[ ! -s $required ]]; then echo "missing cotton method input: $required" >&2; exit 1; fi
done
if [[ -e $result_root || -e $working_root ]]; then echo "refusing to overwrite cotton method trio" >&2; exit 1; fi
mkdir -p "$working_root"/{merged,evaluator} \
    "$working_root/methods"/{miniprot,gemoma,lifton}/{blind,complete_control} \
    "$working_root/consensus"/{union,support2,support3}/{blind,complete_control}
{
    printf 'field\tvalue\ncode_commit\t%s\n' "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'split\texternal_zero_retuning_holdout\ncandidate_truth_access\tfalse\n'
    printf 'evaluator_truth_access\ttrue\nmethod_families\tminiprot,gemoma,lifton\n'
    printf 'references\tGossypium_arboreum_A,Gossypium_raimondii_D\n'
    printf 'within_method_reference_vote_count\t1\nprimary_policy\texact_support_at_least_2_of_3\n'
    printf 'policy_sha256\t%s\n' "$(sha256sum "$policy" | awk '{print $1}')"
    printf 'consensus_module_sha256\t%s\n' "$(sha256sum "$code_root/src/ploidypatch/consensus.py" | awk '{print $1}')"
} > "$working_root/run_contract.tsv"
cd "$code_root"
for method in gemoma lifton; do
    "$python_bin" -m ploidypatch.cli baseline merge-candidate-gffs \
        --candidate "gar_a=${raw[${method}_gar]}" --candidate "gra_d=${raw[${method}_gra]}" \
        --output-gff "$working_root/merged/$method.gff3" \
        --provenance-tsv "$working_root/merged/$method.provenance.tsv" \
        > "$working_root/merged/$method.stdout.json" 2> "$working_root/merged/$method.stderr.log"
done

pids=(); labels=()
for mode in blind complete_control; do
    annotation=$blind_gff; [[ $mode == complete_control ]] && annotation=$source_gff
    (
        "$python_bin" -m ploidypatch.cli baseline adapt-miniprot \
            --perturbed-gff "$annotation" --miniprot-gff "$miniprot_gff" \
            --protein-map "$protein_map" \
            --output-gff "$working_root/methods/miniprot/$mode/candidate.gff3" \
            --decisions-tsv "$working_root/methods/miniprot/$mode/decisions.tsv" \
            > "$working_root/methods/miniprot/$mode/stdout.json" \
            2> "$working_root/methods/miniprot/$mode/stderr.log"
    ) & pids+=("$!"); labels+=("miniprot:$mode")
    for method in gemoma lifton; do
        (
            "$python_bin" -m ploidypatch.cli baseline adapt-gff \
                --perturbed-gff "$annotation" --candidate-gff "$working_root/merged/$method.gff3" \
                --source "$method" --max-existing-cds-overlap 0.2 --max-redundancy-overlap 0.5 \
                --output-gff "$working_root/methods/$method/$mode/candidate.gff3" \
                --decisions-tsv "$working_root/methods/$method/$mode/decisions.tsv" \
                > "$working_root/methods/$method/$mode/stdout.json" \
                2> "$working_root/methods/$method/$mode/stderr.log"
        ) & pids+=("$!"); labels+=("$method:$mode")
    done
done
failed=0
for i in "${!pids[@]}"; do if ! wait "${pids[$i]}"; then echo "failed ${labels[$i]}" >&2; failed=1; fi; done
[[ $failed -eq 0 ]] || exit 1

pids=(); labels=()
for tier in union support2 support3; do
    case $tier in union) support=1 ;; support2) support=2 ;; support3) support=3 ;; esac
    for mode in blind complete_control; do
        base=$blind_gff; [[ $mode == complete_control ]] && base=$source_gff
        (
            "$python_bin" -m ploidypatch.cli baseline select-method-consensus \
                --base-gff "$base" \
                --candidate "miniprot=$working_root/methods/miniprot/$mode/candidate.gff3" \
                --candidate "gemoma=$working_root/methods/gemoma/$mode/candidate.gff3" \
                --candidate "lifton=$working_root/methods/lifton/$mode/candidate.gff3" \
                --min-method-support "$support" --max-redundancy-overlap 0.5 \
                --output-gff "$working_root/consensus/$tier/$mode/candidate.gff3" \
                --decisions-tsv "$working_root/consensus/$tier/$mode/decisions.tsv" \
                > "$working_root/consensus/$tier/$mode/stdout.json" \
                2> "$working_root/consensus/$tier/$mode/stderr.log"
        ) & pids+=("$!"); labels+=("$tier:$mode")
    done
done
failed=0
for i in "${!pids[@]}"; do if ! wait "${pids[$i]}"; then echo "failed ${labels[$i]}" >&2; failed=1; fi; done
[[ $failed -eq 0 ]] || exit 1

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

pids=(); labels=()
for method in miniprot gemoma lifton; do
    (
        "$python_bin" -m ploidypatch.cli benchmark score --source-gff "$source_gff" \
            --perturbed-gff "$blind_gff" \
            --candidate-gff "$working_root/methods/$method/blind/candidate.gff3" \
            --control-candidate-gff "$working_root/methods/$method/complete_control/candidate.gff3" \
            --truth "$truth" --include-event-details > "$working_root/methods/$method/score.json" \
            2> "$working_root/methods/$method/score.stderr.log"
    ) & pids+=("$!"); labels+=("score:$method")
done
for tier in union support2 support3; do
    (
        "$python_bin" -m ploidypatch.cli benchmark score --source-gff "$source_gff" \
            --perturbed-gff "$blind_gff" \
            --candidate-gff "$working_root/consensus/$tier/blind/candidate.gff3" \
            --control-candidate-gff "$working_root/consensus/$tier/complete_control/candidate.gff3" \
            --truth "$truth" --include-event-details > "$working_root/consensus/$tier/score.json" \
            2> "$working_root/consensus/$tier/score.stderr.log"
    ) & pids+=("$!"); labels+=("score:$tier")
done
failed=0
for i in "${!pids[@]}"; do if ! wait "${pids[$i]}"; then echo "failed ${labels[$i]}" >&2; failed=1; fi; done
[[ $failed -eq 0 ]] || exit 1
for score in "$working_root"/methods/*/score.json "$working_root"/consensus/*/score.json; do
    grep -q '"grade": "pass"' "$score" || { echo "cotton score gate failed: $score" >&2; exit 1; }
done
(
    cd "$working_root"
    find . -type f \( -name '*.gff3' -o -name '*.json' -o -name '*.tsv' \) \
        -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS
    sha256sum -c SHA256SUMS >/dev/null
)
du -sb "$working_root" > "$working_root/disk_bytes.txt"
mv "$working_root" "$result_root"
printf 'cotton zero-retuning method trio frozen: %s\n' "$result_root"
