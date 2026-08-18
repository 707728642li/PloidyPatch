#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then echo "usage: $0 PROJECT_ROOT" >&2; exit 2; fi
project_root=$(realpath "$1")
code_root=$project_root/code
python_bin=$project_root/envs/ploidypatch-dev/bin/python
benchmark=$project_root/benchmark/structure/copy_collapse_v0.2/zma_maize1/annotation_copy_collapse_seed20260829
blind_gff=$benchmark/blind/perturbed.gff3
truth=$benchmark/evaluator/truth/hidden_truth.json
source_gff=$project_root/data/derived/holdout_inputs/maize_v2/zea_mays/primary_chromosomes.gff3
upstream=$project_root/results/baselines/maize_v2
blind_root=$project_root/results/copy_collapse/holdout/maize_v2_method_trio
result_root=$project_root/results/copy_collapse/holdout/maize_v2_method_trio_evaluation
working_root=${result_root}.working
declare -A raw=(
    [gemoma_sorghum]="$upstream/gemoma/sorghum_bicolor/upstream/final_annotation.gff"
    [gemoma_setaria]="$upstream/gemoma/setaria_italica/upstream/final_annotation.gff"
    [lifton_sorghum]="$upstream/lifton/sorghum_bicolor/upstream/lifton.gff3"
    [lifton_setaria]="$upstream/lifton/setaria_italica/upstream/lifton.gff3"
)
miniprot_gff=$upstream/miniprot/raw/miniprot.gff3
protein_map=$upstream/miniprot/reference/maize_outgroups.map.tsv

for required in "$python_bin" "$blind_gff" "$truth" "$source_gff" \
                "$blind_root/SHA256SUMS" "$miniprot_gff" "$protein_map" \
                "${raw[@]}"; do
    [[ -s $required ]] || { echo "missing maize evaluator input: $required" >&2; exit 1; }
done
(cd "$blind_root" && sha256sum -c SHA256SUMS >/dev/null)
if [[ -e $result_root || -e $working_root ]]; then
    echo "refusing to overwrite maize method evaluation" >&2; exit 1
fi
mkdir -p "$working_root/merged" "$working_root/controls/methods"/{miniprot,gemoma,lifton} \
    "$working_root/controls/consensus"/{union,support2,support3} \
    "$working_root/scores/methods"/{miniprot,gemoma,lifton} \
    "$working_root/scores/consensus"/{union,support2,support3}
{
    printf 'field\tvalue\ncode_commit\t%s\n' \
        "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'access\tevaluator_only_after_candidate_freeze\n'
    printf 'candidate_freeze_sha256\t%s\n' \
        "$(sha256sum "$blind_root/candidate_freeze.tsv" | awk '{print $1}')"
    printf 'paired_complete_annotation_control\ttrue\n'
    printf 'hidden_event_count\t800\nautomatic_approval\tfalse\n'
} > "$working_root/run_contract.tsv"
{
    printf 'role\tbytes\tsha256\tpath\n'
    for path in "$blind_root/SHA256SUMS" "$blind_root/candidate_freeze.tsv" \
                "$truth" "$source_gff"; do
        printf '%s\t%s\t%s\t%s\n' "$(basename "$path")" \
            "$(stat -Lc %s "$path")" "$(sha256sum "$path" | awk '{print $1}')" "$path"
    done
} > "$working_root/input_manifest.tsv"

cd "$code_root"
for method in gemoma lifton; do
    "$python_bin" -m ploidypatch.cli baseline merge-candidate-gffs \
        --candidate "sorghum=${raw[${method}_sorghum]}" \
        --candidate "setaria=${raw[${method}_setaria]}" \
        --output-gff "$working_root/merged/$method.gff3" \
        --provenance-tsv "$working_root/merged/$method.provenance.tsv" \
        > "$working_root/merged/$method.stdout.json" \
        2> "$working_root/merged/$method.stderr.log"
done

pids=(); labels=()
(
    "$python_bin" -m ploidypatch.cli baseline adapt-miniprot \
        --perturbed-gff "$source_gff" --miniprot-gff "$miniprot_gff" \
        --protein-map "$protein_map" \
        --output-gff "$working_root/controls/methods/miniprot/candidate.gff3" \
        --decisions-tsv "$working_root/controls/methods/miniprot/decisions.tsv" \
        > "$working_root/controls/methods/miniprot/stdout.json" \
        2> "$working_root/controls/methods/miniprot/stderr.log"
) & pids+=("$!"); labels+=("miniprot")
for method in gemoma lifton; do
    (
        "$python_bin" -m ploidypatch.cli baseline adapt-gff \
            --perturbed-gff "$source_gff" \
            --candidate-gff "$working_root/merged/$method.gff3" \
            --source "$method" --max-existing-cds-overlap 0.2 \
            --max-redundancy-overlap 0.5 \
            --output-gff "$working_root/controls/methods/$method/candidate.gff3" \
            --decisions-tsv "$working_root/controls/methods/$method/decisions.tsv" \
            > "$working_root/controls/methods/$method/stdout.json" \
            2> "$working_root/controls/methods/$method/stderr.log"
    ) & pids+=("$!"); labels+=("$method")
done
failed=0
for i in "${!pids[@]}"; do
    if ! wait "${pids[$i]}"; then echo "failed control adapt: ${labels[$i]}" >&2; failed=1; fi
done
[[ $failed -eq 0 ]] || exit 1

pids=(); labels=()
for tier in union support2 support3; do
    case $tier in union) support=1 ;; support2) support=2 ;; support3) support=3 ;; esac
    (
        "$python_bin" -m ploidypatch.cli baseline select-method-consensus \
            --base-gff "$source_gff" \
            --candidate "miniprot=$working_root/controls/methods/miniprot/candidate.gff3" \
            --candidate "gemoma=$working_root/controls/methods/gemoma/candidate.gff3" \
            --candidate "lifton=$working_root/controls/methods/lifton/candidate.gff3" \
            --min-method-support "$support" --max-redundancy-overlap 0.5 \
            --output-gff "$working_root/controls/consensus/$tier/candidate.gff3" \
            --decisions-tsv "$working_root/controls/consensus/$tier/decisions.tsv" \
            > "$working_root/controls/consensus/$tier/stdout.json" \
            2> "$working_root/controls/consensus/$tier/stderr.log"
    ) & pids+=("$!"); labels+=("$tier")
done
failed=0
for i in "${!pids[@]}"; do
    if ! wait "${pids[$i]}"; then echo "failed control consensus: ${labels[$i]}" >&2; failed=1; fi
done
[[ $failed -eq 0 ]] || exit 1

{
    printf 'role\tbytes\tsha256\tpath\n'
    find "$working_root/controls" -type f \
        \( -name candidate.gff3 -o -name decisions.tsv \) -print0 | sort -z \
        | while IFS= read -r -d '' path; do
            role=${path#"$working_root/"}; role=${role//\//_}
            printf '%s\t%s\t%s\t%s\n' "$role" "$(stat -Lc %s "$path")" \
                "$(sha256sum "$path" | awk '{print $1}')" "$path"
        done
} > "$working_root/control_freeze.tsv"

pids=(); labels=()
for method in miniprot gemoma lifton; do
    (
        "$python_bin" -m ploidypatch.cli benchmark score \
            --source-gff "$source_gff" --perturbed-gff "$blind_gff" \
            --candidate-gff "$blind_root/methods/$method/blind/candidate.gff3" \
            --control-candidate-gff "$working_root/controls/methods/$method/candidate.gff3" \
            --truth "$truth" --include-event-details \
            > "$working_root/scores/methods/$method/score.json" \
            2> "$working_root/scores/methods/$method/stderr.log"
    ) & pids+=("$!"); labels+=("method:$method")
done
for tier in union support2 support3; do
    (
        "$python_bin" -m ploidypatch.cli benchmark score \
            --source-gff "$source_gff" --perturbed-gff "$blind_gff" \
            --candidate-gff "$blind_root/consensus/$tier/blind/candidate.gff3" \
            --control-candidate-gff "$working_root/controls/consensus/$tier/candidate.gff3" \
            --truth "$truth" --include-event-details \
            > "$working_root/scores/consensus/$tier/score.json" \
            2> "$working_root/scores/consensus/$tier/stderr.log"
    ) & pids+=("$!"); labels+=("consensus:$tier")
done
failed=0
for i in "${!pids[@]}"; do
    if ! wait "${pids[$i]}"; then echo "failed score: ${labels[$i]}" >&2; failed=1; fi
done
[[ $failed -eq 0 ]] || exit 1
for score in "$working_root"/scores/methods/*/score.json \
             "$working_root"/scores/consensus/*/score.json; do
    grep -q '"grade": "pass"' "$score" || {
        echo "maize score gate failed: $score" >&2; exit 1;
    }
done
(
    cd "$working_root"
    find . -type f \( -name '*.gff3' -o -name '*.json' -o -name '*.tsv' \) \
        -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS
    sha256sum -c SHA256SUMS >/dev/null
)
du -sb "$working_root" > "$working_root/disk_bytes.txt"
mv "$working_root" "$result_root"
printf 'maize method trio evaluation frozen: %s\n' "$result_root"
