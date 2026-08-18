#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then echo "usage: $0 PROJECT_ROOT" >&2; exit 2; fi
project_root=$(realpath "$1")
code_root=$project_root/code
python_bin=$project_root/envs/ploidypatch-dev/bin/python
benchmark_root=$project_root/benchmark/structure/copy_collapse_v0.1/bna_daae/annotation_copy_collapse_seed20260814
blind_gff=$benchmark_root/blind/perturbed.gff3
truth=$benchmark_root/evaluator/truth/hidden_truth.json
source_gff=$project_root/data/derived/normalized_bundles/v0.1/bna_daae_primary/primary_chromosomes.gff3
method_root=$project_root/results/copy_collapse/model_development/brassica_method_trio_v0.1
self_wgd_root=$project_root/results/copy_collapse/model_development/brassica_union_self_wgd_v0.1
miniprot_root=$project_root/results/copy_collapse/miniprot_brassica_v0.1
model_root=$project_root/results/copy_collapse/model_development/glycine_copy_model_v0.1
model_json=$model_root/model.json
result_root=$project_root/results/copy_collapse/model_development/brassica_glycine_model_transfer_v0.1
working_root=${result_root}.working
resume_stage=${PLOIDYPATCH_RESUME_STAGE:-}

for frozen in "$method_root" "$model_root" "$self_wgd_root/blind" "$self_wgd_root/complete_control"; do
    if [[ ! -s $frozen/SHA256SUMS ]]; then echo "unfrozen input: $frozen" >&2; exit 1; fi
    (cd "$frozen" && sha256sum -c SHA256SUMS >/dev/null)
done
for required in "$python_bin" "$blind_gff" "$truth" "$source_gff" "$model_json" \
    "$miniprot_root/blind/decisions.tsv" "$miniprot_root/complete_control/decisions.tsv"; do
    if [[ ! -s $required ]]; then echo "missing Brassica transfer input: $required" >&2; exit 1; fi
done
if [[ -e $result_root ]]; then
    echo "refusing to overwrite Brassica model transfer" >&2; exit 1
fi
if [[ -e $working_root && $resume_stage != evaluator ]]; then
    echo "working Brassica model transfer exists; set PLOIDYPATCH_RESUME_STAGE=evaluator after validating the checkpoint" >&2
    exit 1
fi
mkdir -p "$working_root/blind" "$working_root/complete_control" "$working_root/evaluator"
if [[ $resume_stage == evaluator ]]; then
    for mode in blind complete_control; do
        for checkpoint in features.tsv scores.tsv candidate.gff3 selection.tsv; do
            if [[ ! -s $working_root/$mode/$checkpoint ]]; then
                echo "missing evaluator-resume checkpoint: $working_root/$mode/$checkpoint" >&2
                exit 1
            fi
        done
    done
    {
        printf 'field\tvalue\n'
        printf 'resume_stage\tevaluator\n'
        printf 'resume_code_commit\t%s\n' "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
        printf 'preserved_candidate_scoring\ttrue\n'
        printf 'failed_command_cause\tlegacy_cli_argument_hidden_truth\n'
        for mode in blind complete_control; do
            for checkpoint in features.tsv scores.tsv candidate.gff3 selection.tsv; do
                path=$working_root/$mode/$checkpoint
                printf '%s_%s_sha256\t%s\n' "$mode" "${checkpoint//./_}" \
                    "$(sha256sum "$path" | awk '{print $1}')"
            done
        done
    } > "$working_root/resume_contract.tsv"
fi
if [[ $resume_stage != evaluator ]]; then
{
    printf 'field\tvalue\n'
    printf 'code_commit\t%s\n' "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'source_model_species\tGlycine_max\n'
    printf 'target_species\tBrassica_napus\n'
    printf 'parameter_retuning\tfalse\n'
    printf 'model_sha256\t%s\n' "$(sha256sum "$model_json" | awk '{print $1}')"
    printf 'split\tsecondary_portability_diagnostic_not_pristine_holdout\n'
    printf 'formal_external_holdout_claim_allowed\tfalse\n'
    printf 'candidate_truth_access\tfalse\n'
    printf 'evaluator_truth_access\ttrue\n'
    printf 'paired_complete_annotation_control\ttrue\n'
    printf 'automatic_approval\tfalse\n'
} > "$working_root/run_contract.tsv"
{
    printf 'role\tbytes\tsha256\tpath\n'
    for entry in "model:$model_json" "truth:$truth" "blind_gff:$blind_gff" \
                 "source_gff:$source_gff"; do
        role=${entry%%:*}; path=${entry#*:}
        printf '%s\t%s\t%s\t%s\n' "$role" "$(stat -Lc %s "$path")" \
            "$(sha256sum "$path" | awk '{print $1}')" "$path"
    done
} > "$working_root/input_manifest.tsv"

cd "$code_root"
for mode in blind complete_control; do
    "$python_bin" -m ploidypatch.cli evidence build-copy-features \
        --consensus-decisions "$method_root/consensus/union/$mode/decisions.tsv" \
        --method-decisions "miniprot=$miniprot_root/$mode/decisions.tsv" \
        --method-decisions "gemoma=$method_root/methods/gemoma/$mode/decisions.tsv" \
        --method-decisions "lifton=$method_root/methods/lifton/$mode/decisions.tsv" \
        --wgd-selection "$self_wgd_root/$mode/selected/selection.tsv" \
        --output-tsv "$working_root/$mode/features.tsv" \
        > "$working_root/$mode/features.stdout.json" \
        2> "$working_root/$mode/features.stderr.log"
    "$python_bin" -m ploidypatch.cli evidence score-copy-candidates \
        --features "$working_root/$mode/features.tsv" --model-json "$model_json" \
        --output-tsv "$working_root/$mode/scores.tsv" \
        > "$working_root/$mode/scores.stdout.json" \
        2> "$working_root/$mode/scores.stderr.log"
    if [[ $mode == blind ]]; then base=$blind_gff; else base=$source_gff; fi
    "$python_bin" -m ploidypatch.cli evidence select-scored-copy-candidates \
        --base-gff "$base" \
        --candidate-gff "$method_root/consensus/union/$mode/candidate.gff3" \
        --scores "$working_root/$mode/scores.tsv" --model-json "$model_json" \
        --policy review --output-gff "$working_root/$mode/candidate.gff3" \
        --selection-tsv "$working_root/$mode/selection.tsv" \
        > "$working_root/$mode/selection.stdout.json" \
        2> "$working_root/$mode/selection.stderr.log"
done
else
    cd "$code_root"
fi
"$python_bin" -m ploidypatch.cli benchmark label-copy-features \
    --features "$working_root/blind/features.tsv" --truth "$truth" \
    --output-tsv "$working_root/evaluator/labeled_features.tsv" \
    > "$working_root/evaluator/labels.stdout.json" \
    2> "$working_root/evaluator/labels.stderr.log"
"$python_bin" -m ploidypatch.cli benchmark score-copy-ranking \
    --scores "$working_root/blind/scores.tsv" \
    --labeled-features "$working_root/evaluator/labeled_features.tsv" \
    --output-json "$working_root/evaluator/ranking.json" \
    > "$working_root/evaluator/ranking.stdout.json" \
    2> "$working_root/evaluator/ranking.stderr.log"
/usr/bin/time -v -o "$working_root/evaluator/paired_score.resource.time.txt" \
    "$python_bin" -m ploidypatch.cli benchmark score \
        --source-gff "$source_gff" --perturbed-gff "$blind_gff" \
        --candidate-gff "$working_root/blind/candidate.gff3" \
        --control-candidate-gff "$working_root/complete_control/candidate.gff3" \
        --truth "$truth" --include-event-details \
        > "$working_root/evaluator/paired_score.json" \
        2> "$working_root/evaluator/paired_score.stderr.log"
if ! grep -q '"grade": "pass"' "$working_root/evaluator/paired_score.json"; then
    echo "Brassica transfer paired score gate failed" >&2; exit 1
fi
{
    printf 'role\tbytes\tsha256\tpath\n'
    find "$working_root/blind" "$working_root/complete_control" "$working_root/evaluator" \
        -type f \( -name '*.tsv' -o -name '*.json' -o -name '*.gff3' \) \
        -print0 | sort -z | while IFS= read -r -d '' path; do
            role=${path#"$working_root/"}; role=${role//\//_}
            printf '%s\t%s\t%s\t%s\n' "$role" "$(stat -Lc %s "$path")" \
                "$(sha256sum "$path" | awk '{print $1}')" "$path"
        done
} > "$working_root/output_freeze.tsv"
(
    cd "$working_root"
    find . -type f \( -name '*.gff3' -o -name '*.json' -o -name '*.tsv' \) \
        -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS
    sha256sum -c SHA256SUMS >/dev/null
)
du -sb "$working_root" > "$working_root/disk_bytes.txt"
mv "$working_root" "$result_root"
printf 'Brassica zero-retuning Glycine model diagnostic frozen: %s\n' "$result_root"
