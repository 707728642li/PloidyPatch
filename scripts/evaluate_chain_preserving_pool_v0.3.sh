#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "usage: $0 PROJECT_ROOT {glycine|brassica|cotton|maize}" >&2
    exit 2
fi
project_root=$(realpath "$1")
species=$2
code_root=$project_root/code
python_bin=$project_root/envs/ploidypatch-dev/bin/python

case $species in
glycine)
    split=development_labels_previously_seen
    benchmark=$project_root/benchmark/structure/copy_collapse_v0.1/gma_v21/annotation_copy_collapse_seed20260815
    source_gff=$project_root/data/derived/syngap_glycine_v0.1/gma_complete/primary_chromosomes.gff3
    miniprot_root=$project_root/results/copy_collapse/miniprot_glycine_v0.1
    gemoma_root=$project_root/results/copy_collapse/gemoma_glycine_v0.1
    lifton_root=$project_root/results/copy_collapse/lifton_glycine_v0.1
    blind_miniprot=$miniprot_root/blind/candidate.gff3
    blind_gemoma=$gemoma_root/blind/candidate.gff3
    blind_lifton=$lifton_root/blind/candidate.gff3
    control_miniprot=$miniprot_root/complete_control/candidate.gff3
    control_gemoma=$gemoma_root/complete_control/candidate.gff3
    control_lifton=$lifton_root/complete_control/candidate.gff3
    old_score=$project_root/results/copy_collapse/consensus_union_glycine_dev_v0.1/score.json
    ;;
brassica)
    split=development_labels_previously_seen
    benchmark=$project_root/benchmark/structure/copy_collapse_v0.1/bna_daae/annotation_copy_collapse_seed20260814
    source_gff=$project_root/data/derived/normalized_bundles/v0.1/bna_daae_primary/primary_chromosomes.gff3
    method_root=$project_root/results/copy_collapse/model_development/brassica_method_trio_v0.1
    miniprot_root=$project_root/results/copy_collapse/miniprot_brassica_v0.1
    blind_miniprot=$miniprot_root/blind/candidate.gff3
    blind_gemoma=$method_root/methods/gemoma/blind/candidate.gff3
    blind_lifton=$method_root/methods/lifton/blind/candidate.gff3
    control_miniprot=$miniprot_root/complete_control/candidate.gff3
    control_gemoma=$method_root/methods/gemoma/complete_control/candidate.gff3
    control_lifton=$method_root/methods/lifton/complete_control/candidate.gff3
    old_score=$method_root/consensus/union/score.json
    ;;
cotton)
    split=retrospective_diagnostic_labels_previously_seen
    benchmark=$project_root/benchmark/structure/copy_collapse_v0.1/ghi_ad/annotation_copy_collapse_seed20260817
    source_gff=$project_root/data/derived/holdout_inputs/cotton_v0.1/hirsutum/primary_chromosomes.gff3
    method_root=$project_root/results/copy_collapse/holdout/cotton_method_trio_v0.1
    blind_miniprot=$method_root/methods/miniprot/blind/candidate.gff3
    blind_gemoma=$method_root/methods/gemoma/blind/candidate.gff3
    blind_lifton=$method_root/methods/lifton/blind/candidate.gff3
    control_miniprot=$method_root/methods/miniprot/complete_control/candidate.gff3
    control_gemoma=$method_root/methods/gemoma/complete_control/candidate.gff3
    control_lifton=$method_root/methods/lifton/complete_control/candidate.gff3
    old_score=$method_root/consensus/union/score.json
    ;;
maize)
    split=posthoc_formal_holdout_diagnostic_no_selection
    benchmark=$project_root/benchmark/structure/copy_collapse_v0.2/zma_maize1/annotation_copy_collapse_seed20260829
    source_gff=$project_root/data/derived/holdout_inputs/maize_v2/zea_mays/primary_chromosomes.gff3
    blind_root=$project_root/results/copy_collapse/holdout/maize_v2_method_trio
    control_root=$project_root/results/copy_collapse/holdout/maize_v2_method_trio_evaluation/controls/methods
    blind_miniprot=$blind_root/methods/miniprot/blind/candidate.gff3
    blind_gemoma=$blind_root/methods/gemoma/blind/candidate.gff3
    blind_lifton=$blind_root/methods/lifton/blind/candidate.gff3
    control_miniprot=$control_root/miniprot/candidate.gff3
    control_gemoma=$control_root/gemoma/candidate.gff3
    control_lifton=$control_root/lifton/candidate.gff3
    old_score=$project_root/results/copy_collapse/holdout/maize_v2_method_trio_evaluation/scores/consensus/union/score.json
    ;;
*) echo "unsupported species: $species" >&2; exit 2 ;;
esac

blind_gff=$benchmark/blind/perturbed.gff3
truth=$benchmark/evaluator/truth/hidden_truth.json
result_root=$project_root/results/copy_collapse/candidate_pool_v0.3/$species
working_root=${result_root}.working
for required in "$python_bin" "$source_gff" "$blind_gff" "$truth" "$old_score" \
                "$blind_miniprot" "$blind_gemoma" "$blind_lifton" \
                "$control_miniprot" "$control_gemoma" "$control_lifton"; do
    [[ -s $required ]] || { echo "missing candidate-pool input: $required" >&2; exit 1; }
done
if [[ -e $result_root || -e $working_root ]]; then
    echo "refusing to overwrite candidate-pool result: $result_root" >&2; exit 1
fi
mkdir -p "$working_root"/{blind,complete_control,evaluator}
{
    printf 'field\tvalue\n'
    printf 'code_commit\t%s\n' "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'species\t%s\nsplit\t%s\n' "$species" "$split"
    printf 'candidate_policy\tretain_distinct_phased_cds_chains\n'
    printf 'exact_chain_deduplication\ttrue\n'
    printf 'overlapping_alternative_suppression\tfalse\n'
    printf 'conflict_overlap_threshold\t0.5\n'
    printf 'automatic_approval\tfalse\n'
    printf 'formal_holdout_claim_allowed\tfalse\n'
} > "$working_root/run_contract.tsv"
{
    printf 'role\tbytes\tsha256\tpath\n'
    for entry in "source_gff:$source_gff" "blind_gff:$blind_gff" \
        "truth:$truth" "old_union_score:$old_score" \
        "blind_miniprot:$blind_miniprot" "blind_gemoma:$blind_gemoma" \
        "blind_lifton:$blind_lifton" "control_miniprot:$control_miniprot" \
        "control_gemoma:$control_gemoma" "control_lifton:$control_lifton"; do
        role=${entry%%:*}; path=${entry#*:}
        printf '%s\t%s\t%s\t%s\n' "$role" "$(stat -Lc %s "$path")" \
            "$(sha256sum "$path" | awk '{print $1}')" "$path"
    done
} > "$working_root/input_manifest.tsv"

cd "$code_root"
pids=(); labels=()
for mode in blind complete_control; do
    base=$blind_gff
    miniprot=$blind_miniprot; gemoma=$blind_gemoma; lifton=$blind_lifton
    if [[ $mode == complete_control ]]; then
        base=$source_gff
        miniprot=$control_miniprot; gemoma=$control_gemoma; lifton=$control_lifton
    fi
    (
        /usr/bin/time -v -o "$working_root/$mode/resource.time.txt" \
            "$python_bin" -m ploidypatch.cli baseline select-method-consensus \
            --base-gff "$base" \
            --candidate "miniprot=$miniprot" \
            --candidate "gemoma=$gemoma" \
            --candidate "lifton=$lifton" \
            --min-method-support 1 --max-redundancy-overlap 0.5 \
            --redundancy-policy retain_distinct_chains \
            --output-gff "$working_root/$mode/candidate.gff3" \
            --decisions-tsv "$working_root/$mode/decisions.tsv" \
            > "$working_root/$mode/stdout.json" \
            2> "$working_root/$mode/stderr.log"
    ) & pids+=("$!"); labels+=("$mode")
done
failed=0
for index in "${!pids[@]}"; do
    if ! wait "${pids[$index]}"; then
        echo "candidate-pool build failed: ${labels[$index]}" >&2; failed=1
    fi
done
[[ $failed -eq 0 ]] || exit 1

"$python_bin" -m ploidypatch.cli benchmark score \
    --source-gff "$source_gff" --perturbed-gff "$blind_gff" \
    --candidate-gff "$working_root/blind/candidate.gff3" \
    --control-candidate-gff "$working_root/complete_control/candidate.gff3" \
    --truth "$truth" --include-event-details \
    > "$working_root/evaluator/score.json" \
    2> "$working_root/evaluator/score.stderr.log"
grep -q '"grade": "pass"' "$working_root/evaluator/score.json" || {
    echo "candidate-pool score gate failed: $species" >&2; exit 1;
}
cp "$old_score" "$working_root/evaluator/old_suppressing_union_score.json"
(
    cd "$working_root"
    find . -type f \( -name '*.gff3' -o -name '*.json' -o -name '*.tsv' \) \
        -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS
    sha256sum -c SHA256SUMS >/dev/null
)
du -sb "$working_root" > "$working_root/disk_bytes.txt"
mv "$working_root" "$result_root"
printf 'chain-preserving candidate pool frozen: %s\n' "$result_root"
