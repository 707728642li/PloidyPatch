#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "usage: $0 PROJECT_ROOT gemoma|lifton" >&2
    exit 2
fi

project_root=$(realpath "$1")
method=$2
code_root=$project_root/code
python_bin=$project_root/envs/ploidypatch-dev/bin/python
benchmark_root=$project_root/benchmark/heldout/v0.1/gma_v21_annotation_missing_gene_seed20260809
blind_gff=$benchmark_root/blind/perturbed.gff3
source_gff=$benchmark_root/evaluator/restoration/restored.gff3
truth=$benchmark_root/evaluator/truth/hidden_truth.json
strata=$benchmark_root/evaluator/catalog/gma_v21_primary.annotation_missing_gene.tsv

case $method in
    gemoma)
        baseline_root=$project_root/results/baselines/gemoma_v1.9/glycine_v0.1
        upstream_gff=$baseline_root/raw/upstream/final_annotation.gff
        source_label=gemoma_gso
        ;;
    lifton)
        baseline_root=$project_root/results/baselines/lifton_v1.0.11/glycine_v0.1
        upstream_gff=$baseline_root/raw/upstream/lifton.gff3
        source_label=lifton_gso
        ;;
    *)
        echo "method must be gemoma or lifton" >&2
        exit 2
        ;;
esac

result_root=$baseline_root/evaluation
working_root=${result_root}.working
for required in "$python_bin" "$upstream_gff" "$blind_gff" \
                "$source_gff" "$truth" "$strata"; do
    if [[ ! -s $required ]]; then
        echo "missing or empty evaluation input: $required" >&2
        exit 1
    fi
done
if [[ -e $result_root || -e $working_root ]]; then
    echo "refusing to overwrite evaluation result: $result_root" >&2
    exit 1
fi

mkdir -p "$working_root/blind" "$working_root/complete_control"
{
    printf 'role\tbytes\tsha256\tpath\n'
    for entry in \
        "upstream_gff:$upstream_gff" \
        "blind_gff:$blind_gff" \
        "source_gff:$source_gff" \
        "hidden_truth:$truth" \
        "event_strata:$strata"; do
        role=${entry%%:*}
        path=${entry#*:}
        printf '%s\t%s\t%s\t%s\n' "$role" "$(stat -Lc %s "$path")" \
            "$(sha256sum "$path" | awk '{print $1}')" "$path"
    done
} > "$working_root/input_manifest.tsv"
{
    printf 'field\tvalue\n'
    printf 'method\t%s\n' "$method"
    printf 'candidate_source\t%s\n' "$source_label"
    printf 'max_existing_cds_overlap\t0.2\n'
    printf 'max_redundancy_overlap\t0.5\n'
    printf 'paired_complete_annotation_control\ttrue\n'
} > "$working_root/run_contract.tsv"

cd "$code_root"
/usr/bin/time -v -o "$working_root/blind/resource.time.txt" \
    "$python_bin" -m ploidypatch.cli baseline adapt-gff \
        --perturbed-gff "$blind_gff" \
        --candidate-gff "$upstream_gff" \
        --source "$source_label" \
        --output-gff "$working_root/blind/candidate.gff3" \
        --decisions-tsv "$working_root/blind/decisions.tsv" \
        --max-existing-cds-overlap 0.2 \
        --max-redundancy-overlap 0.5 \
        > "$working_root/blind/stdout.log" \
        2> "$working_root/blind/stderr.log"

/usr/bin/time -v -o "$working_root/complete_control/resource.time.txt" \
    "$python_bin" -m ploidypatch.cli baseline adapt-gff \
        --perturbed-gff "$source_gff" \
        --candidate-gff "$upstream_gff" \
        --source "$source_label" \
        --output-gff "$working_root/complete_control/candidate.gff3" \
        --decisions-tsv "$working_root/complete_control/decisions.tsv" \
        --max-existing-cds-overlap 0.2 \
        --max-redundancy-overlap 0.5 \
        > "$working_root/complete_control/stdout.log" \
        2> "$working_root/complete_control/stderr.log"

/usr/bin/time -v -o "$working_root/score.resource.time.txt" \
    "$python_bin" -m ploidypatch.cli benchmark score \
        --source-gff "$source_gff" \
        --perturbed-gff "$blind_gff" \
        --candidate-gff "$working_root/blind/candidate.gff3" \
        --control-candidate-gff \
            "$working_root/complete_control/candidate.gff3" \
        --truth "$truth" \
        --event-strata "$strata" \
        --stratum-column transcript_count_bin \
        --stratum-column max_exons_bin \
        --include-event-details \
        > "$working_root/score.json" \
        2> "$working_root/score.stderr.log"

for output in \
    "$working_root/blind/candidate.gff3" \
    "$working_root/blind/decisions.tsv" \
    "$working_root/complete_control/candidate.gff3" \
    "$working_root/complete_control/decisions.tsv" \
    "$working_root/score.json"; do
    if [[ ! -s $output ]]; then
        echo "evaluation output is missing or empty: $output" >&2
        exit 1
    fi
done
{
    printf 'role\tbytes\tsha256\tpath\n'
    for output in \
        "$working_root/blind/candidate.gff3" \
        "$working_root/blind/decisions.tsv" \
        "$working_root/complete_control/candidate.gff3" \
        "$working_root/complete_control/decisions.tsv" \
        "$working_root/score.json"; do
        printf '%s\t%s\t%s\t%s\n' "$(basename "$output")" \
            "$(stat -Lc %s "$output")" \
            "$(sha256sum "$output" | awk '{print $1}')" "$output"
    done
} > "$working_root/output_manifest.tsv"
du -sb "$working_root" > "$working_root/disk_bytes.txt"
mv "$working_root" "$result_root"
printf '%s evaluation validated: %s\n' "$method" "$result_root"
