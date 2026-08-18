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
source_gff=$project_root/data/derived/normalized_bundles/v0.1/bna_daae_primary/primary_chromosomes.gff3
baseline_root=$project_root/results/copy_collapse/miniprot_brassica_v0.1
result_root=$project_root/results/copy_collapse/miniprot_brassica_support2_v0.1
working_root=${result_root}.working
perturbed_gff=$benchmark_root/blind/perturbed.gff3
truth=$benchmark_root/evaluator/truth/hidden_truth.json

for required in "$python_bin" "$source_gff" "$perturbed_gff" "$truth" \
                "$baseline_root/blind/candidate.gff3" \
                "$baseline_root/blind/projection_support.tsv" \
                "$baseline_root/complete_control/candidate.gff3" \
                "$baseline_root/complete_control/projection_support.tsv"; do
    if [[ ! -s $required ]]; then
        echo "missing or empty support-tier input: $required" >&2
        exit 1
    fi
done
if [[ -e $result_root || -e $working_root ]]; then
    echo "refusing to overwrite support-tier result: $result_root" >&2
    exit 1
fi

mkdir -p "$working_root/blind" "$working_root/complete_control"
{
    printf 'field\tvalue\n'
    printf 'code_commit\t%s\n' "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'candidate_truth_access\tfalse\n'
    printf 'evaluator_truth_access\ttrue\n'
    printf 'pair_tsv_candidate_access\tfalse\n'
    printf 'min_support_group_count\t2\n'
    printf 'support_groups\tbra_a,bol_c\n'
    printf 'support_unit\tindependent_reference_species\n'
    printf 'threshold_provenance\tprescored_disjoint_development_exploration\n'
    printf 'selector_module_sha256\t%s\n' \
        "$(sha256sum "$code_root/src/ploidypatch/projection_select.py" | awk '{print $1}')"
} > "$working_root/run_contract.tsv"

cd "$code_root"
for mode in blind complete_control; do
    /usr/bin/time -v -o "$working_root/$mode/resource.time.txt" \
        "$python_bin" -m ploidypatch.cli evidence select-projection-support \
            --candidate-gff "$baseline_root/$mode/candidate.gff3" \
            --projection-support "$baseline_root/$mode/projection_support.tsv" \
            --min-support-group-count 2 \
            --output-gff "$working_root/$mode/candidate.gff3" \
            --selection-tsv "$working_root/$mode/selection.tsv" \
            > "$working_root/$mode/stdout.json" \
            2> "$working_root/$mode/stderr.log"
done

{
    printf 'role\tbytes\tsha256\tpath\n'
    for mode in blind complete_control; do
        path=$working_root/$mode/candidate.gff3
        printf '%s_candidate\t%s\t%s\t%s\n' "$mode" "$(stat -Lc %s "$path")" \
            "$(sha256sum "$path" | awk '{print $1}')" "$path"
    done
} > "$working_root/candidate_freeze.tsv"

/usr/bin/time -v -o "$working_root/score.resource.time.txt" \
    "$python_bin" -m ploidypatch.cli benchmark score \
        --source-gff "$source_gff" \
        --perturbed-gff "$perturbed_gff" \
        --candidate-gff "$working_root/blind/candidate.gff3" \
        --control-candidate-gff "$working_root/complete_control/candidate.gff3" \
        --truth "$truth" \
        --include-event-details \
        > "$working_root/score.json" \
        2> "$working_root/score.stderr.log"
if ! grep -q '"grade": "pass"' "$working_root/score.json"; then
    echo "copy-collapse support-tier score quality gate failed" >&2
    exit 1
fi

for output in "$working_root/blind/candidate.gff3" \
              "$working_root/blind/selection.tsv" \
              "$working_root/complete_control/candidate.gff3" \
              "$working_root/complete_control/selection.tsv" \
              "$working_root/candidate_freeze.tsv" \
              "$working_root/score.json"; do
    if [[ ! -s $output ]]; then
        echo "missing or empty copy-collapse support-tier output: $output" >&2
        exit 1
    fi
done
(
    cd "$working_root"
    find . -type f \( -name '*.gff3' -o -name '*.json' -o -name '*.tsv' \) \
        -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS
)
du -sb "$working_root" > "$working_root/disk_bytes.txt"
mv "$working_root" "$result_root"
printf 'copy-collapse independent-support tier frozen: %s\n' "$result_root"
