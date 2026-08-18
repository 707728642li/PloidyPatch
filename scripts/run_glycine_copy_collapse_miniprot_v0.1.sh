#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "usage: $0 PROJECT_ROOT" >&2
    exit 2
fi

project_root=$(realpath "$1")
code_root=$project_root/code
python_bin=$project_root/envs/ploidypatch-dev/bin/python
benchmark_root=$project_root/benchmark/structure/copy_collapse_v0.1/gma_v21/annotation_copy_collapse_seed20260815
source_gff=$project_root/data/derived/syngap_glycine_v0.1/gma_complete/primary_chromosomes.gff3
projection_root=$project_root/results/baselines/miniprot_soybean_v0.1
miniprot_gff=$projection_root/raw/miniprot.gff3
protein_map=$projection_root/reference/gso_v2.map.tsv
result_root=$project_root/results/copy_collapse/miniprot_glycine_v0.1
working_root=${result_root}.working
perturbed_gff=$benchmark_root/blind/perturbed.gff3
truth=$benchmark_root/evaluator/truth/hidden_truth.json

for required in "$python_bin" "$source_gff" "$perturbed_gff" "$truth" \
                "$miniprot_gff" "$protein_map"; do
    if [[ ! -s $required ]]; then
        echo "missing or empty Glycine copy-collapse miniprot input: $required" >&2
        exit 1
    fi
done
if [[ -e $result_root || -e $working_root ]]; then
    echo "refusing to overwrite Glycine copy-collapse miniprot result: $result_root" >&2
    exit 1
fi

mkdir -p "$working_root/blind" "$working_root/complete_control"
{
    printf 'field\tvalue\n'
    printf 'code_commit\t%s\n' "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'split\theldout\n'
    printf 'candidate_truth_access\tfalse\n'
    printf 'evaluator_truth_access\ttrue\n'
    printf 'pair_tsv_candidate_access\tfalse\n'
    printf 'projection_reused\ttrue\n'
    printf 'reference_species\tGlycine_soja\n'
    printf 'min_identity\t0.5\n'
    printf 'min_query_coverage\t0.5\n'
    printf 'require_intact\ttrue\n'
    printf 'max_existing_cds_overlap\t0.2\n'
    printf 'max_redundancy_overlap\t0.5\n'
    printf 'adapter_module_sha256\t%s\n' \
        "$(sha256sum "$code_root/src/ploidypatch/baseline.py" | awk '{print $1}')"
    printf 'score_module_sha256\t%s\n' \
        "$(sha256sum "$code_root/src/ploidypatch/score.py" | awk '{print $1}')"
} > "$working_root/run_contract.tsv"
{
    printf 'role\tbytes\tsha256\tpath\n'
    for entry in \
        "source_gff:$source_gff" \
        "perturbed_gff:$perturbed_gff" \
        "miniprot_gff:$miniprot_gff" \
        "protein_map:$protein_map"; do
        role=${entry%%:*}
        path=${entry#*:}
        printf '%s\t%s\t%s\t%s\n' "$role" "$(stat -Lc %s "$path")" \
            "$(sha256sum "$path" | awk '{print $1}')" "$path"
    done
} > "$working_root/candidate_input_manifest.tsv"

cd "$code_root"
for mode in blind complete_control; do
    if [[ $mode == blind ]]; then
        annotation=$perturbed_gff
    else
        annotation=$source_gff
    fi
    /usr/bin/time -v -o "$working_root/$mode/resource.time.txt" \
        "$python_bin" -m ploidypatch.cli baseline adapt-miniprot \
            --perturbed-gff "$annotation" \
            --miniprot-gff "$miniprot_gff" \
            --protein-map "$protein_map" \
            --output-gff "$working_root/$mode/candidate.gff3" \
            --decisions-tsv "$working_root/$mode/decisions.tsv" \
            > "$working_root/$mode/stdout.json" \
            2> "$working_root/$mode/stderr.log"
    "$python_bin" -m ploidypatch.cli evidence summarize-projection-support \
        --decisions "$working_root/$mode/decisions.tsv" \
        --output-tsv "$working_root/$mode/projection_support.tsv" \
        > "$working_root/$mode/projection_support.stdout.json" \
        2> "$working_root/$mode/projection_support.stderr.log"
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
    echo "Glycine copy-collapse miniprot score quality gate failed" >&2
    exit 1
fi

for output in "$working_root/blind/candidate.gff3" \
              "$working_root/blind/decisions.tsv" \
              "$working_root/blind/projection_support.tsv" \
              "$working_root/complete_control/candidate.gff3" \
              "$working_root/complete_control/decisions.tsv" \
              "$working_root/complete_control/projection_support.tsv" \
              "$working_root/candidate_freeze.tsv" "$working_root/score.json"; do
    if [[ ! -s $output ]]; then
        echo "missing or empty Glycine copy-collapse miniprot output: $output" >&2
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
printf 'Glycine copy-collapse miniprot holdout frozen: %s\n' "$result_root"
