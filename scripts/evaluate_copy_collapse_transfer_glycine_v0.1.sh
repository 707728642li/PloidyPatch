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
benchmark_root=$project_root/benchmark/structure/copy_collapse_v0.1/gma_v21/annotation_copy_collapse_seed20260815
blind_gff=$benchmark_root/blind/perturbed.gff3
source_gff=$project_root/data/derived/syngap_glycine_v0.1/gma_complete/primary_chromosomes.gff3
truth=$benchmark_root/evaluator/truth/hidden_truth.json

case $method in
    gemoma)
        upstream_gff=$project_root/results/baselines/gemoma_v1.9/glycine_v0.1/raw/upstream/final_annotation.gff
        source_label=gemoma_gso
        result_root=$project_root/results/copy_collapse/gemoma_glycine_v0.1
        ;;
    lifton)
        upstream_gff=$project_root/results/baselines/lifton_v1.0.11/glycine_v0.1/raw/upstream/lifton.gff3
        source_label=lifton_gso
        result_root=$project_root/results/copy_collapse/lifton_glycine_v0.1
        ;;
    *)
        echo "method must be gemoma or lifton" >&2
        exit 2
        ;;
esac
working_root=${result_root}.working

for required in "$python_bin" "$upstream_gff" "$blind_gff" \
                "$source_gff" "$truth"; do
    if [[ ! -s $required ]]; then
        echo "missing or empty Glycine copy-collapse input: $required" >&2
        exit 1
    fi
done
if [[ -e $result_root || -e $working_root ]]; then
    echo "refusing to overwrite Glycine copy-collapse result: $result_root" >&2
    exit 1
fi

mkdir -p "$working_root/blind" "$working_root/complete_control"
{
    printf 'field\tvalue\n'
    printf 'code_commit\t%s\n' "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'method\t%s\n' "$method"
    printf 'split\theldout\n'
    printf 'candidate_source\t%s\n' "$source_label"
    printf 'candidate_truth_access\tfalse\n'
    printf 'evaluator_truth_access\ttrue\n'
    printf 'pair_tsv_candidate_access\tfalse\n'
    printf 'upstream_projection_reused\ttrue\n'
    printf 'reference_species\tGlycine_soja\n'
    printf 'max_existing_cds_overlap\t0.2\n'
    printf 'max_redundancy_overlap\t0.5\n'
    printf 'paired_complete_annotation_control\ttrue\n'
    printf 'adapter_module_sha256\t%s\n' \
        "$(sha256sum "$code_root/src/ploidypatch/baseline.py" | awk '{print $1}')"
    printf 'score_module_sha256\t%s\n' \
        "$(sha256sum "$code_root/src/ploidypatch/score.py" | awk '{print $1}')"
} > "$working_root/run_contract.tsv"
{
    printf 'role\tbytes\tsha256\tpath\n'
    for entry in \
        "upstream_gff:$upstream_gff" \
        "blind_gff:$blind_gff" \
        "source_gff:$source_gff"; do
        role=${entry%%:*}
        path=${entry#*:}
        printf '%s\t%s\t%s\t%s\n' "$role" "$(stat -Lc %s "$path")" \
            "$(sha256sum "$path" | awk '{print $1}')" "$path"
    done
} > "$working_root/candidate_input_manifest.tsv"

cd "$code_root"
for mode in blind complete_control; do
    if [[ $mode == blind ]]; then
        annotation=$blind_gff
    else
        annotation=$source_gff
    fi
    /usr/bin/time -v -o "$working_root/$mode/resource.time.txt" \
        "$python_bin" -m ploidypatch.cli baseline adapt-gff \
            --perturbed-gff "$annotation" \
            --candidate-gff "$upstream_gff" \
            --source "$source_label" \
            --output-gff "$working_root/$mode/candidate.gff3" \
            --decisions-tsv "$working_root/$mode/decisions.tsv" \
            --max-existing-cds-overlap 0.2 \
            --max-redundancy-overlap 0.5 \
            > "$working_root/$mode/stdout.json" \
            2> "$working_root/$mode/stderr.log"
done

{
    printf 'role\tbytes\tsha256\tpath\n'
    for mode in blind complete_control; do
        for name in candidate.gff3 decisions.tsv; do
            path=$working_root/$mode/$name
            printf '%s_%s\t%s\t%s\t%s\n' "$mode" "$name" \
                "$(stat -Lc %s "$path")" \
                "$(sha256sum "$path" | awk '{print $1}')" "$path"
        done
    done
} > "$working_root/candidate_freeze.tsv"
{
    printf 'role\tbytes\tsha256\tpath\n'
    printf 'hidden_truth\t%s\t%s\t%s\n' "$(stat -Lc %s "$truth")" \
        "$(sha256sum "$truth" | awk '{print $1}')" "$truth"
} > "$working_root/evaluator_input_manifest.tsv"

/usr/bin/time -v -o "$working_root/score.resource.time.txt" \
    "$python_bin" -m ploidypatch.cli benchmark score \
        --source-gff "$source_gff" \
        --perturbed-gff "$blind_gff" \
        --candidate-gff "$working_root/blind/candidate.gff3" \
        --control-candidate-gff "$working_root/complete_control/candidate.gff3" \
        --truth "$truth" \
        --include-event-details \
        > "$working_root/score.json" \
        2> "$working_root/score.stderr.log"
if ! grep -q '"grade": "pass"' "$working_root/score.json"; then
    echo "Glycine copy-collapse $method score quality gate failed" >&2
    exit 1
fi

for output in "$working_root/blind/candidate.gff3" \
              "$working_root/blind/decisions.tsv" \
              "$working_root/complete_control/candidate.gff3" \
              "$working_root/complete_control/decisions.tsv" \
              "$working_root/candidate_freeze.tsv" \
              "$working_root/score.json"; do
    if [[ ! -s $output ]]; then
        echo "missing or empty Glycine copy-collapse output: $output" >&2
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
printf 'Glycine copy-collapse %s evaluation frozen: %s\n' "$method" "$result_root"
