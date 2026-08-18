#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "usage: $0 PROJECT_ROOT" >&2
    exit 2
fi

project_root=$(realpath "$1")
code_root=$project_root/code
python_bin=$project_root/envs/ploidypatch-dev/bin/python
benchmark_root=$project_root/benchmark/heldout/v0.1/gma_v21_annotation_missing_gene_seed20260809
blind_gff=$benchmark_root/blind/perturbed.gff3
source_gff=$benchmark_root/evaluator/restoration/restored.gff3
truth=$benchmark_root/evaluator/truth/hidden_truth.json
strata=$benchmark_root/evaluator/catalog/gma_v21_primary.annotation_missing_gene.tsv
baseline_root=$project_root/results/baselines/syngap_v1.2.5/glycine_v0.1/genblastg
result_root=$baseline_root/evaluation_phase_normalized_v0.2
working_root=${result_root}.working

resolve_full_gff() {
    local mode=$1
    local manifest=$baseline_root/$mode/output_manifest.tsv
    if [[ ! -s $manifest ]]; then
        echo "validated SynGAP output manifest is absent: $manifest" >&2
        return 1
    fi
    local paths=()
    mapfile -t paths < <(
        awk -F '\t' '$1 == "Gma.SynGAP.gff3" { print $4 }' "$manifest"
    )
    if [[ ${#paths[@]} -ne 1 || ! -s ${paths[0]} ]]; then
        echo "expected one full Gma SynGAP GFF for $mode" >&2
        return 1
    fi
    local expected
    expected=$(awk -F '\t' '$1 == "Gma.SynGAP.gff3" { print $3 }' "$manifest")
    if [[ $(sha256sum "${paths[0]}" | awk '{print $1}') != "$expected" ]]; then
        echo "SynGAP full GFF checksum mismatch for $mode" >&2
        return 1
    fi
    printf '%s\n' "${paths[0]}"
}

blind_upstream=$(resolve_full_gff blind)
control_upstream=$(resolve_full_gff complete_control)
for required in "$python_bin" "$blind_upstream" "$control_upstream" \
                "$blind_gff" "$source_gff" "$truth" "$strata"; do
    if [[ ! -s $required ]]; then
        echo "missing or empty SynGAP evaluation input: $required" >&2
        exit 1
    fi
done
if [[ -e $result_root || -e $working_root ]]; then
    echo "refusing to overwrite SynGAP evaluation: $result_root" >&2
    exit 1
fi
mkdir -p "$working_root/blind" "$working_root/complete_control"
{
    printf 'field\tvalue\n'
    printf 'code_commit\t%s\n' "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'method\tSynGAP-1.2.5-genblastg\n'
    printf 'candidate_source\tsyngap_gso\n'
    printf 'max_existing_cds_overlap\t0.2\n'
    printf 'max_redundancy_overlap\t0.5\n'
    printf 'missing_cds_phase_policy\tinfer_only_when_all_missing_full_cds_first_phase_zero\n'
    printf 'paired_complete_annotation_control\ttrue\n'
    printf 'upstream_completion_manifest_required\ttrue\n'
} > "$working_root/run_contract.tsv"
{
    printf 'role\tbytes\tsha256\tpath\n'
    for entry in \
        "blind_upstream:$blind_upstream" \
        "control_upstream:$control_upstream" \
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

cd "$code_root"
for mode in blind complete_control; do
    if [[ $mode == blind ]]; then
        perturbed=$blind_gff
        upstream=$blind_upstream
    else
        perturbed=$source_gff
        upstream=$control_upstream
    fi
    /usr/bin/time -v -o "$working_root/$mode/resource.time.txt" \
        "$python_bin" -m ploidypatch.cli baseline adapt-gff \
            --perturbed-gff "$perturbed" \
            --candidate-gff "$upstream" \
            --source syngap_gso \
            --output-gff "$working_root/$mode/candidate.gff3" \
            --decisions-tsv "$working_root/$mode/decisions.tsv" \
            --max-existing-cds-overlap 0.2 \
            --max-redundancy-overlap 0.5 \
            --infer-missing-cds-phase \
            > "$working_root/$mode/stdout.log" \
            2> "$working_root/$mode/stderr.log"
done

/usr/bin/time -v -o "$working_root/score.resource.time.txt" \
    "$python_bin" -m ploidypatch.cli benchmark score \
        --source-gff "$source_gff" \
        --perturbed-gff "$blind_gff" \
        --candidate-gff "$working_root/blind/candidate.gff3" \
        --control-candidate-gff "$working_root/complete_control/candidate.gff3" \
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
        echo "SynGAP evaluation output is missing or empty: $output" >&2
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
printf 'SynGAP phase-normalized matched evaluation completed: %s\n' "$result_root"
