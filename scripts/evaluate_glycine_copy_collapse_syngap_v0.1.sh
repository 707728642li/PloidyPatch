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
blind_gff=$benchmark_root/blind/perturbed.gff3
source_gff=$project_root/data/derived/syngap_glycine_v0.1/gma_complete/primary_chromosomes.gff3
truth=$benchmark_root/evaluator/truth/hidden_truth.json
copy_root=$project_root/results/copy_collapse/syngap_glycine_v0.1
blind_root=$copy_root/upstream/blind
control_root=$project_root/results/baselines/syngap_v1.2.5/glycine_v0.1/genblastg/complete_control
result_root=$copy_root/evaluation
working_root=${result_root}.working

resolve_full_gff() {
    local run_root=$1
    local expected_input=$2
    local manifest=$run_root/output_manifest.tsv
    local contract=$run_root/run_contract.tsv
    if [[ ! -s $manifest || ! -s $contract ]]; then
        echo "validated SynGAP manifests are absent: $run_root" >&2
        return 1
    fi
    local configured_input
    configured_input=$(awk -F '\t' '$1 == "sp1_gff" { print $2 }' "$contract")
    if [[ ! -s $configured_input ]] ||
       [[ $(sha256sum "$configured_input" | awk '{print $1}') != \
          $(sha256sum "$expected_input" | awk '{print $1}') ]]; then
        echo "SynGAP target annotation does not match expected input: $run_root" >&2
        return 1
    fi
    local paths=()
    mapfile -t paths < <(
        awk -F '\t' '$1 == "Gma.SynGAP.gff3" { print $4 }' "$manifest"
    )
    if [[ ${#paths[@]} -ne 1 || ! -s ${paths[0]} ]]; then
        echo "expected one full Gma SynGAP GFF: $run_root" >&2
        return 1
    fi
    local expected_sha
    expected_sha=$(awk -F '\t' '$1 == "Gma.SynGAP.gff3" { print $3 }' "$manifest")
    if [[ $(sha256sum "${paths[0]}" | awk '{print $1}') != "$expected_sha" ]]; then
        echo "SynGAP full GFF checksum mismatch: $run_root" >&2
        return 1
    fi
    printf '%s\n' "${paths[0]}"
}

blind_upstream=$(resolve_full_gff "$blind_root" "$blind_gff")
control_upstream=$(resolve_full_gff "$control_root" "$source_gff")
for required in "$python_bin" "$blind_upstream" "$control_upstream" \
                "$blind_gff" "$source_gff" "$truth"; do
    if [[ ! -s $required ]]; then
        echo "missing or empty SynGAP copy-collapse input: $required" >&2
        exit 1
    fi
done
if [[ -e $result_root || -e $working_root ]]; then
    echo "refusing to overwrite SynGAP copy-collapse evaluation: $result_root" >&2
    exit 1
fi

mkdir -p "$working_root/blind" "$working_root/complete_control"
{
    printf 'field\tvalue\n'
    printf 'code_commit\t%s\n' "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'method\tSynGAP-1.2.5-genblastg\n'
    printf 'split\theldout\n'
    printf 'candidate_source\tsyngap_gso\n'
    printf 'candidate_truth_access\tfalse\n'
    printf 'evaluator_truth_access\ttrue\n'
    printf 'pair_tsv_candidate_access\tfalse\n'
    printf 'reference_species\tGlycine_soja\n'
    printf 'blind_upstream_rerun_for_copy_collapse\ttrue\n'
    printf 'complete_control_upstream_reused\ttrue\n'
    printf 'max_existing_cds_overlap\t0.2\n'
    printf 'max_redundancy_overlap\t0.5\n'
    printf 'missing_cds_phase_policy\tinfer_only_when_all_missing_full_cds_first_phase_zero\n'
    printf 'paired_complete_annotation_control\ttrue\n'
    printf 'adapter_module_sha256\t%s\n' \
        "$(sha256sum "$code_root/src/ploidypatch/baseline.py" | awk '{print $1}')"
    printf 'score_module_sha256\t%s\n' \
        "$(sha256sum "$code_root/src/ploidypatch/score.py" | awk '{print $1}')"
} > "$working_root/run_contract.tsv"
{
    printf 'role\tbytes\tsha256\tpath\n'
    for entry in \
        "blind_upstream:$blind_upstream" \
        "control_upstream:$control_upstream" \
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
        upstream=$blind_upstream
    else
        annotation=$source_gff
        upstream=$control_upstream
    fi
    /usr/bin/time -v -o "$working_root/$mode/resource.time.txt" \
        "$python_bin" -m ploidypatch.cli baseline adapt-gff \
            --perturbed-gff "$annotation" \
            --candidate-gff "$upstream" \
            --source syngap_gso \
            --output-gff "$working_root/$mode/candidate.gff3" \
            --decisions-tsv "$working_root/$mode/decisions.tsv" \
            --max-existing-cds-overlap 0.2 \
            --max-redundancy-overlap 0.5 \
            --infer-missing-cds-phase \
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
    echo "SynGAP copy-collapse score quality gate failed" >&2
    exit 1
fi

for output in "$working_root/blind/candidate.gff3" \
              "$working_root/blind/decisions.tsv" \
              "$working_root/complete_control/candidate.gff3" \
              "$working_root/complete_control/decisions.tsv" \
              "$working_root/candidate_freeze.tsv" \
              "$working_root/score.json"; do
    if [[ ! -s $output ]]; then
        echo "missing or empty SynGAP copy-collapse output: $output" >&2
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
printf 'SynGAP Glycine copy-collapse evaluation frozen: %s\n' "$result_root"
