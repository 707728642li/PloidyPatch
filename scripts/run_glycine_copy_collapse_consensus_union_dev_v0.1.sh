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
miniprot_root=$project_root/results/copy_collapse/miniprot_glycine_v0.1
gemoma_root=$project_root/results/copy_collapse/gemoma_glycine_v0.1
lifton_root=$project_root/results/copy_collapse/lifton_glycine_v0.1
result_root=$project_root/results/copy_collapse/consensus_union_glycine_dev_v0.1
working_root=${result_root}.working

for required in "$python_bin" "$blind_gff" "$source_gff" "$truth" \
                "$miniprot_root/blind/candidate.gff3" \
                "$miniprot_root/complete_control/candidate.gff3" \
                "$gemoma_root/blind/candidate.gff3" \
                "$gemoma_root/complete_control/candidate.gff3" \
                "$lifton_root/blind/candidate.gff3" \
                "$lifton_root/complete_control/candidate.gff3"; do
    if [[ ! -s $required ]]; then
        echo "missing or empty consensus-union input: $required" >&2
        exit 1
    fi
done
if [[ -e $result_root || -e $working_root ]]; then
    echo "refusing to overwrite Glycine consensus-union result: $result_root" >&2
    exit 1
fi
mkdir -p "$working_root/blind" "$working_root/complete_control"

{
    printf 'field\tvalue\n'
    printf 'code_commit\t%s\n' "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'split\tpost_holdout_exploratory_development\n'
    printf 'formal_holdout_claim_allowed\tfalse\n'
    printf 'reason\tsoybean_truth_was_already_inspected_before_this_arm\n'
    printf 'candidate_truth_access\tfalse\n'
    printf 'pair_tsv_candidate_access\tfalse\n'
    printf 'method_families\tminiprot,gemoma,lifton\n'
    printf 'min_method_support\t1\n'
    printf 'consensus_unit\texact_seqid_strand_phased_cds_chain\n'
    printf 'max_redundancy_overlap\t0.5\n'
    printf 'output_scope\tcoding_model_only\n'
    printf 'consensus_module_sha256\t%s\n' \
        "$(sha256sum "$code_root/src/ploidypatch/consensus.py" | awk '{print $1}')"
} > "$working_root/run_contract.tsv"
{
    printf 'role\tbytes\tsha256\tpath\n'
    for entry in \
        "blind_gff:$blind_gff" \
        "source_gff:$source_gff" \
        "miniprot_blind:$miniprot_root/blind/candidate.gff3" \
        "miniprot_control:$miniprot_root/complete_control/candidate.gff3" \
        "gemoma_blind:$gemoma_root/blind/candidate.gff3" \
        "gemoma_control:$gemoma_root/complete_control/candidate.gff3" \
        "lifton_blind:$lifton_root/blind/candidate.gff3" \
        "lifton_control:$lifton_root/complete_control/candidate.gff3"; do
        role=${entry%%:*}
        path=${entry#*:}
        printf '%s\t%s\t%s\t%s\n' "$role" "$(stat -Lc %s "$path")" \
            "$(sha256sum "$path" | awk '{print $1}')" "$path"
    done
} > "$working_root/candidate_input_manifest.tsv"

cd "$code_root"
for mode in blind complete_control; do
    if [[ $mode == blind ]]; then
        base=$blind_gff
    else
        base=$source_gff
    fi
    /usr/bin/time -v -o "$working_root/$mode/resource.time.txt" \
        "$python_bin" -m ploidypatch.cli baseline select-method-consensus \
            --base-gff "$base" \
            --candidate "miniprot=$miniprot_root/$mode/candidate.gff3" \
            --candidate "gemoma=$gemoma_root/$mode/candidate.gff3" \
            --candidate "lifton=$lifton_root/$mode/candidate.gff3" \
            --output-gff "$working_root/$mode/candidate.gff3" \
            --decisions-tsv "$working_root/$mode/decisions.tsv" \
            --min-method-support 1 \
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
    echo "Glycine consensus-union score quality gate failed" >&2
    exit 1
fi
(
    cd "$working_root"
    find . -type f \( -name '*.gff3' -o -name '*.json' -o -name '*.tsv' \) \
        -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS
)
du -sb "$working_root" > "$working_root/disk_bytes.txt"
mv "$working_root" "$result_root"
printf 'Glycine post-holdout consensus union frozen: %s\n' "$result_root"
