#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 || $# -gt 2 ]]; then
    echo "usage: $0 PROJECT_ROOT [support2|union]" >&2
    exit 2
fi

project_root=$(realpath "$1")
candidate_arm=${2:-support2}
code_root=$project_root/code
python_bin=$project_root/envs/ploidypatch-dev/bin/python
benchmark_root=$project_root/benchmark/structure/copy_collapse_v0.1/gma_v21/annotation_copy_collapse_seed20260815
blind_gff=$benchmark_root/blind/perturbed.gff3
source_gff=$project_root/data/derived/syngap_glycine_v0.1/gma_complete/primary_chromosomes.gff3
truth=$benchmark_root/evaluator/truth/hidden_truth.json
case $candidate_arm in
    support2)
        result_family=wgd_reanchor_glycine_v0.1
        candidate_source=method_consensus_support2_plus_blind_self_wgd_reanchor
        ;;
    union)
        result_family=wgd_reanchor_union_glycine_dev_v0.1
        candidate_source=method_consensus_support1_union_plus_blind_self_wgd_reanchor
        ;;
    *)
        echo "candidate arm must be support2 or union" >&2
        exit 2
        ;;
esac
self_wgdi_root=$project_root/results/copy_collapse/$result_family/self_wgdi
blind_candidate=$self_wgdi_root/blind/selected/candidate.gff3
control_candidate=$self_wgdi_root/complete_control/selected/candidate.gff3
result_root=$project_root/results/copy_collapse/$result_family/evaluation
working_root=${result_root}.working

for mode in blind complete_control; do
    arm=$self_wgdi_root/$mode
    if [[ ! -s $arm/SHA256SUMS ]]; then
        echo "candidate self-WGD arm is not frozen: $arm" >&2
        exit 1
    fi
    (cd "$arm" && sha256sum -c SHA256SUMS >/dev/null)
done
for required in "$python_bin" "$blind_gff" "$source_gff" "$truth" \
                "$blind_candidate" "$control_candidate"; do
    if [[ ! -s $required ]]; then
        echo "missing or empty WGD reanchor evaluation input: $required" >&2
        exit 1
    fi
done
if [[ -e $result_root || -e $working_root ]]; then
    echo "refusing to overwrite WGD reanchor evaluation: $result_root" >&2
    exit 1
fi
mkdir -p "$working_root"

{
    printf 'field\tvalue\n'
    printf 'code_commit\t%s\n' "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'candidate_arm\t%s\n' "$candidate_arm"
    printf 'split\tpost_holdout_exploratory_development\n'
    printf 'formal_holdout_claim_allowed\tfalse\n'
    printf 'candidate_source\t%s\n' "$candidate_source"
    printf 'candidate_truth_access\tfalse\n'
    printf 'preperturbation_pair_access\tfalse\n'
    printf 'evaluator_truth_access\ttrue\n'
    printf 'paired_complete_annotation_control\ttrue\n'
    printf 'claim_boundary\tcandidate_named_wgd_syntenic_duplicate_not_gene_tree_homeology\n'
    printf 'score_module_sha256\t%s\n' \
        "$(sha256sum "$code_root/src/ploidypatch/score.py" | awk '{print $1}')"
} > "$working_root/run_contract.tsv"
{
    printf 'role\tbytes\tsha256\tpath\n'
    for entry in \
        "blind_candidate:$blind_candidate" \
        "control_candidate:$control_candidate" \
        "blind_gff:$blind_gff" \
        "source_gff:$source_gff"; do
        role=${entry%%:*}
        path=${entry#*:}
        printf '%s\t%s\t%s\t%s\n' "$role" "$(stat -Lc %s "$path")" \
            "$(sha256sum "$path" | awk '{print $1}')" "$path"
    done
} > "$working_root/candidate_freeze.tsv"
{
    printf 'role\tbytes\tsha256\tpath\n'
    printf 'hidden_truth\t%s\t%s\t%s\n' "$(stat -Lc %s "$truth")" \
        "$(sha256sum "$truth" | awk '{print $1}')" "$truth"
} > "$working_root/evaluator_input_manifest.tsv"

cd "$code_root"
/usr/bin/time -v -o "$working_root/score.resource.time.txt" \
    "$python_bin" -m ploidypatch.cli benchmark score \
        --source-gff "$source_gff" \
        --perturbed-gff "$blind_gff" \
        --candidate-gff "$blind_candidate" \
        --control-candidate-gff "$control_candidate" \
        --truth "$truth" \
        --include-event-details \
        > "$working_root/score.json" \
        2> "$working_root/score.stderr.log"
if ! grep -q '"grade": "pass"' "$working_root/score.json"; then
    echo "WGD reanchor score quality gate failed" >&2
    exit 1
fi
(
    cd "$working_root"
    find . -type f \( -name '*.json' -o -name '*.tsv' \) \
        -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS
)
du -sb "$working_root" > "$working_root/disk_bytes.txt"
mv "$working_root" "$result_root"
printf 'Glycine blind candidate self-WGD evaluation frozen: %s\n' "$result_root"
