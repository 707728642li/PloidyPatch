#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "usage: $0 PROJECT_ROOT" >&2
    exit 2
fi

project_root=$(realpath "$1")
code_root=$project_root/code
python_bin=$project_root/envs/ploidypatch-dev/bin/python
primary_validation=$project_root/results/hawthorn/natural_validation_v0.2/candidates.primary_rna.tsv
grouped_junctions=$project_root/results/hawthorn/rna_junctions_grouped_v0.1/all_samples.grouped_aggregate.tsv
result_root=$project_root/results/hawthorn/natural_secondary_rna_v0.1
working_root=${result_root}.working

for required in "$python_bin" "$primary_validation" \
                "${primary_validation}.manifest.json" "$grouped_junctions" \
                "${grouped_junctions}.manifest.json"; do
    if [[ ! -s $required ]]; then
        echo "missing secondary-validation prerequisite: $required" >&2
        exit 1
    fi
done
if [[ -e $result_root || -e $working_root ]]; then
    echo "refusing to overwrite secondary natural validation" >&2
    exit 1
fi
mkdir -p "$working_root"
{
    printf 'field\tvalue\n'
    printf 'code_commit\t%s\n' "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'candidate_discovery_used_secondary_rna\tfalse\n'
    printf 'secondary_group_support_is_context_not_truth\ttrue\n'
    printf 'negative_evidence_policy\tabsence_is_missing_not_contradiction\n'
    printf 'automatic_patch_policy\treview_required\n'
} > "$working_root/run_contract.tsv"

cd "$code_root"
/usr/bin/time -v -o "$working_root/resource.time.txt" \
    "$python_bin" -m ploidypatch.cli evidence validate-natural-secondary-rna \
        --primary-validation "$primary_validation" \
        --grouped-junctions "$grouped_junctions" \
        --output-tsv "$working_root/candidates.secondary_group_rna.tsv" \
        > "$working_root/stdout.log" \
        2> "$working_root/stderr.log"
for output in \
    "$working_root/candidates.secondary_group_rna.tsv" \
    "$working_root/candidates.secondary_group_rna.tsv.manifest.json"; do
    if [[ ! -s $output ]]; then
        echo "secondary natural validation is missing or empty: $output" >&2
        exit 1
    fi
done
(
    cd "$working_root"
    sha256sum ./*.tsv ./*.json > SHA256SUMS
)
du -sb "$working_root" > "$working_root/disk_bytes.txt"
mv "$working_root" "$result_root"
printf 'hawthorn secondary natural validation completed: %s\n' "$result_root"
