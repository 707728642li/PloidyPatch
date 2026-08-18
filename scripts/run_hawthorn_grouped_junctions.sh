#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "usage: $0 PROJECT_ROOT" >&2
    exit 2
fi

project_root=$(realpath "$1")
code_root=$project_root/code
python_bin=$project_root/envs/ploidypatch-dev/bin/python
primary_root=$project_root/results/hawthorn/rna_junctions_v0.1
secondary_root=$project_root/results/hawthorn/rna_junctions_secondary_v0.1
groups=$code_root/config/hawthorn_rna_sample_groups_v0.1.tsv
result_root=$project_root/results/hawthorn/rna_junctions_grouped_v0.1
working_root=${result_root}.working

for required in "$python_bin" "$primary_root/primary_aggregate.tsv" \
                "$secondary_root/all_samples_aggregate.tsv" "$groups"; do
    if [[ ! -s $required ]]; then
        echo "missing grouped-junction prerequisite: $required" >&2
        exit 1
    fi
done
if [[ -e $result_root || -e $working_root ]]; then
    echo "refusing to overwrite grouped junction result" >&2
    exit 1
fi
mkdir -p "$working_root"
{
    printf 'field\tvalue\n'
    printf 'code_commit\t%s\n' "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'sample_group_interpretation\tfilename_stem_only_not_biological_metadata\n'
    printf 'min_reads_per_sample\t2\n'
    printf 'min_samples_per_group\t2\n'
    printf 'min_secondary_groups\t2\n'
    printf 'negative_evidence_policy\tabsence_is_missing_not_contradiction\n'
} > "$working_root/run_contract.tsv"

cd "$code_root"
/usr/bin/time -v -o "$working_root/resource.time.txt" \
    "$python_bin" -m ploidypatch.cli evidence aggregate-junctions \
        --input-dir "$primary_root/samples" \
        --input-dir "$secondary_root/samples" \
        --primary-sample Black-1 \
        --primary-sample Black-2 \
        --primary-sample Black-3 \
        --sample-groups "$groups" \
        --min-reads-per-sample 2 \
        --min-supporting-samples 2 \
        --min-samples-per-group 2 \
        --min-secondary-groups 2 \
        --output-tsv "$working_root/all_samples.grouped_aggregate.tsv" \
        > "$working_root/stdout.log" \
        2> "$working_root/stderr.log"

for output in \
    "$working_root/all_samples.grouped_aggregate.tsv" \
    "$working_root/all_samples.grouped_aggregate.tsv.manifest.json"; do
    if [[ ! -s $output ]]; then
        echo "grouped-junction output is missing or empty: $output" >&2
        exit 1
    fi
done
(
    cd "$working_root"
    sha256sum ./*.tsv ./*.json > SHA256SUMS
)
du -sb "$working_root" > "$working_root/disk_bytes.txt"
mv "$working_root" "$result_root"
printf 'hawthorn grouped junction evidence completed: %s\n' "$result_root"
