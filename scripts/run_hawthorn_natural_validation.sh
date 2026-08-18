#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "usage: $0 PROJECT_ROOT" >&2
    exit 2
fi

project_root=$(realpath "$1")
code_root=$project_root/code
python_bin=$project_root/envs/ploidypatch-dev/bin/python
data_root=$project_root/data/derived/hawthorn_black_projection_v0.1
projection_root=$project_root/results/hawthorn/haplotype_projection_v0.1
junction_root=$project_root/results/hawthorn/rna_junctions_v0.1
result_root=$project_root/results/hawthorn/natural_validation_v0.2
working_root=${result_root}.working
target_gff=$data_root/Black_Primary.gff3
miniprot_gff=$projection_root/raw/miniprot.gff3
protein_map=$projection_root/reference/haplotypes.map.tsv
junctions=$junction_root/primary_aggregate.tsv

for required in "$python_bin" "$target_gff" "$miniprot_gff" \
                "$protein_map" "$junctions" \
                "${junctions}.manifest.json"; do
    if [[ ! -s $required ]]; then
        echo "missing or empty natural-validation input: $required" >&2
        exit 1
    fi
done
if [[ -e $result_root || -e $working_root ]]; then
    echo "refusing to overwrite natural validation: $result_root" >&2
    exit 1
fi

mkdir -p "$working_root"
{
    printf 'field\tvalue\n'
    if git -C "$code_root" rev-parse HEAD >/dev/null 2>&1; then
        printf 'code_commit\t%s\n' "$(git -C "$code_root" rev-parse HEAD)"
    elif [[ -n "${PLOIDYPATCH_CODE_COMMIT:-}" ]]; then
        printf 'code_commit\t%s\n' "$PLOIDYPATCH_CODE_COMMIT"
    else
        printf 'code_commit\tunavailable_server_mirror\n'
    fi
    printf 'natural_module_sha256\t%s\n' \
        "$(sha256sum "$code_root/src/ploidypatch/natural.py" | awk '{print $1}')"
    printf 'cli_module_sha256\t%s\n' \
        "$(sha256sum "$code_root/src/ploidypatch/cli.py" | awk '{print $1}')"
    printf 'candidate_discovery_used_rna\tfalse\n'
    printf 'automatic_patch_policy\treview_required\n'
} > "$working_root/run_contract.tsv"

cd "$code_root"
/usr/bin/time -v -o "$working_root/discovery.resource.time.txt" \
    "$python_bin" -m ploidypatch.cli evidence discover-natural \
        --target-gff "$target_gff" \
        --miniprot-gff "$miniprot_gff" \
        --protein-map "$protein_map" \
        --output-tsv "$working_root/candidates.tsv" \
        --min-identity 0.8 \
        --min-query-coverage 0.8 \
        --max-existing-cds-overlap 0.1 \
        --min-boundary-extension-bp 30 \
        --near-best-score-fraction 0.95 \
        > "$working_root/discovery.stdout.log" \
        2> "$working_root/discovery.stderr.log"

sha256sum "$working_root/candidates.tsv" \
    "$working_root/candidates.tsv.manifest.json" \
    > "$working_root/PRE_RNA_VALIDATION_SHA256SUMS"
/usr/bin/time -v -o "$working_root/rna_validation.resource.time.txt" \
    "$python_bin" -m ploidypatch.cli evidence validate-natural-rna \
        --candidates "$working_root/candidates.tsv" \
        --junctions "$junctions" \
        --output-tsv "$working_root/candidates.primary_rna.tsv" \
        > "$working_root/rna_validation.stdout.log" \
        2> "$working_root/rna_validation.stderr.log"

for output in \
    "$working_root/candidates.tsv" \
    "$working_root/candidates.tsv.manifest.json" \
    "$working_root/candidates.primary_rna.tsv" \
    "$working_root/candidates.primary_rna.tsv.manifest.json"; do
    if [[ ! -s $output ]]; then
        echo "natural-validation output is missing or empty: $output" >&2
        exit 1
    fi
done
sha256sum "$working_root"/*.tsv "$working_root"/*.json \
    > "$working_root/POST_RNA_VALIDATION_SHA256SUMS"
du -sb "$working_root" > "$working_root/disk_bytes.txt"
mv "$working_root" "$result_root"
printf 'hawthorn natural validation completed: %s\n' "$result_root"
