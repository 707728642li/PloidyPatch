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
validation_root=$project_root/results/hawthorn/natural_validation_v0.2
context_root=$project_root/results/hawthorn/natural_assembly_context_v0.1
context_working=${context_root}.working
graph_root=$project_root/results/hawthorn/natural_event_graph_v0.3
graph_working=${graph_root}.working
validation=$validation_root/candidates.primary_rna.tsv
genome=$data_root/Black_Primary.fa
fai=$data_root/Black_Primary.fa.fai

for required in "$python_bin" "$validation" "${validation}.manifest.json" \
                "$genome" "$fai"; do
    if [[ ! -s $required ]]; then
        echo "missing assembly-context prerequisite: $required" >&2
        exit 1
    fi
done

if [[ ! -e $context_root ]]; then
    if [[ -e $context_working ]]; then
        echo "refusing to overwrite assembly-context work: $context_working" >&2
        exit 1
    fi
    mkdir -p "$context_working"
    {
        printf 'field\tvalue\n'
        printf 'code_commit\t%s\n' "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
        printf 'natural_module_sha256\t%s\n' \
            "$(sha256sum "$code_root/src/ploidypatch/natural.py" | awk '{print $1}')"
        printf 'event_graph_module_sha256\t%s\n' \
            "$(sha256sum "$code_root/src/ploidypatch/event_graph.py" | awk '{print $1}')"
        printf 'flank_bp\t5000\n'
        printf 'max_ambiguous_fraction\t0\n'
        printf 'softmask_policy\tcontext_only_uncalibrated\n'
        printf 'automatic_patch_policy\treview_required\n'
    } > "$context_working/run_contract.tsv"
    cd "$code_root"
    /usr/bin/time -v -o "$context_working/resource.time.txt" \
        "$python_bin" -m ploidypatch.cli evidence annotate-natural-assembly \
            --validation-tsv "$validation" \
            --genome "$genome" \
            --fai "$fai" \
            --output-tsv "$context_working/candidates.primary_rna.assembly.tsv" \
            --flank-bp 5000 \
            --max-ambiguous-fraction 0 \
            > "$context_working/stdout.log" \
            2> "$context_working/stderr.log"
    for output in \
        "$context_working/candidates.primary_rna.assembly.tsv" \
        "$context_working/candidates.primary_rna.assembly.tsv.manifest.json"; do
        if [[ ! -s $output ]]; then
            echo "assembly-context output is missing or empty: $output" >&2
            exit 1
        fi
    done
    (
        cd "$context_working"
        sha256sum ./*.tsv ./*.json > SHA256SUMS
    )
    du -sb "$context_working" > "$context_working/disk_bytes.txt"
    mv "$context_working" "$context_root"
fi

context=$context_root/candidates.primary_rna.assembly.tsv
if [[ ! -s $context || ! -s ${context}.manifest.json ]]; then
    echo "published assembly context is incomplete: $context_root" >&2
    exit 1
fi
if [[ -e $graph_root || -e $graph_working ]]; then
    echo "refusing to overwrite natural graph v0.3" >&2
    exit 1
fi
mkdir -p "$graph_working"
cp "$context_root/run_contract.tsv" "$graph_working/run_contract.tsv"
cd "$code_root"
"$python_bin" -m ploidypatch.cli evidence prepare-natural-graph \
    --validation-tsv "$context" \
    --output-candidates "$graph_working/candidates.tsv" \
    --output-evidence "$graph_working/evidence.tsv" \
    > "$graph_working/adapter.stdout.log" \
    2> "$graph_working/adapter.stderr.log"
/usr/bin/time -v -o "$graph_working/resource.time.txt" \
    "$python_bin" -m ploidypatch.cli graph infer \
        --candidates "$graph_working/candidates.tsv" \
        --evidence "$graph_working/evidence.tsv" \
        --output-json "$graph_working/event_graph.json" \
        --decisions-tsv "$graph_working/decisions.tsv" \
        > "$graph_working/infer.stdout.log" \
        2> "$graph_working/infer.stderr.log"
for output in \
    "$graph_working/candidates.tsv" \
    "$graph_working/candidates.tsv.manifest.json" \
    "$graph_working/evidence.tsv" \
    "$graph_working/event_graph.json" \
    "$graph_working/decisions.tsv"; do
    if [[ ! -s $output ]]; then
        echo "natural graph v0.3 output is missing or empty: $output" >&2
        exit 1
    fi
done
(
    cd "$graph_working"
    sha256sum ./*.tsv ./*.json > SHA256SUMS
)
du -sb "$graph_working" > "$graph_working/disk_bytes.txt"
mv "$graph_working" "$graph_root"
printf 'hawthorn assembly context and graph completed: %s\n' "$graph_root"
