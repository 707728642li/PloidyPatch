#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "usage: $0 PROJECT_ROOT MIN_SOURCE_SUPPORT" >&2
    exit 2
fi

project_root=$(realpath "$1")
min_source_support=$2
if [[ $min_source_support != 1 && $min_source_support != 2 ]]; then
    echo "MIN_SOURCE_SUPPORT must be 1 or 2" >&2
    exit 2
fi

code_root=$project_root/code
python_bin=$project_root/envs/ploidypatch-dev/bin/python
source_gff=$project_root/data/derived/structure_sources_v0.1/ath_tair10/source.gff3
benchmark_root=$project_root/benchmark/structure/public_models_v0.1/ath_tair10
candidate_root=$project_root/results/structure_candidates/miniprot_multisource/ath_tair10_v0.1/source_support_${min_source_support}
result_root=$project_root/results/structure_hypotheses/exact_topology/ath_tair10_v0.1/source_support_${min_source_support}
working_root=${result_root}.working
events=(
    annotation_boundary_shift
    annotation_fused_gene
    annotation_missing_internal_exon
    annotation_split_gene
)

for required in \
    "$python_bin" \
    "$source_gff" \
    "$candidate_root/complete_control/candidate.gff3"; do
    if [[ ! -s $required ]]; then
        echo "missing or empty hypothesis input: $required" >&2
        exit 1
    fi
done
for event in "${events[@]}"; do
    event_root=$benchmark_root/${event}_seed20260812
    for required in \
        "$event_root/blind/perturbed.gff3" \
        "$event_root/evaluator/truth/hidden_truth.json" \
        "$candidate_root/$event/candidate.gff3"; do
        if [[ ! -s $required ]]; then
            echo "missing or empty hypothesis input: $required" >&2
            exit 1
        fi
    done
done
if [[ -e $result_root || -e $working_root ]]; then
    echo "refusing to overwrite hypothesis result: $result_root" >&2
    exit 1
fi

mkdir -p "$working_root/complete_control"
{
    printf 'field\tvalue\n'
    printf 'dataset\tath_tair10\n'
    printf 'split\tdevelopment\n'
    printf 'code_commit\t%s\n' "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'candidate_min_source_support\t%s\n' "$min_source_support"
    printf 'topology_ambiguity_policy\temit_only_unique_exact_topology\n'
    printf 'candidate_truth_access\tfalse\n'
    printf 'hypothesis_truth_access\tfalse\n'
    printf 'hypothesis_module_sha256\t%s\n' \
        "$(sha256sum "$code_root/src/ploidypatch/structure_hypothesis.py" | awk '{print $1}')"
} > "$working_root/run_contract.tsv"
{
    printf 'role\tbytes\tsha256\tpath\n'
    for entry in \
        "source_gff:$source_gff" \
        "control_candidate:$candidate_root/complete_control/candidate.gff3"; do
        role=${entry%%:*}
        path=${entry#*:}
        printf '%s\t%s\t%s\t%s\n' "$role" "$(stat -Lc %s "$path")" \
            "$(sha256sum "$path" | awk '{print $1}')" "$path"
    done
    for event in "${events[@]}"; do
        event_root=$benchmark_root/${event}_seed20260812
        for entry in \
            "${event}_blind:$event_root/blind/perturbed.gff3" \
            "${event}_truth:$event_root/evaluator/truth/hidden_truth.json" \
            "${event}_candidate:$candidate_root/$event/candidate.gff3"; do
            role=${entry%%:*}
            path=${entry#*:}
            printf '%s\t%s\t%s\t%s\n' "$role" "$(stat -Lc %s "$path")" \
                "$(sha256sum "$path" | awk '{print $1}')" "$path"
        done
    done
} > "$working_root/input_manifest.tsv"

cd "$code_root"
/usr/bin/time -v -o "$working_root/complete_control/resource.time.txt" \
    "$python_bin" -m ploidypatch.cli baseline infer-structure-hypotheses \
        --annotation-gff "$source_gff" \
        --candidate-gff "$candidate_root/complete_control/candidate.gff3" \
        --output-tsv "$working_root/complete_control/hypotheses.tsv" \
        --candidate-topology-tsv "$working_root/complete_control/candidate_topology.tsv" \
        > "$working_root/complete_control/stdout.json" \
        2> "$working_root/complete_control/stderr.log"

for event in "${events[@]}"; do
    event_root=$benchmark_root/${event}_seed20260812
    event_output=$working_root/$event
    mkdir -p "$event_output"
    /usr/bin/time -v -o "$event_output/resource.time.txt" \
        "$python_bin" -m ploidypatch.cli baseline infer-structure-hypotheses \
            --annotation-gff "$event_root/blind/perturbed.gff3" \
            --candidate-gff "$candidate_root/$event/candidate.gff3" \
            --output-tsv "$event_output/hypotheses.tsv" \
            --candidate-topology-tsv "$event_output/candidate_topology.tsv" \
            > "$event_output/stdout.json" \
            2> "$event_output/stderr.log"
done

for output in \
    "$working_root/complete_control/hypotheses.tsv" \
    "$working_root/complete_control/candidate_topology.tsv"; do
    if [[ ! -s $output ]]; then
        echo "missing or empty hypothesis output: $output" >&2
        exit 1
    fi
done
for event in "${events[@]}"; do
    for output in \
        "$working_root/$event/hypotheses.tsv" \
        "$working_root/$event/candidate_topology.tsv"; do
        if [[ ! -s $output ]]; then
            echo "missing or empty hypothesis output: $output" >&2
            exit 1
        fi
    done
done
(
    cd "$working_root"
    find . -type f \( -name '*.json' -o -name '*.tsv' \) \
        -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS
)
du -sb "$working_root" > "$working_root/disk_bytes.txt"
mv "$working_root" "$result_root"
printf 'TAIR10 exact-topology hypotheses completed: %s\n' "$result_root"
