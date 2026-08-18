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
projection_root=$project_root/results/baselines/miniprot_v0.18/ath_tair10_missing_gene_v0.1/projection
miniprot_gff=$projection_root/miniprot.gff3
protein_map=$projection_root/references.protein_map.tsv
result_root=$project_root/results/structure_candidates/miniprot_multisource/ath_tair10_v0.1/source_support_${min_source_support}
working_root=${result_root}.working
events=(
    annotation_boundary_shift
    annotation_fused_gene
    annotation_missing_internal_exon
    annotation_split_gene
)

for required in "$python_bin" "$source_gff" "$miniprot_gff" "$protein_map"; do
    if [[ ! -s $required ]]; then
        echo "missing or empty structure-candidate input: $required" >&2
        exit 1
    fi
done
for event in "${events[@]}"; do
    event_root=$benchmark_root/${event}_seed20260812
    for required in \
        "$event_root/blind/perturbed.gff3" \
        "$event_root/evaluator/truth/hidden_truth.json"; do
        if [[ ! -s $required ]]; then
            echo "missing or empty frozen benchmark input: $required" >&2
            exit 1
        fi
    done
done
if [[ -e $result_root || -e $working_root ]]; then
    echo "refusing to overwrite structure-candidate result: $result_root" >&2
    exit 1
fi

mkdir -p "$working_root/complete_control"
{
    printf 'field\tvalue\n'
    printf 'dataset\tath_tair10\n'
    printf 'split\tdevelopment\n'
    printf 'code_commit\t%s\n' "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'adapter\tadapt-miniprot-structure\n'
    printf 'min_identity\t0.5\n'
    printf 'min_query_coverage\t0.5\n'
    printf 'min_source_support\t%s\n' "$min_source_support"
    printf 'min_gene_overlap_fraction\t0.1\n'
    printf 'require_intact\ttrue\n'
    printf 'paired_complete_annotation_control\ttrue\n'
    printf 'candidate_truth_access\tfalse\n'
    printf 'evaluator_truth_access\ttrue\n'
    printf 'structure_candidate_module_sha256\t%s\n' \
        "$(sha256sum "$code_root/src/ploidypatch/structure_candidate.py" | awk '{print $1}')"
    printf 'score_module_sha256\t%s\n' \
        "$(sha256sum "$code_root/src/ploidypatch/score.py" | awk '{print $1}')"
} > "$working_root/run_contract.tsv"
{
    printf 'role\tbytes\tsha256\tpath\n'
    for entry in \
        "source_gff:$source_gff" \
        "miniprot_projection:$miniprot_gff" \
        "protein_map:$protein_map"; do
        role=${entry%%:*}
        path=${entry#*:}
        printf '%s\t%s\t%s\t%s\n' "$role" "$(stat -Lc %s "$path")" \
            "$(sha256sum "$path" | awk '{print $1}')" "$path"
    done
    for event in "${events[@]}"; do
        event_root=$benchmark_root/${event}_seed20260812
        for entry in \
            "${event}_blind:$event_root/blind/perturbed.gff3" \
            "${event}_truth:$event_root/evaluator/truth/hidden_truth.json"; do
            role=${entry%%:*}
            path=${entry#*:}
            printf '%s\t%s\t%s\t%s\n' "$role" "$(stat -Lc %s "$path")" \
                "$(sha256sum "$path" | awk '{print $1}')" "$path"
        done
    done
} > "$working_root/input_manifest.tsv"

cd "$code_root"
/usr/bin/time -v -o "$working_root/complete_control/resource.time.txt" \
    "$python_bin" -m ploidypatch.cli baseline adapt-miniprot-structure \
        --annotation-gff "$source_gff" \
        --miniprot-gff "$miniprot_gff" \
        --protein-map "$protein_map" \
        --min-source-support "$min_source_support" \
        --output-gff "$working_root/complete_control/candidate.gff3" \
        --decisions-tsv "$working_root/complete_control/decisions.tsv" \
        > "$working_root/complete_control/stdout.json" \
        2> "$working_root/complete_control/stderr.log"

for event in "${events[@]}"; do
    event_root=$benchmark_root/${event}_seed20260812
    blind_gff=$event_root/blind/perturbed.gff3
    truth=$event_root/evaluator/truth/hidden_truth.json
    event_output=$working_root/$event
    mkdir -p "$event_output"
    /usr/bin/time -v -o "$event_output/adapter.resource.time.txt" \
        "$python_bin" -m ploidypatch.cli baseline adapt-miniprot-structure \
            --annotation-gff "$blind_gff" \
            --miniprot-gff "$miniprot_gff" \
            --protein-map "$protein_map" \
            --min-source-support "$min_source_support" \
            --output-gff "$event_output/candidate.gff3" \
            --decisions-tsv "$event_output/decisions.tsv" \
            > "$event_output/adapter.stdout.json" \
            2> "$event_output/adapter.stderr.log"
    /usr/bin/time -v -o "$event_output/score.resource.time.txt" \
        "$python_bin" -m ploidypatch.cli benchmark score \
            --source-gff "$source_gff" \
            --perturbed-gff "$blind_gff" \
            --candidate-gff "$event_output/candidate.gff3" \
            --control-candidate-gff "$working_root/complete_control/candidate.gff3" \
            --truth "$truth" \
            --include-event-details \
            > "$event_output/score.json" \
            2> "$event_output/score.stderr.log"
done

for output in \
    "$working_root/complete_control/candidate.gff3" \
    "$working_root/complete_control/decisions.tsv"; do
    if [[ ! -s $output ]]; then
        echo "missing or empty structure-candidate output: $output" >&2
        exit 1
    fi
done
for event in "${events[@]}"; do
    for output in \
        "$working_root/$event/candidate.gff3" \
        "$working_root/$event/decisions.tsv" \
        "$working_root/$event/score.json"; do
        if [[ ! -s $output ]]; then
            echo "missing or empty structure-candidate output: $output" >&2
            exit 1
        fi
    done
done
(
    cd "$working_root"
    find . -type f \( -name '*.gff3' -o -name '*.json' -o -name '*.tsv' \) \
        -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS
)
du -sb "$working_root" > "$working_root/disk_bytes.txt"
mv "$working_root" "$result_root"
printf 'TAIR10 structure-candidate evaluation completed: %s\n' "$result_root"
