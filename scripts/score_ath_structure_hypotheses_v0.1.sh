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
hypothesis_root=$project_root/results/structure_hypotheses/exact_topology/ath_tair10_v0.1/source_support_${min_source_support}
result_root=$project_root/results/structure_hypotheses/evaluation/ath_tair10_v0.1/source_support_${min_source_support}
working_root=${result_root}.working
control_hypotheses=$hypothesis_root/complete_control/hypotheses.tsv
events=(
    annotation_boundary_shift
    annotation_fused_gene
    annotation_missing_internal_exon
    annotation_split_gene
)

for required in "$python_bin" "$source_gff" "$control_hypotheses"; do
    if [[ ! -s $required ]]; then
        echo "missing or empty hypothesis-scoring input: $required" >&2
        exit 1
    fi
done
for event in "${events[@]}"; do
    event_root=$benchmark_root/${event}_seed20260812
    for required in \
        "$event_root/blind/perturbed.gff3" \
        "$event_root/evaluator/truth/hidden_truth.json" \
        "$candidate_root/$event/candidate.gff3" \
        "$hypothesis_root/$event/hypotheses.tsv"; do
        if [[ ! -s $required ]]; then
            echo "missing or empty hypothesis-scoring input: $required" >&2
            exit 1
        fi
    done
done
if [[ -e $result_root || -e $working_root ]]; then
    echo "refusing to overwrite hypothesis evaluation: $result_root" >&2
    exit 1
fi

mkdir -p "$working_root"
{
    printf 'field\tvalue\n'
    printf 'dataset\tath_tair10\n'
    printf 'split\tdevelopment\n'
    printf 'code_commit\t%s\n' "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'candidate_min_source_support\t%s\n' "$min_source_support"
    printf 'evaluation_mode\tpaired_complete_annotation_difference\n'
    printf 'score_module_sha256\t%s\n' \
        "$(sha256sum "$code_root/src/ploidypatch/structure_hypothesis_score.py" | awk '{print $1}')"
} > "$working_root/run_contract.tsv"
{
    printf 'role\tbytes\tsha256\tpath\n'
    for entry in \
        "source_gff:$source_gff" \
        "control_hypotheses:$control_hypotheses"; do
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
            "${event}_candidate:$candidate_root/$event/candidate.gff3" \
            "${event}_hypotheses:$hypothesis_root/$event/hypotheses.tsv"; do
            role=${entry%%:*}
            path=${entry#*:}
            printf '%s\t%s\t%s\t%s\n' "$role" "$(stat -Lc %s "$path")" \
                "$(sha256sum "$path" | awk '{print $1}')" "$path"
        done
    done
} > "$working_root/input_manifest.tsv"

cd "$code_root"
for event in "${events[@]}"; do
    event_root=$benchmark_root/${event}_seed20260812
    event_output=$working_root/$event
    mkdir -p "$event_output"
    /usr/bin/time -v -o "$event_output/resource.time.txt" \
        "$python_bin" -m ploidypatch.cli benchmark score-structure-hypotheses \
            --source-gff "$source_gff" \
            --perturbed-gff "$event_root/blind/perturbed.gff3" \
            --candidate-gff "$candidate_root/$event/candidate.gff3" \
            --hypotheses-tsv "$hypothesis_root/$event/hypotheses.tsv" \
            --control-hypotheses-tsv "$control_hypotheses" \
            --truth "$event_root/evaluator/truth/hidden_truth.json" \
            --include-event-details \
            > "$event_output/score.json" \
            2> "$event_output/stderr.log"
    if ! grep -q '"grade": "pass"' "$event_output/score.json"; then
        echo "hypothesis score quality gate failed: $event" >&2
        exit 1
    fi
done

for event in "${events[@]}"; do
    if [[ ! -s $working_root/$event/score.json ]]; then
        echo "missing or empty hypothesis score: $event" >&2
        exit 1
    fi
done
(
    cd "$working_root"
    find . -type f \( -name '*.json' -o -name '*.tsv' \) \
        -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS
)
du -sb "$working_root" > "$working_root/disk_bytes.txt"
mv "$working_root" "$result_root"
printf 'TAIR10 structure-hypothesis scoring completed: %s\n' "$result_root"
