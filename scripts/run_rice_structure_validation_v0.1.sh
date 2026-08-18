#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 2 || $# -gt 3 ]]; then
    echo "usage: $0 PROJECT_ROOT MIN_SOURCE_SUPPORT [cultivar_ungrouped|phylo_grouped]" >&2
    exit 2
fi

project_root=$(realpath "$1")
min_source_support=$2
reference_mode=${3:-cultivar_ungrouped}
if [[ $min_source_support != 1 && $min_source_support != 2 ]]; then
    echo "MIN_SOURCE_SUPPORT must be 1 or 2" >&2
    exit 2
fi

code_root=$project_root/code
python_bin=$project_root/envs/ploidypatch-dev/bin/python
source_gff=$project_root/data/derived/structure_sources_v0.1/osa_irgsp10/source.gff3
benchmark_root=$project_root/benchmark/structure/public_models_v0.1/osa_irgsp10
case $reference_mode in
    cultivar_ungrouped)
        projection_root=$project_root/results/projections/miniprot_v0.18/osa_irgsp10_multisource_v0.1/projection
        result_namespace=osa_irgsp10_v0.1
        reference_label=Oryza_sativa_IR64,Oryza_sativa_MH63,Oryza_sativa_ZS97
        group_args=()
        ;;
    phylo_grouped)
        projection_root=$project_root/results/projections/miniprot_v0.18/osa_irgsp10_phylo_v0.1/projection
        result_namespace=osa_irgsp10_phylo_grouped_v0.1
        reference_label=Oryza_rufipogon,Oryza_glaberrima,Oryza_brachyantha
        source_group_map=$code_root/config/source_groups/osa_phylo_v0.1.tsv
        group_args=(--source-group-map "$source_group_map")
        ;;
    *)
        echo "reference mode must be cultivar_ungrouped or phylo_grouped" >&2
        exit 2
        ;;
esac
miniprot_gff=$projection_root/miniprot.gff3
protein_map=$projection_root/references.protein_map.tsv
result_root=$project_root/results/heldout_structure/$result_namespace/source_support_${min_source_support}
working_root=${result_root}.working
events=(
    annotation_boundary_shift
    annotation_fused_gene
    annotation_missing_internal_exon
    annotation_split_gene
)

for required in "$python_bin" "$source_gff" "$miniprot_gff" "$protein_map"; do
    if [[ ! -s $required ]]; then
        echo "missing or empty rice validation input: $required" >&2
        exit 1
    fi
done
if [[ $reference_mode == phylo_grouped && ! -s $source_group_map ]]; then
    echo "missing or empty source-group map: $source_group_map" >&2
    exit 1
fi
for event in "${events[@]}"; do
    event_root=$benchmark_root/${event}_seed20260813
    for required in \
        "$event_root/blind/perturbed.gff3" \
        "$event_root/evaluator/truth/hidden_truth.json"; do
        if [[ ! -s $required ]]; then
            echo "missing or empty frozen rice benchmark input: $required" >&2
            exit 1
        fi
    done
done
if [[ -e $result_root || -e $working_root ]]; then
    echo "refusing to overwrite rice validation: $result_root" >&2
    exit 1
fi

mkdir -p "$working_root/complete_control"
{
    printf 'field\tvalue\n'
    printf 'dataset\tosa_irgsp10\n'
    printf 'split\theldout_no_retuning\n'
    printf 'code_commit\t%s\n' "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'reference_mode\t%s\n' "$reference_mode"
    printf 'references\t%s\n' "$reference_label"
    printf 'min_identity\t0.5\n'
    printf 'min_query_coverage\t0.5\n'
    printf 'min_source_support\t%s\n' "$min_source_support"
    printf 'min_gene_overlap_fraction\t0.1\n'
    printf 'require_intact\ttrue\n'
    printf 'topology_ambiguity_policy\temit_only_unique_exact_topology\n'
    printf 'parameters_frozen_from\tath_tair10_v0.1\n'
    if [[ $reference_mode == phylo_grouped ]]; then
        printf 'source_group_map_sha256\t%s\n' \
            "$(sha256sum "$source_group_map" | awk '{print $1}')"
    fi
    printf 'candidate_truth_access\tfalse\n'
    printf 'hypothesis_truth_access\tfalse\n'
    printf 'evaluator_truth_access\ttrue\n'
    for module in structure_candidate structure_hypothesis structure_hypothesis_score; do
        printf '%s_module_sha256\t%s\n' "$module" \
            "$(sha256sum "$code_root/src/ploidypatch/$module.py" | awk '{print $1}')"
    done
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
        event_root=$benchmark_root/${event}_seed20260813
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
/usr/bin/time -v -o "$working_root/complete_control/candidate.resource.time.txt" \
    "$python_bin" -m ploidypatch.cli baseline adapt-miniprot-structure \
        --annotation-gff "$source_gff" \
        --miniprot-gff "$miniprot_gff" \
        --protein-map "$protein_map" \
        "${group_args[@]}" \
        --min-source-support "$min_source_support" \
        --output-gff "$working_root/complete_control/candidate.gff3" \
        --decisions-tsv "$working_root/complete_control/decisions.tsv" \
        > "$working_root/complete_control/candidate.stdout.json" \
        2> "$working_root/complete_control/candidate.stderr.log"
/usr/bin/time -v -o "$working_root/complete_control/hypothesis.resource.time.txt" \
    "$python_bin" -m ploidypatch.cli baseline infer-structure-hypotheses \
        --annotation-gff "$source_gff" \
        --candidate-gff "$working_root/complete_control/candidate.gff3" \
        --output-tsv "$working_root/complete_control/hypotheses.tsv" \
        --candidate-topology-tsv "$working_root/complete_control/candidate_topology.tsv" \
        > "$working_root/complete_control/hypothesis.stdout.json" \
        2> "$working_root/complete_control/hypothesis.stderr.log"

for event in "${events[@]}"; do
    event_root=$benchmark_root/${event}_seed20260813
    blind_gff=$event_root/blind/perturbed.gff3
    truth=$event_root/evaluator/truth/hidden_truth.json
    event_output=$working_root/$event
    mkdir -p "$event_output"
    /usr/bin/time -v -o "$event_output/candidate.resource.time.txt" \
        "$python_bin" -m ploidypatch.cli baseline adapt-miniprot-structure \
            --annotation-gff "$blind_gff" \
            --miniprot-gff "$miniprot_gff" \
            --protein-map "$protein_map" \
            "${group_args[@]}" \
            --min-source-support "$min_source_support" \
            --output-gff "$event_output/candidate.gff3" \
            --decisions-tsv "$event_output/decisions.tsv" \
            > "$event_output/candidate.stdout.json" \
            2> "$event_output/candidate.stderr.log"
    /usr/bin/time -v -o "$event_output/hypothesis.resource.time.txt" \
        "$python_bin" -m ploidypatch.cli baseline infer-structure-hypotheses \
            --annotation-gff "$blind_gff" \
            --candidate-gff "$event_output/candidate.gff3" \
            --output-tsv "$event_output/hypotheses.tsv" \
            --candidate-topology-tsv "$event_output/candidate_topology.tsv" \
            > "$event_output/hypothesis.stdout.json" \
            2> "$event_output/hypothesis.stderr.log"
    /usr/bin/time -v -o "$event_output/score.resource.time.txt" \
        "$python_bin" -m ploidypatch.cli benchmark score-structure-hypotheses \
            --source-gff "$source_gff" \
            --perturbed-gff "$blind_gff" \
            --candidate-gff "$event_output/candidate.gff3" \
            --hypotheses-tsv "$event_output/hypotheses.tsv" \
            --control-hypotheses-tsv "$working_root/complete_control/hypotheses.tsv" \
            --truth "$truth" \
            --include-event-details \
            > "$event_output/score.json" \
            2> "$event_output/score.stderr.log"
    if ! grep -q '"grade": "pass"' "$event_output/score.json"; then
        echo "rice hypothesis score quality gate failed: $event" >&2
        exit 1
    fi
done

for output in \
    "$working_root/complete_control/candidate.gff3" \
    "$working_root/complete_control/decisions.tsv" \
    "$working_root/complete_control/hypotheses.tsv" \
    "$working_root/complete_control/candidate_topology.tsv"; do
    if [[ ! -s $output ]]; then
        echo "missing or empty rice validation output: $output" >&2
        exit 1
    fi
done
for event in "${events[@]}"; do
    for output in \
        "$working_root/$event/candidate.gff3" \
        "$working_root/$event/decisions.tsv" \
        "$working_root/$event/hypotheses.tsv" \
        "$working_root/$event/candidate_topology.tsv" \
        "$working_root/$event/score.json"; do
        if [[ ! -s $output ]]; then
            echo "missing or empty rice validation output: $output" >&2
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
printf 'heldout rice structure validation completed: %s\n' "$result_root"
