#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "usage: $0 ISOLATED_PROJECT_ROOT" >&2
    exit 2
fi
project_root=$(realpath "$1")
code_root=$project_root/code
baseline_root=$project_root/results/baselines/populus_v0.4
normalized=$baseline_root/miniprot/normalized
target=$normalized/target/primary_chromosomes.genome.fa

# A second check inside the namespace makes accidental mount broadening fail
# before an upstream program starts.
for forbidden in \
    "$project_root/data/derived/external_inputs/populus_v0.4/evaluator_only" \
    "$project_root/benchmark/structure/copy_collapse_v0.4/ptr_v4/annotation_copy_collapse_seed20260930/evaluator" \
    "$project_root/truth" "$project_root/labels" /nas_data; do
    [[ ! -e $forbidden ]] || {
        echo "forbidden evaluator/NAS path visible in blind namespace: $forbidden" >&2
        exit 1
    }
done

bash "$code_root/scripts/run_populus_miniprot_upstream_v0.4.sh" "$project_root"
[[ -s $target && -s ${target}.fai ]] || {
    echo "candidate-safe target FASTA/index missing after normalization" >&2
    exit 1
}

pids=()
labels=()
for reference in salix_purpurea salix_suchowensis; do
    reference_fasta=$normalized/$reference/primary_chromosomes.genome.fa
    reference_gff=$normalized/$reference/primary_chromosomes.gff3
    lifton_reference_gff=$normalized/$reference/primary_chromosomes.lifton.gff3
    for input in "$reference_fasta" "$reference_gff" "$lifton_reference_gff"; do
        [[ -s $input ]] || { echo "missing normalized candidate reference: $input" >&2; exit 1; }
    done
    (
        bash "$code_root/scripts/run_gemoma_homology.sh" \
            "$project_root/envs/ploidypatch-gemoma" \
            "$baseline_root/gemoma/$reference" \
            "$target" "$reference_fasta" "$reference_gff" "$reference" 32
    ) &
    pids+=("$!")
    labels+=("gemoma:$reference")
    (
        bash "$code_root/scripts/run_lifton_transfer.sh" \
            "$project_root/envs/ploidypatch-lifton" \
            "$baseline_root/lifton/$reference" \
            "$target" "$reference_fasta" "$lifton_reference_gff" "$reference" 16
    ) &
    pids+=("$!")
    labels+=("lifton:$reference")
done
failed=0
for index in "${!pids[@]}"; do
    if ! wait "${pids[$index]}"; then
        printf 'blind upstream failed: %s\n' "${labels[$index]}" >&2
        failed=1
    fi
done
[[ $failed -eq 0 ]] || exit 1

bash "$code_root/scripts/build_populus_method_trio_candidate_pools_v0.4.sh" "$project_root"
bash "$code_root/scripts/run_populus_blind_union_self_wgd_v0.4.sh" "$project_root"
bash "$code_root/scripts/score_populus_candidates_blind_v0.4.sh" "$project_root"

printf 'Populus v0.4 blind pipeline completed inside isolated namespace\n'
