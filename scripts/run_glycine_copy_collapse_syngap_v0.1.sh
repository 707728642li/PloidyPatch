#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "usage: $0 PROJECT_ROOT" >&2
    exit 2
fi

project_root=$(realpath "$1")
derived=$project_root/data/derived/syngap_glycine_v0.1
benchmark_root=$project_root/benchmark/structure/copy_collapse_v0.1/gma_v21/annotation_copy_collapse_seed20260815
reference_dir=$derived/gso_reference
target_dir=$derived/gma_complete
output_root=$project_root/results/copy_collapse/syngap_glycine_v0.1/upstream/blind

bash "$project_root/code/scripts/run_syngap_dual.sh" \
    "$project_root/envs/ploidypatch-syngap" \
    "$output_root" \
    "$target_dir/primary_chromosomes.genome.fa" \
    "$benchmark_root/blind/perturbed.gff3" \
    "$reference_dir/primary_chromosomes.genome.fa" \
    "$reference_dir/primary_chromosomes.gff3" \
    Gma Gso 32 genblastg
