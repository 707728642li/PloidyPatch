#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "usage: $0 PROJECT_ROOT blind|complete_control" >&2
    exit 2
fi

project_root=$(realpath "$1")
mode=$2
derived=$project_root/data/derived/syngap_glycine_v0.1
reference_dir=$derived/gso_reference
output_root=$project_root/results/baselines/syngap_v1.2.5/glycine_v0.1/genblastg

case $mode in
    blind)
        target_dir=$derived/gma_blind
        ;;
    complete_control)
        target_dir=$derived/gma_complete
        ;;
    *)
        echo "mode must be blind or complete_control" >&2
        exit 2
        ;;
esac

bash "$project_root/code/scripts/run_syngap_dual.sh" \
    "$project_root/envs/ploidypatch-syngap" \
    "$output_root/$mode" \
    "$target_dir/primary_chromosomes.genome.fa" \
    "$target_dir/primary_chromosomes.gff3" \
    "$reference_dir/primary_chromosomes.genome.fa" \
    "$reference_dir/primary_chromosomes.gff3" \
    Gma Gso 32 genblastg
