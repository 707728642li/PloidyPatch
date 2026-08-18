#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "usage: $0 PROJECT_ROOT gemoma|lifton" >&2
    exit 2
fi

project_root=$(realpath "$1")
method=$2
derived=$project_root/data/derived/syngap_glycine_v0.1
target_dir=$derived/gma_blind
reference_dir=$derived/gso_reference
target_fasta=$target_dir/primary_chromosomes.genome.fa
reference_fasta=$reference_dir/primary_chromosomes.genome.fa
reference_gff=$reference_dir/primary_chromosomes.gff3

case $method in
    gemoma)
        bash "$project_root/code/scripts/run_gemoma_homology.sh" \
            "$project_root/envs/ploidypatch-gemoma" \
            "$project_root/results/baselines/gemoma_v1.9/glycine_v0.1/raw" \
            "$target_fasta" "$reference_fasta" "$reference_gff" Gso 32
        ;;
    lifton)
        bash "$project_root/code/scripts/run_lifton_transfer.sh" \
            "$project_root/envs/ploidypatch-lifton" \
            "$project_root/results/baselines/lifton_v1.0.11/glycine_v0.1/raw" \
            "$target_fasta" "$reference_fasta" "$reference_gff" Gso 32
        ;;
    *)
        echo "method must be gemoma or lifton" >&2
        exit 2
        ;;
esac
