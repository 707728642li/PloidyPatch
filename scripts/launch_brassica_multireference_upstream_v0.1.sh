#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "usage: $0 PROJECT_ROOT" >&2
    exit 2
fi

project_root=$(realpath "$1")
code_root=$project_root/code
input_root=$project_root/data/derived/multireference_development/brassica_v0.1
target=$project_root/data/derived/normalized_bundles/v0.1/bna_daae_primary/primary_chromosomes.genome.fa
result_root=$project_root/results/baselines/multireference_brassica_v0.1
state_root=$project_root/logs/baseline/multireference_brassica_v0.1

if [[ ! -s $input_root/SHA256SUMS ]]; then
    echo "Brassica multireference inputs are not frozen" >&2
    exit 1
fi
(cd "$input_root" && sha256sum -c SHA256SUMS >/dev/null)
for required in "$target" "${target}.fai"; do
    if [[ ! -s $required ]]; then
        echo "missing Brassica target input: $required" >&2
        exit 1
    fi
done
if [[ -e $result_root || -e $state_root ]]; then
    echo "refusing to reuse Brassica multireference upstream state" >&2
    exit 1
fi
mkdir -p "$result_root" "$state_root"

launch_one() {
    local method=$1
    local reference=$2
    local environment runner reference_id
    case $method in
        gemoma)
            environment=$project_root/envs/ploidypatch-gemoma
            runner=$code_root/scripts/run_gemoma_homology.sh
            ;;
        lifton)
            environment=$project_root/envs/ploidypatch-lifton
            runner=$code_root/scripts/run_lifton_transfer.sh
            ;;
        *) return 2 ;;
    esac
    case $reference in
        brapa) reference_id=Brassica_rapa_A ;;
        bol) reference_id=Brassica_oleracea_C ;;
        *) return 2 ;;
    esac
    local reference_root=$input_root/$reference
    local run_root=$result_root/$method/$reference
    local log=$state_root/$method.$reference.launcher.log
    local pid_file=$state_root/$method.$reference.pid
    mkdir -p "$(dirname "$run_root")"
    nohup setsid bash "$runner" \
        "$environment" "$run_root" "$target" \
        "$reference_root/primary_chromosomes.genome.fa" \
        "$reference_root/primary_chromosomes.gff3" \
        "$reference_id" 32 > "$log" 2>&1 &
    local pid=$!
    printf '%s\n' "$pid" > "$pid_file"
    printf '%s\t%s\t%s\t%s\n' "$method" "$reference" "$pid" "$log"
}

{
    printf 'method\treference\tpid\tlog\n'
    launch_one gemoma brapa
    sleep 10
    launch_one gemoma bol
    sleep 10
    launch_one lifton brapa
    sleep 10
    launch_one lifton bol
} > "$state_root/launch_manifest.tsv"
sleep 5
while IFS=$'\t' read -r method reference pid log; do
    [[ $method == method ]] && continue
    if ! kill -0 "$pid" 2>/dev/null; then
        echo "$method/$reference exited during startup; inspect $log" >&2
        exit 1
    fi
done < "$state_root/launch_manifest.tsv"
printf 'Brassica multireference upstream jobs launched: %s\n' "$state_root/launch_manifest.tsv"
