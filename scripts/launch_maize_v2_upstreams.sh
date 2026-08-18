#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then echo "usage: $0 PROJECT_ROOT" >&2; exit 2; fi
project_root=$(realpath "$1")
code_root=$project_root/code
input_root=$project_root/data/derived/holdout_inputs/maize_v2
target=$input_root/zea_mays/primary_chromosomes.genome.fa
benchmark=$project_root/benchmark/structure/copy_collapse_v0.2/zma_maize1/annotation_copy_collapse_seed20260829
result_root=$project_root/results/baselines/maize_v2
state=$project_root/logs/baseline/maize_v2

[[ -s $benchmark/SHA256SUMS ]] || { echo "maize holdout is not frozen" >&2; exit 1; }
expected_blind_sha=$(awk '$2 == "./blind/perturbed.gff3" {print $1}' "$benchmark/SHA256SUMS")
observed_blind_sha=$(sha256sum "$benchmark/blind/perturbed.gff3" | awk '{print $1}')
[[ -n $expected_blind_sha && $observed_blind_sha == "$expected_blind_sha" ]] || {
    echo "maize blind target checksum failed" >&2; exit 1;
}
if [[ -e $result_root/gemoma || -e $result_root/lifton \
      || -e $result_root/miniprot || -e $state ]]; then
    echo "refusing to reuse maize upstream state" >&2; exit 1
fi
mkdir -p "$result_root/gemoma" "$result_root/lifton" "$state"
{
    printf 'field\tvalue\ncode_commit\t%s\n' \
        "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'candidate_truth_access\tfalse\ntarget_complete_annotation_access\tfalse\n'
    printf 'method_family_reference_vote_policy\tone_vote_after_within_method_merge\n'
    printf 'reference_species\tSorghum_bicolor,Setaria_italica\n'
    printf 'threads_per_transfer\t24\nminiprot_threads\t32\n'
} > "$state/run_contract.tsv"

launch_transfer() {
    local method=$1 ref=$2 reference_id=$3
    local env runner
    case $method in
        gemoma) env=$project_root/envs/ploidypatch-gemoma; runner=$code_root/scripts/run_gemoma_homology.sh ;;
        lifton) env=$project_root/envs/ploidypatch-lifton; runner=$code_root/scripts/run_lifton_transfer.sh ;;
        *) return 2 ;;
    esac
    local run_root=$result_root/$method/$ref
    nohup setsid bash "$runner" "$env" "$run_root" "$target" \
        "$input_root/$ref/primary_chromosomes.genome.fa" \
        "$input_root/$ref/primary_chromosomes.gff3" "$reference_id" 24 \
        > "$state/$method.$ref.log" 2>&1 &
    printf '%s\n' "$!" > "$state/$method.$ref.pid"
}
launch_transfer gemoma sorghum_bicolor Sorghum_bicolor_NCBIv3
sleep 3
launch_transfer gemoma setaria_italica Setaria_italica_v2
sleep 3
launch_transfer lifton sorghum_bicolor Sorghum_bicolor_NCBIv3
sleep 3
launch_transfer lifton setaria_italica Setaria_italica_v2
sleep 3
nohup setsid env PLOIDYPATCH_CODE_COMMIT="${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}" \
    bash "$code_root/scripts/run_maize_v2_miniprot_upstream.sh" "$project_root" \
    > "$state/miniprot.log" 2>&1 &
printf '%s\n' "$!" > "$state/miniprot.pid"
sleep 3
for key in gemoma.sorghum_bicolor gemoma.setaria_italica \
           lifton.sorghum_bicolor lifton.setaria_italica miniprot; do
    pid=$(cat "$state/$key.pid")
    if ! kill -0 "$pid" 2>/dev/null; then
        echo "maize upstream exited during startup: $key" >&2; exit 1
    fi
    printf '%s\t%s\n' "$key" "$pid"
done
