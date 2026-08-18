#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then echo "usage: $0 PROJECT_ROOT" >&2; exit 2; fi
project_root=$(realpath "$1")
code_root=$project_root/code
input_root=$project_root/data/derived/holdout_inputs/cotton_v0.1
target=$input_root/hirsutum/primary_chromosomes.genome.fa
benchmark=$project_root/benchmark/structure/copy_collapse_v0.1/ghi_ad/annotation_copy_collapse_seed20260817
result_root=$project_root/results/baselines/cotton_holdout_v0.1
state=$project_root/logs/baseline/cotton_holdout_v0.1
if [[ ! -s $benchmark/SHA256SUMS ]]; then echo "cotton holdout is not frozen" >&2; exit 1; fi
(cd "$benchmark" && sha256sum -c SHA256SUMS >/dev/null)
if [[ -e $result_root/gemoma || -e $result_root/lifton || -e $result_root/miniprot || -e $state ]]; then
    echo "refusing to reuse cotton upstream state" >&2; exit 1
fi
mkdir -p "$result_root/gemoma" "$result_root/lifton" "$state"
{
    printf 'field\tvalue\ncode_commit\t%s\n' "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'candidate_truth_access\tfalse\ntarget_complete_annotation_access\tfalse\n'
    printf 'method_family_reference_vote_policy\tone_vote_after_within_method_merge\n'
} > "$state/run_contract.tsv"

launch_transfer() {
    local method=$1 ref=$2 species=$3 reference_id=$4
    local env runner
    case $method in
        gemoma) env=$project_root/envs/ploidypatch-gemoma; runner=$code_root/scripts/run_gemoma_homology.sh ;;
        lifton) env=$project_root/envs/ploidypatch-lifton; runner=$code_root/scripts/run_lifton_transfer.sh ;;
        *) return 2 ;;
    esac
    local run_root=$result_root/$method/$ref
    nohup setsid bash "$runner" "$env" "$run_root" "$target" \
        "$input_root/$species/primary_chromosomes.genome.fa" \
        "$input_root/$species/primary_chromosomes.gff3" "$reference_id" 32 \
        > "$state/$method.$ref.log" 2>&1 &
    printf '%s\n' "$!" > "$state/$method.$ref.pid"
}
launch_transfer gemoma gar_a arboreum Gossypium_arboreum_A
sleep 5
launch_transfer gemoma gra_d raimondii Gossypium_raimondii_D
sleep 5
launch_transfer lifton gar_a arboreum Gossypium_arboreum_A
sleep 5
launch_transfer lifton gra_d raimondii Gossypium_raimondii_D
sleep 5
nohup setsid env PLOIDYPATCH_CODE_COMMIT="${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}" \
    bash "$code_root/scripts/run_cotton_miniprot_upstream_v0.1.sh" "$project_root" \
    > "$state/miniprot.log" 2>&1 &
printf '%s\n' "$!" > "$state/miniprot.pid"
sleep 3
for key in gemoma.gar_a gemoma.gra_d lifton.gar_a lifton.gra_d miniprot; do
    pid=$(cat "$state/$key.pid")
    if ! kill -0 "$pid" 2>/dev/null; then echo "cotton upstream exited during startup: $key" >&2; exit 1; fi
    printf '%s\t%s\n' "$key" "$pid"
done
