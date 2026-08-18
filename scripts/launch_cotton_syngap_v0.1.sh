#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then echo "usage: $0 PROJECT_ROOT" >&2; exit 2; fi
project_root=$(realpath "$1")
code_root=$project_root/code
syngap_env=$project_root/envs/ploidypatch-syngap
input_root=$project_root/data/derived/holdout_inputs/cotton_v0.1
prepared=$project_root/data/derived/candidate_inputs/cotton_syngap_v0.1
result_root=$project_root/results/baselines/cotton_holdout_v0.1/syngap
state_root=$project_root/logs/baseline/cotton_holdout_syngap_v0.1
target_genome=$input_root/hirsutum/primary_chromosomes.genome.fa
declare -A reference_genome=(
    [gar_a]="$input_root/arboreum/primary_chromosomes.genome.fa"
    [gra_d]="$input_root/raimondii/primary_chromosomes.genome.fa"
)
declare -A reference_gff=(
    [gar_a]="$project_root/data/derived/candidate_inputs/cotton_holdout_v0.1/gar_a_lifton_compat_v3/reference.compat.gff3"
    [gra_d]="$project_root/data/derived/candidate_inputs/cotton_holdout_v0.1/gra_d_lifton_compat_v2/reference.compat.gff3"
)
declare -A reference_label=([gar_a]=Gar [gra_d]=Gra)

if [[ ! -s $prepared/SHA256SUMS ]]; then echo "cotton SynGAP inputs are not frozen" >&2; exit 1; fi
(cd "$prepared" && sha256sum -c SHA256SUMS >/dev/null)
if [[ -e $result_root || -e $state_root ]]; then
    echo "refusing to reuse cotton SynGAP run state" >&2; exit 1
fi
mkdir -p "$result_root" "$state_root"
{
    printf 'field\tvalue\n'
    printf 'code_commit\t%s\n' "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'process\tgenblastg\nthreads_per_arm\t32\narms\t4\n'
    printf 'candidate_truth_access\tfalse\n'
    printf 'method_role\texternal_comparator_not_primary_policy\n'
} > "$state_root/run_contract.tsv"
for mode in blind complete_control; do
    target_gff=$prepared/target/$mode/annotation.compat.gff3
    for ref in gar_a gra_d; do
        run_root=$result_root/$ref/$mode
        mkdir -p "$(dirname "$run_root")"
        nohup setsid bash "$code_root/scripts/run_syngap_dual.sh" \
            "$syngap_env" "$run_root" "$target_genome" "$target_gff" \
            "${reference_genome[$ref]}" "${reference_gff[$ref]}" \
            Ghi "${reference_label[$ref]}" 32 genblastg \
            > "$state_root/$ref.$mode.launcher.log" 2>&1 &
        printf '%s\n' "$!" > "$state_root/$ref.$mode.pid"
        sleep 15
    done
done
sleep 3
for mode in blind complete_control; do
    for ref in gar_a gra_d; do
        pid=$(cat "$state_root/$ref.$mode.pid")
        if ! kill -0 "$pid" 2>/dev/null; then
            echo "cotton SynGAP arm exited during startup: $ref $mode" >&2; exit 1
        fi
        printf '%s_%s_pid\t%s\n' "$ref" "$mode" "$pid"
    done
done
