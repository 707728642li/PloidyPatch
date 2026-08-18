#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then echo "usage: $0 PROJECT_ROOT CODE_COMMIT" >&2; exit 2; fi
project_root=$1
commit=$2
code_root=$project_root/code
upstream_state=$project_root/logs/baseline/multireference_brassica_v0.1
pipeline_state=$project_root/logs/copy_collapse/brassica_accelerated_pipeline_v0.1
upstream_root=$project_root/results/baselines/multireference_brassica_v0.1
method_root=$project_root/results/copy_collapse/model_development/brassica_method_trio_v0.1
self_root=$project_root/results/copy_collapse/model_development/brassica_union_self_wgd_v0.1
transfer_root=$project_root/results/copy_collapse/model_development/brassica_glycine_model_transfer_v0.1

record_failure() {
    rc=$?
    printf 'pipeline_failed\t%s\texit_%s\n' "$(date --iso-8601=seconds)" "$rc" \
        >> "$pipeline_state/status.tsv"
    exit "$rc"
}
trap record_failure ERR

printf 'watch_started\t%s\ncode_commit\t%s\n' "$(date --iso-8601=seconds)" "$commit" > "$pipeline_state/status.tsv"
declare -A expected=(
    [gemoma.brapa]="$upstream_root/gemoma/brapa/upstream/final_annotation.gff"
    [gemoma.bol]="$upstream_root/gemoma/bol/upstream/final_annotation.gff"
    [lifton.brapa]="$upstream_root/lifton/brapa/upstream/lifton.gff3"
    [lifton.bol]="$upstream_root/lifton/bol/upstream/lifton.gff3"
)
while true; do
    remaining=0
    for key in gemoma.brapa gemoma.bol lifton.brapa lifton.bol; do
        if [[ -s ${expected[$key]} ]]; then continue; fi
        remaining=$((remaining + 1))
        pid=$(cat "$upstream_state/$key.pid")
        if ! kill -0 "$pid" 2>/dev/null; then
            printf 'upstream_failed\t%s\t%s\n' "$key" "$(date --iso-8601=seconds)" >> "$pipeline_state/status.tsv"
            exit 1
        fi
    done
    [[ $remaining -gt 0 ]] || break
    sleep 30
done
printf 'upstream_complete\t%s\n' "$(date --iso-8601=seconds)" >> "$pipeline_state/status.tsv"
if [[ -s $method_root/SHA256SUMS ]]; then
    (cd "$method_root" && sha256sum -c SHA256SUMS >/dev/null)
    printf 'method_trio_reused\t%s\n' "$(date --iso-8601=seconds)" >> "$pipeline_state/status.tsv"
else
    PLOIDYPATCH_CODE_COMMIT="$commit" bash "$code_root/scripts/run_brassica_multireference_method_trio_v0.1.sh" "$project_root"
    printf 'method_trio_complete\t%s\n' "$(date --iso-8601=seconds)" >> "$pipeline_state/status.tsv"
fi
if [[ -s $self_root/blind/SHA256SUMS && -s $self_root/complete_control/SHA256SUMS ]]; then
    for mode in blind complete_control; do
        (cd "$self_root/$mode" && sha256sum -c SHA256SUMS >/dev/null)
    done
    printf 'self_wgd_reused\t%s\n' "$(date --iso-8601=seconds)" >> "$pipeline_state/status.tsv"
else
    PLOIDYPATCH_CODE_COMMIT="$commit" bash "$code_root/scripts/launch_brassica_union_self_wgd_v0.1.sh" "$project_root"
fi
for mode in blind complete_control; do
    while [[ ! -s $self_root/$mode/SHA256SUMS ]]; do
        pid=$(cat "$project_root/logs/copy_collapse/brassica_union_self_wgd_v0.1/$mode.pid")
        if ! kill -0 "$pid" 2>/dev/null; then
            printf 'self_wgd_failed\t%s\t%s\n' "$mode" "$(date --iso-8601=seconds)" >> "$pipeline_state/status.tsv"
            exit 1
        fi
        sleep 30
    done
done
printf 'self_wgd_complete\t%s\n' "$(date --iso-8601=seconds)" >> "$pipeline_state/status.tsv"
if [[ -s $transfer_root/SHA256SUMS ]]; then
    (cd "$transfer_root" && sha256sum -c SHA256SUMS >/dev/null)
    printf 'transfer_reused\t%s\n' "$(date --iso-8601=seconds)" >> "$pipeline_state/status.tsv"
elif [[ -d ${transfer_root}.working ]]; then
    PLOIDYPATCH_CODE_COMMIT="$commit" PLOIDYPATCH_RESUME_STAGE=evaluator \
        bash "$code_root/scripts/evaluate_brassica_glycine_copy_model_transfer_v0.1.sh" "$project_root"
    printf 'transfer_resumed\t%s\n' "$(date --iso-8601=seconds)" >> "$pipeline_state/status.tsv"
else
    PLOIDYPATCH_CODE_COMMIT="$commit" \
        bash "$code_root/scripts/evaluate_brassica_glycine_copy_model_transfer_v0.1.sh" "$project_root"
    printf 'transfer_complete\t%s\n' "$(date --iso-8601=seconds)" >> "$pipeline_state/status.tsv"
fi
printf 'pipeline_complete\t%s\nresult\t%s\n' "$(date --iso-8601=seconds)" "$transfer_root" >> "$pipeline_state/status.tsv"
