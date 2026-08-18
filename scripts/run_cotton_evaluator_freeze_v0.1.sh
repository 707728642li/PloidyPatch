#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then echo "usage: $0 PROJECT_ROOT" >&2; exit 2; fi
project_root=$(realpath "$1")
code_root=$project_root/code
commit=${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}
state=$project_root/logs/copy_collapse/cotton_evaluator_freeze_v0.1
if [[ -e $state ]]; then echo "cotton evaluator state exists" >&2; exit 1; fi
mkdir -p "$state"
record_failure() {
    rc=$?
    printf 'pipeline_failed\t%s\texit_%s\n' "$(date --iso-8601=seconds)" "$rc" >> "$state/status.tsv"
    exit "$rc"
}
trap record_failure ERR
printf 'started\t%s\ncode_commit\t%s\n' "$(date --iso-8601=seconds)" "$commit" > "$state/status.tsv"
declare -A artifact=(
    [prepare_cotton_holdout_wgdi_inputs_v0.1]="$project_root/data/derived/holdout_evaluator/cotton_wgdi_inputs_v0.1"
    [run_cotton_holdout_homeolog_wgdi_v0.1]="$project_root/results/evaluator/cotton_holdout_v0.1/wgdi"
    [infer_cotton_holdout_homeolog_pairs_v0.1]="$project_root/results/evaluator/cotton_holdout_v0.1/homeolog_pairs"
    [run_cotton_copy_collapse_holdout_v0.1]="$project_root/benchmark/structure/copy_collapse_v0.1/ghi_ad/annotation_copy_collapse_seed20260817"
)
for stage in prepare_cotton_holdout_wgdi_inputs_v0.1 \
             run_cotton_holdout_homeolog_wgdi_v0.1 \
             infer_cotton_holdout_homeolog_pairs_v0.1 \
             run_cotton_copy_collapse_holdout_v0.1; do
    root=${artifact[$stage]}
    if [[ -s $root/SHA256SUMS ]]; then
        (cd "$root" && sha256sum -c SHA256SUMS >/dev/null)
        printf '%s_reused\t%s\n' "$stage" "$(date --iso-8601=seconds)" >> "$state/status.tsv"
    else
        PLOIDYPATCH_CODE_COMMIT="$commit" bash "$code_root/scripts/$stage.sh" "$project_root" \
            > "$state/$stage.log" 2>&1
        printf '%s_complete\t%s\n' "$stage" "$(date --iso-8601=seconds)" >> "$state/status.tsv"
    fi
done
printf 'pipeline_complete\t%s\n' "$(date --iso-8601=seconds)" >> "$state/status.tsv"
