#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then echo "usage: $0 PROJECT_ROOT CODE_COMMIT" >&2; exit 2; fi
project_root=$(realpath "$1"); commit=$2; code_root=$project_root/code
upstream=$project_root/results/baselines/cotton_holdout_v0.1
launch_state=$project_root/logs/baseline/cotton_holdout_v0.1
state=$project_root/logs/copy_collapse/cotton_candidate_pipeline_v0.1
result=$project_root/results/copy_collapse/holdout/cotton_method_trio_v0.1
if [[ -e $state ]]; then echo "cotton candidate watcher state exists" >&2; exit 1; fi
mkdir -p "$state"
record_failure() { rc=$?; printf 'pipeline_failed\t%s\texit_%s\n' "$(date --iso-8601=seconds)" "$rc" >> "$state/status.tsv"; exit "$rc"; }
trap record_failure ERR
printf 'started\t%s\ncode_commit\t%s\n' "$(date --iso-8601=seconds)" "$commit" > "$state/status.tsv"
declare -A output=(
    [gemoma.gar_a]="$upstream/gemoma/gar_a/upstream/final_annotation.gff"
    [gemoma.gra_d]="$upstream/gemoma/gra_d/upstream/final_annotation.gff"
    [lifton.gar_a]="$upstream/lifton/gar_a/upstream/lifton.gff3"
    [lifton.gra_d]="$upstream/lifton/gra_d/upstream/lifton.gff3"
    [miniprot]="$upstream/miniprot/SHA256SUMS"
)
while true; do
    remaining=0
    for key in gemoma.gar_a gemoma.gra_d lifton.gar_a lifton.gra_d miniprot; do
        [[ -s ${output[$key]} ]] && continue
        remaining=$((remaining + 1)); pid=$(cat "$launch_state/$key.pid")
        if ! kill -0 "$pid" 2>/dev/null; then echo "upstream failed: $key" >&2; exit 1; fi
    done
    (( remaining > 0 )) || break
    sleep 30
done
printf 'upstreams_complete\t%s\n' "$(date --iso-8601=seconds)" >> "$state/status.tsv"
PLOIDYPATCH_CODE_COMMIT="$commit" bash "$code_root/scripts/run_cotton_holdout_method_trio_v0.1.sh" "$project_root" \
    > "$state/method_trio.log" 2>&1
printf 'method_trio_complete\t%s\nresult\t%s\n' "$(date --iso-8601=seconds)" "$result" >> "$state/status.tsv"
