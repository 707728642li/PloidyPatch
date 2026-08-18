#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then echo "usage: $0 PROJECT_ROOT" >&2; exit 2; fi
project_root=$(realpath "$1")
code_root=$project_root/code
upstream=$project_root/results/baselines/maize_v2
state=$project_root/logs/baseline/maize_v2
status=$state/post_pipeline.tsv
commit=${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}
[[ -d $state ]] || { echo "missing maize upstream state" >&2; exit 1; }
[[ ! -e $status ]] || { echo "maize post-pipeline status exists" >&2; exit 1; }
printf 'watcher_started\t%s\ncode_commit\t%s\n' \
    "$(date --iso-8601=seconds)" "$commit" > "$status"

expected=(
    "$upstream/gemoma/sorghum_bicolor/upstream/final_annotation.gff"
    "$upstream/gemoma/setaria_italica/upstream/final_annotation.gff"
    "$upstream/lifton/sorghum_bicolor/upstream/lifton.gff3"
    "$upstream/lifton/setaria_italica/upstream/lifton.gff3"
    "$upstream/miniprot/raw/miniprot.gff3"
    "$upstream/miniprot/reference/maize_outgroups.map.tsv"
)
keys=(
    gemoma.sorghum_bicolor gemoma.setaria_italica
    lifton.sorghum_bicolor lifton.setaria_italica miniprot
)
while true; do
    ready=true
    for path in "${expected[@]}"; do [[ -s $path ]] || ready=false; done
    if [[ $ready == true ]]; then break; fi
    active=0
    for key in "${keys[@]}"; do
        pid=$(cat "$state/$key.pid")
        if kill -0 "$pid" 2>/dev/null; then active=$((active + 1)); fi
    done
    if [[ $active -eq 0 ]]; then
        printf 'upstream_failed\t%s\tmissing_expected_output\n' \
            "$(date --iso-8601=seconds)" >> "$status"
        exit 1
    fi
    sleep 30
done
printf 'upstreams_complete\t%s\n' "$(date --iso-8601=seconds)" >> "$status"

PLOIDYPATCH_CODE_COMMIT="$commit" \
    bash "$code_root/scripts/build_maize_v2_method_trio_blind.sh" "$project_root" \
    > "$state/build_maize_v2_method_trio_blind.log" 2>&1
printf 'blind_candidates_frozen\t%s\n' "$(date --iso-8601=seconds)" >> "$status"
PLOIDYPATCH_CODE_COMMIT="$commit" \
    bash "$code_root/scripts/evaluate_maize_v2_method_trio.sh" "$project_root" \
    > "$state/evaluate_maize_v2_method_trio.log" 2>&1
printf 'method_evaluation_complete\t%s\n' "$(date --iso-8601=seconds)" >> "$status"
