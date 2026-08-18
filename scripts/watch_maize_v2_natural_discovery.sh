#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then echo "usage: $0 PROJECT_ROOT" >&2; exit 2; fi
project_root=$(realpath "$1")
code_root=$project_root/code
formal=$project_root/results/copy_collapse/holdout/maize_v2_homeolog_ranker_evaluation
formal_state=$project_root/logs/copy_collapse/maize_v2_homeolog_pipeline/status.tsv
state=$project_root/logs/natural/maize_v2_discovery
status=$state/status.tsv
mkdir -p "$state"
[[ ! -e $status ]] || { echo "maize natural discovery status exists" >&2; exit 1; }
printf 'watcher_started\t%s\ncode_commit\t%s\n' "$(date --iso-8601=seconds)" \
    "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}" > "$status"

while [[ ! -s $formal/SHA256SUMS ]]; do
    if [[ -s $formal_state ]] && grep -Eq '^(pipeline_failed|dependency_failed|dependency_stopped)' "$formal_state"; then
        printf 'dependency_failed\t%s\tformal_maize_homeolog_evaluation\n' \
            "$(date --iso-8601=seconds)" >> "$status"
        exit 1
    fi
    sleep 30
done
(cd "$formal" && sha256sum -c SHA256SUMS >/dev/null)
printf 'formal_evaluation_verified\t%s\n' "$(date --iso-8601=seconds)" >> "$status"

if ! PLOIDYPATCH_CODE_COMMIT="${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}" \
    bash "$code_root/scripts/build_maize_v2_natural_method_trio.sh" "$project_root" \
    > "$state/build_method_trio.log" 2>&1; then
    printf 'pipeline_failed\t%s\tbuild_method_trio\n' "$(date --iso-8601=seconds)" >> "$status"
    exit 1
fi
printf 'method_trio_frozen\t%s\n' "$(date --iso-8601=seconds)" >> "$status"

if ! PLOIDYPATCH_CODE_COMMIT="${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}" \
    bash "$code_root/scripts/run_maize_v2_natural_self_wgd.sh" "$project_root" \
    > "$state/run_self_wgd.log" 2>&1; then
    printf 'pipeline_failed\t%s\trun_self_wgd\n' "$(date --iso-8601=seconds)" >> "$status"
    exit 1
fi
printf 'self_wgd_frozen\t%s\n' "$(date --iso-8601=seconds)" >> "$status"

if ! PLOIDYPATCH_CODE_COMMIT="${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}" \
    bash "$code_root/scripts/score_maize_v2_natural_homeolog_ranker.sh" "$project_root" \
    > "$state/score_homeolog_ranker.log" 2>&1; then
    printf 'pipeline_failed\t%s\tscore_homeolog_ranker\n' "$(date --iso-8601=seconds)" >> "$status"
    exit 1
fi
printf 'candidate_and_review_freeze_complete\t%s\n' \
    "$(date --iso-8601=seconds)" >> "$status"
