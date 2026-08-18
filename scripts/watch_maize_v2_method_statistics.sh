#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then echo "usage: $0 PROJECT_ROOT" >&2; exit 2; fi
project_root=$(realpath "$1")
code_root=$project_root/code
evaluation=$project_root/results/copy_collapse/holdout/maize_v2_method_trio_evaluation
dependency_state=$project_root/logs/baseline/maize_v2/post_pipeline.tsv
state=$project_root/logs/copy_collapse/maize_v2_method_statistics
status=$state/status.tsv
commit=${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}
mkdir -p "$state"
[[ ! -e $status ]] || { echo "maize statistics watcher status exists" >&2; exit 1; }
printf 'watcher_started\t%s\ncode_commit\t%s\n' \
    "$(date --iso-8601=seconds)" "$commit" > "$status"
while [[ ! -s $evaluation/SHA256SUMS ]]; do
    if [[ -s $dependency_state ]] && grep -q '^upstream_failed' "$dependency_state"; then
        printf 'dependency_failed\t%s\n' "$(date --iso-8601=seconds)" >> "$status"
        exit 1
    fi
    if ! pgrep -f 'watch_maize_v2_method_pipeline.sh' >/dev/null; then
        printf 'dependency_stopped\t%s\n' "$(date --iso-8601=seconds)" >> "$status"
        exit 1
    fi
    sleep 30
done
PLOIDYPATCH_CODE_COMMIT="$commit" \
    bash "$code_root/scripts/run_maize_v2_method_statistics.sh" "$project_root" \
    > "$state/run_maize_v2_method_statistics.log" 2>&1
printf 'statistics_complete\t%s\n' "$(date --iso-8601=seconds)" >> "$status"
