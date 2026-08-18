#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then echo "usage: $0 PROJECT_ROOT" >&2; exit 2; fi
project_root=$(realpath "$1")
code_root=$project_root/code
method_eval=$project_root/results/copy_collapse/holdout/maize_v2_method_trio_evaluation
method_state=$project_root/logs/baseline/maize_v2/post_pipeline.tsv
state=$project_root/logs/copy_collapse/maize_v2_homeolog_pipeline
status=$state/status.tsv
commit=${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}
mkdir -p "$state"
[[ ! -e $status ]] || { echo "maize homeolog pipeline status exists" >&2; exit 1; }
printf 'watcher_started\t%s\ncode_commit\t%s\n' \
    "$(date --iso-8601=seconds)" "$commit" > "$status"

while [[ ! -s $method_eval/SHA256SUMS ]]; do
    if [[ -s $method_state ]] && grep -q '^upstream_failed' "$method_state"; then
        printf 'dependency_failed\t%s\tmaize_method_pipeline\n' \
            "$(date --iso-8601=seconds)" >> "$status"
        exit 1
    fi
    if ! pgrep -f 'watch_maize_v2_method_pipeline.sh' >/dev/null; then
        printf 'dependency_stopped\t%s\tmaize_method_pipeline\n' \
            "$(date --iso-8601=seconds)" >> "$status"
        exit 1
    fi
    sleep 30
done
(cd "$method_eval" && sha256sum -c SHA256SUMS >/dev/null)
printf 'method_evaluation_verified\t%s\n' "$(date --iso-8601=seconds)" >> "$status"

PLOIDYPATCH_CODE_COMMIT="$commit" \
    bash "$code_root/scripts/run_maize_v2_union_self_wgd.sh" "$project_root" \
    > "$state/run_maize_v2_union_self_wgd.log" 2>&1
printf 'blind_self_wgd_frozen\t%s\n' "$(date --iso-8601=seconds)" >> "$status"
PLOIDYPATCH_CODE_COMMIT="$commit" \
    bash "$code_root/scripts/score_maize_v2_homeolog_ranker_blind.sh" "$project_root" \
    > "$state/score_maize_v2_homeolog_ranker_blind.log" 2>&1
printf 'blind_rank_scores_frozen\t%s\n' "$(date --iso-8601=seconds)" >> "$status"
PLOIDYPATCH_CODE_COMMIT="$commit" \
    bash "$code_root/scripts/evaluate_maize_v2_homeolog_ranker.sh" "$project_root" \
    > "$state/evaluate_maize_v2_homeolog_ranker.log" 2>&1
printf 'formal_rank_evaluation_complete\t%s\n' "$(date --iso-8601=seconds)" >> "$status"
