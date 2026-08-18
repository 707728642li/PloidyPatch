#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then echo "usage: $0 PROJECT_ROOT" >&2; exit 2; fi
project_root=$(realpath "$1")
code_root=$project_root/code
state=$project_root/logs/copy_collapse/maize_v2_evaluator_pairs
commit=${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}
if [[ -e $state ]]; then echo "maize evaluator-pair state exists" >&2; exit 1; fi
mkdir -p "$state"
record_failure() {
    rc=$?
    printf 'pipeline_failed\t%s\texit_%s\n' "$(date --iso-8601=seconds)" "$rc" \
        >> "$state/status.tsv"
    exit "$rc"
}
trap record_failure ERR
printf 'started\t%s\ncode_commit\t%s\ntruth_generated\tfalse\n' \
    "$(date --iso-8601=seconds)" "$commit" > "$state/status.tsv"

for stage in prepare_maize_v2_wgdi_inputs run_maize_v2_outgroup_wgdi \
             infer_maize_v2_outgroup_pairs; do
    PLOIDYPATCH_CODE_COMMIT="$commit" bash "$code_root/scripts/$stage.sh" \
        "$project_root" > "$state/$stage.log" 2>&1
    printf '%s_complete\t%s\n' "$stage" "$(date --iso-8601=seconds)" \
        >> "$state/status.tsv"
done
printf 'pipeline_complete\t%s\ntruth_generated\tfalse\n' \
    "$(date --iso-8601=seconds)" >> "$state/status.tsv"
