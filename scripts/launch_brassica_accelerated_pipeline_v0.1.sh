#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then echo "usage: $0 PROJECT_ROOT" >&2; exit 2; fi
project_root=$(realpath "$1")
code_root=$project_root/code
pipeline_state=$project_root/logs/copy_collapse/brassica_accelerated_pipeline_v0.1
watcher=$code_root/scripts/watch_brassica_accelerated_pipeline_v0.1.sh
if [[ -e $pipeline_state ]]; then echo "pipeline watcher state already exists" >&2; exit 1; fi
if [[ ! -s $watcher ]]; then echo "pipeline watcher script is absent" >&2; exit 1; fi
mkdir -p "$pipeline_state"
commit=${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}
nohup setsid bash "$watcher" "$project_root" "$commit" \
    > "$pipeline_state/watcher.log" 2>&1 &
printf '%s\n' "$!" > "$pipeline_state/watcher.pid"
sleep 2
pid=$(cat "$pipeline_state/watcher.pid")
if ! kill -0 "$pid" 2>/dev/null; then
    echo "accelerated pipeline watcher exited during startup" >&2; exit 1
fi
printf 'Brassica accelerated pipeline watcher PID: %s\n' "$pid"
