#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then echo "usage: $0 PROJECT_ROOT" >&2; exit 2; fi
project_root=$(realpath "$1")
code_root=$project_root/code
python_bin=$project_root/envs/ploidypatch-model/bin/python
rank_root=$project_root/results/copy_collapse/holdout/maize_v2_homeolog_ranker
evaluation_root=$project_root/results/copy_collapse/holdout/maize_v2_homeolog_ranker_evaluation
model_root=$project_root/results/copy_collapse/model_development/homeolog_ranker_v0.2
scores=$rank_root/blind/scores.tsv
topology=$rank_root/blind/topology_features.tsv
labels=$evaluation_root/evaluator/labeled_features.tsv
model=$model_root/model.json
script=$code_root/scripts/diagnose_maize_v2_homeolog_ranker.py
result_root=$project_root/results/copy_collapse/diagnostics/maize_v2_homeolog_ranker_v0.2
working_root=${result_root}.working

for frozen in "$rank_root" "$evaluation_root" "$model_root"; do
    [[ -s $frozen/SHA256SUMS ]] || { echo "unfrozen diagnostic input: $frozen" >&2; exit 1; }
    (cd "$frozen" && sha256sum -c SHA256SUMS >/dev/null)
done
for required in "$python_bin" "$scores" "$topology" "$labels" "$model" "$script"; do
    [[ -s $required ]] || { echo "missing maize diagnostic input: $required" >&2; exit 1; }
done
[[ ! -e $result_root ]] || { echo "refusing to overwrite maize diagnostic" >&2; exit 1; }
if [[ -e $working_root && ${PLOIDYPATCH_RESUME_DIAGNOSTIC:-false} != true ]]; then
    echo "maize diagnostic working directory exists; explicit resume required" >&2
    exit 1
fi
mkdir -p "$working_root"
{
    printf 'field\tvalue\ncode_commit\t%s\n' \
        "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'analysis_role\tpost_hoc_failure_diagnosis_not_model_selection\n'
    printf 'maize_refit\tfalse\nmaize_threshold_tuning\tfalse\n'
    printf 'changes_formal_gate\tfalse\nautomatic_approval\tfalse\n'
} > "$working_root/run_contract.tsv"
{
    printf 'role\tbytes\tsha256\tpath\n'
    for entry in "scores:$scores" "topology:$topology" "labels:$labels" "model:$model"; do
        role=${entry%%:*}; path=${entry#*:}
        printf '%s\t%s\t%s\t%s\n' "$role" "$(stat -Lc %s "$path")" \
            "$(sha256sum "$path" | awk '{print $1}')" "$path"
    done
} > "$working_root/input_manifest.tsv"

if [[ ! -s $working_root/diagnostic.json ]]; then
    cd "$code_root"
    /usr/bin/time -v -o "$working_root/resource.time.txt" \
        env PYTHONPATH="$code_root/src" "$python_bin" "$script" \
        --scores "$scores" --labels "$labels" --topology "$topology" \
        --model "$model" --output-json "$working_root/diagnostic.json" \
        > "$working_root/stdout.log" 2> "$working_root/stderr.log"
else
    printf 'resumed_existing_diagnostic\t%s\n' "$(date --iso-8601=seconds)" \
        >> "$working_root/run_contract.tsv"
fi
"$python_bin" - "$working_root/diagnostic.json" <<'PY'
import json, sys
report = json.load(open(sys.argv[1], encoding="utf-8"))
if report.get("schema_version") != "ploidypatch.maize_homeolog_ranker_diagnostic.v1":
    raise SystemExit("unexpected maize diagnostic schema")
if report.get("policy") != {
    "automatic_approval": False,
    "changes_formal_gate": False,
    "maize_refit": False,
    "maize_threshold_tuning": False,
}:
    raise SystemExit("maize diagnostic policy gate failed")
if report.get("counts", {}).get("candidates") != 15636:
    raise SystemExit("maize diagnostic candidate universe changed")
PY
(
    cd "$working_root"
    find . -type f \( -name '*.json' -o -name '*.tsv' \) -print0 \
        | sort -z | xargs -0 sha256sum > SHA256SUMS
    sha256sum -c SHA256SUMS >/dev/null
)
du -sb "$working_root" > "$working_root/disk_bytes.txt"
mv "$working_root" "$result_root"
printf 'maize post-hoc homeolog diagnostic frozen: %s\n' "$result_root"
