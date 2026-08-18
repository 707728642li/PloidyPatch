#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then echo "usage: $0 PROJECT_ROOT" >&2; exit 2; fi
project_root=$(realpath "$1")
code_root=$project_root/code
dev_python=$project_root/envs/ploidypatch-dev/bin/python
model_python=$project_root/envs/ploidypatch-model/bin/python
benchmark=$project_root/benchmark/structure/copy_collapse_v0.2/zma_maize1/annotation_copy_collapse_seed20260829
truth=$benchmark/evaluator/truth/hidden_truth.json
rank_root=$project_root/results/copy_collapse/holdout/maize_v2_homeolog_ranker
features=$rank_root/blind/copy_features.tsv
scores=$rank_root/blind/scores.tsv
policy=$code_root/config/maize_v2_zero_retuning_policy.tsv
result_root=$project_root/results/copy_collapse/holdout/maize_v2_homeolog_ranker_evaluation
working_root=${result_root}.working

[[ -s $rank_root/SHA256SUMS ]] || { echo "unfrozen maize rank scores" >&2; exit 1; }
(cd "$rank_root" && sha256sum -c SHA256SUMS >/dev/null)
for required in "$dev_python" "$model_python" "$truth" "$features" "$scores" \
                "$policy" "$code_root/scripts/evaluate_frozen_homeolog_ranker_v0.2.py"; do
    [[ -s $required ]] || { echo "missing maize rank evaluator input: $required" >&2; exit 1; }
done
if [[ -e $result_root || -e $working_root ]]; then
    echo "refusing to overwrite maize rank evaluation" >&2; exit 1
fi
mkdir -p "$working_root/evaluator" "$working_root/logs"
{
    printf 'field\tvalue\ncode_commit\t%s\n' \
        "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'split\texternal_zero_retuning_holdout_v0.2\n'
    printf 'candidate_freeze_verified_before_truth_access\ttrue\n'
    printf 'evaluator_truth_access\ttrue\nmodel_refit\tfalse\nthreshold_tuning\tfalse\n'
    printf 'bootstrap_replicates\t20000\nbootstrap_seed\t20260829\n'
    printf 'bootstrap_unit\ttarget_seqid_chromosome\n'
    printf 'minimum_topology_positive_coverage\t0.70\n'
    printf 'automatic_approval\tfalse\n'
} > "$working_root/run_contract.tsv"
{
    printf 'role\tbytes\tsha256\tpath\n'
    for entry in "rank_freeze:$rank_root/SHA256SUMS" \
                 "features:$features" "scores:$scores" \
                 "hidden_truth:$truth" "policy:$policy"; do
        role=${entry%%:*}; path=${entry#*:}
        printf '%s\t%s\t%s\t%s\n' "$role" "$(stat -Lc %s "$path")" \
            "$(sha256sum "$path" | awk '{print $1}')" "$path"
    done
} > "$working_root/evaluator_input_manifest.tsv"

cd "$code_root"
/usr/bin/time -v -o "$working_root/logs/label_features.time.txt" \
    "$dev_python" -m ploidypatch.cli benchmark label-copy-features \
        --features "$features" --truth "$truth" \
        --output-tsv "$working_root/evaluator/labeled_features.tsv" \
        > "$working_root/logs/label_features.stdout.json" \
        2> "$working_root/logs/label_features.stderr.log"
/usr/bin/time -v -o "$working_root/logs/evaluate_ranker.time.txt" \
    "$model_python" "$code_root/scripts/evaluate_frozen_homeolog_ranker_v0.2.py" \
        --scores "$scores" \
        --labels "$working_root/evaluator/labeled_features.tsv" \
        --truth-event-count 800 --bootstrap-replicates 20000 \
        --seed 20260829 --minimum-topology-coverage 0.70 \
        --output-json "$working_root/evaluator/evaluation.json" \
        > "$working_root/logs/evaluate_ranker.stdout.log" \
        2> "$working_root/logs/evaluate_ranker.stderr.log"
for output in "$working_root/evaluator/labeled_features.tsv" \
              "$working_root/evaluator/labeled_features.tsv.manifest.json" \
              "$working_root/evaluator/evaluation.json"; do
    [[ -s $output ]] || { echo "missing maize rank evaluation output: $output" >&2; exit 1; }
done
"$model_python" - "$working_root/evaluator/evaluation.json" <<'PY'
import json, sys
report = json.load(open(sys.argv[1], encoding="utf-8"))
if report["schema_version"] != "ploidypatch.frozen_homeolog_ranker_evaluation.v2":
    raise SystemExit("unexpected maize rank evaluation schema")
if report["policy"]["model_refit_on_external_species"] is not False:
    raise SystemExit("external species model refit detected")
print(json.dumps(report["gates"], sort_keys=True))
PY
(
    cd "$working_root"
    find . -type f \( -name '*.json' -o -name '*.tsv' \) -print0 \
        | sort -z | xargs -0 sha256sum > SHA256SUMS
    sha256sum -c SHA256SUMS >/dev/null
)
du -sb "$working_root" > "$working_root/disk_bytes.txt"
mv "$working_root" "$result_root"
printf 'maize frozen homeolog ranker evaluation: %s\n' "$result_root"
