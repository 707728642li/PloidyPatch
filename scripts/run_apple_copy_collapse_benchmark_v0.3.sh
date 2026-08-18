#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then echo "usage: $0 PROJECT_ROOT" >&2; exit 2; fi
project_root=$(realpath "$1")
code_root=$project_root/code
python_bin=$project_root/envs/ploidypatch-dev/bin/python
source_root=$project_root/data/derived/external_inputs/apple_v0.3
source_gff=$source_root/target_apple/primary_chromosomes.gff3
pair_root=$project_root/results/evaluator/apple_v0.3/truth_pairs
pair_tsv=$pair_root/intersection/pairs.tsv
protocol_root=$project_root/results/protocol_freezes/apple_external_v0.3
execution_root=$project_root/results/protocol_freezes/apple_external_v0.3_execution
policy=$protocol_root/policy.tsv
model=$project_root/results/models/support_conditioned_ranker_v0.3/model.json
seed=20260831
maximum_count=800
result_root=$project_root/benchmark/structure/copy_collapse_v0.3/mdx_gddh13/annotation_copy_collapse_seed${seed}
working_root=${result_root}.working

for required in "$python_bin" "$source_root/SHA256SUMS" "$source_gff" \
    "$pair_root/SHA256SUMS" "$pair_tsv" "$protocol_root/SHA256SUMS" \
    "$execution_root/SHA256SUMS" \
    "$policy" "$model" "$code_root/scripts/audit_copy_pair_selection_truth.py"; do
    [[ -s $required ]] || { echo "missing apple benchmark prerequisite: $required" >&2; exit 1; }
done
(cd "$source_root" && sha256sum -c SHA256SUMS >/dev/null)
(cd "$pair_root" && sha256sum -c SHA256SUMS >/dev/null)
(cd "$protocol_root" && sha256sum -c SHA256SUMS >/dev/null)
(cd "$execution_root" && sha256sum -c SHA256SUMS >/dev/null)
expected_self=$(awk -F '\t' '$1 == "scripts/run_apple_copy_collapse_benchmark_v0.3.sh" {print $3}' \
    "$execution_root/implementation_manifest.tsv")
[[ -n $expected_self && $(sha256sum "$code_root/scripts/run_apple_copy_collapse_benchmark_v0.3.sh" | awk '{print $1}') == "$expected_self" ]] || {
    echo "apple benchmark script differs from execution freeze" >&2; exit 1;
}
expected_model_sha=$(awk -F '\t' '$1 == "model_sha256" {print $2}' "$policy")
observed_model_sha=$(sha256sum "$model" | awk '{print $1}')
[[ -n $expected_model_sha && $observed_model_sha == "$expected_model_sha" ]] || {
    echo "frozen model hash disagrees with apple policy" >&2; exit 1;
}
for relative in src/ploidypatch/copy_pair_sampling.py \
                src/ploidypatch/structure_perturb.py src/ploidypatch/score.py; do
    expected=$(awk -F '\t' -v path="$relative" '$1 == path {print $2}' "$protocol_root/code_manifest.tsv")
    observed=$(sha256sum "$code_root/$relative" | awk '{print $1}')
    [[ -n $expected && $observed == "$expected" ]] || {
        echo "post-freeze module change detected: $relative" >&2; exit 1;
    }
done
[[ ! -e $result_root && ! -e $working_root ]] || {
    echo "refusing to overwrite apple v0.3 benchmark" >&2; exit 1;
}
mkdir -p "$working_root/blind" "$working_root/evaluator/pair_selection" \
    "$working_root/evaluator/truth" "$working_root/evaluator/validation"
{
    printf 'field\tvalue\n'
    printf 'code_commit\t%s\n' "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'dataset_id\tmdx_gddh13\nevent_type\tannotation_copy_collapse\n'
    printf 'requested_count\t%s\nseed\t%s\nsplit\tuntouched_external_v0.3\n' \
        "$maximum_count" "$seed"
    printf 'pair_access\tevaluator_only\npolicy_frozen_before_truth\ttrue\n'
    printf 'candidate_model_frozen_before_truth\ttrue\nautomatic_approval\tfalse\n'
    printf 'candidate_reference_access\tfalse\n'
    printf 'protocol_freeze_sha256sums_sha256\t%s\n' \
        "$(sha256sum "$protocol_root/SHA256SUMS" | awk '{print $1}')"
    printf 'model_sha256\t%s\npair_tsv_sha256\t%s\n' "$observed_model_sha" \
        "$(sha256sum "$pair_tsv" | awk '{print $1}')"
    printf 'pair_sampler_module_sha256\t%s\n' \
        "$(sha256sum "$code_root/src/ploidypatch/copy_pair_sampling.py" | awk '{print $1}')"
    printf 'structure_module_sha256\t%s\n' \
        "$(sha256sum "$code_root/src/ploidypatch/structure_perturb.py" | awk '{print $1}')"
    printf 'score_module_sha256\t%s\n' \
        "$(sha256sum "$code_root/src/ploidypatch/score.py" | awk '{print $1}')"
} > "$working_root/run_contract.tsv"

selected=$working_root/evaluator/pair_selection/selected_pairs.tsv
selection_decisions=$working_root/evaluator/pair_selection/decisions.tsv
cd "$code_root"
"$python_bin" -m ploidypatch.cli benchmark sample-copy-pairs \
    --source-gff "$source_gff" --pairs "$pair_tsv" \
    --count "$maximum_count" --seed "$seed" \
    --output-pairs "$selected" --decisions-tsv "$selection_decisions" \
    > "$working_root/evaluator/pair_selection/stdout.json" \
    2> "$working_root/evaluator/pair_selection/stderr.log"
selected_count=$(( $(wc -l < "$selected") - 1 ))
[[ $selected_count -gt 0 ]] || { echo "apple selected pair set is empty" >&2; exit 1; }
printf 'selected_count\t%s\n' "$selected_count" >> "$working_root/run_contract.tsv"

"$python_bin" - "$selected" "$working_root/evaluator/pair_selection/evaluability.json" <<'PY'
import csv
import json
import sys
from collections import Counter

source, output = sys.argv[1:]
with open(source, encoding="utf-8", newline="") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))
chromosomes = Counter(row["target_seqid"] for row in rows)
complexities = Counter(row["coding_complexity_bin"] for row in rows)
expected_bins = ("one", "two_to_three", "four_to_six", "seven_plus")
report = {
    "schema_version": "ploidypatch.apple_external_evaluability.v1",
    "events": len(rows),
    "target_chromosomes": len(chromosomes),
    "events_by_target_chromosome": dict(sorted(chromosomes.items())),
    "events_by_complexity_bin": {
        name: complexities.get(name, 0) for name in expected_bins
    },
}
report["gates"] = {
    "minimum_500_events": len(rows) >= 500,
    "minimum_15_target_chromosomes": len(chromosomes) >= 15,
    "minimum_20_each_complexity_bin": all(
        complexities.get(name, 0) >= 20 for name in expected_bins
    ),
}
report["formal_evaluable"] = all(report["gates"].values())
with open(output, "x", encoding="utf-8") as handle:
    json.dump(report, handle, indent=2, sort_keys=True)
    handle.write("\n")
PY

"$python_bin" -m ploidypatch.cli benchmark perturb \
    --gff "$source_gff" --output-dir "$working_root/blind" \
    --truth-dir "$working_root/evaluator/truth" \
    --event-type annotation_copy_collapse --pair-tsv "$selected" \
    --count "$selected_count" --seed "$seed" \
    > "$working_root/perturb.stdout.json" 2> "$working_root/perturb.stderr.log"
perturbed=$working_root/blind/perturbed.gff3
truth=$working_root/evaluator/truth/hidden_truth.json
restored=$working_root/evaluator/restored.gff3
"$python_bin" scripts/audit_copy_pair_selection_truth.py \
    --selected-pairs "$selected" --truth "$truth" \
    --output-json "$working_root/evaluator/validation/pair_truth_audit.json"
"$python_bin" -m ploidypatch.cli benchmark restore \
    --perturbed-gff "$perturbed" --truth "$truth" --output-gff "$restored" \
    > "$working_root/evaluator/restoration_report.json" \
    2> "$working_root/evaluator/restoration.stderr.log"
cmp -s "$source_gff" "$restored" || {
    echo "apple restoration is not byte-identical" >&2; exit 1;
}
for mode in noop oracle; do
    candidate=$perturbed; [[ $mode == oracle ]] && candidate=$source_gff
    "$python_bin" -m ploidypatch.cli benchmark score \
        --source-gff "$source_gff" --perturbed-gff "$perturbed" \
        --candidate-gff "$candidate" --truth "$truth" --include-event-details \
        > "$working_root/evaluator/validation/score_$mode.json" \
        2> "$working_root/evaluator/validation/score_$mode.stderr.log"
done
grep -q '"complete_cds_chain_recovery": 0' \
    "$working_root/evaluator/validation/score_noop.json"
grep -q "\"complete_cds_chain_recovery\": $selected_count" \
    "$working_root/evaluator/validation/score_oracle.json"
grep -q '"grade": "pass"' "$working_root/evaluator/validation/score_noop.json"
grep -q '"grade": "pass"' "$working_root/evaluator/validation/score_oracle.json"
grep -q '"grade": "pass"' \
    "$working_root/evaluator/validation/pair_truth_audit.json"
(
    cd "$working_root"
    find . -type f \( -name '*.gff3' -o -name '*.json' -o -name '*.tsv' \) \
        -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS
    sha256sum -c SHA256SUMS >/dev/null
)
du -sb "$working_root" > "$working_root/disk_bytes.txt"
mv "$working_root" "$result_root"
printf 'apple v0.3 copy-collapse benchmark frozen: %s\n' "$result_root"
