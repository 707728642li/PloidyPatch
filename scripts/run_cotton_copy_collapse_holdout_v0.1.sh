#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then echo "usage: $0 PROJECT_ROOT" >&2; exit 2; fi
project_root=$(realpath "$1")
code_root=$project_root/code
python_bin=$project_root/envs/ploidypatch-dev/bin/python
source_gff=$project_root/data/derived/holdout_inputs/cotton_v0.1/hirsutum/primary_chromosomes.gff3
pair_root=$project_root/results/evaluator/cotton_holdout_v0.1/homeolog_pairs
pair_tsv=$pair_root/ghi_ad.AD.homeolog_pairs.tsv
pair_manifest=${pair_tsv}.manifest.json
policy=$code_root/config/copy_collapse_zero_retuning_policy_v0.1.tsv
seed=20260817
maximum_count=800
result_root=$project_root/benchmark/structure/copy_collapse_v0.1/ghi_ad/annotation_copy_collapse_seed${seed}
working_root=${result_root}.working

for required in "$python_bin" "$source_gff" "$pair_root/SHA256SUMS" "$pair_tsv" "$pair_manifest" "$policy"; do
    if [[ ! -s $required ]]; then echo "missing cotton holdout prerequisite: $required" >&2; exit 1; fi
done
(cd "$pair_root" && sha256sum -c SHA256SUMS >/dev/null)
pair_count=$(( $(wc -l < "$pair_tsv") - 1 ))
count=$maximum_count
if (( pair_count < count )); then count=$pair_count; fi
if (( count < 1 )); then echo "cotton pair set is empty" >&2; exit 1; fi
if [[ -e $result_root || -e $working_root ]]; then echo "refusing to overwrite cotton holdout" >&2; exit 1; fi
mkdir -p "$working_root/blind" "$working_root/evaluator/truth" "$working_root/evaluator/validation"
{
    printf 'field\tvalue\ncode_commit\t%s\n' "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'dataset_id\tghi_ad\nevent_type\tannotation_copy_collapse\n'
    printf 'count\t%s\nseed\t%s\nsplit\texternal_zero_retuning_holdout\n' "$count" "$seed"
    printf 'pair_access\tevaluator_only\npolicy_frozen_before_truth\ttrue\n'
    printf 'policy_sha256\t%s\n' "$(sha256sum "$policy" | awk '{print $1}')"
    printf 'pair_tsv_sha256\t%s\n' "$(sha256sum "$pair_tsv" | awk '{print $1}')"
    printf 'structure_module_sha256\t%s\n' \
        "$(sha256sum "$code_root/src/ploidypatch/structure_perturb.py" | awk '{print $1}')"
    printf 'score_module_sha256\t%s\n' \
        "$(sha256sum "$code_root/src/ploidypatch/score.py" | awk '{print $1}')"
} > "$working_root/run_contract.tsv"
cd "$code_root"
"$python_bin" -m ploidypatch.cli benchmark perturb \
    --gff "$source_gff" --output-dir "$working_root/blind" \
    --truth-dir "$working_root/evaluator/truth" \
    --event-type annotation_copy_collapse --pair-tsv "$pair_tsv" \
    --count "$count" --seed "$seed" \
    > "$working_root/perturb.stdout.json" 2> "$working_root/perturb.stderr.log"
perturbed=$working_root/blind/perturbed.gff3
truth=$working_root/evaluator/truth/hidden_truth.json
restored=$working_root/evaluator/restored.gff3
"$python_bin" -m ploidypatch.cli benchmark restore \
    --perturbed-gff "$perturbed" --truth "$truth" --output-gff "$restored" \
    > "$working_root/evaluator/restoration_report.json" \
    2> "$working_root/evaluator/restoration.stderr.log"
cmp -s "$source_gff" "$restored" || { echo "cotton restoration is not byte-identical" >&2; exit 1; }
for mode in noop oracle; do
    candidate=$perturbed; [[ $mode == oracle ]] && candidate=$source_gff
    "$python_bin" -m ploidypatch.cli benchmark score \
        --source-gff "$source_gff" --perturbed-gff "$perturbed" \
        --candidate-gff "$candidate" --truth "$truth" --include-event-details \
        > "$working_root/evaluator/validation/score_$mode.json" \
        2> "$working_root/evaluator/validation/score_$mode.stderr.log"
done
grep -q '"complete_cds_chain_recovery": 0' "$working_root/evaluator/validation/score_noop.json"
grep -q "\"complete_cds_chain_recovery\": $count" "$working_root/evaluator/validation/score_oracle.json"
grep -q '"grade": "pass"' "$working_root/evaluator/validation/score_noop.json"
grep -q '"grade": "pass"' "$working_root/evaluator/validation/score_oracle.json"
(
    cd "$working_root"
    find . -type f \( -name '*.gff3' -o -name '*.json' -o -name '*.tsv' \) \
        -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS
    sha256sum -c SHA256SUMS >/dev/null
)
du -sb "$working_root" > "$working_root/disk_bytes.txt"
mv "$working_root" "$result_root"
printf 'cotton zero-retuning holdout frozen: %s\n' "$result_root"
