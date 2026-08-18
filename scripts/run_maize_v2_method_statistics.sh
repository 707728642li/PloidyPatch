#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then echo "usage: $0 PROJECT_ROOT" >&2; exit 2; fi
project_root=$(realpath "$1")
code_root=$project_root/code
python_bin=$project_root/envs/ploidypatch-dev/bin/python
evaluation=$project_root/results/copy_collapse/holdout/maize_v2_method_trio_evaluation
result_root=$project_root/results/copy_collapse/statistics/maize_v2_method_trio
working_root=${result_root}.working
declare -A score=(
    [miniprot]="$evaluation/scores/methods/miniprot/score.json"
    [gemoma]="$evaluation/scores/methods/gemoma/score.json"
    [lifton]="$evaluation/scores/methods/lifton/score.json"
    [union]="$evaluation/scores/consensus/union/score.json"
    [consensus2]="$evaluation/scores/consensus/support2/score.json"
    [consensus3]="$evaluation/scores/consensus/support3/score.json"
)

[[ -s $evaluation/SHA256SUMS ]] || { echo "unfrozen maize method evaluation" >&2; exit 1; }
(cd "$evaluation" && sha256sum -c SHA256SUMS >/dev/null)
for required in "$python_bin" "${score[@]}"; do
    [[ -s $required ]] || { echo "missing maize statistical input: $required" >&2; exit 1; }
done
if [[ -e $result_root || -e $working_root ]]; then
    echo "refusing to overwrite maize method statistics" >&2; exit 1
fi
mkdir -p "$working_root"
{
    printf 'field\tvalue\ncode_commit\t%s\n' \
        "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'split\texternal_zero_retuning_holdout_v0.2\n'
    printf 'comparison_design\tpaired_hidden_events\n'
    printf 'replicates\t20000\nseed\t20260829\nalpha\t0.05\n'
} > "$working_root/run_contract.tsv"
{
    printf 'label\tbytes\tsha256\tpath\n'
    for label in miniprot gemoma lifton union consensus2 consensus3; do
        path=${score[$label]}
        printf '%s\t%s\t%s\t%s\n' "$label" "$(stat -Lc %s "$path")" \
            "$(sha256sum "$path" | awk '{print $1}')" "$path"
    done
} > "$working_root/input_manifest.tsv"

cd "$code_root"
common_scores=(
    --score "miniprot=${score[miniprot]}"
    --score "gemoma=${score[gemoma]}"
    --score "lifton=${score[lifton]}"
    --score "union=${score[union]}"
    --score "consensus2=${score[consensus2]}"
    --score "consensus3=${score[consensus3]}"
)
"$python_bin" -m ploidypatch.cli benchmark bootstrap-events \
    "${common_scores[@]}" --output-json "$working_root/event_recovery.json" \
    --metric complete_cds_chain_recovery --replicates 20000 \
    --seed 20260829 --alpha 0.05 \
    > "$working_root/event_recovery.stdout.json"
"$python_bin" -m ploidypatch.cli benchmark bootstrap-confusion \
    "${common_scores[@]}" --output-json "$working_root/strict_cds_confusion.json" \
    --section strict_cds_chain --replicates 20000 \
    --seed 20260829 --alpha 0.05 \
    > "$working_root/strict_cds_confusion.stdout.json"
for output in "$working_root/event_recovery.json" \
              "$working_root/strict_cds_confusion.json"; do
    [[ -s $output ]] && grep -q schema_version "$output" || {
        echo "maize statistical output validation failed: $output" >&2; exit 1;
    }
done
(
    cd "$working_root"
    find . -type f \( -name '*.json' -o -name '*.tsv' \) -print0 \
        | sort -z | xargs -0 sha256sum > SHA256SUMS
    sha256sum -c SHA256SUMS >/dev/null
)
du -sb "$working_root" > "$working_root/disk_bytes.txt"
mv "$working_root" "$result_root"
printf 'maize paired method statistics frozen: %s\n' "$result_root"
