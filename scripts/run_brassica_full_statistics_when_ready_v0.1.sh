#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "usage: $0 PROJECT_ROOT" >&2
    exit 2
fi

project_root=$(realpath "$1")
code_root=$project_root/code
python_bin=$project_root/envs/ploidypatch-dev/bin/python
pipeline_state=$project_root/logs/copy_collapse/brassica_accelerated_pipeline_v0.1
method_root=$project_root/results/copy_collapse/model_development/brassica_method_trio_v0.1
transfer_root=$project_root/results/copy_collapse/model_development/brassica_glycine_model_transfer_v0.1
result_root=$project_root/results/copy_collapse/statistics/brassica_full_method_comparison_v0.1
working_root=${result_root}.working

while [[ ! -s $transfer_root/SHA256SUMS ]]; do
    watcher_pid=$(cat "$pipeline_state/watcher.pid")
    if ! kill -0 "$watcher_pid" 2>/dev/null; then
        echo "Brassica accelerated pipeline stopped before transfer result" >&2
        exit 1
    fi
    sleep 30
done

for frozen in "$method_root" "$transfer_root"; do
    (cd "$frozen" && sha256sum -c SHA256SUMS >/dev/null)
done

declare -A scores=(
    [miniprot]="$project_root/results/copy_collapse/miniprot_brassica_v0.1/score.json"
    [gemoma]="$method_root/methods/gemoma/score.json"
    [lifton]="$method_root/methods/lifton/score.json"
    [union]="$method_root/consensus/union/score.json"
    [consensus2]="$method_root/consensus/support2/score.json"
    [consensus3]="$method_root/consensus/support3/score.json"
    [glycine_rank_model]="$transfer_root/evaluator/paired_score.json"
)
for required in "$python_bin" "${scores[@]}"; do
    if [[ ! -s $required ]]; then
        echo "missing or empty Brassica statistical input: $required" >&2
        exit 1
    fi
done
if [[ -e $result_root || -e $working_root ]]; then
    echo "refusing to overwrite Brassica statistics: $result_root" >&2
    exit 1
fi
mkdir -p "$working_root"

{
    printf 'field\tvalue\n'
    printf 'code_commit\t%s\n' "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'benchmark\tbrassica_annotation_copy_collapse_v0.1\n'
    printf 'split\tpost_holdout_secondary_development\n'
    printf 'formal_holdout_claim_allowed\tfalse\n'
    printf 'replicates\t20000\nseed\t20260807\nalpha\t0.05\n'
    printf 'design\tpaired_events_stratified_by_event_type\n'
    printf 'primary_metric\tcomplete_cds_chain_recovery\n'
    printf 'glycine_model_parameter_retuning\tfalse\n'
} > "$working_root/run_contract.tsv"
{
    printf 'label\tbytes\tsha256\tpath\n'
    for label in miniprot gemoma lifton union consensus2 consensus3 \
                 glycine_rank_model; do
        path=${scores[$label]}
        printf '%s\t%s\t%s\t%s\n' "$label" "$(stat -Lc %s "$path")" \
            "$(sha256sum "$path" | awk '{print $1}')" "$path"
    done
} > "$working_root/input_manifest.tsv"

cd "$code_root"
"$python_bin" -m ploidypatch.cli benchmark bootstrap-events \
    --score "miniprot=${scores[miniprot]}" \
    --score "gemoma=${scores[gemoma]}" \
    --score "lifton=${scores[lifton]}" \
    --score "union=${scores[union]}" \
    --score "consensus2=${scores[consensus2]}" \
    --score "consensus3=${scores[consensus3]}" \
    --score "glycine_rank_model=${scores[glycine_rank_model]}" \
    --output-json "$working_root/event_recovery.json" \
    --metric complete_cds_chain_recovery \
    --replicates 20000 --seed 20260807 --alpha 0.05 \
    > "$working_root/event_recovery.stdout.json"
if [[ ! -s $working_root/event_recovery.json ]] || \
   ! grep -q '"schema_version"' "$working_root/event_recovery.json"; then
    echo "Brassica paired event bootstrap failed validation" >&2
    exit 1
fi
(
    cd "$working_root"
    find . -type f \( -name '*.json' -o -name '*.tsv' \) \
        -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS
    sha256sum -c SHA256SUMS >/dev/null
)
du -sb "$working_root" > "$working_root/disk_bytes.txt"
mv "$working_root" "$result_root"
printf 'full Brassica method statistics frozen: %s\n' "$result_root"
