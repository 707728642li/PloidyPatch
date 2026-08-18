#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "usage: $0 PROJECT_ROOT" >&2
    exit 2
fi

project_root=$(realpath "$1")
code_root=$project_root/code
python_bin=$project_root/envs/ploidypatch-dev/bin/python
result_root=$project_root/results/copy_collapse/statistics/glycine_full_method_comparison_v0.1
working_root=${result_root}.working

declare -A scores=(
    [miniprot]="$project_root/results/copy_collapse/miniprot_glycine_v0.1/score.json"
    [gemoma]="$project_root/results/copy_collapse/gemoma_glycine_v0.1/score.json"
    [lifton]="$project_root/results/copy_collapse/lifton_glycine_v0.1/score.json"
    [syngap]="$project_root/results/copy_collapse/syngap_glycine_v0.1/evaluation/score.json"
    [consensus2]="$project_root/results/copy_collapse/consensus_glycine_v0.1/support2/score.json"
    [consensus3]="$project_root/results/copy_collapse/consensus_glycine_v0.1/support3/score.json"
    [consensus2_wgd]="$project_root/results/copy_collapse/wgd_reanchor_glycine_v0.1/evaluation/score.json"
    [rank_model]="$project_root/results/copy_collapse/model_development/glycine_copy_model_evaluation_v0.1/evaluator/score.json"
)

for required in "$python_bin" "${scores[@]}"; do
    if [[ ! -s $required ]]; then
        echo "missing or empty statistical input: $required" >&2
        exit 1
    fi
done
if [[ -e $result_root || -e $working_root ]]; then
    echo "refusing to overwrite statistics: $result_root" >&2
    exit 1
fi
mkdir -p "$working_root"

{
    printf 'field\tvalue\n'
    printf 'code_commit\t%s\n' "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'benchmark\tglycine_annotation_copy_collapse_v0.1\n'
    printf 'replicates\t20000\n'
    printf 'seed\t20260807\n'
    printf 'alpha\t0.05\n'
    printf 'design\tpaired_events_stratified_by_event_type\n'
    printf 'primary_metric\tcomplete_cds_chain_recovery\n'
    printf 'rank_model_status\treview_only_frozen_threshold\n'
} > "$working_root/run_contract.tsv"

{
    printf 'label\tbytes\tsha256\tpath\n'
    for label in miniprot gemoma lifton syngap consensus2 consensus3 \
                 consensus2_wgd rank_model; do
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
    --score "syngap=${scores[syngap]}" \
    --score "consensus2=${scores[consensus2]}" \
    --score "consensus3=${scores[consensus3]}" \
    --score "consensus2_wgd=${scores[consensus2_wgd]}" \
    --score "rank_model=${scores[rank_model]}" \
    --output-json "$working_root/event_recovery.json" \
    --metric complete_cds_chain_recovery \
    --replicates 20000 --seed 20260807 --alpha 0.05 \
    > "$working_root/event_recovery.stdout.json"

if [[ ! -s $working_root/event_recovery.json ]] || \
   ! grep -q '"schema_version"' "$working_root/event_recovery.json"; then
    echo "paired event bootstrap failed validation" >&2
    exit 1
fi
(
    cd "$working_root"
    find . -type f \( -name '*.json' -o -name '*.tsv' \) \
        -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS
)
du -sb "$working_root" > "$working_root/disk_bytes.txt"
mv "$working_root" "$result_root"
printf 'full Glycine method statistics frozen: %s\n' "$result_root"
