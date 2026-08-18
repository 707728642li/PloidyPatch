#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "usage: $0 PROJECT_ROOT" >&2
    exit 2
fi

project_root=$(realpath "$1")
code_root=$project_root/code
python_bin=$project_root/envs/ploidypatch-dev/bin/python
gemoma_score=$project_root/results/baselines/gemoma_v1.9/glycine_v0.1/evaluation/score.json
lifton_score=$project_root/results/baselines/lifton_v1.0.11/glycine_v0.1/evaluation/score.json
syngap_score=$project_root/results/baselines/syngap_v1.2.5/glycine_v0.1/genblastg/evaluation_phase_normalized_v0.2/score.json
result_root=$project_root/results/baselines/matched_bootstrap/glycine_v0.1
working_root=${result_root}.working

for required in "$python_bin" "$gemoma_score" "$lifton_score" "$syngap_score"; do
    if [[ ! -s $required ]]; then
        echo "missing matched-bootstrap prerequisite: $required" >&2
        exit 1
    fi
done
if [[ -e $result_root || -e $working_root ]]; then
    echo "refusing to overwrite matched-bootstrap result: $result_root" >&2
    exit 1
fi
mkdir -p "$working_root"
{
    printf 'field\tvalue\n'
    printf 'code_commit\t%s\n' "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'benchmark\tgma_v21_annotation_missing_gene_seed20260809\n'
    printf 'methods\tGeMoMa-1.9,LiftOn-1.0.11,SynGAP-1.2.5-genblastg\n'
    printf 'bootstrap_replicates\t10000\n'
    printf 'bootstrap_seed\t20260807\n'
    printf 'bootstrap_alpha\t0.05\n'
    printf 'resampling\tpaired_event_type_stratified\n'
} > "$working_root/run_contract.tsv"
{
    printf 'method\tbytes\tsha256\tpath\n'
    for entry in \
        "gemoma:$gemoma_score" \
        "lifton:$lifton_score" \
        "syngap:$syngap_score"; do
        method=${entry%%:*}
        path=${entry#*:}
        printf '%s\t%s\t%s\t%s\n' "$method" "$(stat -Lc %s "$path")" \
            "$(sha256sum "$path" | awk '{print $1}')" "$path"
    done
} > "$working_root/input_manifest.tsv"

metrics=(
    complete_cds_chain_recovery
    complete_transcript_recovery
    exact_cds_gene_grouping
    exact_gene_grouping
)
cd "$code_root"
for metric in "${metrics[@]}"; do
    /usr/bin/time -v -o "$working_root/${metric}.resource.time.txt" \
        "$python_bin" -m ploidypatch.cli benchmark bootstrap-events \
            --score "gemoma=$gemoma_score" \
            --score "lifton=$lifton_score" \
            --score "syngap=$syngap_score" \
            --metric "$metric" \
            --replicates 10000 \
            --seed 20260807 \
            --alpha 0.05 \
            --output-json "$working_root/${metric}.bootstrap.json" \
            > "$working_root/${metric}.stdout.json" \
            2> "$working_root/${metric}.stderr.log"
done

for metric in "${metrics[@]}"; do
    if [[ ! -s $working_root/${metric}.bootstrap.json ]]; then
        echo "matched-bootstrap output is missing: $metric" >&2
        exit 1
    fi
done
(
    cd "$working_root"
    find . -type f \( -name '*.json' -o -name '*.tsv' \) -print0 \
        | sort -z | xargs -0 sha256sum > SHA256SUMS
)
du -sb "$working_root" > "$working_root/disk_bytes.txt"
mv "$working_root" "$result_root"
printf 'matched baseline bootstrap completed: %s\n' "$result_root"
