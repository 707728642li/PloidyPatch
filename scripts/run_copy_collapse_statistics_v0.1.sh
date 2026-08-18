#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "usage: $0 PROJECT_ROOT" >&2
    exit 2
fi

project_root=$(realpath "$1")
code_root=$project_root/code
python_bin=$project_root/envs/ploidypatch-dev/bin/python
brassica=$project_root/results/copy_collapse/miniprot_brassica_v0.1/score.json
glycine_miniprot=$project_root/results/copy_collapse/miniprot_glycine_v0.1/score.json
glycine_gemoma=$project_root/results/copy_collapse/gemoma_glycine_v0.1/score.json
glycine_lifton=$project_root/results/copy_collapse/lifton_glycine_v0.1/score.json
result_root=$project_root/results/copy_collapse/statistics/portability_and_transfer_trio_v0.1
working_root=${result_root}.working

for required in "$python_bin" "$brassica" "$glycine_miniprot" \
                "$glycine_gemoma" "$glycine_lifton"; do
    if [[ ! -s $required ]]; then
        echo "missing or empty copy-collapse statistical input: $required" >&2
        exit 1
    fi
done
if [[ -e $result_root || -e $working_root ]]; then
    echo "refusing to overwrite copy-collapse statistics: $result_root" >&2
    exit 1
fi
mkdir -p "$working_root"

{
    printf 'field\tvalue\n'
    printf 'code_commit\t%s\n' "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'replicates\t20000\n'
    printf 'seed\t20260807\n'
    printf 'alpha\t0.05\n'
    printf 'portability_design\tindependent_benchmark_sets\n'
    printf 'transfer_comparison_design\tpaired_events\n'
    printf 'syngap_included\tfalse_pending_matched_upstream_rerun\n'
} > "$working_root/run_contract.tsv"
{
    printf 'label\tbytes\tsha256\tpath\n'
    for entry in \
        "brassica_miniprot:$brassica" \
        "glycine_miniprot:$glycine_miniprot" \
        "glycine_gemoma:$glycine_gemoma" \
        "glycine_lifton:$glycine_lifton"; do
        label=${entry%%:*}
        path=${entry#*:}
        printf '%s\t%s\t%s\t%s\n' "$label" "$(stat -Lc %s "$path")" \
            "$(sha256sum "$path" | awk '{print $1}')" "$path"
    done
} > "$working_root/input_manifest.tsv"

cd "$code_root"
"$python_bin" -m ploidypatch.cli benchmark bootstrap-independent-events \
    --score "brassica_miniprot=$brassica" \
    --score "glycine_miniprot=$glycine_miniprot" \
    --output-json "$working_root/portability_event_recovery.json" \
    --metric complete_cds_chain_recovery \
    --replicates 20000 --seed 20260807 --alpha 0.05 \
    > "$working_root/portability_event_recovery.stdout.json"
"$python_bin" -m ploidypatch.cli benchmark bootstrap-confusion \
    --score "brassica_miniprot=$brassica" \
    --score "glycine_miniprot=$glycine_miniprot" \
    --output-json "$working_root/portability_strict_cds_confusion.json" \
    --section strict_cds_chain \
    --replicates 20000 --seed 20260807 --alpha 0.05 \
    > "$working_root/portability_strict_cds_confusion.stdout.json"
"$python_bin" -m ploidypatch.cli benchmark bootstrap-events \
    --score "miniprot=$glycine_miniprot" \
    --score "gemoma=$glycine_gemoma" \
    --score "lifton=$glycine_lifton" \
    --output-json "$working_root/glycine_transfer_trio_event_recovery.json" \
    --metric complete_cds_chain_recovery \
    --replicates 20000 --seed 20260807 --alpha 0.05 \
    > "$working_root/glycine_transfer_trio_event_recovery.stdout.json"

for output in "$working_root/portability_event_recovery.json" \
              "$working_root/portability_strict_cds_confusion.json" \
              "$working_root/glycine_transfer_trio_event_recovery.json"; do
    if [[ ! -s $output ]] || ! grep -q '"schema_version"' "$output"; then
        echo "copy-collapse statistical output failed validation: $output" >&2
        exit 1
    fi
done
(
    cd "$working_root"
    find . -type f \( -name '*.json' -o -name '*.tsv' \) \
        -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS
)
du -sb "$working_root" > "$working_root/disk_bytes.txt"
mv "$working_root" "$result_root"
printf 'copy-collapse statistics frozen: %s\n' "$result_root"
