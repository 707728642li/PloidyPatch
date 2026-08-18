#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then echo "usage: $0 PROJECT_ROOT" >&2; exit 2; fi
project_root=$(realpath "$1")
code_root=$project_root/code
model_env=$project_root/envs/ploidypatch-model
feature_root=$project_root/results/copy_collapse/candidate_pool_v0.3_features
evaluation=$project_root/results/copy_collapse/model_development/candidate_pool_ranker_v0.3_evaluation/evaluation.json
protocol=$code_root/docs/RANKER_V0.3_DEVELOPMENT_PROTOCOL.md
result_root=$project_root/results/models/support_conditioned_ranker_v0.3
working_root=${result_root}.working

for required in "$model_env/bin/python" "$evaluation" "$protocol" \
    "$feature_root/glycine/labeled_features.tsv" \
    "$feature_root/glycine/topology_features.tsv" \
    "$feature_root/brassica/labeled_features.tsv" \
    "$feature_root/brassica/topology_features.tsv"; do
    [[ -s $required ]] || { echo "missing v0.3 model-freeze input: $required" >&2; exit 1; }
done
[[ ! -e $result_root && ! -e $working_root ]] || {
    echo "refusing to overwrite v0.3 model freeze" >&2; exit 1;
}
mkdir -p "$working_root"
{
    printf 'field\tvalue\n'
    printf 'code_commit\t%s\n' "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'model_selection\talready_completed_in_canonical_development_evaluation\n'
    printf 'external_species_access\tfalse\n'
    printf 'external_labels_access\tfalse\n'
    printf 'automatic_approval\tfalse\n'
} > "$working_root/run_contract.tsv"

cd "$code_root"
/usr/bin/time -v -o "$working_root/resource.time.txt" \
    env PYTHONPATH="$code_root/src" \
    PLOIDYPATCH_CODE_COMMIT="${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}" \
    conda run -p "$model_env" --no-capture-output python \
    "$code_root/scripts/freeze_support_conditioned_ranker_v0.3.py" \
    --development "glycine=$feature_root/glycine/labeled_features.tsv,$feature_root/glycine/topology_features.tsv" \
    --development "brassica=$feature_root/brassica/labeled_features.tsv,$feature_root/brassica/topology_features.tsv" \
    --canonical-evaluation "$evaluation" --protocol "$protocol" \
    --output-json "$working_root/model.json" \
    > "$working_root/stdout.json" 2> "$working_root/stderr.log"

[[ -s $working_root/model.json ]] || { echo "missing frozen v0.3 model" >&2; exit 1; }
(
    cd "$working_root"
    sha256sum model.json run_contract.tsv stdout.json > SHA256SUMS
    sha256sum -c SHA256SUMS >/dev/null
)
mv "$working_root" "$result_root"
printf 'support-conditioned v0.3 ranker frozen: %s\n' "$result_root"

