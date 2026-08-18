#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then echo "usage: $0 PROJECT_ROOT" >&2; exit 2; fi
project_root=$(realpath "$1")
code_root=$project_root/code
model_env=$project_root/envs/ploidypatch-model
feature_root=$project_root/results/copy_collapse/candidate_pool_v0.3_features
pool_root=$project_root/results/copy_collapse/candidate_pool_v0.3
protocol=$code_root/docs/RANKER_V0.3_DEVELOPMENT_PROTOCOL.md
output=$project_root/results/copy_collapse/model_development/candidate_pool_ranker_v0.3_evaluation
[[ -x $model_env/bin/python ]] || { echo "missing PloidyPatch model environment" >&2; exit 1; }
[[ -s $protocol ]] || { echo "missing v0.3 ranker protocol" >&2; exit 1; }
[[ ! -e $output && ! -e ${output}.partial ]] || {
    echo "refusing to overwrite v0.3 ranker evaluation" >&2; exit 1;
}
for species in glycine brassica cotton maize; do
    [[ -s $feature_root/$species/SHA256SUMS ]] || { echo "unfrozen features: $species" >&2; exit 1; }
    [[ -s $pool_root/$species/SHA256SUMS ]] || { echo "unfrozen pool: $species" >&2; exit 1; }
    (cd "$feature_root/$species" && sha256sum -c SHA256SUMS >/dev/null)
    (cd "$pool_root/$species" && sha256sum -c SHA256SUMS >/dev/null)
done
mkdir -p "$project_root/logs/candidate_pool_ranker_v0.3"
cd "$code_root"
/usr/bin/time -v -o "$project_root/logs/candidate_pool_ranker_v0.3/evaluation.resource.time.txt" \
    env PYTHONPATH="$code_root/src" \
    PLOIDYPATCH_CODE_COMMIT="${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}" \
    conda run -p "$model_env" --no-capture-output python \
    "$code_root/scripts/evaluate_candidate_pool_rankers_v0.3.py" \
    --protocol "$protocol" \
    --dataset "glycine=development,$feature_root/glycine/labeled_features.tsv,$feature_root/glycine/topology_features.tsv,$pool_root/glycine/blind/decisions.tsv" \
    --dataset "brassica=development,$feature_root/brassica/labeled_features.tsv,$feature_root/brassica/topology_features.tsv,$pool_root/brassica/blind/decisions.tsv" \
    --dataset "cotton=retrospective_diagnostic,$feature_root/cotton/labeled_features.tsv,$feature_root/cotton/topology_features.tsv,$pool_root/cotton/blind/decisions.tsv" \
    --dataset "maize=retrospective_diagnostic,$feature_root/maize/labeled_features.tsv,$feature_root/maize/topology_features.tsv,$pool_root/maize/blind/decisions.tsv" \
    --bootstrap-replicates 20000 --bootstrap-seed 20260808 \
    --output-dir "$output"
(cd "$output" && sha256sum -c SHA256SUMS >/dev/null)
printf 'v0.3 ranker evaluation frozen: %s\n' "$output"
