#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then echo "usage: $0 PROJECT_ROOT" >&2; exit 2; fi
project_root=$(realpath "$1")
code_root=$project_root/code
model_env=$project_root/envs/ploidypatch-model
feature_root=$project_root/results/copy_collapse/model_development/homeolog_topology_v0.2
output_root=$project_root/results/copy_collapse/model_development/homeolog_ranker_v0.2

[[ -x $model_env/bin/python ]] || { echo "missing model environment" >&2; exit 1; }
[[ -s $feature_root/SHA256SUMS ]] || { echo "unfrozen topology features" >&2; exit 1; }
(cd "$feature_root" && sha256sum -c SHA256SUMS >/dev/null)
[[ ! -e $output_root && ! -e ${output_root}.partial ]] || {
    echo "refusing to overwrite frozen homeolog ranker" >&2
    exit 1
}

gly_labels=$project_root/results/copy_collapse/model_development/glycine_feature_matrix_v0.1/evaluator/labeled_features.tsv
bra_labels=$project_root/results/copy_collapse/model_development/brassica_glycine_model_transfer_v0.1/evaluator/labeled_features.tsv

cd "$code_root"
/usr/bin/time -v -o "$project_root/logs/copy_collapse/homeolog_ranker_freeze_v0.2.resource.time.txt" \
    env PYTHONPATH="$code_root/src" \
    PLOIDYPATCH_CODE_COMMIT="${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}" \
    conda run -p "$model_env" --no-capture-output python \
    scripts/train_homeolog_ranker_v0.2.py \
    --dataset "glycine=$gly_labels,$feature_root/glycine/features.tsv" \
    --dataset "brassica=$bra_labels,$feature_root/brassica/features.tsv" \
    --output-dir "$output_root"

(cd "$output_root" && sha256sum -c SHA256SUMS >/dev/null)
printf 'homeolog ranker frozen: %s\n' "$output_root"
