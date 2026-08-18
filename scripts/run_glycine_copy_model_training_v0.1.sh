#!/usr/bin/env bash
set -euo pipefail

project_root=${1:-/data/codexli/projects/PloidyPatch}
code_root=$project_root/code
environment=$project_root/envs/ploidypatch-model
input_root=$project_root/results/copy_collapse/model_development/glycine_feature_matrix_v0.1
result_root=$project_root/results/copy_collapse/model_development/glycine_copy_model_v0.1
temporary_time=$project_root/results/copy_collapse/model_development/.glycine_copy_model_v0.1.resource.time.txt

for required in \
    "$environment/bin/python" \
    "$code_root/scripts/train_copy_model_v0.1.py" \
    "$input_root/evaluator/labeled_features.tsv"; do
    if [[ ! -s "$required" ]]; then
        echo "missing or empty Glycine model input: $required" >&2
        exit 1
    fi
done
if [[ -e "$result_root" || -e "$result_root.partial" || -e "$temporary_time" ]]; then
    echo "refusing to overwrite Glycine model artifact" >&2
    exit 1
fi

export PYTHONPATH=$code_root/src
/usr/bin/time -v -o "$temporary_time" \
    conda run -p "$environment" --no-capture-output \
    python "$code_root/scripts/train_copy_model_v0.1.py" \
        --labeled-features "$input_root/evaluator/labeled_features.tsv" \
        --output-dir "$result_root" \
        --seed 20260807 \
        --outer-folds 5 \
        --inner-folds 4 \
        --c-grid 0.03 0.1 0.3 1 3 \
        --bootstrap-replicates 10000 \
        --high-precision-lower 0.90 \
        --high-confidence 0.95 \
        --high-minimum-selected 30

mv "$temporary_time" "$result_root/resource.time.txt"
conda list -p "$environment" --explicit > "$result_root/environment.explicit.txt"
(
    cd "$result_root"
    sha256sum environment.explicit.txt resource.time.txt >> SHA256SUMS
    sha256sum -c SHA256SUMS
)
printf 'Glycine nested grouped copy model frozen: %s\n' "$result_root"
