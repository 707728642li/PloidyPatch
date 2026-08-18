#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then echo "usage: $0 PROJECT_ROOT SPECIES" >&2; exit 2; fi
project_root=$(realpath "$1")
species=$2
code_root=$project_root/code
python_bin=$project_root/envs/ploidypatch-dev/bin/python
model_python=$project_root/envs/ploidypatch-model/bin/python
features=$project_root/results/copy_collapse/candidate_pool_v0.3_features/$species
model=$project_root/results/copy_collapse/model_development/homeolog_ranker_v0.2/model.json
result_root=$project_root/results/copy_collapse/candidate_pool_v0.3_v0.2_scores/$species
working_root=${result_root}.working
for required in "$python_bin" "$model_python" "$features/SHA256SUMS" "$model" \
    "$features/copy_features.tsv" "$features/topology_features.tsv" \
    "$features/labeled_features.tsv"; do
    [[ -s $required ]] || { echo "missing v0.3 score input: $required" >&2; exit 1; }
done
(cd "$features" && sha256sum -c SHA256SUMS >/dev/null)
if [[ -e $result_root || -e $working_root ]]; then
    echo "refusing to overwrite v0.3 v0.2 score: $result_root" >&2; exit 1
fi
mkdir -p "$working_root"
cd "$code_root"
"$python_bin" -m ploidypatch.cli evidence score-homeolog-copy-candidates \
    --copy-features "$features/copy_features.tsv" \
    --topology-features "$features/topology_features.tsv" \
    --model-json "$model" --output-tsv "$working_root/scores.tsv" \
    > "$working_root/scoring.stdout.json" 2> "$working_root/scoring.stderr.log"
env PYTHONPATH="$code_root/src" "$model_python" \
    "$code_root/scripts/evaluate_frozen_homeolog_ranker_v0.2.py" \
    --scores "$working_root/scores.tsv" \
    --labels "$features/labeled_features.tsv" \
    --output-json "$working_root/evaluation.json" \
    --truth-event-count 800 --bootstrap-replicates 20000 --seed 20260808 \
    > "$working_root/evaluation.stdout.json" \
    2> "$working_root/evaluation.stderr.log"
{
    printf 'field\tvalue\ncode_commit\t%s\nspecies\t%s\n' \
        "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}" "$species"
    printf 'model_sha256\t%s\n' "$(sha256sum "$model" | awk '{print $1}')"
    printf 'candidate_policy\tretain_distinct_phased_cds_chains\n'
    printf 'model_policy\tfrozen_v0.2_zero_refit_diagnostic\n'
    printf 'formal_holdout_claim_allowed\tfalse\n'
} > "$working_root/run_contract.tsv"
(
    cd "$working_root"
    find . -type f \( -name '*.json' -o -name '*.tsv' \) -print0 \
        | sort -z | xargs -0 sha256sum > SHA256SUMS
    sha256sum -c SHA256SUMS >/dev/null
)
mv "$working_root" "$result_root"
printf 'v0.3 candidate pool scored by frozen v0.2: %s\n' "$result_root"
