#!/usr/bin/env bash
set -euo pipefail
umask 077

[[ $# -eq 1 ]] || { echo "usage: $0 PROJECT_ROOT" >&2; exit 64; }
project_root=$(realpath "$1")
code_root=$project_root/code
python_bin=$project_root/envs/ploidypatch-dev/bin/python
staged=${PLOIDYPATCH_STAGED_INPUT_ROOT:?}
benchmark=${PLOIDYPATCH_BLIND_BENCHMARK_ROOT:?}
protocol=${PLOIDYPATCH_PROTOCOL_FREEZE:?}
execution=${PLOIDYPATCH_EXECUTION_FREEZE:?}
contract=${PLOIDYPATCH_HOLDOUT_CONTRACT:?}
blind_output=${PLOIDYPATCH_BLIND_OUTPUT_ROOT:?}

[[ ${PLOIDYPATCH_BLIND_RUNNER:-} == 1 ]] || exit 65
[[ ${PLOIDYPATCH_NETWORK_ACCESS:-} == none ]] || exit 66
[[ -z ${PLOIDYPATCH_COMPOSITE_MODEL_FREEZE:-} && -z ${PLOIDYPATCH_RANKER_FREEZE:-} ]] || exit 67
[[ $blind_output == /run/blind-run && $project_root == /run/blind-run/project ]] || exit 68
for required in "$python_bin" "$contract" "$staged/role_manifest.tsv" \
    "$staged/role_manifest.json" "$staged/candidate_only/protein_universes/SHA256SUMS" \
    "$benchmark/perturbed.gff3" "$benchmark/blind_manifest.json" \
    "$benchmark/SHA256SUMS" "$protocol/SHA256SUMS" "$execution/SHA256SUMS"; do
    [[ -s $required ]] || { echo "missing Coffea blind prerequisite: $required" >&2; exit 69; }
done
for forbidden in /nas_data "$staged/evaluator_only" "$staged/target_complete" \
    "$staged/truth" "$staged/labels"; do
    [[ ! -e $forbidden ]] || { echo "forbidden Coffea role visible: $forbidden" >&2; exit 70; }
done
for root in "$staged" "$staged/candidate_only/protein_universes" "$benchmark" "$protocol" "$execution"; do
    (cd "$root" && sha256sum -c SHA256SUMS >/dev/null)
done
PYTHONPATH="$code_root/src" "$python_bin" - "$protocol" "$execution" <<'PY'
from pathlib import Path
import sys
from ploidypatch.coffea_h1_framework import verify_execution
verify_execution(Path(sys.argv[2]), Path(sys.argv[1]))
PY

command_log=$project_root/pipeline_commands.tsv
working_log=${command_log}.working
[[ ! -e $command_log && ! -L $command_log && ! -e $working_log && ! -L $working_log ]] || exit 71
{
    printf 'stage_order\tentry\timplementation_sha256\n'
    for row in \
        '10 scripts/prepare_coffea_blind_candidate_inputs_v1.0.py' \
        '20 scripts/run_coffea_candidate_methods_v1.0.sh' \
        '30 scripts/build_coffea_h1_candidate_pools_v1.0.py'; do
        order=${row%% *}; relative=${row#* }
        printf '%s\t%s\t%s\n' "$order" "$relative" \
            "$(sha256sum "$code_root/$relative" | awk '{print $1}')"
    done
} > "$working_log"
mv "$working_log" "$command_log"

export PYTHONPATH="$code_root/src"
"$python_bin" "$code_root/scripts/prepare_coffea_blind_candidate_inputs_v1.0.py" "$project_root"
bash "$code_root/scripts/run_coffea_candidate_methods_v1.0.sh" "$project_root"
"$python_bin" "$code_root/scripts/build_coffea_h1_candidate_pools_v1.0.py" "$project_root"

root=$project_root/results/copy_collapse/external/coffea_v1.0_h1
[[ -s $root/raw_predictions.manifest.json ]] || exit 72
for scope in combined bua_only mauritiana_only; do
    for arm in retain_distinct suppress_overlap; do
        for relative in candidate.gff3 decisions.tsv candidate.gff3.manifest.json; do
            [[ -s $root/$scope/$arm/blind/$relative ]] || {
                echo "missing Coffea blind output: $scope/$arm/$relative" >&2
                exit 72
            }
        done
    done
done
