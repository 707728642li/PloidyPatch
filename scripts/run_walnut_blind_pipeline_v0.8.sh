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
    "$staged/blind_role_manifest.json" \
    "$benchmark/perturbed.gff3" "$benchmark/blind_manifest.json" "$benchmark/SHA256SUMS" \
    "$protocol/SHA256SUMS" "$execution/SHA256SUMS"; do
    [[ -s $required ]] || { echo "missing Walnut blind prerequisite: $required" >&2; exit 69; }
done
for forbidden in /nas_data "$staged/evaluator_only" "$staged/target_complete" \
    "$staged/truth" "$staged/labels"; do
    [[ ! -e $forbidden ]] || { echo "forbidden Walnut role visible: $forbidden" >&2; exit 70; }
done
(cd "$benchmark" && sha256sum -c SHA256SUMS >/dev/null)
(cd "$protocol" && sha256sum -c SHA256SUMS >/dev/null)
(cd "$execution" && sha256sum -c SHA256SUMS >/dev/null)
PYTHONPATH="$code_root/src" "$python_bin" - "$protocol" "$execution" <<'PY'
from pathlib import Path
import sys
from ploidypatch.walnut_h1_framework import verify_execution
verify_execution(Path(sys.argv[2]), Path(sys.argv[1]))
PY

command_log=$project_root/pipeline_commands.tsv
working_log=${command_log}.working
[[ ! -e $command_log && ! -L $command_log && ! -e $working_log && ! -L $working_log ]] || exit 71
{
    printf 'stage_order\tentry\timplementation_sha256\n'
    for row in \
        '10 scripts/prepare_walnut_blind_candidate_inputs_v0.8.py' \
        '20 scripts/run_walnut_candidate_methods_v0.8.sh' \
        '30 scripts/build_walnut_h1_candidate_pools_v0.8.py'; do
        order=${row%% *}; relative=${row#* }
        printf '%s\t%s\t%s\n' "$order" "$relative" \
            "$(sha256sum "$code_root/$relative" | awk '{print $1}')"
    done
} > "$working_log"
mv "$working_log" "$command_log"

export PYTHONPATH="$code_root/src"
"$python_bin" "$code_root/scripts/prepare_walnut_blind_candidate_inputs_v0.8.py" "$project_root"
bash "$code_root/scripts/run_walnut_candidate_methods_v0.8.sh" "$project_root"
"$python_bin" "$code_root/scripts/build_walnut_h1_candidate_pools_v0.8.py" "$project_root"

root=$project_root/results/copy_collapse/external/walnut_v0.8_h1
for relative in raw_predictions.manifest.json \
    retain_distinct/blind/candidate.gff3 retain_distinct/blind/decisions.tsv \
    retain_distinct/blind/candidate.gff3.manifest.json \
    suppress_overlap/blind/candidate.gff3 suppress_overlap/blind/decisions.tsv \
    suppress_overlap/blind/candidate.gff3.manifest.json; do
    [[ -s $root/$relative ]] || { echo "missing exact Walnut blind output: $relative" >&2; exit 72; }
done
