#!/usr/bin/env bash
set -euo pipefail
umask 077

[[ $# -eq 5 ]] || {
    echo "usage: $0 PROJECT_ROOT EXECUTION_FREEZE EVALUATOR_INPUT_ROOT EVALUATOR_RUN_ROOT BLIND_BENCHMARK_OUTPUT" >&2
    exit 64
}
project_root=$(realpath "$1")
execution=$(realpath "$2")
input_root=$(realpath "$3")
output=$(realpath -m "$4")
blind_output=$(realpath -m "$5")
working=${output}.working
blind_working=${blind_output}.working
code_root=$execution/source
python=$project_root/envs/ploidypatch-dev/bin/python

[[ -x $python && -s $execution/SHA256SUMS && -s $input_root/SHA256SUMS ]] || exit 65
[[ ! -e $output && ! -L $output && ! -e $working && ! -L $working ]] || exit 66
[[ ! -e $blind_output && ! -L $blind_output && ! -e $blind_working && ! -L $blind_working ]] || exit 67
for path in "$execution" "$input_root" "$output" "$blind_output"; do
    [[ $path != /nas_data && $path != /nas_data/* ]] || exit 68
done
(cd "$execution" && sha256sum -c SHA256SUMS >/dev/null)
(cd "$input_root" && sha256sum -c SHA256SUMS >/dev/null)
mkdir -p "$working"
export PYTHONPATH=$code_root/src
export PLOIDYPATCH_EXECUTION_FREEZE=$execution
export PLOIDYPATCH_EVALUATOR_INPUT_ROOT=$input_root
export PLOIDYPATCH_EVALUATOR_ONLY_ROOT=$working
export PLOIDYPATCH_BLIND_BENCHMARK_OUTPUT=$blind_output

command_log=$working/evaluator_pipeline_commands.tsv
{
    printf 'stage_order\tentry\timplementation_sha256\n'
    for row in \
        '10 scripts/prepare_coffea_evaluator_wgdi_inputs_v1.0.py' \
        '20 scripts/run_coffea_evaluator_wgdi_v1.0.sh' \
        '30 scripts/infer_coffea_external_pairs_v1.0.py' \
        '40 scripts/build_coffea_structure_holdout_v1.0.py'; do
        order=${row%% *}; relative=${row#* }
        printf '%s\t%s\t%s\n' "$order" "$relative" \
            "$(sha256sum "$code_root/$relative" | awk '{print $1}')"
    done
} > "$command_log"

"$python" "$code_root/scripts/prepare_coffea_evaluator_wgdi_inputs_v1.0.py" "$project_root"
bash "$code_root/scripts/run_coffea_evaluator_wgdi_v1.0.sh" "$project_root"
"$python" "$code_root/scripts/infer_coffea_external_pairs_v1.0.py" "$project_root"
"$python" "$code_root/scripts/build_coffea_structure_holdout_v1.0.py" "$project_root"

for required in \
    "$working/wgdi/input/SHA256SUMS" "$working/wgdi/evidence/SHA256SUMS" \
    "$working/truth_pairs/SHA256SUMS" "$working/benchmark/SHA256SUMS" \
    "$blind_output/SHA256SUMS" "$blind_output/blind_manifest.json" \
    "$blind_output/perturbed.gff3"; do
    [[ -s $required ]] || { echo "missing Coffea evaluator output: $required" >&2; exit 69; }
done
for root in "$working/wgdi/input" "$working/wgdi/evidence" \
    "$working/truth_pairs" "$working/benchmark" "$blind_output"; do
    (cd "$root" && sha256sum -c SHA256SUMS >/dev/null)
done
PYTHONPATH="$code_root/src" "$python" - "$working" <<'PY'
from pathlib import Path
import sys
from ploidypatch.artifact_manifest import write_sha256sums, verify_sha256sums
root=Path(sys.argv[1]); write_sha256sums(root); verify_sha256sums(root,ignore_checksum_file=True)
PY
chmod -R a-w "$working" "$blind_output"
mv "$working" "$output"
