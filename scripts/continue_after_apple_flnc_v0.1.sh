#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "usage: $0 PROJECT_ROOT APPLE_PIPELINE_PID" >&2
    exit 2
fi
project_root=$(realpath "$1")
apple_pid=$2
[[ $apple_pid =~ ^[1-9][0-9]*$ ]] || {
    echo "APPLE_PIPELINE_PID must be a positive integer" >&2; exit 2;
}

code_root=$project_root/code
dev_python=$project_root/envs/ploidypatch-dev/bin/python
publication_python=$project_root/envs/ploidypatch-publication/bin/python
conda_bin=/data/codexli/software/conda/miniforge3/bin/conda
apple_result=$project_root/results/natural/apple_gddh13_v0.4/validation/ont_flnc_pychopper_v0.1
verification=$project_root/results/verification/apple_flnc_v0.1/full_verification.json
manuscript_source=$project_root/results/manuscript_source/apple_flnc_v0.1
reference=$project_root/results/scaling/core_v0.1_references_3a32ba5
benchmark=$project_root/results/scaling/core_v0.1_benchmark_3a32ba5
summary=$project_root/results/scaling/core_v0.1_summary_3a32ba5
figure=$project_root/results/scaling/core_v0.1_figure_3a32ba5
publication_environment=$project_root/results/environments/ploidypatch_publication_v0.1
frozen_code=$reference/frozen_code

for required in "$dev_python" "$publication_python" "$conda_bin" \
    "$code_root/scripts/verify_apple_flnc_result_v0.1.py" \
    "$code_root/scripts/build_apple_flnc_manuscript_tables_v0.1.py" \
    "$code_root/scripts/freeze_publication_environment_v0.1.py" \
    "$reference/SHA256SUMS" \
    "$frozen_code/scripts/run_core_scaling_benchmark_v0.1.py" \
    "$frozen_code/scripts/summarize_core_scaling_benchmark_v0.1.py" \
    "$frozen_code/scripts/plot_core_scaling_benchmark_v0.1.py"; do
    [[ -s $required ]] || { echo "missing continuation input: $required" >&2; exit 1; }
done
[[ $(sha256sum "$code_root/scripts/verify_apple_flnc_result_v0.1.py" | awk '{print $1}') == \
    84a0431dab4202748749707df60b67287411390a5042296fae57322c87309ef4 ]] || {
    echo "apple verifier differs from audited implementation" >&2; exit 1;
}
[[ $(sha256sum "$code_root/scripts/build_apple_flnc_manuscript_tables_v0.1.py" | awk '{print $1}') == \
    1e2e6a624ccff304abc7007986c0874e05ba45e7fb5f3ed48be9887baa8ef826 ]] || {
    echo "apple source-table builder differs from audited implementation" >&2; exit 1;
}
[[ ! -e ${verification}.working ]] || {
    echo "refusing stale verification working file: ${verification}.working" >&2; exit 1;
}
[[ ! -e $verification || -s $verification ]] || {
    echo "existing verification is empty: $verification" >&2; exit 1;
}
for target in "$manuscript_source" "$benchmark" "$summary" \
              "$figure" "$publication_environment"; do
    [[ ! -e $target && ! -e ${target}.working ]] || {
        echo "refusing continuation overwrite: $target" >&2; exit 1;
    }
done
(cd "$reference" && sha256sum -c SHA256SUMS >/dev/null)

printf 'waiting_for_apple_pid\t%s\n' "$apple_pid"
while kill -0 "$apple_pid" 2>/dev/null; do
    sleep 30
done
[[ -d $apple_result && ! -e ${apple_result}.working ]] || {
    echo "apple pipeline exited without a finalized result" >&2; exit 1;
}

verification_candidate=$verification
verification_created=0
if [[ ! -e $verification ]]; then
    "$dev_python" "$code_root/scripts/verify_apple_flnc_result_v0.1.py" \
        "$apple_result" --project-root "$project_root" \
        --expected-code-commit c1be36270808dc5dab89f337e151c5a4adc3d3fd \
        --output-json "${verification}.working"
    verification_candidate=${verification}.working
    verification_created=1
else
    printf 'reusing_full_apple_verification\t%s\n' "$verification"
fi
"$dev_python" - "$verification_candidate" "$apple_result" <<'PY'
import hashlib
import json
import pathlib
import sys

report = json.loads(pathlib.Path(sys.argv[1]).read_text(encoding="utf-8"))
if report.get("valid") is not True or report.get("verification_mode") != "full":
    raise SystemExit("apple verification is not a valid full verification")
result = pathlib.Path(sys.argv[2]).resolve()
if pathlib.Path(report.get("result_root", "")).resolve() != result:
    raise SystemExit("apple verification refers to another result")
digest = hashlib.sha256((result / "SHA256SUMS").read_bytes()).hexdigest()
if report.get("sha256sums_sha256") != digest:
    raise SystemExit("apple verification/result freeze mismatch")
PY
if [[ $verification_created -eq 1 ]]; then
    mv -- "${verification}.working" "$verification"
fi
"$dev_python" "$code_root/scripts/build_apple_flnc_manuscript_tables_v0.1.py" \
    "$apple_result" "$verification" "$manuscript_source"

"$dev_python" "$frozen_code/scripts/run_core_scaling_benchmark_v0.1.py" \
    --project-root "$project_root" --python-bin "$dev_python" \
    --config "$frozen_code/config/core_scaling_workloads_v0.1.json" \
    --reference-root "$reference" \
    --code-commit 3a32ba53f0fb85f894568a17063e35e9fa8f1597 \
    --output-dir "$benchmark" --conda "$conda_bin"
"$dev_python" "$frozen_code/scripts/summarize_core_scaling_benchmark_v0.1.py" \
    --benchmark-root "$benchmark" --reference-root "$reference" \
    --output-dir "$summary"

"$dev_python" "$code_root/scripts/freeze_publication_environment_v0.1.py" \
    --environment-prefix "$project_root/envs/ploidypatch-publication" \
    --conda "$conda_bin" --output-dir "$publication_environment"
"$publication_python" "$frozen_code/scripts/plot_core_scaling_benchmark_v0.1.py" \
    --summary-root "$summary" --output-dir "$figure"
printf 'post_apple_continuation_complete\n'
