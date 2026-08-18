#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 || $# -gt 2 ]]; then
    echo "usage: $0 PROJECT_ROOT [EXECUTION_FREEZE]" >&2; exit 2
fi
project_root=$(realpath "$1")
execution=${2:-$project_root/results/protocol_freezes/populus_external_v0.4_execution}
execution=$(realpath "$execution")
case "$execution" in
    "$project_root/results/protocol_freezes/populus_external_v0.4_execution"*) ;;
    *) echo "execution freeze must be a Populus v0.4 project freeze" >&2; exit 2 ;;
esac
protocol=$project_root/results/protocol_freezes/populus_external_v0.4
model=$project_root/results/models/ploidypatch_ranker_v0.4
blind=$project_root/results/blind_runs/populus_external_v0.4
blind_project=$blind/project
ranking=$blind_project/results/copy_collapse/external/populus_v0.4_blind_rankings
method=$blind_project/results/copy_collapse/external/populus_v0.4_method_trio
scores=$ranking/scores/v04.tsv
score_manifest=$ranking/scores/v04.tsv.manifest.json
decisions=$method/consensus/primary_union/blind/decisions.tsv
pool_manifest=$method/consensus/primary_union/blind/candidate.gff3.manifest.json
custody=$blind/custody_manifest.json
benchmark=$project_root/benchmark/structure/copy_collapse_v0.4/ptr_v4/annotation_copy_collapse_seed20260930
evaluability=$benchmark/evaluator/pair_selection/evaluability.json
reveal_inputs=$project_root/results/evaluator/populus/v0.4/reveal_inputs
result=$project_root/results/copy_collapse/external/populus_v0.4_reveal
working=${result}.working

for required in "$execution/SHA256SUMS" "$execution/implementation_manifest.tsv" \
    "$execution/environment_bindings.tsv" "$protocol/SHA256SUMS" "$protocol/policy.tsv" \
    "$model/SHA256SUMS" "$blind/SHA256SUMS" "$custody" "$scores" \
    "$score_manifest" "$decisions" "$pool_manifest" \
    "$execution/source/scripts/build_populus_complete_control_reveal_inputs_v0.4.sh" \
    "$execution/source/scripts/evaluate_external_v0.4.py"; do
    [[ -s $required ]] || { echo "missing Populus reveal prerequisite: $required" >&2; exit 1; }
done
for root in "$execution" "$protocol" "$model" "$blind"; do
    (cd "$root" && sha256sum -c SHA256SUMS >/dev/null)
done
expected_runner=$(awk -F '\t' '$1 == "scripts/run_populus_external_reveal_v0.4.sh" {print $3}' \
    "$execution/implementation_manifest.tsv")
[[ $expected_runner =~ ^[0-9a-f]{64}$ && $(sha256sum "$0" | awk '{print $1}') == "$expected_runner" ]] || {
    echo "reveal launcher bytes differ from execution freeze" >&2; exit 1;
}
if find "$blind" -perm /222 -print -quit | grep -q .; then
    echo "blind run is not a read-only freeze" >&2; exit 1
fi
[[ ! -e $reveal_inputs && ! -e ${reveal_inputs}.working ]] || {
    echo "refusing a second or partial Populus evaluator reveal-input build" >&2; exit 1;
}
[[ ! -e $result && ! -e $working ]] || {
    echo "refusing to overwrite Populus external reveal" >&2; exit 1;
}

dev_prefix=$(awk -F '\t' '$1 == "ploidypatch-dev" {print $2}' \
    "$execution/environment_bindings.tsv")
python_bin=$dev_prefix/bin/python
[[ -x $python_bin ]] || { echo "frozen reveal Python is absent" >&2; exit 1; }
model_prefix=$(awk -F '\t' '$1 == "ploidypatch-model" {print $2}' \
    "$execution/environment_bindings.tsv")
model_python=$model_prefix/bin/python
[[ -x $model_python ]] || { echo "frozen model/evaluator Python is absent" >&2; exit 1; }
[[ -d $execution/source/src/ploidypatch ]] || { echo "frozen PloidyPatch source is absent" >&2; exit 1; }
export PYTHONPATH=$execution/source/src

# This is the reveal barrier.  It reads only already-sealed blind artifacts and
# exits before the builder is allowed to open truth or complete annotation.
mkdir -p "$working"
"$python_bin" - "$blind" "$custody" "$scores" "$score_manifest" \
    "$decisions" "$pool_manifest" "$working/reveal_authorization.json" <<'PY'
import hashlib
import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

blind, custody_path, scores, score_manifest, decisions, pool_manifest, output = map(
    Path, sys.argv[1:]
)

def sha(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(8 * 1024 * 1024):
            digest.update(block)
    return digest.hexdigest()

custody = json.loads(custody_path.read_text(encoding="utf-8"))
forbidden = (
    "truth_mounted",
    "complete_target_annotation_mounted",
    "evaluator_references_mounted",
    "nas_data_mounted",
    "network_access",
)
blind_outputs = custody.get("blind_outputs", {})
if (
    custody.get("schema_version") != "ploidypatch.blind_run_custody.v1"
    or any(custody.get(field) is not False for field in forbidden)
    or custody.get("runner_identity") != "bubblewrap_populus_v0.4_blind_runner"
    or custody.get("bubblewrap", {}).get("required_flags")
    != ["--unshare-all", "--unshare-net"]
    or blind_outputs.get("scores_sha256") != sha(scores)
    or blind_outputs.get("score_manifest_sha256") != sha(score_manifest)
    or blind_outputs.get("pool_decisions_sha256") != sha(decisions)
    or blind_outputs.get("pool_manifest_sha256") != sha(pool_manifest)
):
    raise SystemExit("blind custody barrier failed; truth reveal forbidden")
authorization = {
    "schema_version": "ploidypatch.populus_reveal_authorization.v0.4",
    "authorized_at": datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z"),
    "truth_opened": False,
    "blind_run_SHA256SUMS_sha256": sha(blind / "SHA256SUMS"),
    "custody_manifest_sha256": sha(custody_path),
    "blind_scores_sha256": sha(scores),
    "blind_score_manifest_sha256": sha(score_manifest),
    "pool_decisions_sha256": sha(decisions),
    "pool_manifest_sha256": sha(pool_manifest),
}
flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL
descriptor = os.open(output, flags, 0o440)
with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
    json.dump(authorization, handle, indent=2, sort_keys=True)
    handle.write("\n")
PY

# Only after the barrier has been durably written may evaluator-owned code see
# complete annotation/truth and build controls and exact-chain labels.
set +e
PLOIDYPATCH_BLIND_RUN_ROOT="$blind" \
PLOIDYPATCH_REVEAL_AUTHORIZATION="$working/reveal_authorization.json" \
PLOIDYPATCH_EXECUTION_FREEZE_OVERRIDE="$execution" \
    bash "$execution/source/scripts/build_populus_complete_control_reveal_inputs_v0.4.sh" \
    "$project_root"
builder_status=$?
set -e
if [[ $builder_status -ne 0 && ! -s $reveal_inputs/status.json ]]; then
    echo "Populus reveal-input builder failed before freezing a scientific status" >&2
    exit "$builder_status"
fi
for required in "$reveal_inputs/SHA256SUMS" "$reveal_inputs/status.json" \
    "$reveal_inputs/reveal_input_manifest.json"; do
    [[ -s $required ]] || { echo "incomplete evaluator reveal inputs: $required" >&2; exit 1; }
done
(cd "$reveal_inputs" && sha256sum -c SHA256SUMS >/dev/null)
while IFS= read -r -d '' path; do chmod 0440 "$path"; done < <(find "$reveal_inputs" -type f -print0)
while IFS= read -r -d '' path; do chmod 0550 "$path"; done < <(find "$reveal_inputs" -depth -type d -print0)
reveal_status=$("$python_bin" - "$reveal_inputs/status.json" "$reveal_inputs/reveal_input_manifest.json" \
    "$blind" "$custody" "$scores" "$score_manifest" "$decisions" "$pool_manifest" <<'PY'
import hashlib
import json
import sys
from pathlib import Path

status_path, manifest_path, blind, custody, scores, score_manifest, decisions, pool_manifest = map(Path, sys.argv[1:])
def sha(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(8 * 1024 * 1024): digest.update(block)
    return digest.hexdigest()
status = json.loads(status_path.read_text(encoding="utf-8"))
manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
expected = {
    "blind_run_SHA256SUMS_sha256": sha(blind / "SHA256SUMS"),
    "custody_manifest_sha256": sha(custody),
    "blind_scores_sha256": sha(scores),
    "blind_score_manifest_sha256": sha(score_manifest),
    "pool_decisions_sha256": sha(decisions),
    "pool_manifest_sha256": sha(pool_manifest),
}
if (
    manifest.get("schema_version") != "ploidypatch.populus_reveal_inputs.v0.4"
    or manifest.get("generated_after_blind_freeze") is not True
    or manifest.get("evaluator_only") is not True
    or any(manifest.get(key) != value for key, value in expected.items())
):
    raise SystemExit("reveal inputs do not bind the frozen blind run")
value = status.get("status")
if value not in {
    "ready_for_evaluation",
    "not_evaluable_without_rule_relaxation",
    "invalid_run",
}:
    raise SystemExit("unknown Populus reveal-input status")
print(value)
PY
)
if [[ $builder_status -ne 0 && $reveal_status != invalid_run ]]; then
    echo "failed reveal-input builder did not freeze invalid_run status" >&2
    exit "$builder_status"
fi

if [[ $reveal_status != ready_for_evaluation ]]; then
    mkdir -p "$working/evaluation"
    "$python_bin" - "$reveal_status" "$reveal_inputs/status.json" \
        "$reveal_inputs/reveal_input_manifest.json" \
        "$working/reveal_authorization.json" "$protocol/SHA256SUMS" \
        "$model/SHA256SUMS" "$blind/SHA256SUMS" "$custody" \
        "$working/evaluation/evaluation.json" <<'PY'
import hashlib
import json
import os
import sys
from pathlib import Path

(
    formal_outcome, status_path, manifest_path, authorization_path,
    protocol_sums, model_sums, blind_sums, custody_path, output,
) = sys.argv[1:]
def sha(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        while block := handle.read(8 * 1024 * 1024): digest.update(block)
    return digest.hexdigest()
status = json.loads(Path(status_path).read_text(encoding="utf-8"))
report = {
    "schema_version": "ploidypatch.external_ranker_evaluation.v0.4",
    "evaluation_role": "untouched_confirmatory_external_species",
    "formal_outcome": formal_outcome,
    "confirmatory_pass": False,
    "reason": status.get("reason", "fixed preregistered data or sentinel gate did not permit ranking evaluation"),
    "evaluator_invoked": False,
    "fixed_rules_relaxed": False,
    "inputs": {
        "status_sha256": sha(status_path),
        "reveal_input_manifest_sha256": sha(manifest_path),
        "reveal_authorization_sha256": sha(authorization_path),
        "protocol_SHA256SUMS_sha256": sha(protocol_sums),
        "composite_model_SHA256SUMS_sha256": sha(model_sums),
        "blind_run_SHA256SUMS_sha256": sha(blind_sums),
        "custody_manifest_sha256": sha(custody_path),
    },
}
flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL
descriptor = os.open(output, flags, 0o440)
with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
    json.dump(report, handle, indent=2, sort_keys=True)
    handle.write("\n")
PY
    (
        cd "$working"
        find . -type f ! -name SHA256SUMS -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS
        sha256sum -c SHA256SUMS >/dev/null
    )
    while IFS= read -r -d '' path; do chmod 0440 "$path"; done < <(find "$working" -type f -print0)
    while IFS= read -r -d '' path; do chmod 0550 "$path"; done < <(find "$working" -depth -type d -print0)
    if [[ $reveal_status == invalid_run ]]; then
        echo "Populus evaluator sentinel/custody run is invalid; retained at $working" >&2
        exit 1
    fi
    mv "$working" "$result"
    printf 'Populus v0.4 frozen outcome: not_evaluable_without_rule_relaxation (%s)\n' "$result"
    exit 0
fi

for required in "$reveal_inputs/labels/candidate_labels.tsv" \
    "$reveal_inputs/labels/candidate_labels.tsv.manifest.json" "$evaluability" \
    "$reveal_inputs/scores/consensus/primary_union.json" \
    "$reveal_inputs/scores/consensus/legacy_union.json" \
    "$reveal_inputs/scores/consensus/support2.json" \
    "$reveal_inputs/scores/consensus/support3.json" \
    "$reveal_inputs/scores/methods/miniprot.json" \
    "$reveal_inputs/scores/methods/gemoma.json" \
    "$reveal_inputs/scores/methods/lifton.json"; do
    [[ -s $required ]] || { echo "ready reveal lacks evaluator input: $required" >&2; exit 1; }
done

secondary=(
    --secondary-score "miniprot=$reveal_inputs/scores/methods/miniprot.json"
    --secondary-score "gemoma=$reveal_inputs/scores/methods/gemoma.json"
    --secondary-score "lifton=$reveal_inputs/scores/methods/lifton.json"
    --secondary-score "support2=$reveal_inputs/scores/consensus/support2.json"
    --secondary-score "support3=$reveal_inputs/scores/consensus/support3.json"
)
"$model_python" "$execution/source/scripts/evaluate_external_v0.4.py" \
    --scores "$scores" \
    --labels "$reveal_inputs/labels/candidate_labels.tsv" \
    --pool-decisions "$decisions" --pool-manifest "$pool_manifest" \
    --primary-pool-score "$reveal_inputs/scores/consensus/primary_union.json" \
    --legacy-pool-score "$reveal_inputs/scores/consensus/legacy_union.json" \
    --evaluability "$evaluability" --custody-manifest "$custody" \
    --protocol-freeze "$protocol" --composite-model-freeze "$model" \
    --policy "$protocol/policy.tsv" "${secondary[@]}" \
    --output-dir "$working/evaluation"
(
    cd "$working"
    find . -type f ! -name SHA256SUMS -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS
    sha256sum -c SHA256SUMS >/dev/null
)
while IFS= read -r -d '' path; do chmod 0440 "$path"; done < <(find "$working" -type f -print0)
while IFS= read -r -d '' path; do chmod 0550 "$path"; done < <(find "$working" -depth -type d -print0)
mv "$working" "$result"
printf 'Populus v0.4 external result revealed and frozen: %s\n' "$result"
