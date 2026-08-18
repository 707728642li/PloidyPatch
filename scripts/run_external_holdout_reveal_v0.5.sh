#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 7 ]]; then
    echo "usage: $0 PROJECT_ROOT EXECUTION_FREEZE PROTOCOL_FREEZE MODEL_FREEZE BLIND_RUN EVALUATOR_ONLY_ROOT NEW_REVEAL_RESULT_DIR" >&2
    exit 2
fi

project_root=$(realpath "$1")
execution=$(realpath "$2")
protocol=$(realpath "$3")
model=$(realpath "$4")
blind_run=$(realpath "$5")
# Deliberately retain this as an opaque string until custody is frozen and
# validated.  No evaluator-owned path is stat'ed or resolved before the gate.
evaluator_only_argument=$6
output=$(realpath -m "$7")
working=${output}.working
allowed_root=$project_root/results/evaluator
case "$output/" in
    "$allowed_root"/*/) ;;
    *) echo "reveal result must be a new child of $allowed_root" >&2; exit 2 ;;
esac
for path in "$execution" "$protocol" "$model" "$blind_run"; do
    [[ $path != /nas_data && $path != /nas_data/* ]] || {
        echo "frozen reveal prerequisite may not be read from /nas_data: $path" >&2
        exit 2
    }
done
[[ ! -e $output && ! -e $working ]] || {
    echo "refusing to overwrite reveal result or retained failed attempt" >&2
    exit 1
}
for required in \
    "$execution/SHA256SUMS" "$execution/execution_manifest.json" \
    "$execution/environment_bindings.tsv" "$execution/source" \
    "$protocol/SHA256SUMS" "$protocol/protocol_manifest.json" \
    "$protocol/contract.json" "$model/SHA256SUMS" \
    "$blind_run/SHA256SUMS" "$blind_run/custody_manifest.json" \
    "$blind_run/project"; do
    [[ -s $required || -d $required ]] || {
        echo "missing frozen reveal prerequisite: $required" >&2; exit 1;
    }
done

command -v conda >/dev/null || { echo "conda is required" >&2; exit 1; }
mkdir -p "$allowed_root" "$working/environment_checks"

dev_python=
model_python=
environment_count=0
while IFS=$'\t' read -r name prefix explicit_relative explicit_sha pip_relative pip_sha; do
    [[ $name == name ]] && continue
    [[ $name =~ ^[a-z0-9][a-z0-9_.-]*$ && -d $prefix ]] || {
        echo "invalid frozen reveal environment binding: $name=$prefix" >&2; exit 1;
    }
    conda list --explicit -p "$prefix" > "$working/environment_checks/$name.explicit.txt"
    if [[ ! -x $prefix/bin/python ]]; then
        printf '# python unavailable; pip lock not applicable\n' \
            > "$working/environment_checks/$name.pip-freeze.txt"
    elif ! "$prefix/bin/python" -m pip freeze --all \
        > "$working/environment_checks/$name.pip-freeze.txt" 2>/dev/null; then
        printf '# pip unavailable; explicit conda lock is authoritative\n' \
            > "$working/environment_checks/$name.pip-freeze.txt"
    fi
    [[ $(sha256sum "$working/environment_checks/$name.explicit.txt" | awk '{print $1}') == "$explicit_sha" \
       && $(sha256sum "$working/environment_checks/$name.pip-freeze.txt" | awk '{print $1}') == "$pip_sha" \
       && $(sha256sum "$execution/$explicit_relative" | awk '{print $1}') == "$explicit_sha" \
       && $(sha256sum "$execution/$pip_relative" | awk '{print $1}') == "$pip_sha" ]] || {
        echo "reveal environment differs from execution freeze: $name" >&2; exit 1;
    }
    [[ $name != ploidypatch-dev ]] || dev_python=$prefix/bin/python
    [[ $name != ploidypatch-model ]] || model_python=$prefix/bin/python
    environment_count=$((environment_count + 1))
done < "$execution/environment_bindings.tsv"
[[ $environment_count -eq 7 && -x $dev_python && -x $model_python ]] || {
    echo "reveal requires seven frozen environments including dev and model" >&2
    exit 1
}

# This is the hard truth-access barrier.  The evaluator-only argument remains
# unresolved until every frozen tree, custody field and blind output digest has
# been validated against an exact SHA256SUMS universe.
PYTHONPATH="$execution/source/src" "$dev_python" - \
    "$execution" "$protocol" "$model" "$blind_run" <<'PY'
from pathlib import Path
import json
import os
import sys

from ploidypatch.artifact_manifest import sha256_file, verify_sha256sums
from ploidypatch.holdout_contract import load_holdout_contract, safe_relative_path

execution, protocol, model, blind = map(Path, sys.argv[1:])
for root in (execution, protocol, model, blind):
    verify_sha256sums(root, ignore_checksum_file=True)
e = json.loads((execution / "execution_manifest.json").read_text(encoding="utf-8"))
p = json.loads((protocol / "protocol_manifest.json").read_text(encoding="utf-8"))
c_path = blind / "custody_manifest.json"
c = json.loads(c_path.read_text(encoding="utf-8"))
contract = load_holdout_contract(protocol / "contract.json")
patch = e.get("execution_patch")
if patch is None:
    patch_valid = (
        e.get("freeze_stage") == "post_metadata_pre_pair_pre_candidate_pre_label"
        and e.get("code_commit") == p.get("code_commit")
    )
else:
    changed = patch.get("changed_files") if isinstance(patch, dict) else None
    failed = patch.get("failed_attempt") if isinstance(patch, dict) else None
    patch_valid = (
        isinstance(patch, dict)
        and patch.get("schema_version") == "ploidypatch.external_holdout_execution_patch.v0.5"
        and patch.get("freeze_stage") == "post_evaluator_truth_failed_blind_pre_candidate_pre_score_pre_label_execution_patch"
        and e.get("freeze_stage") == patch.get("freeze_stage")
        and patch.get("base_code_commit") == p.get("code_commit")
        and patch.get("patch_code_commit") == e.get("code_commit")
        and patch.get("base_protocol_SHA256SUMS_sha256") == sha256_file(protocol / "SHA256SUMS")
        and isinstance(patch.get("superseded_execution_SHA256SUMS_sha256"), str)
        and len(patch["superseded_execution_SHA256SUMS_sha256"]) == 64
        and patch.get("contract_sha256") == sha256_file(protocol / "contract.json")
        and patch.get("composite_model_SHA256SUMS_sha256") == sha256_file(model / "SHA256SUMS")
        and all(
            patch.get(field) is False
            for field in (
                "scientific_protocol_changed", "contract_or_policy_changed",
                "model_or_threshold_changed", "staged_inputs_changed",
                "truth_or_benchmark_regenerated",
                "candidate_generation_completed_before_patch",
                "formal_scores_generated_before_patch",
                "truth_labels_accessed_before_patch", "automatic_approval",
            )
        )
        and patch.get("evaluator_truth_construction_completed_before_patch") is True
        and patch.get("blind_candidate_wgd_completed_before_patch") is False
        and isinstance(changed, list) and bool(changed)
        and isinstance(failed, dict)
        and isinstance(failed.get("exit_status"), int) and failed["exit_status"] != 0
        and (execution / "superseded_failed_attempt_manifest.tsv").is_file()
        and (execution / "patch_reason.md").is_file()
    )
expected_custody_patch = (
    {
        "active": True,
        "freeze_stage": patch["freeze_stage"],
        "base_code_commit": patch["base_code_commit"],
        "patch_code_commit": patch["patch_code_commit"],
        "superseded_execution_SHA256SUMS_sha256": patch[
            "superseded_execution_SHA256SUMS_sha256"
        ],
        "failed_attempt_tree_sha256": patch["failed_attempt"]["tree_sha256"],
    }
    if patch_valid and isinstance(patch, dict)
    else {"active": False}
)
if (
    e.get("schema_version") != "ploidypatch.external_holdout_execution_freeze.v0.5"
    or p.get("schema_version") != "ploidypatch.external_holdout_protocol_freeze.v0.5"
    or c.get("schema_version") != "ploidypatch.external_holdout_blind_custody.v0.5"
    or any(value.get("holdout_id") != contract.holdout_id for value in (e, p, c))
    or any(value.get("policy_id") != contract.policy_id for value in (p, c))
    or any(value.get("model_version") != contract.model_version for value in (p, c))
    or e.get("contract_sha256") != sha256_file(protocol / "contract.json")
    or e.get("protocol_SHA256SUMS_sha256") != sha256_file(protocol / "SHA256SUMS")
    or e.get("composite_model_SHA256SUMS_sha256") != sha256_file(model / "SHA256SUMS")
    or not patch_valid
    or c.get("execution_patch") != expected_custody_patch
):
    raise SystemExit("execution/protocol/model/custody bindings differ")
for field in (
    "truth_mounted", "complete_target_annotation_mounted",
    "evaluator_references_mounted", "nas_data_mounted", "network_access",
):
    if c.get(field) is not False:
        raise SystemExit(f"custody negative-access field is not false: {field}")
frozen = c.get("frozen_inputs", {})
if (
    frozen.get("execution_SHA256SUMS_sha256") != sha256_file(execution / "SHA256SUMS")
    or frozen.get("protocol_SHA256SUMS_sha256") != sha256_file(protocol / "SHA256SUMS")
    or frozen.get("composite_model_SHA256SUMS_sha256") != sha256_file(model / "SHA256SUMS")
    or not c.get("frozen_before_truth_reveal_at")
):
    raise SystemExit("custody does not bind the frozen execution inputs")
outputs = c.get("blind_outputs")
required = {"scores", "score_manifest", "pool_decisions", "pool_manifest", "command_log"}
if not isinstance(outputs, dict) or not required <= set(outputs):
    raise SystemExit("custody lacks exact blind output bindings")
project = blind / "project"
for name in required:
    item = outputs.get(name)
    if not isinstance(item, dict) or set(item) != {"relative_path", "sha256"}:
        raise SystemExit(f"malformed custody output: {name}")
    relative = safe_relative_path(item["relative_path"], f"custody output {name}")
    path = project.joinpath(*relative.parts)
    try:
        path.resolve().relative_to(project.resolve())
    except ValueError:
        raise SystemExit(f"custody output escapes blind project: {name}")
    if not path.is_file() or path.is_symlink() or sha256_file(path) != item["sha256"]:
        raise SystemExit(f"custody output digest differs: {name}")
if not c.get("mount_manifest", {}).get("sha256") or not c.get("commands"):
    raise SystemExit("custody lacks mount or command evidence")
for path in blind.rglob("*"):
    if path.is_file() and os.stat(path).st_mode & 0o222:
        raise SystemExit(f"blind run is writable after custody freeze: {path}")
PY

authorization=$working/reveal_authorization.json
PYTHONPATH="$execution/source/src" "$dev_python" - \
    "$execution" "$protocol" "$model" "$blind_run" "$authorization" <<'PY'
from datetime import datetime, timezone
from pathlib import Path
import json
import os
import sys
from ploidypatch.artifact_manifest import sha256_file
from ploidypatch.holdout_contract import load_holdout_contract

execution, protocol, model, blind, output = map(Path, sys.argv[1:])
contract = load_holdout_contract(protocol / "contract.json")
custody_path = blind / "custody_manifest.json"
custody = json.loads(custody_path.read_text(encoding="utf-8"))
payload = {
    "schema_version": "ploidypatch.external_holdout_reveal_authorization.v0.5",
    "holdout_id": contract.holdout_id,
    "policy_id": contract.policy_id,
    "model_version": contract.model_version,
    "authorized_at": datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z"),
    "truth_opened": False,
    "truth_reveal_authorized": True,
    "authorized_statuses": [
        "ready_for_evaluation",
        "not_evaluable_without_rule_relaxation",
        "invalid_run",
    ],
    "custody_manifest_sha256": sha256_file(custody_path),
    "blind_run_SHA256SUMS_sha256": sha256_file(blind / "SHA256SUMS"),
    "execution_SHA256SUMS_sha256": sha256_file(execution / "SHA256SUMS"),
    "protocol_SHA256SUMS_sha256": sha256_file(protocol / "SHA256SUMS"),
    "composite_model_SHA256SUMS_sha256": sha256_file(model / "SHA256SUMS"),
    "blind_outputs": {
        name: item["sha256"]
        for name, item in custody["blind_outputs"].items()
        if isinstance(item, dict) and "sha256" in item
    },
}
descriptor = os.open(output, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o440)
with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as handle:
    json.dump(payload, handle, indent=2, sort_keys=True)
    handle.write("\n")
PY

# Custody is now immutable and authorization exists.  Only now may the reveal
# process resolve, validate or expose evaluator-owned bytes to the builder.
evaluator_only=$(realpath "$evaluator_only_argument")
[[ -d $evaluator_only && ! -L $evaluator_only ]] || {
    echo "evaluator-only reveal root is missing or symlinked" >&2; exit 1;
}
execution_manifest=$execution/execution_manifest.json
readarray -t pipeline_values < <(PYTHONPATH="$execution/source/src" "$dev_python" - \
    "$execution_manifest" <<'PY'
import json, sys
value = json.load(open(sys.argv[1], encoding="utf-8"))
print(value["pipeline_entries"]["reveal_input_builder"])
print(value["pipeline_entries"]["evaluator"])
print(value["holdout_id"])
print(value["policy_id"])
PY
)
builder_relative=${pipeline_values[0]}
evaluator_relative=${pipeline_values[1]}
holdout_id=${pipeline_values[2]}
policy_id=${pipeline_values[3]}
builder=$execution/source/$builder_relative
evaluator=$execution/source/$evaluator_relative
[[ -f $builder && -f $evaluator ]] || {
    echo "frozen reveal builder or evaluator is missing" >&2; exit 1;
}

reveal_inputs=$working/reveal_inputs
printf '%q ' env \
    PLOIDYPATCH_HOLDOUT_CONTRACT="$protocol/contract.json" \
    PLOIDYPATCH_PROTOCOL_FREEZE="$protocol" \
    PLOIDYPATCH_EXECUTION_FREEZE="$execution" \
    PLOIDYPATCH_COMPOSITE_MODEL_FREEZE="$model" \
    PLOIDYPATCH_BLIND_RUN_ROOT="$blind_run" \
    PLOIDYPATCH_CUSTODY_MANIFEST="$blind_run/custody_manifest.json" \
    PLOIDYPATCH_REVEAL_AUTHORIZATION="$authorization" \
    PLOIDYPATCH_EVALUATOR_ONLY_ROOT="$evaluator_only" \
    PLOIDYPATCH_REVEAL_INPUTS_OUTPUT="$reveal_inputs" \
    /usr/bin/bash "$builder" "$project_root" > "$working/reveal_builder_command.txt"
printf '\n' >> "$working/reveal_builder_command.txt"

set +e
env \
    PLOIDYPATCH_HOLDOUT_CONTRACT="$protocol/contract.json" \
    PLOIDYPATCH_PROTOCOL_FREEZE="$protocol" \
    PLOIDYPATCH_EXECUTION_FREEZE="$execution" \
    PLOIDYPATCH_EXECUTION_FREEZE_OVERRIDE="$execution" \
    PLOIDYPATCH_COMPOSITE_MODEL_FREEZE="$model" \
    PLOIDYPATCH_BLIND_RUN_ROOT="$blind_run" \
    PLOIDYPATCH_CUSTODY_MANIFEST="$blind_run/custody_manifest.json" \
    PLOIDYPATCH_REVEAL_AUTHORIZATION="$authorization" \
    PLOIDYPATCH_EVALUATOR_ONLY_ROOT="$evaluator_only" \
    PLOIDYPATCH_REVEAL_INPUTS_OUTPUT="$reveal_inputs" \
    PYTHONPATH="$execution/source/src" \
    /usr/bin/bash "$builder" "$project_root" \
    > "$working/reveal_builder.stdout.log" \
    2> "$working/reveal_builder.stderr.log"
builder_status=$?
set -e
[[ -d $reveal_inputs && -s $reveal_inputs/SHA256SUMS \
    && -s $reveal_inputs/status.json \
    && -s $reveal_inputs/reveal_input_manifest.json ]] || {
    echo "reveal builder failed without a sealed formal status" >&2; exit 1;
}
PYTHONPATH="$execution/source/src" "$dev_python" - \
    "$reveal_inputs" "$authorization" "$blind_run/custody_manifest.json" \
    "$blind_run/SHA256SUMS" "$holdout_id" "$policy_id" \
    > "$working/formal_status.txt" <<'PY'
from pathlib import Path
import json
import sys
from ploidypatch.artifact_manifest import sha256_file, verify_sha256sums

root, authorization, custody, blind_sums = map(Path, sys.argv[1:5])
holdout, policy = sys.argv[5:]
verify_sha256sums(root, ignore_checksum_file=True)
status_record = json.loads((root / "status.json").read_text(encoding="utf-8"))
manifest = json.loads((root / "reveal_input_manifest.json").read_text(encoding="utf-8"))
allowed = {
    "ready_for_evaluation",
    "not_evaluable_without_rule_relaxation",
    "invalid_run",
}
status = status_record.get("formal_status", status_record.get("status"))
if (
    status not in allowed
    or manifest.get("formal_status", manifest.get("status")) != status
    or status_record.get("holdout_id") != holdout
    or manifest.get("holdout_id") != holdout
    or status_record.get("policy_id") != policy
    or manifest.get("policy_id") != policy
    or manifest.get("custody_manifest_sha256") != sha256_file(custody)
    or manifest.get("reveal_authorization_sha256") != sha256_file(authorization)
    or manifest.get("blind_run_SHA256SUMS_sha256") != sha256_file(blind_sums)
    or manifest.get("generated_after_blind_freeze") is not True
    or manifest.get("evaluator_only") is not True
):
    raise SystemExit("reveal status or frozen bindings differ")
print(status)
PY
formal_status=$(tr -d '\r\n' < "$working/formal_status.txt")
if [[ $formal_status == invalid_run ]]; then
    [[ $builder_status -ne 0 ]] || {
        echo "invalid reveal status must be returned as failure" >&2; exit 1;
    }
elif [[ $builder_status -ne 0 ]]; then
    echo "reveal builder failed for a non-invalid formal status" >&2
    exit "$builder_status"
fi

mkdir "$working/evaluation"
evaluator_invoked=false
if [[ $formal_status == ready_for_evaluation ]]; then
    evaluator_invoked=true
    rmdir "$working/evaluation"
    printf '%q ' "$model_python" "$evaluator" \
        --execution-freeze "$execution" \
        --protocol-freeze "$protocol" \
        --composite-model-freeze "$model" \
        --blind-run "$blind_run/project" \
        --custody-manifest "$blind_run/custody_manifest.json" \
        --reveal-inputs "$reveal_inputs" \
        --output-dir "$working/evaluation" \
        > "$working/evaluator_command.txt"
    printf '\n' >> "$working/evaluator_command.txt"
    set +e
    PYTHONPATH="$execution/source/src" "$model_python" "$evaluator" \
        --execution-freeze "$execution" \
        --protocol-freeze "$protocol" \
        --composite-model-freeze "$model" \
        --blind-run "$blind_run/project" \
        --custody-manifest "$blind_run/custody_manifest.json" \
        --reveal-inputs "$reveal_inputs" \
        --output-dir "$working/evaluation" \
        > "$working/evaluator.stdout.log" \
        2> "$working/evaluator.stderr.log"
    evaluator_status=$?
    set -e
    if [[ $evaluator_status -ne 0 \
          || ! -s $working/evaluation/SHA256SUMS \
          || ! -s $working/evaluation/evaluation.json ]] \
       || ! (cd "$working/evaluation" && sha256sum -c SHA256SUMS >/dev/null); then
        formal_status=invalid_run
        rm -rf "$working/evaluation"
        mkdir "$working/evaluation"
        "$dev_python" - "$working/evaluation/evaluation.json" "$holdout_id" \
            "$policy_id" "$evaluator_status" <<'PY'
from pathlib import Path
import json, sys
output, holdout, policy, exit_status = sys.argv[1:]
payload = {
    "schema_version": "ploidypatch.external_holdout_non_evaluation.v0.5",
    "holdout_id": holdout,
    "policy_id": policy,
    "formal_outcome": "invalid_run",
    "stage": "formal_evaluator",
    "evaluator_invoked": True,
    "evaluator_exit_status": int(exit_status),
    "rule_relaxation": False,
}
Path(output).write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
PY
        (cd "$working/evaluation" && sha256sum evaluation.json > SHA256SUMS)
        printf '%s\n' "$formal_status" > "$working/formal_status.txt"
    fi
else
    "$dev_python" - "$working/evaluation/evaluation.json" "$holdout_id" \
        "$policy_id" "$formal_status" <<'PY'
from pathlib import Path
import json, sys
output, holdout, policy, status = sys.argv[1:]
payload = {
    "schema_version": "ploidypatch.external_holdout_non_evaluation.v0.5",
    "holdout_id": holdout,
    "policy_id": policy,
    "formal_outcome": status,
    "evaluator_invoked": False,
    "rule_relaxation": False,
}
Path(output).write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
PY
    (cd "$working/evaluation" && sha256sum evaluation.json > SHA256SUMS)
fi

"$dev_python" - "$working/reveal_status.json" "$holdout_id" "$policy_id" \
    "$formal_status" "$evaluator_invoked" <<'PY'
from pathlib import Path
import json, sys
output, holdout, policy, status, invoked = sys.argv[1:]
Path(output).write_text(
    json.dumps(
        {
            "schema_version": "ploidypatch.external_holdout_reveal_status.v0.5",
            "holdout_id": holdout,
            "policy_id": policy,
            "status": status,
            "evaluator_invoked": invoked == "true",
            "rule_relaxation": False,
        },
        indent=2,
        sort_keys=True,
    ) + "\n",
    encoding="utf-8",
)
PY

(
    cd "$working"
    find . -type f ! -path ./SHA256SUMS -printf '%P\0' \
        | sort -z | xargs -0 sha256sum > SHA256SUMS
    sha256sum -c SHA256SUMS >/dev/null
)
while IFS= read -r -d '' path; do chmod 0440 "$path"; done < <(find "$working" -type f -print0)
while IFS= read -r -d '' path; do chmod 0550 "$path"; done < <(find "$working" -depth -type d -print0)
mv "$working" "$output"
printf 'External holdout reveal frozen (%s): %s\n' "$formal_status" "$output"
[[ $formal_status != invalid_run ]]
