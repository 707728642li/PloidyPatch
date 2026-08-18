#!/usr/bin/env bash
set -euo pipefail

namespace_entry() {
    [[ $# -eq 2 ]] || { echo "namespace entry requires PROJECT_ROOT PIPELINE" >&2; exit 2; }
    local project_root=$1
    local pipeline=$2
    local python=/frozen/envs/ploidypatch-dev/bin/python
    [[ -x $python && -f $pipeline ]] || { echo "missing frozen namespace entry" >&2; exit 1; }
    "$python" - \
        /frozen/protocol/contract.json \
        /holdout/blind_role_manifest.json \
        /run/blind-run/mount_manifest.json \
        /run/blind-run/namespace_role_validation.json <<'PY'
from pathlib import Path
import hashlib
import json
import sys

contract_path, role_path, mount_path, output_path = map(Path, sys.argv[1:])
contract = json.loads(contract_path.read_text(encoding="utf-8"))
role = json.loads(role_path.read_text(encoding="utf-8"))
mount = json.loads(mount_path.read_text(encoding="utf-8"))
holdout = contract["holdout_id"]
if role.get("holdout_id") != holdout or mount.get("holdout_id") != holdout:
    raise SystemExit("namespace holdout identity differs")
if role.get("roles") != ["shared_target", "candidate_only", "blind_benchmark"]:
    raise SystemExit("namespace role manifest differs")
required = [
    Path("/holdout/shared_target"), Path("/holdout/candidate_only"),
    Path("/holdout/blind_benchmark"),
]
if not all(path.is_dir() for path in required):
    raise SystemExit("namespace lacks candidate-safe roles")
for forbidden in (
    Path("/holdout/evaluator_only"), Path("/holdout/target_complete"),
    Path("/holdout/truth"), Path("/holdout/labels"), Path("/nas_data")
):
    if forbidden.exists():
        raise SystemExit(f"forbidden namespace path is visible: {forbidden}")
def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()
payload = {
    "schema_version": "ploidypatch.blind_namespace_validation.v0.5",
    "holdout_id": holdout,
    "mount_manifest_sha256": sha(mount_path),
    "host_role_manifest_sha256": sha(role_path),
    "shared_target_visible": True,
    "candidate_only_visible": True,
    "blind_benchmark_visible": True,
    "blind_benchmark_manifest_sha256": sha(Path("/holdout/blind_benchmark/blind_manifest.json")),
    "evaluator_only_visible": False,
    "truth_visible": False,
    "complete_target_annotation_visible": False,
    "nas_data_visible": False,
}
with output_path.open("x", encoding="utf-8", newline="") as handle:
    json.dump(payload, handle, indent=2, sort_keys=True)
    handle.write("\n")
PY
    export PLOIDYPATCH_BLIND_RUNNER=1
    export PLOIDYPATCH_NETWORK_ACCESS=none
    export PLOIDYPATCH_STAGED_INPUT_ROOT=/holdout
    export PLOIDYPATCH_BLIND_BENCHMARK_ROOT=/holdout/blind_benchmark
    export PLOIDYPATCH_HOLDOUT_CONTRACT=/frozen/protocol/contract.json
    export PLOIDYPATCH_PROTOCOL_FREEZE=/frozen/protocol
    export PLOIDYPATCH_EXECUTION_FREEZE=/frozen/execution
    export PLOIDYPATCH_COMPOSITE_MODEL_FREEZE=/frozen/model
    export PLOIDYPATCH_BLIND_OUTPUT_ROOT=/run/blind-run
    /usr/bin/bash "$pipeline" "$project_root"
}

if [[ ${1:-} == --namespace-entry ]]; then
    shift
    namespace_entry "$@"
    exit $?
fi

if [[ $# -ne 6 ]]; then
    echo "usage: $0 PROJECT_ROOT EXECUTION_FREEZE PROTOCOL_FREEZE MODEL_FREEZE BLIND_ROLE_ROOT NEW_BLIND_RUN_DIR" >&2
    exit 2
fi

project_root=$(realpath "$1")
execution=$(realpath "$2")
protocol=$(realpath "$3")
model=$(realpath "$4")
role_root=$(realpath "$5")
output=$(realpath -m "$6")
working=${output}.working
allowed_root=$project_root/results/blind_runs
case "$output/" in
    "$allowed_root"/*/) ;;
    *) echo "blind output must be a new child of $allowed_root" >&2; exit 2 ;;
esac
for path in "$execution" "$protocol" "$model" "$role_root" "$output"; do
    [[ $path != /nas_data && $path != /nas_data/* ]] || {
        echo "blind execution path may not be on /nas_data: $path" >&2; exit 2;
    }
done
[[ ! -e $output && ! -e $working ]] || {
    echo "refusing to overwrite blind run or retained failed attempt" >&2; exit 1;
}

for required in \
    "$execution/SHA256SUMS" "$execution/execution_manifest.json" \
    "$execution/implementation_manifest.tsv" "$execution/environment_bindings.tsv" \
    "$protocol/SHA256SUMS" "$protocol/protocol_manifest.json" \
    "$protocol/contract.json" "$protocol/role_manifest.tsv" \
    "$protocol/role_contract.json" "$model/SHA256SUMS" \
    "$role_root/SHA256SUMS" "$role_root/role_manifest.json" \
    "$role_root/shared_target" "$role_root/candidate_only" \
    "$role_root/blind_benchmark" "$role_root/blind_benchmark/SHA256SUMS" \
    "$role_root/blind_benchmark/blind_manifest.json" \
    "$role_root/blind_benchmark/perturbed.gff3"; do
    [[ -s $required || -d $required ]] || { echo "missing blind input: $required" >&2; exit 1; }
done
for root in "$execution" "$protocol" "$model" "$role_root"; do
    (cd "$root" && sha256sum -c SHA256SUMS >/dev/null)
done
for root in "$role_root/shared_target" "$role_root/candidate_only" "$role_root/blind_benchmark"; do
    if find "$root" -type l -print -quit | grep -q .; then
        echo "symlinks are forbidden in blind data roles: $root" >&2
        exit 1
    fi
done
if find "$role_root/shared_target" -type f \
        \( -iname '*complete*' -o -iname '*truth*' -o -iname '*label*' \) \
        -print -quit | grep -q .; then
    echo "shared target role contains a forbidden file name" >&2
    exit 1
fi

command -v bwrap >/dev/null || { echo "bubblewrap is required" >&2; exit 1; }
command -v conda >/dev/null || { echo "conda is required" >&2; exit 1; }

mkdir -p "$allowed_root" "$working/project" "$working/environment_checks"

declare -a environment_names=()
declare -a environment_prefixes=()
lock_python=
while IFS=$'\t' read -r name prefix explicit_relative explicit_sha pip_relative pip_sha; do
    [[ $name == name ]] && continue
    [[ $name =~ ^[a-z0-9][a-z0-9_.-]*$ && -d $prefix \
       && $prefix != /nas_data && $prefix != /nas_data/* ]] || {
        echo "invalid frozen environment binding: $name=$prefix" >&2; exit 1;
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
        echo "environment differs from execution freeze: $name" >&2; exit 1;
    }
    environment_names+=("$name")
    environment_prefixes+=("$(realpath "$prefix")")
    [[ $name != ploidypatch-dev ]] || lock_python=$prefix/bin/python
done < "$execution/environment_bindings.tsv"
[[ ${#environment_names[@]} -eq 7 ]] || { echo "seven frozen environments are required" >&2; exit 1; }
[[ -x $lock_python ]] || { echo "frozen ploidypatch-dev Python is required" >&2; exit 1; }

PYTHONPATH="$execution/source/src" "$lock_python" - \
    "$execution" "$protocol" "$model" "$role_root" <<'PY'
from pathlib import Path
import json
import sys
from ploidypatch.artifact_manifest import sha256_file, verify_sha256sums

execution, protocol, model, role = map(Path, sys.argv[1:])
for root in (execution, protocol, model, role):
    verify_sha256sums(root, ignore_checksum_file=True)
e = json.loads((execution / "execution_manifest.json").read_text(encoding="utf-8"))
p = json.loads((protocol / "protocol_manifest.json").read_text(encoding="utf-8"))
r = json.loads((role / "role_manifest.json").read_text(encoding="utf-8"))
if e.get("schema_version") != "ploidypatch.external_holdout_execution_freeze.v0.5":
    raise SystemExit("wrong execution schema")
if p.get("schema_version") != "ploidypatch.external_holdout_protocol_freeze.v0.5":
    raise SystemExit("wrong protocol schema")
patch = e.get("execution_patch")
if patch is None:
    if (
        e.get("freeze_stage") != "post_metadata_pre_pair_pre_candidate_pre_label"
        or e.get("code_commit") != p.get("code_commit")
    ):
        raise SystemExit("base execution stage or commit differs")
elif (
    not isinstance(patch, dict)
    or patch.get("schema_version") != "ploidypatch.external_holdout_execution_patch.v0.5"
    or patch.get("freeze_stage") != "post_evaluator_truth_failed_blind_pre_candidate_pre_score_pre_label_execution_patch"
    or e.get("freeze_stage") != patch.get("freeze_stage")
    or patch.get("base_code_commit") != p.get("code_commit")
    or patch.get("patch_code_commit") != e.get("code_commit")
    or patch.get("base_protocol_SHA256SUMS_sha256") != sha256_file(protocol / "SHA256SUMS")
    or not isinstance(patch.get("superseded_execution_SHA256SUMS_sha256"), str)
    or len(patch["superseded_execution_SHA256SUMS_sha256"]) != 64
    or patch.get("scientific_protocol_changed") is not False
    or patch.get("contract_or_policy_changed") is not False
    or patch.get("model_or_threshold_changed") is not False
    or patch.get("staged_inputs_changed") is not False
    or patch.get("truth_or_benchmark_regenerated") is not False
    or patch.get("evaluator_truth_construction_completed_before_patch") is not True
    or patch.get("blind_candidate_wgd_completed_before_patch") is not False
    or patch.get("candidate_generation_completed_before_patch") is not False
    or patch.get("formal_scores_generated_before_patch") is not False
    or patch.get("truth_labels_accessed_before_patch") is not False
    or not (execution / "superseded_failed_attempt_manifest.tsv").is_file()
    or not (execution / "patch_reason.md").is_file()
):
    raise SystemExit("execution patch provenance or scientific firewall differs")
if not e.get("holdout_id") or e.get("holdout_id") != p.get("holdout_id") or r.get("holdout_id") != e.get("holdout_id"):
    raise SystemExit("holdout identities differ")
if r.get("schema_version") != "ploidypatch.blind_role_manifest.v0.5":
    raise SystemExit("wrong blind role schema")
if r.get("roles") != ["shared_target", "candidate_only", "blind_benchmark"]:
    raise SystemExit("blind role list differs")
if any(r.get(key) is not False for key in (
    "truth_access", "complete_target_annotation_present", "evaluator_references_present"
)):
    raise SystemExit("blind role manifest violates firewall")
benchmark = role / "blind_benchmark"
verify_sha256sums(benchmark, ignore_checksum_file=True)
blind = json.loads((benchmark / "blind_manifest.json").read_text(encoding="utf-8"))
if (
    r.get("contract_sha256") != sha256_file(protocol / "contract.json")
    or r.get("protocol_SHA256SUMS_sha256") != sha256_file(protocol / "SHA256SUMS")
    or r.get("blind_benchmark_SHA256SUMS_sha256") != sha256_file(benchmark / "SHA256SUMS")
    or r.get("blind_benchmark_manifest_sha256") != sha256_file(benchmark / "blind_manifest.json")
    or r.get("network_access") is not False
    or blind.get("schema_version") != "ploidypatch.blind_benchmark_input.v0.5"
    or blind.get("truth_access") is not False
    or blind.get("complete_target_annotation_access") is not False
    or blind.get("perturbed_annotation", {}).get("file_name") != "perturbed.gff3"
    or blind.get("perturbed_annotation", {}).get("sha256") != sha256_file(benchmark / "perturbed.gff3")
):
    raise SystemExit("blind benchmark role or frozen bindings differ")
PY

readarray -t frozen_values < <(PYTHONPATH="$execution/source/src" "$lock_python" - \
    "$execution/execution_manifest.json" <<'PY'
import json, sys
value = json.load(open(sys.argv[1], encoding="utf-8"))
print(value["holdout_id"])
print(value["code_commit"])
print(value["pipeline_entries"]["blind_pipeline"])
PY
)
holdout_id=${frozen_values[0]}
code_commit=${frozen_values[1]}
blind_pipeline=${frozen_values[2]}
[[ $holdout_id =~ ^[a-z0-9][a-z0-9._-]+$ && $code_commit =~ ^[0-9a-f]{40}$ ]] || {
    echo "malformed frozen holdout identity or commit" >&2; exit 1;
}
[[ -f $execution/source/$blind_pipeline ]] || { echo "missing frozen blind pipeline" >&2; exit 1; }

expected_runner=$(awk -F '\t' '$1 == "scripts/run_external_holdout_blind_isolated_v0.5.sh" {print $3}' \
    "$execution/implementation_manifest.tsv")
[[ $expected_runner =~ ^[0-9a-f]{64}$ && $(sha256sum "$0" | awk '{print $1}') == "$expected_runner" ]] || {
    echo "launcher bytes differ from execution freeze" >&2; exit 1;
}

internal=/run/blind-run/project
mkdir -p "$working/project/envs"
bwrap_version=$(bwrap --version | tr '\n' ' ')
conda_base=$(realpath "$(conda info --base)")
[[ $conda_base != /nas_data && $conda_base != /nas_data/* ]] || {
    echo "conda base may not be mounted from /nas_data" >&2; exit 1;
}

mount_tsv=$working/mounts.tsv
printf 'role\thost_path\tnamespace_path\tread_only\n' > "$mount_tsv"
record_mount() {
    printf '%s\t%s\t%s\t%s\n' "$1" "$2" "$3" "$4" >> "$mount_tsv"
}

declare -a bwrap_args=(
    --unshare-all --unshare-net --die-with-parent --new-session --cap-drop ALL
    --clearenv
    --setenv HOME /home/blind
    --setenv TMPDIR /tmp
    --setenv LANG C.UTF-8
    --setenv LC_ALL C.UTF-8
    --setenv PLOIDYPATCH_BLIND_RUNNER 1
    --setenv PLOIDYPATCH_NETWORK_ACCESS none
    --setenv PLOIDYPATCH_HOLDOUT_ID "$holdout_id"
    --setenv PLOIDYPATCH_CODE_COMMIT "$code_commit"
    --setenv PLOIDYPATCH_STAGED_INPUT_ROOT /holdout
    --setenv PLOIDYPATCH_BLIND_BENCHMARK_ROOT /holdout/blind_benchmark
    --setenv PLOIDYPATCH_HOLDOUT_CONTRACT /frozen/protocol/contract.json
    --setenv PLOIDYPATCH_PROTOCOL_FREEZE /frozen/protocol
    --setenv PLOIDYPATCH_EXECUTION_FREEZE /frozen/execution
    --setenv PLOIDYPATCH_COMPOSITE_MODEL_FREEZE /frozen/model
    --setenv PLOIDYPATCH_BLIND_OUTPUT_ROOT /run/blind-run
    --setenv PYTHONPATH /frozen/source/src
    --setenv PATH "$conda_base/condabin:/usr/bin:/bin"
    --proc /proc --dev /dev --tmpfs /tmp --dir /home --dir /home/blind --dir /run
    --dir /frozen --dir /frozen/envs --dir /holdout
    --ro-bind /usr /usr
    --ro-bind /etc /etc
)
record_mount system_usr /usr /usr true
record_mount system_etc /etc /etc true
if [[ $(realpath /bin) == /usr/bin ]]; then
    bwrap_args+=(--symlink usr/bin /bin)
elif [[ -d /bin ]]; then
    bwrap_args+=(--ro-bind /bin /bin)
    record_mount system_bin /bin /bin true
else
    echo "host /bin cannot be represented in blind namespace" >&2
    exit 1
fi
if [[ -d /lib ]]; then
    bwrap_args+=(--ro-bind /lib /lib)
    record_mount system_lib /lib /lib true
fi
if [[ -d /lib64 ]]; then
    bwrap_args+=(--ro-bind /lib64 /lib64)
    record_mount system_lib64 /lib64 /lib64 true
fi

declare -A made_parent=()
add_parent_dir() {
    local parent=$1
    local current= part
    local -a parts=()
    IFS='/' read -r -a parts <<< "${parent#/}"
    for part in "${parts[@]}"; do
        [[ -n $part ]] || continue
        current=$current/$part
        if [[ -z ${made_parent[$current]+x} ]]; then
            bwrap_args+=(--dir "$current")
            made_parent[$current]=1
        fi
    done
}
add_parent_dir "$(dirname "$conda_base")"
bwrap_args+=(--ro-bind "$conda_base" "$conda_base")
record_mount system_conda "$conda_base" "$conda_base" true
bwrap_args+=(--bind "$working" /run/blind-run)
record_mount blind_output "$working" /run/blind-run false
for index in "${!environment_names[@]}"; do
    name=${environment_names[$index]}
    prefix=${environment_prefixes[$index]}
    add_parent_dir "$(dirname "$prefix")"
    bwrap_args+=(--ro-bind "$prefix" "$prefix")
    bwrap_args+=(--ro-bind "$prefix" "/frozen/envs/$name")
    bwrap_args+=(--ro-bind "$prefix" "$internal/envs/$name")
    record_mount "frozen_environment:$name" "$prefix" "/frozen/envs/$name" true
    record_mount "frozen_environment:$name" "$prefix" "$prefix" true
    record_mount "frozen_environment:$name" "$prefix" "$internal/envs/$name" true
done

bwrap_args+=(
    --ro-bind "$execution" /frozen/execution
    --ro-bind "$execution/source" /frozen/source
    --ro-bind "$execution/source" "$internal/code"
    --ro-bind "$protocol" /frozen/protocol
    --ro-bind "$model" /frozen/model
    --ro-bind "$role_root/shared_target" /holdout/shared_target
    --ro-bind "$role_root/candidate_only" /holdout/candidate_only
    --ro-bind "$role_root/blind_benchmark" /holdout/blind_benchmark
    --ro-bind "$role_root/role_manifest.json" /holdout/blind_role_manifest.json
    --ro-bind "$protocol/role_manifest.tsv" /holdout/role_manifest.tsv
    --ro-bind "$protocol/role_contract.json" /holdout/role_contract.json
    --chdir "$internal"
    /usr/bin/bash /frozen/source/scripts/run_external_holdout_blind_isolated_v0.5.sh
    --namespace-entry "$internal" "/frozen/source/$blind_pipeline"
)
record_mount frozen_execution "$execution" /frozen/execution true
record_mount frozen_source "$execution/source" /frozen/source true
record_mount frozen_source "$execution/source" "$internal/code" true
record_mount frozen_protocol "$protocol" /frozen/protocol true
record_mount frozen_model "$model" /frozen/model true
record_mount shared_target "$role_root/shared_target" /holdout/shared_target true
record_mount candidate_only "$role_root/candidate_only" /holdout/candidate_only true
record_mount blind_benchmark "$role_root/blind_benchmark" /holdout/blind_benchmark true
record_mount system_role_metadata "$role_root/role_manifest.json" /holdout/blind_role_manifest.json true
record_mount system_role_metadata "$protocol/role_manifest.tsv" /holdout/role_manifest.tsv true
record_mount system_role_metadata "$protocol/role_contract.json" /holdout/role_contract.json true

"$lock_python" - "$mount_tsv" "$holdout_id" "$working/mount_manifest.json" <<'PY'
import csv, json, sys
from pathlib import Path
source, holdout, output = sys.argv[1:]
with open(source, encoding="utf-8", newline="") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))
mounts = [
    {
        "role": row["role"],
        "host_path": row["host_path"],
        "namespace_path": row["namespace_path"],
        "read_only": row["read_only"] == "true",
    }
    for row in rows
]
with open(output, "x", encoding="utf-8", newline="") as handle:
    json.dump(
        {"schema_version": "ploidypatch.blind_mount_manifest.v0.5", "holdout_id": holdout, "mounts": mounts},
        handle, indent=2, sort_keys=True,
    )
    handle.write("\n")
PY
rm "$mount_tsv"

printf '%q ' bwrap "${bwrap_args[@]}" > "$working/bwrap_command.txt"
printf '\n' >> "$working/bwrap_command.txt"

set +e
bwrap "${bwrap_args[@]}" > "$working/stdout.log" 2> "$working/stderr.log"
status=$?
set -e
printf '%s\n' "$status" > "$working/exit_status.txt"
if [[ $status -ne 0 ]]; then
    echo "isolated blind run failed; retained without reveal at $working" >&2
    exit "$status"
fi

for required in "$working/namespace_role_validation.json"; do
    [[ -s $required ]] || { echo "missing namespace validation: $required" >&2; exit 1; }
done

PYTHONPATH="$execution/source/src" "$lock_python" \
    "$execution/source/scripts/finalize_external_holdout_blind_custody_v0.5.py" \
    --blind-project-root "$working/project" \
    --blind-output-root "$working" \
    --execution-freeze "$execution" \
    --protocol-freeze "$protocol" \
    --composite-model-freeze "$model" \
    --blind-role-root "$role_root" \
    --runner-command "$working/bwrap_command.txt" \
    --mount-manifest "$working/mount_manifest.json" \
    --namespace-validation "$working/namespace_role_validation.json" \
    --bwrap-version "$bwrap_version" \
    --output "$working/custody_manifest.json"

(
    cd "$working"
    find . -type f ! -path ./SHA256SUMS -printf '%P\0' \
        | sort -z | xargs -0 sha256sum > SHA256SUMS
    sha256sum -c SHA256SUMS >/dev/null
)
while IFS= read -r -d '' path; do chmod 0440 "$path"; done < <(find "$working" -type f -print0)
while IFS= read -r -d '' path; do chmod 0550 "$path"; done < <(find "$working" -depth -type d -print0)
mv "$working" "$output"
printf 'External holdout isolated blind run frozen: %s\n' "$output"
