#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 || $# -gt 3 ]]; then
    echo "usage: $0 PROJECT_ROOT [NEW_BLIND_RUN_DIR] [EXECUTION_FREEZE]" >&2
    exit 2
fi
project_root=$(realpath "$1")
output=${2:-$project_root/results/blind_runs/populus_external_v0.4}
output=$(realpath -m "$output")
allowed_root=$project_root/results/blind_runs
case "$output/" in
    "$allowed_root"/*/) ;;
    *) echo "blind run output must be a new child of $allowed_root" >&2; exit 2 ;;
esac
[[ $output != /nas_data/* ]] || { echo "blind output may not be on /nas_data" >&2; exit 2; }

execution=${3:-$project_root/results/protocol_freezes/populus_external_v0.4_execution}
execution=$(realpath "$execution")
case "$execution" in
    "$project_root/results/protocol_freezes/populus_external_v0.4_execution"*) ;;
    *) echo "execution freeze must be a Populus v0.4 project freeze" >&2; exit 2 ;;
esac
protocol=$project_root/results/protocol_freezes/populus_external_v0.4
model=$project_root/results/models/ploidypatch_ranker_v0.4
stage=$project_root/data/derived/external_inputs/populus_v0.4
shared_target=$stage/shared_target
candidate_only=$stage/candidate_only
blind_annotation=$project_root/benchmark/structure/copy_collapse_v0.4/ptr_v4/annotation_copy_collapse_seed20260930/blind
working=${output}.working

for required in \
    "$execution/SHA256SUMS" "$execution/execution_manifest.json" \
    "$execution/implementation_manifest.tsv" "$execution/environment_bindings.tsv" \
    "$execution/source/scripts/run_populus_blind_pipeline_v0.4.sh" \
    "$execution/source/scripts/finalize_populus_blind_custody_v0.4.py" \
    "$protocol/SHA256SUMS" "$model/SHA256SUMS" \
    "$stage/SHA256SUMS" "$stage/BLIND_SHA256SUMS" "$shared_target" "$candidate_only" \
    "$blind_annotation/perturbed.gff3" "$blind_annotation/blind_manifest.json" \
    "$blind_annotation/SHA256SUMS"; do
    [[ -s $required || -d $required ]] || { echo "missing blind-run input: $required" >&2; exit 1; }
done
for root in "$execution" "$protocol" "$model"; do
    (cd "$root" && sha256sum -c SHA256SUMS >/dev/null)
done
(cd "$stage" && sha256sum -c SHA256SUMS >/dev/null)
if ! awk -F '  ' '
    NF != 2 || $2 ~ /^\// || $2 ~ /(^|\/)\.\.($|\/)/ ||
    $2 !~ /^(shared_target|candidate_only)\// { exit 1 }
' "$stage/BLIND_SHA256SUMS"; then
    echo "blind-role checksum list escapes shared_target/candidate_only" >&2
    exit 1
fi
(cd "$stage" && sha256sum -c BLIND_SHA256SUMS >/dev/null)
if ! awk -F '  ' '
    NF != 2 || $2 ~ /^\// || $2 ~ /(^|\/)\.\.($|\/)/ ||
    tolower($2) ~ /(truth|label|evaluator|complete)/ { exit 1 }
' "$blind_annotation/SHA256SUMS"; then
    echo "blind annotation checksum list contains a forbidden path" >&2
    exit 1
fi
(cd "$blind_annotation" && sha256sum -c SHA256SUMS >/dev/null)
for root in "$shared_target" "$candidate_only" "$blind_annotation"; do
    if find "$root" -type l -print -quit | grep -q .; then
        echo "symlinks are forbidden in isolated blind input: $root" >&2
        exit 1
    fi
done
if find "$shared_target" -type f \( -iname '*gff*' -o -iname '*gtf*' \
        -o -iname '*protein*' -o -iname '*.pep*' -o -iname '*.faa*' \) \
        -print -quit | grep -q .; then
    echo "shared target mount contains annotation/protein bytes" >&2
    exit 1
fi
mapfile -t candidate_species < <(find "$candidate_only" -mindepth 1 -maxdepth 1 \
    -type d -printf '%f\n' | sort)
[[ ${candidate_species[*]} == "Salix_purpurea Salix_suchowensis" ]] || {
    echo "candidate-only mount species differ from frozen Salix roles" >&2
    exit 1
}

expected_runner=$(awk -F '\t' '$1 == "scripts/run_populus_blind_isolated_v0.4.sh" {print $3}' \
    "$execution/implementation_manifest.tsv")
[[ $expected_runner =~ ^[0-9a-f]{64}$ && $(sha256sum "$0" | awk '{print $1}') == "$expected_runner" ]] || {
    echo "launcher bytes differ from the execution freeze" >&2
    exit 1
}
command -v bwrap >/dev/null || { echo "bubblewrap is required" >&2; exit 1; }
command -v conda >/dev/null || { echo "conda is required for frozen-environment verification" >&2; exit 1; }
[[ ! -e $output && ! -e $working ]] || {
    echo "refusing to overwrite blind run or retained failed attempt" >&2
    exit 1
}
mkdir -p "$allowed_root" "$working/project" "$working/environment_checks"

# Re-resolve every environment from the freeze, reproduce both locks, and only
# then mount those exact prefixes read-only.  A changed conda or pip package set
# invalidates the attempt before target code starts.
declare -a environment_names=()
declare -a environment_prefixes=()
lock_python=
while IFS=$'\t' read -r name prefix explicit_relative explicit_sha pip_relative pip_sha; do
    [[ $name == name ]] && continue
    [[ $name =~ ^[a-z0-9][a-z0-9_.-]*$ && -d $prefix ]] || {
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
code_commit=$("$lock_python" -c \
    'import json,sys; print(json.load(open(sys.argv[1]))["code_commit"])' \
    "$execution/execution_manifest.json")
[[ $code_commit =~ ^[0-9a-f]{40}$ ]] || { echo "invalid execution-freeze commit" >&2; exit 1; }

internal=/run/project
mkdir -p \
    "$working/project/code" \
    "$working/project/data/derived/external_inputs/populus_v0.4/shared_target" \
    "$working/project/data/derived/external_inputs/populus_v0.4/candidate_only" \
    "$working/project/benchmark/structure/copy_collapse_v0.4/ptr_v4/annotation_copy_collapse_seed20260930/blind" \
    "$working/project/results/protocol_freezes/populus_external_v0.4" \
    "$working/project/results/protocol_freezes/populus_external_v0.4_execution" \
    "$working/project/results/models/ploidypatch_ranker_v0.4" \
    "$working/project/envs"

bwrap_version=$(bwrap --version | tr '\n' ' ')
conda_base=$(realpath "$(conda info --base)")
declare -a bwrap_args=(
    --unshare-all --unshare-net --die-with-parent --new-session --cap-drop ALL
    --clearenv
    --setenv HOME /home/blind
    --setenv TMPDIR /tmp
    --setenv LANG C.UTF-8
    --setenv LC_ALL C.UTF-8
    --setenv PLOIDYPATCH_BLIND_RUNNER 1
    --setenv PLOIDYPATCH_CODE_COMMIT "$code_commit"
    --setenv PYTHONPATH "$internal/code/src"
    --setenv PATH "$conda_base/condabin:/usr/bin:/bin"
    --proc /proc --dev /dev --tmpfs /tmp --dir /home --dir /home/blind --dir /run
    --ro-bind /usr /usr
    --ro-bind /etc /etc
)
if [[ $(realpath /bin) == /usr/bin ]]; then
    # Ubuntu exposes /bin as a host symlink into /usr.  Mounting /usr alone
    # does not recreate that symlink in an empty bubblewrap root, while MMseqs
    # generated workflows use #!/bin/sh.  Recreate the same immutable alias.
    bwrap_args+=(--symlink usr/bin /bin)
elif [[ -d /bin ]]; then
    bwrap_args+=(--ro-bind /bin /bin)
else
    echo "host /bin cannot be represented in the blind namespace" >&2
    exit 1
fi
[[ ! -d /lib ]] || bwrap_args+=(--ro-bind /lib /lib)
[[ ! -d /lib64 ]] || bwrap_args+=(--ro-bind /lib64 /lib64)

# Conda entry points contain absolute shebangs.  Bind the frozen prefixes (and
# the read-only conda launcher) at their original absolute locations and again
# at the canonical isolated project paths; both mounts are read-only.
declare -A made_parent=()
add_parent_dir() {
    local parent=$1
    if [[ -z ${made_parent[$parent]+x} ]]; then
        bwrap_args+=(--dir "$parent")
        made_parent[$parent]=1
    fi
}
add_parent_dir "$(dirname "$conda_base")"
bwrap_args+=(--ro-bind "$conda_base" "$conda_base")
bwrap_args+=(--bind "$working/project" "$internal")
for index in "${!environment_names[@]}"; do
    name=${environment_names[$index]}
    prefix=${environment_prefixes[$index]}
    add_parent_dir "$(dirname "$prefix")"
    bwrap_args+=(--ro-bind "$prefix" "$prefix")
    mkdir -p "$working/project/envs/$name"
    bwrap_args+=(--ro-bind "$prefix" "$internal/envs/$name")
done

bwrap_args+=(
    --ro-bind "$execution/source" "$internal/code"
    --ro-bind "$execution" "$internal/results/protocol_freezes/populus_external_v0.4_execution"
    --ro-bind "$protocol" "$internal/results/protocol_freezes/populus_external_v0.4"
    --ro-bind "$model" "$internal/results/models/ploidypatch_ranker_v0.4"
    --ro-bind "$shared_target" "$internal/data/derived/external_inputs/populus_v0.4/shared_target"
    --ro-bind "$candidate_only" "$internal/data/derived/external_inputs/populus_v0.4/candidate_only"
    --ro-bind "$blind_annotation" "$internal/benchmark/structure/copy_collapse_v0.4/ptr_v4/annotation_copy_collapse_seed20260930/blind"
    --chdir "$internal"
    /usr/bin/bash "$internal/code/scripts/run_populus_blind_pipeline_v0.4.sh" "$internal"
)
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

"$lock_python" \
    "$execution/source/scripts/finalize_populus_blind_custody_v0.4.py" \
    --blind-project-root "$working/project" \
    --execution-freeze "$execution" --protocol-freeze "$protocol" \
    --composite-model-freeze "$model" \
    --shared-target "$shared_target" --candidate-only "$candidate_only" \
    --blind-annotation-root "$blind_annotation" \
    --blind-role-checksums "$stage/BLIND_SHA256SUMS" \
    --runner-command "$working/bwrap_command.txt" --bwrap-version "$bwrap_version" \
    --output "$working/custody_manifest.json"
(
    cd "$working"
    find . -type f ! -name SHA256SUMS -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS
    sha256sum -c SHA256SUMS >/dev/null
)
while IFS= read -r -d '' path; do chmod 0440 "$path"; done < <(find "$working" -type f -print0)
while IFS= read -r -d '' path; do chmod 0550 "$path"; done < <(find "$working" -depth -type d -print0)
mv "$working" "$output"
printf 'Populus v0.4 isolated blind run frozen: %s\n' "$output"
