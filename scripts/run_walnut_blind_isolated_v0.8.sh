#!/usr/bin/env bash
set -euo pipefail
umask 077

namespace_entry() {
    local project_root=$1 pipeline=$2 python=/frozen/envs/ploidypatch-dev/bin/python
    [[ -d $project_root && -f $pipeline && -x $python ]] || return 80
    [[ ${PLOIDYPATCH_NETWORK_ACCESS:-} == none ]] || return 81
    [[ ${PLOIDYPATCH_STAGED_INPUT_ROOT:-} == /holdout ]] || return 82
    [[ ${PLOIDYPATCH_BLIND_BENCHMARK_ROOT:-} == /holdout/blind_benchmark ]] || return 83
    [[ ${PLOIDYPATCH_PROTOCOL_FREEZE:-} == /frozen/protocol ]] || return 84
    [[ ${PLOIDYPATCH_EXECUTION_FREEZE:-} == /frozen/execution ]] || return 85
    [[ -z ${PLOIDYPATCH_COMPOSITE_MODEL_FREEZE:-} && -z ${PLOIDYPATCH_RANKER_FREEZE:-} ]] || return 86
    for required in /holdout/shared_target /holdout/candidate_only /holdout/blind_benchmark \
        /holdout/blind_role_manifest.json /holdout/role_manifest.tsv /holdout/role_contract.json; do
        [[ -e $required ]] || return 87
    done
    for forbidden in /holdout/evaluator_only /holdout/target_complete /holdout/truth \
        /holdout/labels /frozen/model /nas_data; do
        [[ ! -e $forbidden ]] || return 88
    done
    "$python" - /run/blind-run/namespace_role_validation.json <<'PY'
from pathlib import Path
import hashlib, json, os, sys
out=Path(sys.argv[1]); sha=lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
payload={'schema_version':'ploidypatch.walnut_core_h1_namespace_validation.v0.8',
 'holdout_id':os.environ['PLOIDYPATCH_HOLDOUT_ID'],
 'mount_manifest_sha256':sha('/run/blind-run/mount_manifest.json'),
 'host_role_manifest_sha256':sha('/holdout/blind_role_manifest.json'),
 'blind_benchmark_manifest_sha256':sha('/holdout/blind_benchmark/blind_manifest.json'),
 'shared_target_visible':True,'candidate_only_visible':True,'blind_benchmark_visible':True,
 'evaluator_only_visible':False,'complete_target_annotation_visible':False,
 'truth_visible':False,'labels_visible':False,'nas_data_visible':False,
 'model_visible':False,'ranker_visible':False}
fd=os.open(out,os.O_WRONLY|os.O_CREAT|os.O_EXCL,0o440)
with os.fdopen(fd,'w',encoding='utf-8') as h: json.dump(payload,h,indent=2,sort_keys=True);h.write('\n')
PY
    /usr/bin/bash "$pipeline" "$project_root"
}

if [[ ${1:-} == --namespace-entry ]]; then
    [[ $# -eq 3 ]] || exit 64
    namespace_entry "$2" "$3"
    exit $?
fi

[[ $# -eq 4 ]] || { echo "usage: $0 PROTOCOL EXECUTION ROLE_ROOT OUTPUT" >&2; exit 64; }
protocol=$(realpath "$1"); execution=$(realpath "$2"); role_root=$(realpath "$3")
output=$4; working=${output}.working
[[ -d $protocol && -d $execution && -d $role_root ]] || exit 65
[[ ! -e $output && ! -L $output && ! -e $working && ! -L $working ]] || exit 66
command -v bwrap >/dev/null; command -v conda >/dev/null
mkdir -p "$working/project/envs"

lock_python=$(command -v python3)
readarray -t frozen < <(PYTHONPATH="$execution/source/src" "$lock_python" - "$protocol" "$execution" <<'PY'
from pathlib import Path
import sys
from ploidypatch.walnut_h1_framework import verify_execution
manifest,_,contract=verify_execution(Path(sys.argv[2]),Path(sys.argv[1]))
print(contract.holdout_id);print(manifest['code_commit']);print(manifest['pipeline_entries']['blind_pipeline'])
for item in manifest['environments']: print(f"ENV\t{item['name']}\t{item['host_prefix']}")
PY
)
holdout_id=${frozen[0]}; code_commit=${frozen[1]}; pipeline_relative=${frozen[2]}
pipeline=$execution/source/$pipeline_relative
[[ -f $pipeline && $holdout_id == walnut_walnut2_v0.8 ]] || exit 67
declare -a environment_names=() environment_prefixes=()
for line in "${frozen[@]:3}"; do
    IFS=$'\t' read -r marker name prefix <<< "$line"; [[ $marker == ENV ]] || exit 68
    environment_names+=("$name"); environment_prefixes+=("$(realpath "$prefix")")
done
[[ ${#environment_names[@]} -eq 6 ]] || exit 69

mount_tsv=$working/mounts.tsv
printf 'role\thost_path\tnamespace_path\tread_only\n' > "$mount_tsv"
record_mount() { printf '%s\t%s\t%s\t%s\n' "$1" "$2" "$3" "$4" >> "$mount_tsv"; }
conda_base=$(realpath "$(conda info --base)"); [[ $conda_base != /nas_data* ]] || exit 70
internal=/run/blind-run/project
declare -a args=(--unshare-all --unshare-net --die-with-parent --new-session --cap-drop ALL --clearenv
 --setenv HOME /home/blind --setenv TMPDIR /tmp --setenv LANG C.UTF-8 --setenv LC_ALL C.UTF-8
 --setenv PLOIDYPATCH_BLIND_RUNNER 1 --setenv PLOIDYPATCH_NETWORK_ACCESS none
 --setenv PLOIDYPATCH_HOLDOUT_ID "$holdout_id" --setenv PLOIDYPATCH_CODE_COMMIT "$code_commit"
 --setenv PLOIDYPATCH_STAGED_INPUT_ROOT /holdout
 --setenv PLOIDYPATCH_BLIND_BENCHMARK_ROOT /holdout/blind_benchmark
 --setenv PLOIDYPATCH_HOLDOUT_CONTRACT /frozen/protocol/contract.json
 --setenv PLOIDYPATCH_PROTOCOL_FREEZE /frozen/protocol
 --setenv PLOIDYPATCH_EXECUTION_FREEZE /frozen/execution
 --setenv PLOIDYPATCH_BLIND_OUTPUT_ROOT /run/blind-run
 --setenv PYTHONPATH /frozen/source/src --setenv PATH "$conda_base/condabin:/usr/bin:/bin"
 --proc /proc --dev /dev --tmpfs /tmp --dir /home --dir /home/blind --dir /run
 --dir /frozen --dir /frozen/envs --dir /holdout --ro-bind /usr /usr --ro-bind /etc /etc)
record_mount system_usr /usr /usr true; record_mount system_etc /etc /etc true
declare -A made_parent=()
add_parent_dir() {
    local current= part
    local -a parts=()
    IFS='/' read -r -a parts <<< "${1#/}"
    for part in "${parts[@]}"; do
        [[ -n $part ]] || continue; current=$current/$part
        if [[ -z ${made_parent[$current]+x} ]]; then
            args+=(--dir "$current"); made_parent[$current]=1
        fi
    done
}
if [[ $(realpath /bin) == /usr/bin ]]; then args+=(--symlink usr/bin /bin); else
    args+=(--ro-bind /bin /bin); record_mount system_bin /bin /bin true; fi
for library in /lib /lib64; do if [[ -d $library ]]; then args+=(--ro-bind "$library" "$library");
    record_mount "system_${library#/}" "$library" "$library" true; fi; done
add_parent_dir "$(dirname "$conda_base")"
args+=(--ro-bind "$conda_base" "$conda_base" --bind "$working" /run/blind-run)
record_mount system_conda "$conda_base" "$conda_base" true
record_mount blind_output "$working" /run/blind-run false
for index in "${!environment_names[@]}"; do
    name=${environment_names[$index]}; prefix=${environment_prefixes[$index]}
    [[ $prefix != /nas_data* ]] || exit 71
    add_parent_dir "$(dirname "$prefix")"
    args+=(--ro-bind "$prefix" "$prefix" --ro-bind "$prefix" "/frozen/envs/$name" --ro-bind "$prefix" "$internal/envs/$name")
    record_mount "frozen_environment:$name" "$prefix" "$prefix" true
    record_mount "frozen_environment:$name" "$prefix" "/frozen/envs/$name" true
    record_mount "frozen_environment:$name" "$prefix" "$internal/envs/$name" true
done
args+=(--ro-bind "$execution" /frozen/execution --ro-bind "$execution/source" /frozen/source
 --ro-bind "$execution/source" "$internal/code" --ro-bind "$protocol" /frozen/protocol
 --ro-bind "$role_root/shared_target" /holdout/shared_target
 --ro-bind "$role_root/candidate_only" /holdout/candidate_only
 --ro-bind "$role_root/blind_benchmark" /holdout/blind_benchmark
 --ro-bind "$role_root/role_manifest.json" /holdout/blind_role_manifest.json
 --ro-bind "$protocol/role_manifest.tsv" /holdout/role_manifest.tsv
 --ro-bind "$protocol/role_contract.json" /holdout/role_contract.json --chdir "$internal"
 /usr/bin/bash /frozen/source/scripts/run_walnut_blind_isolated_v0.8.sh
 --namespace-entry "$internal" "/frozen/source/$pipeline_relative")
record_mount frozen_execution "$execution" /frozen/execution true
record_mount frozen_source "$execution/source" /frozen/source true
record_mount frozen_source "$execution/source" "$internal/code" true
record_mount frozen_protocol "$protocol" /frozen/protocol true
for role in shared_target candidate_only blind_benchmark; do record_mount "$role" "$role_root/$role" "/holdout/$role" true; done
record_mount system_role_metadata "$role_root/role_manifest.json" /holdout/blind_role_manifest.json true
record_mount system_role_metadata "$protocol/role_manifest.tsv" /holdout/role_manifest.tsv true
record_mount system_role_metadata "$protocol/role_contract.json" /holdout/role_contract.json true

PYTHONPATH="$execution/source/src" "$lock_python" - "$mount_tsv" "$working/mount_manifest.json" "$holdout_id" <<'PY'
import csv,json,os,sys
with open(sys.argv[1],encoding='utf-8',newline='') as h: rows=list(csv.DictReader(h,delimiter='\t'))
for row in rows: row['read_only']=row['read_only']=='true'
payload={'schema_version':'ploidypatch.walnut_core_h1_mount_manifest.v0.8','holdout_id':sys.argv[3],
 'mounts':rows,'ranker_or_model_mounted':False,'network_access':False}
fd=os.open(sys.argv[2],os.O_WRONLY|os.O_CREAT|os.O_EXCL,0o440)
with os.fdopen(fd,'w',encoding='utf-8') as h: json.dump(payload,h,indent=2,sort_keys=True);h.write('\n')
PY
rm "$mount_tsv"
bwrap_version=$(bwrap --version | tr '\n' ' ')
printf '%q ' bwrap "${args[@]}" > "$working/bwrap_command.txt"; printf '\n' >> "$working/bwrap_command.txt"
set +e; bwrap "${args[@]}" > "$working/stdout.log" 2> "$working/stderr.log"; status=$?; set -e
printf '%s\n' "$status" > "$working/exit_status.txt"
[[ $status -eq 0 ]] || { echo "blind pipeline failed; working tree retained" >&2; exit "$status"; }
dev_python=
for index in "${!environment_names[@]}"; do
    [[ ${environment_names[$index]} == ploidypatch-dev ]] && dev_python=${environment_prefixes[$index]}/bin/python
done
[[ -x $dev_python ]] || exit 72
PYTHONPATH="$execution/source/src" "$dev_python" "$execution/source/scripts/finalize_walnut_blind_custody_v0.8.py" \
 --blind-project-root "$working/project" --blind-output-root "$working" --execution-freeze "$execution" \
 --protocol-freeze "$protocol" --blind-role-root "$role_root" --runner-command "$working/bwrap_command.txt" \
 --mount-manifest "$working/mount_manifest.json" --namespace-validation "$working/namespace_role_validation.json" \
 --bwrap-version "$bwrap_version" --output "$working/custody_manifest.json"
mv "$working" "$output"
