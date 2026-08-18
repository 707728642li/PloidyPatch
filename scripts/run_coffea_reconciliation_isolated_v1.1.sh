#!/usr/bin/env bash
set -euo pipefail
umask 077

namespace_entry() {
    local project_root=$1 python=/frozen/envs/ploidypatch-dev/bin/python
    [[ -d $project_root && -x $python ]] || return 80
    [[ ${PLOIDYPATCH_NETWORK_ACCESS:-} == none ]] || return 81
    for required in /holdout/role_manifest.json /holdout/role_manifest.tsv \
        /holdout/SHA256SUMS /holdout/blind_benchmark/perturbed.gff3 \
        /blind-run-a/project /blind-run-b/project /frozen/protocol/SHA256SUMS \
        /frozen/execution/SHA256SUMS; do
        [[ -e $required ]] || return 82
    done
    for forbidden in /holdout/evaluator_only /holdout/target_complete \
        /holdout/truth /holdout/labels /frozen/model /nas_data; do
        [[ ! -e $forbidden ]] || return 83
    done
    "$python" - /run/blind-run/namespace_role_validation.json <<'PY'
from pathlib import Path
import hashlib,json,os,sys
sha=lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
payload={'schema_version':'ploidypatch.coffea_core_h1_namespace_validation.v1.0',
 'holdout_id':os.environ['PLOIDYPATCH_HOLDOUT_ID'],
 'mount_manifest_sha256':sha('/run/blind-run/mount_manifest.json'),
 'host_role_manifest_sha256':sha('/holdout/role_manifest.json'),
 'role_root_SHA256SUMS_sha256':sha('/holdout/SHA256SUMS'),
 'blind_benchmark_manifest_sha256':sha('/holdout/blind_benchmark/blind_manifest.json'),
 'shared_target_visible':True,'candidate_only_visible':True,'blind_benchmark_visible':True,
 'reproducibility_run_a_visible':True,'reproducibility_run_b_visible':True,
 'evaluator_only_visible':False,'complete_target_annotation_visible':False,
 'truth_visible':False,'labels_visible':False,'nas_data_visible':False,
 'model_visible':False,'ranker_visible':False}
out=Path(sys.argv[1]);out.write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n')
PY
    cat > "$project_root/pipeline_commands.tsv" <<'EOF'
stage_order	entry	purpose
10	scripts/reconcile_coffea_blind_replicates_v1.1.py	exact dual-run projection stability intersection and six pool rebuilds
EOF
    "$python" /frozen/source/scripts/reconcile_coffea_blind_replicates_v1.1.py \
        --pools-a /blind-run-a/project/results/copy_collapse/external/coffea_v1.0_h1.invalid_run \
        --pools-b /blind-run-b/project/results/copy_collapse/external/coffea_v1.0_h1 \
        --base-gff /holdout/blind_benchmark/perturbed.gff3 \
        --formal-project-root "$project_root" --role-root /holdout \
        --protocol-freeze /frozen/protocol --execution-freeze /frozen/execution \
        --output-dir "$project_root/results/copy_collapse/external/coffea_v1.0_h1"
}

if [[ ${1:-} == --namespace-entry ]]; then
    [[ $# -eq 2 ]] || exit 64
    namespace_entry "$2"
    exit $?
fi

[[ $# -eq 6 ]] || {
    echo "usage: $0 PROTOCOL EXECUTION ROLE_ROOT RUN_A RUN_B OUTPUT" >&2
    exit 64
}
protocol=$(realpath "$1"); execution=$(realpath "$2"); role_root=$(realpath "$3")
run_a=$(realpath "$4"); run_b=$(realpath "$5"); output=$6; working=${output}.working
for root in "$protocol" "$execution" "$role_root" "$run_a" "$run_b"; do
    [[ -d $root && ! -L $root ]] || exit 65
done
[[ ! -e $output && ! -L $output && ! -e $working && ! -L $working ]] || exit 66
command -v bwrap >/dev/null; command -v conda >/dev/null
(cd "$role_root" && sha256sum -c SHA256SUMS >/dev/null)

host_python=$(command -v python3)
readarray -t frozen < <(PYTHONPATH="$execution/source/src" "$host_python" - \
    "$protocol" "$execution" "$run_a" "$run_b" <<'PY'
from pathlib import Path
import sys
from ploidypatch.coffea_h1_framework import verify_execution
from ploidypatch.reproducible_projection import verify_tree_manifest
protocol,execution,run_a,run_b=map(Path,sys.argv[1:])
manifest,_,contract=verify_execution(execution,protocol)
patch=manifest['execution_patch']
if patch.get('patch_sequence') != 3:
    raise SystemExit('Coffea reconciliation requires patch sequence 3')
for root,key in ((run_a,'reproducibility_run_a_manifest_relative_path'),
                 (run_b,'reproducibility_run_b_manifest_relative_path')):
    relative=patch[key]
    verify_tree_manifest(root=root,manifest=execution/relative)
print(contract.holdout_id);print(manifest['code_commit'])
for item in manifest['environments']:
    print(f"ENV\t{item['name']}\t{item['host_prefix']}")
PY
)
holdout_id=${frozen[0]}; code_commit=${frozen[1]}
[[ $holdout_id == coffea_et39_v1.0 ]] || exit 67
declare -a environment_names=() environment_prefixes=()
for line in "${frozen[@]:2}"; do
    IFS=$'\t' read -r marker name prefix <<< "$line"; [[ $marker == ENV ]] || exit 68
    environment_names+=("$name"); environment_prefixes+=("$(realpath "$prefix")")
done
[[ ${#environment_names[@]} -eq 6 ]] || exit 69

mkdir -p "$working/project/results/baselines/coffea"
cp -al "$run_b/project/results/baselines/coffea/v1.0" \
    "$working/project/results/baselines/coffea/v1.0"

mount_tsv=$working/mounts.tsv
printf 'role\thost_path\tnamespace_path\tread_only\n' > "$mount_tsv"
record_mount() { printf '%s\t%s\t%s\t%s\n' "$1" "$2" "$3" "$4" >> "$mount_tsv"; }
conda_base=$(realpath "$(conda info --base)"); [[ $conda_base != /nas_data* ]] || exit 70
internal=/run/blind-run/project
declare -a args=(--unshare-all --unshare-net --die-with-parent --new-session
 --cap-drop ALL --clearenv --setenv HOME /home/blind --setenv TMPDIR /tmp
 --setenv LANG C.UTF-8 --setenv LC_ALL C.UTF-8
 --setenv PLOIDYPATCH_BLIND_RUNNER 1 --setenv PLOIDYPATCH_NETWORK_ACCESS none
 --setenv PLOIDYPATCH_HOLDOUT_ID "$holdout_id" --setenv PLOIDYPATCH_CODE_COMMIT "$code_commit"
 --setenv PYTHONPATH /frozen/source/src --setenv PATH "$conda_base/condabin:/usr/bin:/bin"
 --proc /proc --dev /dev --tmpfs /tmp --dir /home --dir /home/blind --dir /run
 --dir /frozen --dir /frozen/envs --dir /holdout --ro-bind /usr /usr --ro-bind /etc /etc)
record_mount system_usr /usr /usr true; record_mount system_etc /etc /etc true
declare -A made_parent=()
add_parent_dir() {
    local current= part; local -a parts=()
    IFS='/' read -r -a parts <<< "${1#/}"
    for part in "${parts[@]}"; do
        [[ -n $part ]] || continue; current=$current/$part
        if [[ -z ${made_parent[$current]+x} ]]; then args+=(--dir "$current"); made_parent[$current]=1; fi
    done
}
if [[ $(realpath /bin) == /usr/bin ]]; then args+=(--symlink usr/bin /bin); else
    args+=(--ro-bind /bin /bin); record_mount system_bin /bin /bin true; fi
for library in /lib /lib64; do
    if [[ -d $library ]]; then args+=(--ro-bind "$library" "$library");
        record_mount "system_${library#/}" "$library" "$library" true; fi
done
add_parent_dir "$(dirname "$conda_base")"
args+=(--ro-bind "$conda_base" "$conda_base" --bind "$working" /run/blind-run)
record_mount system_conda "$conda_base" "$conda_base" true
record_mount blind_output "$working" /run/blind-run false
for index in "${!environment_names[@]}"; do
    name=${environment_names[$index]}; prefix=${environment_prefixes[$index]}
    [[ $prefix != /nas_data* ]] || exit 71
    add_parent_dir "$(dirname "$prefix")"
    args+=(--ro-bind "$prefix" "$prefix" --ro-bind "$prefix" "/frozen/envs/$name")
    record_mount "frozen_environment:$name" "$prefix" "$prefix" true
    record_mount "frozen_environment:$name" "$prefix" "/frozen/envs/$name" true
done
args+=(--ro-bind "$execution" /frozen/execution --ro-bind "$execution/source" /frozen/source
 --ro-bind "$protocol" /frozen/protocol --ro-bind "$role_root/shared_target" /holdout/shared_target
 --ro-bind "$role_root/candidate_only" /holdout/candidate_only
 --ro-bind "$role_root/blind_benchmark" /holdout/blind_benchmark
 --ro-bind "$role_root/role_manifest.json" /holdout/role_manifest.json
 --ro-bind "$role_root/role_manifest.tsv" /holdout/role_manifest.tsv
 --ro-bind "$role_root/SHA256SUMS" /holdout/SHA256SUMS
 --ro-bind "$run_a" /blind-run-a --ro-bind "$run_b" /blind-run-b
 --chdir "$internal" /usr/bin/bash
 /frozen/source/scripts/run_coffea_reconciliation_isolated_v1.1.sh
 --namespace-entry "$internal")
record_mount frozen_execution "$execution" /frozen/execution true
record_mount frozen_source "$execution/source" /frozen/source true
record_mount frozen_protocol "$protocol" /frozen/protocol true
for role in shared_target candidate_only blind_benchmark; do
    record_mount "$role" "$role_root/$role" "/holdout/$role" true
done
record_mount system_role_metadata "$role_root/role_manifest.json" /holdout/role_manifest.json true
record_mount system_role_metadata "$role_root/role_manifest.tsv" /holdout/role_manifest.tsv true
record_mount system_role_metadata "$role_root/SHA256SUMS" /holdout/SHA256SUMS true
record_mount reproducibility_run_a "$run_a" /blind-run-a true
record_mount reproducibility_run_b "$run_b" /blind-run-b true

PYTHONPATH="$execution/source/src" "$host_python" - "$mount_tsv" \
    "$working/mount_manifest.json" "$holdout_id" <<'PY'
import csv,json,sys
with open(sys.argv[1],encoding='utf-8',newline='') as h:
    rows=list(csv.DictReader(h,delimiter='\t'))
for row in rows: row['read_only']=row['read_only']=='true'
payload={'schema_version':'ploidypatch.coffea_core_h1_mount_manifest.v1.0',
 'holdout_id':sys.argv[3],'mounts':rows,'ranker_or_model_mounted':False,
 'network_access':False,'reproducibility_patch':True}
with open(sys.argv[2],'x',encoding='utf-8') as h:
    json.dump(payload,h,indent=2,sort_keys=True);h.write('\n')
PY
rm "$mount_tsv"
printf '%q ' bwrap "${args[@]}" > "$working/bwrap_command.txt"; printf '\n' >> "$working/bwrap_command.txt"
set +e; bwrap "${args[@]}" > "$working/stdout.log" 2> "$working/stderr.log"; status=$?; set -e
printf '%s\n' "$status" > "$working/exit_status.txt"
[[ $status -eq 0 ]] || { echo "reconciliation pipeline failed; working tree retained" >&2; exit "$status"; }
dev_python=
for index in "${!environment_names[@]}"; do
    [[ ${environment_names[$index]} == ploidypatch-dev ]] && dev_python=${environment_prefixes[$index]}/bin/python
done
[[ -x $dev_python ]] || exit 72
PYTHONPATH="$execution/source/src" "$dev_python" \
 "$execution/source/scripts/finalize_coffea_blind_custody_v1.0.py" \
 --blind-project "$working/project" --blind-run "$working" --execution-freeze "$execution" \
 --protocol-freeze "$protocol" --blind-role-root "$role_root" \
 --runner-command "$working/bwrap_command.txt" --mount-manifest "$working/mount_manifest.json" \
 --namespace-validation "$working/namespace_role_validation.json" \
 --output "$working/custody_manifest.json"
PYTHONPATH="$execution/source/src" "$dev_python" - "$working" <<'PY'
from pathlib import Path
import sys
from ploidypatch.artifact_manifest import write_sha256sums, verify_sha256sums
root=Path(sys.argv[1]); write_sha256sums(root); verify_sha256sums(root,ignore_checksum_file=True)
PY
chmod -R a-w "$working"
mv "$working" "$output"
