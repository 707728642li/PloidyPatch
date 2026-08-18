#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then echo "usage: $0 PROJECT_ROOT" >&2; exit 2; fi
project_root=$(realpath "$1")
code_root=$project_root/code
policy=$code_root/config/apple_external_validation_policy_v0.3.tsv
protocol=$code_root/docs/APPLE_EXTERNAL_VALIDATION_PROTOCOL_v0.3.md
public_root=$project_root/data/public/apple_external_v0.3
normalized_root=$project_root/data/derived/external_inputs/apple_v0.3
model_root=$project_root/results/models/support_conditioned_ranker_v0.3
model=$model_root/model.json
result_root=$project_root/results/protocol_freezes/apple_external_v0.3
working_root=${result_root}.working

synteny_env=$project_root/envs/ploidypatch-synteny
syngap_env=$project_root/envs/ploidypatch-syngap
baseline_env=$project_root/envs/ploidypatch-baseline
gemoma_env=$project_root/envs/ploidypatch-gemoma
lifton_env=$project_root/envs/ploidypatch-lifton
diamond_bin=$synteny_env/bin/diamond
wgdi_bin=$synteny_env/bin/wgdi
gffread_bin=$syngap_env/bin/gffread
miniprot_bin=$baseline_env/bin/miniprot
gemoma_bin=$gemoma_env/bin/GeMoMa
gemoma_jar=$gemoma_env/share/gemoma-1.9-0/GeMoMa-1.9.jar
lifton_bin=$lifton_env/bin/lifton

modules=(
    src/ploidypatch/structure_perturb.py
    src/ploidypatch/score.py
    src/ploidypatch/copy_pair_sampling.py
    src/ploidypatch/homeolog_pairs.py
    src/ploidypatch/self_wgd_pairs.py
    src/ploidypatch/consensus.py
    src/ploidypatch/copy_features.py
    src/ploidypatch/homeolog_topology.py
    src/ploidypatch/wgd_candidate_select.py
    src/ploidypatch/support_ranker.py
    scripts/freeze_apple_external_protocol_v0.3.sh
)
for required in "$policy" "$protocol" "$public_root/SHA256SUMS" \
    "$normalized_root/SHA256SUMS" "$model_root/SHA256SUMS" "$model" \
    "$diamond_bin" "$wgdi_bin" "$gffread_bin" "$miniprot_bin" \
    "$gemoma_bin" "$gemoma_jar" "$lifton_bin"; do
    [[ -s $required ]] || { echo "missing apple protocol-freeze input: $required" >&2; exit 1; }
done
for relative in "${modules[@]}"; do
    [[ -s $code_root/$relative ]] || { echo "missing frozen module: $relative" >&2; exit 1; }
done
[[ $("$diamond_bin" version | awk '{print $3}') == 2.2.2 ]]
[[ $("$wgdi_bin" --version 2>&1 | tail -1) == 0.75 ]]
[[ $("$gffread_bin" --version 2>&1 | head -1) == 0.12.9 ]]
[[ $("$miniprot_bin" --version 2>&1 | head -1) == 0.18-r281 ]]
[[ $("$lifton_env/bin/python" -c 'import importlib.metadata; print(importlib.metadata.version("lifton"))') == 1.0.11 ]]
grep -q '"version": "1.9"' "$gemoma_env"/conda-meta/gemoma-1.9-*.json
[[ ! -e $result_root && ! -e $working_root ]] || {
    echo "refusing to overwrite apple protocol freeze" >&2; exit 1;
}
(cd "$public_root" && sha256sum -c SHA256SUMS >/dev/null)
(cd "$normalized_root" && sha256sum -c SHA256SUMS >/dev/null)
(cd "$model_root" && sha256sum -c SHA256SUMS >/dev/null)

policy_value() {
    awk -F '\t' -v wanted="$1" '$1 == wanted { print $2; found=1 } END { if (!found) exit 1 }' "$policy"
}
[[ $(policy_value public_bundle_sha256sums_sha256) == $(sha256sum "$public_root/SHA256SUMS" | awk '{print $1}') ]]
[[ $(policy_value normalized_bundle_sha256sums_sha256) == $(sha256sum "$normalized_root/SHA256SUMS" | awk '{print $1}') ]]
[[ $(policy_value model_sha256) == $(sha256sum "$model" | awk '{print $1}') ]]

mkdir -p "$working_root"
cp "$policy" "$working_root/policy.tsv"
cp "$protocol" "$working_root/protocol.md"
{
    printf 'field\tvalue\n'
    printf 'code_commit\t%s\n' "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'stage\tpre_pair_pre_event_pre_candidate_protocol_freeze\n'
    printf 'hidden_pair_enumeration\tfalse\n'
    printf 'hidden_event_generation\tfalse\n'
    printf 'candidate_generation\tfalse\n'
    printf 'candidate_score_generation\tfalse\n'
    printf 'external_label_access\tfalse\n'
    printf 'reference_roles_mutable\tfalse\n'
    printf 'automatic_approval\tfalse\n'
} > "$working_root/run_contract.tsv"
{
    printf 'role\tbytes\tsha256\tpath\n'
    for entry in \
        "public_SHA256SUMS:$public_root/SHA256SUMS" \
        "normalized_SHA256SUMS:$normalized_root/SHA256SUMS" \
        "model:$model" "model_SHA256SUMS:$model_root/SHA256SUMS" \
        "policy:$policy" "protocol:$protocol"; do
        role=${entry%%:*}; path=${entry#*:}
        printf '%s\t%s\t%s\t%s\n' "$role" "$(stat -Lc %s "$path")" \
            "$(sha256sum "$path" | awk '{print $1}')" "$path"
    done
} > "$working_root/input_manifest.tsv"
{
    printf 'tool\tversion\tsha256\tpath\n'
    printf 'DIAMOND\t2.2.2\t%s\t%s\n' "$(sha256sum "$diamond_bin" | awk '{print $1}')" "$diamond_bin"
    printf 'WGDI\t0.75\t%s\t%s\n' "$(sha256sum "$wgdi_bin" | awk '{print $1}')" "$wgdi_bin"
    printf 'gffread\t0.12.9\t%s\t%s\n' "$(sha256sum "$gffread_bin" | awk '{print $1}')" "$gffread_bin"
    printf 'miniprot\t0.18-r281\t%s\t%s\n' "$(sha256sum "$miniprot_bin" | awk '{print $1}')" "$miniprot_bin"
    printf 'GeMoMa-wrapper\t1.9\t%s\t%s\n' "$(sha256sum "$gemoma_bin" | awk '{print $1}')" "$gemoma_bin"
    printf 'GeMoMa-jar\t1.9\t%s\t%s\n' "$(sha256sum "$gemoma_jar" | awk '{print $1}')" "$gemoma_jar"
    printf 'LiftOn\t1.0.11\t%s\t%s\n' "$(sha256sum "$lifton_bin" | awk '{print $1}')" "$lifton_bin"
} > "$working_root/software.tsv"
{
    printf 'path\tsha256\n'
    for relative in "${modules[@]}"; do
        printf '%s\t%s\n' "$relative" "$(sha256sum "$code_root/$relative" | awk '{print $1}')"
    done
} > "$working_root/code_manifest.tsv"
(
    cd "$working_root"
    find . -type f ! -name SHA256SUMS -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS
    sha256sum -c SHA256SUMS >/dev/null
)
mv "$working_root" "$result_root"
printf 'apple v0.3 external protocol frozen: %s\n' "$result_root"
