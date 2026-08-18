#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then echo "usage: $0 PROJECT_ROOT" >&2; exit 2; fi
project_root=$(realpath "$1")
code_root=$project_root/code
python_bin=$project_root/envs/ploidypatch-dev/bin/python
input_root=${PLOIDYPATCH_STAGED_INPUT_ROOT:?PLOIDYPATCH_STAGED_INPUT_ROOT is required}
blind_benchmark_root=${PLOIDYPATCH_BLIND_BENCHMARK_ROOT:?PLOIDYPATCH_BLIND_BENCHMARK_ROOT is required}
blind_gff=$blind_benchmark_root/perturbed.gff3
upstream=$project_root/results/baselines/actinidia_v0.5
protocol_root=${PLOIDYPATCH_PROTOCOL_FREEZE:?PLOIDYPATCH_PROTOCOL_FREEZE is required}
execution_root=${PLOIDYPATCH_EXECUTION_FREEZE:?PLOIDYPATCH_EXECUTION_FREEZE is required}
contract_path=${PLOIDYPATCH_HOLDOUT_CONTRACT:?PLOIDYPATCH_HOLDOUT_CONTRACT is required}
context_verifier=$code_root/scripts/verify_external_holdout_blind_context_v0.5.py
declare -A raw=(
    [gemoma_actinidia_eriantha]="$upstream/gemoma/actinidia_eriantha/upstream/final_annotation.gff"
    [gemoma_actinidia_rufa]="$upstream/gemoma/actinidia_rufa/upstream/final_annotation.gff"
    [lifton_actinidia_eriantha]="$upstream/lifton/actinidia_eriantha/upstream/lifton.gff3"
    [lifton_actinidia_rufa]="$upstream/lifton/actinidia_rufa/upstream/lifton.gff3"
)
miniprot_root=$upstream/miniprot
miniprot_gff=$miniprot_root/raw/miniprot.gff3
protein_map=$miniprot_root/reference/actinidia_candidate_refs.map.tsv
result_root=$project_root/results/copy_collapse/external/actinidia_v0.5_method_trio
working_root=${result_root}.working

[[ ${PLOIDYPATCH_BLIND_RUNNER:-} == 1 ]] || {
    echo "Actinidia candidate pools must run inside the frozen blind runner" >&2; exit 1;
}
[[ ${PLOIDYPATCH_NETWORK_ACCESS:-} == none ]] || {
    echo "Actinidia candidate pools require a network-disabled namespace" >&2; exit 1;
}
for forbidden in /nas_data "$input_root/evaluator_only" "$blind_benchmark_root/evaluator" \
    "$blind_benchmark_root/truth" "$blind_benchmark_root/complete"; do
    [[ ! -e $forbidden ]] || { echo "forbidden blind-runner path is visible: $forbidden" >&2; exit 1; }
done
if [[ -r /proc/self/mountinfo ]] && grep -Eq '/nas_data|/evaluator_only|/target_complete|/truth_references' /proc/self/mountinfo; then
    echo "forbidden evaluator or NAS mount detected in blind namespace" >&2; exit 1
fi

verify_tree() { (cd "$1" && sha256sum -c SHA256SUMS >/dev/null); }
verify_implementation() {
    local relative=$1 expected observed
    expected=$(awk -F '\t' -v path="$relative" 'NR > 1 && $1 == path {print $3}' \
        "$execution_root/implementation_manifest.tsv")
    observed=$(sha256sum "$code_root/$relative" | awk '{print $1}')
    [[ $expected =~ ^[0-9a-f]{64}$ && $observed == "$expected" ]] || {
        echo "implementation differs from execution freeze: $relative" >&2; exit 1;
    }
}

for required in "$python_bin" "$blind_gff" "$protocol_root/SHA256SUMS" \
    "$blind_benchmark_root/blind_manifest.json" "$blind_benchmark_root/SHA256SUMS" \
    "$execution_root/SHA256SUMS" "$miniprot_root/SHA256SUMS" "$miniprot_gff" \
    "$protein_map" "$contract_path" "$input_root/role_manifest.tsv" \
    "$input_root/role_contract.json" "$context_verifier" "${raw[@]}"; do
    [[ -s $required ]] || { echo "missing Actinidia method-pool input: $required" >&2; exit 1; }
done
verify_tree "$protocol_root"
verify_tree "$execution_root"
verify_tree "$miniprot_root"
verify_implementation scripts/build_actinidia_method_trio_candidate_pools_v0.5.sh
verify_implementation scripts/verify_external_holdout_blind_context_v0.5.py
verify_implementation src/ploidypatch/consensus.py
verify_implementation src/ploidypatch/baseline.py
[[ ! -e $result_root && ! -e $working_root ]] || {
    echo "refusing to overwrite Actinidia method trio" >&2; exit 1;
}
mkdir -p "$working_root"/{merged,freeze} \
    "$working_root/methods"/{miniprot,gemoma,lifton}/blind \
    "$working_root/consensus"/{primary_union,legacy_union,support2,support3}/blind
cp /proc/self/mountinfo "$working_root/freeze/blind_runner.mountinfo"
PYTHONPATH="$code_root/src" "$python_bin" "$context_verifier" \
    --input-root "$input_root" --contract "$contract_path" \
    --protocol-freeze "$protocol_root" --execution-freeze "$execution_root" \
    --blind-benchmark-root "$blind_benchmark_root" \
    --expected-holdout-id actinidia_red5_v0.5 --expected-primary-chromosomes 29 \
    --output-json "$working_root/freeze/blind_context.json"
{
    printf 'field\tvalue\ncode_commit\t%s\n' "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'split\ttarget_level_predeclared_untouched_secondary_replication\ncandidate_truth_access\tfalse\n'
    printf 'complete_target_annotation_access\tfalse\nevaluator_reference_access\tfalse\n'
    printf 'method_families\tminiprot,gemoma,lifton\n'
    printf 'candidate_references\tActinidia_eriantha_White,Actinidia_rufa_ARU\n'
    printf 'within_method_reference_vote_count\t1\n'
    printf 'adapter_min_identity\t0.5\nadapter_min_query_coverage\t0.5\n'
    printf 'adapter_max_existing_cds_overlap\t0.2\nadapter_max_redundancy_overlap\t0.5\n'
    printf 'primary_candidate_policy\tretain_distinct_phased_CDS_chains\n'
    printf 'legacy_candidate_policy\tsuppress_overlapping\nautomatic_approval\tfalse\n'
    printf 'protocol_freeze_sha256sums_sha256\t%s\n' \
        "$(sha256sum "$protocol_root/SHA256SUMS" | awk '{print $1}')"
    printf 'execution_freeze_sha256sums_sha256\t%s\n' \
        "$(sha256sum "$execution_root/SHA256SUMS" | awk '{print $1}')"
} > "$working_root/run_contract.tsv"
{
    printf 'role\tbytes\tsha256\tpath\n'
    for entry in "blind_gff:$blind_gff" "miniprot_gff:$miniprot_gff" \
        "protein_map:$protein_map" \
        "gemoma_actinidia_eriantha:${raw[gemoma_actinidia_eriantha]}" \
        "gemoma_actinidia_rufa:${raw[gemoma_actinidia_rufa]}" \
        "lifton_actinidia_eriantha:${raw[lifton_actinidia_eriantha]}" \
        "lifton_actinidia_rufa:${raw[lifton_actinidia_rufa]}"; do
        role=${entry%%:*}; path=${entry#*:}
        printf '%s\t%s\t%s\t%s\n' "$role" "$(stat -Lc %s "$path")" \
            "$(sha256sum "$path" | awk '{print $1}')" "$path"
    done
} > "$working_root/input_manifest.tsv"

cd "$code_root"
for method in gemoma lifton; do
    "$python_bin" -m ploidypatch.cli baseline merge-candidate-gffs \
        --candidate "actinidia_eriantha=${raw[${method}_actinidia_eriantha]}" \
        --candidate "actinidia_rufa=${raw[${method}_actinidia_rufa]}" \
        --output-gff "$working_root/merged/$method.gff3" \
        --provenance-tsv "$working_root/merged/$method.provenance.tsv" \
        > "$working_root/merged/$method.stdout.json" \
        2> "$working_root/merged/$method.stderr.log"
done

pids=(); labels=()
(
    "$python_bin" -m ploidypatch.cli baseline adapt-miniprot \
        --perturbed-gff "$blind_gff" --miniprot-gff "$miniprot_gff" \
        --protein-map "$protein_map" --min-identity 0.5 \
        --min-query-coverage 0.5 --max-existing-cds-overlap 0.2 \
        --max-redundancy-overlap 0.5 \
        --output-gff "$working_root/methods/miniprot/blind/candidate.gff3" \
        --decisions-tsv "$working_root/methods/miniprot/blind/decisions.tsv" \
        > "$working_root/methods/miniprot/blind/stdout.json" \
        2> "$working_root/methods/miniprot/blind/stderr.log"
) & pids+=("$!"); labels+=(miniprot)
for method in gemoma lifton; do
    (
        "$python_bin" -m ploidypatch.cli baseline adapt-gff \
            --perturbed-gff "$blind_gff" \
            --candidate-gff "$working_root/merged/$method.gff3" \
            --source "$method" --max-existing-cds-overlap 0.2 \
            --max-redundancy-overlap 0.5 \
            --output-gff "$working_root/methods/$method/blind/candidate.gff3" \
            --decisions-tsv "$working_root/methods/$method/blind/decisions.tsv" \
            > "$working_root/methods/$method/blind/stdout.json" \
            2> "$working_root/methods/$method/blind/stderr.log"
    ) & pids+=("$!"); labels+=("$method")
done
failed=0
for i in "${!pids[@]}"; do
    if ! wait "${pids[$i]}"; then echo "failed Actinidia adapt: ${labels[$i]}" >&2; failed=1; fi
done
[[ $failed -eq 0 ]] || exit 1

pids=(); labels=()
for pool in primary_union legacy_union support2 support3; do
    case $pool in
        primary_union) support=1; redundancy=retain_distinct_chains ;;
        legacy_union) support=1; redundancy=suppress_overlapping ;;
        support2) support=2; redundancy=suppress_overlapping ;;
        support3) support=3; redundancy=suppress_overlapping ;;
    esac
    (
        "$python_bin" -m ploidypatch.cli baseline select-method-consensus \
            --base-gff "$blind_gff" \
            --candidate "miniprot=$working_root/methods/miniprot/blind/candidate.gff3" \
            --candidate "gemoma=$working_root/methods/gemoma/blind/candidate.gff3" \
            --candidate "lifton=$working_root/methods/lifton/blind/candidate.gff3" \
            --min-method-support "$support" --max-redundancy-overlap 0.5 \
            --redundancy-policy "$redundancy" \
            --output-gff "$working_root/consensus/$pool/blind/candidate.gff3" \
            --decisions-tsv "$working_root/consensus/$pool/blind/decisions.tsv" \
            > "$working_root/consensus/$pool/blind/stdout.json" \
            2> "$working_root/consensus/$pool/blind/stderr.log"
    ) & pids+=("$!"); labels+=("$pool")
done
failed=0
for i in "${!pids[@]}"; do
    if ! wait "${pids[$i]}"; then echo "failed Actinidia pool: ${labels[$i]}" >&2; failed=1; fi
done
[[ $failed -eq 0 ]] || exit 1

for pool in primary_union legacy_union support2 support3; do
    for output in candidate.gff3 decisions.tsv candidate.gff3.manifest.json; do
        [[ -s $working_root/consensus/$pool/blind/$output ]] || {
            echo "missing Actinidia pool artifact: $pool/$output" >&2; exit 1;
        }
    done
done
{
    printf 'role\tbytes\tsha256\tpath\n'
    find "$working_root/methods" "$working_root/consensus" -type f \
        \( -name candidate.gff3 -o -name decisions.tsv -o -name '*.manifest.json' \) \
        -print0 | sort -z | while IFS= read -r -d '' path; do
            role=${path#"$working_root/"}; role=${role//\//_}
            printf '%s\t%s\t%s\t%s\n' "$role" "$(stat -Lc %s "$path")" \
                "$(sha256sum "$path" | awk '{print $1}')" "$path"
        done
} > "$working_root/freeze/candidate_freeze.tsv"
(
    cd "$working_root"
    find . -type f \( -name '*.gff3' -o -name '*.json' -o -name '*.tsv' -o -name '*.mountinfo' \) \
        -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS
    sha256sum -c SHA256SUMS >/dev/null
)
du -sb "$working_root" > "$working_root/disk_bytes.txt"
mv "$working_root" "$result_root"
printf 'Actinidia truth-blind method pools frozen: %s\n' "$result_root"
