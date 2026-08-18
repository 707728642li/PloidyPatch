#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then echo "usage: $0 PROJECT_ROOT" >&2; exit 2; fi
project_root=$(realpath "$1")
code_root=$project_root/code
python_bin=$project_root/envs/ploidypatch-dev/bin/python
input_root=$project_root/data/derived/external_inputs/populus_v0.4
benchmark=$project_root/benchmark/structure/copy_collapse_v0.4/ptr_v4/annotation_copy_collapse_seed20260930
blind_gff=$benchmark/blind/perturbed.gff3
upstream=$project_root/results/baselines/populus_v0.4
protocol_root=$project_root/results/protocol_freezes/populus_external_v0.4
execution_root=$project_root/results/protocol_freezes/populus_external_v0.4_execution
declare -A raw=(
    [gemoma_salix_purpurea]="$upstream/gemoma/salix_purpurea/upstream/final_annotation.gff"
    [gemoma_salix_suchowensis]="$upstream/gemoma/salix_suchowensis/upstream/final_annotation.gff"
    [lifton_salix_purpurea]="$upstream/lifton/salix_purpurea/upstream/lifton.gff3"
    [lifton_salix_suchowensis]="$upstream/lifton/salix_suchowensis/upstream/lifton.gff3"
)
miniprot_root=$upstream/miniprot
miniprot_gff=$miniprot_root/raw/miniprot.gff3
protein_map=$miniprot_root/reference/populus_candidate_refs.map.tsv
result_root=$project_root/results/copy_collapse/external/populus_v0.4_method_trio
working_root=${result_root}.working

[[ ${PLOIDYPATCH_BLIND_RUNNER:-} == 1 ]] || {
    echo "Populus candidate pools must run inside the frozen blind runner" >&2; exit 1;
}
for forbidden in /nas_data "$input_root/evaluator_only" "$benchmark/evaluator" \
    "$benchmark/truth" "$benchmark/complete"; do
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
    "$execution_root/SHA256SUMS" "$miniprot_root/SHA256SUMS" "$miniprot_gff" \
    "$protein_map" "${raw[@]}"; do
    [[ -s $required ]] || { echo "missing Populus method-pool input: $required" >&2; exit 1; }
done
verify_tree "$protocol_root"
verify_tree "$execution_root"
verify_tree "$miniprot_root"
verify_implementation scripts/build_populus_method_trio_candidate_pools_v0.4.sh
verify_implementation src/ploidypatch/consensus.py
verify_implementation src/ploidypatch/baseline.py
[[ ! -e $result_root && ! -e $working_root ]] || {
    echo "refusing to overwrite Populus method trio" >&2; exit 1;
}
mkdir -p "$working_root"/{merged,freeze} \
    "$working_root/methods"/{miniprot,gemoma,lifton}/blind \
    "$working_root/consensus"/{primary_union,legacy_union,support2,support3}/blind
cp /proc/self/mountinfo "$working_root/freeze/blind_runner.mountinfo"
{
    printf 'field\tvalue\ncode_commit\t%s\n' "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'split\tuntouched_confirmatory_external_v0.4\ncandidate_truth_access\tfalse\n'
    printf 'complete_target_annotation_access\tfalse\nevaluator_reference_access\tfalse\n'
    printf 'method_families\tminiprot,gemoma,lifton\n'
    printf 'candidate_references\tSalix_purpurea,Salix_suchowensis\n'
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
        "gemoma_salix_purpurea:${raw[gemoma_salix_purpurea]}" \
        "gemoma_salix_suchowensis:${raw[gemoma_salix_suchowensis]}" \
        "lifton_salix_purpurea:${raw[lifton_salix_purpurea]}" \
        "lifton_salix_suchowensis:${raw[lifton_salix_suchowensis]}"; do
        role=${entry%%:*}; path=${entry#*:}
        printf '%s\t%s\t%s\t%s\n' "$role" "$(stat -Lc %s "$path")" \
            "$(sha256sum "$path" | awk '{print $1}')" "$path"
    done
} > "$working_root/input_manifest.tsv"

cd "$code_root"
for method in gemoma lifton; do
    "$python_bin" -m ploidypatch.cli baseline merge-candidate-gffs \
        --candidate "salix_purpurea=${raw[${method}_salix_purpurea]}" \
        --candidate "salix_suchowensis=${raw[${method}_salix_suchowensis]}" \
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
    if ! wait "${pids[$i]}"; then echo "failed Populus adapt: ${labels[$i]}" >&2; failed=1; fi
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
    if ! wait "${pids[$i]}"; then echo "failed Populus pool: ${labels[$i]}" >&2; failed=1; fi
done
[[ $failed -eq 0 ]] || exit 1

for pool in primary_union legacy_union support2 support3; do
    for output in candidate.gff3 decisions.tsv candidate.gff3.manifest.json; do
        [[ -s $working_root/consensus/$pool/blind/$output ]] || {
            echo "missing Populus pool artifact: $pool/$output" >&2; exit 1;
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
printf 'Populus truth-blind method pools frozen: %s\n' "$result_root"
