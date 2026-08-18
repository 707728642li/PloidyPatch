#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then echo "usage: $0 PROJECT_ROOT" >&2; exit 2; fi
project_root=$(realpath "$1")
code_root=$project_root/code
python_bin=$project_root/envs/ploidypatch-dev/bin/python
benchmark=$project_root/benchmark/structure/copy_collapse_v0.3/mdx_gddh13/annotation_copy_collapse_seed20260831
blind_gff=$benchmark/blind/perturbed.gff3
source_gff=$project_root/data/derived/external_inputs/apple_v0.3/target_apple/primary_chromosomes.gff3
upstream=$project_root/results/baselines/apple_v0.3
protocol_root=$project_root/results/protocol_freezes/apple_external_v0.3
execution_root=$project_root/results/protocol_freezes/apple_external_v0.3_execution
result_root=$project_root/results/copy_collapse/external/apple_v0.3_method_trio
working_root=${result_root}.working
declare -A raw=(
    [gemoma_pear]="$upstream/gemoma/pear/upstream/final_annotation.gff"
    [gemoma_peach]="$upstream/gemoma/peach/upstream/final_annotation.gff"
    [lifton_pear]="$upstream/lifton/pear/upstream/lifton.gff3"
    [lifton_peach]="$upstream/lifton/peach/upstream/lifton.gff3"
)
miniprot_gff=$upstream/miniprot/raw/miniprot.gff3
protein_map=$upstream/miniprot/reference/apple_candidate_refs.map.tsv

for required in "$python_bin" "$benchmark/SHA256SUMS" "$blind_gff" \
    "$source_gff" "$protocol_root/SHA256SUMS" "$execution_root/SHA256SUMS" "$miniprot_gff" \
    "$protein_map" "${raw[@]}"; do
    [[ -s $required ]] || { echo "missing apple method input: $required" >&2; exit 1; }
done
(cd "$benchmark" && sha256sum -c SHA256SUMS >/dev/null)
(cd "$protocol_root" && sha256sum -c SHA256SUMS >/dev/null)
(cd "$execution_root" && sha256sum -c SHA256SUMS >/dev/null)
expected_self=$(awk -F '\t' '$1 == "scripts/build_apple_method_trio_candidate_pools_v0.3.sh" {print $3}' \
    "$execution_root/implementation_manifest.tsv")
[[ -n $expected_self && $(sha256sum "$code_root/scripts/build_apple_method_trio_candidate_pools_v0.3.sh" | awk '{print $1}') == "$expected_self" ]] || {
    echo "apple candidate-pool script differs from execution freeze" >&2; exit 1;
}
expected_consensus=$(awk -F '\t' '$1 == "src/ploidypatch/consensus.py" {print $2}' "$protocol_root/code_manifest.tsv")
observed_consensus=$(sha256sum "$code_root/src/ploidypatch/consensus.py" | awk '{print $1}')
[[ -n $expected_consensus && $observed_consensus == "$expected_consensus" ]] || {
    echo "post-freeze consensus implementation change detected" >&2; exit 1;
}
[[ ! -e $result_root && ! -e $working_root ]] || {
    echo "refusing to overwrite apple method trio" >&2; exit 1;
}
mkdir -p "$working_root"/{merged,freeze} \
    "$working_root/methods"/{miniprot,gemoma,lifton}/{blind,complete_control} \
    "$working_root/consensus"/{primary_union,legacy_union,support2,support3}/{blind,complete_control}
{
    printf 'field\tvalue\ncode_commit\t%s\n' \
        "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'split\tuntouched_external_v0.3\ncandidate_truth_access\tfalse\n'
    printf 'hidden_pair_access\tfalse\nexternal_label_access\tfalse\n'
    printf 'blind_and_control_adaptation\tindependent_same_raw_predictions\n'
    printf 'method_families\tminiprot,gemoma,lifton\n'
    printf 'candidate_references\tPyrus_communis,Prunus_persica\n'
    printf 'within_method_reference_vote_count\t1\n'
    printf 'primary_candidate_policy\tretain_distinct_phased_CDS_chains\n'
    printf 'legacy_candidate_policy\tsuppress_overlapping\n'
    printf 'automatic_approval\tfalse\n'
    printf 'protocol_freeze_sha256sums_sha256\t%s\n' \
        "$(sha256sum "$protocol_root/SHA256SUMS" | awk '{print $1}')"
    printf 'consensus_module_sha256\t%s\n' "$observed_consensus"
} > "$working_root/run_contract.tsv"
{
    printf 'role\tbytes\tsha256\tpath\n'
    for entry in "blind_gff:$blind_gff" "complete_control_gff:$source_gff" \
        "miniprot_gff:$miniprot_gff" "protein_map:$protein_map" \
        "gemoma_pear:${raw[gemoma_pear]}" "gemoma_peach:${raw[gemoma_peach]}" \
        "lifton_pear:${raw[lifton_pear]}" "lifton_peach:${raw[lifton_peach]}"; do
        role=${entry%%:*}; path=${entry#*:}
        printf '%s\t%s\t%s\t%s\n' "$role" "$(stat -Lc %s "$path")" \
            "$(sha256sum "$path" | awk '{print $1}')" "$path"
    done
} > "$working_root/input_manifest.tsv"

cd "$code_root"
for method in gemoma lifton; do
    "$python_bin" -m ploidypatch.cli baseline merge-candidate-gffs \
        --candidate "pear=${raw[${method}_pear]}" \
        --candidate "peach=${raw[${method}_peach]}" \
        --output-gff "$working_root/merged/$method.gff3" \
        --provenance-tsv "$working_root/merged/$method.provenance.tsv" \
        > "$working_root/merged/$method.stdout.json" \
        2> "$working_root/merged/$method.stderr.log"
done

pids=(); labels=()
for mode in blind complete_control; do
    annotation=$blind_gff; [[ $mode == complete_control ]] && annotation=$source_gff
    (
        "$python_bin" -m ploidypatch.cli baseline adapt-miniprot \
            --perturbed-gff "$annotation" --miniprot-gff "$miniprot_gff" \
            --protein-map "$protein_map" --min-identity 0.5 \
            --min-query-coverage 0.5 --max-existing-cds-overlap 0.2 \
            --max-redundancy-overlap 0.5 \
            --output-gff "$working_root/methods/miniprot/$mode/candidate.gff3" \
            --decisions-tsv "$working_root/methods/miniprot/$mode/decisions.tsv" \
            > "$working_root/methods/miniprot/$mode/stdout.json" \
            2> "$working_root/methods/miniprot/$mode/stderr.log"
    ) & pids+=("$!"); labels+=("miniprot:$mode")
    for method in gemoma lifton; do
        (
            "$python_bin" -m ploidypatch.cli baseline adapt-gff \
                --perturbed-gff "$annotation" \
                --candidate-gff "$working_root/merged/$method.gff3" \
                --source "$method" --max-existing-cds-overlap 0.2 \
                --max-redundancy-overlap 0.5 \
                --output-gff "$working_root/methods/$method/$mode/candidate.gff3" \
                --decisions-tsv "$working_root/methods/$method/$mode/decisions.tsv" \
                > "$working_root/methods/$method/$mode/stdout.json" \
                2> "$working_root/methods/$method/$mode/stderr.log"
        ) & pids+=("$!"); labels+=("$method:$mode")
    done
done
failed=0
for i in "${!pids[@]}"; do
    if ! wait "${pids[$i]}"; then echo "failed apple adapt: ${labels[$i]}" >&2; failed=1; fi
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
    for mode in blind complete_control; do
        base=$blind_gff; [[ $mode == complete_control ]] && base=$source_gff
        (
            "$python_bin" -m ploidypatch.cli baseline select-method-consensus \
                --base-gff "$base" \
                --candidate "miniprot=$working_root/methods/miniprot/$mode/candidate.gff3" \
                --candidate "gemoma=$working_root/methods/gemoma/$mode/candidate.gff3" \
                --candidate "lifton=$working_root/methods/lifton/$mode/candidate.gff3" \
                --min-method-support "$support" --max-redundancy-overlap 0.5 \
                --redundancy-policy "$redundancy" \
                --output-gff "$working_root/consensus/$pool/$mode/candidate.gff3" \
                --decisions-tsv "$working_root/consensus/$pool/$mode/decisions.tsv" \
                > "$working_root/consensus/$pool/$mode/stdout.json" \
                2> "$working_root/consensus/$pool/$mode/stderr.log"
        ) & pids+=("$!"); labels+=("$pool:$mode")
    done
done
failed=0
for i in "${!pids[@]}"; do
    if ! wait "${pids[$i]}"; then echo "failed apple pool: ${labels[$i]}" >&2; failed=1; fi
done
[[ $failed -eq 0 ]] || exit 1

{
    printf 'role\tbytes\tsha256\tpath\n'
    find "$working_root/methods" "$working_root/consensus" -type f \
        \( -name candidate.gff3 -o -name decisions.tsv \) -print0 | sort -z \
        | while IFS= read -r -d '' path; do
            role=${path#"$working_root/"}; role=${role//\//_}
            printf '%s\t%s\t%s\t%s\n' "$role" "$(stat -Lc %s "$path")" \
                "$(sha256sum "$path" | awk '{print $1}')" "$path"
        done
} > "$working_root/freeze/candidate_freeze.tsv"
(
    cd "$working_root"
    find . -type f \( -name '*.gff3' -o -name '*.json' -o -name '*.tsv' \) \
        -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS
    sha256sum -c SHA256SUMS >/dev/null
)
du -sb "$working_root" > "$working_root/disk_bytes.txt"
mv "$working_root" "$result_root"
printf 'apple truth-blind method pools frozen: %s\n' "$result_root"
