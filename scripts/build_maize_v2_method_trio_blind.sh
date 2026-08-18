#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then echo "usage: $0 PROJECT_ROOT" >&2; exit 2; fi
project_root=$(realpath "$1")
code_root=$project_root/code
python_bin=$project_root/envs/ploidypatch-dev/bin/python
benchmark=$project_root/benchmark/structure/copy_collapse_v0.2/zma_maize1/annotation_copy_collapse_seed20260829
blind_gff=$benchmark/blind/perturbed.gff3
upstream=$project_root/results/baselines/maize_v2
policy=$code_root/config/maize_v2_zero_retuning_policy.tsv
result_root=$project_root/results/copy_collapse/holdout/maize_v2_method_trio
working_root=${result_root}.working
declare -A raw=(
    [gemoma_sorghum]="$upstream/gemoma/sorghum_bicolor/upstream/final_annotation.gff"
    [gemoma_setaria]="$upstream/gemoma/setaria_italica/upstream/final_annotation.gff"
    [lifton_sorghum]="$upstream/lifton/sorghum_bicolor/upstream/lifton.gff3"
    [lifton_setaria]="$upstream/lifton/setaria_italica/upstream/lifton.gff3"
)
miniprot_gff=$upstream/miniprot/raw/miniprot.gff3
protein_map=$upstream/miniprot/reference/maize_outgroups.map.tsv

for required in "$python_bin" "$blind_gff" "$policy" "$miniprot_gff" \
                "$protein_map" "${raw[@]}"; do
    [[ -s $required ]] || { echo "missing maize blind method input: $required" >&2; exit 1; }
done
if [[ -e $result_root || -e $working_root ]]; then
    echo "refusing to overwrite maize blind method trio" >&2; exit 1
fi
mkdir -p "$working_root/merged" "$working_root/methods"/{miniprot,gemoma,lifton}/blind \
    "$working_root/consensus"/{union,support2,support3}/blind
{
    printf 'field\tvalue\ncode_commit\t%s\n' \
        "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'split\texternal_zero_retuning_holdout_v0.2\n'
    printf 'candidate_truth_access\tfalse\ntarget_complete_annotation_access\tfalse\n'
    printf 'method_families\tminiprot,gemoma,lifton\n'
    printf 'references\tSorghum_bicolor,Setaria_italica\n'
    printf 'within_method_reference_vote_count\t1\n'
    printf 'candidate_universe\texact_phased_CDS_union\n'
    printf 'policy_sha256\t%s\n' "$(sha256sum "$policy" | awk '{print $1}')"
    printf 'consensus_module_sha256\t%s\n' \
        "$(sha256sum "$code_root/src/ploidypatch/consensus.py" | awk '{print $1}')"
} > "$working_root/run_contract.tsv"
{
    printf 'role\tbytes\tsha256\tpath\n'
    printf 'blind_gff\t%s\t%s\t%s\n' "$(stat -Lc %s "$blind_gff")" \
        "$(sha256sum "$blind_gff" | awk '{print $1}')" "$blind_gff"
    for key in gemoma_sorghum gemoma_setaria lifton_sorghum lifton_setaria; do
        path=${raw[$key]}
        printf '%s\t%s\t%s\t%s\n' "$key" "$(stat -Lc %s "$path")" \
            "$(sha256sum "$path" | awk '{print $1}')" "$path"
    done
    for entry in "miniprot:$miniprot_gff" "protein_map:$protein_map"; do
        role=${entry%%:*}; path=${entry#*:}
        printf '%s\t%s\t%s\t%s\n' "$role" "$(stat -Lc %s "$path")" \
            "$(sha256sum "$path" | awk '{print $1}')" "$path"
    done
} > "$working_root/input_manifest.tsv"

cd "$code_root"
for method in gemoma lifton; do
    "$python_bin" -m ploidypatch.cli baseline merge-candidate-gffs \
        --candidate "sorghum=${raw[${method}_sorghum]}" \
        --candidate "setaria=${raw[${method}_setaria]}" \
        --output-gff "$working_root/merged/$method.gff3" \
        --provenance-tsv "$working_root/merged/$method.provenance.tsv" \
        > "$working_root/merged/$method.stdout.json" \
        2> "$working_root/merged/$method.stderr.log"
done

pids=(); labels=()
(
    "$python_bin" -m ploidypatch.cli baseline adapt-miniprot \
        --perturbed-gff "$blind_gff" --miniprot-gff "$miniprot_gff" \
        --protein-map "$protein_map" \
        --output-gff "$working_root/methods/miniprot/blind/candidate.gff3" \
        --decisions-tsv "$working_root/methods/miniprot/blind/decisions.tsv" \
        > "$working_root/methods/miniprot/blind/stdout.json" \
        2> "$working_root/methods/miniprot/blind/stderr.log"
) & pids+=("$!"); labels+=("miniprot")
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
    if ! wait "${pids[$i]}"; then echo "failed blind adapt: ${labels[$i]}" >&2; failed=1; fi
done
[[ $failed -eq 0 ]] || exit 1

pids=(); labels=()
for tier in union support2 support3; do
    case $tier in union) support=1 ;; support2) support=2 ;; support3) support=3 ;; esac
    (
        "$python_bin" -m ploidypatch.cli baseline select-method-consensus \
            --base-gff "$blind_gff" \
            --candidate "miniprot=$working_root/methods/miniprot/blind/candidate.gff3" \
            --candidate "gemoma=$working_root/methods/gemoma/blind/candidate.gff3" \
            --candidate "lifton=$working_root/methods/lifton/blind/candidate.gff3" \
            --min-method-support "$support" --max-redundancy-overlap 0.5 \
            --output-gff "$working_root/consensus/$tier/blind/candidate.gff3" \
            --decisions-tsv "$working_root/consensus/$tier/blind/decisions.tsv" \
            > "$working_root/consensus/$tier/blind/stdout.json" \
            2> "$working_root/consensus/$tier/blind/stderr.log"
    ) & pids+=("$!"); labels+=("$tier")
done
failed=0
for i in "${!pids[@]}"; do
    if ! wait "${pids[$i]}"; then echo "failed blind consensus: ${labels[$i]}" >&2; failed=1; fi
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
} > "$working_root/candidate_freeze.tsv"
(
    cd "$working_root"
    find . -type f \( -name '*.gff3' -o -name '*.json' -o -name '*.tsv' \) \
        -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS
    sha256sum -c SHA256SUMS >/dev/null
)
du -sb "$working_root" > "$working_root/disk_bytes.txt"
mv "$working_root" "$result_root"
printf 'maize blind method trio frozen: %s\n' "$result_root"
