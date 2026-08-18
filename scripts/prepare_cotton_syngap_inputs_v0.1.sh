#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then echo "usage: $0 PROJECT_ROOT" >&2; exit 2; fi
project_root=$(realpath "$1")
code_root=$project_root/code
python_bin=$project_root/envs/ploidypatch-dev/bin/python
syngap_env=$project_root/envs/ploidypatch-syngap
input_root=$project_root/data/derived/holdout_inputs/cotton_v0.1
candidate_root=$project_root/data/derived/candidate_inputs/cotton_holdout_v0.1
benchmark=$project_root/benchmark/structure/copy_collapse_v0.1/ghi_ad/annotation_copy_collapse_seed20260817
result_root=$project_root/data/derived/candidate_inputs/cotton_syngap_v0.1
working_root=${result_root}.working

declare -A source_gff=(
    [blind]="$benchmark/blind/perturbed.gff3"
    [complete_control]="$input_root/hirsutum/primary_chromosomes.gff3"
)
declare -A reference_gff=(
    [gar_a]="$candidate_root/gar_a_lifton_compat_v3/reference.compat.gff3"
    [gra_d]="$candidate_root/gra_d_lifton_compat_v2/reference.compat.gff3"
)
declare -A reference_genome=(
    [gar_a]="$input_root/arboreum/primary_chromosomes.genome.fa"
    [gra_d]="$input_root/raimondii/primary_chromosomes.genome.fa"
)
target_genome=$input_root/hirsutum/primary_chromosomes.genome.fa
for required in "$python_bin" "$syngap_env/bin/syngap" "$target_genome" \
                "${source_gff[@]}" "${reference_gff[@]}" \
                "${reference_genome[@]}"; do
    if [[ ! -s $required ]]; then echo "missing cotton SynGAP input: $required" >&2; exit 1; fi
done
if [[ -e $result_root || -e $working_root ]]; then
    echo "refusing to overwrite cotton SynGAP candidate inputs" >&2; exit 1
fi
mkdir -p "$working_root/target" "$working_root/preflight"
cd "$code_root"
for mode in blind complete_control; do
    mkdir -p "$working_root/target/$mode"
    gene_records=$(awk -F '\t' '!/^#/ && NF == 9 && $3 == "gene" { count++ } END { print count + 0 }' "${source_gff[$mode]}")
    if [[ $gene_records -gt 0 ]]; then
        ln -s "$(realpath "${source_gff[$mode]}")" \
            "$working_root/target/$mode/annotation.compat.gff3"
        printf 'transcript_id\tsynthesized_gene_id\n' \
            > "$working_root/target/$mode/transcript_gene_map.tsv"
        printf '{"input_gene_records": %s, "policy": "existing_gene_hierarchy_preserved"}\n' \
            "$gene_records" > "$working_root/target/$mode/stdout.json"
        : > "$working_root/target/$mode/stderr.log"
    else
        "$python_bin" scripts/synthesize_root_transcript_genes.py \
            --input-gff "${source_gff[$mode]}" \
            --output-gff "$working_root/target/$mode/annotation.compat.gff3" \
            --mapping-tsv "$working_root/target/$mode/transcript_gene_map.tsv" \
            > "$working_root/target/$mode/stdout.json" \
            2> "$working_root/target/$mode/stderr.log"
    fi
done

pids=(); labels=()
for mode in blind complete_control; do
    bash "$code_root/scripts/preflight_syngap_annotation.sh" \
        "$syngap_env" "$working_root/preflight/target_$mode" \
        "Ghi_$mode" "$target_genome" \
        "$working_root/target/$mode/annotation.compat.gff3" \
        > "$working_root/preflight/target_$mode.launcher.log" 2>&1 &
    pids+=("$!"); labels+=("target_$mode")
    sleep 5
done
for ref in gar_a gra_d; do
    bash "$code_root/scripts/preflight_syngap_annotation.sh" \
        "$syngap_env" "$working_root/preflight/$ref" "$ref" \
        "${reference_genome[$ref]}" "${reference_gff[$ref]}" \
        > "$working_root/preflight/$ref.launcher.log" 2>&1 &
    pids+=("$!"); labels+=("$ref")
    sleep 5
done
failed=0
for i in "${!pids[@]}"; do
    if ! wait "${pids[$i]}"; then
        echo "cotton SynGAP preflight failed: ${labels[$i]}" >&2
        failed=1
    fi
done
[[ $failed -eq 0 ]] || exit 1

{
    printf 'field\tvalue\n'
    printf 'code_commit\t%s\n' "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'target\tGossypium_hirsutum\n'
    printf 'references\tGossypium_arboreum,Gossypium_raimondii\n'
    printf 'target_compatibility\tpreserve_existing_genes_or_synthesize_only_when_absent\n'
    printf 'reference_compatibility\tvalidated_LiftOn_compatibility_GFFs\n'
    printf 'truth_access\tfalse\n'
} > "$working_root/run_contract.tsv"
{
    printf 'role\tbytes\tsha256\tpath\n'
    for mode in blind complete_control; do
        path=$working_root/target/$mode/annotation.compat.gff3
        printf 'target_%s\t%s\t%s\t%s\n' "$mode" "$(stat -Lc %s "$path")" \
            "$(sha256sum "$path" | awk '{print $1}')" "$path"
    done
    for ref in gar_a gra_d; do
        path=${reference_gff[$ref]}
        printf 'reference_%s\t%s\t%s\t%s\n' "$ref" "$(stat -Lc %s "$path")" \
            "$(sha256sum "$path" | awk '{print $1}')" "$path"
    done
} > "$working_root/input_manifest.tsv"
(
    cd "$working_root"
    find . -type f \( -name '*.json' -o -name '*.tsv' -o -name '*.gff3' \) \
        -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS
    sha256sum -c SHA256SUMS >/dev/null
)
du -sb "$working_root" > "$working_root/disk_bytes.txt"
mv "$working_root" "$result_root"
printf 'cotton SynGAP candidate inputs frozen: %s\n' "$result_root"
