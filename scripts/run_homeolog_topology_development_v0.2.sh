#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then echo "usage: $0 PROJECT_ROOT" >&2; exit 2; fi
project_root=$(realpath "$1")
code_root=$project_root/code
python_bin=$project_root/envs/ploidypatch-dev/bin/python
result_root=$project_root/results/copy_collapse/model_development/homeolog_topology_v0.2
working_root=${result_root}.working

if [[ -e $result_root || -e $working_root ]]; then
    echo "refusing to overwrite homeolog-topology development result" >&2
    exit 1
fi
mkdir -p "$working_root"

declare -A features selections candidates bases
features[glycine]=$project_root/results/copy_collapse/model_development/glycine_feature_matrix_v0.1/blind/features.tsv
selections[glycine]=$project_root/results/copy_collapse/wgd_reanchor_union_glycine_dev_v0.1/self_wgdi/blind/selected/selection.tsv
candidates[glycine]=$project_root/results/copy_collapse/consensus_union_glycine_dev_v0.1/blind/candidate.gff3
bases[glycine]=$project_root/benchmark/structure/copy_collapse_v0.1/gma_v21/annotation_copy_collapse_seed20260815/blind/perturbed.gff3

features[brassica]=$project_root/results/copy_collapse/model_development/brassica_glycine_model_transfer_v0.1/blind/features.tsv
selections[brassica]=$project_root/results/copy_collapse/model_development/brassica_union_self_wgd_v0.1/blind/selected/selection.tsv
candidates[brassica]=$project_root/results/copy_collapse/model_development/brassica_method_trio_v0.1/consensus/union/blind/candidate.gff3
bases[brassica]=$project_root/benchmark/structure/copy_collapse_v0.1/bna_daae/annotation_copy_collapse_seed20260814/blind/perturbed.gff3

features[cotton]=$project_root/results/copy_collapse/holdout/cotton_glycine_model_transfer_v0.1/blind/features.tsv
selections[cotton]=$project_root/results/copy_collapse/holdout/cotton_union_self_wgd_v0.1/blind/selected/selection.tsv
candidates[cotton]=$project_root/results/copy_collapse/holdout/cotton_method_trio_v0.1/consensus/union/blind/candidate.gff3
bases[cotton]=$project_root/benchmark/structure/copy_collapse_v0.1/ghi_ad/annotation_copy_collapse_seed20260817/blind/perturbed.gff3

{
    printf 'field\tvalue\n'
    printf 'code_commit\t%s\n' "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'truth_access\tfalse\n'
    printf 'development_species\tGlycine_max,Brassica_napus\n'
    printf 'post_holdout_diagnostic_species\tGossypium_hirsutum\n'
    printf 'topology_policy\tbest_existing_WGD_partner_isoform_pair\n'
    printf 'parallel_species_jobs\t3\n'
} > "$working_root/run_contract.tsv"
{
    printf 'species\trole\tbytes\tsha256\tpath\n'
    for species in glycine brassica cotton; do
        for role in features selections candidates bases; do
            case $role in
                features) path=${features[$species]} ;;
                selections) path=${selections[$species]} ;;
                candidates) path=${candidates[$species]} ;;
                bases) path=${bases[$species]} ;;
            esac
            [[ -s $path ]] || { echo "missing topology input: $path" >&2; exit 1; }
            printf '%s\t%s\t%s\t%s\t%s\n' "$species" "$role" \
                "$(stat -Lc %s "$path")" "$(sha256sum "$path" | awk '{print $1}')" \
                "$path"
        done
    done
} > "$working_root/input_manifest.tsv"

cd "$code_root"
pids=(); labels=()
for species in glycine brassica cotton; do
    mkdir -p "$working_root/$species"
    /usr/bin/time -v -o "$working_root/$species/resource.time.txt" \
        "$python_bin" -m ploidypatch.cli evidence build-homeolog-topology-features \
        --copy-features "${features[$species]}" \
        --wgd-selection "${selections[$species]}" \
        --candidate-gff "${candidates[$species]}" \
        --base-gff "${bases[$species]}" \
        --output-tsv "$working_root/$species/features.tsv" \
        > "$working_root/$species/stdout.json" \
        2> "$working_root/$species/stderr.log" &
    pids+=("$!"); labels+=("$species")
done
failed=0
for index in "${!pids[@]}"; do
    if ! wait "${pids[$index]}"; then
        echo "homeolog topology failed: ${labels[$index]}" >&2
        failed=1
    fi
done
[[ $failed -eq 0 ]] || exit 1

(
    cd "$working_root"
    find . -type f \( -name '*.json' -o -name '*.tsv' \) -print0 \
        | sort -z | xargs -0 sha256sum > SHA256SUMS
    sha256sum -c SHA256SUMS >/dev/null
)
du -sb "$working_root" > "$working_root/disk_bytes.txt"
mv "$working_root" "$result_root"
printf 'homeolog topology features frozen: %s\n' "$result_root"
