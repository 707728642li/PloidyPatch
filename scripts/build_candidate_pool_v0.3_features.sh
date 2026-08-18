#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "usage: $0 PROJECT_ROOT {glycine|brassica|cotton|maize}" >&2; exit 2
fi
project_root=$(realpath "$1")
species=$2
code_root=$project_root/code
python_bin=$project_root/envs/ploidypatch-dev/bin/python
pool=$project_root/results/copy_collapse/candidate_pool_v0.3/$species
candidate_gff=$pool/blind/candidate.gff3
pool_decisions=$pool/blind/decisions.tsv

case $species in
glycine)
    split=development_labels_previously_seen
    benchmark=$project_root/benchmark/structure/copy_collapse_v0.1/gma_v21/annotation_copy_collapse_seed20260815
    base_gff=$benchmark/blind/perturbed.gff3
    miniprot_decisions=$project_root/results/copy_collapse/miniprot_glycine_v0.1/blind/decisions.tsv
    gemoma_decisions=$project_root/results/copy_collapse/gemoma_glycine_v0.1/blind/decisions.tsv
    lifton_decisions=$project_root/results/copy_collapse/lifton_glycine_v0.1/blind/decisions.tsv
    prior_wgd=$project_root/results/copy_collapse/wgd_reanchor_union_glycine_dev_v0.1/self_wgdi/blind/selected/selection.tsv
    ;;
brassica)
    split=development_labels_previously_seen
    benchmark=$project_root/benchmark/structure/copy_collapse_v0.1/bna_daae/annotation_copy_collapse_seed20260814
    base_gff=$benchmark/blind/perturbed.gff3
    method_root=$project_root/results/copy_collapse/model_development/brassica_method_trio_v0.1
    miniprot_decisions=$project_root/results/copy_collapse/miniprot_brassica_v0.1/blind/decisions.tsv
    gemoma_decisions=$method_root/methods/gemoma/blind/decisions.tsv
    lifton_decisions=$method_root/methods/lifton/blind/decisions.tsv
    prior_wgd=$project_root/results/copy_collapse/model_development/brassica_union_self_wgd_v0.1/blind/selected/selection.tsv
    ;;
cotton)
    split=retrospective_diagnostic_labels_previously_seen
    benchmark=$project_root/benchmark/structure/copy_collapse_v0.1/ghi_ad/annotation_copy_collapse_seed20260817
    base_gff=$benchmark/blind/perturbed.gff3
    method_root=$project_root/results/copy_collapse/holdout/cotton_method_trio_v0.1
    miniprot_decisions=$method_root/methods/miniprot/blind/decisions.tsv
    gemoma_decisions=$method_root/methods/gemoma/blind/decisions.tsv
    lifton_decisions=$method_root/methods/lifton/blind/decisions.tsv
    prior_wgd=$project_root/results/copy_collapse/holdout/cotton_union_self_wgd_v0.1/blind/selected/selection.tsv
    ;;
maize)
    split=posthoc_formal_holdout_diagnostic_no_selection
    benchmark=$project_root/benchmark/structure/copy_collapse_v0.2/zma_maize1/annotation_copy_collapse_seed20260829
    base_gff=$benchmark/blind/perturbed.gff3
    method_root=$project_root/results/copy_collapse/holdout/maize_v2_method_trio
    miniprot_decisions=$method_root/methods/miniprot/blind/decisions.tsv
    gemoma_decisions=$method_root/methods/gemoma/blind/decisions.tsv
    lifton_decisions=$method_root/methods/lifton/blind/decisions.tsv
    prior_wgd=$project_root/results/copy_collapse/holdout/maize_v2_union_self_wgd/blind/selected/selection.tsv
    ;;
*) echo "unsupported species: $species" >&2; exit 2 ;;
esac
truth=$benchmark/evaluator/truth/hidden_truth.json
result_root=$project_root/results/copy_collapse/candidate_pool_v0.3_features/$species
working_root=${result_root}.working
for required in "$python_bin" "$pool/SHA256SUMS" "$candidate_gff" \
    "$pool_decisions" "$base_gff" "$truth" "$miniprot_decisions" \
    "$gemoma_decisions" "$lifton_decisions" "$prior_wgd"; do
    [[ -s $required ]] || { echo "missing v0.3 feature input: $required" >&2; exit 1; }
done
(cd "$pool" && sha256sum -c SHA256SUMS >/dev/null)
if [[ -e $result_root || -e $working_root ]]; then
    echo "refusing to overwrite v0.3 feature result: $result_root" >&2; exit 1
fi
mkdir -p "$working_root"
{
    printf 'field\tvalue\ncode_commit\t%s\n' \
        "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'species\t%s\nsplit\t%s\n' "$species" "$split"
    printf 'candidate_policy\tretain_distinct_phased_cds_chains\n'
    printf 'wgd_propagation\tunique_existing_partner_within_conflict_set\n'
    printf 'truth_access_phase\tlabeling_only_after_feature_freeze\n'
    printf 'formal_holdout_claim_allowed\tfalse\n'
} > "$working_root/run_contract.tsv"
cd "$code_root"
"$python_bin" -m ploidypatch.cli evidence propagate-wgd-conflict-partners \
    --base-gff "$base_gff" --candidate-gff "$candidate_gff" \
    --pool-decisions "$pool_decisions" --prior-wgd-selection "$prior_wgd" \
    --output-selection "$working_root/wgd_selection.tsv" \
    > "$working_root/wgd_selection.stdout.json" \
    2> "$working_root/wgd_selection.stderr.log"
"$python_bin" -m ploidypatch.cli evidence build-copy-features \
    --consensus-decisions "$pool_decisions" \
    --method-decisions "miniprot=$miniprot_decisions" \
    --method-decisions "gemoma=$gemoma_decisions" \
    --method-decisions "lifton=$lifton_decisions" \
    --wgd-selection "$working_root/wgd_selection.tsv" \
    --output-tsv "$working_root/copy_features.tsv" \
    > "$working_root/copy_features.stdout.json" \
    2> "$working_root/copy_features.stderr.log"
"$python_bin" -m ploidypatch.cli evidence build-homeolog-topology-features \
    --copy-features "$working_root/copy_features.tsv" \
    --wgd-selection "$working_root/wgd_selection.tsv" \
    --candidate-gff "$candidate_gff" --base-gff "$base_gff" \
    --output-tsv "$working_root/topology_features.tsv" \
    > "$working_root/topology_features.stdout.json" \
    2> "$working_root/topology_features.stderr.log"
{
    printf 'feature\tbytes\tsha256\n'
    for path in "$working_root/wgd_selection.tsv" "$working_root/copy_features.tsv" \
                "$working_root/topology_features.tsv"; do
        printf '%s\t%s\t%s\n' "$(basename "$path")" "$(stat -Lc %s "$path")" \
            "$(sha256sum "$path" | awk '{print $1}')"
    done
} > "$working_root/blind_feature_freeze.tsv"
"$python_bin" -m ploidypatch.cli benchmark label-copy-features \
    --features "$working_root/copy_features.tsv" --truth "$truth" \
    --output-tsv "$working_root/labeled_features.tsv" \
    > "$working_root/labeled_features.stdout.json" \
    2> "$working_root/labeled_features.stderr.log"
for output in wgd_selection.tsv copy_features.tsv topology_features.tsv \
              labeled_features.tsv blind_feature_freeze.tsv; do
    [[ -s $working_root/$output ]] || { echo "missing v0.3 feature output: $output" >&2; exit 1; }
done
(
    cd "$working_root"
    find . -type f \( -name '*.json' -o -name '*.tsv' \) -print0 \
        | sort -z | xargs -0 sha256sum > SHA256SUMS
    sha256sum -c SHA256SUMS >/dev/null
)
mv "$working_root" "$result_root"
printf 'v0.3 feature matrix frozen: %s\n' "$result_root"
