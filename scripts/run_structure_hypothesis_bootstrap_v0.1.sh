#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "usage: $0 PROJECT_ROOT" >&2
    exit 2
fi

project_root=$(realpath "$1")
code_root=$project_root/code
python_bin=$project_root/envs/ploidypatch-dev/bin/python
result_root=$project_root/results/structure_hypotheses/bootstrap_v0.1
working_root=${result_root}.working
events=(
    annotation_boundary_shift
    annotation_fused_gene
    annotation_missing_internal_exon
    annotation_split_gene
)
ath_root=$project_root/results/structure_hypotheses/evaluation/ath_tair10_v0.1
rice_cultivar_root=$project_root/results/heldout_structure/osa_irgsp10_v0.1
rice_phylo_root=$project_root/results/heldout_structure/osa_irgsp10_phylo_grouped_v0.1

if [[ -e $result_root || -e $working_root ]]; then
    echo "refusing to overwrite structure bootstrap: $result_root" >&2
    exit 1
fi
for cohort_root in "$ath_root" "$rice_cultivar_root" "$rice_phylo_root"; do
    for support in 1 2; do
        for event in "${events[@]}"; do
            required=$cohort_root/source_support_${support}/$event/score.json
            if [[ ! -s $required ]]; then
                echo "missing or empty structure score: $required" >&2
                exit 1
            fi
        done
    done
done

mkdir -p "$working_root"
{
    printf 'field\tvalue\n'
    printf 'code_commit\t%s\n' "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'replicates\t10000\n'
    printf 'seed\t20260816\n'
    printf 'alpha\t0.05\n'
    printf 'resampling_unit\tevent\n'
    printf 'stratified_by_event_type\ttrue\n'
} > "$working_root/run_contract.tsv"
{
    printf 'cohort\ttier\tevent\tbytes\tsha256\tpath\n'
    for cohort_entry in \
        "ath_tair10:$ath_root" \
        "osa_cultivar_ungrouped:$rice_cultivar_root" \
        "osa_phylo_grouped:$rice_phylo_root"; do
        cohort=${cohort_entry%%:*}
        cohort_root=${cohort_entry#*:}
        for support in 1 2; do
            tier=tier_b
            if [[ $support == 2 ]]; then tier=tier_a; fi
            for event in "${events[@]}"; do
                path=$cohort_root/source_support_${support}/$event/score.json
                printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
                    "$cohort" "$tier" "$event" "$(stat -Lc %s "$path")" \
                    "$(sha256sum "$path" | awk '{print $1}')" "$path"
            done
        done
    done
} > "$working_root/input_manifest.tsv"

cd "$code_root"
for cohort_entry in \
    "ath_tair10:$ath_root" \
    "osa_cultivar_ungrouped:$rice_cultivar_root" \
    "osa_phylo_grouped:$rice_phylo_root"; do
    cohort=${cohort_entry%%:*}
    cohort_root=${cohort_entry#*:}
    score_args=()
    for support in 1 2; do
        tier=tier_b
        if [[ $support == 2 ]]; then tier=tier_a; fi
        for event in "${events[@]}"; do
            score_args+=(
                --score "$tier=$cohort_root/source_support_${support}/$event/score.json"
            )
        done
    done
    "$python_bin" -m ploidypatch.cli benchmark bootstrap-structure-hypotheses \
        "${score_args[@]}" \
        --replicates 10000 \
        --seed 20260816 \
        --alpha 0.05 \
        --output-json "$working_root/$cohort.json" \
        > "$working_root/$cohort.stdout.json" \
        2> "$working_root/$cohort.stderr.log"
done

for output in \
    "$working_root/ath_tair10.json" \
    "$working_root/osa_cultivar_ungrouped.json" \
    "$working_root/osa_phylo_grouped.json"; do
    if [[ ! -s $output ]]; then
        echo "missing or empty structure bootstrap output: $output" >&2
        exit 1
    fi
done
(
    cd "$working_root"
    sha256sum ./*.json ./*.tsv > SHA256SUMS
)
du -sb "$working_root" > "$working_root/disk_bytes.txt"
mv "$working_root" "$result_root"
printf 'structure-hypothesis bootstraps completed: %s\n' "$result_root"
