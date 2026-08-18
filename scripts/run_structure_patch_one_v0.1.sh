#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 3 ]]; then
    echo "usage: $0 PROJECT_ROOT ath_tair10|osa_irgsp10 annotation_fused_gene|annotation_split_gene" >&2
    exit 2
fi

project_root=$(realpath "$1")
cohort=$2
event=$3
code_root=$project_root/code
python_bin=$project_root/envs/ploidypatch-dev/bin/python

case $event in
    annotation_fused_gene|annotation_split_gene) ;;
    *)
        echo "automatic patch validation is frozen to fused/split hypotheses" >&2
        exit 2
        ;;
esac
case $cohort in
    ath_tair10)
        source_gff=$project_root/data/derived/structure_sources_v0.1/ath_tair10/source.gff3
        benchmark_root=$project_root/benchmark/structure/public_models_v0.1/ath_tair10/${event}_seed20260812
        candidate_root=$project_root/results/structure_candidates/miniprot_multisource/ath_tair10_v0.1/source_support_2/$event
        hypothesis_root=$project_root/results/structure_hypotheses/exact_topology/ath_tair10_v0.1
        hypotheses=$hypothesis_root/source_support_2/$event/hypotheses.tsv
        ;;
    osa_irgsp10)
        source_gff=$project_root/data/derived/structure_sources_v0.1/osa_irgsp10/source.gff3
        benchmark_root=$project_root/benchmark/structure/public_models_v0.1/osa_irgsp10/${event}_seed20260813
        candidate_root=$project_root/results/heldout_structure/osa_irgsp10_phylo_grouped_v0.1/source_support_2/$event
        hypotheses=$candidate_root/hypotheses.tsv
        ;;
    *)
        echo "unsupported structure-patch cohort: $cohort" >&2
        exit 2
        ;;
esac

perturbed_gff=$benchmark_root/blind/perturbed.gff3
truth=$benchmark_root/evaluator/truth/hidden_truth.json
candidate_gff=$candidate_root/candidate.gff3
result_root=$project_root/results/structure_repairs/exact_topology_v0.1/$cohort/$event
working_root=${result_root}.working
for required in "$python_bin" "$source_gff" "$perturbed_gff" "$truth" \
                "$candidate_gff" "$hypotheses"; do
    if [[ ! -s $required ]]; then
        echo "missing or empty structure-patch input: $required" >&2
        exit 1
    fi
done
if [[ -e $result_root || -e $working_root ]]; then
    echo "refusing to overwrite structure-patch result: $result_root" >&2
    exit 1
fi

mkdir -p "$working_root"
{
    printf 'field\tvalue\n'
    printf 'cohort\t%s\n' "$cohort"
    printf 'event_type\t%s\n' "$event"
    printf 'code_commit\t%s\n' "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'tier\tA\n'
    printf 'min_support_group_count\t2\n'
    printf 'automatic_patch_event_policy\tfused_and_split_only\n'
    printf 'source_hierarchy_policy\tcomplete_single_transcript_gene_only\n'
    printf 'compiler_truth_access\tfalse\n'
    printf 'evaluator_truth_access\ttrue\n'
    printf 'structure_patch_module_sha256\t%s\n' \
        "$(sha256sum "$code_root/src/ploidypatch/structure_patch.py" | awk '{print $1}')"
} > "$working_root/run_contract.tsv"
{
    printf 'role\tbytes\tsha256\tpath\n'
    for entry in \
        "source_gff:$source_gff" \
        "perturbed_gff:$perturbed_gff" \
        "hidden_truth:$truth" \
        "candidate_gff:$candidate_gff" \
        "hypotheses:$hypotheses"; do
        role=${entry%%:*}
        path=${entry#*:}
        printf '%s\t%s\t%s\t%s\n' "$role" "$(stat -Lc %s "$path")" \
            "$(sha256sum "$path" | awk '{print $1}')" "$path"
    done
} > "$working_root/input_manifest.tsv"

cd "$code_root"
"$python_bin" -m ploidypatch.cli patch compile-structure \
    --annotation-gff "$perturbed_gff" \
    --candidate-gff "$candidate_gff" \
    --hypotheses-tsv "$hypotheses" \
    --output-edits-json "$working_root/edits.json" \
    --event-type "$event" \
    --min-support-group-count 2 \
    > "$working_root/compile.stdout.json" \
    2> "$working_root/compile.stderr.log"
"$python_bin" -m ploidypatch.cli patch create \
    --source-gff "$perturbed_gff" \
    --edits-json "$working_root/edits.json" \
    --output-patch "$working_root/repair.patch.json" \
    > "$working_root/create.stdout.json" \
    2> "$working_root/create.stderr.log"
"$python_bin" -m ploidypatch.cli patch apply \
    --source-gff "$perturbed_gff" \
    --patch "$working_root/repair.patch.json" \
    --output-gff "$working_root/repaired.gff3" \
    > "$working_root/apply.stdout.json" \
    2> "$working_root/apply.stderr.log"
"$python_bin" -m ploidypatch.cli patch revert \
    --patched-gff "$working_root/repaired.gff3" \
    --patch "$working_root/repair.patch.json" \
    --output-gff "$working_root/reverted.gff3" \
    > "$working_root/revert.stdout.json" \
    2> "$working_root/revert.stderr.log"
if ! cmp -s "$perturbed_gff" "$working_root/reverted.gff3"; then
    echo "reverted structure patch is not byte-identical to blind input" >&2
    exit 1
fi
/usr/bin/time -v -o "$working_root/score.resource.time.txt" \
    "$python_bin" -m ploidypatch.cli benchmark score \
        --source-gff "$source_gff" \
        --perturbed-gff "$perturbed_gff" \
        --candidate-gff "$working_root/repaired.gff3" \
        --truth "$truth" \
        --include-event-details \
        > "$working_root/score.json" \
        2> "$working_root/score.stderr.log"
if ! grep -q '"grade": "pass"' "$working_root/score.json"; then
    echo "structure-patch score quality gate failed" >&2
    exit 1
fi

for output in "$working_root/edits.json" "$working_root/repair.patch.json" \
              "$working_root/repaired.gff3" "$working_root/reverted.gff3" \
              "$working_root/score.json"; do
    if [[ ! -s $output ]]; then
        echo "missing or empty structure-patch output: $output" >&2
        exit 1
    fi
done
(
    cd "$working_root"
    find . -type f \( -name '*.gff3' -o -name '*.json' -o -name '*.tsv' \) \
        -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS
)
du -sb "$working_root" > "$working_root/disk_bytes.txt"
mv "$working_root" "$result_root"
printf 'structure patch validated: %s\n' "$result_root"
