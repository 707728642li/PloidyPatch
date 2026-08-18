#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 6 ]]; then
    echo "usage: $0 PROJECT_ROOT DATASET_ID EVENT_TYPE COUNT SEED SPLIT" >&2
    exit 2
fi

project_root=$(realpath "$1")
dataset_id=$2
event_type=$3
count=$4
seed=$5
split=$6
code_root=$project_root/code
python_bin=$project_root/envs/ploidypatch-dev/bin/python
namespace=${PLOIDYPATCH_STRUCTURE_NAMESPACE:-v0.1}
if [[ ! $namespace =~ ^[A-Za-z0-9][A-Za-z0-9_.-]*$ ]]; then
    echo "structure namespace must be a safe identifier" >&2
    exit 2
fi

case $dataset_id in
    bna_daae)
        source_gff=$project_root/data/derived/normalized_bundles/v0.1/bna_daae_primary/primary_chromosomes.gff3
        ;;
    gma_v21)
        source_gff=$project_root/data/derived/syngap_glycine_v0.1/gma_complete/primary_chromosomes.gff3
        ;;
    ath_tair10|osa_irgsp10)
        source_gff=$project_root/data/derived/structure_sources_v0.1/$dataset_id/source.gff3
        ;;
    *)
        echo "unsupported structure-benchmark dataset: $dataset_id" >&2
        exit 2
        ;;
esac
case $event_type in
    annotation_missing_internal_exon|annotation_boundary_shift|annotation_split_gene|annotation_fused_gene)
        ;;
    *)
        echo "unsupported frozen structure event: $event_type" >&2
        exit 2
        ;;
esac
if [[ ! $count =~ ^[1-9][0-9]*$ || ! $seed =~ ^[1-9][0-9]*$ ]]; then
    echo "count and seed must be positive integers" >&2
    exit 2
fi
if [[ $split != development && $split != heldout ]]; then
    echo "split must be development or heldout" >&2
    exit 2
fi

result_root=$project_root/benchmark/structure/$namespace/$dataset_id/${event_type}_seed${seed}
working_root=${result_root}.working
for required in "$python_bin" "$source_gff"; do
    if [[ ! -s $required ]]; then
        echo "missing structure-benchmark prerequisite: $required" >&2
        exit 1
    fi
done
if [[ -e $result_root || -e $working_root ]]; then
    echo "refusing to overwrite structure benchmark: $result_root" >&2
    exit 1
fi
mkdir -p "$working_root/blind" "$working_root/evaluator/truth" \
    "$working_root/evaluator/validation"
{
    printf 'field\tvalue\n'
    printf 'code_commit\t%s\n' "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'dataset_id\t%s\n' "$dataset_id"
    printf 'event_type\t%s\n' "$event_type"
    printf 'count\t%s\n' "$count"
    printf 'seed\t%s\n' "$seed"
    printf 'split\t%s\n' "$split"
    printf 'source_gff_sha256\t%s\n' "$(sha256sum "$source_gff" | awk '{print $1}')"
    printf 'structure_module_sha256\t%s\n' \
        "$(sha256sum "$code_root/src/ploidypatch/structure_perturb.py" | awk '{print $1}')"
    printf 'score_module_sha256\t%s\n' \
        "$(sha256sum "$code_root/src/ploidypatch/score.py" | awk '{print $1}')"
} > "$working_root/run_contract.tsv"

cd "$code_root"
/usr/bin/time -v -o "$working_root/perturb.resource.time.txt" \
    "$python_bin" -m ploidypatch.cli benchmark perturb \
        --gff "$source_gff" \
        --output-dir "$working_root/blind" \
        --truth-dir "$working_root/evaluator/truth" \
        --event-type "$event_type" \
        --count "$count" \
        --seed "$seed" \
        > "$working_root/perturb.stdout.json" \
        2> "$working_root/perturb.stderr.log"

perturbed=$working_root/blind/perturbed.gff3
truth=$working_root/evaluator/truth/hidden_truth.json
restored=$working_root/evaluator/restored.gff3
/usr/bin/time -v -o "$working_root/restore.resource.time.txt" \
    "$python_bin" -m ploidypatch.cli benchmark restore \
        --perturbed-gff "$perturbed" \
        --truth "$truth" \
        --output-gff "$restored" \
        > "$working_root/evaluator/restoration_report.json" \
        2> "$working_root/evaluator/restoration.stderr.log"
if ! cmp -s "$source_gff" "$restored"; then
    echo "restored GFF is not byte-identical to source" >&2
    exit 1
fi

for mode in noop oracle; do
    if [[ $mode == noop ]]; then
        candidate=$perturbed
    else
        candidate=$source_gff
    fi
    /usr/bin/time -v -o "$working_root/evaluator/validation/${mode}.resource.time.txt" \
        "$python_bin" -m ploidypatch.cli benchmark score \
            --source-gff "$source_gff" \
            --perturbed-gff "$perturbed" \
            --candidate-gff "$candidate" \
            --truth "$truth" \
            --include-event-details \
            > "$working_root/evaluator/validation/score_${mode}.json" \
            2> "$working_root/evaluator/validation/score_${mode}.stderr.log"
done

for output in \
    "$working_root/blind/perturbed.gff3" \
    "$working_root/blind/manifest.json" \
    "$working_root/evaluator/truth/hidden_truth.json" \
    "$working_root/evaluator/restored.gff3" \
    "$working_root/evaluator/restoration_report.json" \
    "$working_root/evaluator/validation/score_noop.json" \
    "$working_root/evaluator/validation/score_oracle.json"; do
    if [[ ! -s $output ]]; then
        echo "structure-benchmark output is missing or empty: $output" >&2
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
printf 'structure benchmark completed: %s\n' "$result_root"
