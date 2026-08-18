#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "usage: $0 PROJECT_ROOT" >&2
    exit 2
fi

project_root=$(realpath "$1")
code_root=$project_root/code
python_bin=$project_root/envs/ploidypatch-dev/bin/python
miniprot_bin=$project_root/envs/ploidypatch-baseline/bin/miniprot
public_root=$project_root/data/public/ensembl_plants_62
source_gff=$project_root/data/derived/structure_sources_v0.1/ath_tair10/source.gff3
target_genome=$public_root/arabidopsis_thaliana/genome/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
aly_protein=$public_root/arabidopsis_lyrata/annotation/Arabidopsis_lyrata.v.1.0.pep.all.fa.gz
aha_protein=$public_root/arabidopsis_halleri/annotation/Arabidopsis_halleri.Ahal2.2.pep.all.fa.gz
benchmark_root=$project_root/benchmark/public_models_missing_gene/v0.1/ath_tair10_seed20260814
benchmark_working=${benchmark_root}.working
result_root=$project_root/results/baselines/miniprot_v0.18/ath_tair10_missing_gene_v0.1
result_working=${result_root}.working

for required in "$python_bin" "$miniprot_bin" "$source_gff" "$target_genome" \
                "$aly_protein" "$aha_protein"; do
    if [[ ! -s $required ]]; then
        echo "missing Arabidopsis benchmark prerequisite: $required" >&2
        exit 1
    fi
done
for target in "$benchmark_root" "$benchmark_working" \
              "$result_root" "$result_working"; do
    if [[ -e $target ]]; then
        echo "refusing to overwrite Arabidopsis benchmark artifact: $target" >&2
        exit 1
    fi
done

mkdir -p "$benchmark_working/blind" "$benchmark_working/evaluator/truth"
cd "$code_root"
"$python_bin" -m ploidypatch.cli benchmark perturb \
    --gff "$source_gff" \
    --output-dir "$benchmark_working/blind" \
    --truth-dir "$benchmark_working/evaluator/truth" \
    --event-type annotation_missing_gene \
    --count 800 \
    --seed 20260814 \
    > "$benchmark_working/perturb.stdout.json" \
    2> "$benchmark_working/perturb.stderr.log"
"$python_bin" -m ploidypatch.cli benchmark restore \
    --perturbed-gff "$benchmark_working/blind/perturbed.gff3" \
    --truth "$benchmark_working/evaluator/truth/hidden_truth.json" \
    --output-gff "$benchmark_working/evaluator/restored.gff3" \
    > "$benchmark_working/evaluator/restoration_report.json" \
    2> "$benchmark_working/evaluator/restoration.stderr.log"
if ! cmp -s "$source_gff" "$benchmark_working/evaluator/restored.gff3"; then
    echo "Arabidopsis missing-gene restoration is not byte-identical" >&2
    exit 1
fi
{
    printf 'field\tvalue\n'
    printf 'code_commit\t%s\n' "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'dataset\tath_tair10\n'
    printf 'event_type\tannotation_missing_gene\n'
    printf 'events\t800\n'
    printf 'seed\t20260814\n'
    printf 'role\tpublic_model_development\n'
} > "$benchmark_working/run_contract.tsv"
(
    cd "$benchmark_working"
    find . -type f \( -name '*.gff3' -o -name '*.json' -o -name '*.tsv' \) \
        -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS
)
mv "$benchmark_working" "$benchmark_root"

mkdir -p "$result_working/index" "$result_working/projection" \
    "$result_working/blind" "$result_working/complete_control"
{
    printf 'field\tvalue\n'
    printf 'code_commit\t%s\n' "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'target\tArabidopsis_thaliana_TAIR10\n'
    printf 'references\tArabidopsis_lyrata_v1,Arabidopsis_halleri_Ahal2.2\n'
    printf 'target_derived_reference_proteins\tfalse\n'
    printf 'candidate_generation_used_hidden_truth\tfalse\n'
    printf 'min_identity\t0.5\n'
    printf 'min_query_coverage\t0.5\n'
    printf 'max_existing_cds_overlap\t0.2\n'
    printf 'max_redundancy_overlap\t0.5\n'
} > "$result_working/run_contract.tsv"
{
    printf 'role\tbytes\tsha256\tpath\n'
    for entry in \
        "source_gff:$source_gff" \
        "target_genome:$target_genome" \
        "aly_protein:$aly_protein" \
        "aha_protein:$aha_protein"; do
        role=${entry%%:*}
        path=${entry#*:}
        printf '%s\t%s\t%s\t%s\n' "$role" "$(stat -Lc %s "$path")" \
            "$(sha256sum "$path" | awk '{print $1}')" "$path"
    done
} > "$result_working/input_manifest.tsv"

"$python_bin" -m ploidypatch.cli baseline prepare-proteins \
    --protein "aly=$aly_protein" \
    --protein "aha=$aha_protein" \
    --output-fasta "$result_working/projection/references.protein.fa" \
    --output-map "$result_working/projection/references.protein_map.tsv" \
    > "$result_working/projection/prepare.stdout.json" \
    2> "$result_working/projection/prepare.stderr.log"
/usr/bin/time -v -o "$result_working/index/resource.time.txt" \
    "$miniprot_bin" -t64 -d "$result_working/index/ath_tair10.mpi" \
        "$target_genome" \
        > "$result_working/index/stdout.log" \
        2> "$result_working/index/stderr.log"
/usr/bin/time -v -o "$result_working/projection/resource.time.txt" \
    "$miniprot_bin" -I -t64 --gff-only -N4 --outn 4 --outs 0.8 --outc 0.5 \
        "$result_working/index/ath_tair10.mpi" \
        "$result_working/projection/references.protein.fa" \
        > "$result_working/projection/miniprot.gff3" \
        2> "$result_working/projection/miniprot.stderr.log"
if ! grep -q '^##gff-version' "$result_working/projection/miniprot.gff3"; then
    echo "Arabidopsis miniprot output lacks GFF sentinel" >&2
    exit 1
fi

for mode in blind complete_control; do
    if [[ $mode == blind ]]; then
        annotation=$benchmark_root/blind/perturbed.gff3
    else
        annotation=$benchmark_root/evaluator/restored.gff3
    fi
    /usr/bin/time -v -o "$result_working/$mode/resource.time.txt" \
        "$python_bin" -m ploidypatch.cli baseline adapt-miniprot \
            --perturbed-gff "$annotation" \
            --miniprot-gff "$result_working/projection/miniprot.gff3" \
            --protein-map "$result_working/projection/references.protein_map.tsv" \
            --output-gff "$result_working/$mode/candidate.gff3" \
            --decisions-tsv "$result_working/$mode/decisions.tsv" \
            > "$result_working/$mode/stdout.json" \
            2> "$result_working/$mode/stderr.log"
done
/usr/bin/time -v -o "$result_working/score.resource.time.txt" \
    "$python_bin" -m ploidypatch.cli benchmark score \
        --source-gff "$source_gff" \
        --perturbed-gff "$benchmark_root/blind/perturbed.gff3" \
        --candidate-gff "$result_working/blind/candidate.gff3" \
        --control-candidate-gff "$result_working/complete_control/candidate.gff3" \
        --truth "$benchmark_root/evaluator/truth/hidden_truth.json" \
        --include-event-details \
        > "$result_working/score.json" \
        2> "$result_working/score.stderr.log"

for output in \
    "$result_working/projection/miniprot.gff3" \
    "$result_working/blind/candidate.gff3" \
    "$result_working/complete_control/candidate.gff3" \
    "$result_working/score.json"; do
    if [[ ! -s $output ]]; then
        echo "Arabidopsis miniprot output is missing or empty: $output" >&2
        exit 1
    fi
done
(
    cd "$result_working"
    find . -type f \( -name '*.gff3' -o -name '*.json' -o -name '*.tsv' \) \
        -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS
)
du -sb "$result_working" > "$result_working/disk_bytes.txt"
mv "$result_working" "$result_root"
printf 'Arabidopsis missing-gene miniprot baseline completed: %s\n' "$result_root"
