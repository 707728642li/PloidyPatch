#!/usr/bin/env bash
set -euo pipefail
umask 077

[[ $# -eq 1 ]] || { echo "usage: $0 PROJECT_ROOT" >&2; exit 64; }
project_root=$(realpath "$1")
code_root=$project_root/code
dev_python=$project_root/envs/ploidypatch-dev/bin/python
miniprot=$project_root/envs/ploidypatch-baseline/bin/miniprot
normalized=$project_root/results/baselines/walnut/v0.8/normalized
baseline_root=$project_root/results/baselines/walnut/v0.8
target=$normalized/target_walnut2/primary_chromosomes.genome.fa

[[ ${PLOIDYPATCH_BLIND_RUNNER:-} == 1 ]] || { echo "Walnut candidate methods require blind runner" >&2; exit 65; }
[[ ${PLOIDYPATCH_NETWORK_ACCESS:-} == none ]] || { echo "Walnut candidate methods require network=none" >&2; exit 66; }
[[ -x $dev_python && -x $miniprot && -s $target && -s ${target}.fai ]] || exit 67
for forbidden in /nas_data /holdout/evaluator_only /holdout/target_complete /holdout/truth /holdout/labels; do
    [[ ! -e $forbidden ]] || { echo "forbidden Walnut blind path visible: $forbidden" >&2; exit 68; }
done

declare -A reference_id=(
    [candidate_mandshurica]=juglans_mandshurica
    [candidate_carya]=carya_illinoinensis
)
for bundle in candidate_mandshurica candidate_carya; do
    for item in primary_chromosomes.genome.fa primary_chromosomes.gff3 \
        primary_chromosomes.lifton.gff3 provider.protein.fa; do
        [[ -s $normalized/$bundle/$item ]] || { echo "missing Walnut candidate input: $bundle/$item" >&2; exit 69; }
    done
done

index_root=$baseline_root/miniprot_index
index_working=${index_root}.working
[[ ! -e $index_root && ! -L $index_root && ! -e $index_working && ! -L $index_working ]] || exit 70
mkdir -p "$index_working"
/usr/bin/time -v -o "$index_working/resource.time.txt" \
    "$miniprot" -t32 -d "$index_working/walnut2.mpi" "$target" \
    > "$index_working/stdout.log" 2> "$index_working/stderr.log"
(cd "$index_working" && find . -type f ! -name SHA256SUMS -printf '%P\0' | sort -z | xargs -0 sha256sum > SHA256SUMS && sha256sum -c SHA256SUMS >/dev/null)
mv "$index_working" "$index_root"

run_miniprot() {
    local bundle=$1
    local output=$baseline_root/miniprot/$bundle
    local working=${output}.working
    [[ ! -e $output && ! -L $output && ! -e $working && ! -L $working ]] || return 71
    mkdir -p "$working"/{reference,raw,logs}
    PYTHONPATH="$code_root/src" "$dev_python" -m ploidypatch.cli baseline prepare-proteins \
        --protein "$bundle=$normalized/$bundle/provider.protein.fa" \
        --output-fasta "$working/reference/protein.fa" \
        --output-map "$working/reference/protein.map.tsv" \
        > "$working/reference/prepare.stdout.json" 2> "$working/reference/prepare.stderr.log"
    /usr/bin/time -v -o "$working/logs/projection.time.txt" \
        "$miniprot" -I -t32 --gff-only -N4 --outn 4 --outs 0.8 --outc 0.5 \
        "$index_root/walnut2.mpi" "$working/reference/protein.fa" \
        > "$working/raw/miniprot.gff3" 2> "$working/logs/projection.stderr.log"
    grep -q '^##gff-version' "$working/raw/miniprot.gff3"
    {
        printf 'field\tvalue\nreference_bundle\t%s\nthreads\t32\n' "$bundle"
        printf 'target_genome_sha256\t%s\n' "$(sha256sum "$target" | awk '{print $1}')"
        printf 'reference_protein_sha256\t%s\ntruth_access\tfalse\nranker_access\tfalse\n' \
            "$(sha256sum "$normalized/$bundle/provider.protein.fa" | awk '{print $1}')"
    } > "$working/run_contract.tsv"
    (cd "$working" && find . -type f ! -name SHA256SUMS -printf '%P\0' | sort -z | xargs -0 sha256sum > SHA256SUMS && sha256sum -c SHA256SUMS >/dev/null)
    mv "$working" "$output"
}

pids=(); labels=()
for bundle in candidate_mandshurica candidate_carya; do
    (run_miniprot "$bundle") & pids+=("$!"); labels+=("miniprot:$bundle")
    (
        bash "$code_root/scripts/run_gemoma_homology.sh" \
            "$project_root/envs/ploidypatch-gemoma" "$baseline_root/gemoma/$bundle" \
            "$target" "$normalized/$bundle/primary_chromosomes.genome.fa" \
            "$normalized/$bundle/primary_chromosomes.gff3" "${reference_id[$bundle]}" 32
    ) & pids+=("$!"); labels+=("gemoma:$bundle")
    (
        bash "$code_root/scripts/run_lifton_transfer.sh" \
            "$project_root/envs/ploidypatch-lifton" "$baseline_root/lifton/$bundle" \
            "$target" "$normalized/$bundle/primary_chromosomes.genome.fa" \
            "$normalized/$bundle/primary_chromosomes.lifton.gff3" "${reference_id[$bundle]}" 16
    ) & pids+=("$!"); labels+=("lifton:$bundle")
done
failed=0
for index in "${!pids[@]}"; do
    if ! wait "${pids[$index]}"; then
        printf 'Walnut candidate method failed: %s\n' "${labels[$index]}" >&2
        failed=1
    fi
done
[[ $failed -eq 0 ]] || exit 72
for method in miniprot gemoma lifton; do
    for bundle in candidate_mandshurica candidate_carya; do
        method_root=$baseline_root/$method/$bundle
        if [[ $method == miniprot ]]; then
            [[ -s $method_root/SHA256SUMS ]] || exit 73
            (cd "$method_root" && sha256sum -c SHA256SUMS >/dev/null)
        elif [[ $method == gemoma ]]; then
            PYTHONPATH="$code_root/src" "$dev_python" \
                "$code_root/scripts/verify_published_method_output_v0.8.py" \
                "$method_root" final_annotation.gff \
                unfiltered_predictions_from_species_0.gff \
                reference_gene_table.tabular predicted_proteins.fasta \
                protocol_GeMoMaPipeline.txt
        else
            PYTHONPATH="$code_root/src" "$dev_python" \
                "$code_root/scripts/verify_published_method_output_v0.8.py" \
                "$method_root" lifton.gff3
        fi
    done
done
