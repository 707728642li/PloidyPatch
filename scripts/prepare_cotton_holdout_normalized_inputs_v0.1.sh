#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "usage: $0 PROJECT_ROOT" >&2
    exit 2
fi

project_root=$(realpath "$1")
code_root=$project_root/code
python_bin=$project_root/envs/ploidypatch-dev/bin/python
parallel_bin=/data/codexli/software/conda/miniforge3/bin/parallel
public_root=$project_root/data/public/pgcp_cotton_holdout_v0.1
result_root=$project_root/data/derived/holdout_inputs/cotton_v0.1
working_root=${result_root}.working

declare -A genome=(
    [hirsutum]="$public_root/gossypium_hirsutum/gossypium_hirsutum.genomic.fa.gz"
    [arboreum]="$public_root/gossypium_arboreum/gossypium_arboreum.genomic.fa.gz"
    [raimondii]="$public_root/gossypium_raimondii/gossypium_raimondii.genomic.fa.gz"
    [barbadense]="$public_root/gossypium_barbadense/gossypium_barbadense.genomic.fa.gz"
)
declare -A gff=(
    [hirsutum]="$public_root/gossypium_hirsutum/gossypium_hirsutum.genomic.gff.gz"
    [arboreum]="$public_root/gossypium_arboreum/gossypium_arboreum.genomic.gff.gz"
    [raimondii]="$public_root/gossypium_raimondii/gossypium_raimondii.genomic.gff.gz"
    [barbadense]="$public_root/gossypium_barbadense/gossypium_barbadense.genomic.gff.gz"
)
declare -A primary=(
    [hirsutum]="$code_root/config/primary_seqids/gossypium_hirsutum_pgcp_v245.tsv"
    [arboreum]="$code_root/config/primary_seqids/gossypium_arboreum_pgcp_v240.tsv"
    [raimondii]="$code_root/config/primary_seqids/gossypium_raimondii_pgcp_v247.tsv"
    [barbadense]="$code_root/config/primary_seqids/gossypium_barbadense_pgcp_v242.tsv"
)
declare -A expected_seqids=(
    [hirsutum]=26 [arboreum]=13 [raimondii]=13 [barbadense]=26
)

if [[ ! -x $parallel_bin ]]; then
    echo "missing GNU parallel in the canonical server base environment: $parallel_bin" >&2
    exit 1
fi
for species in hirsutum arboreum raimondii barbadense; do
    for required in "${genome[$species]}" "${gff[$species]}" \
                    "${primary[$species]}"; do
        if [[ ! -s $required ]]; then
            echo "missing or empty cotton input: $required" >&2
            exit 1
        fi
    done
done
if [[ -e $result_root || -e $working_root ]]; then
    echo "refusing to overwrite normalized cotton holdout inputs" >&2
    exit 1
fi
mkdir -p "$working_root"

{
    printf 'field\tvalue\n'
    printf 'code_commit\t%s\n' "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'status\tpretruth_format_normalization_only\n'
    printf 'hidden_benchmark_generated\tfalse\n'
    printf 'target_homeolog_pairs_enumerated\tfalse\n'
    printf 'candidate_threshold_access\tfalse\n'
    printf 'target_species\tGossypium_hirsutum\n'
    printf 'reference_species\tGossypium_arboreum,Gossypium_raimondii,Gossypium_barbadense\n'
    printf 'normalization_module_sha256\t%s\n' \
        "$(sha256sum "$code_root/src/ploidypatch/normalize.py" | awk '{print $1}')"
} > "$working_root/run_contract.tsv"
{
    printf 'species\trole\tbytes\tsha256\tpath\n'
    for species in hirsutum arboreum raimondii barbadense; do
        for role in genome gff primary; do
            case $role in
                genome) path=${genome[$species]} ;;
                gff) path=${gff[$species]} ;;
                primary) path=${primary[$species]} ;;
            esac
            printf '%s\t%s\t%s\t%s\t%s\n' "$species" "$role" \
                "$(stat -Lc %s "$path")" \
                "$(sha256sum "$path" | awk '{print $1}')" "$path"
        done
    done
} > "$working_root/input_manifest.tsv"

commands=$working_root/normalization_commands.tsv
{
    for species in hirsutum arboreum raimondii barbadense; do
        printf '%s\t%s\t%s\t%s\n' "$species" "${gff[$species]}" \
            "${genome[$species]}" "${primary[$species]}"
    done
} > "$commands"
cd "$code_root"
"$parallel_bin" --colsep '\t' --delay 2 -j 4 \
    --joblog "$working_root/parallel.joblog.tsv" \
    '/usr/bin/time -v -o '"$working_root"'/{1}.resource.time.txt '"$python_bin"' -m ploidypatch.cli normalize primary-annotation --gff {2} --genome {3} --primary-seqid-table {4} --output-dir '"$working_root"'/{1} > '"$working_root"'/{1}.stdout.json 2> '"$working_root"'/{1}.stderr.log' \
    :::: "$commands"

for species in hirsutum arboreum raimondii barbadense; do
    for required in \
        "$working_root/$species/primary_chromosomes.genome.fa" \
        "$working_root/$species/primary_chromosomes.genome.fa.fai" \
        "$working_root/$species/primary_chromosomes.gff3" \
        "$working_root/$species/manifest.json"; do
        if [[ ! -s $required ]]; then
            echo "missing normalized cotton artifact: $required" >&2
            exit 1
        fi
    done
    observed=$(wc -l < "$working_root/$species/primary_chromosomes.genome.fa.fai")
    if [[ $observed -ne ${expected_seqids[$species]} ]]; then
        echo "$species primary chromosome count $observed is not ${expected_seqids[$species]}" >&2
        exit 1
    fi
done
(
    cd "$working_root"
    find . -type f \( -name '*.fa' -o -name '*.fai' -o -name '*.gff3' \
        -o -name '*.json' -o -name '*.tsv' \) \
        -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS
    sha256sum -c SHA256SUMS >/dev/null
)
du -sb "$working_root" > "$working_root/disk_bytes.txt"
mv "$working_root" "$result_root"
printf 'normalized cotton holdout inputs frozen: %s\n' "$result_root"
