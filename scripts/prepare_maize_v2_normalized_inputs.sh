#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then echo "usage: $0 PROJECT_ROOT" >&2; exit 2; fi
project_root=$(realpath "$1")
code_root=$project_root/code
python_bin=$project_root/envs/ploidypatch-dev/bin/python
parallel_bin=/data/codexli/software/conda/miniforge3/bin/parallel
public_root=$project_root/data/public/maize_v2_holdout
result_root=$project_root/data/derived/holdout_inputs/maize_v2
working_root=${result_root}.working

species=(zea_mays sorghum_bicolor setaria_italica)
declare -A primary=(
    [zea_mays]="$code_root/config/primary_seqids/zea_mays_ensembl62.tsv"
    [sorghum_bicolor]="$code_root/config/primary_seqids/sorghum_bicolor_ensembl62.tsv"
    [setaria_italica]="$code_root/config/primary_seqids/setaria_italica_ensembl62.tsv"
)
declare -A expected=(
    [zea_mays]=10 [sorghum_bicolor]=10 [setaria_italica]=9
)

[[ -x $python_bin ]] || { echo "missing project Python" >&2; exit 1; }
[[ -x $parallel_bin ]] || { echo "missing GNU parallel" >&2; exit 1; }
for name in "${species[@]}"; do
    for required in "$public_root/$name/genome.fa.gz" \
                    "$public_root/$name/annotation.gff3.gz" \
                    "${primary[$name]}"; do
        [[ -s $required ]] || { echo "missing normalization input: $required" >&2; exit 1; }
    done
done
if [[ -e $result_root || -e $working_root ]]; then
    echo "refusing to overwrite normalized maize holdout inputs" >&2
    exit 1
fi
mkdir -p "$working_root"

{
    printf 'field\tvalue\n'
    printf 'code_commit\t%s\n' "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'status\tpretruth_format_normalization_only\n'
    printf 'hidden_benchmark_generated\tfalse\n'
    printf 'target_homeolog_pairs_enumerated\tfalse\n'
    printf 'maize_label_access\tfalse\n'
    printf 'target_species\tZea_mays_B73_NAM5\n'
    printf 'reference_species\tSorghum_bicolor_NCBIv3,Setaria_italica_v2\n'
    printf 'normalization_module_sha256\t%s\n' \
        "$(sha256sum "$code_root/src/ploidypatch/normalize.py" | awk '{print $1}')"
} > "$working_root/run_contract.tsv"
{
    printf 'species\trole\tbytes\tsha256\tpath\n'
    for name in "${species[@]}"; do
        for role in genome annotation primary; do
            case $role in
                genome) path=$public_root/$name/genome.fa.gz ;;
                annotation) path=$public_root/$name/annotation.gff3.gz ;;
                primary) path=${primary[$name]} ;;
            esac
            printf '%s\t%s\t%s\t%s\t%s\n' "$name" "$role" \
                "$(stat -Lc %s "$path")" \
                "$(sha256sum "$path" | awk '{print $1}')" "$path"
        done
    done
} > "$working_root/input_manifest.tsv"

commands=$working_root/normalization_commands.tsv
{
    for name in "${species[@]}"; do
        printf '%s\t%s\t%s\t%s\n' "$name" \
            "$public_root/$name/annotation.gff3.gz" \
            "$public_root/$name/genome.fa.gz" "${primary[$name]}"
    done
} > "$commands"

cd "$code_root"
"$parallel_bin" --colsep '\t' --delay 1 -j 3 \
    --joblog "$working_root/parallel.joblog.tsv" \
    '/usr/bin/time -v -o '"$working_root"'/{1}.resource.time.txt '"$python_bin"' -m ploidypatch.cli normalize primary-annotation --gff {2} --genome {3} --primary-seqid-table {4} --output-dir '"$working_root"'/{1} > '"$working_root"'/{1}.stdout.json 2> '"$working_root"'/{1}.stderr.log' \
    :::: "$commands"

for name in "${species[@]}"; do
    for required in \
        "$working_root/$name/primary_chromosomes.genome.fa" \
        "$working_root/$name/primary_chromosomes.genome.fa.fai" \
        "$working_root/$name/primary_chromosomes.gff3" \
        "$working_root/$name/manifest.json"; do
        [[ -s $required ]] || { echo "missing normalized artifact: $required" >&2; exit 1; }
    done
    observed=$(wc -l < "$working_root/$name/primary_chromosomes.genome.fa.fai")
    [[ $observed -eq ${expected[$name]} ]] || {
        echo "$name primary chromosome count $observed is not ${expected[$name]}" >&2
        exit 1
    }
done
(
    cd "$working_root"
    find . -type f \( -name '*.fa' -o -name '*.fai' -o -name '*.gff3' \
        -o -name '*.json' -o -name '*.tsv' \) -print0 \
        | sort -z | xargs -0 sha256sum > SHA256SUMS
    sha256sum -c SHA256SUMS >/dev/null
)
du -sb "$working_root" > "$working_root/disk_bytes.txt"
mv "$working_root" "$result_root"
printf 'normalized maize v0.2 holdout inputs frozen: %s\n' "$result_root"
