#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "usage: $0 PROJECT_ROOT" >&2
    exit 2
fi

project_root=$(realpath "$1")
code_root=$project_root/code
python_bin=$project_root/envs/ploidypatch-dev/bin/python
public_root=$project_root/data/public/ensembl_plants_62
result_root=$project_root/data/derived/multireference_development/brassica_v0.1
working_root=${result_root}.working

declare -A genome=(
    [brapa]="$public_root/brassica_rapa/genome/Brassica_rapa.Brapa_1.0.dna.toplevel.fa.gz"
    [bol]="$public_root/brassica_oleracea/genome/Brassica_oleracea.BOL.dna.toplevel.fa.gz"
)
declare -A gff=(
    [brapa]="$public_root/brassica_rapa/annotation/Brassica_rapa.Brapa_1.0.62.gff3.gz"
    [bol]="$public_root/brassica_oleracea/annotation/Brassica_oleracea.BOL.62.gff3.gz"
)
declare -A primary=(
    [brapa]="$code_root/config/primary_seqids/brassica_rapa_ensembl62.tsv"
    [bol]="$code_root/config/primary_seqids/brassica_oleracea_ensembl62.tsv"
)

for reference in brapa bol; do
    for required in "${genome[$reference]}" "${gff[$reference]}" "${primary[$reference]}"; do
        if [[ ! -s $required ]]; then
            echo "missing or empty Brassica reference input: $required" >&2
            exit 1
        fi
    done
done
if [[ -e $result_root || -e $working_root ]]; then
    echo "refusing to overwrite Brassica multireference input bundle" >&2
    exit 1
fi
mkdir -p "$working_root"

{
    printf 'field\tvalue\n'
    printf 'code_commit\t%s\n' "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'development_species\tBrassica_napus\n'
    printf 'reference_lineages\tBrassica_rapa_A,Brassica_oleracea_C\n'
    printf 'method_family_vote_policy\tmultiple_references_per_method_count_as_one_vote\n'
    printf 'target_truth_access\tfalse\n'
    printf 'normalization_module_sha256\t%s\n' \
        "$(sha256sum "$code_root/src/ploidypatch/normalize.py" | awk '{print $1}')"
} > "$working_root/run_contract.tsv"
{
    printf 'reference\trole\tbytes\tsha256\tpath\n'
    for reference in brapa bol; do
        for role in genome gff primary; do
            case $role in
                genome) path=${genome[$reference]} ;;
                gff) path=${gff[$reference]} ;;
                primary) path=${primary[$reference]} ;;
            esac
            printf '%s\t%s\t%s\t%s\t%s\n' "$reference" "$role" \
                "$(stat -Lc %s "$path")" \
                "$(sha256sum "$path" | awk '{print $1}')" "$path"
        done
    done
} > "$working_root/input_manifest.tsv"

cd "$code_root"
for reference in brapa bol; do
    /usr/bin/time -v -o "$working_root/$reference.resource.time.txt" \
        "$python_bin" -m ploidypatch.cli normalize primary-annotation \
            --gff "${gff[$reference]}" \
            --genome "${genome[$reference]}" \
            --primary-seqid-table "${primary[$reference]}" \
            --output-dir "$working_root/$reference" \
            > "$working_root/$reference.stdout.json" \
            2> "$working_root/$reference.stderr.log"
done

for reference in brapa bol; do
    for required in \
        "$working_root/$reference/primary_chromosomes.genome.fa" \
        "$working_root/$reference/primary_chromosomes.genome.fa.fai" \
        "$working_root/$reference/primary_chromosomes.gff3" \
        "$working_root/$reference/manifest.json"; do
        if [[ ! -s $required ]]; then
            echo "missing normalized Brassica reference artifact: $required" >&2
            exit 1
        fi
    done
done
if [[ $(wc -l < "$working_root/brapa/primary_chromosomes.genome.fa.fai") -ne 10 ]]; then
    echo "B. rapa primary chromosome count is not ten" >&2
    exit 1
fi
if [[ $(wc -l < "$working_root/bol/primary_chromosomes.genome.fa.fai") -ne 9 ]]; then
    echo "B. oleracea primary chromosome count is not nine" >&2
    exit 1
fi
(
    cd "$working_root"
    find . -type f \( -name '*.fa' -o -name '*.fai' -o -name '*.gff3' -o -name '*.json' -o -name '*.tsv' \) \
        -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS
    sha256sum -c SHA256SUMS >/dev/null
)
du -sb "$working_root" > "$working_root/disk_bytes.txt"
mv "$working_root" "$result_root"
printf 'Brassica multireference inputs frozen: %s\n' "$result_root"
