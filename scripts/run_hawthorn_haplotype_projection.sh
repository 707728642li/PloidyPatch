#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "usage: $0 PROJECT_ROOT HAWTHORN_SOURCE_ROOT" >&2
    exit 2
fi

project_root=$(realpath "$1")
source_root=$(realpath "$2")
code_root=$project_root/code
python_bin=$project_root/envs/ploidypatch-dev/bin/python
miniprot_bin=$project_root/envs/ploidypatch-baseline/bin/miniprot
data_root=$project_root/data/derived/hawthorn_black_projection_v0.1
data_working=${data_root}.working
result_root=$project_root/results/hawthorn/haplotype_projection_v0.1
result_working=${result_root}.working

inputs=(
    Black_Primary.fa
    Black_Primary.fa.fai
    Black_Primary.gff3
    Black_Hap1.gff3
    Black_Hap1.protein.fa
    Black_Hap2.gff3
    Black_Hap2.protein.fa
)
for name in "${inputs[@]}"; do
    if [[ ! -s $source_root/$name ]]; then
        echo "missing or empty hawthorn input: $source_root/$name" >&2
        exit 1
    fi
done
for required in "$python_bin" "$miniprot_bin" "$code_root/pyproject.toml"; do
    if [[ ! -e $required ]]; then
        echo "missing project executable or code input: $required" >&2
        exit 1
    fi
done
if [[ -e $data_root || -e $data_working || \
      -e $result_root || -e $result_working ]]; then
    echo "refusing to overwrite hawthorn projection data or result" >&2
    exit 1
fi

mkdir -p "$data_working" "$result_working/reference" \
    "$result_working/index" "$result_working/raw"
for name in "${inputs[@]}"; do
    cp --reflink=auto "$source_root/$name" "$data_working/$name"
done
{
    printf 'file_name\tbytes\tsha256\tsource_path\n'
    for name in "${inputs[@]}"; do
        source_hash=$(sha256sum "$source_root/$name" | awk '{print $1}')
        copied_hash=$(sha256sum "$data_working/$name" | awk '{print $1}')
        if [[ $source_hash != "$copied_hash" ]]; then
            echo "copied hawthorn input checksum mismatch: $name" >&2
            exit 1
        fi
        printf '%s\t%s\t%s\t%s\n' "$name" \
            "$(stat -Lc %s "$data_working/$name")" "$copied_hash" \
            "$source_root/$name"
    done
} > "$data_working/input_manifest.tsv"
mv "$data_working" "$data_root"

cd "$code_root"
"$python_bin" -m ploidypatch.cli baseline prepare-proteins \
    --protein "hap1=$data_root/Black_Hap1.protein.fa" \
    --protein "hap2=$data_root/Black_Hap2.protein.fa" \
    --output-fasta "$result_working/reference/haplotypes.protein.fa" \
    --output-map "$result_working/reference/haplotypes.map.tsv" \
    > "$result_working/reference/prepare.stdout.log" \
    2> "$result_working/reference/prepare.stderr.log"

/usr/bin/time -v -o "$result_working/index/resource.time.txt" \
    "$miniprot_bin" -t64 -d "$result_working/index/Black_Primary.mpi" \
        "$data_root/Black_Primary.fa" \
        > "$result_working/index/stdout.log" \
        2> "$result_working/index/stderr.log"
/usr/bin/time -v -o "$result_working/raw/resource.time.txt" \
    "$miniprot_bin" -I -t64 --gff-only -N4 --outn 4 --outs 0.8 --outc 0.5 \
        "$result_working/index/Black_Primary.mpi" \
        "$result_working/reference/haplotypes.protein.fa" \
        > "$result_working/raw/miniprot.gff3" \
        2> "$result_working/raw/stderr.log"

if [[ ! -s $result_working/index/Black_Primary.mpi || \
      ! -s $result_working/raw/miniprot.gff3 ]]; then
    echo "miniprot projection output is absent or empty" >&2
    exit 1
fi
if ! grep -q '^##gff-version' "$result_working/raw/miniprot.gff3"; then
    echo "miniprot output lacks the GFF version sentinel" >&2
    exit 1
fi
{
    printf 'field\tvalue\n'
    printf 'miniprot_version\t%s\n' \
        "$($miniprot_bin --version 2>&1 | head -1)"
    printf 'miniprot_sha256\t%s\n' \
        "$(sha256sum "$miniprot_bin" | awk '{print $1}')"
    printf 'threads\t64\n'
    printf 'parameters\t-I -N4 --outn 4 --outs 0.8 --outc 0.5\n'
    printf 'discovery_rna_used\tfalse\n'
} > "$result_working/run_contract.tsv"
{
    printf 'role\tbytes\tsha256\tpath\n'
    for artifact in \
        "$result_working/reference/haplotypes.protein.fa" \
        "$result_working/reference/haplotypes.map.tsv" \
        "$result_working/index/Black_Primary.mpi" \
        "$result_working/raw/miniprot.gff3"; do
        printf '%s\t%s\t%s\t%s\n' "$(basename "$artifact")" \
            "$(stat -Lc %s "$artifact")" \
            "$(sha256sum "$artifact" | awk '{print $1}')" "$artifact"
    done
} > "$result_working/output_manifest.tsv"
du -sb "$result_working" > "$result_working/disk_bytes.txt"
mv "$result_working" "$result_root"
printf 'hawthorn haplotype projection validated: %s\n' "$result_root"
