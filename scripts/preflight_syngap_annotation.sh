#!/usr/bin/env bash
set -euo pipefail
export LC_ALL=C

if [[ $# -ne 5 ]]; then
    echo "usage: $0 ENV_PREFIX OUTPUT_DIR LABEL GENOME_FASTA ANNOTATION_GFF" >&2
    exit 2
fi

env_prefix=$(realpath "$1")
output_dir=$(realpath -m "$2")
label=$3
genome_fasta=$(realpath "$4")
annotation_gff=$(realpath "$5")

if [[ ! $label =~ ^[A-Za-z0-9][A-Za-z0-9_.-]*$ ]]; then
    echo "label must be a shell-safe identifier" >&2
    exit 2
fi
for input_path in "$env_prefix" "$genome_fasta" "$annotation_gff"; do
    if [[ ! $input_path =~ ^/[A-Za-z0-9._/-]+$ ]]; then
        echo "path is unsafe for SynGAP 1.2.5: $input_path" >&2
        exit 2
    fi
done
if [[ ! -x $env_prefix/bin/syngap ]]; then
    echo "SynGAP executable is absent: $env_prefix/bin/syngap" >&2
    exit 2
fi
if [[ ! -s $genome_fasta || ! -s $annotation_gff ]]; then
    echo "genome or annotation input is absent/empty" >&2
    exit 2
fi
if [[ -e $output_dir ]]; then
    echo "refusing to overwrite preflight output: $output_dir" >&2
    exit 2
fi
mkdir -p "$output_dir"

bed=$output_dir/$label.primary.bed
primary_ids=$output_dir/$label.primary.ids
all_cds=$output_dir/$label.all.cds.fa
all_pep=$output_dir/$label.all.pep.fa
primary_cds=$output_dir/$label.primary.cds.fa
primary_pep=$output_dir/$label.primary.pep.fa

/usr/bin/time -v -o "$output_dir/jcvi_bed.time.txt" \
    conda run -p "$env_prefix" --no-capture-output \
    python -m jcvi.formats.gff bed \
    --type=mRNA --key=ID --parent_key=Parent --primary_only \
    "$annotation_gff" -o "$bed" \
    > "$output_dir/jcvi_bed.stdout.log" \
    2> "$output_dir/jcvi_bed.stderr.log"

cut -f4 "$bed" > "$primary_ids"
sort "$primary_ids" | uniq -d > "$output_dir/duplicate_primary.ids"
if [[ -s $output_dir/duplicate_primary.ids ]]; then
    echo "JCVI primary IDs are not unique" >&2
    exit 1
fi

/usr/bin/time -v -o "$output_dir/gffread_cds.time.txt" \
    conda run -p "$env_prefix" --no-capture-output \
    gffread "$annotation_gff" -g "$genome_fasta" -x "$all_cds" \
    > "$output_dir/gffread_cds.stdout.log" \
    2> "$output_dir/gffread_cds.stderr.log"
/usr/bin/time -v -o "$output_dir/gffread_pep.time.txt" \
    conda run -p "$env_prefix" --no-capture-output \
    gffread "$annotation_gff" -g "$genome_fasta" -y "$all_pep" -S \
    > "$output_dir/gffread_pep.stdout.log" \
    2> "$output_dir/gffread_pep.stderr.log"

conda run -p "$env_prefix" --no-capture-output \
    seqkit grep -f "$primary_ids" "$all_cds" -o "$primary_cds" \
    > "$output_dir/seqkit_cds.stdout.log" \
    2> "$output_dir/seqkit_cds.stderr.log"
conda run -p "$env_prefix" --no-capture-output \
    seqkit grep -f "$primary_ids" "$all_pep" -o "$primary_pep" \
    > "$output_dir/seqkit_pep.stdout.log" \
    2> "$output_dir/seqkit_pep.stderr.log"

conda run -p "$env_prefix" --no-capture-output \
    seqkit seq -n -i "$primary_cds" \
    | sort -u > "$output_dir/primary_cds.ids"
conda run -p "$env_prefix" --no-capture-output \
    seqkit seq -n -i "$primary_pep" \
    | sort -u > "$output_dir/primary_pep.ids"
sort -u "$primary_ids" > "$output_dir/primary.ids.sorted"
comm -23 "$output_dir/primary.ids.sorted" "$output_dir/primary_cds.ids" \
    > "$output_dir/missing_primary_cds.ids"
comm -23 "$output_dir/primary.ids.sorted" "$output_dir/primary_pep.ids" \
    > "$output_dir/missing_primary_pep.ids"
if [[ -s $output_dir/missing_primary_cds.ids ||
      -s $output_dir/missing_primary_pep.ids ]]; then
    echo "one or more JCVI primary transcripts lack extracted CDS/protein" >&2
    exit 1
fi

primary_count=$(wc -l < "$primary_ids")
cds_count=$(grep -c '^>' "$primary_cds")
pep_count=$(grep -c '^>' "$primary_pep")
if [[ $primary_count -ne $cds_count || $primary_count -ne $pep_count ]]; then
    echo "primary ID/CDS/protein record counts disagree" >&2
    exit 1
fi

{
    printf 'field\tvalue\n'
    printf 'label\t%s\n' "$label"
    printf 'genome_path\t%s\n' "$genome_fasta"
    printf 'genome_sha256\t%s\n' "$(sha256sum "$genome_fasta" | awk '{print $1}')"
    printf 'annotation_path\t%s\n' "$annotation_gff"
    printf 'annotation_sha256\t%s\n' "$(sha256sum "$annotation_gff" | awk '{print $1}')"
    printf 'jcvi_primary_ids\t%s\n' "$primary_count"
    printf 'primary_cds_records\t%s\n' "$cds_count"
    printf 'primary_pep_records\t%s\n' "$pep_count"
    printf 'missing_primary_cds\t0\n'
    printf 'missing_primary_pep\t0\n'
} > "$output_dir/preflight_summary.tsv"

{
    printf 'artifact\tbytes\tsha256\n'
    for artifact in "$bed" "$primary_ids" "$all_cds" "$all_pep" \
                    "$primary_cds" "$primary_pep"; do
        artifact_sha=$(sha256sum "$artifact" | awk '{print $1}')
        printf '%s\t%s\t%s\n' \
            "$(basename "$artifact")" "$(stat -c %s "$artifact")" "$artifact_sha"
    done
} > "$output_dir/output_manifest.tsv"

printf 'SynGAP annotation preflight passed: %s (%s primary transcripts)\n' \
    "$label" "$primary_count"
