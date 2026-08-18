#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 3 ]]; then
    echo "usage: $0 HAWTHORN_ROOT OUTPUT_DIR SAMTOOLS_ENV_PREFIX" >&2
    exit 2
fi

source_root=$(realpath "$1")
output_dir=$(realpath -m "$2")
env_prefix=$(realpath "$3")
rna_root=$source_root/RNASeq
bam_dir=$rna_root/bam_folder
raw_dir=$rna_root/raw_data
clean_dir=$rna_root/clean_fastq
rna_reference=$rna_root/genome/Black_Primary.genome.fa
rna_gff=$rna_root/genome/Black_Primary.gff3
catalog_reference=$source_root/Black_Primary.fa
catalog_gff=$source_root/Black_Primary.gff3
catalog_fai=$source_root/Black_Primary.fa.fai
samtools=$env_prefix/bin/samtools

for required in "$rna_root" "$bam_dir" "$raw_dir" "$clean_dir" \
                "$rna_reference" "$rna_gff" "$catalog_reference" \
                "$catalog_gff" "$catalog_fai" "$samtools"; do
    if [[ ! -e $required ]]; then
        echo "required hawthorn RNA audit input is absent: $required" >&2
        exit 2
    fi
done
if [[ -e $output_dir ]]; then
    echo "refusing to overwrite RNA audit: $output_dir" >&2
    exit 2
fi
mkdir -p "$output_dir/headers"

"$samtools" --version > "$output_dir/samtools_version.txt"
{
    printf 'role\tbytes\tsha256\tpath\n'
    for entry in \
        "rna_reference:$rna_reference" \
        "catalog_reference:$catalog_reference" \
        "rna_gff:$rna_gff" \
        "catalog_gff:$catalog_gff"; do
        role=${entry%%:*}
        path=${entry#*:}
        digest=$(sha256sum "$path" | awk '{print $1}')
        printf '%s\t%s\t%s\t%s\n' \
            "$role" "$(stat -Lc %s "$path")" "$digest" "$path"
    done
} > "$output_dir/reference_identity.tsv"

awk 'BEGIN { OFS="\t" } { print $1, $2 }' "$catalog_fai" \
    > "$output_dir/catalog_reference.sq.tsv"
catalog_sq_sha=$(sha256sum "$output_dir/catalog_reference.sq.tsv" | awk '{print $1}')

shopt -s nullglob
bams=("$bam_dir"/*.sorted.bam)
if [[ ${#bams[@]} -eq 0 ]]; then
    echo "no sorted BAM files found" >&2
    exit 1
fi
printf 'sample\tbam_bytes\tbai_present\tquickcheck\tcoordinate_sorted\t' \
    > "$output_dir/bam_audit.tsv"
printf 'sq_count\tsq_sha256\treference_match\tmapped_alignments\t' \
    >> "$output_dir/bam_audit.tsv"
printf 'unmapped_alignments\tread_groups\tprograms\tbam_path\n' \
    >> "$output_dir/bam_audit.tsv"

for bam in "${bams[@]}"; do
    sample=$(basename "$bam" .sorted.bam)
    header=$output_dir/headers/$sample.header.sam
    sq=$output_dir/headers/$sample.sq.tsv
    "$samtools" view -H "$bam" > "$header"
    awk -F '\t' '
        BEGIN { OFS="\t" }
        $1 == "@SQ" {
            sn=""; ln=""
            for (i=2; i<=NF; i++) {
                if ($i ~ /^SN:/) sn=substr($i, 4)
                if ($i ~ /^LN:/) ln=substr($i, 4)
            }
            if (sn == "" || ln == "") exit 1
            print sn, ln
        }
    ' "$header" > "$sq"
    sq_sha=$(sha256sum "$sq" | awk '{print $1}')
    if [[ $sq_sha == "$catalog_sq_sha" ]]; then
        reference_match=true
    else
        reference_match=false
    fi
    if "$samtools" quickcheck "$bam"; then
        quickcheck=pass
    else
        quickcheck=fail
    fi
    if grep -Fq 'SO:coordinate' "$header"; then
        coordinate_sorted=true
    else
        coordinate_sorted=false
    fi
    if [[ -s ${bam}.bai ]]; then
        bai_present=true
    else
        bai_present=false
    fi
    read mapped unmapped < <(
        "$samtools" idxstats "$bam" | awk \
            '{ mapped += $3; unmapped += $4 } END { print mapped+0, unmapped+0 }'
    )
    read_groups=$(awk -F '\t' '$1 == "@RG" { n++ } END { print n+0 }' "$header")
    programs=$(awk -F '\t' '$1 == "@PG" { n++ } END { print n+0 }' "$header")
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$sample" "$(stat -Lc %s "$bam")" "$bai_present" "$quickcheck" \
        "$coordinate_sorted" "$(wc -l < "$sq")" "$sq_sha" \
        "$reference_match" "$mapped" "$unmapped" "$read_groups" \
        "$programs" "$bam" >> "$output_dir/bam_audit.tsv"
done

printf 'sample\traw_read1_bytes\traw_read2_bytes\tclean_read1_bytes\tclean_read2_bytes\n' \
    > "$output_dir/fastq_pairs.tsv"
raw_read1=("$raw_dir"/*_1.fq.gz)
for read1 in "${raw_read1[@]}"; do
    sample=$(basename "$read1" _1.fq.gz)
    read2=$raw_dir/${sample}_2.fq.gz
    clean1=$clean_dir/${sample}_1.clean.fq.gz
    clean2=$clean_dir/${sample}_2.clean.fq.gz
    for mate in "$read2" "$clean1" "$clean2"; do
        if [[ ! -s $mate ]]; then
            echo "missing FASTQ mate or clean file for $sample: $mate" >&2
            exit 1
        fi
    done
    printf '%s\t%s\t%s\t%s\t%s\n' "$sample" \
        "$(stat -Lc %s "$read1")" "$(stat -Lc %s "$read2")" \
        "$(stat -Lc %s "$clean1")" "$(stat -Lc %s "$clean2")" \
        >> "$output_dir/fastq_pairs.tsv"
done

find "$source_root" -type f \
    \( -iname '*isoseq*' -o -iname '*iso-seq*' -o -iname '*flnc*' \
       -o -iname '*pacbio*' -o -iname '*nanopore*' -o -iname '*hifi*' \
       -o -iname '*ccs*' \) -print > "$output_dir/long_read_name_hits.txt"

{
    printf 'metric\tvalue\n'
    printf 'bam_files\t%s\n' "${#bams[@]}"
    printf 'bam_indexes\t%s\n' "$(find "$bam_dir" -maxdepth 1 -type f -name '*.sorted.bam.bai' | wc -l)"
    printf 'raw_fastq_pairs\t%s\n' "${#raw_read1[@]}"
    printf 'bam_quickcheck_pass\t%s\n' "$(awk -F '\t' 'NR>1 && $4=="pass" {n++} END {print n+0}' "$output_dir/bam_audit.tsv")"
    printf 'bam_reference_match\t%s\n' "$(awk -F '\t' 'NR>1 && $8=="true" {n++} END {print n+0}' "$output_dir/bam_audit.tsv")"
    printf 'coordinate_sorted_bams\t%s\n' "$(awk -F '\t' 'NR>1 && $5=="true" {n++} END {print n+0}' "$output_dir/bam_audit.tsv")"
    printf 'bam_with_read_groups\t%s\n' "$(awk -F '\t' 'NR>1 && $11>0 {n++} END {print n+0}' "$output_dir/bam_audit.tsv")"
    printf 'long_read_name_hits\t%s\n' "$(wc -l < "$output_dir/long_read_name_hits.txt")"
} > "$output_dir/summary.tsv"

{
    printf 'artifact\tbytes\tsha256\n'
    for artifact in "$output_dir"/*.tsv "$output_dir"/*.txt; do
        [[ -f $artifact ]] || continue
        [[ $(basename "$artifact") != audit_artifacts.tsv ]] || continue
        printf '%s\t%s\t%s\n' "$(basename "$artifact")" \
            "$(stat -Lc %s "$artifact")" \
            "$(sha256sum "$artifact" | awk '{print $1}')"
    done
} > "$output_dir/audit_artifacts.tsv"

printf 'Hawthorn RNA audit passed: %s\n' "$output_dir"
