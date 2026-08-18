#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 3 ]]; then
    echo "usage: $0 RUN_ROOT TARGET_FASTA PUBLICATION_MODE" >&2
    exit 2
fi

run_root=$(realpath -m "$1")
target_fasta=$(realpath "$2")
publication_mode=$3
working_root=${run_root}.working
if [[ ! $publication_mode =~ ^[A-Za-z0-9][A-Za-z0-9_.-]*$ ]]; then
    echo "publication mode must be a safe identifier" >&2
    exit 2
fi
if [[ -e $run_root || ! -d $working_root ]]; then
    echo "expected one unpublished GeMoMa working directory: $working_root" >&2
    exit 2
fi
if [[ ! -s $working_root/resource.time.txt ]] || ! awk '
    $1 == "Exit" && $2 == "status:" && $3 == "0" { success = 1 }
    END { exit !success }
' "$working_root/resource.time.txt"; then
    echo "GeMoMa resource record does not contain exit status 0" >&2
    exit 1
fi

upstream=$working_root/upstream
final_gff=$upstream/final_annotation.gff
unfiltered_gff=$upstream/unfiltered_predictions_from_species_0.gff
reference_table=$upstream/reference_gene_table.tabular
predicted_proteins=$upstream/predicted_proteins.fasta
protocol=$upstream/protocol_GeMoMaPipeline.txt
sentinels=(
    "$final_gff"
    "$unfiltered_gff"
    "$reference_table"
    "$predicted_proteins"
    "$protocol"
)
for sentinel in "${sentinels[@]}"; do
    if [[ ! -s $sentinel ]]; then
        echo "missing or empty GeMoMa sentinel: $sentinel" >&2
        exit 1
    fi
done
if grep -Eqi \
    'Exception in thread|OutOfMemoryError|Segmentation fault|(^|[[:space:]])Killed([[:space:]]|$)|command not found' \
    "$working_root/stdout.log" "$working_root/stderr.log"; then
    echo "GeMoMa log contains a fatal signature" >&2
    exit 1
fi
if ! grep -Fq 'No errors detected.' "$working_root/stdout.log"; then
    echo "GeMoMa success signature is absent" >&2
    exit 1
fi

target_fai=${target_fasta}.fai
if [[ ! -s $target_fai ]]; then
    echo "target FASTA index is required for seqid validation: $target_fai" >&2
    exit 1
fi
if ! awk -F '\t' '
    NR == FNR { valid[$1] = 1; next }
    /^#/ || NF == 0 { next }
    NF != 9 || !($1 in valid) || $4 !~ /^[0-9]+$/ ||
        $5 !~ /^[0-9]+$/ || $4 > $5 { exit 1 }
' "$target_fai" "$final_gff"; then
    echo "GeMoMa final GFF failed target-coordinate validation" >&2
    exit 1
fi

awk -F '\t' '
    BEGIN { OFS="\t"; print "feature", "count" }
    /^#/ || NF != 9 { next }
    $3 == "gene" { gene++ }
    $3 == "mRNA" || $3 == "transcript" { transcript++ }
    $3 == "CDS" { cds++ }
    END {
        print "gene", gene + 0
        print "transcript", transcript + 0
        print "CDS", cds + 0
        if (gene == 0 || transcript == 0 || cds == 0) exit 1
    }
' "$final_gff" > "$working_root/feature_counts.tsv"

{
    printf 'role\tbytes\tsha256\tpath\n'
    for artifact in "${sentinels[@]}"; do
        digest=$(sha256sum "$artifact" | awk '{print $1}')
        printf '%s\t%s\t%s\t%s\n' \
            "$(basename "$artifact")" "$(stat -Lc %s "$artifact")" \
            "$digest" "$artifact"
    done
} > "$working_root/output_manifest.tsv"
{
    printf 'field\tvalue\n'
    printf 'publication_mode\t%s\n' "$publication_mode"
    printf 'publisher_script_sha256\t%s\n' \
        "$(sha256sum "$0" | awk '{print $1}')"
    printf 'published_at\t%s\n' "$(date --iso-8601=seconds)"
} > "$working_root/publication.tsv"
du -sb "$working_root" > "$working_root/disk_bytes.txt"

mv "$working_root" "$run_root"
printf 'GeMoMa homology run validated: %s\n' "$run_root"
