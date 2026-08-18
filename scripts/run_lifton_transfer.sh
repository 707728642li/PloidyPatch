#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 7 || $# -gt 8 ]]; then
    echo "usage: $0 ENV_PREFIX RUN_ROOT TARGET_FASTA REFERENCE_FASTA REFERENCE_GFF REFERENCE_ID THREADS [--resume-completed]" >&2
    exit 2
fi

env_prefix=$1
run_root=$2
target_fasta=$3
reference_fasta=$4
reference_gff=$5
reference_id=$6
threads=$7
resume_completed=false
if [[ $# -eq 8 ]]; then
    [[ $8 == --resume-completed ]] || { echo "unknown LiftOn option: $8" >&2; exit 2; }
    resume_completed=true
fi

if [[ ! $reference_id =~ ^[A-Za-z0-9][A-Za-z0-9_.-]*$ ]]; then
    echo "reference ID must be a safe identifier" >&2
    exit 2
fi
if [[ ! $threads =~ ^[1-9][0-9]*$ ]]; then
    echo "threads must be a positive integer" >&2
    exit 2
fi
for executable in lifton minimap2 miniprot; do
    if [[ ! -x $env_prefix/bin/$executable ]]; then
        echo "required executable is absent: $env_prefix/bin/$executable" >&2
        exit 2
    fi
done
for input_path in "$target_fasta" "$reference_fasta" "$reference_gff"; do
    if [[ ! -s $input_path ]]; then
        echo "missing or empty input: $input_path" >&2
        exit 2
    fi
done

env_prefix=$(realpath "$env_prefix")
run_root=$(realpath -m "$run_root")
target_fasta=$(realpath "$target_fasta")
reference_fasta=$(realpath "$reference_fasta")
reference_gff=$(realpath "$reference_gff")
working_root=${run_root}.working
for safe_path in "$env_prefix" "$run_root" "$target_fasta" \
                 "$reference_fasta" "$reference_gff"; do
    if [[ ! $safe_path =~ ^/[A-Za-z0-9._/-]+$ ]]; then
        echo "path contains characters unsafe for an upstream child command: $safe_path" >&2
        exit 2
    fi
done
if [[ -e $run_root ]]; then
    echo "refusing to overwrite a completed run: $run_root" >&2
    exit 2
fi
if [[ -e $working_root && $resume_completed == false ]]; then
    echo "refusing to overwrite a working run: $working_root" >&2
    exit 2
fi

if [[ $resume_completed == true ]]; then
    for artifact in "$working_root/run_contract.tsv" \
                    "$working_root/input_manifest.tsv" \
                    "$working_root/command.sh" "$working_root/stdout.log" \
                    "$working_root/stderr.log" "$working_root/resource.time.txt"; do
        [[ -s $artifact ]] || { echo "resume artifact missing: $artifact" >&2; exit 1; }
    done
else
    mkdir -p "$working_root/upstream"
fi
upstream=$working_root/upstream
output_gff=$upstream/lifton.gff3
unmapped=$upstream/unmapped_features.txt
if [[ $resume_completed == false ]]; then
{
    printf 'field\tvalue\n'
    printf 'environment\t%s\n' "$env_prefix"
    printf 'target_fasta\t%s\n' "$target_fasta"
    printf 'reference_fasta\t%s\n' "$reference_fasta"
    printf 'reference_gff\t%s\n' "$reference_gff"
    printf 'reference_id\t%s\n' "$reference_id"
    printf 'threads\t%s\n' "$threads"
    printf 'feature_scope\tgene_only\n'
    printf 'locus_pipeline\ttrue\n'
    printf 'validate_output\ttrue\n'
    printf 'reference_proteins\ttool_extracted_from_declared_reference\n'
    printf 'rescue_merge_alignment\tLiftOn_1.0.11_defaults\n'
} > "$working_root/run_contract.tsv"

lifton_command=(
    conda run -p "$env_prefix" --no-capture-output
    lifton "$target_fasta" "$reference_fasta"
    -g "$reference_gff" -o "$output_gff" -u "$unmapped"
    -t "$threads" --gene-only --locus-pipeline --validate-output -time
)
printf '%q ' "${lifton_command[@]}" > "$working_root/command.sh"
printf '\n' >> "$working_root/command.sh"

{
    printf 'role\tbytes\tsha256\tpath\n'
    for entry in \
        "target_fasta:$target_fasta" \
        "reference_fasta:$reference_fasta" \
        "reference_gff:$reference_gff"; do
        role=${entry%%:*}
        path=${entry#*:}
        digest=$(sha256sum "$path" | awk '{print $1}')
        printf '%s\t%s\t%s\t%s\n' \
            "$role" "$(stat -Lc %s "$path")" "$digest" "$path"
    done
} > "$working_root/input_manifest.tsv"

cd "$working_root"
/usr/bin/time -v -o resource.time.txt \
    "${lifton_command[@]}" > stdout.log 2> stderr.log
else
    cd "$working_root"
fi

if [[ ! -s $output_gff ]]; then
    echo "missing or empty LiftOn output GFF: $output_gff" >&2
    exit 1
fi
if grep -Eqi 'Segmentation fault|Killed|command not found' \
    stdout.log stderr.log; then
    echo "LiftOn log contains a fatal signature" >&2
    exit 1
fi
caught_gene_errors=$(awk '
    /^\[ERROR\] Error during (Liftoff|Miniprot) gene processing/ { count++ }
    END { print count + 0 }
' stdout.log stderr.log)
tracebacks=$(awk '
    index($0, "Traceback (most recent call last)") { count++ }
    END { print count + 0 }
' stdout.log stderr.log)
if [[ $tracebacks -ne $caught_gene_errors ]]; then
    echo "LiftOn log contains an unclassified traceback" >&2
    exit 1
fi
{
    printf 'warning\tcount\n'
    printf 'caught_per_gene_processing_error\t%s\n' "$caught_gene_errors"
    printf 'classified_traceback\t%s\n' "$tracebacks"
} > upstream_warning_counts.tsv

target_fai=${target_fasta}.fai
if [[ ! -s $target_fai ]]; then
    echo "target FASTA index is required for seqid validation: $target_fai" >&2
    exit 1
fi
if ! awk -F '\t' '
    NR == FNR { valid[$1] = 1; next }
    /^#/ || NF == 0 { next }
    NF != 9 || !($1 in valid) || $4 !~ /^[0-9]+$/ || $5 !~ /^[0-9]+$/ || $4 > $5 { exit 1 }
' "$target_fai" "$output_gff"; then
    echo "LiftOn output GFF failed target-coordinate validation" >&2
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
' "$output_gff" > feature_counts.tsv

{
    printf 'role\tbytes\tsha256\tpath\n'
    digest=$(sha256sum "$output_gff" | awk '{print $1}')
    printf 'lifton.gff3\t%s\t%s\t%s\n' \
        "$(stat -Lc %s "$output_gff")" "$digest" "$output_gff"
    if [[ -e $unmapped ]]; then
        digest=$(sha256sum "$unmapped" | awk '{print $1}')
        printf 'unmapped_features.txt\t%s\t%s\t%s\n' \
            "$(stat -Lc %s "$unmapped")" "$digest" "$unmapped"
    fi
} > output_manifest.tsv
du -sb "$working_root" > disk_bytes.txt

cd "$(dirname "$run_root")"
mv "$working_root" "$run_root"
printf 'LiftOn transfer run validated: %s\n' "$run_root"
