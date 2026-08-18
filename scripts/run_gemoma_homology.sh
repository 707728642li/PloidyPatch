#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 7 ]]; then
    echo "usage: $0 ENV_PREFIX RUN_ROOT TARGET_FASTA REFERENCE_FASTA REFERENCE_GFF REFERENCE_ID THREADS" >&2
    exit 2
fi

env_prefix=$1
run_root=$2
target_fasta=$3
reference_fasta=$4
reference_gff=$5
reference_id=$6
threads=$7

if [[ ! $reference_id =~ ^[A-Za-z0-9][A-Za-z0-9_.-]*$ ]]; then
    echo "reference ID must be a safe identifier" >&2
    exit 2
fi
if [[ ! $threads =~ ^[1-9][0-9]*$ ]]; then
    echo "threads must be a positive integer" >&2
    exit 2
fi
if [[ ! -x $env_prefix/bin/GeMoMa ]]; then
    echo "GeMoMa executable is absent: $env_prefix/bin/GeMoMa" >&2
    exit 2
fi
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
if [[ -e $run_root || -e $working_root ]]; then
    echo "refusing to overwrite a completed or working run: $run_root" >&2
    exit 2
fi

mkdir -p "$working_root/upstream"
upstream=$working_root/upstream
{
    printf 'field\tvalue\n'
    printf 'environment\t%s\n' "$env_prefix"
    printf 'target_fasta\t%s\n' "$target_fasta"
    printf 'reference_fasta\t%s\n' "$reference_fasta"
    printf 'reference_gff\t%s\n' "$reference_gff"
    printf 'reference_id\t%s\n' "$reference_id"
    printf 'threads\t%s\n' "$threads"
    printf 'jvm_initial_heap\t8g\n'
    printf 'jvm_max_heap\t128g\n'
    printf 'rna_evidence\tNO\n'
    printf 'search\tmmseqs\n'
    printf 'GeMoMa.Score\tReAlign\n'
    printf 'AnnotationFinalizer.r\tNO\n'
    printf 'output_individual_predictions\ttrue\n'
} > "$working_root/run_contract.tsv"

gemoma_command=(
    conda run -p "$env_prefix" --no-capture-output
    GeMoMa -Xms8g -Xmx128g GeMoMaPipeline
    "threads=$threads" r=NO tblastn=false GeMoMa.Score=ReAlign
    AnnotationFinalizer.r=NO o=true "t=$target_fasta"
    s=own "i=$reference_id" "a=$reference_gff" "g=$reference_fasta"
    "outdir=$upstream"
)
printf '%q ' "${gemoma_command[@]}" > "$working_root/command.sh"
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

/usr/bin/time -v -o "$working_root/resource.time.txt" \
    "${gemoma_command[@]}" \
    > "$working_root/stdout.log" 2> "$working_root/stderr.log"

bash "$(dirname "$0")/publish_gemoma_working.sh" \
    "$run_root" "$target_fasta" fresh_run
