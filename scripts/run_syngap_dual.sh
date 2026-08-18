#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 10 ]]; then
    echo "usage: $0 ENV_PREFIX RUN_ROOT SP1_FASTA SP1_GFF SP2_FASTA SP2_GFF SP1 SP2 THREADS PROCESS" >&2
    exit 2
fi

env_prefix=$1
run_root=$2
sp1_fasta=$3
sp1_gff=$4
sp2_fasta=$5
sp2_gff=$6
sp1=$7
sp2=$8
threads=$9
process=${10}

if [[ ! $sp1 =~ ^[A-Za-z0-9][A-Za-z0-9_.-]*$ ]] ||
   [[ ! $sp2 =~ ^[A-Za-z0-9][A-Za-z0-9_.-]*$ ]]; then
    echo "species names must be shell-safe identifiers" >&2
    exit 2
fi
if [[ ! $threads =~ ^[1-9][0-9]*$ ]]; then
    echo "threads must be a positive integer" >&2
    exit 2
fi
if [[ $process != genblastg && $process != miniprot ]]; then
    echo "process must be genblastg or miniprot" >&2
    exit 2
fi
if [[ -e $run_root ]]; then
    echo "refusing to overwrite run root: $run_root" >&2
    exit 2
fi
if [[ ! -x $env_prefix/bin/syngap ]]; then
    echo "SynGAP executable is absent from environment: $env_prefix" >&2
    exit 2
fi

env_prefix=$(realpath "$env_prefix")
run_root=$(realpath -m "$run_root")
if [[ ! $env_prefix =~ ^/[A-Za-z0-9._/-]+$ ]] ||
   [[ ! $run_root =~ ^/[A-Za-z0-9._/-]+$ ]]; then
    echo "environment and run-root paths must be shell-safe for SynGAP 1.2.5" >&2
    exit 2
fi

for input_path in "$sp1_fasta" "$sp1_gff" "$sp2_fasta" "$sp2_gff"; do
    if [[ ! -s $input_path ]]; then
        echo "missing or empty input: $input_path" >&2
        exit 2
    fi
    resolved=$(realpath "$input_path")
    # SynGAP 1.2.5 constructs multiple unquoted os.system commands internally.
    # Reject paths that could be split or interpreted by its child shell.
    if [[ ! $resolved =~ ^/[A-Za-z0-9._/-]+$ ]]; then
        echo "input path is unsafe for SynGAP 1.2.5: $resolved" >&2
        exit 2
    fi
done
sp1_fasta=$(realpath "$sp1_fasta")
sp1_gff=$(realpath "$sp1_gff")
sp2_fasta=$(realpath "$sp2_fasta")
sp2_gff=$(realpath "$sp2_gff")

mkdir -p "$run_root"
cd "$run_root"
{
    printf 'field\tvalue\n'
    printf 'environment\t%s\n' "$(realpath "$env_prefix")"
    printf 'sp1\t%s\n' "$sp1"
    printf 'sp1_fasta\t%s\n' "$(realpath "$sp1_fasta")"
    printf 'sp1_gff\t%s\n' "$(realpath "$sp1_gff")"
    printf 'sp2\t%s\n' "$sp2"
    printf 'sp2_fasta\t%s\n' "$(realpath "$sp2_fasta")"
    printf 'sp2_gff\t%s\n' "$(realpath "$sp2_gff")"
    printf 'threads\t%s\n' "$threads"
    printf 'process\t%s\n' "$process"
    printf 'datatype\tnucl\n'
    printf 'cscore\t0.7\n'
    printf 'evalue\t1e-5\n'
    printf 'rank\t5\n'
    printf 'coverage\t0.5\n'
} > run_contract.tsv

printf '%q ' \
    conda run -p "$env_prefix" --no-capture-output \
    syngap dual \
    --sp1fa "$sp1_fasta" --sp1gff "$sp1_gff" \
    --sp2fa "$sp2_fasta" --sp2gff "$sp2_gff" \
    --sp1 "$sp1" --sp2 "$sp2" --threads "$threads" --process "$process" \
    > command.sh
printf '\n' >> command.sh

/usr/bin/time -v -o resource.time.txt \
    conda run -p "$env_prefix" --no-capture-output \
    syngap dual \
    --sp1fa "$sp1_fasta" --sp1gff "$sp1_gff" \
    --sp2fa "$sp2_fasta" --sp2gff "$sp2_gff" \
    --sp1 "$sp1" --sp2 "$sp2" --threads "$threads" --process "$process" \
    > stdout.log 2> stderr.log

shopt -s nullglob
outer_dirs=(SynGAP_dual_*)
if [[ ${#outer_dirs[@]} -ne 1 || ! -d ${outer_dirs[0]} ]]; then
    echo "expected exactly one SynGAP dual directory" >&2
    exit 1
fi
inner_dirs=("${outer_dirs[0]}"/SynGAP_"${sp1}"_"${sp2}"_*)
if [[ ${#inner_dirs[@]} -ne 1 || ! -d ${inner_dirs[0]} ]]; then
    echo "expected exactly one nested SynGAP analysis directory" >&2
    exit 1
fi
results_dir=${inner_dirs[0]}/results
full_gff=$results_dir/${sp1}.SynGAP.gff3
clean_gff=$results_dir/${sp1}.SynGAP.clean.gff3
miss_gff=$results_dir/${sp1}.SynGAP.clean.miss_annotated.gff3
mis_gff=$results_dir/${sp1}.SynGAP.clean.mis_annotated.gff3
anchors=$results_dir/${sp1}.${sp2}.anchors
gap_anchors=$results_dir/${sp1}.${sp2}.anchors.gap

if [[ ! -s $full_gff ]]; then
    echo "missing or empty full polished annotation: $full_gff" >&2
    exit 1
fi
for sentinel in "$clean_gff" "$miss_gff" "$mis_gff" "$anchors" "$gap_anchors"; do
    if [[ ! -e $sentinel ]]; then
        echo "missing SynGAP sentinel: $sentinel" >&2
        exit 1
    fi
done
if ! grep -Fq "SynGAP analysis for" stdout.log ||
   ! grep -Fq "Please check the result files" stdout.log; then
    echo "SynGAP completion messages are absent" >&2
    exit 1
fi
if grep -Eqi 'Traceback|command not found|No such file or directory|Segmentation fault|Killed' \
    stdout.log stderr.log; then
    echo "SynGAP logs contain a fatal child-command signature" >&2
    exit 1
fi

source_bytes=$(stat -c %s "$sp1_gff")
source_sha=$(sha256sum "$sp1_gff" | awk '{print $1}')
prefix_sha=$(head -c "$source_bytes" "$full_gff" | sha256sum | awk '{print $1}')
if [[ $source_sha != "$prefix_sha" ]]; then
    echo "full SynGAP output does not preserve the target input GFF prefix" >&2
    exit 1
fi

{
    printf 'role\tbytes\tsha256\tpath\n'
    for artifact in "$full_gff" "$clean_gff" "$miss_gff" "$mis_gff" \
                    "$anchors" "$gap_anchors"; do
        artifact_sha=$(sha256sum "$artifact" | awk '{print $1}')
        printf '%s\t%s\t%s\t%s\n' \
            "$(basename "$artifact")" "$(stat -Lc %s "$artifact")" \
            "$artifact_sha" "$(realpath "$artifact")"
    done
} > output_manifest.tsv

printf 'SynGAP dual run validated: %s\n' "$(realpath "$results_dir")"
