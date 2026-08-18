#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "usage: $0 PROJECT_ROOT HAWTHORN_SOURCE_ROOT" >&2
    exit 2
fi

project_root=$(realpath "$1")
hawthorn_root=$(realpath "$2")
code_root="$project_root/code"
python_bin="$project_root/envs/ploidypatch-dev/bin/python"
samtools_bin="$project_root/envs/ploidypatch-pav/bin/samtools"
bam_root="$hawthorn_root/RNASeq/bam_folder"
result_root="$project_root/results/hawthorn/rna_junctions_v0.1"
working_root="${result_root}.working"
log_root="$project_root/logs/hawthorn/rna_junctions_v0.1"

for required in "$code_root/pyproject.toml" "$python_bin" "$samtools_bin"; do
    if [[ ! -e "$required" ]]; then
        echo "required project artifact is absent: $required" >&2
        exit 1
    fi
done
if [[ -e "$result_root" || -e "$working_root" ]]; then
    echo "refusing to overwrite existing result or working directory" >&2
    exit 1
fi

if git -C "$code_root" rev-parse HEAD >/dev/null 2>&1; then
    code_commit=$(git -C "$code_root" rev-parse HEAD)
elif [[ -n "${PLOIDYPATCH_CODE_COMMIT:-}" ]]; then
    code_commit="$PLOIDYPATCH_CODE_COMMIT"
elif [[ -s "$code_root/.snapshot_commit" ]]; then
    code_commit=$(<"$code_root/.snapshot_commit")
else
    echo "code snapshot identity is unavailable" >&2
    exit 1
fi

mkdir -p "$working_root/samples" "$log_root"
{
    printf 'field\tvalue\n'
    printf 'project_root\t%s\n' "$project_root"
    printf 'hawthorn_source_root\t%s\n' "$hawthorn_root"
    printf 'code_commit\t%s\n' "$code_commit"
    if [[ -n "${PLOIDYPATCH_CODE_ARCHIVE_SHA256:-}" ]]; then
        printf 'code_archive_sha256\t%s\n' \
            "$PLOIDYPATCH_CODE_ARCHIVE_SHA256"
    elif [[ -s "$code_root/.snapshot_archive_sha256" ]]; then
        printf 'code_archive_sha256\t%s\n' \
            "$(<"$code_root/.snapshot_archive_sha256")"
    fi
    printf 'python\t%s\n' "$($python_bin --version 2>&1)"
    printf 'samtools_path\t%s\n' "$samtools_bin"
    printf 'samtools_sha256\t%s\n' "$(sha256sum "$samtools_bin" | cut -d ' ' -f 1)"
    printf 'min_mapq\t20\n'
    printf 'excluded_flag_mask\t2308\n'
    printf 'strand_policy\tunstranded\n'
} > "$working_root/environment.tsv"

pids=()
samples=(Black-1 Black-2 Black-3)
for sample in "${samples[@]}"; do
    bam="$bam_root/${sample}.sorted.bam"
    if [[ ! -f "$bam" || ! -f "${bam}.bai" ]]; then
        echo "BAM or adjacent BAI is absent: $bam" >&2
        exit 1
    fi
    (
        cd "$code_root"
        "$python_bin" -m ploidypatch.cli evidence extract-bam-junctions \
            --bam "$bam" \
            --sample "$sample" \
            --samtools "$samtools_bin" \
            --output-tsv "$working_root/samples/${sample}.junctions.tsv" \
            --threads 4 \
            --min-mapq 20
    ) > "$log_root/${sample}.stdout.log" \
      2> "$log_root/${sample}.stderr.log" &
    pids+=("$!")
done

status=0
for pid in "${pids[@]}"; do
    if ! wait "$pid"; then
        status=1
    fi
done
if [[ $status -ne 0 ]]; then
    echo "one or more primary junction extractions failed; preserving $working_root" >&2
    exit 1
fi

(
    cd "$code_root"
    "$python_bin" -m ploidypatch.cli evidence aggregate-junctions \
        --input-dir "$working_root/samples" \
        --primary-sample Black-1 \
        --primary-sample Black-2 \
        --primary-sample Black-3 \
        --output-tsv "$working_root/primary_aggregate.tsv" \
        --min-reads-per-sample 2 \
        --min-supporting-samples 2
) > "$log_root/aggregate.stdout.log" \
  2> "$log_root/aggregate.stderr.log"

sha256sum "$working_root"/samples/*.junctions.tsv \
    "$working_root"/samples/*.manifest.json \
    "$working_root"/primary_aggregate.tsv \
    "$working_root"/primary_aggregate.tsv.manifest.json \
    > "$working_root/SHA256SUMS"
mv "$working_root" "$result_root"
printf 'completed\t%s\n' "$(date --iso-8601=seconds)" > "$log_root/completed.tsv"
