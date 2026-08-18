#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "usage: $0 PROJECT_ROOT" >&2
    exit 2
fi

project_root=$(realpath "$1")
code_root=$project_root/code
python_bin=$project_root/envs/ploidypatch-dev/bin/python
samtools_bin=$project_root/envs/ploidypatch-pav/bin/samtools
audit=$project_root/results/hawthorn/rna_audit_v0.1/bam_audit.tsv
primary_root=$project_root/results/hawthorn/rna_junctions_v0.1
result_root=$project_root/results/hawthorn/rna_junctions_secondary_v0.1
working_root=${result_root}.working
log_root=$project_root/logs/hawthorn/rna_junctions_secondary_v0.1

for required in "$python_bin" "$samtools_bin" "$audit" \
                "$primary_root/primary_aggregate.tsv"; do
    if [[ ! -s $required ]]; then
        echo "missing secondary-junction prerequisite: $required" >&2
        exit 1
    fi
done
if [[ -e $result_root || -e $working_root ]]; then
    echo "refusing to overwrite secondary junction result" >&2
    exit 1
fi
mkdir -p "$working_root/samples" "$log_root"

awk -F '\t' 'BEGIN { OFS="\t"; print "sample", "bam_path" }
    NR > 1 && $1 != "Black-1" && $1 != "Black-2" && $1 != "Black-3" {
        print $1, $13
    }
' "$audit" > "$working_root/sample_bams.tsv"
sample_count=$(($(wc -l < "$working_root/sample_bams.tsv") - 1))
if [[ $sample_count -ne 57 ]]; then
    echo "expected 57 secondary samples, found $sample_count" >&2
    exit 1
fi

{
    printf 'field\tvalue\n'
    printf 'samples\t%s\n' "$sample_count"
    printf 'concurrent_samples\t8\n'
    printf 'samtools_threads_per_sample\t4\n'
    printf 'min_mapq\t20\n'
    printf 'excluded_flag_mask\t2308\n'
    printf 'strand_policy\tunstranded\n'
    printf 'negative_evidence_policy\tabsence_is_missing_not_contradiction\n'
    printf 'samtools_sha256\t%s\n' \
        "$(sha256sum "$samtools_bin" | awk '{print $1}')"
} > "$working_root/run_contract.tsv"

pids=()
pid_samples=()
status=0
while IFS=$'\t' read -r sample bam; do
    if [[ $sample == sample ]]; then
        continue
    fi
    if [[ ! $sample =~ ^[A-Za-z0-9][A-Za-z0-9_.-]*$ ]]; then
        echo "unsafe sample identifier in RNA audit: $sample" >&2
        exit 1
    fi
    if [[ ! -s $bam || ! -s ${bam}.bai ]]; then
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
    pid_samples+=("$sample")
    if [[ ${#pids[@]} -ge 8 ]]; then
        if ! wait "${pids[0]}"; then
            echo "junction extraction failed: ${pid_samples[0]}" >&2
            status=1
        fi
        pids=("${pids[@]:1}")
        pid_samples=("${pid_samples[@]:1}")
    fi
done < "$working_root/sample_bams.tsv"
for index in "${!pids[@]}"; do
    if ! wait "${pids[$index]}"; then
        echo "junction extraction failed: ${pid_samples[$index]}" >&2
        status=1
    fi
done
if [[ $status -ne 0 ]]; then
    echo "preserving failed secondary junction work: $working_root" >&2
    exit 1
fi

cd "$code_root"
/usr/bin/time -v -o "$working_root/aggregate.resource.time.txt" \
    "$python_bin" -m ploidypatch.cli evidence aggregate-junctions \
        --input-dir "$primary_root/samples" \
        --input-dir "$working_root/samples" \
        --primary-sample Black-1 \
        --primary-sample Black-2 \
        --primary-sample Black-3 \
        --output-tsv "$working_root/all_samples_aggregate.tsv" \
        --min-reads-per-sample 2 \
        --min-supporting-samples 2 \
        > "$working_root/aggregate.stdout.log" \
        2> "$working_root/aggregate.stderr.log"

secondary_outputs=$(find "$working_root/samples" -type f \
    -name '*.junctions.tsv' | wc -l)
if [[ $secondary_outputs -ne 57 || \
      ! -s $working_root/all_samples_aggregate.tsv ]]; then
    echo "secondary or aggregate junction outputs are incomplete" >&2
    exit 1
fi
sha256sum "$working_root"/samples/*.junctions.tsv \
    "$working_root"/samples/*.manifest.json \
    "$working_root"/all_samples_aggregate.tsv \
    "$working_root"/all_samples_aggregate.tsv.manifest.json \
    > "$working_root/SHA256SUMS"
du -sb "$working_root" > "$working_root/disk_bytes.txt"
mv "$working_root" "$result_root"
printf 'hawthorn secondary junction evidence completed: %s\n' "$result_root"
