#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 3 ]]; then
    echo "usage: $0 ENV_PREFIX OUTPUT_DIR MAMBA_TIME_LOG" >&2
    exit 2
fi

env_prefix=$(realpath "$1")
output_dir=$(realpath -m "$2")
mamba_time_log=$(realpath "$3")

if [[ ! -x $env_prefix/bin/syngap ]]; then
    echo "SynGAP executable is absent: $env_prefix/bin/syngap" >&2
    exit 2
fi
if [[ -e $output_dir ]]; then
    echo "refusing to overwrite environment audit: $output_dir" >&2
    exit 2
fi
mkdir -p "$output_dir"

conda list -p "$env_prefix" --explicit > "$output_dir/conda_explicit.txt"
conda list -p "$env_prefix" > "$output_dir/conda_list.txt"
conda list -p "$env_prefix" --json > "$output_dir/conda_list.json"

shopt -s nullglob
metadata=("$env_prefix"/conda-meta/syngap-1.2.5-*.json)
if [[ ${#metadata[@]} -ne 1 ]]; then
    echo "expected exactly one SynGAP 1.2.5 conda metadata record" >&2
    exit 1
fi
cp "${metadata[0]}" "$output_dir/syngap_conda_metadata.json"

conda run -p "$env_prefix" --no-capture-output syngap --help \
    > "$output_dir/syngap_help.txt" \
    2> "$output_dir/syngap_help.stderr.txt"
conda run -p "$env_prefix" --no-capture-output python --version \
    > "$output_dir/python_version.txt" 2>&1
conda run -p "$env_prefix" --no-capture-output bedtools --version \
    > "$output_dir/bedtools_version.txt" 2>&1
conda run -p "$env_prefix" --no-capture-output gffread --version \
    > "$output_dir/gffread_version.txt" 2>&1
conda run -p "$env_prefix" --no-capture-output seqkit version \
    > "$output_dir/seqkit_version.txt" 2>&1
conda run -p "$env_prefix" --no-capture-output diamond version \
    > "$output_dir/diamond_version.txt" 2>&1
conda run -p "$env_prefix" --no-capture-output lastal --version \
    > "$output_dir/lastal_version.txt" 2>&1
conda run -p "$env_prefix" --no-capture-output needle -version \
    > "$output_dir/needle_version.txt" 2>&1

sha256sum "$env_prefix/bin/syngap" "${metadata[0]}" \
    > "$output_dir/entrypoint_sha256.txt"
find -L "$env_prefix" \
    \( -path '*/genBlast_v138_linux_x86_64/genblast_v138_linux_x86_64' \
       -o -path '*/CPC2.py' \
       -o -path '*/SwissProt.tar.gz' \) \
    -type f | sort > "$output_dir/bundled_resource_paths.txt"

for required_pattern in \
    genBlast_v138_linux_x86_64/genblast_v138_linux_x86_64 \
    CPC2.py SwissProt.tar.gz; do
    if ! grep -Fq "$required_pattern" "$output_dir/bundled_resource_paths.txt"; then
        echo "missing bundled SynGAP resource: $required_pattern" >&2
        exit 1
    fi
done
while IFS= read -r resource; do
    sha256sum "$resource"
done < "$output_dir/bundled_resource_paths.txt" \
    > "$output_dir/bundled_resource_sha256.txt"

cp "$mamba_time_log" "$output_dir/mamba_create.time.txt"
{
    printf 'artifact\tbytes\tsha256\n'
    for artifact in "$output_dir"/*; do
        [[ -f $artifact ]] || continue
        [[ $(basename "$artifact") != audit_artifacts.tsv ]] || continue
        artifact_sha=$(sha256sum "$artifact" | awk '{print $1}')
        printf '%s\t%s\t%s\n' \
            "$(basename "$artifact")" "$(stat -c %s "$artifact")" "$artifact_sha"
    done
} > "$output_dir/audit_artifacts.tsv"

printf 'SynGAP environment audit passed: %s\n' "$output_dir"
