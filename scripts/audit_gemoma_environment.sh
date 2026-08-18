#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 3 ]]; then
    echo "usage: $0 ENV_PREFIX OUTPUT_DIR MAMBA_TIME_LOG" >&2
    exit 2
fi

env_prefix=$(realpath "$1")
output_dir=$(realpath -m "$2")
mamba_time_log=$(realpath "$3")
jar=$env_prefix/share/gemoma-1.9-0/GeMoMa-1.9.jar

if [[ ! -x $env_prefix/bin/GeMoMa || ! -s $jar ]]; then
    echo "GeMoMa executable or jar is absent" >&2
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
metadata=("$env_prefix"/conda-meta/gemoma-1.9-*.json)
if [[ ${#metadata[@]} -ne 1 ]]; then
    echo "expected exactly one GeMoMa 1.9 conda metadata record" >&2
    exit 1
fi
cp "${metadata[0]}" "$output_dir/gemoma_conda_metadata.json"
cp "$mamba_time_log" "$output_dir/mamba_create.time.txt"

conda run -p "$env_prefix" --no-capture-output GeMoMa -h \
    > "$output_dir/gemoma_help.txt" 2> "$output_dir/gemoma_help.stderr.txt"
conda run -p "$env_prefix" --no-capture-output java -version \
    > "$output_dir/java_version.txt" 2>&1
conda run -p "$env_prefix" --no-capture-output tblastn -version \
    > "$output_dir/tblastn_version.txt" 2>&1
conda run -p "$env_prefix" --no-capture-output mmseqs version \
    > "$output_dir/mmseqs_version.txt" 2>&1
conda run -p "$env_prefix" --no-capture-output python --version \
    > "$output_dir/python_version.txt" 2>&1

sha256sum "$env_prefix/bin/GeMoMa" "$jar" "${metadata[0]}" \
    > "$output_dir/entrypoint_sha256.txt"
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

printf 'GeMoMa environment audit passed: %s\n' "$output_dir"
