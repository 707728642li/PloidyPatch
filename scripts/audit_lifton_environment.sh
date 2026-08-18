#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 5 ]]; then
    echo "usage: $0 ENV_PREFIX OUTPUT_DIR SDIST MAMBA_TIME_LOG PIP_TIME_LOG" >&2
    exit 2
fi

env_prefix=$(realpath "$1")
output_dir=$(realpath -m "$2")
sdist=$(realpath "$3")
mamba_time_log=$(realpath "$4")
pip_time_log=$(realpath "$5")
expected_sdist_sha=c2125db9bede3640e13c6aa6a0672887aaf2611710442f9bf69017197d188f88

if [[ ! -x $env_prefix/bin/lifton || ! -x $env_prefix/bin/minimap2 ||
      ! -x $env_prefix/bin/miniprot ]]; then
    echo "LiftOn or a required aligner executable is absent" >&2
    exit 2
fi
if ! sha256sum "$sdist" | grep -Fq "$expected_sdist_sha"; then
    echo "LiftOn source archive hash mismatch" >&2
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
conda run -p "$env_prefix" --no-capture-output pip freeze --all \
    > "$output_dir/pip_freeze.txt"
conda run -p "$env_prefix" --no-capture-output pip inspect --local \
    > "$output_dir/pip_inspect.json"
conda run -p "$env_prefix" --no-capture-output pip check \
    > "$output_dir/pip_check.txt"

conda run -p "$env_prefix" --no-capture-output lifton -V \
    > "$output_dir/lifton_version.txt" 2>&1
conda run -p "$env_prefix" --no-capture-output lifton -h \
    > "$output_dir/lifton_help.txt" 2> "$output_dir/lifton_help.stderr.txt"
conda run -p "$env_prefix" --no-capture-output python -c \
    'import lifton; print(lifton.__version__)' \
    > "$output_dir/lifton_import_version.txt" 2>&1
conda run -p "$env_prefix" --no-capture-output python --version \
    > "$output_dir/python_version.txt" 2>&1
conda run -p "$env_prefix" --no-capture-output minimap2 --version \
    > "$output_dir/minimap2_version.txt" 2>&1
conda run -p "$env_prefix" --no-capture-output miniprot --version \
    > "$output_dir/miniprot_version.txt" 2>&1

for required_flag in \
    --gene-only --no-miniprot-candidate --no-adaptive-rescue-floor \
    --no-miniprot-rescue --legacy-merge --full-dp-align \
    --locus-pipeline --validate-output; do
    if ! grep -Fq -- "$required_flag" "$output_dir/lifton_help.txt"; then
        echo "LiftOn help lacks expected v1.0.11 flag: $required_flag" >&2
        exit 1
    fi
done
if ! grep -Fq 'No broken requirements found' "$output_dir/pip_check.txt"; then
    echo "pip dependency audit did not pass" >&2
    exit 1
fi
if ! grep -Fq '1.0.11' "$output_dir/lifton_version.txt" ||
   ! grep -Fq '1.0.11' "$output_dir/lifton_import_version.txt"; then
    echo "LiftOn CLI/import versions disagree with the frozen release" >&2
    exit 1
fi

shopt -s nullglob
dist_info=("$env_prefix"/lib/python3.11/site-packages/lifton-1.0.11.dist-info)
if [[ ${#dist_info[@]} -ne 1 ]]; then
    echo "expected exactly one LiftOn 1.0.11 dist-info directory" >&2
    exit 1
fi
cp "${dist_info[0]}/METADATA" "$output_dir/lifton_METADATA.txt"
cp "${dist_info[0]}/RECORD" "$output_dir/lifton_RECORD.csv"
if [[ -s ${dist_info[0]}/direct_url.json ]]; then
    cp "${dist_info[0]}/direct_url.json" "$output_dir/lifton_direct_url.json"
fi

minimap_metadata=("$env_prefix"/conda-meta/minimap2-2.30-*.json)
miniprot_metadata=("$env_prefix"/conda-meta/miniprot-0.18-*.json)
if [[ ${#minimap_metadata[@]} -ne 1 || ${#miniprot_metadata[@]} -ne 1 ]]; then
    echo "expected one frozen minimap2 and miniprot conda record" >&2
    exit 1
fi
cp "${minimap_metadata[0]}" "$output_dir/minimap2_conda_metadata.json"
cp "${miniprot_metadata[0]}" "$output_dir/miniprot_conda_metadata.json"
cp "$mamba_time_log" "$output_dir/mamba_create.time.txt"
cp "$pip_time_log" "$output_dir/pip_install.time.txt"

sha256sum "$sdist" "$env_prefix/bin/lifton" \
    "$env_prefix/lib/python3.11/site-packages/lifton/lifton.py" \
    "${minimap_metadata[0]}" "${miniprot_metadata[0]}" \
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

printf 'LiftOn environment audit passed: %s\n' "$output_dir"
