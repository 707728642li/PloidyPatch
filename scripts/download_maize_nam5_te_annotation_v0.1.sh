#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then echo "usage: $0 PROJECT_ROOT" >&2; exit 2; fi
project_root=$(realpath "$1")
result_root=$project_root/data/validation/maize_nam5_te_v0.1
working_root=${result_root}.working
url=https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0.TE.gff3.gz
file_name=Zm-B73-REFERENCE-NAM-5.0.TE.gff3.gz

if [[ -e $result_root || -e $working_root ]]; then
    echo "refusing to overwrite maize NAM-5.0 TE annotation" >&2; exit 1
fi
mkdir -p "$working_root/files"
axel -n 8 -o "$working_root/files/$file_name" "$url"
gzip -t "$working_root/files/$file_name"
{
    printf 'field\tvalue\n'
    printf 'source_url\t%s\n' "$url"
    printf 'download_tool\taxel_-n_8\n'
    printf 'downloaded_bytes\t%s\n' "$(stat -Lc %s "$working_root/files/$file_name")"
    printf 'downloaded_sha256\t%s\n' \
        "$(sha256sum "$working_root/files/$file_name" | awk '{print $1}')"
    printf 'gzip_integrity\tpass\n'
    printf 'role\tofficial_repeat_context_not_candidate_generation_or_truth\n'
} > "$working_root/source.tsv"
(
    cd "$working_root"
    sha256sum "files/$file_name" source.tsv > SHA256SUMS
    sha256sum -c SHA256SUMS >/dev/null
)
du -sb "$working_root" > "$working_root/disk_bytes.txt"
mv "$working_root" "$result_root"
printf 'maize NAM-5.0 TE annotation frozen: %s\n' "$result_root"
