#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then echo "usage: $0 PROJECT_ROOT" >&2; exit 2; fi
project_root=$(realpath "$1")
rank_root=$project_root/results/natural/apple_gddh13_v0.4/discovery/rankings
result_root=$project_root/data/validation/apple_gddh13_te_v0.1
working_root=${result_root}.working
url=https://iris.angers.inra.fr/gddh13/downloads/GDDH13_1-1_TE.gff3.bz2
expected_bytes=34251728

[[ -s $rank_root/SHA256SUMS && -s $rank_root/run_contract.tsv ]] || {
    echo "apple natural rankings are not frozen" >&2; exit 1;
}
(cd "$rank_root" && sha256sum -c SHA256SUMS >/dev/null)
grep -q $'^RNA_access\tfalse$' "$rank_root/run_contract.tsv" || {
    echo "candidate ranks were not declared RNA-blind" >&2; exit 1;
}
[[ ! -e $result_root && ! -e $working_root ]] || {
    echo "refusing to overwrite apple TE validation bundle" >&2; exit 1;
}
mkdir -p "$working_root"/{files,logs}
axel -n 8 -o "$working_root/files/GDDH13_1-1_TE.gff3.bz2" "$url" \
    > "$working_root/logs/axel.log" 2>&1
[[ $(stat -Lc %s "$working_root/files/GDDH13_1-1_TE.gff3.bz2") == $expected_bytes ]] || {
    echo "apple TE download byte count differs" >&2; exit 1;
}
bzip2 -t "$working_root/files/GDDH13_1-1_TE.gff3.bz2"
bzip2 -dc "$working_root/files/GDDH13_1-1_TE.gff3.bz2" \
    > "$working_root/files/GDDH13_1-1_TE.gff3"
[[ -s $working_root/files/GDDH13_1-1_TE.gff3 ]] || {
    echo "decompressed apple TE GFF is empty" >&2; exit 1;
}
{
    printf 'field\tvalue\nsource_url\t%s\n' "$url"
    printf 'source_last_modified\tTue,_17_Oct_2017_09:09:14_GMT\n'
    printf 'source_etag\t59e5c8ba-20aa3d0\nsource_content_length\t%s\n' \
        "$expected_bytes"
    printf 'candidate_and_rank_freeze_precedes_validation_access\ttrue\n'
    printf 'rank_freeze_sha256sums_sha256\t%s\n' \
        "$(sha256sum "$rank_root/SHA256SUMS" | awk '{print $1}')"
} > "$working_root/run_contract.tsv"
{
    printf 'role\tbytes\tsha256\tpath\n'
    for entry in \
        "compressed_te_gff:$working_root/files/GDDH13_1-1_TE.gff3.bz2" \
        "te_gff:$working_root/files/GDDH13_1-1_TE.gff3"; do
        role=${entry%%:*}; path=${entry#*:}
        printf '%s\t%s\t%s\t%s\n' "$role" "$(stat -Lc %s "$path")" \
            "$(sha256sum "$path" | awk '{print $1}')" "$path"
    done
} > "$working_root/input_manifest.tsv"
(
    cd "$working_root"
    find . -type f ! -name SHA256SUMS -print0 | sort -z | xargs -0 sha256sum \
        > SHA256SUMS
    sha256sum -c SHA256SUMS >/dev/null
)
chmod -R a-w "$working_root"
mv "$working_root" "$result_root"
printf 'apple GDDH13 TE annotation frozen: %s\n' "$result_root"

