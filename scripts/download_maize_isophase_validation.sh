#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then echo "usage: $0 PROJECT_ROOT" >&2; exit 2; fi
project_root=$(realpath "$1")
freeze=$project_root/results/natural/maize_v2/discovery/homeolog_ranker
out=$project_root/data/validation/maize_isophase_zenodo_2611319
working=${out}.working

[[ -s $freeze/SHA256SUMS ]] || {
    echo "maize natural candidate/rank freeze is missing" >&2; exit 1;
}
(cd "$freeze" && sha256sum -c SHA256SUMS >/dev/null)
grep -q $'^candidate_freeze_precedes_validation_access\ttrue$' \
    "$freeze/run_contract.tsv" || {
    echo "candidate/rank freeze does not declare the evidence firewall" >&2
    exit 1
}
if [[ -e $out || -e $working ]]; then
    echo "refusing to overwrite maize Iso-Seq validation data" >&2; exit 1
fi
mkdir -p "$working/files"

cat > "$working/source_contract.tsv" <<'EOF'
filename	bytes	zenodo_md5	url
F1maize.FINAL.fasta	187170721	3011fb13e8bb06f03e17971a5538a476	https://zenodo.org/api/records/2611319/files/F1maize.FINAL.fasta/content
F1maize.FINAL.demux_FL_count.txt	2653173	b894b7e18834525ae7e5e6fe7d9c59b8	https://zenodo.org/api/records/2611319/files/F1maize.FINAL.demux_FL_count.txt/content
F1maize.FINAL.gff	59650292	8c3816625c69e9da526616cc0170e8f1	https://zenodo.org/api/records/2611319/files/F1maize.FINAL.gff/content
EOF

printf 'field\tvalue\n' > "$working/download_status.tsv"
printf 'candidate_freeze_sha256\t%s\n' \
    "$(sha256sum "$freeze/SHA256SUMS" | awk '{print $1}')" \
    >> "$working/download_status.tsv"
printf 'download_started\t%s\ncode_commit\t%s\n' \
    "$(date --iso-8601=seconds)" \
    "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}" \
    >> "$working/download_status.tsv"

tail -n +2 "$working/source_contract.tsv" |
while IFS=$'\t' read -r name expected_bytes expected_md5 url; do
    partial=$working/files/$name.partial
    axel -n 8 -o "$partial" "$url"
    observed_bytes=$(stat -Lc %s "$partial")
    observed_md5=$(md5sum "$partial" | awk '{print $1}')
    [[ $observed_bytes == "$expected_bytes" ]] || {
        echo "size mismatch for $name" >&2; exit 1;
    }
    [[ $observed_md5 == "$expected_md5" ]] || {
        echo "MD5 mismatch for $name" >&2; exit 1;
    }
    mv "$partial" "$working/files/$name"
done

{
    printf 'filename\tbytes\tmd5\tsha256\n'
    tail -n +2 "$working/source_contract.tsv" |
    while IFS=$'\t' read -r name expected_bytes expected_md5 url; do
        path=$working/files/$name
        printf '%s\t%s\t%s\t%s\n' "$name" "$(stat -Lc %s "$path")" \
            "$(md5sum "$path" | awk '{print $1}')" \
            "$(sha256sum "$path" | awk '{print $1}')"
    done
} > "$working/file_manifest.tsv"
printf 'download_complete\t%s\n' "$(date --iso-8601=seconds)" \
    >> "$working/download_status.tsv"
(
    cd "$working"
    find . -type f ! -name SHA256SUMS -print0 | sort -z | xargs -0 sha256sum \
        > SHA256SUMS
    sha256sum -c SHA256SUMS >/dev/null
)
mv "$working" "$out"
printf 'maize Iso-Seq validation data frozen: %s\n' "$out"
