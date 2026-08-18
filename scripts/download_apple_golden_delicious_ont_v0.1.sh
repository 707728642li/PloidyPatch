#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then echo "usage: $0 PROJECT_ROOT" >&2; exit 2; fi
project_root=$(realpath "$1")
rank_root=$project_root/results/natural/apple_gddh13_v0.4/discovery/rankings
result_root=$project_root/data/validation/apple_golden_delicious_ont_cra021523_v0.1
working_root=${result_root}.working

[[ -s $rank_root/SHA256SUMS && -s $rank_root/run_contract.tsv ]] || {
    echo "apple natural candidates and ranks are not frozen" >&2; exit 1;
}
(cd "$rank_root" && sha256sum -c SHA256SUMS >/dev/null)
grep -q $'^candidate_and_rank_freeze_precedes_validation_access\ttrue$' \
    "$rank_root/run_contract.tsv" || {
    echo "candidate/RNA evidence firewall is not declared" >&2; exit 1;
}
grep -q $'^RNA_access\tfalse$' "$rank_root/run_contract.tsv" || {
    echo "candidate ranks were not declared RNA-blind" >&2; exit 1;
}
[[ ! -e $result_root && ! -e $working_root ]] || {
    echo "refusing to overwrite apple ONT validation bundle" >&2; exit 1;
}
mkdir -p "$working_root"/{files,metadata,logs}
cat > "$working_root/metadata/file_contract.tsv" <<'EOF'
accession	tissue	md5	url
CRR1429911	root	9aa0ed42d22c1d243108451a61bd3f2f	https://download.cncb.ac.cn/gsa4/CRA021523/CRR1429911/CRR1429911.fastq.gz
CRR1429912	stem	58a9a32399c531cdcc8977f5cf240c4b	https://download.cncb.ac.cn/gsa4/CRA021523/CRR1429912/CRR1429912.fastq.gz
CRR1429913	leaf	75d1775b932bcdb8b7a3b79456757c15	https://download.cncb.ac.cn/gsa4/CRA021523/CRR1429913/CRR1429913.fastq.gz
CRR1429914	flower	c54b598838f64d6096e63512e0cf4bd0	https://download.cncb.ac.cn/gsa4/CRA021523/CRR1429914/CRR1429914.fastq.gz
CRR1429915	young_fruit_45DAB	083cd900dff4c7cc6e7ba30639d71b20	https://download.cncb.ac.cn/gsa4/CRA021523/CRR1429915/CRR1429915.fastq.gz
CRR1429916	expanding_fruit_95DAB	16f391099c84c8f8178100f935967c03	https://download.cncb.ac.cn/gsa4/CRA021523/CRR1429916/CRR1429916.fastq.gz
CRR1429917	mature_fruit_145DAB	319e7135b89fbc2acf22ef3e589bf684	https://download.cncb.ac.cn/gsa4/CRA021523/CRR1429917/CRR1429917.fastq.gz
EOF
curl -L -sS https://download.cncb.ac.cn/gsa4/CRA021523/md5sum.txt \
    > "$working_root/metadata/official_md5sum.txt"
for accession in CRR1429911 CRR1429912 CRR1429913 CRR1429914 \
                 CRR1429915 CRR1429916 CRR1429917; do
    grep -q "$accession/$accession.fastq.gz" \
        "$working_root/metadata/official_md5sum.txt" || {
        echo "official GSA checksum list lacks $accession" >&2; exit 1;
    }
done
while IFS=$'\t' read -r accession tissue expected_md5 url; do
    [[ $accession == accession ]] && continue
    official_md5=$(awk -v accession="$accession" \
        '$2 ~ ("/" accession "/" accession "[.]fastq[.]gz$") {print $1}' \
        "$working_root/metadata/official_md5sum.txt")
    [[ $official_md5 =~ ^[0-9a-f]{32}$ && $official_md5 == "$expected_md5" ]] || {
        echo "file contract differs from official GSA MD5 for $accession" >&2; exit 1;
    }
done < "$working_root/metadata/file_contract.tsv"
{
    printf 'field\tvalue\n'
    printf 'code_commit\t%s\n' \
        "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'download_script_sha256\t%s\n' \
        "$(sha256sum "$project_root/code/scripts/download_apple_golden_delicious_ont_v0.1.sh" | awk '{print $1}')"
    printf 'source\tGSA_CRA021523\nplatform\tOxford_Nanopore\n'
    printf 'cultivar\tGolden_Delicious\ntissues\t7\n'
    printf 'candidate_and_rank_freeze_precedes_download\ttrue\n'
    printf 'rank_freeze_sha256sums_sha256\t%s\n' \
        "$(sha256sum "$rank_root/SHA256SUMS" | awk '{print $1}')"
    printf 'download_connections_per_file\t8\nparallel_files\t7\n'
} > "$working_root/run_contract.tsv"

pids=()
while IFS=$'\t' read -r accession tissue expected_md5 url; do
    [[ $accession == accession ]] && continue
    (
        axel -n 8 -o "$working_root/files/$accession.fastq.gz" "$url" \
            > "$working_root/logs/$accession.axel.log" 2>&1
        observed=$(md5sum "$working_root/files/$accession.fastq.gz" | awk '{print $1}')
        [[ $observed == "$expected_md5" ]] || {
            echo "$accession MD5 mismatch" >&2; exit 1;
        }
        gzip -t "$working_root/files/$accession.fastq.gz"
    ) &
    pids+=("$!")
done < "$working_root/metadata/file_contract.tsv"
failed=0
for pid in "${pids[@]}"; do
    wait "$pid" || failed=1
done
[[ $failed == 0 ]] || { echo "one or more apple ONT downloads failed" >&2; exit 1; }

{
    printf 'accession\ttissue\tbytes\tmd5\tsha256\tpath\n'
    while IFS=$'\t' read -r accession tissue expected_md5 url; do
        [[ $accession == accession ]] && continue
        path=$working_root/files/$accession.fastq.gz
        printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$accession" "$tissue" \
            "$(stat -Lc %s "$path")" "$expected_md5" \
            "$(sha256sum "$path" | awk '{print $1}')" "$path"
    done < "$working_root/metadata/file_contract.tsv"
} > "$working_root/input_manifest.tsv"
(
    cd "$working_root"
    find . -type f ! -name SHA256SUMS -print0 | sort -z | xargs -0 sha256sum \
        > SHA256SUMS
    sha256sum -c SHA256SUMS >/dev/null
)
du -sb "$working_root" > "$working_root/disk_bytes.txt"
chmod -R a-w "$working_root"
mv "$working_root" "$result_root"
printf 'apple Golden Delicious ONT validation bundle frozen: %s\n' "$result_root"
