#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then echo "usage: $0 PROJECT_ROOT" >&2; exit 2; fi
project_root=$(realpath "$1")
rank_root=$project_root/results/natural/apple_gddh13_v0.4/discovery/rankings
result_root=$project_root/data/validation/apple_golden_delicious_ont_cra021523_v0.1
working_root=${result_root}.working
contract=$working_root/metadata/file_contract.tsv

[[ ! -e $result_root && -d $working_root ]] || {
    echo "refusing to overwrite or ambiguously resume apple ONT bundle" >&2; exit 1;
}
for required in "$rank_root/SHA256SUMS" "$rank_root/run_contract.tsv" \
    "$contract" "$working_root/metadata/official_md5sum.txt"; do
    [[ -s $required ]] || { echo "missing resume input: $required" >&2; exit 1; }
done
(cd "$rank_root" && sha256sum -c SHA256SUMS >/dev/null)
grep -q $'^candidate_and_rank_freeze_precedes_validation_access\ttrue$' \
    "$rank_root/run_contract.tsv" || {
    echo "candidate/RNA evidence firewall is not declared" >&2; exit 1;
}
if [[ ! -e $working_root/resume_contract.tsv ]]; then
    {
        printf 'field\tvalue\n'
        printf 'code_commit\t%s\n' "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
        printf 'reason\tparallel_axel_TLS_connections_closed\n'
        printf 'source_changed\tfalse\nchecksums_changed\tfalse\n'
        printf 'resume_transport\taxel_st_state\n'
        printf 'download_connections_per_file\t8\nparallel_files\t1\n'
        printf 'candidate_coordinates_or_ranks_modified\tfalse\n'
        printf 'resume_script_sha256\t%s\n' \
            "$(sha256sum "$project_root/code/scripts/resume_apple_golden_delicious_ont_v0.1.sh" | awk '{print $1}')"
    } > "$working_root/resume_contract.tsv"
else
    grep -q $'^source_changed\tfalse$' "$working_root/resume_contract.tsv"
    grep -q $'^checksums_changed\tfalse$' "$working_root/resume_contract.tsv"
    grep -q $'^candidate_coordinates_or_ranks_modified\tfalse$' \
        "$working_root/resume_contract.tsv"
fi
{
    printf '%s\t%s\t%s\n' "$(date -u +%FT%TZ)" \
        "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}" \
        "$(sha256sum "$project_root/code/scripts/resume_apple_golden_delicious_ont_v0.1.sh" | awk '{print $1}')"
} >> "$working_root/resume_invocations.tsv"
[[ -e $working_root/transport_attempts.tsv ]] || \
    printf 'accession\tattempt\texit_code\n' > "$working_root/transport_attempts.tsv"

while IFS=$'\t' read -r accession tissue expected_md5 url; do
    [[ $accession == accession ]] && continue
    output=$working_root/files/$accession.fastq.gz
    state=${output}.st
    complete=false
    if [[ -s $output && ! -e $state ]]; then
        observed=$(md5sum "$output" | awk '{print $1}')
        if [[ $observed == "$expected_md5" ]] && gzip -t "$output"; then
            complete=true
        fi
    fi
    if [[ $complete == false ]]; then
        attempt=0
        while true; do
            attempt=$((attempt + 1))
            if axel -n 8 -o "$output" "$url" \
                >> "$working_root/logs/$accession.resume.axel.log" 2>&1; then
                printf '%s\t%s\t0\n' "$accession" "$attempt" \
                    >> "$working_root/transport_attempts.tsv"
                break
            else
                status=$?
            fi
            printf '%s\t%s\t%s\n' "$accession" "$attempt" "$status" \
                >> "$working_root/transport_attempts.tsv"
            [[ $attempt -lt 20 ]] || {
                echo "$accession exhausted transport retries" >&2; exit 1;
            }
            sleep 2
        done
    fi
    observed=$(md5sum "$output" | awk '{print $1}')
    [[ $observed == "$expected_md5" ]] || {
        echo "$accession MD5 mismatch after resume" >&2; exit 1;
    }
    gzip -t "$output"
done < "$contract"

{
    printf 'accession\ttissue\tbytes\tmd5\tsha256\tpath\n'
    while IFS=$'\t' read -r accession tissue expected_md5 url; do
        [[ $accession == accession ]] && continue
        path=$working_root/files/$accession.fastq.gz
        printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$accession" "$tissue" \
            "$(stat -Lc %s "$path")" "$expected_md5" \
            "$(sha256sum "$path" | awk '{print $1}')" "$path"
    done < "$contract"
} > "$working_root/input_manifest.tsv"
du -sb "$working_root" > "$working_root/disk_bytes.txt"
(
    cd "$working_root"
    find . -type f ! -name SHA256SUMS -print0 | sort -z | xargs -0 sha256sum \
        > SHA256SUMS
    sha256sum -c SHA256SUMS >/dev/null
)
chmod -R a-w "$working_root"
mv "$working_root" "$result_root"
printf 'apple Golden Delicious ONT validation bundle resumed and frozen: %s\n' \
    "$result_root"
