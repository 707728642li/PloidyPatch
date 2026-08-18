#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then echo "usage: $0 PROJECT_ROOT" >&2; exit 2; fi
project_root=$(realpath "$1")
result_root=$project_root/data/validation/apple_golden_delicious_ont_cra021523_v0.1
working_root=${result_root}.working
contract=$working_root/metadata/file_contract.tsv
accession=CRR1429917
output=$working_root/files/$accession.fastq.gz
audit=$working_root/parallel_resume_contract.tsv
attempts=$working_root/parallel_transport_attempts.tsv

[[ ! -e $result_root && -d $working_root && -s $contract ]] || {
    echo "refusing to overwrite or accelerate an ambiguous ONT bundle" >&2; exit 1;
}
[[ ! -e $audit && ! -e $attempts ]] || {
    echo "refusing to overwrite parallel transport audit" >&2; exit 1;
}
row=$(awk -F '\t' -v accession="$accession" '$1 == accession {print}' "$contract")
[[ -n $row && $(printf '%s\n' "$row" | wc -l) -eq 1 ]] || {
    echo "final accession is not unique in the frozen contract" >&2; exit 1;
}
IFS=$'\t' read -r observed_accession tissue expected_md5 url <<< "$row"
[[ $observed_accession == "$accession" && $expected_md5 =~ ^[0-9a-f]{32}$ ]] || {
    echo "malformed final accession contract" >&2; exit 1;
}
[[ $url == "https://download.cncb.ac.cn/gsa4/CRA021523/$accession/$accession.fastq.gz" ]] || {
    echo "final accession source changed" >&2; exit 1;
}
{
    printf 'field\tvalue\n'
    printf 'code_commit\t%s\n' "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'accession\t%s\n' "$accession"
    printf 'reason\toverlap_final_accession_with_rate_limited_penultimate_accession\n'
    printf 'source_changed\tfalse\nchecksums_changed\tfalse\n'
    printf 'candidate_coordinates_or_ranks_modified\tfalse\n'
    printf 'parallel_files_total\t2\nconnections\t8\n'
    printf 'helper_script_sha256\t%s\n' \
        "$(sha256sum "$project_root/code/scripts/parallel_resume_apple_ont_final_accession_v0.1.sh" | awk '{print $1}')"
} > "$audit"
printf 'accession\tattempt\texit_code\n' > "$attempts"

attempt=0
while true; do
    attempt=$((attempt + 1))
    if axel -n 8 -o "$output" "$url" \
        >> "$working_root/logs/$accession.parallel.axel.log" 2>&1; then
        printf '%s\t%s\t0\n' "$accession" "$attempt" >> "$attempts"
        break
    else
        status=$?
    fi
    printf '%s\t%s\t%s\n' "$accession" "$attempt" "$status" >> "$attempts"
    [[ $attempt -lt 20 ]] || {
        echo "$accession exhausted parallel transport retries" >&2; exit 1;
    }
    sleep 2
done
observed_md5=$(md5sum "$output" | awk '{print $1}')
[[ $observed_md5 == "$expected_md5" ]] || {
    echo "$accession MD5 mismatch after parallel resume" >&2; exit 1;
}
gzip -t "$output"
printf 'parallel final accession complete: %s\n' "$accession"
