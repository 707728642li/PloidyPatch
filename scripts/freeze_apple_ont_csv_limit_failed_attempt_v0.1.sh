#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then echo "usage: $0 PROJECT_ROOT" >&2; exit 2; fi
project_root=$(realpath "$1")
validation_root=$project_root/results/natural/apple_gddh13_v0.4/validation
working_root=$validation_root/ont_raw_v0.1.working
failed_root=$validation_root/ont_raw_v0.1_failed_csv_field_limit
evidence=$working_root/evidence/candidate_ont_evidence.tsv
manifest=${evidence}.manifest.json
stderr=$working_root/evidence/review_join.stderr.log

for required in "$evidence" "$manifest" "$stderr" \
    "$working_root/alignment/combined/golden_delicious_ont.candidate_query_universe.paf" \
    "$working_root/alignment/combined/query_filter_manifest.tsv"; do
    [[ -s $required ]] || { echo "missing failed-attempt artifact: $required" >&2; exit 1; }
done
[[ ! -e $failed_root && ! -e $validation_root/ont_raw_v0.1 ]] || {
    echo "refusing to overwrite apple ONT failed-attempt freeze" >&2; exit 1;
}
grep -q 'field larger than field limit (131072)' "$stderr" || {
    echo "failed attempt does not contain the expected csv limit error" >&2; exit 1;
}
[[ $(wc -l < "$evidence") -eq 29145 ]] || {
    echo "failed attempt evidence is incomplete" >&2; exit 1;
}
{
    printf 'field\tvalue\n'
    printf 'formal_outcome\texecution_failure_before_rank_summary\n'
    printf 'failure_stage\tjoin_isoseq_review_rankings\n'
    printf 'failure_reason\tpython_csv_default_field_limit_131072\n'
    printf 'raw_alignment_complete\ttrue\nstrict_chain_evidence_complete\ttrue\n'
    printf 'labels_accessed\ttrue\nscientific_thresholds_changed\tfalse\n'
    printf 'candidate_coordinates_or_ranks_modified\tfalse\n'
    printf 'automatic_annotation_patch\tfalse\n'
    printf 'required_retry\tparser_only_patch_new_nonoverwriting_result_root\n'
    printf 'original_code_commit\t1736241f04a7e2456c8ecd11f55cbf82a13fb287\n'
    printf 'candidate_evidence_sha256\t%s\n' "$(sha256sum "$evidence" | awk '{print $1}')"
    printf 'candidate_evidence_manifest_sha256\t%s\n' "$(sha256sum "$manifest" | awk '{print $1}')"
    printf 'freeze_script_sha256\t%s\n' \
        "$(sha256sum "$project_root/code/scripts/freeze_apple_ont_csv_limit_failed_attempt_v0.1.sh" | awk '{print $1}')"
} > "$working_root/failure_manifest.tsv"
(
    cd "$working_root"
    find . -type f ! -name SHA256SUMS -print0 | sort -z | xargs -0 sha256sum \
        > SHA256SUMS
    sha256sum -c SHA256SUMS >/dev/null
)
chmod -R a-w "$working_root"
mv "$working_root" "$failed_root"
printf 'apple ONT csv-limit failed attempt frozen: %s\n' "$failed_root"
