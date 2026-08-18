#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then echo "usage: $0 PROJECT_ROOT" >&2; exit 2; fi
project_root=$(realpath "$1")
code_root=$project_root/code
python_bin=$project_root/envs/ploidypatch-dev/bin/python
minimap2=$project_root/envs/ploidypatch-pav/bin/minimap2
bundle=$project_root/data/derived/external_inputs/apple_v0.3/target_apple
genome=$bundle/primary_chromosomes.genome.fa
base_gff=$bundle/primary_chromosomes.gff3
method_root=$project_root/results/copy_collapse/external/apple_v0.3_method_trio
candidate_gff=$method_root/consensus/primary_union/complete_control/candidate.gff3
rank_root=$project_root/results/natural/apple_gddh13_v0.4/discovery/rankings
rankings=$rank_root/natural/review_rankings.tsv
te_root=$project_root/data/validation/apple_gddh13_te_v0.1
repeat_gff=$te_root/files/GDDH13_1-1_TE.gff3
failed_root=$project_root/results/natural/apple_gddh13_v0.4/validation/ont_raw_v0.1_failed_csv_field_limit
evidence=$failed_root/evidence/candidate_ont_evidence.tsv
result_root=$project_root/results/natural/apple_gddh13_v0.4/validation/ont_raw_v0.1_patch1
working_root=${result_root}.working

for required in "$python_bin" "$minimap2" "$genome" "$genome.fai" \
    "$base_gff" "$candidate_gff" "$rankings" "$repeat_gff" \
    "$failed_root/SHA256SUMS" "$failed_root/failure_manifest.tsv" \
    "$evidence" "$evidence.manifest.json"; do
    [[ -s $required ]] || { echo "missing patch1 input: $required" >&2; exit 1; }
done
(cd "$failed_root" && sha256sum -c SHA256SUMS >/dev/null)
grep -q $'^required_retry\tparser_only_patch_new_nonoverwriting_result_root$' \
    "$failed_root/failure_manifest.tsv"
grep -q $'^labels_accessed\ttrue$' "$failed_root/failure_manifest.tsv"
[[ ! -e $result_root && ! -e $working_root ]] || {
    echo "refusing to overwrite apple ONT patch1 result" >&2; exit 1;
}
mkdir -p "$working_root"/{evidence,queries,self_map,audit,logs,freeze}
{
    printf 'field\tvalue\n'
    printf 'code_commit\t%s\n' "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'retry_type\tparser_only_post_alignment_patch\n'
    printf 'labels_seen_before_patch\ttrue\n'
    printf 'candidate_coordinates_or_ranks_modified\tfalse\n'
    printf 'alignment_recomputed\tfalse\nstrict_chain_evidence_recomputed\tfalse\n'
    printf 'scientific_thresholds_changed\tfalse\nreview_budgets_changed\tfalse\n'
    printf 'minimum_case_study_full_chain_reads\t2\n'
    printf 'review_budgets\t100,146,250,292,500,583\n'
    printf 'bootstrap_replicates\t20000\nbootstrap_seed\t20261004\n'
    printf 'primary_estimator\tv04_guard\ncomparator_estimator\tbaseline\n'
    printf 'automatic_annotation_patch\tfalse\n'
    printf 'failed_attempt_sha256sums_sha256\t%s\n' \
        "$(sha256sum "$failed_root/SHA256SUMS" | awk '{print $1}')"
} > "$working_root/run_contract.tsv"
{
    printf 'role\tbytes\tsha256\tpath\n'
    for entry in "failed_freeze:$failed_root/SHA256SUMS" \
        "candidate_evidence:$evidence" "evidence_manifest:$evidence.manifest.json" \
        "genome:$genome" "genome_fai:$genome.fai" "base_gff:$base_gff" \
        "candidate_gff:$candidate_gff" "review_rankings:$rankings" \
        "repeat_gff:$repeat_gff"; do
        role=${entry%%:*}; path=${entry#*:}
        printf '%s\t%s\t%s\t%s\n' "$role" "$(stat -Lc %s "$path")" \
            "$(sha256sum "$path" | awk '{print $1}')" "$path"
    done
} > "$working_root/input_manifest.tsv"
{
    printf 'path\tsha256\n'
    for relative in scripts/resume_apple_ont_post_alignment_patch1_v0.1.sh \
        src/ploidypatch/isoseq_validation.py src/ploidypatch/natural_audit.py \
        src/ploidypatch/cli.py; do
        printf '%s\t%s\n' "$relative" \
            "$(sha256sum "$code_root/$relative" | awk '{print $1}')"
    done
} > "$working_root/freeze/code_manifest.tsv"

cd "$code_root"
"$python_bin" -m ploidypatch.cli evidence join-isoseq-review-rankings \
    --evidence "$evidence" --review-rankings "$rankings" \
    --review-budget 100 --review-budget 146 --review-budget 250 \
    --review-budget 292 --review-budget 500 --review-budget 583 \
    --comparator-estimator baseline --primary-estimator v04_guard \
    --output-tsv "$working_root/evidence/ranked_ont_evidence.tsv" \
    --output-summary "$working_root/evidence/review_yield.json" \
    > "$working_root/evidence/review_join.stdout.json" \
    2> "$working_root/evidence/review_join.stderr.log"
"$python_bin" -m ploidypatch.cli evidence bootstrap-isoseq-review-yield \
    --evidence "$evidence" --review-rankings "$rankings" \
    --review-budget 100 --review-budget 146 --review-budget 250 \
    --review-budget 292 --review-budget 500 --review-budget 583 \
    --comparator-estimator baseline --primary-estimator v04_guard \
    --replicates 20000 --seed 20261004 \
    --output-json "$working_root/evidence/bootstrap.json" \
    > "$working_root/evidence/bootstrap.stdout.json" \
    2> "$working_root/evidence/bootstrap.stderr.log"
"$python_bin" -m ploidypatch.cli evidence export-natural-candidate-cds \
    --candidate-gff "$candidate_gff" --genome-fasta "$genome" \
    --output-fasta "$working_root/queries/candidate_cds.fa" \
    > "$working_root/queries/stdout.json" 2> "$working_root/queries/stderr.log"
/usr/bin/time -v -o "$working_root/logs/candidate_self_map.time.txt" \
    "$minimap2" -x splice:hq -t 64 --secondary=yes -N 20 -c \
    "$genome" "$working_root/queries/candidate_cds.fa" \
    > "$working_root/self_map/candidate_cds_to_genome.paf" \
    2> "$working_root/self_map/minimap2.stderr.log"
"$python_bin" -m ploidypatch.cli evidence audit-natural-candidates \
    --candidate-gff "$candidate_gff" --base-gff "$base_gff" \
    --genome-fasta "$genome" --review-rankings "$rankings" \
    --isoseq-evidence "$evidence" \
    --self-map-paf "$working_root/self_map/candidate_cds_to_genome.paf" \
    --repeat-gff "$repeat_gff" --repeat-flank-bp 2000 \
    --minimum-full-length-read-support 2 \
    --review-budget 100 --review-budget 146 --review-budget 250 \
    --review-budget 292 --review-budget 500 --review-budget 583 \
    --minimum-query-coverage 0.90 --minimum-identity 0.98 \
    --near-equal-score-fraction 0.95 \
    --output-tsv "$working_root/audit/candidate_biological_audit.tsv" \
    --output-summary "$working_root/audit/summary.json" \
    > "$working_root/audit/stdout.json" 2> "$working_root/audit/stderr.log"

for required in "$working_root/evidence/review_yield.json" \
    "$working_root/evidence/bootstrap.json" \
    "$working_root/self_map/candidate_cds_to_genome.paf" \
    "$working_root/audit/candidate_biological_audit.tsv" \
    "$working_root/audit/summary.json"; do
    [[ -s $required ]] || { echo "missing apple ONT patch1 output: $required" >&2; exit 1; }
done
"$python_bin" - "$working_root/audit/summary.json" \
    "$working_root/evidence/review_yield.json" <<'PY'
import json
import sys

audit = json.load(open(sys.argv[1], encoding="utf-8"))
review = json.load(open(sys.argv[2], encoding="utf-8"))
if audit["counts"]["candidates"] != 29144:
    raise SystemExit("unexpected apple patch1 candidate count")
if audit["parameters"]["minimum_full_length_read_support"] != 2:
    raise SystemExit("raw-read case threshold changed")
delta = review["primary"]["estimator_delta"]
if delta["primary_estimator"] != "v04_guard" or delta["comparator_estimator"] != "baseline":
    raise SystemExit("patch1 review comparison changed")
PY
(
    cd "$working_root"
    find . -type f ! -name SHA256SUMS -print0 | sort -z | xargs -0 sha256sum \
        > SHA256SUMS
    sha256sum -c SHA256SUMS >/dev/null
)
chmod -R a-w "$working_root"
mv "$working_root" "$result_root"
printf 'apple Golden Delicious raw ONT patch1 validation frozen: %s\n' "$result_root"
