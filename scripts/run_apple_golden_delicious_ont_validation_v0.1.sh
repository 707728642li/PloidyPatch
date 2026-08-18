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
ont_root=$project_root/data/validation/apple_golden_delicious_ont_cra021523_v0.1
te_root=$project_root/data/validation/apple_gddh13_te_v0.1
repeat_gff=$te_root/files/GDDH13_1-1_TE.gff3
result_root=$project_root/results/natural/apple_gddh13_v0.4/validation/ont_raw_v0.1
working_root=${result_root}.working

verify_tree() { (cd "$1" && sha256sum -c SHA256SUMS >/dev/null); }
for required in "$python_bin" "$minimap2" "$genome" "$genome.fai" \
    "$base_gff" "$candidate_gff" "$rankings" "$rank_root/SHA256SUMS" \
    "$ont_root/SHA256SUMS" "$ont_root/metadata/file_contract.tsv" \
    "$te_root/SHA256SUMS" "$repeat_gff" "$method_root/SHA256SUMS"; do
    [[ -s $required ]] || { echo "missing apple ONT validation input: $required" >&2; exit 1; }
done
for root in "$rank_root" "$ont_root" "$te_root" "$method_root"; do
    verify_tree "$root"
done
grep -q $'^candidate_and_rank_freeze_precedes_validation_access\ttrue$' \
    "$rank_root/run_contract.tsv" || {
    echo "candidate/RNA evidence firewall is absent" >&2; exit 1;
}
[[ ! -e $result_root && ! -e $working_root ]] || {
    echo "refusing to overwrite apple raw ONT validation" >&2; exit 1;
}
mkdir -p "$working_root"/{alignment/by_run,alignment/combined,evidence,queries,self_map,audit,logs,freeze}
{
    printf 'field\tvalue\ncode_commit\t%s\n' \
        "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'target\tMalus_domestica_GDDH13_v1.1\n'
    printf 'validation_source\tGolden_Delicious_GSA_CRA021523_raw_ONT\n'
    printf 'candidate_and_rank_freeze_verified_before_read_access\ttrue\n'
    printf 'candidate_coordinates_or_ranks_modified\tfalse\n'
    printf 'alignment_preset\tsplice-k14\nsecondary_alignments\t10\n'
    printf 'shared_minimap2_index\ttrue\nindex_builds\t1\n'
    printf 'threads_per_run\t16\nparallel_runs\t7\n'
    printf 'minimum_query_coverage\t0.85\nminimum_identity\t0.90\n'
    printf 'minimum_mapq\t20\nmaximum_secondary_score_fraction\t0.95\n'
    printf 'minimum_candidate_cds_coverage\t0.90\nflank_bp\t5000\n'
    printf 'alignment_strand_source\tminimap2_ts\n'
    printf 'alignment_strand_formula\treference_transcript_strand=paf_query_target_strand*ts_query_transcript_relation\n'
    printf 'unspliced_ts_missing_policy\texclude_from_strand_specific_support\n'
    printf 'candidate_query_prefilter\tlossless_all_alignments_for_any_query_with_transcript_strand_overlapping_a_candidate_span\n'
    printf 'minimum_case_study_full_chain_reads\t2\n'
    printf 'single_exon_is_not_full_chain_positive\ttrue\n'
    printf 'repeat_annotation\tGDDH13_v1.1_official_TE_GFF\n'
    printf 'review_budgets\t100,146,250,292,500,583\n'
    printf 'automatic_annotation_patch\tfalse\n'
    printf 'minimap2_version\t%s\n' "$("$minimap2" --version)"
    printf 'rank_freeze_sha256sums_sha256\t%s\n' \
        "$(sha256sum "$rank_root/SHA256SUMS" | awk '{print $1}')"
    printf 'ont_freeze_sha256sums_sha256\t%s\n' \
        "$(sha256sum "$ont_root/SHA256SUMS" | awk '{print $1}')"
    printf 'te_freeze_sha256sums_sha256\t%s\n' \
        "$(sha256sum "$te_root/SHA256SUMS" | awk '{print $1}')"
} > "$working_root/run_contract.tsv"
{
    printf 'role\tbytes\tsha256\tpath\n'
    for entry in "genome:$genome" "genome_fai:$genome.fai" \
        "base_gff:$base_gff" "candidate_gff:$candidate_gff" \
        "review_rankings:$rankings" "repeat_gff:$repeat_gff" \
        "rank_freeze:$rank_root/SHA256SUMS" "ont_freeze:$ont_root/SHA256SUMS" \
        "te_freeze:$te_root/SHA256SUMS"; do
        role=${entry%%:*}; path=${entry#*:}
        printf '%s\t%s\t%s\t%s\n' "$role" "$(stat -Lc %s "$path")" \
            "$(sha256sum "$path" | awk '{print $1}')" "$path"
    done
} > "$working_root/input_manifest.tsv"
{
    printf 'path\tsha256\n'
    for relative in scripts/run_apple_golden_delicious_ont_validation_v0.1.sh \
        src/ploidypatch/isoseq_validation.py src/ploidypatch/natural_audit.py \
        src/ploidypatch/cli.py; do
        printf '%s\t%s\n' "$relative" \
            "$(sha256sum "$code_root/$relative" | awk '{print $1}')"
    done
} > "$working_root/freeze/code_manifest.tsv"

genome_index=$working_root/alignment/gddh13_primary.splice_k14.mmi
/usr/bin/time -v -o "$working_root/logs/genome_index.time.txt" \
    "$minimap2" -x splice -k14 -d "$genome_index" "$genome" \
    > "$working_root/logs/genome_index.stdout.log" \
    2> "$working_root/logs/genome_index.stderr.log"
[[ -s $genome_index ]] || { echo "apple minimap2 index is empty" >&2; exit 1; }

pids=()
while IFS=$'\t' read -r accession tissue expected_md5 url; do
    [[ $accession == accession ]] && continue
    input=$ont_root/files/$accession.fastq.gz
    output=$working_root/alignment/by_run/$accession.paf
    [[ -s $input ]] || { echo "missing ONT FASTQ: $input" >&2; exit 1; }
    (
        /usr/bin/time -v -o "$working_root/logs/$accession.minimap2.time.txt" \
            "$minimap2" -t 16 -x splice -k14 --secondary=yes -N 10 \
            -c --cs=long "$genome_index" "$input" > "$output" \
            2> "$working_root/logs/$accession.minimap2.stderr.log"
        [[ -s $output ]] || { echo "$accession produced an empty PAF" >&2; exit 1; }
    ) &
    pids+=("$!")
done < "$ont_root/metadata/file_contract.tsv"
failed=0
for pid in "${pids[@]}"; do
    wait "$pid" || failed=1
done
[[ $failed == 0 ]] || { echo "one or more ONT alignments failed" >&2; exit 1; }

combined=$working_root/alignment/combined/golden_delicious_ont.candidate_query_universe.paf
cd "$code_root"
filter_args=()
while IFS=$'\t' read -r accession _; do
    [[ $accession == accession ]] && continue
    filter_args+=(--paf-input "$accession=$working_root/alignment/by_run/$accession.paf")
done < "$ont_root/metadata/file_contract.tsv"
"$python_bin" -m ploidypatch.cli evidence filter-candidate-query-paf \
    --candidate-gff "$candidate_gff" "${filter_args[@]}" \
    --alignment-strand-source minimap2_ts \
    --output-paf "$combined" \
    --output-counts "$working_root/evidence/raw_read_counts.tsv" \
    --output-summary "$working_root/alignment/combined/query_filter_summary.tsv" \
    --output-manifest "$working_root/alignment/combined/query_filter_manifest.json" \
    > "$working_root/alignment/combined/query_filter.stdout.json" \
    2> "$working_root/alignment/combined/query_filter.stderr.log"
"$python_bin" -m ploidypatch.cli evidence validate-isoseq-candidates \
    --candidate-gff "$candidate_gff" --paf "$combined" \
    --selected-counts "$working_root/evidence/raw_read_counts.tsv" \
    --genome-fasta "$genome" --minimum-query-coverage 0.85 \
    --minimum-identity 0.90 --minimum-mapq 20 \
    --maximum-secondary-score-fraction 0.95 \
    --minimum-candidate-cds-coverage 0.90 --flank-bp 5000 \
    --alignment-strand-source minimap2_ts \
    --output-evidence "$working_root/evidence/candidate_ont_evidence.tsv" \
    > "$working_root/evidence/validation.stdout.json" \
    2> "$working_root/evidence/validation.stderr.log"
"$python_bin" -m ploidypatch.cli evidence join-isoseq-review-rankings \
    --evidence "$working_root/evidence/candidate_ont_evidence.tsv" \
    --review-rankings "$rankings" \
    --review-budget 100 --review-budget 146 --review-budget 250 \
    --review-budget 292 --review-budget 500 --review-budget 583 \
    --comparator-estimator baseline --primary-estimator v04_guard \
    --output-tsv "$working_root/evidence/ranked_ont_evidence.tsv" \
    --output-summary "$working_root/evidence/review_yield.json" \
    > "$working_root/evidence/review_join.stdout.json" \
    2> "$working_root/evidence/review_join.stderr.log"
"$python_bin" -m ploidypatch.cli evidence bootstrap-isoseq-review-yield \
    --evidence "$working_root/evidence/candidate_ont_evidence.tsv" \
    --review-rankings "$rankings" \
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
    --isoseq-evidence "$working_root/evidence/candidate_ont_evidence.tsv" \
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

for required in "$working_root/evidence/candidate_ont_evidence.tsv" \
    "$working_root/evidence/candidate_ont_evidence.tsv.manifest.json" \
    "$working_root/evidence/review_yield.json" \
    "$working_root/evidence/bootstrap.json" \
    "$working_root/alignment/combined/query_filter_manifest.json" \
    "$genome_index" \
    "$working_root/self_map/candidate_cds_to_genome.paf" \
    "$working_root/audit/candidate_biological_audit.tsv" \
    "$working_root/audit/summary.json"; do
    [[ -s $required ]] || { echo "missing apple ONT validation output: $required" >&2; exit 1; }
done
"$python_bin" - "$working_root/audit/summary.json" \
    "$working_root/evidence/candidate_ont_evidence.tsv.manifest.json" <<'PY'
import json
import sys

audit = json.load(open(sys.argv[1], encoding="utf-8"))
evidence = json.load(open(sys.argv[2], encoding="utf-8"))
if audit["counts"]["candidates"] != 29144:
    raise SystemExit("unexpected apple natural candidate count")
if audit["parameters"]["minimum_full_length_read_support"] != 2:
    raise SystemExit("raw-read case threshold changed")
if evidence["parameters"]["selected_count_field"] != "full_length_reads":
    raise SystemExit("raw-read evidence did not use generic per-read counts")
if evidence["parameters"]["alignment_strand_source"] != "minimap2_ts":
    raise SystemExit("raw unstranded cDNA evidence did not use minimap2 ts")
if evidence["parameters"]["single_exon_is_not_full_chain_positive"] is not True:
    raise SystemExit("single-exon support boundary changed")
PY
(
    cd "$working_root"
    find . -type f ! -name SHA256SUMS -print0 | sort -z | xargs -0 sha256sum \
        > SHA256SUMS
    sha256sum -c SHA256SUMS >/dev/null
)
du -sb "$working_root" > "$working_root/disk_bytes.txt"
chmod -R a-w "$working_root"
mv "$working_root" "$result_root"
printf 'apple Golden Delicious raw ONT validation frozen: %s\n' "$result_root"
