#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then echo "usage: $0 PROJECT_ROOT" >&2; exit 2; fi
project_root=$(realpath "$1")
code_root=$project_root/code
python_bin=$project_root/envs/ploidypatch-dev/bin/python
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
validation_parent=$project_root/results/natural/apple_gddh13_v0.4/validation
failed_root=$validation_parent/ont_raw_v0.1_failed_csv_field_limit
patch1_root=$validation_parent/ont_raw_v0.1_patch1
result_root=$validation_parent/ont_raw_v0.1_patch2_transcript_strand
working_root=${result_root}.working
code_commit=${PLOIDYPATCH_CODE_COMMIT:-}

[[ $code_commit =~ ^[0-9a-f]{40}$ ]] || {
    echo "PLOIDYPATCH_CODE_COMMIT must be a full git SHA" >&2; exit 1;
}
verify_tree() { (cd "$1" && sha256sum -c SHA256SUMS >/dev/null); }
for required in "$python_bin" "$genome" "$genome.fai" "$base_gff" \
    "$candidate_gff" "$rankings" "$repeat_gff" \
    "$ont_root/metadata/file_contract.tsv" \
    "$failed_root/SHA256SUMS" "$patch1_root/SHA256SUMS" \
    "$patch1_root/self_map/candidate_cds_to_genome.paf" \
    "$rank_root/SHA256SUMS" "$ont_root/SHA256SUMS" \
    "$te_root/SHA256SUMS" "$method_root/SHA256SUMS"; do
    [[ -s $required ]] || { echo "missing patch2 input: $required" >&2; exit 1; }
done
[[ ! -e $result_root && ! -e $working_root ]] || {
    echo "refusing to overwrite apple transcript-strand patch2" >&2; exit 1;
}
for root in "$failed_root" "$patch1_root" "$rank_root" "$ont_root" \
    "$te_root" "$method_root"; do
    verify_tree "$root"
done

mkdir -p "$working_root"/{alignment/combined,evidence,audit,self_map,logs,freeze}
{
    printf 'field\tvalue\n'
    printf 'code_commit\t%s\n' "$code_commit"
    printf 'target\tMalus_domestica_GDDH13_v1.1\n'
    printf 'validation_source\tGolden_Delicious_GSA_CRA021523_raw_ONT\n'
    printf 'patch_reason\tcorrect_PAF_query_orientation_vs_transcript_strand_semantics\n'
    printf 'labels_seen\ttrue\n'
    printf 'posthoc_semantic_correction\ttrue\n'
    printf 'candidate_coordinates_modified\tfalse\n'
    printf 'candidate_ranks_modified\tfalse\n'
    printf 'evidence_thresholds_modified\tfalse\n'
    printf 'raw_reads_remapped\tfalse\n'
    printf 'frozen_raw_PAF_reused\ttrue\n'
    printf 'alignment_strand_source\tminimap2_ts\n'
    printf 'alignment_strand_formula\treference_transcript_strand=paf_query_target_strand*ts_query_transcript_relation\n'
    printf 'unspliced_ts_missing_policy\texclude_from_strand_specific_support\n'
    printf 'minimum_query_coverage\t0.85\nminimum_identity\t0.90\n'
    printf 'minimum_mapq\t20\nmaximum_secondary_score_fraction\t0.95\n'
    printf 'minimum_candidate_cds_coverage\t0.90\nflank_bp\t5000\n'
    printf 'minimum_case_study_full_chain_reads\t2\n'
    printf 'review_budgets\t100,146,250,292,500,583\n'
    printf 'bootstrap_replicates\t20000\nbootstrap_seed\t20261004\n'
    printf 'automatic_annotation_patch\tfalse\n'
    printf 'patch1_status\tdescriptive_invalid_transcript_strand_interpretation\n'
    printf 'upstream_reference\thttps://github.com/lh3/minimap2/blob/master/cookbook.md\n'
} > "$working_root/run_contract.tsv"
{
    printf 'role\tbytes\tsha256\tpath\n'
    for entry in "genome:$genome" "genome_fai:$genome.fai" \
        "base_gff:$base_gff" "candidate_gff:$candidate_gff" \
        "review_rankings:$rankings" "repeat_gff:$repeat_gff" \
        "failed_attempt_freeze:$failed_root/SHA256SUMS" \
        "patch1_freeze:$patch1_root/SHA256SUMS" \
        "reused_self_map:$patch1_root/self_map/candidate_cds_to_genome.paf" \
        "rank_freeze:$rank_root/SHA256SUMS" \
        "ont_freeze:$ont_root/SHA256SUMS" "te_freeze:$te_root/SHA256SUMS"; do
        role=${entry%%:*}; path=${entry#*:}
        printf '%s\t%s\t%s\t%s\n' "$role" "$(stat -Lc %s "$path")" \
            "$(sha256sum "$path" | awk '{print $1}')" "$path"
    done
} > "$working_root/input_manifest.tsv"
{
    printf 'path\tsha256\n'
    for relative in scripts/resume_apple_ont_transcript_strand_patch2_v0.1.sh \
        scripts/run_apple_golden_delicious_ont_validation_v0.1.sh \
        src/ploidypatch/isoseq_validation.py src/ploidypatch/natural_audit.py \
        src/ploidypatch/cli.py; do
        printf '%s\t%s\n' "$relative" \
            "$(sha256sum "$code_root/$relative" | awk '{print $1}')"
    done
} > "$working_root/freeze/code_manifest.tsv"

cd "$code_root"
combined=$working_root/alignment/combined/golden_delicious_ont.transcript_strand_candidate_query_universe.paf
filter_args=()
while IFS=$'\t' read -r accession _; do
    [[ $accession == accession ]] && continue
    filter_args+=(--paf-input "$accession=$failed_root/alignment/by_run/$accession.paf")
done < "$ont_root/metadata/file_contract.tsv"
/usr/bin/time -v -o "$working_root/logs/query_filter.time.txt" \
    "$python_bin" -m ploidypatch.cli evidence filter-candidate-query-paf \
    --candidate-gff "$candidate_gff" "${filter_args[@]}" \
    --alignment-strand-source minimap2_ts \
    --output-paf "$combined" \
    --output-counts "$working_root/evidence/raw_read_counts.tsv" \
    --output-summary "$working_root/alignment/combined/query_filter_summary.tsv" \
    --output-manifest "$working_root/alignment/combined/query_filter_manifest.json" \
    > "$working_root/alignment/combined/query_filter.stdout.json" \
    2> "$working_root/alignment/combined/query_filter.stderr.log"

/usr/bin/time -v -o "$working_root/logs/validation.time.txt" \
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

cp --reflink=auto "$patch1_root/self_map/candidate_cds_to_genome.paf" \
    "$working_root/self_map/candidate_cds_to_genome.paf"
[[ $(sha256sum "$working_root/self_map/candidate_cds_to_genome.paf" | awk '{print $1}') == \
   $(sha256sum "$patch1_root/self_map/candidate_cds_to_genome.paf" | awk '{print $1}') ]] || {
    echo "reused candidate self-map changed" >&2; exit 1;
}
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

"$python_bin" - "$failed_root/evidence/candidate_ont_evidence.tsv" \
    "$working_root/evidence/candidate_ont_evidence.tsv" \
    "$patch1_root/evidence/review_yield.json" \
    "$working_root/evidence/review_yield.json" \
    "$patch1_root/audit/summary.json" "$working_root/audit/summary.json" \
    "$working_root/evidence/patch1_vs_patch2.json" <<'PY'
import csv
import hashlib
import json
import sys
from collections import Counter

csv.field_size_limit(2**31 - 1)
old_path, new_path, old_review_path, new_review_path, old_audit_path, new_audit_path, output = sys.argv[1:]

def sha(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()

def read_evidence(path):
    with open(path, encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return {row["candidate_digest"]: row for row in rows}

old = read_evidence(old_path)
new = read_evidence(new_path)
if set(old) != set(new):
    raise SystemExit("patch2 candidate universe differs from patch1")
transitions = Counter((old[digest]["evidence_state"], new[digest]["evidence_state"]) for digest in old)
old_full = {digest for digest, row in old.items() if row["evidence_state"] == "full_chain_supported"}
new_full = {digest for digest, row in new.items() if row["evidence_state"] == "full_chain_supported"}
old_reads2 = {digest for digest in old_full if int(old[digest]["supporting_full_length_reads"]) >= 2}
new_reads2 = {digest for digest in new_full if int(new[digest]["supporting_full_length_reads"]) >= 2}
old_review = json.load(open(old_review_path, encoding="utf-8"))
new_review = json.load(open(new_review_path, encoding="utf-8"))
old_audit = json.load(open(old_audit_path, encoding="utf-8"))
new_audit = json.load(open(new_audit_path, encoding="utf-8"))
report = {
    "schema_version": "ploidypatch.apple_ont_strand_patch_comparison.v1",
    "scope": "descriptive_posthoc_semantic_correction",
    "inputs": {
        "patch1_evidence_sha256": sha(old_path),
        "patch2_evidence_sha256": sha(new_path),
        "patch1_review_sha256": sha(old_review_path),
        "patch2_review_sha256": sha(new_review_path),
        "patch1_audit_sha256": sha(old_audit_path),
        "patch2_audit_sha256": sha(new_audit_path),
    },
    "counts": {
        "candidate_universe": len(old),
        "patch1_full_chain": len(old_full),
        "patch2_full_chain": len(new_full),
        "full_chain_intersection": len(old_full & new_full),
        "patch1_full_chain_reads_ge_2": len(old_reads2),
        "patch2_full_chain_reads_ge_2": len(new_reads2),
        "reads_ge_2_intersection": len(old_reads2 & new_reads2),
        "patch1_case_study_ready": old_audit["counts"]["case_study_ready"],
        "patch2_case_study_ready": new_audit["counts"]["case_study_ready"],
    },
    "state_transitions": {
        f"{left}->{right}": count
        for (left, right), count in sorted(transitions.items())
    },
    "review_primary": {
        "patch1": old_review["primary"],
        "patch2": new_review["primary"],
    },
}
with open(output, "x", encoding="utf-8", newline="") as handle:
    json.dump(report, handle, indent=2, sort_keys=True)
    handle.write("\n")
PY

"$python_bin" - "$failed_root/SHA256SUMS" \
    "$working_root/alignment/combined/query_filter_manifest.json" \
    "$working_root/evidence/candidate_ont_evidence.tsv.manifest.json" \
    "$working_root/audit/summary.json" <<'PY'
import json
import sys

source_sums, filter_path, evidence_path, audit_path = sys.argv[1:]
expected = {}
with open(source_sums, encoding="utf-8") as handle:
    for line in handle:
        digest, relative = line.rstrip("\n").split("  ", 1)
        if relative.startswith("./alignment/by_run/") and relative.endswith(".paf"):
            expected[relative.rsplit("/", 1)[-1][:-4]] = digest
query_filter = json.load(open(filter_path, encoding="utf-8"))
observed = {row["accession"]: row["sha256"] for row in query_filter["inputs"]["paf"]}
if expected != observed or len(observed) != 7:
    raise SystemExit("patch2 raw PAF hashes do not match frozen alignment source")
if query_filter["parameters"]["alignment_strand_source"] != "minimap2_ts":
    raise SystemExit("patch2 filter did not use minimap2 ts semantics")
if query_filter["counts"]["strand_unavailable_alignments"] < 1:
    raise SystemExit("patch2 did not expose unstranded alignments")
evidence = json.load(open(evidence_path, encoding="utf-8"))
if evidence["schema_version"] != "ploidypatch.isoseq_candidate_validation.v2":
    raise SystemExit("patch2 evidence schema is not v2")
if evidence["parameters"]["alignment_strand_source"] != "minimap2_ts":
    raise SystemExit("patch2 validator did not use minimap2 ts semantics")
if evidence["counts"]["candidate_models"] != 29144:
    raise SystemExit("patch2 candidate universe changed")
if evidence["counts"]["selected_transcripts"] != query_filter["counts"]["retained_queries"]:
    raise SystemExit("patch2 filter/validator read universes differ")
audit = json.load(open(audit_path, encoding="utf-8"))
if audit["counts"]["candidates"] != 29144:
    raise SystemExit("patch2 audit candidate universe changed")
if audit["parameters"]["minimum_full_length_read_support"] != 2:
    raise SystemExit("patch2 case threshold changed")
PY

for required in \
    "$working_root/alignment/combined/query_filter_manifest.json" \
    "$working_root/evidence/candidate_ont_evidence.tsv" \
    "$working_root/evidence/candidate_ont_evidence.tsv.manifest.json" \
    "$working_root/evidence/review_yield.json" \
    "$working_root/evidence/bootstrap.json" \
    "$working_root/evidence/patch1_vs_patch2.json" \
    "$working_root/self_map/candidate_cds_to_genome.paf" \
    "$working_root/audit/candidate_biological_audit.tsv" \
    "$working_root/audit/summary.json"; do
    [[ -s $required ]] || { echo "missing patch2 output: $required" >&2; exit 1; }
done
du -sb "$working_root" > "$working_root/disk_bytes.txt"
(
    cd "$working_root"
    find . -type f ! -name SHA256SUMS -print0 | sort -z | \
        xargs -0 sha256sum > SHA256SUMS
    sha256sum -c SHA256SUMS >/dev/null
)
chmod -R a-w "$working_root"
mv "$working_root" "$result_root"
printf 'apple ONT transcript-strand patch2 frozen: %s\n' "$result_root"
