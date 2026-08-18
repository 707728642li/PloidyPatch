#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then echo "usage: $0 PROJECT_ROOT" >&2; exit 2; fi
project_root=$(realpath "$1")
code_root=$project_root/code
python_bin=$project_root/envs/ploidypatch-dev/bin/python
minimap2_bin=$project_root/envs/ploidypatch-pav/bin/minimap2
bundle=$project_root/data/derived/holdout_inputs/maize_v2/zea_mays
bundle_root=$project_root/data/derived/holdout_inputs/maize_v2
genome=$bundle/primary_chromosomes.genome.fa
base_gff=$bundle/primary_chromosomes.gff3
method_root=$project_root/results/natural/maize_v2/discovery/method_trio
candidate_gff=$method_root/consensus/union/natural/candidate.gff3
rank_root=$project_root/results/natural/maize_v2/discovery/homeolog_ranker
rankings=$rank_root/natural/review_rankings.tsv
isoseq_root=$project_root/results/natural/maize_v2/validation/isoseq_v0.1
isoseq=$isoseq_root/evidence/candidate_isoseq_evidence.tsv
repeat_root=$project_root/data/validation/maize_nam5_te_v0.1
repeat_gff=$repeat_root/files/Zm-B73-REFERENCE-NAM-5.0.TE.gff3.gz
repeat_seqid_map=$code_root/config/maize_nam5_te_seqid_map_v0.1.tsv
result_root=$project_root/results/natural/maize_v2/validation/biological_audit_v0.1
working_root=${result_root}.working

for frozen in "$bundle_root" "$method_root" "$rank_root" "$isoseq_root" "$repeat_root"; do
    [[ -s $frozen/SHA256SUMS ]] || { echo "unfrozen biological-audit input: $frozen" >&2; exit 1; }
    (cd "$frozen" && sha256sum -c SHA256SUMS >/dev/null)
done
for required in "$python_bin" "$minimap2_bin" "$genome" "$genome.fai" \
                "$base_gff" "$candidate_gff" "$rankings" "$isoseq" \
                "$repeat_gff" "$repeat_seqid_map"; do
    [[ -s $required ]] || { echo "missing biological-audit input: $required" >&2; exit 1; }
done
if [[ -e $result_root || -e $working_root ]]; then
    echo "refusing to overwrite maize natural biological audit" >&2; exit 1
fi
mkdir -p "$working_root"/{queries,self_map,audit,logs}
{
    printf 'field\tvalue\n'
    printf 'code_commit\t%s\n' "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'split\tnatural_current_annotation_post_isoseq_audit_v0.1\n'
    printf 'candidate_freeze_precedes_validation_access\ttrue\n'
    printf 'isoseq_truth_access\ttrue\nself_map_truth_access\tfalse\n'
    printf 'minimum_query_coverage\t0.90\nminimum_identity\t0.98\n'
    printf 'near_equal_score_fraction\t0.95\n'
    printf 'repeat_annotation\tNAM-5.0_official_TE_GFF\nrepeat_flank_bp\t2000\n'
    printf 'review_budgets\t25,50,100,200\n'
    printf 'minimap2_preset\tsplice:hq\nminimap2_threads\t64\n'
    printf 'minimap2_version\t%s\n' "$("$minimap2_bin" --version)"
    printf 'automatic_approval\tfalse\n'
} > "$working_root/run_contract.tsv"
{
    printf 'role\tbytes\tsha256\tpath\n'
    for entry in "genome:$genome" "genome_fai:$genome.fai" \
                 "base_gff:$base_gff" "candidate_gff:$candidate_gff" \
                 "review_rankings:$rankings" "isoseq_evidence:$isoseq" \
                 "repeat_gff:$repeat_gff" "repeat_seqid_map:$repeat_seqid_map"; do
        role=${entry%%:*}; path=${entry#*:}
        printf '%s\t%s\t%s\t%s\n' "$role" "$(stat -Lc %s "$path")" \
            "$(sha256sum "$path" | awk '{print $1}')" "$path"
    done
} > "$working_root/input_manifest.tsv"

cd "$code_root"
"$python_bin" -m ploidypatch.cli evidence export-natural-candidate-cds \
    --candidate-gff "$candidate_gff" --genome-fasta "$genome" \
    --output-fasta "$working_root/queries/candidate_cds.fa" \
    > "$working_root/queries/stdout.json" \
    2> "$working_root/queries/stderr.log"
/usr/bin/time -v -o "$working_root/logs/minimap2.time.txt" \
    "$minimap2_bin" -x splice:hq -t 64 --secondary=yes -N 20 -c \
        "$genome" "$working_root/queries/candidate_cds.fa" \
        > "$working_root/self_map/candidate_cds_to_genome.paf" \
        2> "$working_root/self_map/minimap2.stderr.log"
"$python_bin" -m ploidypatch.cli evidence audit-natural-candidates \
    --candidate-gff "$candidate_gff" --base-gff "$base_gff" \
    --genome-fasta "$genome" --review-rankings "$rankings" \
    --isoseq-evidence "$isoseq" \
    --self-map-paf "$working_root/self_map/candidate_cds_to_genome.paf" \
    --repeat-gff "$repeat_gff" --repeat-seqid-map "$repeat_seqid_map" \
    --repeat-flank-bp 2000 \
    --review-budget 25 --review-budget 50 --review-budget 100 --review-budget 200 \
    --minimum-query-coverage 0.90 --minimum-identity 0.98 \
    --near-equal-score-fraction 0.95 \
    --output-tsv "$working_root/audit/candidate_biological_audit.tsv" \
    --output-summary "$working_root/audit/summary.json" \
    > "$working_root/audit/stdout.json" \
    2> "$working_root/audit/stderr.log"

for required in "$working_root/queries/candidate_cds.fa" \
                "$working_root/self_map/candidate_cds_to_genome.paf" \
                "$working_root/audit/candidate_biological_audit.tsv" \
                "$working_root/audit/summary.json"; do
    [[ -s $required ]] || { echo "missing biological-audit output: $required" >&2; exit 1; }
done
"$python_bin" - "$working_root/audit/summary.json" <<'PY'
import json, sys
report = json.load(open(sys.argv[1], encoding="utf-8"))
if report["counts"]["candidates"] != 14853:
    raise SystemExit("unexpected maize natural candidate count")
if report["counts"]["ranking_rows"] != 29706:
    raise SystemExit("unexpected maize natural ranking row count")
if report["interpretation"]["automatic_approval"] is not False:
    raise SystemExit("automatic approval unexpectedly enabled")
PY
(
    cd "$working_root"
    find . -type f \( -name '*.json' -o -name '*.tsv' -o -name '*.fa' \
        -o -name '*.paf' \) -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS
    sha256sum -c SHA256SUMS >/dev/null
)
du -sb "$working_root" > "$working_root/disk_bytes.txt"
mv "$working_root" "$result_root"
printf 'maize natural biological audit frozen: %s\n' "$result_root"
