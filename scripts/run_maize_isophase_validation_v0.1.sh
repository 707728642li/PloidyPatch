#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then echo "usage: $0 PROJECT_ROOT" >&2; exit 2; fi
project_root=$(realpath "$1")
code_root=$project_root/code
python_bin=$project_root/envs/ploidypatch-dev/bin/python
minimap2=$project_root/envs/ploidypatch-pav/bin/minimap2
evidence=$project_root/data/validation/maize_isophase_zenodo_2611319
discovery=$project_root/results/natural/maize_v2/discovery
rank_root=$discovery/homeolog_ranker
candidate_gff=$discovery/method_trio/consensus/union/natural/candidate.gff3
rankings=$rank_root/natural/review_rankings.tsv
genome=$project_root/data/derived/holdout_inputs/maize_v2/zea_mays/primary_chromosomes.genome.fa
protocol=$code_root/docs/MAIZE_NATURAL_VALIDATION_PROTOCOL_v0.1.md
result_root=$project_root/results/natural/maize_v2/validation/isoseq_v0.1
working_root=${result_root}.working

for required in "$python_bin" "$minimap2" "$evidence/SHA256SUMS" \
    "$rank_root/SHA256SUMS" "$candidate_gff" "$rankings" "$genome" \
    "$genome.fai" "$protocol"; do
    [[ -s $required ]] || { echo "missing maize Iso-Seq validation input: $required" >&2; exit 1; }
done
(cd "$evidence" && sha256sum -c SHA256SUMS >/dev/null)
(cd "$rank_root" && sha256sum -c SHA256SUMS >/dev/null)
grep -q $'^candidate_freeze_precedes_validation_access\ttrue$' \
    "$rank_root/run_contract.tsv" || {
    echo "candidate/rank evidence firewall is not declared" >&2; exit 1;
}
if [[ -e $result_root || -e $working_root ]]; then
    echo "refusing to overwrite maize Iso-Seq validation" >&2; exit 1
fi
mkdir -p "$working_root"/{prepared,alignment,evidence}
{
    printf 'field\tvalue\n'
    printf 'code_commit\t%s\n' "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'protocol\tmaize_natural_validation_v0.1\n'
    printf 'candidate_and_rank_freeze_verified\ttrue\n'
    printf 'candidate_coordinates_or_ranks_modified\tfalse\n'
    printf 'primary_b73_count_columns\tEM1,R1,END1\n'
    printf 'minimum_b73_full_length_reads\t2\n'
    printf 'alignment_preset\tsplice:hq,-uf\n'
    printf 'secondary_alignments_retained_for_ambiguity_filter\ttrue\n'
    printf 'automatic_annotation_patch\tfalse\n'
    printf 'minimap2_version\t%s\n' "$($minimap2 --version)"
    printf 'protocol_sha256\t%s\n' "$(sha256sum "$protocol" | awk '{print $1}')"
} > "$working_root/run_contract.tsv"
{
    printf 'role\tbytes\tsha256\tpath\n'
    for entry in \
        "candidate_gff:$candidate_gff" "review_rankings:$rankings" \
        "candidate_rank_freeze:$rank_root/SHA256SUMS" \
        "validation_freeze:$evidence/SHA256SUMS" \
        "transcript_fasta:$evidence/files/F1maize.FINAL.fasta" \
        "transcript_counts:$evidence/files/F1maize.FINAL.demux_FL_count.txt" \
        "target_genome:$genome" "target_genome_fai:$genome.fai"; do
        role=${entry%%:*}; path=${entry#*:}
        printf '%s\t%s\t%s\t%s\n' "$role" "$(stat -Lc %s "$path")" \
            "$(sha256sum "$path" | awk '{print $1}')" "$path"
    done
} > "$working_root/input_manifest.tsv"

cd "$code_root"
"$python_bin" -m ploidypatch.cli evidence prepare-b73-isoseq \
    --fasta "$evidence/files/F1maize.FINAL.fasta" \
    --counts-csv "$evidence/files/F1maize.FINAL.demux_FL_count.txt" \
    --minimum-b73-reads 2 \
    --output-fasta "$working_root/prepared/b73_observed.fasta" \
    --output-counts "$working_root/prepared/b73_observed_counts.tsv" \
    > "$working_root/prepared/stdout.json" \
    2> "$working_root/prepared/stderr.log"

{
    printf '%q ' "$minimap2" -t 64 -x splice:hq -uf --secondary=yes -N 10 \
        -c --cs=long "$genome" "$working_root/prepared/b73_observed.fasta"
    printf '\n'
} > "$working_root/alignment/command.txt"
/usr/bin/time -v -o "$working_root/alignment/resource.time.txt" \
    "$minimap2" -t 64 -x splice:hq -uf --secondary=yes -N 10 -c --cs=long \
    "$genome" "$working_root/prepared/b73_observed.fasta" \
    > "$working_root/alignment/b73_observed.nam5.paf" \
    2> "$working_root/alignment/minimap2.stderr.log"
[[ -s $working_root/alignment/b73_observed.nam5.paf ]] || {
    echo "minimap2 produced an empty PAF" >&2; exit 1;
}

"$python_bin" -m ploidypatch.cli evidence validate-isoseq-candidates \
    --candidate-gff "$candidate_gff" \
    --paf "$working_root/alignment/b73_observed.nam5.paf" \
    --selected-counts "$working_root/prepared/b73_observed_counts.tsv" \
    --genome-fasta "$genome" \
    --minimum-query-coverage 0.90 --minimum-identity 0.98 \
    --minimum-mapq 20 --maximum-secondary-score-fraction 0.95 \
    --minimum-candidate-cds-coverage 0.90 --flank-bp 5000 \
    --output-evidence "$working_root/evidence/candidate_isoseq_evidence.tsv" \
    > "$working_root/evidence/validation.stdout.json" \
    2> "$working_root/evidence/validation.stderr.log"
"$python_bin" -m ploidypatch.cli evidence join-isoseq-review-rankings \
    --evidence "$working_root/evidence/candidate_isoseq_evidence.tsv" \
    --review-rankings "$rankings" \
    --review-budget 25 --review-budget 50 --review-budget 100 --review-budget 200 \
    --output-tsv "$working_root/evidence/ranked_isoseq_evidence.tsv" \
    --output-summary "$working_root/evidence/review_yield.json" \
    > "$working_root/evidence/review_join.stdout.json" \
    2> "$working_root/evidence/review_join.stderr.log"

for required in \
    "$working_root/prepared/b73_observed.fasta.manifest.json" \
    "$working_root/evidence/candidate_isoseq_evidence.tsv.manifest.json" \
    "$working_root/evidence/review_yield.json"; do
    [[ -s $required ]] || { echo "missing validation output: $required" >&2; exit 1; }
done
(
    cd "$working_root"
    find . -type f ! -name SHA256SUMS -print0 | sort -z | xargs -0 sha256sum \
        > SHA256SUMS
    sha256sum -c SHA256SUMS >/dev/null
)
du -sb "$working_root" > "$working_root/disk_bytes.txt"
mv "$working_root" "$result_root"
printf 'maize Iso-Seq validation frozen: %s\n' "$result_root"
