#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then echo "usage: $0 PROJECT_ROOT" >&2; exit 2; fi
project_root=$(realpath "$1")
code_root=$project_root/code
dev_python=$project_root/envs/ploidypatch-dev/bin/python
pychopper_env=$project_root/envs/ploidypatch-pychopper
pychopper_python=$pychopper_env/bin/python
minimap2=$project_root/envs/ploidypatch-pav/bin/minimap2
conda_exe=/data/codexli/software/conda/miniforge3/bin/conda
filter_script=$code_root/scripts/filter_fastq_length_quality.py
pychopper_wrapper=$code_root/scripts/run_pychopper_seeded.py
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
raw_alignment_root=$validation_parent/ont_raw_v0.1_failed_csv_field_limit
raw_patch2_root=$validation_parent/ont_raw_v0.1_patch2_transcript_strand
source_index=$raw_alignment_root/alignment/gddh13_primary.splice_k14.mmi
result_root=$validation_parent/ont_flnc_pychopper_v0.1
working_root=${result_root}.working
code_commit=${PLOIDYPATCH_CODE_COMMIT:-}
autotune_seed=20261005

[[ $code_commit =~ ^[0-9a-f]{40}$ ]] || {
    echo "PLOIDYPATCH_CODE_COMMIT must be a full git SHA" >&2; exit 1;
}
verify_tree() { (cd "$1" && sha256sum -c SHA256SUMS >/dev/null); }
for required in "$dev_python" "$pychopper_python" \
    "$pychopper_env/bin/pychopper" "$minimap2" "$conda_exe" \
    "$filter_script" "$pychopper_wrapper" "$genome" "$genome.fai" \
    "$base_gff" "$candidate_gff" "$rankings" "$repeat_gff" \
    "$ont_root/metadata/file_contract.tsv" "$source_index" \
    "$raw_alignment_root/SHA256SUMS" "$raw_patch2_root/SHA256SUMS" \
    "$raw_patch2_root/self_map/candidate_cds_to_genome.paf" \
    "$rank_root/SHA256SUMS" "$ont_root/SHA256SUMS" \
    "$te_root/SHA256SUMS" "$method_root/SHA256SUMS"; do
    [[ -s $required ]] || { echo "missing FLNC validation input: $required" >&2; exit 1; }
done
[[ ! -e $result_root && ! -e $working_root ]] || {
    echo "refusing to overwrite apple Pychopper FLNC validation" >&2; exit 1;
}
for root in "$raw_alignment_root" "$raw_patch2_root" "$rank_root" \
    "$ont_root" "$te_root" "$method_root"; do
    verify_tree "$root"
done
grep -q $'^candidate_and_rank_freeze_precedes_validation_access\ttrue$' \
    "$rank_root/run_contract.tsv" || {
    echo "candidate/RNA evidence firewall is absent" >&2; exit 1;
}

pychopper_version=$(
    "$pychopper_python" -c \
        'from importlib.metadata import version; print(version("pychopper"))'
)
[[ $pychopper_version == 2.7.10 ]] || {
    echo "expected Pychopper 2.7.10, observed $pychopper_version" >&2; exit 1;
}
pandas_version=$(
    "$pychopper_python" -c \
        'from importlib.metadata import version; import pandas, pytz; print(version("pandas"))'
)
[[ ${pandas_version%%.*} == 2 ]] || {
    echo "Pychopper 2.7.10 reporting requires pandas 2.x; observed $pandas_version" >&2
    exit 1
}
pytz_version=$(
    "$pychopper_python" -c \
        'from importlib.metadata import version; print(version("pytz"))'
)
tqdm_version=$(
    "$pychopper_python" -c \
        'from importlib.metadata import version; print(version("tqdm"))'
)
mkdir -p "$working_root"/{preprocessing/by_run,flnc/by_run,alignment/by_run,alignment/combined,evidence,audit,self_map,logs,freeze}
"$conda_exe" list -p "$pychopper_env" --explicit \
    > "$working_root/freeze/pychopper_conda_explicit.txt"
"$conda_exe" list -p "$pychopper_env" --json \
    > "$working_root/freeze/pychopper_conda_packages.json"
if "$pychopper_python" -m pip --version >/dev/null 2>&1; then
    "$pychopper_python" -m pip freeze --all \
        > "$working_root/freeze/pychopper_pip_freeze.txt"
else
    printf '# pip is absent; the conda explicit lock is authoritative.\n' \
        > "$working_root/freeze/pychopper_pip_freeze.txt"
fi

{
    printf 'field\tvalue\n'
    printf 'code_commit\t%s\n' "$code_commit"
    printf 'target\tMalus_domestica_GDDH13_v1.1\n'
    printf 'validation_source\tGolden_Delicious_GSA_CRA021523_raw_ONT\n'
    printf 'source_protocol\tSQK-PCS109_cDNA-PCR\n'
    printf 'source_publication_doi\t10.3389/fpls.2026.1819201\n'
    printf 'analysis_role\tdescriptive_posthoc_natural_validation\n'
    printf 'labels_seen\ttrue\n'
    printf 'candidate_coordinates_modified\tfalse\n'
    printf 'candidate_ranks_modified\tfalse\n'
    printf 'evidence_thresholds_modified\tfalse\n'
    printf 'input_read_filter\tlength_ge_500_before_Pychopper\n'
    printf 'mean_quality_filter\tPychopper_mean_error_probability_PHRED_ge_7\n'
    printf 'pychopper_version\t%s\n' "$pychopper_version"
    printf 'pandas_version\t%s\n' "$pandas_version"
    printf 'pandas_compatibility\tPychopper_2.7.10_reporting_requires_pandas_2.x\n'
    printf 'pytz_version\t%s\n' "$pytz_version"
    printf 'tqdm_version\t%s\n' "$tqdm_version"
    printf 'pychopper_upstream_tqdm_metadata_pin\t4.26.0\n'
    printf 'tqdm_runtime_policy\tBioconda_curated_modern_runtime_metadata_exception\n'
    printf 'pychopper_kit\tPCS109\n'
    printf 'pychopper_method\tphmm\n'
    printf 'pychopper_cutoff\tautotuned\n'
    printf 'pychopper_autotune_numpy_seed\t%s\n' "$autotune_seed"
    printf 'source_date_epoch\t0\n'
    printf 'pychopper_minimum_mean_quality\t7\n'
    printf 'pychopper_minimum_output_segment_length\t50\n'
    printf 'pychopper_primary_full_length_only\ttrue\n'
    printf 'pychopper_rescued_segments_used\tfalse\n'
    printf 'threads_per_tissue\t16\nparallel_tissues\t7\n'
    printf 'alignment_preset\tsplice-uf-k14-G1000000\n'
    printf 'secondary_alignments\t10\n'
    printf 'frozen_minimap2_index_reused\ttrue\n'
    printf 'alignment_strand_source\tquery_orientation\n'
    printf 'strand_rationale\tPychopper_orients_FLNC_reads_5prime_to_3prime_before_alignment\n'
    printf 'minimum_query_coverage\t0.85\nminimum_identity\t0.90\n'
    printf 'minimum_mapq\t20\nmaximum_secondary_score_fraction\t0.95\n'
    printf 'minimum_candidate_cds_coverage\t0.90\nflank_bp\t5000\n'
    printf 'minimum_case_study_full_chain_reads\t2\n'
    printf 'single_exon_is_not_full_chain_positive\ttrue\n'
    printf 'review_budgets\t100,146,250,292,500,583\n'
    printf 'bootstrap_replicates\t20000\nbootstrap_seed\t20261004\n'
    printf 'repeat_annotation\tGDDH13_v1.1_official_TE_GFF\n'
    printf 'automatic_annotation_patch\tfalse\n'
    printf 'raw_ts_patch2_comparator\t%s\n' "$raw_patch2_root"
    printf 'minimap2_version\t%s\n' "$("$minimap2" --version)"
} > "$working_root/run_contract.tsv"
{
    printf 'role\tbytes\tsha256\tpath\n'
    for entry in "genome:$genome" "genome_fai:$genome.fai" \
        "base_gff:$base_gff" "candidate_gff:$candidate_gff" \
        "review_rankings:$rankings" "repeat_gff:$repeat_gff" \
        "minimap2_index:$source_index" \
        "raw_alignment_freeze:$raw_alignment_root/SHA256SUMS" \
        "raw_ts_patch2_freeze:$raw_patch2_root/SHA256SUMS" \
        "reused_self_map:$raw_patch2_root/self_map/candidate_cds_to_genome.paf" \
        "rank_freeze:$rank_root/SHA256SUMS" \
        "ont_freeze:$ont_root/SHA256SUMS" "te_freeze:$te_root/SHA256SUMS"; do
        role=${entry%%:*}; path=${entry#*:}
        printf '%s\t%s\t%s\t%s\n' "$role" "$(stat -Lc %s "$path")" \
            "$(sha256sum "$path" | awk '{print $1}')" "$path"
    done
} > "$working_root/input_manifest.tsv"
{
    printf 'path\tsha256\n'
    for relative in scripts/run_apple_golden_delicious_ont_flnc_validation_v0.1.sh \
        scripts/filter_fastq_length_quality.py scripts/run_pychopper_seeded.py \
        src/ploidypatch/isoseq_validation.py src/ploidypatch/natural_audit.py \
        src/ploidypatch/cli.py; do
        printf '%s\t%s\n' "$relative" \
            "$(sha256sum "$code_root/$relative" | awk '{print $1}')"
    done
    printf '%s\t%s\n' "envs/ploidypatch-pychopper/bin/pychopper" \
        "$(sha256sum "$pychopper_env/bin/pychopper" | awk '{print $1}')"
} > "$working_root/freeze/code_manifest.tsv"

# Remove reads shorter than 500 nt before Pychopper.  Pychopper's own -Q 7
# implementation then applies the paper's mean-error-probability quality gate;
# keeping that computation in Pychopper avoids an expensive duplicate pass.
# The length-filtered copy is materialized because default cutoff autotuning
# makes multiple passes over the input.
pids=()
while IFS=$'\t' read -r accession tissue expected_md5 url; do
    [[ $accession == accession ]] && continue
    input=$ont_root/files/$accession.fastq.gz
    filtered=$working_root/preprocessing/by_run/$accession.len500.fastq.gz
    summary=$working_root/preprocessing/by_run/$accession.filter_summary.json
    [[ -s $input ]] || { echo "missing ONT FASTQ: $input" >&2; exit 1; }
    (
        set -o pipefail
        /usr/bin/time -v -o "$working_root/logs/$accession.prefilter.time.txt" \
            "$dev_python" "$filter_script" "$input" \
            --minimum-length 500 \
            --summary-json "$summary" | gzip -1 > "$filtered"
        [[ -s $filtered && -s $summary ]] || exit 1
    ) > "$working_root/logs/$accession.prefilter.stdout.log" \
      2> "$working_root/logs/$accession.prefilter.stderr.log" &
    pids+=("$!")
done < "$ont_root/metadata/file_contract.tsv"
failed=0
for pid in "${pids[@]}"; do wait "$pid" || failed=1; done
[[ $failed == 0 ]] || { echo "one or more ONT preprocessing jobs failed" >&2; exit 1; }

# Pychopper's upstream CLI uses NumPy sampling for cutoff autotuning without a
# seed argument.  The frozen adapter supplies a literal seed and otherwise
# invokes the unmodified Pychopper 2.7.10 main routine.
pids=()
while IFS=$'\t' read -r accession tissue expected_md5 url; do
    [[ $accession == accession ]] && continue
    filtered=$working_root/preprocessing/by_run/$accession.len500.fastq.gz
    flnc=$working_root/flnc/by_run/$accession.flnc.fastq.gz
    (
        set -o pipefail
        mkdir -p "$working_root/flnc/by_run/$accession.mplconfig"
        MPLCONFIGDIR="$working_root/flnc/by_run/$accession.mplconfig" \
        OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONHASHSEED=0 \
        SOURCE_DATE_EPOCH=0 \
        PATH="$pychopper_env/bin:/usr/bin:/bin" \
        /usr/bin/time -v -o "$working_root/logs/$accession.pychopper.time.txt" \
            "$pychopper_python" "$pychopper_wrapper" "$autotune_seed" \
            -k PCS109 -m phmm -Q 7 -z 50 -t 16 \
            -r "$working_root/flnc/by_run/$accession.report.pdf" \
            -S "$working_root/flnc/by_run/$accession.stats.tsv" \
            "$filtered" - | gzip -1 > "$flnc"
        [[ -s $flnc ]] || exit 1
        [[ -s $working_root/flnc/by_run/$accession.stats.tsv ]] || exit 1
        [[ -s $working_root/flnc/by_run/$accession.report.pdf ]] || exit 1
    ) > "$working_root/logs/$accession.pychopper.stdout.log" \
      2> "$working_root/logs/$accession.pychopper.stderr.log" &
    pids+=("$!")
done < "$ont_root/metadata/file_contract.tsv"
failed=0
for pid in "${pids[@]}"; do wait "$pid" || failed=1; done
[[ $failed == 0 ]] || { echo "one or more Pychopper jobs failed" >&2; exit 1; }

"$dev_python" - "$ont_root/metadata/file_contract.tsv" \
    "$working_root/preprocessing/by_run" "$working_root/flnc/by_run" \
    "$working_root/flnc/classification_summary.tsv" \
    "$working_root/flnc/classification_summary.json" <<'PY'
import csv
import gzip
import json
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

contract, filter_root, stats_root, output_tsv, output_json = sys.argv[1:]

def count_fastq_records(path):
    count = 0
    with gzip.open(path, "rt", encoding="ascii", newline="") as handle:
        while True:
            name = handle.readline()
            if not name:
                break
            sequence = handle.readline()
            plus = handle.readline()
            quality = handle.readline()
            if not sequence or not plus or not quality:
                raise SystemExit(f"{path}: truncated FASTQ record {count + 1}")
            sequence = sequence.rstrip("\r\n")
            quality = quality.rstrip("\r\n")
            if not name.startswith("@") or not plus.startswith("+"):
                raise SystemExit(f"{path}: malformed FASTQ record {count + 1}")
            if len(sequence) != len(quality):
                raise SystemExit(f"{path}: sequence/quality length mismatch at record {count + 1}")
            count += 1
    return count

with open(contract, encoding="utf-8", newline="") as handle:
    contract_rows = list(csv.DictReader(handle, delimiter="\t"))
with ThreadPoolExecutor(max_workers=7) as executor:
    output_counts = dict(zip(
        (row["accession"] for row in contract_rows),
        executor.map(
            count_fastq_records,
            (Path(stats_root) / f'{row["accession"]}.flnc.fastq.gz' for row in contract_rows),
        ),
        strict=True,
    ))
rows = []
for row in contract_rows:
        accession = row["accession"]
        qc = json.load(open(Path(filter_root) / f"{accession}.filter_summary.json", encoding="utf-8"))
        values = {}
        with open(Path(stats_root) / f"{accession}.stats.tsv", encoding="utf-8", newline="") as stats_handle:
            for stat in csv.DictReader(stats_handle, delimiter="\t"):
                values[(stat["Category"], stat["Name"])] = int(float(stat["Value"]))
        pass_reads = values[("ReadStats", "PassReads")]
        qc_fail = values[("ReadStats", "QcFail")]
        primary_classified = values[("Classification", "Primers_found")]
        primary_flnc = output_counts[accession]
        output_length_failed = values[("ReadStats", "LenFail")]
        rescued_segments = values[("Classification", "Rescue")]
        rescued_reads = sum(
            value
            for (category, _), value in values.items()
            if category == "RescueSegmentNr"
        )
        unusable = values[("Classification", "Unusable")]
        if pass_reads + qc_fail != qc["passed_records"]:
            raise SystemExit(f"{accession}: preprocessing and Pychopper Q7 counts disagree")
        if primary_classified + rescued_reads + unusable != pass_reads:
            raise SystemExit(f"{accession}: Pychopper read classes do not conserve Q7-pass reads")
        if not 0 < primary_flnc <= primary_classified <= pass_reads:
            raise SystemExit(f"{accession}: invalid primary FLNC count")
        if primary_classified - primary_flnc > output_length_failed:
            raise SystemExit(f"{accession}: primary-classified/output discrepancy is unexplained")
        record = {
            "accession": accession,
            "tissue": row["tissue"],
            "input_reads": qc["input_records"],
            "prefilter_pass_reads": qc["passed_records"],
            "prefilter_length_failed": qc["length_failed_records"],
            "pychopper_quality_failed": qc_fail,
            "pychopper_primary_classified": primary_classified,
            "pychopper_primary_flnc": primary_flnc,
            "pychopper_output_segments_below_50": output_length_failed,
            "pychopper_rescued_reads_excluded": rescued_reads,
            "pychopper_rescued_segments_excluded": rescued_segments,
            "pychopper_unusable_reads": unusable,
            "primary_flnc_fraction_of_q7_pass": primary_flnc / pass_reads,
            "within_publication_reported_count_range": 3_587_082 <= primary_flnc <= 4_693_774,
            "within_publication_reported_fraction_range": 0.8584 <= primary_flnc / pass_reads <= 0.8786,
        }
        rows.append(record)
fields = list(rows[0])
with open(output_tsv, "x", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
report = {
    "schema_version": "ploidypatch.apple_pychopper_flnc_summary.v1",
    "paper_reported_full_length_reads_per_tissue": [3_587_082, 4_693_774],
    "paper_reported_full_length_fraction_per_tissue": [0.8584, 0.8786],
    "reported_ranges_are_descriptive_not_validity_gates": True,
    "tissues": rows,
    "totals": {
        key: sum(row[key] for row in rows)
        for key in (
            "input_reads", "prefilter_pass_reads", "prefilter_length_failed",
            "pychopper_quality_failed", "pychopper_primary_classified",
            "pychopper_primary_flnc", "pychopper_output_segments_below_50",
            "pychopper_rescued_reads_excluded",
            "pychopper_rescued_segments_excluded", "pychopper_unusable_reads",
        )
    },
}
with open(output_json, "x", encoding="utf-8", newline="") as handle:
    json.dump(report, handle, indent=2, sort_keys=True)
    handle.write("\n")
PY

# Align the oriented primary FLNC reads exactly with the publication's
# splice/-uf/k14/1-Mb intron semantics, retaining secondary alignments so the
# validator can reject ambiguous loci.
pids=()
while IFS=$'\t' read -r accession tissue expected_md5 url; do
    [[ $accession == accession ]] && continue
    flnc=$working_root/flnc/by_run/$accession.flnc.fastq.gz
    output=$working_root/alignment/by_run/$accession.paf
    (
        /usr/bin/time -v -o "$working_root/logs/$accession.minimap2.time.txt" \
            "$minimap2" -t 16 -x splice -uf -k14 -G 1000000 \
            --secondary=yes -N 10 -c --cs=long "$source_index" "$flnc" \
            > "$output" 2> "$working_root/logs/$accession.minimap2.stderr.log"
        [[ -s $output ]] || exit 1
    ) &
    pids+=("$!")
done < "$ont_root/metadata/file_contract.tsv"
failed=0
for pid in "${pids[@]}"; do wait "$pid" || failed=1; done
[[ $failed == 0 ]] || { echo "one or more FLNC alignments failed" >&2; exit 1; }

cd "$code_root"
combined=$working_root/alignment/combined/golden_delicious_flnc.candidate_query_universe.paf
filter_args=()
while IFS=$'\t' read -r accession _; do
    [[ $accession == accession ]] && continue
    filter_args+=(--paf-input "$accession=$working_root/alignment/by_run/$accession.paf")
done < "$ont_root/metadata/file_contract.tsv"
/usr/bin/time -v -o "$working_root/logs/query_filter.time.txt" \
    "$dev_python" -m ploidypatch.cli evidence filter-candidate-query-paf \
    --candidate-gff "$candidate_gff" "${filter_args[@]}" \
    --alignment-strand-source query_orientation \
    --output-paf "$combined" \
    --output-counts "$working_root/evidence/flnc_read_counts.tsv" \
    --output-summary "$working_root/alignment/combined/query_filter_summary.tsv" \
    --output-manifest "$working_root/alignment/combined/query_filter_manifest.json" \
    > "$working_root/alignment/combined/query_filter.stdout.json" \
    2> "$working_root/alignment/combined/query_filter.stderr.log"
/usr/bin/time -v -o "$working_root/logs/validation.time.txt" \
    "$dev_python" -m ploidypatch.cli evidence validate-isoseq-candidates \
    --candidate-gff "$candidate_gff" --paf "$combined" \
    --selected-counts "$working_root/evidence/flnc_read_counts.tsv" \
    --genome-fasta "$genome" --minimum-query-coverage 0.85 \
    --minimum-identity 0.90 --minimum-mapq 20 \
    --maximum-secondary-score-fraction 0.95 \
    --minimum-candidate-cds-coverage 0.90 --flank-bp 5000 \
    --alignment-strand-source query_orientation \
    --output-evidence "$working_root/evidence/candidate_flnc_evidence.tsv" \
    > "$working_root/evidence/validation.stdout.json" \
    2> "$working_root/evidence/validation.stderr.log"
"$dev_python" -m ploidypatch.cli evidence join-isoseq-review-rankings \
    --evidence "$working_root/evidence/candidate_flnc_evidence.tsv" \
    --review-rankings "$rankings" \
    --review-budget 100 --review-budget 146 --review-budget 250 \
    --review-budget 292 --review-budget 500 --review-budget 583 \
    --comparator-estimator baseline --primary-estimator v04_guard \
    --output-tsv "$working_root/evidence/ranked_flnc_evidence.tsv" \
    --output-summary "$working_root/evidence/review_yield.json" \
    > "$working_root/evidence/review_join.stdout.json" \
    2> "$working_root/evidence/review_join.stderr.log"
"$dev_python" -m ploidypatch.cli evidence bootstrap-isoseq-review-yield \
    --evidence "$working_root/evidence/candidate_flnc_evidence.tsv" \
    --review-rankings "$rankings" \
    --review-budget 100 --review-budget 146 --review-budget 250 \
    --review-budget 292 --review-budget 500 --review-budget 583 \
    --comparator-estimator baseline --primary-estimator v04_guard \
    --replicates 20000 --seed 20261004 \
    --output-json "$working_root/evidence/bootstrap.json" \
    > "$working_root/evidence/bootstrap.stdout.json" \
    2> "$working_root/evidence/bootstrap.stderr.log"

cp --reflink=auto "$raw_patch2_root/self_map/candidate_cds_to_genome.paf" \
    "$working_root/self_map/candidate_cds_to_genome.paf"
[[ $(sha256sum "$working_root/self_map/candidate_cds_to_genome.paf" | awk '{print $1}') == \
   $(sha256sum "$raw_patch2_root/self_map/candidate_cds_to_genome.paf" | awk '{print $1}') ]] || {
    echo "reused candidate self-map changed" >&2; exit 1;
}
"$dev_python" -m ploidypatch.cli evidence audit-natural-candidates \
    --candidate-gff "$candidate_gff" --base-gff "$base_gff" \
    --genome-fasta "$genome" --review-rankings "$rankings" \
    --isoseq-evidence "$working_root/evidence/candidate_flnc_evidence.tsv" \
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

"$dev_python" - "$raw_patch2_root/evidence/candidate_ont_evidence.tsv" \
    "$working_root/evidence/candidate_flnc_evidence.tsv" \
    "$raw_patch2_root/evidence/review_yield.json" \
    "$working_root/evidence/review_yield.json" \
    "$raw_patch2_root/audit/summary.json" "$working_root/audit/summary.json" \
    "$working_root/evidence/raw_ts_patch2_vs_flnc.json" <<'PY'
import csv
import hashlib
import json
import sys
from collections import Counter

csv.field_size_limit(2**31 - 1)
raw_path, flnc_path, raw_review_path, flnc_review_path, raw_audit_path, flnc_audit_path, output = sys.argv[1:]

def sha(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()

def read_evidence(path):
    with open(path, encoding="utf-8", newline="") as handle:
        return {row["candidate_digest"]: row for row in csv.DictReader(handle, delimiter="\t")}

raw = read_evidence(raw_path)
flnc = read_evidence(flnc_path)
if set(raw) != set(flnc):
    raise SystemExit("FLNC candidate universe differs from raw-ts patch2")
raw_full = {d for d, row in raw.items() if row["evidence_state"] == "full_chain_supported"}
flnc_full = {d for d, row in flnc.items() if row["evidence_state"] == "full_chain_supported"}
raw_reads2 = {d for d in raw_full if int(raw[d]["supporting_full_length_reads"]) >= 2}
flnc_reads2 = {d for d in flnc_full if int(flnc[d]["supporting_full_length_reads"]) >= 2}
transitions = Counter((raw[d]["evidence_state"], flnc[d]["evidence_state"]) for d in raw)
raw_audit = json.load(open(raw_audit_path, encoding="utf-8"))
flnc_audit = json.load(open(flnc_audit_path, encoding="utf-8"))
report = {
    "schema_version": "ploidypatch.apple_raw_ts_vs_pychopper_flnc.v1",
    "scope": "descriptive_labels_seen_candidates_and_ranks_unchanged",
    "inputs": {
        "raw_ts_evidence_sha256": sha(raw_path),
        "flnc_evidence_sha256": sha(flnc_path),
        "raw_ts_review_sha256": sha(raw_review_path),
        "flnc_review_sha256": sha(flnc_review_path),
        "raw_ts_audit_sha256": sha(raw_audit_path),
        "flnc_audit_sha256": sha(flnc_audit_path),
    },
    "counts": {
        "candidate_universe": len(raw),
        "raw_ts_full_chain": len(raw_full),
        "pychopper_flnc_full_chain": len(flnc_full),
        "full_chain_intersection": len(raw_full & flnc_full),
        "raw_ts_full_chain_reads_ge_2": len(raw_reads2),
        "pychopper_flnc_full_chain_reads_ge_2": len(flnc_reads2),
        "reads_ge_2_intersection": len(raw_reads2 & flnc_reads2),
        "raw_ts_case_study_ready": raw_audit["counts"]["case_study_ready"],
        "pychopper_flnc_case_study_ready": flnc_audit["counts"]["case_study_ready"],
    },
    "state_transitions": {
        f"{left}->{right}": count
        for (left, right), count in sorted(transitions.items())
    },
    "review_primary": {
        "raw_ts": json.load(open(raw_review_path, encoding="utf-8"))["primary"],
        "pychopper_flnc": json.load(open(flnc_review_path, encoding="utf-8"))["primary"],
    },
}
with open(output, "x", encoding="utf-8", newline="") as handle:
    json.dump(report, handle, indent=2, sort_keys=True)
    handle.write("\n")
PY

"$dev_python" - "$working_root/flnc/classification_summary.json" \
    "$working_root/alignment/combined/query_filter_manifest.json" \
    "$working_root/evidence/candidate_flnc_evidence.tsv.manifest.json" \
    "$working_root/audit/summary.json" <<'PY'
import json
import sys

classification_path, filter_path, evidence_path, audit_path = sys.argv[1:]
classification = json.load(open(classification_path, encoding="utf-8"))
query_filter = json.load(open(filter_path, encoding="utf-8"))
if len(classification["tissues"]) != 7 or len(query_filter["inputs"]["paf"]) != 7:
    raise SystemExit("FLNC validation did not contain exactly seven tissues")
if query_filter["parameters"]["alignment_strand_source"] != "query_orientation":
    raise SystemExit("oriented FLNC filter did not use PAF query orientation")
if query_filter["counts"]["strand_unavailable_alignments"] != 0:
    raise SystemExit("oriented FLNC alignments unexpectedly lacked strand")
evidence = json.load(open(evidence_path, encoding="utf-8"))
if evidence["schema_version"] != "ploidypatch.isoseq_candidate_validation.v2":
    raise SystemExit("FLNC evidence schema is not v2")
if evidence["parameters"]["alignment_strand_source"] != "query_orientation":
    raise SystemExit("oriented FLNC validator used the wrong strand source")
if evidence["counts"]["candidate_models"] != 29144:
    raise SystemExit("FLNC candidate universe changed")
if evidence["counts"]["selected_transcripts"] != query_filter["counts"]["retained_queries"]:
    raise SystemExit("FLNC filter/validator read universes differ")
audit = json.load(open(audit_path, encoding="utf-8"))
if audit["counts"]["candidates"] != 29144:
    raise SystemExit("FLNC audit candidate universe changed")
if audit["parameters"]["minimum_full_length_read_support"] != 2:
    raise SystemExit("FLNC case threshold changed")
PY

for required in \
    "$working_root/flnc/classification_summary.json" \
    "$working_root/alignment/combined/query_filter_manifest.json" \
    "$working_root/evidence/candidate_flnc_evidence.tsv" \
    "$working_root/evidence/candidate_flnc_evidence.tsv.manifest.json" \
    "$working_root/evidence/review_yield.json" \
    "$working_root/evidence/bootstrap.json" \
    "$working_root/evidence/raw_ts_patch2_vs_flnc.json" \
    "$working_root/self_map/candidate_cds_to_genome.paf" \
    "$working_root/audit/candidate_biological_audit.tsv" \
    "$working_root/audit/summary.json"; do
    [[ -s $required ]] || { echo "missing FLNC validation output: $required" >&2; exit 1; }
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
printf 'apple Golden Delicious Pychopper FLNC validation frozen: %s\n' "$result_root"
