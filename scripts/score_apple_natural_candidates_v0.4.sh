#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then echo "usage: $0 PROJECT_ROOT" >&2; exit 2; fi
project_root=$(realpath "$1")
code_root=$project_root/code
python_bin=$project_root/envs/ploidypatch-dev/bin/python
bundle=$project_root/data/derived/external_inputs/apple_v0.3/target_apple
base_gff=$bundle/primary_chromosomes.gff3
method_root=$project_root/results/copy_collapse/external/apple_v0.3_method_trio
candidate_root=$method_root/consensus/primary_union/complete_control
candidate_gff=$candidate_root/candidate.gff3
pool_decisions=$candidate_root/decisions.tsv
pool_manifest=$candidate_root/candidate.gff3.manifest.json
self_wgd_root=$project_root/results/natural/apple_gddh13_v0.4/discovery/self_wgd
prior_wgd=$self_wgd_root/selected/selection.tsv
composite_root=$project_root/results/models/ploidypatch_ranker_v0.4
model_v03=$composite_root/model_v0.3.json
guard_policy=$composite_root/guard_policy.json
result_root=$project_root/results/natural/apple_gddh13_v0.4/discovery/rankings
working_root=${result_root}.working

verify_tree() { (cd "$1" && sha256sum -c SHA256SUMS >/dev/null); }
for required in "$python_bin" "$base_gff" "$method_root/SHA256SUMS" \
    "$self_wgd_root/SHA256SUMS" "$composite_root/SHA256SUMS" \
    "$candidate_gff" "$pool_decisions" "$pool_manifest" "$prior_wgd" \
    "$model_v03" "$guard_policy" "$composite_root/composite_manifest.json" \
    "$method_root/methods/miniprot/complete_control/decisions.tsv" \
    "$method_root/methods/gemoma/complete_control/decisions.tsv" \
    "$method_root/methods/lifton/complete_control/decisions.tsv"; do
    [[ -s $required ]] || { echo "missing apple natural ranking input: $required" >&2; exit 1; }
done
verify_tree "$method_root"
verify_tree "$self_wgd_root"
verify_tree "$composite_root"
grep -q $'^RNA_access\tfalse$' "$self_wgd_root/run_contract.tsv" || {
    echo "self-WGD discovery did not declare RNA blindness" >&2; exit 1;
}
[[ ! -e $result_root && ! -e $working_root ]] || {
    echo "refusing to overwrite apple natural rankings" >&2; exit 1;
}
mkdir -p "$working_root"/{features,scores,freeze,natural}
{
    printf 'field\tvalue\ncode_commit\t%s\n' \
        "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'split\tnatural_current_annotation_discovery_v0.4\n'
    printf 'target\tMalus_domestica_GDDH13_v1.1\n'
    printf 'truth_access\tfalse\nRNA_access\tfalse\nvalidation_sequence_access\tfalse\n'
    printf 'candidate_and_rank_freeze_precedes_validation_access\ttrue\n'
    printf 'candidate_policy\tretain_distinct_phased_CDS_chains\n'
    printf 'rank_order\tdescending_score_then_candidate_digest\n'
    printf 'interpretation\treview_priority_only\ncalibration\tnone\n'
    printf 'automatic_approval\tfalse\n'
    printf 'v03_model_sha256\t%s\n' "$(sha256sum "$model_v03" | awk '{print $1}')"
    printf 'v04_guard_policy_sha256\t%s\n' \
        "$(sha256sum "$guard_policy" | awk '{print $1}')"
    printf 'composite_sha256sums_sha256\t%s\n' \
        "$(sha256sum "$composite_root/SHA256SUMS" | awk '{print $1}')"
} > "$working_root/run_contract.tsv"
{
    printf 'role\tbytes\tsha256\tpath\n'
    for entry in "base_gff:$base_gff" "candidate_gff:$candidate_gff" \
        "pool_decisions:$pool_decisions" "pool_manifest:$pool_manifest" \
        "prior_wgd:$prior_wgd" "model_v03:$model_v03" \
        "guard_policy:$guard_policy" \
        "self_wgd_freeze:$self_wgd_root/SHA256SUMS" \
        "method_trio_freeze:$method_root/SHA256SUMS"; do
        role=${entry%%:*}; path=${entry#*:}
        printf '%s\t%s\t%s\t%s\n' "$role" "$(stat -Lc %s "$path")" \
            "$(sha256sum "$path" | awk '{print $1}')" "$path"
    done
} > "$working_root/input_manifest.tsv"
{
    printf 'path\tsha256\n'
    for relative in scripts/score_apple_natural_candidates_v0.4.sh \
        src/ploidypatch/copy_features.py src/ploidypatch/homeolog_topology.py \
        src/ploidypatch/wgd_candidate_select.py src/ploidypatch/support_ranker.py \
        src/ploidypatch/conflict_guard.py src/ploidypatch/cli.py; do
        printf '%s\t%s\n' "$relative" \
            "$(sha256sum "$code_root/$relative" | awk '{print $1}')"
    done
} > "$working_root/freeze/code_manifest.tsv"

cd "$code_root"
"$python_bin" -m ploidypatch.cli evidence propagate-wgd-conflict-partners \
    --base-gff "$base_gff" --candidate-gff "$candidate_gff" \
    --pool-decisions "$pool_decisions" --prior-wgd-selection "$prior_wgd" \
    --output-selection "$working_root/features/wgd_selection.tsv" \
    > "$working_root/features/wgd_selection.stdout.json" \
    2> "$working_root/features/wgd_selection.stderr.log"
"$python_bin" -m ploidypatch.cli evidence build-copy-features \
    --consensus-decisions "$pool_decisions" \
    --method-decisions "miniprot=$method_root/methods/miniprot/complete_control/decisions.tsv" \
    --method-decisions "gemoma=$method_root/methods/gemoma/complete_control/decisions.tsv" \
    --method-decisions "lifton=$method_root/methods/lifton/complete_control/decisions.tsv" \
    --wgd-selection "$working_root/features/wgd_selection.tsv" \
    --output-tsv "$working_root/features/copy_features.tsv" \
    > "$working_root/features/copy_features.stdout.json" \
    2> "$working_root/features/copy_features.stderr.log"
"$python_bin" -m ploidypatch.cli evidence build-homeolog-topology-features \
    --copy-features "$working_root/features/copy_features.tsv" \
    --wgd-selection "$working_root/features/wgd_selection.tsv" \
    --candidate-gff "$candidate_gff" --base-gff "$base_gff" \
    --output-tsv "$working_root/features/topology_features.tsv" \
    > "$working_root/features/topology_features.stdout.json" \
    2> "$working_root/features/topology_features.stderr.log"
"$python_bin" -m ploidypatch.cli evidence score-support-conditioned-candidates \
    --copy-features "$working_root/features/copy_features.tsv" \
    --topology-features "$working_root/features/topology_features.tsv" \
    --model-json "$model_v03" --output-tsv "$working_root/scores/v03.tsv" \
    > "$working_root/scores/v03.stdout.json" \
    2> "$working_root/scores/v03.stderr.log"
"$python_bin" -m ploidypatch.cli evidence apply-conflict-winner-guard \
    --v03-scores "$working_root/scores/v03.tsv" \
    --pool-decisions "$pool_decisions" --pool-manifest "$pool_manifest" \
    --output-tsv "$working_root/scores/v04.tsv" \
    > "$working_root/scores/v04.stdout.json" \
    2> "$working_root/scores/v04.stderr.log"

"$python_bin" - "$working_root/scores/v04.tsv" \
    "$working_root/scores/v04.tsv.manifest.json" \
    "$working_root/natural/review_rankings.tsv" <<'PY'
import csv
import json
import math
import sys
from pathlib import Path

scores_path, manifest_path, output_path = map(Path, sys.argv[1:])
manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
with scores_path.open(encoding="utf-8", newline="") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))
if not rows or manifest.get("truth_access") is not False:
    raise SystemExit("natural score table is empty or not truth-blind")
audit = manifest.get("winner_audit", {})
if (
    manifest.get("counts", {}).get("automatic_approved") != 0
    or audit.get("mismatch_count") != 0
    or audit.get("baseline_mapping_sha256") != audit.get("v04_guard_mapping_sha256")
):
    raise SystemExit("natural v0.4 guard safety invariant failed")
estimators = (
    ("baseline", "v03_baseline_logit"),
    ("v03_primary", "v03_primary_rank_score"),
    ("v04_guard", "v04_primary_rank_score"),
)
output_path.parent.mkdir(parents=True, exist_ok=True)
with output_path.open("x", encoding="utf-8", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(("estimator", "review_rank", "candidate_digest", "rank_score"))
    for estimator, field in estimators:
        ordered = []
        for row in rows:
            score = float(row[field])
            if not math.isfinite(score):
                raise SystemExit(f"non-finite {field}")
            ordered.append((score, row["candidate_digest"]))
        ordered.sort(key=lambda item: (-item[0], item[1]))
        for rank, (score, digest) in enumerate(ordered, start=1):
            writer.writerow((estimator, rank, digest, format(score, ".17g")))
PY

for output in features/wgd_selection.tsv features/copy_features.tsv \
    features/topology_features.tsv scores/v03.tsv scores/v03.tsv.manifest.json \
    scores/v04.tsv scores/v04.tsv.manifest.json natural/review_rankings.tsv; do
    [[ -s $working_root/$output ]] || { echo "missing apple natural score output: $output" >&2; exit 1; }
done
(
    cd "$working_root"
    find . -type f ! -name SHA256SUMS -print0 | sort -z | xargs -0 sha256sum \
        > SHA256SUMS
    sha256sum -c SHA256SUMS >/dev/null
)
du -sb "$working_root" > "$working_root/disk_bytes.txt"
chmod -R a-w "$working_root"
mv "$working_root" "$result_root"
printf 'apple natural RNA-blind rankings frozen: %s\n' "$result_root"

