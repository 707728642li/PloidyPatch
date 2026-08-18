#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then echo "usage: $0 PROJECT_ROOT" >&2; exit 2; fi
project_root=$(realpath "$1")
code_root=$project_root/code
python_bin=$project_root/envs/ploidypatch-dev/bin/python
input_root=${PLOIDYPATCH_STAGED_INPUT_ROOT:?PLOIDYPATCH_STAGED_INPUT_ROOT is required}
blind_benchmark_root=${PLOIDYPATCH_BLIND_BENCHMARK_ROOT:?PLOIDYPATCH_BLIND_BENCHMARK_ROOT is required}
base_gff=$blind_benchmark_root/perturbed.gff3
method_root=$project_root/results/copy_collapse/external/actinidia_v0.5_method_trio
candidate_gff=$method_root/consensus/primary_union/blind/candidate.gff3
pool_decisions=$method_root/consensus/primary_union/blind/decisions.tsv
pool_manifest=$method_root/consensus/primary_union/blind/candidate.gff3.manifest.json
self_wgd_root=$project_root/results/copy_collapse/external/actinidia_v0.5_blind_self_wgd
prior_wgd=$self_wgd_root/selected/selection.tsv
protocol_root=${PLOIDYPATCH_PROTOCOL_FREEZE:?PLOIDYPATCH_PROTOCOL_FREEZE is required}
execution_root=${PLOIDYPATCH_EXECUTION_FREEZE:?PLOIDYPATCH_EXECUTION_FREEZE is required}
contract_path=${PLOIDYPATCH_HOLDOUT_CONTRACT:?PLOIDYPATCH_HOLDOUT_CONTRACT is required}
context_verifier=$code_root/scripts/verify_external_holdout_blind_context_v0.5.py
policy=$protocol_root/protocol_artifacts/config/actinidia_external_validation_policy_v0.5.tsv
composite_root=${PLOIDYPATCH_COMPOSITE_MODEL_FREEZE:?PLOIDYPATCH_COMPOSITE_MODEL_FREEZE is required}
model_v03=$composite_root/model_v0.3.json
guard_policy=$composite_root/guard_policy.json
result_root=$project_root/results/copy_collapse/external/actinidia_v0.5_blind_rankings
working_root=${result_root}.working

[[ ${PLOIDYPATCH_BLIND_RUNNER:-} == 1 ]] || {
    echo "Actinidia scoring must run inside the frozen blind runner" >&2; exit 1;
}
[[ ${PLOIDYPATCH_NETWORK_ACCESS:-} == none ]] || {
    echo "Actinidia scoring requires a network-disabled namespace" >&2; exit 1;
}
for forbidden in /nas_data "$input_root/evaluator_only" "$blind_benchmark_root/evaluator" \
    "$blind_benchmark_root/truth" "$blind_benchmark_root/complete"; do
    [[ ! -e $forbidden ]] || { echo "forbidden blind-runner path is visible: $forbidden" >&2; exit 1; }
done
if [[ -r /proc/self/mountinfo ]] && grep -Eq '/nas_data|/evaluator_only|/target_complete|/truth_references' /proc/self/mountinfo; then
    echo "forbidden evaluator or NAS mount detected in blind namespace" >&2; exit 1
fi

verify_tree() { (cd "$1" && sha256sum -c SHA256SUMS >/dev/null); }
verify_implementation() {
    local relative=$1 expected observed
    expected=$(awk -F '\t' -v path="$relative" 'NR > 1 && $1 == path {print $3}' \
        "$execution_root/implementation_manifest.tsv")
    observed=$(sha256sum "$code_root/$relative" | awk '{print $1}')
    [[ $expected =~ ^[0-9a-f]{64}$ && $observed == "$expected" ]] || {
        echo "implementation differs from execution freeze: $relative" >&2; exit 1;
    }
}

for required in "$python_bin" "$base_gff" "$method_root/SHA256SUMS" \
    "$blind_benchmark_root/blind_manifest.json" "$blind_benchmark_root/SHA256SUMS" \
    "$self_wgd_root/SHA256SUMS" "$protocol_root/SHA256SUMS" \
    "$execution_root/SHA256SUMS" "$composite_root/SHA256SUMS" "$policy" \
    "$candidate_gff" "$pool_decisions" "$pool_manifest" "$prior_wgd" \
    "$model_v03" "$guard_policy" "$composite_root/composite_manifest.json" \
    "$contract_path" "$input_root/role_manifest.tsv" \
    "$input_root/role_contract.json" "$context_verifier" \
    "$method_root/methods/miniprot/blind/decisions.tsv" \
    "$method_root/methods/gemoma/blind/decisions.tsv" \
    "$method_root/methods/lifton/blind/decisions.tsv"; do
    [[ -s $required ]] || { echo "missing Actinidia blind ranking input: $required" >&2; exit 1; }
done
for root in "$method_root" "$self_wgd_root" "$protocol_root" "$execution_root" \
            "$composite_root"; do
    verify_tree "$root"
done
verify_implementation scripts/score_actinidia_candidates_blind_v0.5.sh
verify_implementation scripts/verify_external_holdout_blind_context_v0.5.py
for relative in src/ploidypatch/copy_features.py \
    src/ploidypatch/homeolog_topology.py src/ploidypatch/wgd_candidate_select.py \
    src/ploidypatch/support_ranker.py src/ploidypatch/conflict_guard.py \
    src/ploidypatch/cli.py; do
    verify_implementation "$relative"
done

"$python_bin" - "$composite_root" "$policy" <<'PY'
import csv
import hashlib
import json
import sys
from pathlib import Path

root = Path(sys.argv[1])
policy_path = Path(sys.argv[2])

def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()

with policy_path.open(encoding="utf-8", newline="") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))
policy = {row["field"]: row["value"] for row in rows}
if len(policy) != len(rows):
    raise SystemExit("duplicate Actinidia policy field")
if (
    policy.get("policy_id") != "ploidypatch_actinidia_external_validation_v0.5"
    or policy.get("model_version") != "PloidyPatch_ranker_v0.4"
    or policy.get("estimator_version") != "PloidyPatch_ranker_v0.4"
    or policy.get("multiple_references_per_method_vote") != "one_method_family_vote"
    or policy.get("candidate_truth_access")
    != "false_until_blind_scores_pool_manifests_and_custody_hashes_are_frozen"
    or policy.get("automatic_copy_addition_approval") != "false"
):
    raise SystemExit("Actinidia policy/composite binding differs")
manifest = json.loads((root / "composite_manifest.json").read_text(encoding="utf-8"))
guard = json.loads((root / "guard_policy.json").read_text(encoding="utf-8"))
model_sha = sha(root / "model_v0.3.json")
if (
    manifest.get("schema_version") != "ploidypatch.composite_ranker.v0.4"
    or manifest.get("truth_access") is not False
    or manifest.get("automatic_approval") is not False
    or manifest.get("components", {}).get("model_v0.3") != model_sha
    or manifest.get("components", {}).get("guard_policy") != sha(root / "guard_policy.json")
):
    raise SystemExit("composite ranker manifest fails identity or blind-use gate")
required_pool = guard.get("required_pool", {})
guard_rule = guard.get("policy", {})
if (
    guard.get("schema_version") != "ploidypatch.conflict_winner_guard_policy.v1"
    or guard.get("base_model", {}).get("sha256") != model_sha
    or required_pool.get("schema_version") != "ploidypatch.method_candidate_pool.v2"
    or required_pool.get("conflict_action") != "retain_all_for_ranking_and_review"
    or required_pool.get("max_redundancy_overlap") != 0.5
    or required_pool.get("min_method_support") != 1
    or required_pool.get("redundancy_policy") != "retain_distinct_chains"
    or guard_rule.get("automatic_approval") is not False
    or guard_rule.get("calibrated_probability") is not False
):
    raise SystemExit("frozen v0.4 guard policy differs")
PY

[[ ! -e $result_root && ! -e $working_root ]] || {
    echo "refusing to overwrite Actinidia blind rankings" >&2; exit 1;
}
mkdir -p "$working_root"/{features,scores,freeze}
cp /proc/self/mountinfo "$working_root/freeze/blind_runner.mountinfo"
PYTHONPATH="$code_root/src" "$python_bin" "$context_verifier" \
    --input-root "$input_root" --contract "$contract_path" \
    --protocol-freeze "$protocol_root" --execution-freeze "$execution_root" \
    --blind-benchmark-root "$blind_benchmark_root" \
    --expected-holdout-id actinidia_red5_v0.5 --expected-primary-chromosomes 29 \
    --model-freeze "$composite_root" \
    --output-json "$working_root/freeze/blind_context.json"
{
    printf 'field\tvalue\ncode_commit\t%s\n' "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'split\ttarget_level_predeclared_untouched_secondary_replication\ntruth_access\tfalse\n'
    printf 'hidden_pair_access\tfalse\nexternal_label_access\tfalse\n'
    printf 'complete_target_annotation_access\tfalse\nevaluator_reference_access\tfalse\n'
    printf 'candidate_policy\tretain_distinct_phased_CDS_chains\n'
    printf 'wgd_source\tblind_candidate_recomputation_only\n'
    printf 'v03_model_sha256\t%s\n' "$(sha256sum "$model_v03" | awk '{print $1}')"
    printf 'v04_guard_policy_sha256\t%s\n' "$(sha256sum "$guard_policy" | awk '{print $1}')"
    printf 'composite_sha256sums_sha256\t%s\n' \
        "$(sha256sum "$composite_root/SHA256SUMS" | awk '{print $1}')"
    printf 'automatic_approval\tfalse\ncalibration\tnone_rank_scores_only\n'
    printf 'protocol_freeze_sha256sums_sha256\t%s\n' \
        "$(sha256sum "$protocol_root/SHA256SUMS" | awk '{print $1}')"
    printf 'execution_freeze_sha256sums_sha256\t%s\n' \
        "$(sha256sum "$execution_root/SHA256SUMS" | awk '{print $1}')"
} > "$working_root/run_contract.tsv"
{
    printf 'role\tbytes\tsha256\tpath\n'
    for entry in "base_gff:$base_gff" "candidate_gff:$candidate_gff" \
        "pool_decisions:$pool_decisions" "pool_manifest:$pool_manifest" \
        "prior_wgd:$prior_wgd" "composite_manifest:$composite_root/composite_manifest.json" \
        "model_v03:$model_v03" "guard_policy:$guard_policy"; do
        role=${entry%%:*}; path=${entry#*:}
        printf '%s\t%s\t%s\t%s\n' "$role" "$(stat -Lc %s "$path")" \
            "$(sha256sum "$path" | awk '{print $1}')" "$path"
    done
} > "$working_root/input_manifest.tsv"

cd "$code_root"
"$python_bin" -m ploidypatch.cli evidence propagate-wgd-conflict-partners \
    --base-gff "$base_gff" --candidate-gff "$candidate_gff" \
    --pool-decisions "$pool_decisions" --prior-wgd-selection "$prior_wgd" \
    --output-selection "$working_root/features/wgd_selection.tsv" \
    > "$working_root/features/wgd_selection.stdout.json" \
    2> "$working_root/features/wgd_selection.stderr.log"
"$python_bin" -m ploidypatch.cli evidence build-copy-features \
    --consensus-decisions "$pool_decisions" \
    --method-decisions "miniprot=$method_root/methods/miniprot/blind/decisions.tsv" \
    --method-decisions "gemoma=$method_root/methods/gemoma/blind/decisions.tsv" \
    --method-decisions "lifton=$method_root/methods/lifton/blind/decisions.tsv" \
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
{
    printf 'feature\tbytes\tsha256\n'
    for path in "$working_root/features/wgd_selection.tsv" \
        "$working_root/features/copy_features.tsv" \
        "$working_root/features/copy_features.tsv.manifest.json" \
        "$working_root/features/topology_features.tsv" \
        "$working_root/features/topology_features.tsv.manifest.json"; do
        [[ -s $path ]] || { echo "missing blind feature artifact: $path" >&2; exit 1; }
        printf '%s\t%s\t%s\n' "$(basename "$path")" "$(stat -Lc %s "$path")" \
            "$(sha256sum "$path" | awk '{print $1}')"
    done
} > "$working_root/freeze/blind_feature_freeze.tsv"

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
    "$working_root/scores/v04.tsv.manifest.json" <<'PY'
import csv
import hashlib
import json
import sys
from pathlib import Path

scores = Path(sys.argv[1])
manifest_path = Path(sys.argv[2])
manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
with scores.open(encoding="utf-8", newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    fields = reader.fieldnames or []
    rows = list(reader)
if not rows or "v04_automatic_approval" not in fields:
    raise SystemExit("frozen v0.4 score table is empty or lacks approval boundary")
if any("truth" in field.lower() or "label" in field.lower() for field in fields):
    raise SystemExit("truth/label field leaked into blind frozen v0.4 scores")
if any(row["v04_automatic_approval"] != "0" for row in rows):
    raise SystemExit("frozen v0.4 guard emitted an automatic approval")
score_sha = hashlib.sha256(scores.read_bytes()).hexdigest()
audit = manifest.get("winner_audit", {})
counts = manifest.get("counts", {})
if (
    manifest.get("schema_version") != "ploidypatch.conflict_winner_guard_scores.v1"
    or manifest.get("truth_access") is not False
    or manifest.get("policy", {}).get("automatic_approval") is not False
    or counts.get("automatic_approved") != 0
    or counts.get("winner_mismatch_count") != 0
    or audit.get("mismatch_count") != 0
    or audit.get("baseline_mapping_sha256") != audit.get("v04_guard_mapping_sha256")
    or manifest.get("outputs", {}).get("scores", {}).get("sha256") != score_sha
    or manifest.get("outputs", {}).get("scores", {}).get("rows") != len(rows)
):
    raise SystemExit("frozen v0.4 production guard failed its blind safety contract")
PY

for output in scores/v03.tsv scores/v03.tsv.manifest.json scores/v04.tsv \
    scores/v04.tsv.manifest.json freeze/blind_feature_freeze.tsv; do
    [[ -s $working_root/$output ]] || { echo "missing Actinidia blind score output: $output" >&2; exit 1; }
done
{
    printf 'score\tbytes\tsha256\n'
    for path in "$working_root"/scores/*.tsv "$working_root"/scores/*.manifest.json; do
        printf '%s\t%s\t%s\n' "$(basename "$path")" "$(stat -Lc %s "$path")" \
            "$(sha256sum "$path" | awk '{print $1}')"
    done
} > "$working_root/freeze/blind_score_freeze.tsv"
(
    cd "$working_root"
    find . -type f \( -name '*.json' -o -name '*.tsv' -o -name '*.mountinfo' \) \
        -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS
    sha256sum -c SHA256SUMS >/dev/null
)
du -sb "$working_root" > "$working_root/disk_bytes.txt"
mv "$working_root" "$result_root"
printf 'Actinidia truth-blind v0.5 rankings frozen: %s\n' "$result_root"
