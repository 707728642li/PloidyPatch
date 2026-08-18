#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then echo "usage: $0 PROJECT_ROOT" >&2; exit 2; fi
project_root=$(realpath "$1")
source_root=$project_root/data/derived/external_inputs/populus_v0.4
source_gff=$source_root/normalized/target_populus/primary_chromosomes.gff3
source_genome=$source_root/normalized/target_populus/primary_chromosomes.genome.fa
benchmark=$project_root/benchmark/structure/copy_collapse_v0.4/ptr_v4/annotation_copy_collapse_seed20260930
blind_gff=$benchmark/blind/perturbed.gff3
truth=$benchmark/evaluator/truth/hidden_truth.json
evaluability=$benchmark/evaluator/pair_selection/evaluability.json
blind_run_root=$project_root/results/blind_runs/populus_external_v0.4
blind_project=$blind_run_root/project
upstream=$blind_project/results/baselines/populus_v0.4
method_root=$blind_project/results/copy_collapse/external/populus_v0.4_method_trio
ranking_root=$blind_project/results/copy_collapse/external/populus_v0.4_blind_rankings
blind_target_genome=$upstream/miniprot/normalized/target/primary_chromosomes.genome.fa
custody=$blind_run_root/custody_manifest.json
protocol_root=$project_root/results/protocol_freezes/populus_external_v0.4
execution_root=${PLOIDYPATCH_EXECUTION_FREEZE_OVERRIDE:-$project_root/results/protocol_freezes/populus_external_v0.4_execution}
execution_root=$(realpath "$execution_root")
case "$execution_root" in
    "$project_root/results/protocol_freezes/populus_external_v0.4_execution"*) ;;
    *) echo "execution freeze must be a Populus v0.4 project freeze" >&2; exit 2 ;;
esac
code_root=$execution_root/source
export PYTHONPATH="$code_root/src${PYTHONPATH:+:$PYTHONPATH}"
environment_bindings=$execution_root/environment_bindings.tsv
[[ -s $environment_bindings ]] || { echo "missing frozen environment bindings" >&2; exit 1; }
dev_prefix=$(awk -F '\t' '$1 == "ploidypatch-dev" {print $2}' "$environment_bindings")
[[ $dev_prefix == /* ]] || { echo "invalid frozen ploidypatch-dev binding" >&2; exit 1; }
python_bin=$dev_prefix/bin/python
policy=$protocol_root/policy.tsv
blind_scores=$ranking_root/scores/v04.tsv
blind_score_manifest=$ranking_root/scores/v04.tsv.manifest.json
pool_decisions=$method_root/consensus/primary_union/blind/decisions.tsv
pool_manifest=$method_root/consensus/primary_union/blind/candidate.gff3.manifest.json
result_root=$project_root/results/evaluator/populus/v0.4/reveal_inputs
working_root=${result_root}.working
self_relative=scripts/build_populus_complete_control_reveal_inputs_v0.4.sh

[[ -n ${PLOIDYPATCH_BLIND_RUN_ROOT:-} && -n ${PLOIDYPATCH_REVEAL_AUTHORIZATION:-} ]] || {
    echo "builder requires frozen blind-run root and reveal authorization" >&2; exit 1;
}
authorized_blind_root=$(realpath "$PLOIDYPATCH_BLIND_RUN_ROOT")
authorization=$(realpath "$PLOIDYPATCH_REVEAL_AUTHORIZATION")
[[ $authorized_blind_root == "$blind_run_root" && -s $authorization ]] || {
    echo "reveal authorization does not name the canonical frozen blind run" >&2; exit 1;
}

verify_implementation() {
    local relative=$1 manifest=$execution_root/implementation_manifest.tsv
    local rows=()
    mapfile -t rows < <(awk -F '\t' -v path="$relative" '$1 == path {print $2 "\t" $3}' "$manifest")
    [[ ${#rows[@]} -eq 1 ]] || { echo "execution freeze has no unique row for $relative" >&2; return 1; }
    local expected_bytes expected_sha
    IFS=$'\t' read -r expected_bytes expected_sha <<< "${rows[0]}"
    [[ $expected_bytes =~ ^[0-9]+$ && $expected_sha =~ ^[0-9a-f]{64}$ ]] || {
        echo "malformed execution implementation row for $relative" >&2; return 1;
    }
    [[ $(stat -Lc %s "$code_root/$relative") == "$expected_bytes" \
        && $(sha256sum "$code_root/$relative" | awk '{print $1}') == "$expected_sha" ]] || {
        echo "implementation differs from execution freeze: $relative" >&2; return 1;
    }
}
verify_tree() { (cd "$1" && sha256sum -c SHA256SUMS >/dev/null); }

implementation_dependencies=(
    "$self_relative"
    src/ploidypatch/baseline.py
    src/ploidypatch/cli.py
    src/ploidypatch/consensus.py
    src/ploidypatch/copy_features.py
    src/ploidypatch/score.py
)
# This first prerequisite block is deliberately truth-free.  Complete target
# annotation and hidden truth are not even stat'ed until custody validation.
for required in "$python_bin" "$protocol_root/SHA256SUMS" \
                "$policy" \
                "$execution_root/SHA256SUMS" "$execution_root/implementation_manifest.tsv" \
                "$blind_run_root/SHA256SUMS" "$custody" "$method_root/SHA256SUMS" \
                "$ranking_root/SHA256SUMS" "$benchmark/blind/SHA256SUMS" \
                "$blind_scores" "$blind_score_manifest" "$pool_decisions" "$pool_manifest" \
                "${implementation_dependencies[@]/#/$code_root/}"; do
    [[ -s $required ]] || { echo "missing pre-reveal custody prerequisite: $required" >&2; exit 1; }
done
verify_tree "$protocol_root"
verify_tree "$execution_root"
verify_tree "$blind_run_root"
verify_tree "$method_root"
verify_tree "$ranking_root"
verify_tree "$benchmark/blind"
for relative in "${implementation_dependencies[@]}"; do verify_implementation "$relative"; done
policy_value() {
    local key=$1
    awk -F '\t' -v key="$key" '$1 == key {print $2}' "$policy"
}
[[ $(policy_value policy_id) == ploidypatch_populus_external_validation_v0.4 \
    && $(policy_value adapter_min_identity) == 0.5 \
    && $(policy_value adapter_min_query_coverage) == 0.5 \
    && $(policy_value adapter_require_intact) == true \
    && $(policy_value adapter_max_existing_cds_overlap) == 0.2 \
    && $(policy_value adapter_max_redundancy_overlap) == 0.5 \
    && $(policy_value primary_candidate_policy) == retain_distinct_phased_CDS_chains \
    && $(policy_value legacy_candidate_comparator) == suppress_strongly_overlapping_alternative_chains \
    && $(policy_value automatic_copy_addition_approval) == false ]] || {
    echo "complete-control adaptation differs from the frozen Populus policy" >&2; exit 1;
}

"$python_bin" - "$authorization" "$blind_run_root/SHA256SUMS" "$custody" \
    "$blind_scores" "$blind_score_manifest" "$pool_decisions" "$pool_manifest" <<'PY'
import hashlib
import json
import sys
from pathlib import Path

authorization_path, blind_sums, custody_path, scores, score_manifest, decisions, pool_manifest = map(Path, sys.argv[1:])
def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()
with custody_path.open(encoding="utf-8") as handle:
    custody = json.load(handle)
with authorization_path.open(encoding="utf-8") as handle:
    authorization = json.load(handle)
forbidden = (
    "truth_mounted", "complete_target_annotation_mounted",
    "evaluator_references_mounted", "nas_data_mounted", "network_access",
)
outputs = custody.get("blind_outputs", {})
expected_authorization = {
    "blind_run_SHA256SUMS_sha256": sha256(blind_sums),
    "custody_manifest_sha256": sha256(custody_path),
    "blind_scores_sha256": sha256(scores),
    "blind_score_manifest_sha256": sha256(score_manifest),
    "pool_decisions_sha256": sha256(decisions),
    "pool_manifest_sha256": sha256(pool_manifest),
}
if (
    custody.get("schema_version") != "ploidypatch.blind_run_custody.v1"
    or any(custody.get(field) is not False for field in forbidden)
    or not custody.get("frozen_before_truth_reveal_at")
    or outputs.get("scores_sha256") != sha256(scores)
    or outputs.get("score_manifest_sha256") != sha256(score_manifest)
    or outputs.get("pool_decisions_sha256") != sha256(decisions)
    or outputs.get("pool_manifest_sha256") != sha256(pool_manifest)
    or authorization.get("schema_version") != "ploidypatch.populus_reveal_authorization.v0.4"
    or authorization.get("truth_opened") is not False
    or any(authorization.get(key) != value for key, value in expected_authorization.items())
):
    raise SystemExit("blind custody is absent, late, or disagrees with frozen blind artifacts")
with score_manifest.open(encoding="utf-8") as handle:
    score = json.load(handle)
with pool_manifest.open(encoding="utf-8") as handle:
    pool = json.load(handle)
if (
    score.get("schema_version") != "ploidypatch.conflict_winner_guard_scores.v1"
    or score.get("truth_access") is not False
    or score.get("outputs", {}).get("scores", {}).get("sha256") != sha256(scores)
    or score.get("inputs", {}).get("pool_decisions") != sha256(decisions)
    or score.get("inputs", {}).get("pool_manifest") != sha256(pool_manifest)
    or score.get("winner_audit", {}).get("mismatch_count") != 0
    or pool.get("schema_version") != "ploidypatch.method_candidate_pool.v2"
    or pool.get("outputs", {}).get("decisions", {}).get("sha256") != sha256(decisions)
):
    raise SystemExit("blind scores or candidate pool fail the frozen pre-reveal contract")
PY

# Custody has now passed.  Evaluator-owned source, truth and controls may be read.
for required in "$source_root/EVALUATOR_SHA256SUMS" "$source_gff" "$source_genome" \
                "$benchmark/SHA256SUMS" "$benchmark/blind/blind_manifest.json" \
                "$blind_gff" "$truth" "$evaluability" \
                "$upstream/miniprot/SHA256SUMS" \
                "$blind_target_genome" \
                "$upstream/miniprot/raw/miniprot.gff3" \
                "$upstream/miniprot/reference/populus_candidate_refs.map.tsv" \
                "$method_root/merged/gemoma.gff3" "$method_root/merged/lifton.gff3" \
                "$ranking_root/features/copy_features.tsv"; do
    [[ -s $required ]] || { echo "missing evaluator-owned reveal prerequisite: $required" >&2; exit 1; }
done
(cd "$source_root" && sha256sum -c EVALUATOR_SHA256SUMS >/dev/null)
verify_tree "$benchmark"
verify_tree "$upstream/miniprot"
[[ ! -e $result_root && ! -e $working_root ]] || {
    echo "refusing to overwrite Populus reveal inputs" >&2; exit 1;
}
mkdir -p "$working_root"/{complete_control/methods,complete_control/consensus,labels,scores/methods,scores/consensus}
for method in miniprot gemoma lifton; do mkdir -p "$working_root/complete_control/methods/$method"; done
for pool in primary_union legacy_union support2 support3; do mkdir -p "$working_root/complete_control/consensus/$pool"; done
record_invalid() {
    local status=$?
    if [[ -d $working_root ]]; then
        printf '{"status":"invalid_run","stage":"complete_control_reveal_builder","exit_status":%s}\n' \
            "$status" > "$working_root/invalid_run.json" || true
    fi
    exit "$status"
}
trap record_invalid ERR

formal_outcome=$("$python_bin" - "$evaluability" <<'PY'
import json
import sys
with open(sys.argv[1], encoding="utf-8") as handle:
    report = json.load(handle)
print(report.get("formal_outcome", "invalid_run"))
PY
)
genome_sentinel=$("$python_bin" - "$source_genome" "$blind_target_genome" \
    "$benchmark/blind/blind_manifest.json" <<'PY'
import hashlib
import json
import sys
from pathlib import Path
source, blind, manifest_path = map(Path, sys.argv[1:])
def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()
with manifest_path.open(encoding="utf-8") as handle:
    expected = json.load(handle).get("target_genome", {}).get("sha256")
print("pass" if expected == sha256(source) == sha256(blind) else "fail")
PY
)
[[ $genome_sentinel == pass ]] || formal_outcome=invalid_run
{
    printf 'field\tvalue\n'
    printf 'code_commit\t%s\n' "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'access\tevaluator_only_after_blind_score_and_custody_freeze\n'
    printf 'generated_after_blind_freeze\ttrue\nevaluator_only\ttrue\n'
    printf 'blind_raw_predictions_reused_without_rerun\ttrue\n'
    printf 'blind_and_complete_control_adaptation\tindependent_same_frozen_raw_predictions\n'
    printf 'model_refit\tfalse\nthreshold_tuning\tfalse\nreference_change\tfalse\n'
    printf 'formal_evaluator_environment\tploidypatch-model\n'
    printf 'automatic_approval\tfalse\nbenchmark_formal_outcome\t%s\n' "$formal_outcome"
    printf 'blind_complete_target_genome_sentinel\t%s\n' "$genome_sentinel"
    printf 'custody_manifest_sha256\t%s\n' "$(sha256sum "$custody" | awk '{print $1}')"
    printf 'reveal_authorization_sha256\t%s\n' "$(sha256sum "$authorization" | awk '{print $1}')"
    printf 'blind_run_sha256sums_sha256\t%s\n' \
        "$(sha256sum "$blind_run_root/SHA256SUMS" | awk '{print $1}')"
} > "$working_root/run_contract.tsv"

write_manifest() {
    local status=$1
    local reason=$2
    "$python_bin" - "$working_root" "$status" "$reason" \
        "$blind_run_root/SHA256SUMS" "$custody" \
        "$blind_scores" "$blind_score_manifest" "$pool_decisions" "$pool_manifest" \
        "$benchmark/SHA256SUMS" "$truth" "$evaluability" "$authorization" \
        "$method_root/SHA256SUMS" "$ranking_root/SHA256SUMS" <<'PY'
import hashlib
import json
import sys
from pathlib import Path

(
    root, status, reason, blind_sums, custody, scores, score_manifest, decisions,
    pool_manifest, benchmark_sums, truth, evaluability, authorization,
    method_sums, ranking_sums,
) = sys.argv[1:]
root = Path(root)
def sha256(value: str | Path) -> str:
    path = Path(value)
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()
manifest = {
    "schema_version": "ploidypatch.populus_reveal_inputs.v0.4",
    "generated_after_blind_freeze": True,
    "evaluator_only": True,
    "status": status,
    "reason": reason,
    "blind_run_SHA256SUMS_sha256": sha256(blind_sums),
    "custody_manifest_sha256": sha256(custody),
    "reveal_authorization_sha256": sha256(authorization),
    "blind_scores_sha256": sha256(scores),
    "blind_score_manifest_sha256": sha256(score_manifest),
    "pool_decisions_sha256": sha256(decisions),
    "pool_manifest_sha256": sha256(pool_manifest),
    "benchmark_SHA256SUMS_sha256": sha256(benchmark_sums),
    "hidden_truth_sha256": sha256(truth),
    "evaluability_sha256": sha256(evaluability),
    "blind_method_pool_SHA256SUMS_sha256": sha256(method_sums),
    "blind_ranking_SHA256SUMS_sha256": sha256(ranking_sums),
    "automatic_approval": False,
    "formal_evaluator_environment": "ploidypatch-model",
}
if status == "ready_for_evaluation":
    paths = {
        "labels": "labels/candidate_labels.tsv",
        "labels_manifest": "labels/candidate_labels.tsv.manifest.json",
        "primary_pool_score": "scores/consensus/primary_union.json",
        "legacy_pool_score": "scores/consensus/legacy_union.json",
        "secondary:miniprot": "scores/methods/miniprot.json",
        "secondary:gemoma": "scores/methods/gemoma.json",
        "secondary:lifton": "scores/methods/lifton.json",
        "secondary:support2": "scores/consensus/support2.json",
        "secondary:support3": "scores/consensus/support3.json",
    }
    manifest["evaluation_inputs"] = {
        role: {"relative_path": relative, "sha256": sha256(root / relative)}
        for role, relative in paths.items()
    }
    control_paths = sorted(
        path for path in (root / "complete_control").rglob("*")
        if path.is_file() and (
            path.name in {"candidate.gff3", "decisions.tsv"}
            or path.name.endswith(".manifest.json")
        )
    )
    manifest["complete_control_artifacts"] = {
        path.relative_to(root).as_posix(): sha256(path) for path in control_paths
    }
with (root / "reveal_input_manifest.json").open("x", encoding="utf-8") as handle:
    json.dump(manifest, handle, indent=2, sort_keys=True)
    handle.write("\n")
with (root / "status.json").open("x", encoding="utf-8") as handle:
    json.dump(
        {
            "schema_version": "ploidypatch.populus_reveal_input_status.v0.4",
            "status": status,
            "reason": reason,
            "generated_after_blind_freeze": True,
            "evaluator_only": True,
        },
        handle,
        indent=2,
        sort_keys=True,
    )
    handle.write("\n")
PY
}
finalize_tree() {
    du -sb "$working_root" > "$working_root/disk_bytes.txt"
    (
        cd "$working_root"
        find . -type f \( -name '*.gff3' -o -name '*.json' -o -name '*.tsv' \
            -o -name disk_bytes.txt \) \
            ! -name SHA256SUMS -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS
        sha256sum -c SHA256SUMS >/dev/null
    )
    trap - ERR
    mv "$working_root" "$result_root"
}

if [[ $formal_outcome == invalid_run ]]; then
    write_manifest invalid_run benchmark_or_genome_sentinel_invalid
    finalize_tree
    echo "benchmark is invalid; frozen reveal builder refusal: $result_root" >&2
    exit 1
elif [[ $formal_outcome == not_evaluable_without_rule_relaxation ]]; then
    write_manifest not_evaluable_without_rule_relaxation benchmark_fixed_data_gates_failed
    finalize_tree
    printf 'Populus reveal builder retained fixed-rule non-evaluable outcome: %s\n' "$result_root"
    exit 0
elif [[ $formal_outcome != formally_evaluable_pending_blind_and_complete_control_reveal ]]; then
    echo "unrecognized Populus benchmark outcome: $formal_outcome" >&2
    exit 1
fi

cd "$code_root"
miniprot_gff=$upstream/miniprot/raw/miniprot.gff3
protein_map=$upstream/miniprot/reference/populus_candidate_refs.map.tsv
"$python_bin" -m ploidypatch.cli baseline adapt-miniprot \
    --perturbed-gff "$source_gff" --miniprot-gff "$miniprot_gff" \
    --protein-map "$protein_map" --min-identity 0.5 --min-query-coverage 0.5 \
    --max-existing-cds-overlap 0.2 --max-redundancy-overlap 0.5 \
    --output-gff "$working_root/complete_control/methods/miniprot/candidate.gff3" \
    --decisions-tsv "$working_root/complete_control/methods/miniprot/decisions.tsv" \
    > "$working_root/complete_control/methods/miniprot/stdout.json" \
    2> "$working_root/complete_control/methods/miniprot/stderr.log"
for method in gemoma lifton; do
    "$python_bin" -m ploidypatch.cli baseline adapt-gff \
        --perturbed-gff "$source_gff" --candidate-gff "$method_root/merged/$method.gff3" \
        --source "$method" --max-existing-cds-overlap 0.2 --max-redundancy-overlap 0.5 \
        --output-gff "$working_root/complete_control/methods/$method/candidate.gff3" \
        --decisions-tsv "$working_root/complete_control/methods/$method/decisions.tsv" \
        > "$working_root/complete_control/methods/$method/stdout.json" \
        2> "$working_root/complete_control/methods/$method/stderr.log"
done

pids=(); labels=()
for pool in primary_union legacy_union support2 support3; do
    case $pool in
        primary_union) support=1; redundancy=retain_distinct_chains ;;
        legacy_union) support=1; redundancy=suppress_overlapping ;;
        support2) support=2; redundancy=suppress_overlapping ;;
        support3) support=3; redundancy=suppress_overlapping ;;
    esac
    (
        "$python_bin" -m ploidypatch.cli baseline select-method-consensus \
            --base-gff "$source_gff" \
            --candidate "miniprot=$working_root/complete_control/methods/miniprot/candidate.gff3" \
            --candidate "gemoma=$working_root/complete_control/methods/gemoma/candidate.gff3" \
            --candidate "lifton=$working_root/complete_control/methods/lifton/candidate.gff3" \
            --min-method-support "$support" --max-redundancy-overlap 0.5 \
            --redundancy-policy "$redundancy" \
            --output-gff "$working_root/complete_control/consensus/$pool/candidate.gff3" \
            --decisions-tsv "$working_root/complete_control/consensus/$pool/decisions.tsv" \
            > "$working_root/complete_control/consensus/$pool/stdout.json" \
            2> "$working_root/complete_control/consensus/$pool/stderr.log"
    ) & pids+=("$!"); labels+=("pool:$pool")
done
failed=0
for index in "${!pids[@]}"; do
    if ! wait "${pids[$index]}"; then echo "failed complete-control ${labels[$index]}" >&2; failed=1; fi
done
[[ $failed -eq 0 ]] || exit 1

raw_labels=$working_root/labels/candidate_labels.raw.tsv
"$python_bin" -m ploidypatch.cli benchmark label-copy-features \
    --features "$ranking_root/features/copy_features.tsv" --truth "$truth" \
    --output-tsv "$raw_labels" \
    > "$working_root/labels/label.stdout.json" 2> "$working_root/labels/label.stderr.log"
"$python_bin" - "$raw_labels" "$blind_scores" "$pool_decisions" "$pool_manifest" \
    "$working_root/labels/candidate_labels.tsv" <<'PY'
import csv
import hashlib
import json
import sys
from pathlib import Path

raw_path, score_path, decisions_path, pool_manifest_path, output = map(Path, sys.argv[1:])
def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()
def read(path: Path):
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return reader.fieldnames or [], list(reader)
raw_fields, raw = read(raw_path)
score_fields, scores = read(score_path)
decision_fields, decisions = read(decisions_path)
if not {"candidate_digest", "label_exact_cds"} <= set(raw_fields):
    raise SystemExit("primitive candidate labels lack exact-CDS fields")
raw_by_digest = {row["candidate_digest"]: row for row in raw}
score_digests = {row["candidate_digest"] for row in scores}
accepted = {
    row.get("candidate_digest") or row.get("consensus_digest")
    for row in decisions if row.get("status") == "accepted"
}
if (
    len(raw_by_digest) != len(raw) or len(score_digests) != len(scores)
    or "" in accepted or score_digests != set(raw_by_digest) or score_digests != accepted
):
    raise SystemExit("label, blind-score and accepted-pool candidate universes differ")
with output.open("x", encoding="utf-8", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(("candidate_digest", "label_exact_cds"))
    for digest in sorted(score_digests):
        label = raw_by_digest[digest]["label_exact_cds"]
        if label not in {"0", "1"}:
            raise SystemExit("non-binary exact-CDS label")
        writer.writerow((digest, label))
with pool_manifest_path.open(encoding="utf-8") as handle:
    pool = json.load(handle)
manifest = {
    "schema_version": "ploidypatch.external_candidate_labels.v1",
    "evaluator_only": True,
    "generated_after_blind_freeze": True,
    "blind_scores_sha256": sha256(score_path),
    "pool_manifest_sha256": sha256(pool_manifest_path),
    "pool_decisions_sha256": sha256(decisions_path),
    "primitive_label_manifest_sha256": sha256(Path(str(raw_path) + ".manifest.json")),
    "counts": {
        "candidates": len(raw),
        "positive_exact_cds": sum(row["label_exact_cds"] == "1" for row in raw),
        "negative_candidates": sum(row["label_exact_cds"] == "0" for row in raw),
    },
    "outputs": {"labels": {"file_name": output.name, "rows": len(raw), "sha256": sha256(output)}},
}
with Path(str(output) + ".manifest.json").open("x", encoding="utf-8") as handle:
    json.dump(manifest, handle, indent=2, sort_keys=True)
    handle.write("\n")
PY

pids=(); labels=()
for method in miniprot gemoma lifton; do
    (
        "$python_bin" -m ploidypatch.cli benchmark score \
            --source-gff "$source_gff" --perturbed-gff "$blind_gff" \
            --candidate-gff "$method_root/methods/$method/blind/candidate.gff3" \
            --control-candidate-gff "$working_root/complete_control/methods/$method/candidate.gff3" \
            --truth "$truth" --include-event-details \
            > "$working_root/scores/methods/$method.json" \
            2> "$working_root/scores/methods/$method.stderr.log"
    ) & pids+=("$!"); labels+=("method:$method")
done
for pool in primary_union legacy_union support2 support3; do
    (
        "$python_bin" -m ploidypatch.cli benchmark score \
            --source-gff "$source_gff" --perturbed-gff "$blind_gff" \
            --candidate-gff "$method_root/consensus/$pool/blind/candidate.gff3" \
            --control-candidate-gff "$working_root/complete_control/consensus/$pool/candidate.gff3" \
            --truth "$truth" --include-event-details \
            > "$working_root/scores/consensus/$pool.json" \
            2> "$working_root/scores/consensus/$pool.stderr.log"
    ) & pids+=("$!"); labels+=("pool-score:$pool")
done
failed=0
for index in "${!pids[@]}"; do
    if ! wait "${pids[$index]}"; then echo "failed evaluator score ${labels[$index]}" >&2; failed=1; fi
done
[[ $failed -eq 0 ]] || exit 1
"$python_bin" - "$working_root" <<'PY'
import json
import sys
from pathlib import Path
root = Path(sys.argv[1])
for path in sorted((root / "scores").rglob("*.json")):
    with path.open(encoding="utf-8") as handle:
        report = json.load(handle)
    if report.get("quality_gate", {}).get("grade") != "pass":
        raise SystemExit(f"score quality gate failed: {path}")
PY
label_status=$("$python_bin" - "$working_root/labels/candidate_labels.tsv.manifest.json" <<'PY'
import json
import sys
with open(sys.argv[1], encoding="utf-8") as handle:
    counts = json.load(handle)["counts"]
print("ready_for_evaluation" if counts["positive_exact_cds"] > 0 and counts["negative_candidates"] > 0 else "not_evaluable_without_rule_relaxation")
PY
)
label_reason=all_reveal_inputs_valid
[[ $label_status == ready_for_evaluation ]] || label_reason=candidate_labels_lack_both_classes
write_manifest "$label_status" "$label_reason"
finalize_tree
printf 'Populus complete-control reveal inputs frozen (%s): %s\n' "$label_status" "$result_root"
