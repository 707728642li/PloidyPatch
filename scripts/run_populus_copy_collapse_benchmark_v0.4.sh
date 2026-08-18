#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then echo "usage: $0 PROJECT_ROOT" >&2; exit 2; fi
project_root=$(realpath "$1")
source_root=$project_root/data/derived/external_inputs/populus_v0.4
source_gff=$source_root/normalized/target_populus/primary_chromosomes.gff3
source_genome=$source_root/normalized/target_populus/primary_chromosomes.genome.fa
pair_root=$project_root/results/evaluator/populus/v0.4/truth_pairs
pair_tsv=$pair_root/intersection/pairs.tsv
protocol_root=$project_root/results/protocol_freezes/populus_external_v0.4
execution_root=$project_root/results/protocol_freezes/populus_external_v0.4_execution
code_root=$execution_root/source
export PYTHONPATH="$code_root/src${PYTHONPATH:+:$PYTHONPATH}"
environment_bindings=$execution_root/environment_bindings.tsv
[[ -s $environment_bindings ]] || { echo "missing frozen environment bindings" >&2; exit 1; }
dev_prefix=$(awk -F '\t' '$1 == "ploidypatch-dev" {print $2}' "$environment_bindings")
[[ $dev_prefix == /* ]] || { echo "invalid frozen ploidypatch-dev binding" >&2; exit 1; }
python_bin=$dev_prefix/bin/python
policy=$protocol_root/policy.tsv
seed=20260930
maximum_count=800
result_root=$project_root/benchmark/structure/copy_collapse_v0.4/ptr_v4/annotation_copy_collapse_seed${seed}
working_root=${result_root}.working
self_relative=scripts/run_populus_copy_collapse_benchmark_v0.4.sh

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

implementation_dependencies=(
    "$self_relative"
    scripts/audit_copy_pair_selection_truth.py
    src/ploidypatch/cli.py
    src/ploidypatch/copy_pair_sampling.py
    src/ploidypatch/perturb.py
    src/ploidypatch/structure_perturb.py
    src/ploidypatch/score.py
)
for required in "$python_bin" "$source_root/EVALUATOR_SHA256SUMS" "$source_gff" "$source_genome" \
                "$pair_root/SHA256SUMS" "$pair_tsv" "$protocol_root/SHA256SUMS" \
                "$execution_root/SHA256SUMS" "$execution_root/implementation_manifest.tsv" \
                "$policy" "${implementation_dependencies[@]/#/$code_root/}"; do
    [[ -s $required ]] || { echo "missing Populus benchmark prerequisite: $required" >&2; exit 1; }
done
(cd "$source_root" && sha256sum -c EVALUATOR_SHA256SUMS >/dev/null)
(cd "$pair_root" && sha256sum -c SHA256SUMS >/dev/null)
(cd "$protocol_root" && sha256sum -c SHA256SUMS >/dev/null)
(cd "$execution_root" && sha256sum -c SHA256SUMS >/dev/null)
for relative in "${implementation_dependencies[@]}"; do verify_implementation "$relative"; done
policy_value() {
    local key=$1
    awk -F '\t' -v key="$key" '$1 == key {print $2}' "$policy"
}
[[ $(policy_value policy_id) == ploidypatch_populus_external_validation_v0.4 \
    && $(policy_value truth_sampler_seed) == "$seed" \
    && $(policy_value truth_event_count) == 800_if_eligible_else_all_eligible_with_shortfall_reported \
    && $(policy_value minimum_formal_event_count) == 500 \
    && $(policy_value minimum_target_chromosomes) == 15 \
    && $(policy_value minimum_events_per_complexity_bin) == 20 ]] || {
    echo "Populus benchmark policy identity, seed, count, or evaluability gates differ" >&2; exit 1;
}
[[ ! -e $result_root && ! -e $working_root ]] || {
    echo "refusing to overwrite Populus v0.4 benchmark" >&2; exit 1;
}
mkdir -p "$working_root/blind" "$working_root/evaluator/pair_selection" \
    "$working_root/evaluator/truth" "$working_root/evaluator/validation"
record_invalid() {
    local status=$?
    if [[ -d $working_root ]]; then
        printf 'field\tvalue\nformal_status\tinvalid_run\nstage\thidden_event_benchmark\nexit_status\t%s\n' \
            "$status" > "$working_root/invalid_run.tsv" || true
    fi
    exit "$status"
}
trap record_invalid ERR
{
    printf 'field\tvalue\n'
    printf 'code_commit\t%s\n' "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'dataset_id\tptr_v4\nevent_type\tannotation_copy_collapse\n'
    printf 'requested_count\t%s\nseed\t%s\nsplit\tuntouched_external_v0.4\n' \
        "$maximum_count" "$seed"
    printf 'pair_access\tevaluator_only\npolicy_frozen_before_truth\ttrue\n'
    printf 'candidate_model_frozen_before_truth\ttrue\nautomatic_approval\tfalse\n'
    printf 'candidate_reference_access\tfalse\n'
    printf 'fixed_pair_rules_relaxed\tfalse\n'
    printf 'protocol_freeze_sha256sums_sha256\t%s\n' \
        "$(sha256sum "$protocol_root/SHA256SUMS" | awk '{print $1}')"
    printf 'execution_freeze_sha256sums_sha256\t%s\n' \
        "$(sha256sum "$execution_root/SHA256SUMS" | awk '{print $1}')"
    printf 'pair_tsv_sha256\t%s\ntarget_genome_sha256\t%s\n' \
        "$(sha256sum "$pair_tsv" | awk '{print $1}')" \
        "$(sha256sum "$source_genome" | awk '{print $1}')"
} > "$working_root/run_contract.tsv"

finalize_tree() {
    du -sb "$working_root" > "$working_root/disk_bytes.txt"
    (
        cd "$working_root"
        find . -type f \( -name '*.gff3' -o -name '*.fa' -o -name '*.json' \
            -o -name '*.tsv' -o -name SHA256SUMS -o -name disk_bytes.txt \) \
            ! -path './SHA256SUMS' -print0 \
            | sort -z | xargs -0 sha256sum > SHA256SUMS
        sha256sum -c SHA256SUMS >/dev/null
    )
    trap - ERR
    mv "$working_root" "$result_root"
}

available_pairs=$(( $(wc -l < "$pair_tsv") - 1 ))
[[ $available_pairs -ge 0 ]] || { echo "Populus truth-pair table lacks a header" >&2; exit 1; }
printf 'available_intersection_pairs\t%s\n' "$available_pairs" >> "$working_root/run_contract.tsv"
if [[ $available_pairs -eq 0 ]]; then
    "$python_bin" - "$working_root/evaluator/pair_selection/evaluability.json" <<'PY'
import json
import sys
from pathlib import Path

report = {
    "schema_version": "ploidypatch.populus_external_evaluability.v0.4",
    "events": 0,
    "target_chromosomes": 0,
    "events_by_target_chromosome": {},
    "complexity_bins": {"one": 0, "two_to_three": 0, "four_to_six": 0, "seven_plus": 0},
    "sentinels_executed": False,
    "sentinels": {
        "blind_noop_exact_recovery": None,
        "complete_oracle_exact_recovery": None,
        "restoration_byte_identical": None,
        "blind_complete_genome_sha256_identical": None,
    },
    "data_gates": {
        "minimum_events": False,
        "minimum_target_chromosomes": False,
        "four_complexity_bins_present": True,
        "minimum_events_each_complexity_bin": False,
    },
    "formal_evaluable": False,
    "formal_outcome": "not_evaluable_without_rule_relaxation",
    "reason": "zero_pairs_pass_frozen_self_and_two_outgroup_intersection",
    "fixed_rules_relaxed": False,
}
with Path(sys.argv[1]).open("x", encoding="utf-8") as handle:
    json.dump(report, handle, indent=2, sort_keys=True)
    handle.write("\n")
PY
    printf 'field\tvalue\nformal_status\tnot_evaluable_without_rule_relaxation\nstage\thidden_event_benchmark\n' \
        > "$working_root/stage_status.tsv"
    finalize_tree
    printf 'Populus v0.4 benchmark is not evaluable under frozen pair rules: %s\n' "$result_root"
    exit 0
fi

selected=$working_root/evaluator/pair_selection/selected_pairs.tsv
selection_decisions=$working_root/evaluator/pair_selection/decisions.tsv
cd "$code_root"
"$python_bin" -m ploidypatch.cli benchmark sample-copy-pairs \
    --source-gff "$source_gff" --pairs "$pair_tsv" \
    --count "$maximum_count" --seed "$seed" \
    --output-pairs "$selected" --decisions-tsv "$selection_decisions" \
    > "$working_root/evaluator/pair_selection/stdout.json" \
    2> "$working_root/evaluator/pair_selection/stderr.log"
for required in "$selected" "$selected.manifest.json" "$selection_decisions"; do
    [[ -s $required ]] || { echo "missing frozen Populus pair-sampling artifact: $required" >&2; exit 1; }
done
selected_count=$(( $(wc -l < "$selected") - 1 ))
[[ $selected_count -gt 0 && $selected_count -le $maximum_count ]] || {
    echo "Populus selected pair count is outside 1..$maximum_count" >&2; exit 1;
}
printf 'selected_count\t%s\n' "$selected_count" >> "$working_root/run_contract.tsv"

"$python_bin" -m ploidypatch.cli benchmark perturb \
    --gff "$source_gff" --output-dir "$working_root/blind" \
    --truth-dir "$working_root/evaluator/truth" \
    --event-type annotation_copy_collapse --pair-tsv "$selected" \
    --count "$selected_count" --seed "$seed" \
    > "$working_root/perturb.stdout.json" 2> "$working_root/perturb.stderr.log"
perturbed=$working_root/blind/perturbed.gff3
truth=$working_root/evaluator/truth/hidden_truth.json
restored=$working_root/evaluator/restored.gff3
[[ -s $perturbed && -s $truth && -s $working_root/blind/manifest.json ]] || {
    echo "Populus perturbation lacks blind or evaluator output" >&2; exit 1;
}
mv "$working_root/blind/manifest.json" "$working_root/evaluator/perturbation_manifest.json"
"$python_bin" - "$perturbed" "$source_genome" "$working_root/blind/blind_manifest.json" <<'PY'
import hashlib
import json
import sys
from pathlib import Path

perturbed, genome, output = map(Path, sys.argv[1:])
def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()
manifest = {
    "schema_version": "ploidypatch.blind_benchmark_input.v0.4",
    "truth_access": False,
    "complete_target_annotation_access": False,
    "perturbed_annotation": {"file_name": perturbed.name, "sha256": sha256(perturbed)},
    "target_genome": {"mount_role": "shared_target_genome", "sha256": sha256(genome)},
}
with output.open("x", encoding="utf-8") as handle:
    json.dump(manifest, handle, indent=2, sort_keys=True)
    handle.write("\n")
PY
(
    cd "$working_root/blind"
    sha256sum perturbed.gff3 blind_manifest.json > SHA256SUMS
    sha256sum -c SHA256SUMS >/dev/null
)

"$python_bin" scripts/audit_copy_pair_selection_truth.py \
    --selected-pairs "$selected" --truth "$truth" \
    --output-json "$working_root/evaluator/validation/pair_truth_audit.json"
"$python_bin" -m ploidypatch.cli benchmark restore \
    --perturbed-gff "$perturbed" --truth "$truth" --output-gff "$restored" \
    > "$working_root/evaluator/restoration_report.json" \
    2> "$working_root/evaluator/restoration.stderr.log"
if cmp -s "$source_gff" "$restored"; then restoration_identical=true; else restoration_identical=false; fi
for mode in noop oracle; do
    candidate=$perturbed; [[ $mode == oracle ]] && candidate=$source_gff
    "$python_bin" -m ploidypatch.cli benchmark score \
        --source-gff "$source_gff" --perturbed-gff "$perturbed" \
        --candidate-gff "$candidate" --truth "$truth" --include-event-details \
        > "$working_root/evaluator/validation/score_$mode.json" \
        2> "$working_root/evaluator/validation/score_$mode.stderr.log"
done
"$python_bin" - "$selected" \
    "$working_root/evaluator/validation/score_noop.json" \
    "$working_root/evaluator/validation/score_oracle.json" \
    "$working_root/evaluator/validation/pair_truth_audit.json" \
    "$source_genome" "$source_genome" "$restoration_identical" \
    "$(policy_value minimum_formal_event_count)" \
    "$(policy_value minimum_target_chromosomes)" \
    "$(policy_value minimum_events_per_complexity_bin)" \
    "$working_root/evaluator/pair_selection/evaluability.json" <<'PY'
import csv
import hashlib
import json
import sys
from collections import Counter
from pathlib import Path

(
    selected_path, noop_path, oracle_path, audit_path, blind_genome,
    complete_genome, restoration_value, minimum_events,
    minimum_chromosomes, minimum_per_bin, output,
) = sys.argv[1:]
def load(path: str):
    with open(path, encoding="utf-8") as handle:
        return json.load(handle)
def sha256(path: str) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()
with open(selected_path, encoding="utf-8", newline="") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))
chromosomes = Counter(row["target_seqid"] for row in rows)
complexities = Counter(row["coding_complexity_bin"] for row in rows)
expected_bins = ("one", "two_to_three", "four_to_six", "seven_plus")
noop, oracle, audit = load(noop_path), load(oracle_path), load(audit_path)
events = len(rows)
sentinels = {
    "blind_noop_exact_recovery": noop.get("event_recovery", {}).get("complete_cds_chain_recovery"),
    "complete_oracle_exact_recovery": oracle.get("event_recovery", {}).get("complete_cds_chain_recovery"),
    "restoration_byte_identical": restoration_value == "true",
    "blind_complete_genome_sha256_identical": sha256(blind_genome) == sha256(complete_genome),
    "noop_quality_grade_pass": noop.get("quality_gate", {}).get("grade") == "pass",
    "oracle_quality_grade_pass": oracle.get("quality_gate", {}).get("grade") == "pass",
    "pair_truth_audit_pass": audit.get("grade") == "pass",
}
sentinel_valid = (
    sentinels["blind_noop_exact_recovery"] == 0
    and sentinels["complete_oracle_exact_recovery"] == events
    and all(value is True for key, value in sentinels.items() if key not in {
        "blind_noop_exact_recovery", "complete_oracle_exact_recovery"
    })
)
complexity = {name: complexities.get(name, 0) for name in expected_bins}
data_gates = {
    "minimum_events": events >= int(minimum_events),
    "minimum_target_chromosomes": len(chromosomes) >= int(minimum_chromosomes),
    "four_complexity_bins_present": len(complexity) == 4,
    "minimum_events_each_complexity_bin": min(complexity.values()) >= int(minimum_per_bin),
}
formal_evaluable = sentinel_valid and all(data_gates.values())
if not sentinel_valid:
    outcome = "invalid_run"
elif not formal_evaluable:
    outcome = "not_evaluable_without_rule_relaxation"
else:
    outcome = "formally_evaluable_pending_blind_and_complete_control_reveal"
report = {
    "schema_version": "ploidypatch.populus_external_evaluability.v0.4",
    "events": events,
    "target_chromosomes": len(chromosomes),
    "events_by_target_chromosome": dict(sorted(chromosomes.items())),
    "complexity_bins": complexity,
    "sentinels_executed": True,
    "sentinels": sentinels,
    "data_gates": data_gates,
    "formal_evaluable": formal_evaluable,
    "formal_outcome": outcome,
    "fixed_rules_relaxed": False,
}
with Path(output).open("x", encoding="utf-8") as handle:
    json.dump(report, handle, indent=2, sort_keys=True)
    handle.write("\n")
PY
formal_outcome=$("$python_bin" - "$working_root/evaluator/pair_selection/evaluability.json" <<'PY'
import json
import sys
with open(sys.argv[1], encoding="utf-8") as handle:
    print(json.load(handle)["formal_outcome"])
PY
)
printf 'field\tvalue\nformal_status\t%s\nstage\thidden_event_benchmark\n' "$formal_outcome" \
    > "$working_root/stage_status.tsv"
finalize_tree
if [[ $formal_outcome == invalid_run ]]; then
    echo "Populus benchmark sentinel violation was frozen as invalid_run: $result_root" >&2
    exit 1
fi
printf 'Populus v0.4 copy-collapse benchmark frozen (%s): %s\n' "$formal_outcome" "$result_root"
