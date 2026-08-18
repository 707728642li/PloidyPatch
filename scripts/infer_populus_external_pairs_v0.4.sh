#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then echo "usage: $0 PROJECT_ROOT" >&2; exit 2; fi
project_root=$(realpath "$1")
input_root=$project_root/data/derived/external_evaluator/populus_v0.4_wgdi_inputs
normalized_root=$project_root/data/derived/external_inputs/populus_v0.4
wgdi_root=$project_root/results/evaluator/populus/v0.4/wgdi
protocol_root=$project_root/results/protocol_freezes/populus_external_v0.4
execution_root=$project_root/results/protocol_freezes/populus_external_v0.4_execution
code_root=$execution_root/source
export PYTHONPATH="$code_root/src${PYTHONPATH:+:$PYTHONPATH}"
environment_bindings=$execution_root/environment_bindings.tsv
[[ -s $environment_bindings ]] || { echo "missing frozen environment bindings" >&2; exit 1; }
dev_prefix=$(awk -F '\t' '$1 == "ploidypatch-dev" {print $2}' "$environment_bindings")
[[ $dev_prefix == /* ]] || { echo "invalid frozen ploidypatch-dev binding" >&2; exit 1; }
python_bin=$dev_prefix/bin/python
result_root=$project_root/results/evaluator/populus/v0.4/truth_pairs
working_root=${result_root}.working
query_gff=$input_root/ptr/ptr.wgdi.gff
representatives=$input_root/ptr/ptr.representatives.tsv
source_gff=$normalized_root/normalized/target_populus/primary_chromosomes.gff3
self_collinearity=$wgdi_root/collinearity/ptr_self.tsv
manihot_collinearity=$wgdi_root/collinearity/ptr_vs_mes.tsv
ricinus_collinearity=$wgdi_root/collinearity/ptr_vs_rco.tsv
self_relative=scripts/infer_populus_external_pairs_v0.4.sh

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
    scripts/build_wgdi_source_alias_gff.py
    src/ploidypatch/cli.py
    src/ploidypatch/pair_consensus.py
)
for required in "$python_bin" "$input_root/SHA256SUMS" "$normalized_root/EVALUATOR_SHA256SUMS" \
                "$wgdi_root/SHA256SUMS" "$protocol_root/SHA256SUMS" \
                "$execution_root/SHA256SUMS" "$execution_root/implementation_manifest.tsv" \
                "$query_gff" "$representatives" "$source_gff" "$self_collinearity" \
                "$manihot_collinearity" "$ricinus_collinearity" \
                "${implementation_dependencies[@]/#/$code_root/}"; do
    [[ -s $required ]] || { echo "missing Populus pair prerequisite: $required" >&2; exit 1; }
done
(cd "$input_root" && sha256sum -c SHA256SUMS >/dev/null)
(cd "$normalized_root" && sha256sum -c EVALUATOR_SHA256SUMS >/dev/null)
(cd "$wgdi_root" && sha256sum -c SHA256SUMS >/dev/null)
(cd "$protocol_root" && sha256sum -c SHA256SUMS >/dev/null)
(cd "$execution_root" && sha256sum -c SHA256SUMS >/dev/null)
for relative in "${implementation_dependencies[@]}"; do verify_implementation "$relative"; done
[[ ! -e $result_root && ! -e $working_root ]] || {
    echo "refusing to overwrite Populus evaluator truth pairs" >&2; exit 1;
}
mkdir -p "$working_root"/{mapping,self_wgdi,two_outgroups,intersection}
record_invalid() {
    local status=$?
    if [[ -d $working_root ]]; then
        printf 'field\tvalue\nformal_status\tinvalid_run\nstage\ttruth_pair_inference\nexit_status\t%s\n' \
            "$status" > "$working_root/invalid_run.tsv" || true
    fi
    exit "$status"
}
trap record_invalid ERR
{
    printf 'field\tvalue\n'
    printf 'code_commit\t%s\n' "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'access\tevaluator_only_after_protocol_and_execution_freeze\n'
    printf 'candidate_reference_access\tfalse\ntruth_label_generation\tfalse\n'
    printf 'hidden_event_generation\tfalse\n'
    printf 'final_rule\texact_unordered_pair_intersection_of_self_wgdi_and_two_outgroup_support\n'
    printf 'min_block_pairs\t20\noutgroup_min_support_groups\t2\n'
    printf 'require_cross_seqid\ttrue\nrequire_reciprocal_unique\ttrue\n'
    printf 'identifier_mapping\texact_unique_alias_only\n'
    printf 'protocol_freeze_sha256sums_sha256\t%s\n' \
        "$(sha256sum "$protocol_root/SHA256SUMS" | awk '{print $1}')"
    printf 'execution_freeze_sha256sums_sha256\t%s\n' \
        "$(sha256sum "$execution_root/SHA256SUMS" | awk '{print $1}')"
    printf 'pair_intersection_module_sha256\t%s\n' \
        "$(sha256sum "$code_root/src/ploidypatch/pair_consensus.py" | awk '{print $1}')"
    printf 'identifier_alias_adapter_sha256\t%s\n' \
        "$(sha256sum "$code_root/scripts/build_wgdi_source_alias_gff.py" | awk '{print $1}')"
} > "$working_root/run_contract.tsv"

cd "$code_root"
source_alias_gff=$working_root/mapping/ptr.source_alias.gff3
"$python_bin" scripts/build_wgdi_source_alias_gff.py \
    --source-gff "$source_gff" --representatives "$representatives" \
    --output-gff "$source_alias_gff" \
    > "$working_root/mapping/stdout.json" 2> "$working_root/mapping/stderr.log"
"$python_bin" -m ploidypatch.cli evidence infer-self-wgd-pairs \
    --query-wgdi-gff "$query_gff" --collinearity "$self_collinearity" \
    --source-gff "$source_alias_gff" --wgd-event salicoid_wgd_self \
    --min-block-pairs 20 \
    --output-pairs "$working_root/self_wgdi/pairs.tsv" \
    --decisions-tsv "$working_root/self_wgdi/decisions.tsv" \
    > "$working_root/self_wgdi/stdout.json" 2> "$working_root/self_wgdi/stderr.log"
"$python_bin" -m ploidypatch.cli evidence infer-outgroup-duplicated-pairs \
    --query-wgdi-gff "$query_gff" --source-gff "$source_alias_gff" \
    --collinearity "manihot=$manihot_collinearity" \
    --collinearity "ricinus=$ricinus_collinearity" \
    --wgd-event salicoid_wgd_two_euphorbiaceae_outgroups \
    --min-support-group-count 2 --min-block-pairs 20 \
    --output-pairs "$working_root/two_outgroups/pairs.tsv" \
    --decisions-tsv "$working_root/two_outgroups/decisions.tsv" \
    > "$working_root/two_outgroups/stdout.json" \
    2> "$working_root/two_outgroups/stderr.log"
"$python_bin" -m ploidypatch.cli evidence intersect-copy-pair-evidence \
    --pairs "self_wgdi=$working_root/self_wgdi/pairs.tsv" \
    --pairs "two_outgroups=$working_root/two_outgroups/pairs.tsv" \
    --pair-set-label salicoid_wgd_self_and_two_euphorbiaceae_outgroups \
    --output-pairs "$working_root/intersection/pairs.tsv" \
    --decisions-tsv "$working_root/intersection/decisions.tsv" \
    > "$working_root/intersection/stdout.json" \
    2> "$working_root/intersection/stderr.log"

for path in self_wgdi/pairs.tsv two_outgroups/pairs.tsv intersection/pairs.tsv; do
    [[ -s $working_root/$path ]] || { echo "missing Populus pair output: $path" >&2; exit 1; }
done
{
    printf 'pair_set\taccepted_pairs\n'
    for label in self_wgdi two_outgroups intersection; do
        count=$(( $(wc -l < "$working_root/$label/pairs.tsv") - 1 ))
        [[ $count -ge 0 ]] || { echo "pair table lacks header: $label" >&2; exit 1; }
        printf '%s\t%s\n' "$label" "$count"
    done
} > "$working_root/pair_counts.tsv"
"$python_bin" - "$working_root/pair_counts.tsv" "$working_root/pair_status.json" <<'PY'
import csv
import json
import sys
from pathlib import Path

counts_path, output = map(Path, sys.argv[1:])
with counts_path.open(encoding="utf-8", newline="") as handle:
    counts = {row["pair_set"]: int(row["accepted_pairs"]) for row in csv.DictReader(handle, delimiter="\t")}
final = counts["intersection"]
status = {
    "schema_version": "ploidypatch.populus_external_pair_status.v0.4",
    "pair_counts": counts,
    "fixed_pair_rules_relaxed": False,
    "truth_labels_generated": False,
    "formal_status": (
        "pair_evidence_frozen_pending_event_sampling"
        if final > 0 else "not_evaluable_without_rule_relaxation"
    ),
    "reason": None if final > 0 else "zero_pairs_pass_frozen_self_and_two_outgroup_intersection",
}
with output.open("x", encoding="utf-8") as handle:
    json.dump(status, handle, indent=2, sort_keys=True)
    handle.write("\n")
PY
printf 'field\tvalue\nformal_status\tpair_evidence_frozen\nstage\ttruth_pair_inference\n' \
    > "$working_root/stage_status.tsv"
du -sb "$working_root" > "$working_root/disk_bytes.txt"
(
    cd "$working_root"
    find . -type f \( -name '*.tsv' -o -name '*.json' -o -name '*.gff3' \
        -o -name disk_bytes.txt \) \
        ! -name SHA256SUMS -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS
    sha256sum -c SHA256SUMS >/dev/null
)
trap - ERR
mv "$working_root" "$result_root"
printf 'Populus evaluator truth pairs frozen: %s\n' "$result_root"
