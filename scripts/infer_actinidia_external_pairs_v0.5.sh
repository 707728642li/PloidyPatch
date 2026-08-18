
#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then echo "usage: $0 PROJECT_ROOT" >&2; exit 2; fi
project_root=$(realpath "$1")
input_root=$project_root/data/derived/external_evaluator/actinidia_v0.5_wgdi_inputs
normalized_root=$project_root/data/derived/external_inputs/actinidia/v0.5
wgdi_root=$project_root/results/evaluator/actinidia/v0.5/wgdi
protocol_root=$project_root/results/protocol_freezes/actinidia_external_v0.5
execution_root=$project_root/results/protocol_freezes/actinidia_external_v0.5_execution
code_root=$execution_root/source
export PYTHONPATH="$code_root/src${PYTHONPATH:+:$PYTHONPATH}"
environment_bindings=$execution_root/environment_bindings.tsv
[[ -s $environment_bindings ]] || { echo "missing frozen environment bindings" >&2; exit 1; }
dev_prefix=$(awk -F '\t' '$1 == "ploidypatch-dev" {print $2}' "$environment_bindings")
[[ $dev_prefix == /* ]] || { echo "invalid frozen ploidypatch-dev binding" >&2; exit 1; }
python_bin=$dev_prefix/bin/python
result_root=$project_root/results/evaluator/actinidia/v0.5/truth_pairs
working_root=${result_root}.working
query_gff=$input_root/red5/red5.wgdi.gff
representatives=$input_root/red5/red5.representatives.tsv
source_gff=$normalized_root/normalized/target_red5/primary_chromosomes.gff3
self_collinearity=$wgdi_root/collinearity/red5_self.tsv
rhododendron_collinearity=$wgdi_root/collinearity/red5_vs_rhs.tsv
diospyros_collinearity=$wgdi_root/collinearity/red5_vs_dol.tsv
self_relative=scripts/infer_actinidia_external_pairs_v0.5.sh

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
                "$rhododendron_collinearity" "$diospyros_collinearity" \
                "${implementation_dependencies[@]/#/$code_root/}"; do
    [[ -s $required ]] || { echo "missing Actinidia pair prerequisite: $required" >&2; exit 1; }
done
(cd "$input_root" && sha256sum -c SHA256SUMS >/dev/null)
(cd "$normalized_root" && sha256sum -c EVALUATOR_SHA256SUMS >/dev/null)
(cd "$wgdi_root" && sha256sum -c SHA256SUMS >/dev/null)
(cd "$protocol_root" && sha256sum -c SHA256SUMS >/dev/null)
(cd "$execution_root" && sha256sum -c SHA256SUMS >/dev/null)
for relative in "${implementation_dependencies[@]}"; do verify_implementation "$relative"; done
[[ ! -e $result_root && ! -e $working_root ]] || {
    echo "refusing to overwrite Actinidia evaluator truth pairs" >&2; exit 1;
}
mkdir -p "$working_root"/{mapping,self_wgdi,evaluator_groups/rhododendron,evaluator_groups/diospyros,two_evaluators,intersection}
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
    printf 'target_event\tActinidia_specific_Ad_alpha\n'
    printf 'final_rule\texact_unordered_pair_intersection_of_self_wgdi_and_two_evaluator_groups\n'
    printf 'min_block_pairs\t20\nevaluator_groups_required\t2\n'
    printf 'counterpart_target_multiplicity\texactly_two_no_truncation\n'
    printf 'per_group_pair_consistency\treject_if_either_target_has_any_other_exact_1_to_2_pair_in_that_group\n'
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
source_alias_gff=$working_root/mapping/red5.source_alias.gff3
"$python_bin" scripts/build_wgdi_source_alias_gff.py \
    --source-gff "$source_gff" --representatives "$representatives" \
    --output-gff "$source_alias_gff" \
    > "$working_root/mapping/stdout.json" 2> "$working_root/mapping/stderr.log"
"$python_bin" -m ploidypatch.cli evidence infer-self-wgd-pairs \
    --query-wgdi-gff "$query_gff" --collinearity "$self_collinearity" \
    --source-gff "$source_alias_gff" --wgd-event Actinidia_specific_Ad_alpha_self \
    --min-block-pairs 20 \
    --output-pairs "$working_root/self_wgdi/pairs.tsv" \
    --decisions-tsv "$working_root/self_wgdi/decisions.tsv" \
    > "$working_root/self_wgdi/stdout.json" 2> "$working_root/self_wgdi/stderr.log"
"$python_bin" -m ploidypatch.cli evidence infer-outgroup-duplicated-pairs \
    --query-wgdi-gff "$query_gff" --source-gff "$source_alias_gff" \
    --collinearity "rhododendron=$rhododendron_collinearity" \
    --wgd-event Actinidia_specific_Ad_alpha_rhododendron \
    --min-support-group-count 1 --min-block-pairs 20 \
    --output-pairs "$working_root/evaluator_groups/rhododendron/pairs.tsv" \
    --decisions-tsv "$working_root/evaluator_groups/rhododendron/decisions.tsv" \
    > "$working_root/evaluator_groups/rhododendron/stdout.json" \
    2> "$working_root/evaluator_groups/rhododendron/stderr.log"
"$python_bin" -m ploidypatch.cli evidence infer-outgroup-duplicated-pairs \
    --query-wgdi-gff "$query_gff" --source-gff "$source_alias_gff" \
    --collinearity "diospyros=$diospyros_collinearity" \
    --wgd-event Actinidia_specific_Ad_alpha_diospyros \
    --min-support-group-count 1 --min-block-pairs 20 \
    --output-pairs "$working_root/evaluator_groups/diospyros/pairs.tsv" \
    --decisions-tsv "$working_root/evaluator_groups/diospyros/decisions.tsv" \
    > "$working_root/evaluator_groups/diospyros/stdout.json" \
    2> "$working_root/evaluator_groups/diospyros/stderr.log"

# Do not rely on the v0.4 combined-group reciprocal gate: it only examines
# pairs already meeting the group threshold.  Each evaluator is first gated
# independently, so an alternate exact-1:2 pair in either group rejects every
# incident main pair, including a two-serial-WGD pattern seen by one group only.
"$python_bin" -m ploidypatch.cli evidence intersect-copy-pair-evidence \
    --pairs "rhododendron=$working_root/evaluator_groups/rhododendron/pairs.tsv" \
    --pairs "diospyros=$working_root/evaluator_groups/diospyros/pairs.tsv" \
    --pair-set-label Actinidia_specific_Ad_alpha_two_evaluator_groups \
    --output-pairs "$working_root/two_evaluators/pairs.tsv" \
    --decisions-tsv "$working_root/two_evaluators/decisions.tsv" \
    > "$working_root/two_evaluators/stdout.json" \
    2> "$working_root/two_evaluators/stderr.log"

"$python_bin" - \
    "$working_root/evaluator_groups/rhododendron/decisions.tsv" \
    "$working_root/evaluator_groups/diospyros/decisions.tsv" \
    "$working_root/two_evaluators/pairs.tsv" \
    "$working_root/pair_consistency_audit.tsv" \
    "$working_root/pair_consistency_audit.json" <<'PY'
import csv
import hashlib
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

rh_decisions, dol_decisions, accepted_path, audit_tsv, audit_json = map(Path, sys.argv[1:])

def read(path):
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not {"gene_id_a", "gene_id_b", "status", "reason"} <= set(reader.fieldnames or ()):
            raise ValueError(f"pair-consistency input schema is incomplete: {path}")
        return list(reader)

def pair(row):
    a, b = row["gene_id_a"], row["gene_id_b"]
    if not a or not b or a == b:
        raise ValueError("invalid exact unordered pair in evaluator decisions")
    return tuple(sorted((a, b)))

summaries = []
accepted_by_group = {}
all_pairs_by_group = {}
for label, path in (("rhododendron", rh_decisions), ("diospyros", dol_decisions)):
    rows = read(path)
    pairs = {pair(row) for row in rows}
    accepted = {pair(row) for row in rows if row["status"] == "accepted"}
    partners = defaultdict(set)
    for a, b in pairs:
        partners[a].add(b)
        partners[b].add(a)
    conflicts = {gene for gene, values in partners.items() if len(values) > 1}
    if any(a in conflicts or b in conflicts for a, b in accepted):
        raise ValueError(f"{label} accepted a pair incident to an alternate qualifying pair")
    reasons = Counter(row["reason"] for row in rows)
    summaries.append(
        {
            "evaluator_group": label,
            "qualifying_exact_1_to_2_pairs": len(pairs),
            "accepted_pair_consistent_pairs": len(accepted),
            "discordant_target_members": len(conflicts),
            "rejected_nonreciprocal_multiple_partners": reasons.get("nonreciprocal_multiple_partners", 0),
            "rejected_support_below_threshold": reasons.get("support_below_threshold", 0),
        }
    )
    accepted_by_group[label] = accepted
    all_pairs_by_group[label] = pairs

with accepted_path.open(encoding="utf-8", newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    observed = {pair(row) for row in reader}
expected = accepted_by_group["rhododendron"] & accepted_by_group["diospyros"]
if observed != expected:
    raise ValueError("two-evaluator output is not the exact unordered intersection of pair-consistent groups")

fields = tuple(summaries[0])
with audit_tsv.open("x", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(summaries)
payload = {
    "schema_version": "ploidypatch.actinidia_evaluator_pair_consistency.v0.5",
    "target_event": "Actinidia_specific_Ad_alpha",
    "rule": "reject_pair_if_either_target_member_has_any_alternate_exact_1_to_2_pair_in_either_evaluator_group",
    "evaluator_groups": summaries,
    "accepted_exact_two_group_intersection_pairs": len(observed),
    "fixed_rules_relaxed": False,
    "truth_labels_generated": False,
    "inputs": {
        path.name + ":" + label: hashlib.sha256(path.read_bytes()).hexdigest()
        for label, path in (("rhododendron", rh_decisions), ("diospyros", dol_decisions))
    },
}
with audit_json.open("x", encoding="utf-8") as handle:
    json.dump(payload, handle, indent=2, sort_keys=True)
    handle.write("\n")
PY
"$python_bin" -m ploidypatch.cli evidence intersect-copy-pair-evidence \
    --pairs "self_wgdi=$working_root/self_wgdi/pairs.tsv" \
    --pairs "two_evaluators=$working_root/two_evaluators/pairs.tsv" \
    --pair-set-label Actinidia_specific_Ad_alpha_self_and_two_evaluator_groups \
    --output-pairs "$working_root/intersection/pairs.tsv" \
    --decisions-tsv "$working_root/intersection/decisions.tsv" \
    > "$working_root/intersection/stdout.json" \
    2> "$working_root/intersection/stderr.log"

for path in self_wgdi/pairs.tsv evaluator_groups/rhododendron/pairs.tsv \
            evaluator_groups/diospyros/pairs.tsv two_evaluators/pairs.tsv \
            intersection/pairs.tsv pair_consistency_audit.tsv pair_consistency_audit.json; do
    [[ -s $working_root/$path ]] || { echo "missing Actinidia pair output: $path" >&2; exit 1; }
done
{
    printf 'pair_set\taccepted_pairs\n'
    for label in self_wgdi evaluator_groups/rhododendron \
                 evaluator_groups/diospyros two_evaluators intersection; do
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
    "schema_version": "ploidypatch.actinidia_external_pair_status.v0.5",
    "pair_counts": counts,
    "fixed_pair_rules_relaxed": False,
    "truth_labels_generated": False,
    "formal_status": (
        "pair_evidence_frozen_pending_event_sampling"
        if final > 0 else "not_evaluable_without_rule_relaxation"
    ),
    "reason": None if final > 0 else "zero_pairs_pass_frozen_self_and_two_evaluator_group_intersection",
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
printf 'Actinidia evaluator truth pairs frozen: %s\n' "$result_root"
