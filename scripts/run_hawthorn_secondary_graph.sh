#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "usage: $0 PROJECT_ROOT" >&2
    exit 2
fi

project_root=$(realpath "$1")
code_root=$project_root/code
python_bin=$project_root/envs/ploidypatch-dev/bin/python
assembly_validation=$project_root/results/hawthorn/natural_assembly_context_v0.1/candidates.primary_rna.assembly.tsv
grouped_junctions=$project_root/results/hawthorn/rna_junctions_grouped_v0.1/all_samples.grouped_aggregate.tsv
prior_decisions=$project_root/results/hawthorn/natural_event_graph_v0.3/decisions.tsv
secondary_root=$project_root/results/hawthorn/natural_secondary_rna_v0.2
secondary_working=${secondary_root}.working
graph_root=$project_root/results/hawthorn/natural_event_graph_v0.4
graph_working=${graph_root}.working

for required in "$python_bin" "$assembly_validation" \
                "${assembly_validation}.manifest.json" "$grouped_junctions" \
                "${grouped_junctions}.manifest.json" "$prior_decisions"; do
    if [[ ! -s $required ]]; then
        echo "missing secondary-graph prerequisite: $required" >&2
        exit 1
    fi
done
for target in "$secondary_root" "$secondary_working" \
              "$graph_root" "$graph_working"; do
    if [[ -e $target ]]; then
        echo "refusing to overwrite secondary-graph artifact: $target" >&2
        exit 1
    fi
done

mkdir -p "$secondary_working"
{
    printf 'field\tvalue\n'
    printf 'code_commit\t%s\n' "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'assembly_context_preserved\ttrue\n'
    printf 'secondary_group_edge_direction\tcontext\n'
    printf 'secondary_group_support_changes_score\tfalse\n'
    printf 'negative_evidence_policy\tabsence_is_missing_not_contradiction\n'
    printf 'automatic_patch_policy\treview_required\n'
} > "$secondary_working/run_contract.tsv"
cd "$code_root"
/usr/bin/time -v -o "$secondary_working/resource.time.txt" \
    "$python_bin" -m ploidypatch.cli evidence validate-natural-secondary-rna \
        --primary-validation "$assembly_validation" \
        --grouped-junctions "$grouped_junctions" \
        --output-tsv "$secondary_working/candidates.assembly.secondary_rna.tsv" \
        > "$secondary_working/stdout.log" \
        2> "$secondary_working/stderr.log"
for output in \
    "$secondary_working/candidates.assembly.secondary_rna.tsv" \
    "$secondary_working/candidates.assembly.secondary_rna.tsv.manifest.json"; do
    if [[ ! -s $output ]]; then
        echo "secondary assembly/RNA output is missing or empty: $output" >&2
        exit 1
    fi
done
(
    cd "$secondary_working"
    sha256sum ./*.tsv ./*.json > SHA256SUMS
)
du -sb "$secondary_working" > "$secondary_working/disk_bytes.txt"
mv "$secondary_working" "$secondary_root"

mkdir -p "$graph_working"
cp "$secondary_root/run_contract.tsv" "$graph_working/run_contract.tsv"
"$python_bin" -m ploidypatch.cli evidence prepare-natural-graph \
    --validation-tsv "$secondary_root/candidates.assembly.secondary_rna.tsv" \
    --output-candidates "$graph_working/candidates.tsv" \
    --output-evidence "$graph_working/evidence.tsv" \
    > "$graph_working/adapter.stdout.log" \
    2> "$graph_working/adapter.stderr.log"
/usr/bin/time -v -o "$graph_working/resource.time.txt" \
    "$python_bin" -m ploidypatch.cli graph infer \
        --candidates "$graph_working/candidates.tsv" \
        --evidence "$graph_working/evidence.tsv" \
        --output-json "$graph_working/event_graph.json" \
        --decisions-tsv "$graph_working/decisions.tsv" \
        > "$graph_working/infer.stdout.log" \
        2> "$graph_working/infer.stderr.log"

{
    printf 'candidate_id\tprior_decision\tsecondary_context_decision\tchanged\n'
    join -t $'\t' -1 1 -2 1 \
        <(tail -n +2 "$prior_decisions" | cut -f2,11 | sort -k1,1) \
        <(tail -n +2 "$graph_working/decisions.tsv" | cut -f2,11 | sort -k1,1) \
        | awk -F '\t' 'BEGIN { OFS="\t" } { print $1, $2, $3, ($2 == $3 ? "false" : "true") }'
} > "$graph_working/decision_transition.tsv"
prior_count=$(tail -n +2 "$prior_decisions" | wc -l)
new_count=$(tail -n +2 "$graph_working/decisions.tsv" | wc -l)
transition_count=$(tail -n +2 "$graph_working/decision_transition.tsv" | wc -l)
decision_changes=$(awk -F '\t' '$4 == "true" { count++ } END { print count + 0 }' \
    "$graph_working/decision_transition.tsv")
if [[ $prior_count -ne $new_count || $prior_count -ne $transition_count ]]; then
    echo "secondary graph candidate set differs from v0.3" >&2
    exit 1
fi
if [[ $decision_changes -ne 0 ]]; then
    echo "context-only secondary RNA unexpectedly changed decisions" >&2
    exit 1
fi
{
    printf 'metric\tvalue\n'
    printf 'candidates_compared\t%s\n' "$transition_count"
    printf 'decision_changes\t%s\n' "$decision_changes"
} > "$graph_working/decision_transition.summary.tsv"

for output in \
    "$graph_working/candidates.tsv" \
    "$graph_working/candidates.tsv.manifest.json" \
    "$graph_working/evidence.tsv" \
    "$graph_working/event_graph.json" \
    "$graph_working/decisions.tsv" \
    "$graph_working/decision_transition.tsv" \
    "$graph_working/decision_transition.summary.tsv"; do
    if [[ ! -s $output ]]; then
        echo "secondary natural graph output is missing or empty: $output" >&2
        exit 1
    fi
done
(
    cd "$graph_working"
    sha256sum ./*.tsv ./*.json > SHA256SUMS
)
du -sb "$graph_working" > "$graph_working/disk_bytes.txt"
mv "$graph_working" "$graph_root"
printf 'hawthorn secondary context graph completed: %s\n' "$graph_root"
