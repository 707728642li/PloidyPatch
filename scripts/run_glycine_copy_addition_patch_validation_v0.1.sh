#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "usage: $0 PROJECT_ROOT" >&2
    exit 2
fi

project_root=$(realpath "$1")
code_root=$project_root/code
python_bin=$project_root/envs/ploidypatch-dev/bin/python
benchmark_root=$project_root/benchmark/structure/copy_collapse_v0.1/gma_v21/annotation_copy_collapse_seed20260815
source_gff=$project_root/data/derived/syngap_glycine_v0.1/gma_complete/primary_chromosomes.gff3
blind_gff=$benchmark_root/blind/perturbed.gff3
truth=$benchmark_root/evaluator/truth/hidden_truth.json
wgd_root=$project_root/results/copy_collapse/wgd_reanchor_glycine_v0.1/self_wgdi
review_candidate=$wgd_root/blind/selected/candidate.gff3
control_candidate=$wgd_root/complete_control/selected/candidate.gff3
result_root=$project_root/results/copy_collapse/copy_addition_patch_validation_glycine_dev_v0.1
working_root=${result_root}.working

for required in "$python_bin" "$source_gff" "$blind_gff" "$truth" \
                "$review_candidate" "$control_candidate"; do
    if [[ ! -s $required ]]; then
        echo "missing or empty copy-addition patch input: $required" >&2
        exit 1
    fi
done
if [[ -e $result_root || -e $working_root ]]; then
    echo "refusing to overwrite copy-addition patch validation: $result_root" >&2
    exit 1
fi
mkdir -p "$working_root"

{
    printf 'field\tvalue\n'
    printf 'code_commit\t%s\n' "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'split\tpost_holdout_engineering_validation\n'
    printf 'formal_holdout_claim_allowed\tfalse\n'
    printf 'automatic_approval\tfalse\n'
    printf 'candidate_set\tblind_support2_consensus_with_mode_specific_self_wgd_reanchor\n'
    printf 'patch_operation\tEOF_only_copy_addition\n'
    printf 'existing_feature_modifications\t0\n'
    printf 'existing_feature_deletions\t0\n'
    printf 'copy_patch_module_sha256\t%s\n' \
        "$(sha256sum "$code_root/src/ploidypatch/copy_patch.py" | awk '{print $1}')"
    printf 'generic_patch_module_sha256\t%s\n' \
        "$(sha256sum "$code_root/src/ploidypatch/patch.py" | awk '{print $1}')"
} > "$working_root/run_contract.tsv"
{
    printf 'role\tbytes\tsha256\tpath\n'
    for entry in \
        "blind_gff:$blind_gff" \
        "review_candidate:$review_candidate" \
        "control_candidate:$control_candidate" \
        "source_gff:$source_gff"; do
        role=${entry%%:*}
        path=${entry#*:}
        printf '%s\t%s\t%s\t%s\n' "$role" "$(stat -Lc %s "$path")" \
            "$(sha256sum "$path" | awk '{print $1}')" "$path"
    done
} > "$working_root/candidate_input_manifest.tsv"

cd "$code_root"
"$python_bin" -m ploidypatch.cli patch compile-copy-additions \
    --annotation-gff "$blind_gff" \
    --candidate-gff "$review_candidate" \
    --output-edits-json "$working_root/copy_additions.edits.json" \
    > "$working_root/compile.stdout.json" \
    2> "$working_root/compile.stderr.log"
"$python_bin" -m ploidypatch.cli patch create \
    --source-gff "$blind_gff" \
    --edits-json "$working_root/copy_additions.edits.json" \
    --output-patch "$working_root/copy_additions.patch.json" \
    > "$working_root/create.stdout.json" \
    2> "$working_root/create.stderr.log"
"$python_bin" -m ploidypatch.cli patch apply \
    --source-gff "$blind_gff" \
    --patch "$working_root/copy_additions.patch.json" \
    --output-gff "$working_root/patched.gff3" \
    > "$working_root/apply.stdout.json" \
    2> "$working_root/apply.stderr.log"
if ! cmp -s "$working_root/patched.gff3" "$review_candidate"; then
    echo "applied copy-addition patch differs from reviewed candidate GFF" >&2
    exit 1
fi
"$python_bin" -m ploidypatch.cli patch revert \
    --patched-gff "$working_root/patched.gff3" \
    --patch "$working_root/copy_additions.patch.json" \
    --output-gff "$working_root/reverted.gff3" \
    > "$working_root/revert.stdout.json" \
    2> "$working_root/revert.stderr.log"
if ! cmp -s "$working_root/reverted.gff3" "$blind_gff"; then
    echo "reverted copy-addition patch is not byte-identical to blind GFF" >&2
    exit 1
fi

{
    printf 'role\tbytes\tsha256\tpath\n'
    for entry in \
        "edits:$working_root/copy_additions.edits.json" \
        "patch:$working_root/copy_additions.patch.json" \
        "patched:$working_root/patched.gff3" \
        "reverted:$working_root/reverted.gff3"; do
        role=${entry%%:*}
        path=${entry#*:}
        printf '%s\t%s\t%s\t%s\n' "$role" "$(stat -Lc %s "$path")" \
            "$(sha256sum "$path" | awk '{print $1}')" "$path"
    done
} > "$working_root/patch_freeze.tsv"
{
    printf 'role\tbytes\tsha256\tpath\n'
    printf 'hidden_truth\t%s\t%s\t%s\n' "$(stat -Lc %s "$truth")" \
        "$(sha256sum "$truth" | awk '{print $1}')" "$truth"
} > "$working_root/evaluator_input_manifest.tsv"

/usr/bin/time -v -o "$working_root/score.resource.time.txt" \
    "$python_bin" -m ploidypatch.cli benchmark score \
        --source-gff "$source_gff" \
        --perturbed-gff "$blind_gff" \
        --candidate-gff "$working_root/patched.gff3" \
        --control-candidate-gff "$control_candidate" \
        --truth "$truth" \
        --include-event-details \
        > "$working_root/score.json" \
        2> "$working_root/score.stderr.log"
if ! grep -q '"grade": "pass"' "$working_root/score.json"; then
    echo "copy-addition patch score quality gate failed" >&2
    exit 1
fi
(
    cd "$working_root"
    find . -type f \( -name '*.gff3' -o -name '*.json' -o -name '*.tsv' \) \
        -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS
)
du -sb "$working_root" > "$working_root/disk_bytes.txt"
mv "$working_root" "$result_root"
printf 'Glycine copy-addition patch validation frozen: %s\n' "$result_root"
