#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then echo "usage: $0 PROJECT_ROOT" >&2; exit 2; fi
project_root=$(realpath "$1")
code_root=$project_root/code
python_bin=$project_root/envs/ploidypatch-dev/bin/python
miniprot_bin=$project_root/envs/ploidypatch-baseline/bin/miniprot
input_root=$project_root/data/derived/holdout_inputs/maize_v2
protein_root=$project_root/data/public/maize_v2_holdout
target=$input_root/zea_mays/primary_chromosomes.genome.fa
result_root=$project_root/results/baselines/maize_v2/miniprot
working_root=${result_root}.working
policy=$code_root/config/maize_v2_zero_retuning_policy.tsv

for required in "$python_bin" "$miniprot_bin" "$target" "${target}.fai" \
                "$protein_root/sorghum_bicolor/protein.fa.gz" \
                "$protein_root/setaria_italica/protein.fa.gz" "$policy"; do
    [[ -s $required ]] || { echo "missing maize miniprot input: $required" >&2; exit 1; }
done
if [[ -e $result_root || -e $working_root ]]; then
    echo "refusing to overwrite maize miniprot" >&2; exit 1
fi
mkdir -p "$working_root"/{reference,index,raw,logs}
{
    printf 'field\tvalue\ncode_commit\t%s\n' \
        "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'candidate_truth_access\tfalse\ntarget_complete_annotation_access\tfalse\n'
    printf 'reference_species\tSorghum_bicolor,Setaria_italica\nthreads\t32\n'
    printf 'within_method_reference_vote_count\t1\n'
    printf 'policy_sha256\t%s\n' "$(sha256sum "$policy" | awk '{print $1}')"
} > "$working_root/run_contract.tsv"

cd "$code_root"
"$python_bin" -m ploidypatch.cli baseline prepare-proteins \
    --protein "sorghum=$protein_root/sorghum_bicolor/protein.fa.gz" \
    --protein "setaria=$protein_root/setaria_italica/protein.fa.gz" \
    --output-fasta "$working_root/reference/maize_outgroups.pep.fa" \
    --output-map "$working_root/reference/maize_outgroups.map.tsv" \
    > "$working_root/reference/prepare.stdout.json" \
    2> "$working_root/reference/prepare.stderr.log"
/usr/bin/time -v -o "$working_root/logs/index.time.txt" \
    "$miniprot_bin" -t32 -d "$working_root/index/zma.mpi" "$target" \
    > "$working_root/logs/index.stdout.log" 2> "$working_root/logs/index.stderr.log"
/usr/bin/time -v -o "$working_root/logs/projection.time.txt" \
    "$miniprot_bin" -I -t32 --gff-only -N4 --outn 4 --outs 0.8 --outc 0.5 \
    "$working_root/index/zma.mpi" "$working_root/reference/maize_outgroups.pep.fa" \
    > "$working_root/raw/miniprot.gff3" \
    2> "$working_root/logs/projection.stderr.log"
grep -q '^##gff-version' "$working_root/raw/miniprot.gff3"
for output in reference/maize_outgroups.pep.fa \
              reference/maize_outgroups.map.tsv raw/miniprot.gff3; do
    [[ -s $working_root/$output ]] || { echo "missing maize miniprot output: $output" >&2; exit 1; }
done
(
    cd "$working_root"
    find reference raw -type f \( -name '*.fa' -o -name '*.tsv' \
        -o -name '*.json' -o -name '*.gff3' \) -print0 \
        | sort -z | xargs -0 sha256sum > SHA256SUMS
    sha256sum -c SHA256SUMS >/dev/null
)
du -sb "$working_root" > "$working_root/disk_bytes.txt"
mv "$working_root" "$result_root"
printf 'maize miniprot upstream frozen: %s\n' "$result_root"
