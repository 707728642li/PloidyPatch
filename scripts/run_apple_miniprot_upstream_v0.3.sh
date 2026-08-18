#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then echo "usage: $0 PROJECT_ROOT" >&2; exit 2; fi
project_root=$(realpath "$1")
code_root=$project_root/code
python_bin=$project_root/envs/ploidypatch-dev/bin/python
miniprot_bin=$project_root/envs/ploidypatch-baseline/bin/miniprot
source_root=$project_root/data/derived/external_inputs/apple_v0.3
public_root=$project_root/data/public/apple_external_v0.3
protocol_root=$project_root/results/protocol_freezes/apple_external_v0.3
target=$source_root/target_apple/primary_chromosomes.genome.fa
result_root=$project_root/results/baselines/apple_v0.3/miniprot
working_root=${result_root}.working

for required in "$python_bin" "$miniprot_bin" "$target" "${target}.fai" \
    "$public_root/candidate_pear/protein.fa.gz" \
    "$public_root/candidate_peach/protein.fa.gz" "$protocol_root/SHA256SUMS"; do
    [[ -s $required ]] || { echo "missing apple miniprot input: $required" >&2; exit 1; }
done
(cd "$protocol_root" && sha256sum -c SHA256SUMS >/dev/null)
[[ ! -e $result_root && ! -e $working_root ]] || {
    echo "refusing to overwrite apple miniprot upstream" >&2; exit 1;
}
mkdir -p "$working_root"/{reference,index,raw,logs}
{
    printf 'field\tvalue\ncode_commit\t%s\n' \
        "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'candidate_truth_access\tfalse\ntarget_complete_annotation_access\tfalse\n'
    printf 'reference_species\tPyrus_communis,Prunus_persica\nthreads\t32\n'
    printf 'within_method_reference_vote_count\t1\n'
    printf 'protocol_freeze_sha256sums_sha256\t%s\n' \
        "$(sha256sum "$protocol_root/SHA256SUMS" | awk '{print $1}')"
} > "$working_root/run_contract.tsv"

cd "$code_root"
"$python_bin" -m ploidypatch.cli baseline prepare-proteins \
    --protein "pear=$public_root/candidate_pear/protein.fa.gz" \
    --protein "peach=$public_root/candidate_peach/protein.fa.gz" \
    --output-fasta "$working_root/reference/apple_candidate_refs.pep.fa" \
    --output-map "$working_root/reference/apple_candidate_refs.map.tsv" \
    > "$working_root/reference/prepare.stdout.json" \
    2> "$working_root/reference/prepare.stderr.log"
/usr/bin/time -v -o "$working_root/logs/index.time.txt" \
    "$miniprot_bin" -t32 -d "$working_root/index/mdx.mpi" "$target" \
    > "$working_root/logs/index.stdout.log" 2> "$working_root/logs/index.stderr.log"
/usr/bin/time -v -o "$working_root/logs/projection.time.txt" \
    "$miniprot_bin" -I -t32 --gff-only -N4 --outn 4 --outs 0.8 --outc 0.5 \
    "$working_root/index/mdx.mpi" "$working_root/reference/apple_candidate_refs.pep.fa" \
    > "$working_root/raw/miniprot.gff3" \
    2> "$working_root/logs/projection.stderr.log"
grep -q '^##gff-version' "$working_root/raw/miniprot.gff3"
for output in reference/apple_candidate_refs.pep.fa \
              reference/apple_candidate_refs.map.tsv raw/miniprot.gff3; do
    [[ -s $working_root/$output ]] || { echo "missing apple miniprot output: $output" >&2; exit 1; }
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
printf 'apple miniprot upstream frozen: %s\n' "$result_root"
