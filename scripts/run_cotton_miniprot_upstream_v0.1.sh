#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then echo "usage: $0 PROJECT_ROOT" >&2; exit 2; fi
project_root=$(realpath "$1")
code_root=$project_root/code
python_bin=$project_root/envs/ploidypatch-dev/bin/python
gffread_bin=$project_root/envs/ploidypatch-syngap/bin/gffread
miniprot_bin=$project_root/envs/ploidypatch-baseline/bin/miniprot
input_root=$project_root/data/derived/holdout_inputs/cotton_v0.1
target=$input_root/hirsutum/primary_chromosomes.genome.fa
result_root=$project_root/results/baselines/cotton_holdout_v0.1/miniprot
working_root=${result_root}.working
policy=$code_root/config/copy_collapse_zero_retuning_policy_v0.1.tsv

for required in "$python_bin" "$gffread_bin" "$miniprot_bin" "$target" "${target}.fai" \
                "$input_root/arboreum/primary_chromosomes.gff3" \
                "$input_root/arboreum/primary_chromosomes.genome.fa" \
                "$input_root/raimondii/primary_chromosomes.gff3" \
                "$input_root/raimondii/primary_chromosomes.genome.fa" "$policy"; do
    if [[ ! -s $required ]]; then echo "missing cotton miniprot input: $required" >&2; exit 1; fi
done
if [[ -e $result_root || -e $working_root ]]; then echo "refusing to overwrite cotton miniprot" >&2; exit 1; fi
mkdir -p "$working_root"/{reference,index,raw,logs}
{
    printf 'field\tvalue\ncode_commit\t%s\n' "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'candidate_truth_access\tfalse\ntarget_complete_annotation_access\tfalse\n'
    printf 'reference_species\tGossypium_arboreum,Gossypium_raimondii\nthreads\t64\n'
    printf 'policy_sha256\t%s\n' "$(sha256sum "$policy" | awk '{print $1}')"
} > "$working_root/run_contract.tsv"
for entry in gar_a:arboreum gra_d:raimondii; do
    label=${entry%%:*}; species=${entry#*:}
    "$gffread_bin" "$input_root/$species/primary_chromosomes.gff3" \
        -g "$input_root/$species/primary_chromosomes.genome.fa" \
        -y "$working_root/reference/$label.pep.fa" -S \
        > "$working_root/logs/$label.gffread.stdout.log" \
        2> "$working_root/logs/$label.gffread.stderr.log" &
done
wait
cd "$code_root"
"$python_bin" -m ploidypatch.cli baseline prepare-proteins \
    --protein "gar_a=$working_root/reference/gar_a.pep.fa" \
    --protein "gra_d=$working_root/reference/gra_d.pep.fa" \
    --output-fasta "$working_root/reference/cotton_diploids.pep.fa" \
    --output-map "$working_root/reference/cotton_diploids.map.tsv" \
    > "$working_root/reference/prepare.stdout.json" \
    2> "$working_root/reference/prepare.stderr.log"
/usr/bin/time -v -o "$working_root/logs/index.time.txt" \
    "$miniprot_bin" -t64 -d "$working_root/index/ghi_ad.mpi" "$target" \
    > "$working_root/logs/index.stdout.log" 2> "$working_root/logs/index.stderr.log"
/usr/bin/time -v -o "$working_root/logs/projection.time.txt" \
    "$miniprot_bin" -I -t64 --gff-only -N4 --outn 4 --outs 0.8 --outc 0.5 \
        "$working_root/index/ghi_ad.mpi" "$working_root/reference/cotton_diploids.pep.fa" \
        > "$working_root/raw/miniprot.gff3" 2> "$working_root/logs/projection.stderr.log"
grep -q '^##gff-version' "$working_root/raw/miniprot.gff3"
for output in reference/cotton_diploids.pep.fa reference/cotton_diploids.map.tsv raw/miniprot.gff3; do
    if [[ ! -s $working_root/$output ]]; then echo "missing cotton miniprot output: $output" >&2; exit 1; fi
done
(
    cd "$working_root"
    find reference raw -type f \( -name '*.fa' -o -name '*.tsv' -o -name '*.json' -o -name '*.gff3' \) \
        -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS
    sha256sum -c SHA256SUMS >/dev/null
)
du -sb "$working_root" > "$working_root/disk_bytes.txt"
mv "$working_root" "$result_root"
printf 'cotton miniprot upstream frozen: %s\n' "$result_root"
