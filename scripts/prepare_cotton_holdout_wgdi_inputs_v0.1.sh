#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then echo "usage: $0 PROJECT_ROOT" >&2; exit 2; fi
project_root=$(realpath "$1")
code_root=$project_root/code
python_bin=$project_root/envs/ploidypatch-dev/bin/python
gffread_bin=$project_root/envs/ploidypatch-syngap/bin/gffread
parallel_bin=/data/codexli/software/conda/miniforge3/bin/parallel
source_root=$project_root/data/derived/holdout_inputs/cotton_v0.1
policy=$code_root/config/copy_collapse_zero_retuning_policy_v0.1.tsv
result_root=$project_root/data/derived/holdout_evaluator/cotton_wgdi_inputs_v0.1
working_root=${result_root}.working

for required in "$python_bin" "$gffread_bin" "$parallel_bin" "$policy" \
                "$source_root/SHA256SUMS"; do
    if [[ ! -s $required ]]; then echo "missing cotton WGDI prerequisite: $required" >&2; exit 1; fi
done
(cd "$source_root" && sha256sum -c SHA256SUMS >/dev/null)
if [[ -e $result_root || -e $working_root ]]; then
    echo "refusing to overwrite cotton WGDI inputs" >&2; exit 1
fi
mkdir -p "$working_root"
{
    printf 'field\tvalue\n'
    printf 'code_commit\t%s\n' "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'access\tevaluator_only_after_rule_freeze\n'
    printf 'hidden_truth_access\tfalse\n'
    printf 'homeolog_pair_enumeration\tnot_in_this_stage\n'
    printf 'policy_sha256\t%s\n' "$(sha256sum "$policy" | awk '{print $1}')"
} > "$working_root/run_contract.tsv"

prepare_one() {
    local species=$1 prefix=$2
    local bundle=$source_root/$species out=$working_root/$prefix
    local gff=$bundle/primary_chromosomes.gff3
    local genome=$bundle/primary_chromosomes.genome.fa
    local fai=$bundle/primary_chromosomes.genome.fa.fai
    mkdir -p "$out"
    /usr/bin/time -v -o "$out/gffread.resource.time.txt" \
        "$gffread_bin" "$gff" -g "$genome" -y "$out/$prefix.all.pep.fa" -S \
        > "$out/gffread.stdout.log" 2> "$out/gffread.stderr.log"
    /usr/bin/time -v -o "$out/prepare.resource.time.txt" \
        "$python_bin" -m ploidypatch.cli evidence prepare-wgdi \
            --gff "$gff" --protein "$out/$prefix.all.pep.fa" --fai "$fai" \
            --output-dir "$out" --prefix "$prefix" --min-genes-per-seqid 100 \
            > "$out/prepare.stdout.json" 2> "$out/prepare.stderr.log"
}
export -f prepare_one
export source_root working_root gffread_bin python_bin
printf '%s\t%s\n' hirsutum ghi_ad arboreum gar_a raimondii gra_d \
    | "$parallel_bin" --jobs 3 --delay 1 --colsep '\t' prepare_one {1} {2}

for prefix in ghi_ad gar_a gra_d; do
    for suffix in wgdi.gff wgdi.lens wgdi.pep.fa wgdi_inputs.manifest.json; do
        if [[ ! -s $working_root/$prefix/$prefix.$suffix ]]; then
            echo "missing cotton WGDI input: $prefix.$suffix" >&2; exit 1
        fi
    done
done
(
    cd "$working_root"
    find . -type f \( -name '*.gff' -o -name '*.fa' -o -name '*.json' -o -name '*.tsv' \) \
        -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS
    sha256sum -c SHA256SUMS >/dev/null
)
du -sb "$working_root" > "$working_root/disk_bytes.txt"
mv "$working_root" "$result_root"
printf 'cotton evaluator WGDI inputs frozen: %s\n' "$result_root"
