#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then echo "usage: $0 PROJECT_ROOT" >&2; exit 2; fi
project_root=$(realpath "$1")
code_root=$project_root/code
python_bin=$project_root/envs/ploidypatch-dev/bin/python
gffread_bin=$project_root/envs/ploidypatch-syngap/bin/gffread
parallel_bin=/data/codexli/software/conda/miniforge3/bin/parallel
source_root=$project_root/data/derived/external_inputs/apple_v0.3
protocol_root=$project_root/results/protocol_freezes/apple_external_v0.3
result_root=$project_root/data/derived/external_evaluator/apple_v0.3_wgdi_inputs
working_root=${result_root}.working

declare -A prefix=(
    [target_apple]=mdx [evaluator_rose]=rch [evaluator_strawberry]=fve
)
species=(target_apple evaluator_rose evaluator_strawberry)
for required in "$python_bin" "$gffread_bin" "$parallel_bin" \
                "$source_root/SHA256SUMS" "$protocol_root/SHA256SUMS"; do
    [[ -s $required ]] || { echo "missing apple evaluator input prerequisite: $required" >&2; exit 1; }
done
(cd "$source_root" && sha256sum -c SHA256SUMS >/dev/null)
(cd "$protocol_root" && sha256sum -c SHA256SUMS >/dev/null)
for name in "${species[@]}"; do
    bundle=$source_root/$name
    for required in "$bundle/primary_chromosomes.gff3" \
                    "$bundle/primary_chromosomes.genome.fa" \
                    "$bundle/primary_chromosomes.genome.fa.fai"; do
        [[ -s $required ]] || { echo "missing apple evaluator source: $required" >&2; exit 1; }
    done
done
[[ ! -e $result_root && ! -e $working_root ]] || {
    echo "refusing to overwrite apple evaluator WGDI inputs" >&2; exit 1;
}
mkdir -p "$working_root"
{
    printf 'field\tvalue\n'
    printf 'code_commit\t%s\n' "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'access\tevaluator_only_after_protocol_freeze\n'
    printf 'candidate_reference_access\tfalse\n'
    printf 'hidden_pair_enumeration\tnot_in_this_stage\n'
    printf 'hidden_event_generation\tfalse\n'
    printf 'external_label_access\tfalse\n'
    printf 'protocol_freeze_sha256sums_sha256\t%s\n' \
        "$(sha256sum "$protocol_root/SHA256SUMS" | awk '{print $1}')"
} > "$working_root/run_contract.tsv"

prepare_one() {
    local name=$1 short=$2
    local bundle=$source_root/$name out=$working_root/$short
    local gff=$bundle/primary_chromosomes.gff3
    local genome=$bundle/primary_chromosomes.genome.fa
    local fai=$bundle/primary_chromosomes.genome.fa.fai
    mkdir -p "$out"
    /usr/bin/time -v -o "$out/gffread.resource.time.txt" \
        "$gffread_bin" "$gff" -g "$genome" -y "$out/$short.all.pep.fa" -S \
        > "$out/gffread.stdout.log" 2> "$out/gffread.stderr.log"
    /usr/bin/time -v -o "$out/prepare.resource.time.txt" \
        "$python_bin" -m ploidypatch.cli evidence prepare-wgdi \
        --gff "$gff" --protein "$out/$short.all.pep.fa" --fai "$fai" \
        --output-dir "$out" --prefix "$short" --min-genes-per-seqid 100 \
        > "$out/prepare.stdout.json" 2> "$out/prepare.stderr.log"
}
export -f prepare_one
export source_root working_root gffread_bin python_bin
printf '%s\t%s\n' target_apple mdx evaluator_rose rch evaluator_strawberry fve \
    | "$parallel_bin" --jobs 3 --delay 1 --colsep '\t' prepare_one {1} {2}

for short in mdx rch fve; do
    for suffix in wgdi.gff wgdi.lens wgdi.pep.fa wgdi_inputs.manifest.json; do
        [[ -s $working_root/$short/$short.$suffix ]] || {
            echo "missing apple evaluator WGDI artifact: $short.$suffix" >&2; exit 1;
        }
    done
done
(
    cd "$working_root"
    find . -type f \( -name '*.gff' -o -name '*.fa' -o -name '*.json' \
        -o -name '*.tsv' \) -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS
    sha256sum -c SHA256SUMS >/dev/null
)
du -sb "$working_root" > "$working_root/disk_bytes.txt"
mv "$working_root" "$result_root"
printf 'apple evaluator WGDI inputs frozen: %s\n' "$result_root"
