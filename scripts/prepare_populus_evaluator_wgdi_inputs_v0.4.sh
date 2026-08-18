#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then echo "usage: $0 PROJECT_ROOT" >&2; exit 2; fi
project_root=$(realpath "$1")
parallel_bin=/data/codexli/software/conda/miniforge3/bin/parallel
source_root=$project_root/data/derived/external_inputs/populus_v0.4
protocol_root=$project_root/results/protocol_freezes/populus_external_v0.4
execution_root=$project_root/results/protocol_freezes/populus_external_v0.4_execution
code_root=$execution_root/source
export PYTHONPATH="$code_root/src${PYTHONPATH:+:$PYTHONPATH}"
environment_bindings=$execution_root/environment_bindings.tsv
[[ -s $environment_bindings ]] || { echo "missing frozen environment bindings" >&2; exit 1; }
dev_prefix=$(awk -F '\t' '$1 == "ploidypatch-dev" {print $2}' "$environment_bindings")
syngap_prefix=$(awk -F '\t' '$1 == "ploidypatch-syngap" {print $2}' "$environment_bindings")
[[ $dev_prefix == /* && $syngap_prefix == /* ]] || {
    echo "invalid frozen dev/syngap environment binding" >&2; exit 1;
}
python_bin=$dev_prefix/bin/python
gffread_bin=$syngap_prefix/bin/gffread
result_root=$project_root/data/derived/external_evaluator/populus_v0.4_wgdi_inputs
working_root=${result_root}.working
self_relative=scripts/prepare_populus_evaluator_wgdi_inputs_v0.4.sh

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
    src/ploidypatch/cli.py
    src/ploidypatch/synteny_io.py
)

for required in "$python_bin" "$gffread_bin" "$parallel_bin" \
                "$source_root/EVALUATOR_SHA256SUMS" "$protocol_root/SHA256SUMS" \
                "$execution_root/SHA256SUMS" "$execution_root/implementation_manifest.tsv" \
                "${implementation_dependencies[@]/#/$code_root/}"; do
    [[ -s $required ]] || { echo "missing Populus evaluator WGDI-input prerequisite: $required" >&2; exit 1; }
done
(cd "$source_root" && sha256sum -c EVALUATOR_SHA256SUMS >/dev/null)
(cd "$protocol_root" && sha256sum -c SHA256SUMS >/dev/null)
(cd "$execution_root" && sha256sum -c SHA256SUMS >/dev/null)
for relative in "${implementation_dependencies[@]}"; do verify_implementation "$relative"; done
declare -A prefix=(
    [target_populus]=ptr [evaluator_manihot]=mes [evaluator_ricinus]=rco
)
species=(target_populus evaluator_manihot evaluator_ricinus)
for name in "${species[@]}"; do
    bundle=$source_root/normalized/$name
    for required in "$bundle/primary_chromosomes.gff3" \
                    "$bundle/primary_chromosomes.genome.fa" \
                    "$bundle/primary_chromosomes.genome.fa.fai" "$bundle/manifest.json"; do
        [[ -s $required ]] || { echo "missing normalized Populus evaluator source: $required" >&2; exit 1; }
    done
done
[[ ! -e $result_root && ! -e $working_root ]] || {
    echo "refusing to overwrite Populus evaluator WGDI inputs" >&2; exit 1;
}
mkdir -p "$working_root"
record_invalid() {
    local status=$?
    if [[ -d $working_root ]]; then
        printf 'field\tvalue\nformal_status\tinvalid_run\nstage\tprepare_wgdi_inputs\nexit_status\t%s\n' \
            "$status" > "$working_root/invalid_run.tsv" || true
    fi
    exit "$status"
}
trap record_invalid ERR
{
    printf 'field\tvalue\n'
    printf 'code_commit\t%s\n' "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'access\tevaluator_only_after_protocol_and_execution_freeze\n'
    printf 'candidate_reference_access\tfalse\n'
    printf 'hidden_pair_enumeration\tfalse\n'
    printf 'hidden_event_generation\tfalse\n'
    printf 'truth_label_generation\tfalse\n'
    printf 'prefix_target\tptr\nprefix_evaluator_manihot\tmes\nprefix_evaluator_ricinus\trco\n'
    printf 'protocol_freeze_sha256sums_sha256\t%s\n' \
        "$(sha256sum "$protocol_root/SHA256SUMS" | awk '{print $1}')"
    printf 'execution_freeze_sha256sums_sha256\t%s\n' \
        "$(sha256sum "$execution_root/SHA256SUMS" | awk '{print $1}')"
} > "$working_root/run_contract.tsv"

prepare_one() {
    local name=$1 short=$2
    local bundle=$source_root/normalized/$name out=$working_root/$short
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
printf '%s\t%s\n' target_populus ptr evaluator_manihot mes evaluator_ricinus rco \
    | "$parallel_bin" --jobs 3 --delay 1 --colsep '\t' prepare_one {1} {2}

printf 'prefix\tgenes\tchromosomes\n' > "$working_root/wgdi_input_counts.tsv"
for short in ptr mes rco; do
    for suffix in wgdi.gff wgdi.lens wgdi.pep.fa wgdi_inputs.manifest.json representatives.tsv; do
        [[ -s $working_root/$short/$short.$suffix ]] || {
            echo "missing Populus evaluator WGDI artifact: $short.$suffix" >&2; exit 1;
        }
    done
    genes=$(wc -l < "$working_root/$short/$short.wgdi.gff")
    chromosomes=$(wc -l < "$working_root/$short/$short.wgdi.lens")
    [[ $genes -ge 10000 && $chromosomes -ge 10 ]] || {
        echo "$short WGDI input content floor failed" >&2; exit 1;
    }
    printf '%s\t%s\t%s\n' "$short" "$genes" "$chromosomes" >> "$working_root/wgdi_input_counts.tsv"
done
printf 'field\tvalue\nformal_status\tvalid_pretruth_artifact\nstage\tprepare_wgdi_inputs\n' \
    > "$working_root/stage_status.tsv"
du -sb "$working_root" > "$working_root/disk_bytes.txt"
(
    cd "$working_root"
    find . -type f \( -name '*.gff' -o -name '*.fa' -o -name '*.json' \
        -o -name '*.tsv' -o -name disk_bytes.txt \) ! -name SHA256SUMS -print0 \
        | sort -z | xargs -0 sha256sum > SHA256SUMS
    sha256sum -c SHA256SUMS >/dev/null
)
trap - ERR
mv "$working_root" "$result_root"
printf 'Populus evaluator WGDI inputs frozen: %s\n' "$result_root"
