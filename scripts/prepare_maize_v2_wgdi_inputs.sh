#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then echo "usage: $0 PROJECT_ROOT" >&2; exit 2; fi
project_root=$(realpath "$1")
code_root=$project_root/code
python_bin=$project_root/envs/ploidypatch-dev/bin/python
parallel_bin=/data/codexli/software/conda/miniforge3/bin/parallel
source_root=$project_root/data/derived/holdout_inputs/maize_v2
protein_root=$project_root/data/public/maize_v2_holdout
policy=$code_root/config/maize_v2_zero_retuning_policy.tsv
result_root=$project_root/data/derived/holdout_evaluator/maize_v2_wgdi_inputs
working_root=${result_root}.working

declare -A prefix=(
    [zea_mays]=zma [sorghum_bicolor]=sbi [setaria_italica]=sit
)
species=(zea_mays sorghum_bicolor setaria_italica)
for required in "$python_bin" "$parallel_bin" "$policy" "$source_root/SHA256SUMS"; do
    [[ -s $required ]] || { echo "missing maize WGDI prerequisite: $required" >&2; exit 1; }
done
(cd "$source_root" && sha256sum -c SHA256SUMS >/dev/null)
for name in "${species[@]}"; do
    for required in "$source_root/$name/primary_chromosomes.gff3" \
                    "$source_root/$name/primary_chromosomes.genome.fa.fai" \
                    "$protein_root/$name/protein.fa.gz"; do
        [[ -s $required ]] || { echo "missing WGDI input: $required" >&2; exit 1; }
    done
done
if [[ -e $result_root || -e $working_root ]]; then
    echo "refusing to overwrite maize evaluator WGDI inputs" >&2; exit 1
fi
mkdir -p "$working_root"
{
    printf 'field\tvalue\n'
    printf 'code_commit\t%s\n' "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'access\tevaluator_only_after_protocol_freeze\n'
    printf 'hidden_truth_access\tfalse\n'
    printf 'homeolog_pair_enumeration\tnot_in_this_stage\n'
    printf 'policy_sha256\t%s\n' "$(sha256sum "$policy" | awk '{print $1}')"
} > "$working_root/run_contract.tsv"
{
    printf 'species\trole\tbytes\tsha256\tpath\n'
    for name in "${species[@]}"; do
        for role in gff fai protein; do
            case $role in
                gff) path=$source_root/$name/primary_chromosomes.gff3 ;;
                fai) path=$source_root/$name/primary_chromosomes.genome.fa.fai ;;
                protein) path=$protein_root/$name/protein.fa.gz ;;
            esac
            printf '%s\t%s\t%s\t%s\t%s\n' "$name" "$role" \
                "$(stat -Lc %s "$path")" "$(sha256sum "$path" | awk '{print $1}')" \
                "$path"
        done
    done
} > "$working_root/input_manifest.tsv"

prepare_one() {
    local name=$1 short=$2
    local bundle=$source_root/$name out=$working_root/$short
    /usr/bin/time -v -o "$working_root/$short.resource.time.txt" \
        "$python_bin" -m ploidypatch.cli evidence prepare-wgdi \
        --gff "$bundle/primary_chromosomes.gff3" \
        --protein "$protein_root/$name/protein.fa.gz" \
        --fai "$bundle/primary_chromosomes.genome.fa.fai" \
        --output-dir "$out" --prefix "$short" --min-genes-per-seqid 100 \
        > "$working_root/$short.stdout.json" \
        2> "$working_root/$short.stderr.log"
}
export -f prepare_one
export source_root protein_root working_root python_bin
printf '%s\t%s\n' zea_mays zma sorghum_bicolor sbi setaria_italica sit \
    | "$parallel_bin" --jobs 3 --delay 1 --colsep '\t' prepare_one {1} {2}

for short in zma sbi sit; do
    for suffix in wgdi.gff wgdi.lens wgdi.pep.fa wgdi_inputs.manifest.json; do
        [[ -s $working_root/$short/$short.$suffix ]] || {
            echo "missing maize WGDI artifact: $short.$suffix" >&2; exit 1;
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
printf 'maize evaluator WGDI inputs frozen: %s\n' "$result_root"
