#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then echo "usage: $0 PROJECT_ROOT" >&2; exit 2; fi
project_root=$(realpath "$1")
code_root=$project_root/code
python_bin=$project_root/envs/ploidypatch-dev/bin/python
parallel_bin=/data/codexli/software/conda/miniforge3/bin/parallel
public_root=$project_root/data/public/apple_external_v0.3
result_root=$project_root/data/derived/external_inputs/apple_v0.3
working_root=${result_root}.working

species=(target_apple candidate_pear candidate_peach evaluator_rose evaluator_strawberry)
declare -A primary=(
    [target_apple]="$code_root/config/primary_seqids/malus_gddh13_v1.1.tsv"
    [candidate_pear]="$code_root/config/primary_seqids/pyrus_bartlettdh_v2.0.tsv"
    [candidate_peach]="$code_root/config/primary_seqids/prunus_persica_ncbiv2.tsv"
    [evaluator_rose]="$code_root/config/primary_seqids/rosa_chinensis_rchiobhm_v2.tsv"
    [evaluator_strawberry]="$code_root/config/primary_seqids/fragaria_vesca_v4.0.tsv"
)
declare -A expected=(
    [target_apple]=17 [candidate_pear]=17 [candidate_peach]=8
    [evaluator_rose]=7 [evaluator_strawberry]=7
)
declare -A expected_note_repairs=(
    [target_apple]=3 [candidate_pear]=0 [candidate_peach]=0
    [evaluator_rose]=0 [evaluator_strawberry]=0
)
declare -A expected_dropped_introns=(
    [target_apple]=0 [candidate_pear]=5 [candidate_peach]=0
    [evaluator_rose]=0 [evaluator_strawberry]=0
)
declare -A expected_fasta_directives=(
    [target_apple]=0 [candidate_pear]=0 [candidate_peach]=0
    [evaluator_rose]=0 [evaluator_strawberry]=1
)

for required in "$python_bin" "$parallel_bin" "$public_root/SHA256SUMS"; do
    [[ -s $required ]] || { echo "missing apple normalization prerequisite: $required" >&2; exit 1; }
done
(cd "$public_root" && sha256sum -c SHA256SUMS >/dev/null)
for name in "${species[@]}"; do
    for required in "$public_root/$name/genome.fa.gz" \
                    "$public_root/$name/annotation.gff3.gz" \
                    "${primary[$name]}"; do
        [[ -s $required ]] || { echo "missing apple normalization input: $required" >&2; exit 1; }
    done
done
[[ ! -e $result_root && ! -e $working_root ]] || {
    echo "refusing to overwrite normalized apple external inputs" >&2; exit 1;
}
cd "$code_root"
mkdir -p "$working_root"
{
    printf 'field\tvalue\n'
    printf 'code_commit\t%s\n' "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'stage\tpretruth_format_normalization_only\n'
    printf 'hidden_pair_enumeration\tfalse\n'
    printf 'hidden_event_generation\tfalse\n'
    printf 'external_label_access\tfalse\n'
    printf 'candidate_evaluator_reference_separation\ttrue\n'
    printf 'provider_gff_compatibility\texplicit_narrow_audited\n'
} > "$working_root/run_contract.tsv"
printf 'species\trepair_unescaped_note_semicolons\tdrop_invalid_intron_intervals\tstrip_embedded_fasta\texpected_note_repairs\texpected_dropped_introns\texpected_fasta_directives\n' > "$working_root/provider_compatibility_policy.tsv"
mkdir -p "$working_root/provider_compatible"
for name in "${species[@]}"; do
    compatibility_args=()
    [[ $name == target_apple ]] && compatibility_args+=(--repair-unescaped-note-semicolons)
    [[ $name == candidate_pear ]] && compatibility_args+=(--drop-invalid-intron-intervals)
    [[ $name == evaluator_strawberry ]] && compatibility_args+=(--strip-embedded-fasta)
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$name" \
        "$([[ $name == target_apple ]] && printf true || printf false)" \
        "$([[ $name == candidate_pear ]] && printf true || printf false)" \
        "$([[ $name == evaluator_strawberry ]] && printf true || printf false)" \
        "${expected_note_repairs[$name]}" "${expected_dropped_introns[$name]}" \
        "${expected_fasta_directives[$name]}" \
        >> "$working_root/provider_compatibility_policy.tsv"
    /usr/bin/time -v -o "$working_root/$name.provider_compatibility.resource.time.txt" \
        "$python_bin" -m ploidypatch.cli normalize provider-gff3 \
        --gff "$public_root/$name/annotation.gff3.gz" \
        --output-dir "$working_root/provider_compatible/$name" \
        "${compatibility_args[@]}" \
        > "$working_root/$name.provider_compatibility.stdout.json" \
        2> "$working_root/$name.provider_compatibility.stderr.log"
    read -r repaired dropped fasta_directives < <(
        "$python_bin" - "$working_root/provider_compatible/$name/manifest.json" <<'PY'
import json
import sys

with open(sys.argv[1], encoding="utf-8") as handle:
    observed = json.load(handle)["observed"]
print(
    observed.get("repaired_note_records", 0),
    observed.get("dropped_invalid_intron_intervals", 0),
    observed.get("stripped_embedded_fasta_directives", 0),
)
PY
    )
    [[ $repaired -eq ${expected_note_repairs[$name]} \
        && $dropped -eq ${expected_dropped_introns[$name]} \
        && $fasta_directives -eq ${expected_fasta_directives[$name]} ]] || {
        echo "$name provider compatibility counts differ from frozen expectations" >&2
        exit 1
    }
done
commands=$working_root/normalization_commands.tsv
for name in "${species[@]}"; do
    printf '%s\t%s\t%s\t%s\n' "$name" \
        "$working_root/provider_compatible/$name/sanitized.gff3" \
        "$public_root/$name/genome.fa.gz" "${primary[$name]}" >> "$commands"
done

cd "$code_root"
"$parallel_bin" --colsep '\t' --delay 1 -j 5 \
    --joblog "$working_root/parallel.joblog.tsv" \
    '/usr/bin/time -v -o '"$working_root"'/{1}.resource.time.txt '"$python_bin"' -m ploidypatch.cli normalize primary-annotation --gff {2} --genome {3} --primary-seqid-table {4} --output-dir '"$working_root"'/{1} > '"$working_root"'/{1}.stdout.json 2> '"$working_root"'/{1}.stderr.log' \
    :::: "$commands"

printf 'species\tchromosomes\tgenes\ttranscripts\tcds_segments\n' > "$working_root/annotation_counts.tsv"
for name in "${species[@]}"; do
    bundle=$working_root/$name
    for required in "$bundle/primary_chromosomes.genome.fa" \
                    "$bundle/primary_chromosomes.genome.fa.fai" \
                    "$bundle/primary_chromosomes.gff3" "$bundle/manifest.json"; do
        [[ -s $required ]] || { echo "missing normalized apple artifact: $required" >&2; exit 1; }
    done
    chromosomes=$(wc -l < "$bundle/primary_chromosomes.genome.fa.fai")
    [[ $chromosomes -eq ${expected[$name]} ]] || {
        echo "$name chromosome count $chromosomes is not ${expected[$name]}" >&2; exit 1;
    }
    genes=$(grep -Pc '\tgene\t' "$bundle/primary_chromosomes.gff3")
    transcripts=$(grep -Pc '\t(mRNA|transcript)\t' "$bundle/primary_chromosomes.gff3")
    cds=$(grep -Pc '\tCDS\t' "$bundle/primary_chromosomes.gff3")
    [[ $genes -ge 10000 && $transcripts -ge 10000 && $cds -ge 20000 ]] || {
        echo "$name annotation content floor failed" >&2; exit 1;
    }
    printf '%s\t%s\t%s\t%s\t%s\n' "$name" "$chromosomes" "$genes" "$transcripts" "$cds" >> "$working_root/annotation_counts.tsv"
done
(
    cd "$working_root"
    find . -type f \( -name '*.fa' -o -name '*.fai' -o -name '*.gff3' \
        -o -name '*.json' -o -name '*.tsv' \) -print0 \
        | sort -z | xargs -0 sha256sum > SHA256SUMS
    sha256sum -c SHA256SUMS >/dev/null
)
du -sb "$working_root" > "$working_root/disk_bytes.txt"
mv "$working_root" "$result_root"
printf 'normalized apple external inputs frozen: %s\n' "$result_root"
