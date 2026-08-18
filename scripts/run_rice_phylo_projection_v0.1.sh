#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "usage: $0 PROJECT_ROOT" >&2
    exit 2
fi

project_root=$(realpath "$1")
code_root=$project_root/code
python_bin=$project_root/envs/ploidypatch-dev/bin/python
miniprot_bin=$project_root/envs/ploidypatch-baseline/bin/miniprot
public_root=$project_root/data/public/ensembl_plants_62
source_gff=$project_root/data/derived/structure_sources_v0.1/osa_irgsp10/source.gff3
target_genome=$public_root/oryza_sativa/genome/Oryza_sativa.IRGSP-1.0.dna.toplevel.fa.gz
reference_root=$public_root/oryza_phylo_references_v0.1
oru_protein=$reference_root/Oryza_rufipogon.OR_W1943.pep.all.fa.gz
ogl_protein=$reference_root/Oryza_glaberrima.Oryza_glaberrima_V1.pep.all.fa.gz
obr_protein=$reference_root/Oryza_brachyantha.Oryza_brachyantha.v1.4b.pep.all.fa.gz
result_root=$project_root/results/projections/miniprot_v0.18/osa_irgsp10_phylo_v0.1
working_root=${result_root}.working

for required in "$python_bin" "$miniprot_bin" "$source_gff" \
                "$target_genome" "$oru_protein" "$ogl_protein" \
                "$obr_protein"; do
    if [[ ! -s $required ]]; then
        echo "missing or empty rice phylogenetic projection input: $required" >&2
        exit 1
    fi
done
if [[ -e $result_root || -e $working_root ]]; then
    echo "refusing to overwrite rice phylogenetic projection: $result_root" >&2
    exit 1
fi

mkdir -p "$working_root/index" "$working_root/projection"
{
    printf 'field\tvalue\n'
    printf 'target\tOryza_sativa_IRGSP-1.0\n'
    printf 'validation_role\theldout_phylogenetically_independent\n'
    printf 'code_commit\t%s\n' "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'miniprot_version\t%s\n' "$($miniprot_bin --version 2>&1 | head -1)"
    printf 'references\tOryza_rufipogon_OR_W1943,Oryza_glaberrima_V1,Oryza_brachyantha_v1.4b\n'
    printf 'reference_genome_groups\tAsian_AA,African_AA,FF\n'
    printf 'target_derived_reference_proteins\tfalse\n'
    printf 'miniprot_parameters\t-I -t64 --gff-only -N4 --outn 4 --outs 0.8 --outc 0.5\n'
    printf 'parameters_frozen_from\tath_tair10_v0.1\n'
} > "$working_root/run_contract.tsv"
{
    printf 'role\tbytes\tsha256\tpath\n'
    for entry in \
        "source_gff:$source_gff" \
        "target_genome:$target_genome" \
        "oru_protein:$oru_protein" \
        "ogl_protein:$ogl_protein" \
        "obr_protein:$obr_protein"; do
        role=${entry%%:*}
        path=${entry#*:}
        printf '%s\t%s\t%s\t%s\n' "$role" "$(stat -Lc %s "$path")" \
            "$(sha256sum "$path" | awk '{print $1}')" "$path"
    done
} > "$working_root/input_manifest.tsv"

cd "$code_root"
"$python_bin" -m ploidypatch.cli baseline prepare-proteins \
    --protein "oru=$oru_protein" \
    --protein "ogl=$ogl_protein" \
    --protein "obr=$obr_protein" \
    --output-fasta "$working_root/projection/references.protein.fa" \
    --output-map "$working_root/projection/references.protein_map.tsv" \
    > "$working_root/projection/prepare.stdout.json" \
    2> "$working_root/projection/prepare.stderr.log"
/usr/bin/time -v -o "$working_root/index/resource.time.txt" \
    "$miniprot_bin" -t64 -d "$working_root/index/osa_irgsp10.mpi" \
        "$target_genome" \
        > "$working_root/index/stdout.log" \
        2> "$working_root/index/stderr.log"
/usr/bin/time -v -o "$working_root/projection/resource.time.txt" \
    "$miniprot_bin" -I -t64 --gff-only -N4 --outn 4 --outs 0.8 --outc 0.5 \
        "$working_root/index/osa_irgsp10.mpi" \
        "$working_root/projection/references.protein.fa" \
        > "$working_root/projection/miniprot.gff3" \
        2> "$working_root/projection/miniprot.stderr.log"

for output in \
    "$working_root/projection/references.protein.fa" \
    "$working_root/projection/references.protein_map.tsv" \
    "$working_root/index/osa_irgsp10.mpi" \
    "$working_root/projection/miniprot.gff3"; do
    if [[ ! -s $output ]]; then
        echo "missing or empty rice phylogenetic projection output: $output" >&2
        exit 1
    fi
done
if ! grep -q '^##gff-version' "$working_root/projection/miniprot.gff3"; then
    echo "rice phylogenetic miniprot output lacks GFF sentinel" >&2
    exit 1
fi
(
    cd "$working_root"
    find . -type f \( -name '*.gff3' -o -name '*.json' -o -name '*.tsv' \) \
        -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS
)
du -sb "$working_root" > "$working_root/disk_bytes.txt"
mv "$working_root" "$result_root"
printf 'rice phylogenetic miniprot projection completed: %s\n' "$result_root"
