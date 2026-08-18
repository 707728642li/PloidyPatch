#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "usage: $0 PROJECT_ROOT" >&2
    exit 2
fi

project_root=$(realpath "$1")
public_source=/nas_data/NFS/Public_genome_data/Ensemble/fasta
result_root=$project_root/data/public/ensembl_plants_62/oryza_phylo_references_v0.1
working_root=${result_root}.working
oru=$public_source/oryza_rufipogon/pep/Oryza_rufipogon.OR_W1943.pep.all.fa.gz
ogl=$public_source/oryza_glaberrima/pep/Oryza_glaberrima.Oryza_glaberrima_V1.pep.all.fa.gz
obr=$public_source/oryza_brachyantha/pep/Oryza_brachyantha.Oryza_brachyantha.v1.4b.pep.all.fa.gz

for required in "$oru" "$ogl" "$obr"; do
    if [[ ! -s $required ]]; then
        echo "missing or empty public Oryza protein input: $required" >&2
        exit 1
    fi
done
if [[ -e $result_root || -e $working_root ]]; then
    echo "refusing to overwrite staged Oryza references: $result_root" >&2
    exit 1
fi

mkdir -p "$working_root"
{
    printf 'source\tspecies\tassembly\tsource_bytes\tsource_sha256\tstaged_file\n'
    for entry in \
        "oru:Oryza_rufipogon:OR_W1943:$oru" \
        "ogl:Oryza_glaberrima:Oryza_glaberrima_V1:$ogl" \
        "obr:Oryza_brachyantha:Oryza_brachyantha.v1.4b:$obr"; do
        source=${entry%%:*}
        remainder=${entry#*:}
        species=${remainder%%:*}
        remainder=${remainder#*:}
        assembly=${remainder%%:*}
        path=${remainder#*:}
        staged=$working_root/$(basename "$path")
        cp "$path" "$staged"
        if ! cmp -s "$path" "$staged"; then
            echo "staged Oryza protein is not byte-identical: $source" >&2
            exit 1
        fi
        printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
            "$source" "$species" "$assembly" "$(stat -Lc %s "$path")" \
            "$(sha256sum "$path" | awk '{print $1}')" "$(basename "$staged")"
    done
} > "$working_root/input_manifest.tsv"
(
    cd "$working_root"
    sha256sum ./*.fa.gz ./input_manifest.tsv > SHA256SUMS
)
du -sb "$working_root" > "$working_root/disk_bytes.txt"
mv "$working_root" "$result_root"
printf 'Oryza phylogenetic reference proteins staged: %s\n' "$result_root"
