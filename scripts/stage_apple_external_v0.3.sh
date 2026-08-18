#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then echo "usage: $0 PROJECT_ROOT" >&2; exit 2; fi
project_root=$(realpath "$1")
result_root=$project_root/data/public/apple_external_v0.3
working_root=${result_root}.working
nfs=/nas_data/NFS/Public_genome_data

declare -A source=(
    [target_apple/genome.fa.gz]="$nfs/GDR_database/raw_data/Malus_x_domestica_Genome_GDDH13_v1.1/GDDH13_1-1_formatted.fasta.gz"
    [target_apple/annotation.gff3.gz]="$nfs/GDR_database/raw_data/Malus_x_domestica_Genome_GDDH13_v1.1/gene_models_20170612.gff3.gz"
    [candidate_pear/genome.fa.gz]="$nfs/GDR_database/raw_data/Pyrus_communis_Bartlett_DH_Genome_v2.0/PyrusCommunis_BartlettDHv2.0.fasta.gz"
    [candidate_pear/annotation.gff3.gz]="$nfs/GDR_database/raw_data/Pyrus_communis_Bartlett_DH_Genome_v2.0/PyrusCommunis_BartlettDHv2.0.gff.gz"
    [candidate_pear/protein.fa.gz]="$nfs/plantgarden_database/raw_data/t23211/PyrusCommunis_BartlettDHv2.0.pep.fasta.gz"
    [candidate_peach/genome.fa.gz]="$nfs/genome_database/prunus_persica/dna/Prunus_persica.Prunus_persica_NCBIv2.dna.toplevel.fa.gz"
    [candidate_peach/annotation.gff3.gz]="$nfs/Ensemble/gff3/prunus_persica/Prunus_persica.Prunus_persica_NCBIv2.62.gff3.gz"
    [candidate_peach/protein.fa.gz]="$nfs/genome_database/prunus_persica/pep/Prunus_persica.Prunus_persica_NCBIv2.pep.all.fa.gz"
    [evaluator_rose/genome.fa.gz]="$nfs/genome_database/rosa_chinensis/dna/Rosa_chinensis.RchiOBHm-V2.dna.toplevel.fa.gz"
    [evaluator_rose/annotation.gff3.gz]="$nfs/Ensemble/gff3/rosa_chinensis/Rosa_chinensis.RchiOBHm-V2.62.gff3.gz"
    [evaluator_strawberry/genome.fa.gz]="$nfs/GDR_database/raw_data/Fragaria_vesca_Genome_v4.0/Fragaria_vesca_v4.0.a1.fasta.gz"
    [evaluator_strawberry/annotation.gff3.gz]="$nfs/GDR_database/raw_data/Fragaria_vesca_Genome_v4.0/Fragaria_vesca_v4.0.a1.transcripts.gff3.gz"
)

[[ ! -e $result_root && ! -e $working_root ]] || {
    echo "refusing to overwrite apple external public bundle" >&2; exit 1;
}
for destination in "${!source[@]}"; do
    [[ -s ${source[$destination]} ]] || {
        echo "missing apple external source: ${source[$destination]}" >&2; exit 1;
    }
done
mkdir -p "$working_root"
{
    printf 'field\tvalue\n'
    printf 'code_commit\t%s\n' "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'stage\tpublic_source_copy_only\n'
    printf 'hidden_pair_enumeration\tfalse\n'
    printf 'hidden_event_generation\tfalse\n'
    printf 'external_label_access\tfalse\n'
    printf 'target\tMalus_x_domestica_GDDH13_v1.1\n'
    printf 'candidate_references\tPyrus_communis_BartlettDH_v2.0,Prunus_persica_NCBIv2_Ensembl62\n'
    printf 'evaluator_only_references\tRosa_chinensis_RchiOBHm-V2_Ensembl62,Fragaria_vesca_v4.0\n'
} > "$working_root/run_contract.tsv"
printf 'destination\tsource_path\tsource_bytes\tsource_sha256\n' > "$working_root/input_manifest.tsv"
for destination in "${!source[@]}"; do
    path=${source[$destination]}
    mkdir -p "$working_root/$(dirname "$destination")"
    cp --reflink=auto --preserve=timestamps "$path" "$working_root/$destination"
    source_hash=$(sha256sum "$path" | awk '{print $1}')
    copied_hash=$(sha256sum "$working_root/$destination" | awk '{print $1}')
    [[ $source_hash == "$copied_hash" ]] || { echo "copy checksum mismatch: $destination" >&2; exit 1; }
    printf '%s\t%s\t%s\t%s\n' "$destination" "$path" \
        "$(stat -Lc %s "$path")" "$source_hash" >> "$working_root/input_manifest.tsv"
done
(
    cd "$working_root"
    find . -type f ! -name SHA256SUMS -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS
    sha256sum -c SHA256SUMS >/dev/null
)
du -sb "$working_root" > "$working_root/disk_bytes.txt"
mv "$working_root" "$result_root"
printf 'apple external public bundle frozen: %s\n' "$result_root"

