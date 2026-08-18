#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 DATASET_DIR" >&2
    exit 2
fi

dataset_dir=$1
if [[ ! -d "$dataset_dir" ]]; then
    echo "Dataset directory not found: $dataset_dir" >&2
    exit 2
fi

cd "$dataset_dir"

printf 'section\tfile\tmetric1\tmetric2\tmetric3\tmetric4\n'

for fai in ./*.fa.fai; do
    [[ -e "$fai" ]] || continue
    awk -v file="${fai#./}" '
        BEGIN { n = 0; total = 0; max_len = 0 }
        {
            n++
            total += $2
            if ($2 > max_len) max_len = $2
        }
        END {
            printf "assembly\t%s\t%d\t%d\t%d\tNA\n", file, n, total, max_len
        }
    ' "$fai"
done

for gff in ./*.gff3; do
    [[ -e "$gff" ]] || continue
    awk -F '\t' -v file="${gff#./}" '
        BEGIN { gene = 0; tx = 0; exon = 0; cds = 0 }
        !/^#/ && NF >= 9 {
            if ($3 == "gene") gene++
            else if ($3 == "mRNA" || $3 == "transcript") tx++
            else if ($3 == "exon") exon++
            else if ($3 == "CDS") cds++
        }
        END {
            printf "gff3\t%s\t%d\t%d\t%d\t%d\n", file, gene, tx, exon, cds
        }
    ' "$gff"
done

for protein in ./*.protein.fa; do
    [[ -e "$protein" ]] || continue
    awk -v file="${protein#./}" '
        /^>/ { n++ }
        END { printf "protein\t%s\t%d\tNA\tNA\tNA\n", file, n }
    ' "$protein"
done

for cds_file in ./*.cds.fa; do
    [[ -e "$cds_file" ]] || continue
    awk -v file="${cds_file#./}" '
        /^>/ { n++ }
        END { printf "cds\t%s\t%d\tNA\tNA\tNA\n", file, n }
    ' "$cds_file"
done
