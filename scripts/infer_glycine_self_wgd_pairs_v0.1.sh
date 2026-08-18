#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "usage: $0 PROJECT_ROOT" >&2
    exit 2
fi

project_root=$(realpath "$1")
code_root=$project_root/code
python_bin=$project_root/envs/ploidypatch-dev/bin/python
synteny_env=$project_root/envs/ploidypatch-synteny
diamond_bin=$synteny_env/bin/diamond
wgdi_bin=$synteny_env/bin/wgdi
bundle_root=$project_root/data/derived/syngap_glycine_v0.1/gma_complete
source_gff=$bundle_root/primary_chromosomes.gff3
source_fai=$bundle_root/primary_chromosomes.genome.fa.fai
source_protein=$project_root/data/public/ensembl_plants_62/glycine_max/annotation/Glycine_max.Glycine_max_v2.1.pep.all.fa.gz
result_root=$project_root/results/evidence/wgdi/glycine_self_v0.1
working_root=${result_root}.working

for required in "$python_bin" "$diamond_bin" "$wgdi_bin" "$source_gff" \
                "$source_fai" "$source_protein"; do
    if [[ ! -s $required ]]; then
        echo "missing or empty Glycine self-WGD input: $required" >&2
        exit 1
    fi
done
if [[ -e $result_root || -e $working_root ]]; then
    echo "refusing to overwrite Glycine self-WGD result: $result_root" >&2
    exit 1
fi

mkdir -p "$working_root/input" "$working_root/db" "$working_root/blast" \
    "$working_root/config" "$working_root/collinearity" \
    "$working_root/pairs" "$working_root/logs"
{
    printf 'field\tvalue\n'
    printf 'code_commit\t%s\n' "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'target\tGlycine_max_v2.1\n'
    printf 'relationship_scope\tcandidate_named_wgd_syntenic_duplicate\n'
    printf 'wgd_event\tglycine_max_paleopolyploid_WGD\n'
    printf 'min_block_pairs\t20\n'
    printf 'require_different_seqids\ttrue\n'
    printf 'require_reciprocal_unique\ttrue\n'
    printf 'benchmark_access\tevaluator_only\n'
    printf 'self_wgd_module_sha256\t%s\n' \
        "$(sha256sum "$code_root/src/ploidypatch/self_wgd_pairs.py" | awk '{print $1}')"
} > "$working_root/run_contract.tsv"

cd "$code_root"
/usr/bin/time -v -o "$working_root/logs/prepare.time.txt" \
    "$python_bin" -m ploidypatch.cli evidence prepare-wgdi \
        --gff "$source_gff" \
        --protein "$source_protein" \
        --fai "$source_fai" \
        --output-dir "$working_root/input" \
        --prefix gma_v21 \
        --min-genes-per-seqid 100 \
        --primary-chromosomes-only \
        > "$working_root/logs/prepare.stdout.json" \
        2> "$working_root/logs/prepare.stderr.log"

/usr/bin/time -v -o "$working_root/logs/diamond_makedb.time.txt" \
    "$diamond_bin" makedb \
        --in "$working_root/input/gma_v21.wgdi.pep.fa" \
        --db "$working_root/db/gma_v21" \
        > "$working_root/logs/diamond_makedb.stdout.log" \
        2> "$working_root/logs/diamond_makedb.stderr.log"
/usr/bin/time -v -o "$working_root/logs/diamond_self.time.txt" \
    "$diamond_bin" blastp \
        --query "$working_root/input/gma_v21.wgdi.pep.fa" \
        --db "$working_root/db/gma_v21" \
        --out "$working_root/blast/gma_v21_self.tsv" \
        --outfmt 6 --evalue 1e-5 --max-target-seqs 20 \
        --more-sensitive --threads 96 \
        > "$working_root/logs/diamond_self.stdout.log" \
        2> "$working_root/logs/diamond_self.stderr.log"

config=$working_root/config/gma_v21_self.conf
{
    printf '[collinearity]\n'
    printf 'gff1 = %s\n' "$working_root/input/gma_v21.wgdi.gff"
    printf 'gff2 = %s\n' "$working_root/input/gma_v21.wgdi.gff"
    printf 'lens1 = %s\n' "$working_root/input/gma_v21.wgdi.lens"
    printf 'lens2 = %s\n' "$working_root/input/gma_v21.wgdi.lens"
    printf 'blast = %s\n' "$working_root/blast/gma_v21_self.tsv"
    printf 'blast_reverse = false\ncomparison = genomes\nmultiple = 2\n'
    printf 'process = 64\nevalue = 1e-5\nscore = 100\n'
    printf 'grading = 50,40,25\nmg = 40,40\npvalue = 0.2\n'
    printf 'repeat_number = 20\nposition = order\n'
    printf 'savefile = %s\n' "$working_root/collinearity/gma_v21_self.tsv"
} > "$config"
/usr/bin/time -v -o "$working_root/logs/wgdi_self.time.txt" \
    "$wgdi_bin" -icl "$config" \
    > "$working_root/logs/wgdi_self.stdout.log" \
    2> "$working_root/logs/wgdi_self.stderr.log"

/usr/bin/time -v -o "$working_root/logs/pair_inference.time.txt" \
    "$python_bin" -m ploidypatch.cli evidence infer-self-wgd-pairs \
        --query-wgdi-gff "$working_root/input/gma_v21.wgdi.gff" \
        --collinearity "$working_root/collinearity/gma_v21_self.tsv" \
        --wgd-event glycine_max_paleopolyploid_WGD \
        --min-block-pairs 20 \
        --output-pairs "$working_root/pairs/gma_v21.self_wgd_pairs.tsv" \
        --decisions-tsv "$working_root/pairs/decisions.tsv" \
        > "$working_root/pairs/stdout.json" \
        2> "$working_root/pairs/stderr.log"

for output in "$working_root/input/gma_v21.wgdi.gff" \
              "$working_root/blast/gma_v21_self.tsv" \
              "$working_root/collinearity/gma_v21_self.tsv" \
              "$working_root/pairs/gma_v21.self_wgd_pairs.tsv" \
              "$working_root/pairs/gma_v21.self_wgd_pairs.tsv.manifest.json" \
              "$working_root/pairs/decisions.tsv" "$working_root/pairs/stdout.json"; do
    if [[ ! -s $output ]]; then
        echo "missing or empty Glycine self-WGD output: $output" >&2
        exit 1
    fi
done
(
    cd "$working_root"
    find . -type f \( -name '*.json' -o -name '*.tsv' -o -name '*.gff' \) \
        -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS
)
du -sb "$working_root" > "$working_root/disk_bytes.txt"
mv "$working_root" "$result_root"
printf 'Glycine self-WGD pairs frozen: %s\n' "$result_root"
