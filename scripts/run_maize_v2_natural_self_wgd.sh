#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then echo "usage: $0 PROJECT_ROOT" >&2; exit 2; fi
project_root=$(realpath "$1")
code_root=$project_root/code
python_bin=$project_root/envs/ploidypatch-dev/bin/python
synteny_env=$project_root/envs/ploidypatch-synteny
syngap_env=$project_root/envs/ploidypatch-syngap
diamond_bin=$synteny_env/bin/diamond
wgdi_bin=$synteny_env/bin/wgdi
bundle=$project_root/data/derived/holdout_inputs/maize_v2/zea_mays
genome=$bundle/primary_chromosomes.genome.fa
fai=$bundle/primary_chromosomes.genome.fa.fai
base_gff=$bundle/primary_chromosomes.gff3
method_root=$project_root/results/natural/maize_v2/discovery/method_trio
candidate_gff=$method_root/consensus/union/natural/candidate.gff3
result_root=$project_root/results/natural/maize_v2/discovery/union_self_wgd
working_root=${result_root}.working
prefix=zma_union_candidate_natural

[[ -s $method_root/SHA256SUMS ]] || { echo "unfrozen maize natural candidates" >&2; exit 1; }
(cd "$method_root" && sha256sum -c SHA256SUMS >/dev/null)
for required in "$python_bin" "$diamond_bin" "$wgdi_bin" \
                "$syngap_env/bin/gffread" "$genome" "$fai" \
                "$base_gff" "$candidate_gff"; do
    [[ -s $required ]] || { echo "missing maize natural self-WGD input: $required" >&2; exit 1; }
done
if [[ -e $result_root || -e $working_root ]]; then
    echo "refusing to overwrite maize natural self-WGD" >&2; exit 1
fi
mkdir -p "$working_root"/{input,db,blast,config,collinearity,pairs,selected,logs}
{
    printf 'field\tvalue\ncode_commit\t%s\n' \
        "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'mode\tnatural_current_annotation_v0.1\n'
    printf 'candidate_truth_access\tfalse\nvalidation_evidence_access\tfalse\n'
    printf 'candidate_source\tmethod_consensus_support1_union\n'
    printf 'wgd_event\tzea_mays_lineage_tetraploidy\nmin_block_pairs\t20\n'
    printf 'require_different_seqids\ttrue\nrequire_reciprocal_unique\ttrue\n'
    printf 'candidate_candidate_pair_policy\treject_as_circular\n'
    printf 'diamond_threads\t64\nwgdi_processes\t64\nautomatic_approval\tfalse\n'
} > "$working_root/run_contract.tsv"
{
    printf 'role\tbytes\tsha256\tpath\n'
    for entry in "base_gff:$base_gff" "candidate_gff:$candidate_gff" \
                 "genome:$genome" "fai:$fai"; do
        role=${entry%%:*}; path=${entry#*:}
        printf '%s\t%s\t%s\t%s\n' "$role" "$(stat -Lc %s "$path")" \
            "$(sha256sum "$path" | awk '{print $1}')" "$path"
    done
} > "$working_root/candidate_input_manifest.tsv"

/usr/bin/time -v -o "$working_root/logs/gffread_pep.time.txt" \
    conda run -p "$syngap_env" --no-capture-output \
    gffread "$candidate_gff" -g "$genome" \
        -y "$working_root/input/$prefix.all.pep.fa" -S \
        > "$working_root/logs/gffread_pep.stdout.log" \
        2> "$working_root/logs/gffread_pep.stderr.log"
cd "$code_root"
/usr/bin/time -v -o "$working_root/logs/prepare.time.txt" \
    "$python_bin" -m ploidypatch.cli evidence prepare-wgdi \
        --gff "$candidate_gff" --protein "$working_root/input/$prefix.all.pep.fa" \
        --fai "$fai" --output-dir "$working_root/input" \
        --prefix "$prefix" --min-genes-per-seqid 100 \
        > "$working_root/logs/prepare.stdout.json" \
        2> "$working_root/logs/prepare.stderr.log"
/usr/bin/time -v -o "$working_root/logs/diamond_makedb.time.txt" \
    "$diamond_bin" makedb --in "$working_root/input/$prefix.wgdi.pep.fa" \
        --db "$working_root/db/$prefix" \
        > "$working_root/logs/diamond_makedb.stdout.log" \
        2> "$working_root/logs/diamond_makedb.stderr.log"
/usr/bin/time -v -o "$working_root/logs/diamond_self.time.txt" \
    "$diamond_bin" blastp --query "$working_root/input/$prefix.wgdi.pep.fa" \
        --db "$working_root/db/$prefix" \
        --out "$working_root/blast/${prefix}_self.tsv" --outfmt 6 --evalue 1e-5 \
        --max-target-seqs 20 --more-sensitive --threads 64 \
        > "$working_root/logs/diamond_self.stdout.log" \
        2> "$working_root/logs/diamond_self.stderr.log"
config=$working_root/config/${prefix}_self.conf
{
    printf '[collinearity]\n'
    printf 'gff1 = %s\ngff2 = %s\n' "$working_root/input/$prefix.wgdi.gff" \
        "$working_root/input/$prefix.wgdi.gff"
    printf 'lens1 = %s\nlens2 = %s\n' "$working_root/input/$prefix.wgdi.lens" \
        "$working_root/input/$prefix.wgdi.lens"
    printf 'blast = %s\n' "$working_root/blast/${prefix}_self.tsv"
    printf 'blast_reverse = false\ncomparison = genomes\nmultiple = 2\nprocess = 64\n'
    printf 'evalue = 1e-5\nscore = 100\ngrading = 50,40,25\nmg = 40,40\n'
    printf 'pvalue = 0.2\nrepeat_number = 20\nposition = order\n'
    printf 'savefile = %s\n' "$working_root/collinearity/${prefix}_self.tsv"
} > "$config"
/usr/bin/time -v -o "$working_root/logs/wgdi_self.time.txt" \
    "$wgdi_bin" -icl "$config" > "$working_root/logs/wgdi_self.stdout.log" \
    2> "$working_root/logs/wgdi_self.stderr.log"
/usr/bin/time -v -o "$working_root/logs/pair_inference.time.txt" \
    "$python_bin" -m ploidypatch.cli evidence infer-self-wgd-pairs \
        --query-wgdi-gff "$working_root/input/$prefix.wgdi.gff" \
        --collinearity "$working_root/collinearity/${prefix}_self.tsv" \
        --source-gff "$candidate_gff" \
        --wgd-event zea_mays_lineage_tetraploidy --min-block-pairs 20 \
        --output-pairs "$working_root/pairs/${prefix}.self_wgd_pairs.tsv" \
        --decisions-tsv "$working_root/pairs/decisions.tsv" \
        > "$working_root/pairs/stdout.json" 2> "$working_root/pairs/stderr.log"
/usr/bin/time -v -o "$working_root/logs/candidate_selection.time.txt" \
    "$python_bin" -m ploidypatch.cli evidence select-wgd-supported-candidates \
        --base-gff "$base_gff" --candidate-gff "$candidate_gff" \
        --pairs "$working_root/pairs/${prefix}.self_wgd_pairs.tsv" \
        --output-gff "$working_root/selected/candidate.gff3" \
        --selection-tsv "$working_root/selected/selection.tsv" \
        > "$working_root/selected/stdout.json" \
        2> "$working_root/selected/stderr.log"
for output in "$working_root/input/$prefix.wgdi.gff" \
              "$working_root/input/$prefix.wgdi.pep.fa" \
              "$working_root/blast/${prefix}_self.tsv" \
              "$working_root/collinearity/${prefix}_self.tsv" \
              "$working_root/pairs/${prefix}.self_wgd_pairs.tsv" \
              "$working_root/selected/candidate.gff3" \
              "$working_root/selected/selection.tsv"; do
    [[ -s $output ]] || { echo "missing maize natural self-WGD output: $output" >&2; exit 1; }
done
(
    cd "$working_root"
    find . -type f \( -name '*.json' -o -name '*.tsv' -o -name '*.gff' \
        -o -name '*.gff3' -o -name '*.fa' \) -print0 \
        | sort -z | xargs -0 sha256sum > SHA256SUMS
    sha256sum -c SHA256SUMS >/dev/null
)
du -sb "$working_root" > "$working_root/disk_bytes.txt"
mv "$working_root" "$result_root"
printf 'maize natural self-WGD frozen: %s\n' "$result_root"

