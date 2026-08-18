#!/usr/bin/env bash
set -euo pipefail

project_root=${PLOIDYPATCH_PROJECT_ROOT:-/data/codexli/projects/PloidyPatch}
synteny_env=$project_root/envs/ploidypatch-synteny
dev_env=$project_root/envs/ploidypatch-dev
freeze_root=$project_root/results/frozen/ploidypatch_v0.6_155f72a
code_root=$freeze_root/source
input_root=$project_root/results/quality/helianthus_v0.6/base_only_wgdi_inputs
primary_root=$project_root/data/derived/quality_preflight/helianthus_v0.6_primary
primary_table=$project_root/config/primary_seqids/helianthus_annuus_hanxrq_r2.0.tsv
result_root=$project_root/results/quality/helianthus_v0.6/fixed_target_backbone
working_root=$result_root.working

python_bin=$dev_env/bin/python
diamond_bin=$synteny_env/bin/diamond
wgdi_bin=$synteny_env/bin/wgdi

for required in \
    "$python_bin" "$diamond_bin" "$wgdi_bin" \
    "$freeze_root/source.tar.gz" "$freeze_root/SHA256SUMS" \
    "$code_root/src/ploidypatch/fixed_target_backbone.py" \
    "$code_root/src/ploidypatch/cli.py" \
    "$input_root/han.wgdi.gff" "$input_root/han.wgdi.lens" \
    "$input_root/han.wgdi.pep.fa" "$input_root/han.representatives.tsv" \
    "$input_root/han.wgdi_inputs.manifest.json" \
    "$primary_root/primary_chromosomes.gff3" "$primary_table"; do
    [[ -s $required ]] || { echo "missing fixed-backbone input: $required" >&2; exit 1; }
done
[[ ! -e $result_root && ! -e $working_root ]] || {
    echo "refusing to overwrite fixed-backbone preflight" >&2
    exit 1
}
(cd "$freeze_root" && sha256sum -c SHA256SUMS >/dev/null)
[[ $($diamond_bin version | awk '{print $3}') == 2.2.2 ]]
[[ $($wgdi_bin --version 2>&1 | awk 'NR==1 {print $NF}') == 0.75 ]]

mkdir -p "$working_root"/{input,db,blast,config,collinearity,backbone,logs}
cp "$input_root"/han.wgdi.gff "$working_root/input/"
cp "$input_root"/han.wgdi.lens "$working_root/input/"
cp "$input_root"/han.wgdi.pep.fa "$working_root/input/"
cp "$input_root"/han.representatives.tsv "$working_root/input/"
cp "$input_root"/han.wgdi_inputs.manifest.json "$working_root/input/"
{
    printf 'schema_version\tploidypatch.helianthus_fixed_backbone_run.v0.6\n'
    printf 'role\ttarget_quality_only\n'
    printf 'candidate_access\tfalse\ntruth_or_label_access\tfalse\n'
    printf 'code_commit\t155f72a\n'
    printf 'runner_sha256\t%s\n' "$(sha256sum "$0" | awk '{print $1}')"
    printf 'diamond_version\t2.2.2\nwgdi_version\t0.75\n'
    printf 'diamond_threads\t64\nwgdi_processes\t64\n'
    printf 'diamond_evalue\t1e-5\ndiamond_max_target_seqs\t20\n'
    printf 'diamond_mode\tmore-sensitive\n'
    printf 'wgdi_multiple\t2\nwgdi_score\t100\nwgdi_grading\t50,40,25\n'
    printf 'wgdi_mg\t40,40\nwgdi_pvalue\t0.2\nwgdi_repeat_number\t20\n'
    printf 'backbone_min_block_pairs\t20\n'
} > "$working_root/config/run_contract.tsv"

/usr/bin/time -v -o "$working_root/logs/diamond_makedb.time.txt" \
    "$diamond_bin" makedb --in "$working_root/input/han.wgdi.pep.fa" \
    --db "$working_root/db/han" \
    > "$working_root/logs/diamond_makedb.stdout.log" \
    2> "$working_root/logs/diamond_makedb.stderr.log"
/usr/bin/time -v -o "$working_root/logs/diamond_self.time.txt" \
    "$diamond_bin" blastp --query "$working_root/input/han.wgdi.pep.fa" \
    --db "$working_root/db/han" --out "$working_root/blast/han_self.tsv" \
    --outfmt 6 --evalue 1e-5 --max-target-seqs 20 --more-sensitive --threads 64 \
    > "$working_root/logs/diamond_self.stdout.log" \
    2> "$working_root/logs/diamond_self.stderr.log"

config=$working_root/config/han_self.conf
{
    printf '[collinearity]\n'
    printf 'gff1 = %s\ngff2 = %s\n' "$working_root/input/han.wgdi.gff" \
        "$working_root/input/han.wgdi.gff"
    printf 'lens1 = %s\nlens2 = %s\n' "$working_root/input/han.wgdi.lens" \
        "$working_root/input/han.wgdi.lens"
    printf 'blast = %s\n' "$working_root/blast/han_self.tsv"
    printf 'blast_reverse = false\ncomparison = genomes\nmultiple = 2\nprocess = 64\n'
    printf 'evalue = 1e-5\nscore = 100\ngrading = 50,40,25\nmg = 40,40\n'
    printf 'pvalue = 0.2\nrepeat_number = 20\nposition = order\n'
    printf 'savefile = %s\n' "$working_root/collinearity/han_self.tsv"
} > "$config"
/usr/bin/time -v -o "$working_root/logs/wgdi_self.time.txt" \
    "$wgdi_bin" -icl "$config" \
    > "$working_root/logs/wgdi_self.stdout.log" \
    2> "$working_root/logs/wgdi_self.stderr.log"

PYTHONPATH=$code_root/src /usr/bin/time -v \
    -o "$working_root/logs/backbone.time.txt" \
    "$python_bin" -m ploidypatch.cli evidence build-fixed-target-backbone \
    --base-gff "$primary_root/primary_chromosomes.gff3" \
    --query-wgdi-gff "$working_root/input/han.wgdi.gff" \
    --base-only-collinearity "$working_root/collinearity/han_self.tsv" \
    --primary-seqid-table "$primary_table" --min-block-pairs 20 \
    --output-dir "$working_root/backbone/fixed" \
    > "$working_root/logs/backbone.stdout.json" \
    2> "$working_root/logs/backbone.stderr.log"

for required in "$working_root/blast/han_self.tsv" \
    "$working_root/collinearity/han_self.tsv" \
    "$working_root/backbone/fixed/genes.tsv" \
    "$working_root/backbone/fixed/edges.tsv" \
    "$working_root/backbone/fixed/cells.tsv" \
    "$working_root/backbone/fixed/manifest.json" \
    "$working_root/backbone/fixed/SHA256SUMS"; do
    [[ -s $required ]] || { echo "missing fixed-backbone output: $required" >&2; exit 1; }
done
(
    cd "$working_root"
    find . -type f ! -name SHA256SUMS -printf '%P\0' | sort -z \
        | while IFS= read -r -d '' relative; do sha256sum "$relative"; done \
        > SHA256SUMS
    sha256sum -c SHA256SUMS >/dev/null
)
chmod -R a-w "$working_root"
mv "$working_root" "$result_root"
printf 'Helianthus fixed target backbone complete: %s\n' "$result_root"
