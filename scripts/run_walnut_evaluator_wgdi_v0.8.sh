#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then echo "usage: $0 PROJECT_ROOT" >&2; exit 2; fi
project_root=$(realpath "$1")
evaluator_root=$(realpath "${PLOIDYPATCH_EVALUATOR_ONLY_ROOT:?PLOIDYPATCH_EVALUATOR_ONLY_ROOT is required}")
execution=$(realpath "${PLOIDYPATCH_EXECUTION_FREEZE:?PLOIDYPATCH_EXECUTION_FREEZE is required}")
code_root=$execution/source
python=$project_root/envs/ploidypatch-dev/bin/python
synteny_env=$project_root/envs/ploidypatch-synteny
diamond=$synteny_env/bin/diamond
wgdi=$synteny_env/bin/wgdi
parallel=${PLOIDYPATCH_PARALLEL_BIN:-/data/codexli/software/conda/miniforge3/bin/parallel}
input=$evaluator_root/wgdi/input
output=$evaluator_root/wgdi/evidence
working=$output.working
ks_helper=$code_root/scripts/prepare_walnut_target_ks_inputs_v0.8.py

for path in "$python" "$diamond" "$wgdi" "$parallel" "$ks_helper" \
    "$input/SHA256SUMS" "$input/target/jre.wgdi.gff" "$input/target/jre.wgdi.pep.fa" \
    "$input/target/jre.wgdi.cds.fa" "$input/target/jre.wgdi.lens" \
    "$input/corylus/cav.wgdi.gff" "$input/corylus/cav.wgdi.pep.fa" "$input/corylus/cav.wgdi.lens" \
    "$input/castanea/cmo.wgdi.gff" "$input/castanea/cmo.wgdi.pep.fa" "$input/castanea/cmo.wgdi.lens"; do
    [[ -s $path ]] || { echo "missing Walnut evaluator WGDI prerequisite: $path" >&2; exit 1; }
done
(cd "$input" && sha256sum -c SHA256SUMS >/dev/null)
[[ ! -e $output && ! -e $working ]] || { echo "refusing to overwrite Walnut WGDI evidence" >&2; exit 1; }
mkdir -p "$working"/{db,blast,collinearity,config,logs,ks/jobs}
trap 'status=$?; printf "status\tinvalid\nexit_status\t%s\n" "$status" > "$working/stage_status.tsv" || true; exit "$status"' ERR

parallel_version=$("$parallel" --version)
parallel_version=${parallel_version%%$'\n'*}
{
    printf 'tool\tpath\tbytes\tsha256\tversion\n'
    printf 'gnu_parallel\t%s\t%s\t%s\t%s\n' \
        "$parallel" "$(stat -Lc %s "$parallel")" \
        "$(sha256sum "$parallel" | awk '{print $1}')" "$parallel_version"
} > "$working/tool_manifest.tsv"

for entry in target:jre corylus:cav castanea:cmo; do
    label=${entry%%:*}; prefix=${entry#*:}
    "$diamond" makedb --in "$input/$label/$prefix.wgdi.pep.fa" --db "$working/db/$prefix" \
        >"$working/logs/$prefix.makedb.stdout.log" 2>"$working/logs/$prefix.makedb.stderr.log"
done
run_blast() {
    local label=$1 query=$2 db=$3
    "$diamond" blastp --query "$query" --db "$db" --out "$working/blast/$label.tsv" \
        --outfmt 6 --evalue 1e-5 --max-target-seqs 20 --more-sensitive --threads 40 \
        >"$working/logs/$label.diamond.stdout.log" 2>"$working/logs/$label.diamond.stderr.log"
}
run_blast jre_self "$input/target/jre.wgdi.pep.fa" "$working/db/jre" & p1=$!
run_blast jre_vs_cav "$input/target/jre.wgdi.pep.fa" "$working/db/cav" & p2=$!
run_blast jre_vs_cmo "$input/target/jre.wgdi.pep.fa" "$working/db/cmo" & p3=$!
wait "$p1"; wait "$p2"; wait "$p3"

write_conf() {
    local name=$1 q=$2 s=$3 qprefix=$4 sprefix=$5 multiple=$6
    cat > "$working/config/$name.conf" <<EOF
[collinearity]
gff1 = $input/$q/$qprefix.wgdi.gff
gff2 = $input/$s/$sprefix.wgdi.gff
lens1 = $input/$q/$qprefix.wgdi.lens
lens2 = $input/$s/$sprefix.wgdi.lens
blast = $working/blast/$name.tsv
blast_reverse = false
comparison = genomes
multiple = $multiple
process = 64
evalue = 1e-5
score = 100
grading = 50,40,25
mg = 40,40
pvalue = 0.2
repeat_number = 20
position = order
savefile = $working/collinearity/$name.tsv
EOF
}
write_conf jre_self target target jre jre 2
write_conf jre_vs_cav target corylus jre cav 2
write_conf jre_vs_cmo target castanea jre cmo 2
for name in jre_self jre_vs_cav jre_vs_cmo; do
    "$wgdi" -icl "$working/config/$name.conf" >"$working/logs/$name.wgdi.stdout.log" 2>"$working/logs/$name.wgdi.stderr.log" &
    eval "pid_$name=$!"
done
wait "$pid_jre_self"; wait "$pid_jre_vs_cav"; wait "$pid_jre_vs_cmo"

PYTHONPATH="$code_root/src" "$python" "$ks_helper" \
    --collinearity "$working/collinearity/jre_self.tsv" \
    --wgdi-gff "$input/target/jre.wgdi.gff" --output-dir "$working/ks/target" --shards 64
export working wgdi input
find "$working/ks/target/shards" -name 'pairs.*.tsv' -type f -print0 | sort -z | \
    "$parallel" -0 --jobs 64 --delay 0.02 --halt now,fail=1 '
      shard={}; stem={/.}; job="$working/ks/jobs/$stem"; mkdir "$job";
      printf "[ks]\ncds_file = %s\npep_file = %s\nalign_software = muscle\npairs_file = %s\nks_file = %s\n" \
        "$input/target/jre.wgdi.cds.fa" "$input/target/jre.wgdi.pep.fa" "$shard" "$job/ks.tsv" > "$job/ks.conf";
      cd "$job"; "$wgdi" -ks ks.conf >stdout.log 2>stderr.log
    '
PYTHONPATH="$code_root/src" "$python" - "$working/ks" <<'PY'
import csv, sys
from pathlib import Path
root=Path(sys.argv[1]); rows=[]; seen=set(); fields=['id1','id2','ka_NG86','ks_NG86','ka_YN00','ks_YN00']
for path in sorted(root.glob('jobs/pairs.*/ks.tsv')):
    with path.open() as h:
        reader=csv.DictReader(h, delimiter='\t')
        if reader.fieldnames != fields: raise SystemExit(f'invalid YN00 header: {path}')
        for row in reader:
            key=tuple(sorted((row['id1'],row['id2'])))
            if key in seen: raise SystemExit(f'duplicate YN00 pair: {key}')
            seen.add(key); rows.append(row)
with (root/'target/ks_merged.tsv').open('x', newline='') as h:
    writer=csv.DictWriter(h, fieldnames=fields, delimiter='\t', lineterminator='\n'); writer.writeheader(); writer.writerows(sorted(rows,key=lambda x:(x['id1'],x['id2'])))
PY
printf 'status\tready_for_pair_inference\n' > "$working/stage_status.tsv"
(
    cd "$working"
    find . -type f ! -name SHA256SUMS -printf '%P\0' | sort -z | xargs -0 sha256sum > SHA256SUMS
    sha256sum -c SHA256SUMS >/dev/null
)
trap - ERR
mv "$working" "$output"
