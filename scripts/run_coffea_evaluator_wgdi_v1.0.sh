#!/usr/bin/env bash
set -euo pipefail
umask 077

[[ $# -eq 1 ]] || { echo "usage: $0 PROJECT_ROOT" >&2; exit 64; }
project_root=$(realpath "$1")
evaluator_root=$(realpath "${PLOIDYPATCH_EVALUATOR_ONLY_ROOT:?}")
execution=$(realpath "${PLOIDYPATCH_EXECUTION_FREEZE:?}")
code_root=$execution/source
python=$project_root/envs/ploidypatch-dev/bin/python
synteny_env=$project_root/envs/ploidypatch-synteny
diamond=$synteny_env/bin/diamond
wgdi=$synteny_env/bin/wgdi
parallel=${PLOIDYPATCH_PARALLEL_BIN:-/data/codexli/software/conda/miniforge3/bin/parallel}
input=$evaluator_root/wgdi/input
output=$evaluator_root/wgdi/evidence
working=$output.working
ks_helper=$code_root/scripts/prepare_coffea_target_ks_inputs_v1.0.py

for executable in "$python" "$diamond" "$wgdi" "$parallel"; do
    [[ -x $executable ]] || { echo "missing Coffea WGDI executable: $executable" >&2; exit 65; }
    resolved=$(realpath "$executable")
    [[ $resolved != /nas_data && $resolved != /nas_data/* ]] || {
        echo "Coffea WGDI executable resolves into forbidden NAS: $executable" >&2; exit 65;
    }
done
for path in "$ks_helper" "$input/SHA256SUMS" "$input/target/car.wgdi.gff" \
    "$input/target/car.wgdi.pep.fa" "$input/target/car.wgdi.cds.fa" \
    "$input/target/car.wgdi.lens" \
    "$input/gardenia/gja.wgdi.gff" "$input/gardenia/gja.wgdi.pep.fa" \
    "$input/gardenia/gja.wgdi.lens" \
    "$input/ophiorrhiza/opu.wgdi.gff" "$input/ophiorrhiza/opu.wgdi.pep.fa" \
    "$input/ophiorrhiza/opu.wgdi.lens"; do
    [[ -s $path && ! -L $path ]] || { echo "missing Coffea WGDI prerequisite: $path" >&2; exit 65; }
done
(cd "$input" && sha256sum -c SHA256SUMS >/dev/null)
[[ ! -e $output && ! -L $output && ! -e $working && ! -L $working ]] || {
    echo "refusing to overwrite Coffea evaluator WGDI evidence" >&2; exit 66;
}
mkdir -p "$working"/{db,blast,collinearity,config,logs,ks/jobs}
trap 'status=$?; printf "status\tinvalid_run\nexit_status\t%s\n" "$status" > "$working/stage_status.tsv" || true; exit "$status"' ERR

for entry in target:car gardenia:gja ophiorrhiza:opu; do
    label=${entry%%:*}; prefix=${entry#*:}
    "$diamond" makedb --in "$input/$label/$prefix.wgdi.pep.fa" \
        --db "$working/db/$prefix" \
        >"$working/logs/$prefix.makedb.stdout.log" \
        2>"$working/logs/$prefix.makedb.stderr.log"
done
run_blast() {
    local label=$1 query=$2 db=$3
    "$diamond" blastp --query "$query" --db "$db" --out "$working/blast/$label.tsv" \
        --outfmt 6 --evalue 1e-5 --max-target-seqs 20 --more-sensitive --threads 40 \
        >"$working/logs/$label.diamond.stdout.log" \
        2>"$working/logs/$label.diamond.stderr.log"
}
run_blast car_self "$input/target/car.wgdi.pep.fa" "$working/db/car" & p1=$!
run_blast car_vs_gja "$input/target/car.wgdi.pep.fa" "$working/db/gja" & p2=$!
run_blast car_vs_opu "$input/target/car.wgdi.pep.fa" "$working/db/opu" & p3=$!
wait "$p1"; wait "$p2"; wait "$p3"

write_conf() {
    local name=$1 query=$2 subject=$3 qprefix=$4 sprefix=$5
    printf '%s\n' \
        '[collinearity]' \
        "gff1 = $input/$query/$qprefix.wgdi.gff" \
        "gff2 = $input/$subject/$sprefix.wgdi.gff" \
        "lens1 = $input/$query/$qprefix.wgdi.lens" \
        "lens2 = $input/$subject/$sprefix.wgdi.lens" \
        "blast = $working/blast/$name.tsv" \
        'blast_reverse = false' 'comparison = genomes' 'multiple = 2' \
        'process = 40' 'evalue = 1e-5' 'score = 100' 'grading = 50,40,25' \
        'mg = 40,40' 'pvalue = 0.2' 'repeat_number = 20' 'position = order' \
        "savefile = $working/collinearity/$name.tsv" \
        > "$working/config/$name.conf"
}
write_conf car_self target target car car
write_conf car_vs_gja target gardenia car gja
write_conf car_vs_opu target ophiorrhiza car opu
for name in car_self car_vs_gja car_vs_opu; do
    "$wgdi" -icl "$working/config/$name.conf" \
        >"$working/logs/$name.wgdi.stdout.log" \
        2>"$working/logs/$name.wgdi.stderr.log" &
    eval "pid_$name=$!"
done
wait "$pid_car_self"; wait "$pid_car_vs_gja"; wait "$pid_car_vs_opu"

PYTHONPATH="$code_root/src" "$python" "$ks_helper" \
    --collinearity "$working/collinearity/car_self.tsv" \
    --wgdi-gff "$input/target/car.wgdi.gff" \
    --output-dir "$working/ks/target" --shards 64
export working wgdi input
find "$working/ks/target/shards" -name 'pairs.*.tsv' -type f -print0 | sort -z | \
    "$parallel" -0 --jobs 64 --delay 0.02 --halt now,fail=1 '
      shard={}; stem={/.}; job="$working/ks/jobs/$stem"; mkdir "$job";
      printf "[ks]\ncds_file = %s\npep_file = %s\nalign_software = muscle\npairs_file = %s\nks_file = %s\n" \
        "$input/target/car.wgdi.cds.fa" "$input/target/car.wgdi.pep.fa" \
        "$shard" "$job/ks.tsv" > "$job/ks.conf";
      cd "$job"; "$wgdi" -ks ks.conf >stdout.log 2>stderr.log
    '
PYTHONPATH="$code_root/src" "$python" - "$working/ks" <<'PY'
import csv, sys
from pathlib import Path
root=Path(sys.argv[1]); rows=[]; seen=set()
fields=['id1','id2','ka_NG86','ks_NG86','ka_YN00','ks_YN00']
for path in sorted(root.glob('jobs/pairs.*/ks.tsv')):
    with path.open() as handle:
        reader=csv.DictReader(handle, delimiter='\t')
        if reader.fieldnames != fields: raise SystemExit(f'invalid YN00 header: {path}')
        for row in reader:
            key=tuple(sorted((row['id1'],row['id2'])))
            if key in seen: raise SystemExit(f'duplicate YN00 pair: {key}')
            seen.add(key); rows.append(row)
with (root/'target/ks_merged.tsv').open('x', newline='') as handle:
    writer=csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
    writer.writeheader(); writer.writerows(sorted(rows,key=lambda row:(row['id1'],row['id2'])))
PY
printf 'status\tready_for_pair_inference\n' > "$working/stage_status.tsv"
(
    cd "$working"
    find . -type f ! -name SHA256SUMS -printf '%P\0' | sort -z | \
        xargs -0 sha256sum > SHA256SUMS
    sha256sum -c SHA256SUMS >/dev/null
)
trap - ERR
mv "$working" "$output"
