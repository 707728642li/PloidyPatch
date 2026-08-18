#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then echo "usage: $0 PROJECT_ROOT" >&2; exit 2; fi
project_root=$(realpath "$1")
code_root=$project_root/code
export PYTHONPATH="$code_root/src${PYTHONPATH:+:$PYTHONPATH}"
python_bin=$project_root/envs/ploidypatch-dev/bin/python
miniprot_bin=$project_root/envs/ploidypatch-baseline/bin/miniprot
samtools_bin=$project_root/envs/ploidypatch-syngap/bin/samtools
input_root=${PLOIDYPATCH_STAGED_INPUT_ROOT:?PLOIDYPATCH_STAGED_INPUT_ROOT is required}
blind_benchmark_root=${PLOIDYPATCH_BLIND_BENCHMARK_ROOT:?PLOIDYPATCH_BLIND_BENCHMARK_ROOT is required}
protocol_root=${PLOIDYPATCH_PROTOCOL_FREEZE:?PLOIDYPATCH_PROTOCOL_FREEZE is required}
execution_root=${PLOIDYPATCH_EXECUTION_FREEZE:?PLOIDYPATCH_EXECUTION_FREEZE is required}
contract_path=${PLOIDYPATCH_HOLDOUT_CONTRACT:?PLOIDYPATCH_HOLDOUT_CONTRACT is required}
context_verifier=$code_root/scripts/verify_external_holdout_blind_context_v0.5.py
target_gz=$input_root/shared_target/Actinidia_chinensis_Red5/Actinidia_chinensis.Red5_PS1_1.69.0.dna.toplevel.fa.gz
declare -A candidate_genome=(
    [actinidia_eriantha]="$input_root/candidate_only/Actinidia_eriantha_White/Actinidia_eriantha.chr.fasta.gz"
    [actinidia_rufa]="$input_root/candidate_only/Actinidia_rufa_ARU/ARU_r1.0.pmol.fasta.gz"
)
declare -A candidate_gff=(
    [actinidia_eriantha]="$input_root/candidate_only/Actinidia_eriantha_White/Actinidia_eriantha.gff3.gz"
    [actinidia_rufa]="$input_root/candidate_only/Actinidia_rufa_ARU/ARU1.0.genes.gff.gz"
)
declare -A candidate_protein=(
    [actinidia_eriantha]="$input_root/candidate_only/Actinidia_eriantha_White/Actinidia_eriantha.pep.fasta.gz"
    [actinidia_rufa]="$input_root/candidate_only/Actinidia_rufa_ARU/ARU1.0.proteins.fasta.gz"
)
declare -A seqid_table=(
    [actinidia_eriantha]="$code_root/config/primary_seqids/actinidia_eriantha_white_v1.0.tsv"
    [actinidia_rufa]="$code_root/config/primary_seqids/actinidia_rufa_aru_r1.0.tsv"
)
result_root=$project_root/results/baselines/actinidia_v0.5/miniprot
working_root=${result_root}.working

[[ ${PLOIDYPATCH_BLIND_RUNNER:-} == 1 ]] || {
    echo "Actinidia candidate generation must run inside the frozen blind runner" >&2; exit 1;
}
[[ ${PLOIDYPATCH_NETWORK_ACCESS:-} == none ]] || {
    echo "Actinidia candidate generation requires a network-disabled namespace" >&2; exit 1;
}
for forbidden in /nas_data "$input_root/evaluator_only"; do
    [[ ! -e $forbidden ]] || { echo "forbidden blind-runner path is visible: $forbidden" >&2; exit 1; }
done
if [[ -r /proc/self/mountinfo ]] && grep -Eq '/nas_data|/evaluator_only|/target_complete|/truth_references' /proc/self/mountinfo; then
    echo "forbidden evaluator or NAS mount detected in blind namespace" >&2; exit 1
fi

verify_tree() { (cd "$1" && sha256sum -c SHA256SUMS >/dev/null); }
verify_implementation() {
    local relative=$1 expected observed
    expected=$(awk -F '\t' -v path="$relative" 'NR > 1 && $1 == path {print $3}' \
        "$execution_root/implementation_manifest.tsv")
    observed=$(sha256sum "$code_root/$relative" | awk '{print $1}')
    [[ $expected =~ ^[0-9a-f]{64}$ && $observed == "$expected" ]] || {
        echo "implementation differs from execution freeze: $relative" >&2; exit 1;
    }
}
for required in "$python_bin" "$miniprot_bin" "$samtools_bin" "$target_gz" \
    "$protocol_root/SHA256SUMS" "$execution_root/SHA256SUMS" "$contract_path" \
    "$input_root/role_manifest.tsv" "$input_root/role_contract.json" "$context_verifier" \
    "$blind_benchmark_root/perturbed.gff3" "$blind_benchmark_root/blind_manifest.json" \
    "$blind_benchmark_root/SHA256SUMS" \
    "${candidate_genome[@]}" "${candidate_gff[@]}" "${candidate_protein[@]}" \
    "${seqid_table[@]}" \
    "$code_root/config/primary_seqids/actinidia_chinensis_red5.tsv"; do
    [[ -s $required ]] || { echo "missing Actinidia blind miniprot input: $required" >&2; exit 1; }
done
verify_tree "$protocol_root"
verify_tree "$execution_root"
verify_implementation scripts/run_actinidia_miniprot_upstream_v0.5.sh
verify_implementation scripts/verify_external_holdout_blind_context_v0.5.py
verify_implementation scripts/normalize_maker_transcript_hierarchy_v0.5.py
verify_implementation scripts/synthesize_missing_transcript_exons.py
verify_implementation src/ploidypatch/normalize.py
verify_implementation src/ploidypatch/baseline.py
verify_implementation src/ploidypatch/gff_compat.py
[[ ! -e $result_root && ! -e $working_root ]] || {
    echo "refusing to overwrite Actinidia miniprot upstream" >&2; exit 1;
}
mkdir -p "$working_root"/{normalized,provider_compatible,target,reference,index,raw,logs,freeze}
cp /proc/self/mountinfo "$working_root/freeze/blind_runner.mountinfo"
PYTHONPATH="$code_root/src" "$python_bin" "$context_verifier" \
    --input-root "$input_root" --contract "$contract_path" \
    --protocol-freeze "$protocol_root" --execution-freeze "$execution_root" \
    --blind-benchmark-root "$blind_benchmark_root" \
    --expected-holdout-id actinidia_red5_v0.5 --expected-primary-chromosomes 29 \
    --output-json "$working_root/freeze/blind_context.json"
{
    printf 'field\tvalue\ncode_commit\t%s\n' "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'stage\ttruth_blind_candidate_upstream_v0.5\nblind_runner\ttrue\n'
    printf 'candidate_truth_access\tfalse\ntarget_complete_annotation_access\tfalse\n'
    printf 'evaluator_reference_access\tfalse\nnas_data_access\tfalse\nnetwork_access\tfalse\n'
    printf 'reference_species\tActinidia_eriantha_White,Actinidia_rufa_ARU\nthreads\t32\n'
    printf 'within_method_reference_vote_count\t1\n'
    printf 'protocol_freeze_sha256sums_sha256\t%s\n' \
        "$(sha256sum "$protocol_root/SHA256SUMS" | awk '{print $1}')"
    printf 'execution_freeze_sha256sums_sha256\t%s\n' \
        "$(sha256sum "$execution_root/SHA256SUMS" | awk '{print $1}')"
} > "$working_root/run_contract.tsv"

target_bundle=$working_root/normalized/target
mkdir -p "$target_bundle"
"$python_bin" - "$target_gz" "$code_root/config/primary_seqids/actinidia_chinensis_red5.tsv" \
    "$target_bundle/primary_chromosomes.genome.fa" <<'PY'
import csv
import sys
from pathlib import Path

from ploidypatch.io import iter_fasta

source, table, output = sys.argv[1:]
with open(table, encoding="utf-8", newline="") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))
allowed = {row["seqid"] for row in rows}
if len(allowed) != 29:
    raise SystemExit("Actinidia primary-seqid freeze must contain exactly 29 seqids")
observed = set()
with Path(output).open("xb") as outgoing:
    for seqid, header, sequence in iter_fasta(source):
        if seqid not in allowed:
            continue
        if seqid in observed:
            raise SystemExit(f"duplicate primary FASTA seqid: {seqid}")
        if not sequence:
            raise SystemExit(f"empty primary FASTA seqid: {seqid}")
        observed.add(seqid)
        outgoing.write(f">{header}\n".encode("utf-8"))
        encoded = sequence.encode("ascii")
        for index in range(0, len(encoded), 60):
            outgoing.write(encoded[index:index + 60] + b"\n")
if observed != allowed:
    raise SystemExit(f"primary FASTA mismatch: missing={sorted(allowed - observed)}")
PY
"$samtools_bin" faidx "$target_bundle/primary_chromosomes.genome.fa"
"$python_bin" - "$target_gz" "$code_root/config/primary_seqids/actinidia_chinensis_red5.tsv" \
    "$target_bundle/primary_chromosomes.genome.fa" \
    "$target_bundle/primary_chromosomes.genome.fa.fai" "$target_bundle/manifest.json" <<'PY'
import csv
import hashlib
import json
import sys
from pathlib import Path

source, table, genome, fai, output = map(Path, sys.argv[1:])
with table.open(encoding="utf-8", newline="") as handle:
    seqids = [row["seqid"] for row in csv.DictReader(handle, delimiter="\t")]

def sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(8 * 1024 * 1024):
            digest.update(block)
    return digest.hexdigest()

manifest = {
    "schema_version": "ploidypatch.blind_target_genome_normalization.v0.5",
    "truth_access": False,
    "complete_target_annotation_access": False,
    "selection": "exact_frozen_primary_seqid_table",
    "primary_seqids": seqids,
    "inputs": {
        "source_genome_sha256": sha(source),
        "primary_seqid_table_sha256": sha(table),
    },
    "outputs": {
        "genome_sha256": sha(genome),
        "fai_sha256": sha(fai),
        "chromosomes": len(seqids),
    },
}
with output.open("x", encoding="utf-8") as handle:
    json.dump(manifest, handle, indent=2, sort_keys=True)
    handle.write("\n")
PY

"$python_bin" - "$target_bundle/primary_chromosomes.genome.fa" \
    "$blind_benchmark_root/blind_manifest.json" <<'PY'
import hashlib
import json
import sys
from pathlib import Path

genome, manifest_path = map(Path, sys.argv[1:])
digest = hashlib.sha256()
with genome.open("rb") as handle:
    while block := handle.read(8 * 1024 * 1024):
        digest.update(block)
manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
if manifest.get("target_genome", {}).get("sha256") != digest.hexdigest():
    raise SystemExit("normalized Red5 primary genome differs from sealed blind benchmark")
PY

cd "$code_root"
for species in actinidia_eriantha actinidia_rufa; do
    hierarchy_mode=shared_gene_transcript_id
    [[ $species != actinidia_rufa ]] || hierarchy_mode=gene_with_direct_children
    hierarchy_gff=$working_root/provider_compatible/$species.gff3
    "$python_bin" "$code_root/scripts/normalize_maker_transcript_hierarchy_v0.5.py" \
        --input-gff "${candidate_gff[$species]}" --output-gff "$hierarchy_gff" \
        --mode "$hierarchy_mode" \
        > "$working_root/logs/hierarchy.$species.stdout.json" \
        2> "$working_root/logs/hierarchy.$species.stderr.log"
    "$python_bin" - "$hierarchy_gff.manifest.json" "$hierarchy_mode" <<'PY'
import json
import sys

report = json.load(open(sys.argv[1], encoding="utf-8"))
cds = report.get("cds_rows_sha256", {})
if (
    report.get("schema_version")
    != "ploidypatch.maker_transcript_hierarchy_compat.v0.5"
    or report.get("mode") != sys.argv[2]
    or report.get("coordinate_or_cds_changes") is not False
    or report.get("source_genes", 0) < 1
    or report.get("source_genes") != report.get("output_transcripts")
    or cds.get("input") != cds.get("output")
):
    raise SystemExit("provider hierarchy adapter violated its frozen CDS invariant")
PY
    "$python_bin" -m ploidypatch.cli normalize primary-annotation \
        --gff "$hierarchy_gff" --genome "${candidate_genome[$species]}" \
        --primary-seqid-table "${seqid_table[$species]}" \
        --output-dir "$working_root/normalized/$species" \
        > "$working_root/logs/normalize.$species.stdout.json" \
        2> "$working_root/logs/normalize.$species.stderr.log"
    "$python_bin" "$code_root/scripts/synthesize_missing_transcript_exons.py" \
        --input-gff "$working_root/normalized/$species/primary_chromosomes.gff3" \
        --output-gff "$working_root/normalized/$species/primary_chromosomes.lifton.gff3" \
        --repair-parent-bounds \
        > "$working_root/logs/lifton_gff_compat.$species.stdout.json" \
        2> "$working_root/logs/lifton_gff_compat.$species.stderr.log"
    "$python_bin" - \
        "$working_root/normalized/$species/primary_chromosomes.lifton.gff3.manifest.json" <<'PY'
import json
import sys

report = json.load(open(sys.argv[1], encoding="utf-8"))
counts = report.get("counts", {})
repairs = report.get("parent_bound_repairs", [])
if (
    report.get("schema_version")
    != "ploidypatch.missing_transcript_exon_compat.v2"
    or report.get("repair_parent_bounds") is not True
    or report.get("child_coordinate_or_cds_changes") is not False
    or counts.get("unresolved_transcripts") != 0
    or counts.get("input_cds_records") != counts.get("output_cds_records")
    or counts.get("parent_bounds_repaired") != len(repairs)
    or counts.get("parent_bounds_repaired")
    != counts.get("transcript_parent_bounds_repaired", 0)
    + counts.get("gene_parent_bounds_repaired", 0)
    or report.get("cds_rows_sha256", {}).get("input")
    != report.get("cds_rows_sha256", {}).get("output")
    or any(item.get("feature_type") not in {"gene", "transcript"} for item in repairs)
    or any(
        item.get("output_start", 0) > item.get("input_start", 0)
        or item.get("output_end", 0) < item.get("input_end", 0)
        for item in repairs
    )
):
    raise SystemExit("LiftOn GFF compatibility adapter violated its frozen invariant")
PY
done

"$python_bin" -m ploidypatch.cli baseline prepare-proteins \
    --protein "actinidia_eriantha=${candidate_protein[actinidia_eriantha]}" \
    --protein "actinidia_rufa=${candidate_protein[actinidia_rufa]}" \
    --output-fasta "$working_root/reference/actinidia_candidate_refs.pep.fa" \
    --output-map "$working_root/reference/actinidia_candidate_refs.map.tsv" \
    > "$working_root/reference/prepare.stdout.json" \
    2> "$working_root/reference/prepare.stderr.log"
/usr/bin/time -v -o "$working_root/logs/index.time.txt" \
    "$miniprot_bin" -t32 -d "$working_root/index/red5.mpi" \
    "$target_bundle/primary_chromosomes.genome.fa" \
    > "$working_root/logs/index.stdout.log" 2> "$working_root/logs/index.stderr.log"
/usr/bin/time -v -o "$working_root/logs/projection.time.txt" \
    "$miniprot_bin" -I -t32 --gff-only -N4 --outn 4 --outs 0.8 --outc 0.5 \
    "$working_root/index/red5.mpi" "$working_root/reference/actinidia_candidate_refs.pep.fa" \
    > "$working_root/raw/miniprot.gff3" 2> "$working_root/logs/projection.stderr.log"
grep -q '^##gff-version' "$working_root/raw/miniprot.gff3"
for output in normalized/target/primary_chromosomes.genome.fa \
              normalized/target/primary_chromosomes.genome.fa.fai \
              provider_compatible/actinidia_eriantha.gff3 \
              provider_compatible/actinidia_eriantha.gff3.manifest.json \
              provider_compatible/actinidia_rufa.gff3 \
              provider_compatible/actinidia_rufa.gff3.manifest.json \
              normalized/actinidia_eriantha/primary_chromosomes.genome.fa \
              normalized/actinidia_eriantha/primary_chromosomes.gff3 \
              normalized/actinidia_eriantha/primary_chromosomes.lifton.gff3 \
              normalized/actinidia_eriantha/primary_chromosomes.lifton.gff3.manifest.json \
              normalized/actinidia_rufa/primary_chromosomes.genome.fa \
              normalized/actinidia_rufa/primary_chromosomes.gff3 \
              normalized/actinidia_rufa/primary_chromosomes.lifton.gff3 \
              normalized/actinidia_rufa/primary_chromosomes.lifton.gff3.manifest.json \
              reference/actinidia_candidate_refs.pep.fa \
              reference/actinidia_candidate_refs.map.tsv raw/miniprot.gff3; do
    [[ -s $working_root/$output ]] || { echo "missing Actinidia miniprot output: $output" >&2; exit 1; }
done
{
    printf 'role\tbytes\tsha256\tpath\n'
    for entry in "target_genome:$target_gz" \
        "actinidia_eriantha_genome:${candidate_genome[actinidia_eriantha]}" \
        "actinidia_eriantha_gff3:${candidate_gff[actinidia_eriantha]}" \
        "actinidia_eriantha_protein:${candidate_protein[actinidia_eriantha]}" \
        "actinidia_rufa_genome:${candidate_genome[actinidia_rufa]}" \
        "actinidia_rufa_gff3:${candidate_gff[actinidia_rufa]}" \
        "actinidia_rufa_protein:${candidate_protein[actinidia_rufa]}"; do
        role=${entry%%:*}; path=${entry#*:}
        printf '%s\t%s\t%s\t%s\n' "$role" "$(stat -Lc %s "$path")" \
            "$(sha256sum "$path" | awk '{print $1}')" "$path"
    done
} > "$working_root/input_manifest.tsv"
(
    cd "$working_root"
    { find normalized provider_compatible reference raw freeze -type f -print0; \
      printf '%s\0' run_contract.tsv input_manifest.tsv; } \
        | sort -z | xargs -0 sha256sum > SHA256SUMS
    sha256sum -c SHA256SUMS >/dev/null
)
du -sb "$working_root" > "$working_root/disk_bytes.txt"
mv "$working_root" "$result_root"
printf 'Actinidia truth-blind miniprot upstream frozen: %s\n' "$result_root"
