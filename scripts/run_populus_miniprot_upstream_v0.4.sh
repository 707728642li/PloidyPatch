#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then echo "usage: $0 PROJECT_ROOT" >&2; exit 2; fi
project_root=$(realpath "$1")
code_root=$project_root/code
export PYTHONPATH="$code_root/src${PYTHONPATH:+:$PYTHONPATH}"
python_bin=$project_root/envs/ploidypatch-dev/bin/python
miniprot_bin=$project_root/envs/ploidypatch-baseline/bin/miniprot
samtools_bin=$project_root/envs/ploidypatch-syngap/bin/samtools
input_root=$project_root/data/derived/external_inputs/populus_v0.4
protocol_root=$project_root/results/protocol_freezes/populus_external_v0.4
execution_root=$project_root/results/protocol_freezes/populus_external_v0.4_execution
preflight=$protocol_root/preflight_input_manifest.tsv
target_gz=$input_root/shared_target/Populus_trichocarpa/Ptrichocarpa_533_v4.0.fa.gz
declare -A candidate_genome=(
    [salix_purpurea]="$input_root/candidate_only/Salix_purpurea/Spurpurea_519_v5.0.fa.gz"
    [salix_suchowensis]="$input_root/candidate_only/Salix_suchowensis/GCA_017552425.1_ASM1755242v1_genomic.fna.gz"
)
declare -A candidate_gff=(
    [salix_purpurea]="$input_root/candidate_only/Salix_purpurea/Spurpurea_519_v5.1.gene.gff3.gz"
    [salix_suchowensis]="$input_root/candidate_only/Salix_suchowensis/GCA_017552425.1_ASM1755242v1_genomic.gff.gz"
)
declare -A candidate_protein=(
    [salix_purpurea]="$input_root/candidate_only/Salix_purpurea/Spurpurea_519_v5.1.protein_primaryTranscriptOnly.fa.gz"
    [salix_suchowensis]="$input_root/candidate_only/Salix_suchowensis/GCA_017552425.1_ASM1755242v1_protein.faa.gz"
)
declare -A seqid_table=(
    [salix_purpurea]="$code_root/config/primary_seqids/salix_purpurea_v5.0.tsv"
    [salix_suchowensis]="$code_root/config/primary_seqids/salix_suchowensis_GCA_017552425.1.tsv"
)
result_root=$project_root/results/baselines/populus_v0.4/miniprot
working_root=${result_root}.working

[[ ${PLOIDYPATCH_BLIND_RUNNER:-} == 1 ]] || {
    echo "Populus candidate generation must run inside the frozen blind runner" >&2; exit 1;
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
verify_preflight() {
    local role=$1 species=$2 artifact=$3 path=$4 expected observed
    expected=$(awk -F '\t' -v role="$role" -v species="$species" -v artifact="$artifact" \
        'NR > 1 && $1 == role && $2 == species && $4 == artifact {print $6}' "$preflight")
    observed=$(sha256sum "$path" | awk '{print $1}')
    [[ $expected =~ ^[0-9a-f]{64}$ && $observed == "$expected" ]] || {
        echo "staged blind input differs from protocol preflight: $species/$artifact" >&2; exit 1;
    }
}

for required in "$python_bin" "$miniprot_bin" "$samtools_bin" "$target_gz" \
    "$protocol_root/SHA256SUMS" "$execution_root/SHA256SUMS" "$preflight" \
    "${candidate_genome[@]}" "${candidate_gff[@]}" "${candidate_protein[@]}" \
    "${seqid_table[@]}" \
    "$code_root/config/primary_seqids/populus_trichocarpa_v4.0.tsv"; do
    [[ -s $required ]] || { echo "missing Populus blind miniprot input: $required" >&2; exit 1; }
done
verify_tree "$protocol_root"
verify_tree "$execution_root"
verify_implementation scripts/run_populus_miniprot_upstream_v0.4.sh
verify_implementation scripts/synthesize_missing_transcript_exons.py
verify_implementation src/ploidypatch/normalize.py
verify_implementation src/ploidypatch/baseline.py
verify_implementation src/ploidypatch/gff_compat.py
verify_preflight target Populus_trichocarpa genome "$target_gz"
for species in salix_purpurea salix_suchowensis; do
    staged_species=Salix_purpurea; [[ $species == salix_suchowensis ]] && staged_species=Salix_suchowensis
    verify_preflight candidate_reference "$staged_species" genome "${candidate_genome[$species]}"
    verify_preflight candidate_reference "$staged_species" gff3 "${candidate_gff[$species]}"
    verify_preflight candidate_reference "$staged_species" protein "${candidate_protein[$species]}"
done
[[ ! -e $result_root && ! -e $working_root ]] || {
    echo "refusing to overwrite Populus miniprot upstream" >&2; exit 1;
}
mkdir -p "$working_root"/{normalized,target,reference,index,raw,logs,freeze}
cp /proc/self/mountinfo "$working_root/freeze/blind_runner.mountinfo"
{
    printf 'field\tvalue\ncode_commit\t%s\n' "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'stage\ttruth_blind_candidate_upstream_v0.4\nblind_runner\ttrue\n'
    printf 'candidate_truth_access\tfalse\ntarget_complete_annotation_access\tfalse\n'
    printf 'evaluator_reference_access\tfalse\nnas_data_access\tfalse\nnetwork_access\tfalse\n'
    printf 'reference_species\tSalix_purpurea,Salix_suchowensis\nthreads\t32\n'
    printf 'within_method_reference_vote_count\t1\n'
    printf 'protocol_freeze_sha256sums_sha256\t%s\n' \
        "$(sha256sum "$protocol_root/SHA256SUMS" | awk '{print $1}')"
    printf 'execution_freeze_sha256sums_sha256\t%s\n' \
        "$(sha256sum "$execution_root/SHA256SUMS" | awk '{print $1}')"
} > "$working_root/run_contract.tsv"

target_bundle=$working_root/normalized/target
mkdir -p "$target_bundle"
"$python_bin" - "$target_gz" "$code_root/config/primary_seqids/populus_trichocarpa_v4.0.tsv" \
    "$target_bundle/primary_chromosomes.genome.fa" <<'PY'
import csv
import sys
from pathlib import Path

from ploidypatch.io import iter_fasta

source, table, output = sys.argv[1:]
with open(table, encoding="utf-8", newline="") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))
allowed = {row["seqid"] for row in rows}
if len(allowed) != 19:
    raise SystemExit("Populus primary-seqid freeze must contain exactly 19 seqids")
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
"$python_bin" - "$target_gz" "$code_root/config/primary_seqids/populus_trichocarpa_v4.0.tsv" \
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
    "schema_version": "ploidypatch.blind_target_genome_normalization.v0.4",
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

cd "$code_root"
for species in salix_purpurea salix_suchowensis; do
    "$python_bin" -m ploidypatch.cli normalize primary-annotation \
        --gff "${candidate_gff[$species]}" --genome "${candidate_genome[$species]}" \
        --primary-seqid-table "${seqid_table[$species]}" \
        --output-dir "$working_root/normalized/$species" \
        > "$working_root/logs/normalize.$species.stdout.json" \
        2> "$working_root/logs/normalize.$species.stderr.log"
    "$python_bin" "$code_root/scripts/synthesize_missing_transcript_exons.py" \
        --input-gff "$working_root/normalized/$species/primary_chromosomes.gff3" \
        --output-gff "$working_root/normalized/$species/primary_chromosomes.lifton.gff3" \
        > "$working_root/logs/lifton_gff_compat.$species.stdout.json" \
        2> "$working_root/logs/lifton_gff_compat.$species.stderr.log"
    "$python_bin" - \
        "$working_root/normalized/$species/primary_chromosomes.lifton.gff3.manifest.json" <<'PY'
import json
import sys

report = json.load(open(sys.argv[1], encoding="utf-8"))
counts = report.get("counts", {})
if (
    report.get("schema_version")
    != "ploidypatch.missing_transcript_exon_compat.v1"
    or report.get("coordinate_or_cds_changes") is not False
    or counts.get("unresolved_transcripts") != 0
    or counts.get("input_cds_records") != counts.get("output_cds_records")
    or report.get("cds_rows_sha256", {}).get("input")
    != report.get("cds_rows_sha256", {}).get("output")
):
    raise SystemExit("LiftOn GFF compatibility adapter violated its frozen invariant")
PY
done

"$python_bin" -m ploidypatch.cli baseline prepare-proteins \
    --protein "salix_purpurea=${candidate_protein[salix_purpurea]}" \
    --protein "salix_suchowensis=${candidate_protein[salix_suchowensis]}" \
    --output-fasta "$working_root/reference/populus_candidate_refs.pep.fa" \
    --output-map "$working_root/reference/populus_candidate_refs.map.tsv" \
    > "$working_root/reference/prepare.stdout.json" \
    2> "$working_root/reference/prepare.stderr.log"
/usr/bin/time -v -o "$working_root/logs/index.time.txt" \
    "$miniprot_bin" -t32 -d "$working_root/index/ptr_v4.mpi" \
    "$target_bundle/primary_chromosomes.genome.fa" \
    > "$working_root/logs/index.stdout.log" 2> "$working_root/logs/index.stderr.log"
/usr/bin/time -v -o "$working_root/logs/projection.time.txt" \
    "$miniprot_bin" -I -t32 --gff-only -N4 --outn 4 --outs 0.8 --outc 0.5 \
    "$working_root/index/ptr_v4.mpi" "$working_root/reference/populus_candidate_refs.pep.fa" \
    > "$working_root/raw/miniprot.gff3" 2> "$working_root/logs/projection.stderr.log"
grep -q '^##gff-version' "$working_root/raw/miniprot.gff3"
for output in normalized/target/primary_chromosomes.genome.fa \
              normalized/target/primary_chromosomes.genome.fa.fai \
              normalized/salix_purpurea/primary_chromosomes.genome.fa \
              normalized/salix_purpurea/primary_chromosomes.gff3 \
              normalized/salix_purpurea/primary_chromosomes.lifton.gff3 \
              normalized/salix_purpurea/primary_chromosomes.lifton.gff3.manifest.json \
              normalized/salix_suchowensis/primary_chromosomes.genome.fa \
              normalized/salix_suchowensis/primary_chromosomes.gff3 \
              normalized/salix_suchowensis/primary_chromosomes.lifton.gff3 \
              normalized/salix_suchowensis/primary_chromosomes.lifton.gff3.manifest.json \
              reference/populus_candidate_refs.pep.fa \
              reference/populus_candidate_refs.map.tsv raw/miniprot.gff3; do
    [[ -s $working_root/$output ]] || { echo "missing Populus miniprot output: $output" >&2; exit 1; }
done
{
    printf 'role\tbytes\tsha256\tpath\n'
    for entry in "target_genome:$target_gz" \
        "salix_purpurea_genome:${candidate_genome[salix_purpurea]}" \
        "salix_purpurea_gff3:${candidate_gff[salix_purpurea]}" \
        "salix_purpurea_protein:${candidate_protein[salix_purpurea]}" \
        "salix_suchowensis_genome:${candidate_genome[salix_suchowensis]}" \
        "salix_suchowensis_gff3:${candidate_gff[salix_suchowensis]}" \
        "salix_suchowensis_protein:${candidate_protein[salix_suchowensis]}"; do
        role=${entry%%:*}; path=${entry#*:}
        printf '%s\t%s\t%s\t%s\n' "$role" "$(stat -Lc %s "$path")" \
            "$(sha256sum "$path" | awk '{print $1}')" "$path"
    done
} > "$working_root/input_manifest.tsv"
(
    cd "$working_root"
    { find normalized reference raw freeze -type f -print0; \
      printf '%s\0' run_contract.tsv input_manifest.tsv; } \
        | sort -z | xargs -0 sha256sum > SHA256SUMS
    sha256sum -c SHA256SUMS >/dev/null
)
du -sb "$working_root" > "$working_root/disk_bytes.txt"
mv "$working_root" "$result_root"
printf 'Populus truth-blind miniprot upstream frozen: %s\n' "$result_root"
