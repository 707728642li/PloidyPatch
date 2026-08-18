#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then echo "usage: $0 PROJECT_ROOT" >&2; exit 2; fi
project_root=$(realpath "$1")
parallel_bin=/data/codexli/software/conda/miniforge3/bin/parallel
input_root=$project_root/data/derived/holdout_inputs/populus_v0.4
protocol_root=$project_root/results/protocol_freezes/populus_external_v0.4
execution_root=$project_root/results/protocol_freezes/populus_external_v0.4_execution
code_root=$execution_root/source
export PYTHONPATH="$code_root/src${PYTHONPATH:+:$PYTHONPATH}"
environment_bindings=$execution_root/environment_bindings.tsv
[[ -s $environment_bindings ]] || { echo "missing frozen environment bindings" >&2; exit 1; }
dev_prefix=$(awk -F '\t' '$1 == "ploidypatch-dev" {print $2}' "$environment_bindings")
[[ $dev_prefix == /* ]] || { echo "invalid frozen ploidypatch-dev binding" >&2; exit 1; }
python_bin=$dev_prefix/bin/python
result_root=$project_root/data/derived/external_inputs/populus_v0.4
working_root=${result_root}.working
self_relative=scripts/prepare_populus_external_normalized_inputs_v0.4.sh

verify_implementation() {
    local relative=$1 manifest=$execution_root/implementation_manifest.tsv
    local rows=()
    mapfile -t rows < <(awk -F '\t' -v path="$relative" '$1 == path {print $2 "\t" $3}' "$manifest")
    [[ ${#rows[@]} -eq 1 ]] || { echo "execution freeze has no unique row for $relative" >&2; return 1; }
    local expected_bytes expected_sha
    IFS=$'\t' read -r expected_bytes expected_sha <<< "${rows[0]}"
    [[ $expected_bytes =~ ^[0-9]+$ && $expected_sha =~ ^[0-9a-f]{64}$ ]] || {
        echo "malformed execution implementation row for $relative" >&2; return 1;
    }
    [[ $(stat -Lc %s "$code_root/$relative") == "$expected_bytes" \
        && $(sha256sum "$code_root/$relative" | awk '{print $1}') == "$expected_sha" ]] || {
        echo "implementation differs from execution freeze: $relative" >&2; return 1;
    }
}

implementation_dependencies=(
    "$self_relative"
    src/ploidypatch/cli.py
    src/ploidypatch/normalize.py
)

for required in "$python_bin" "$parallel_bin" "$input_root/role_manifest.tsv" \
                "$input_root/role_contract.json" "$input_root/SHA256SUMS" \
                "$protocol_root/SHA256SUMS" "$protocol_root/preflight_input_manifest.tsv" \
                "$execution_root/SHA256SUMS" \
                "$execution_root/implementation_manifest.tsv" \
                "${implementation_dependencies[@]/#/$code_root/}" \
                "$code_root/config/primary_seqids/populus_trichocarpa_v4.0.tsv" \
                "$code_root/config/primary_seqids/manihot_esculenta_v6.tsv" \
                "$code_root/config/primary_seqids/ricinus_communis_wild_castor.tsv"; do
    [[ -s $required ]] || { echo "missing Populus normalization prerequisite: $required" >&2; exit 1; }
done
(cd "$protocol_root" && sha256sum -c SHA256SUMS >/dev/null)
(cd "$execution_root" && sha256sum -c SHA256SUMS >/dev/null)
(cd "$input_root" && sha256sum -c SHA256SUMS >/dev/null)
for relative in "${implementation_dependencies[@]}"; do verify_implementation "$relative"; done
[[ ! -e $result_root && ! -e $working_root ]] || {
    echo "refusing to overwrite normalized Populus evaluator inputs" >&2; exit 1;
}
mkdir -p "$working_root"/{provider_compatible,normalized}
record_invalid() {
    local status=$?
    if [[ -d $working_root ]]; then
        printf 'field\tvalue\nformal_status\tinvalid_run\nstage\tnormalization\nexit_status\t%s\n' \
            "$status" > "$working_root/invalid_run.tsv" || true
    fi
    exit "$status"
}
trap record_invalid ERR

"$python_bin" - "$input_root" "$working_root" <<'PY'
import csv
import hashlib
import shutil
import sys
from pathlib import Path

source, destination = map(Path, sys.argv[1:])
def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()
with (source / "role_manifest.tsv").open(encoding="utf-8", newline="") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))
if not rows:
    raise ValueError("staged role manifest is empty")
seen = set()
for row in rows:
    relative = Path(row["staged_relative_path"])
    if relative.is_absolute() or ".." in relative.parts or relative in seen:
        raise ValueError(f"unsafe or duplicate staged role path: {relative}")
    seen.add(relative)
    incoming = source / relative
    outgoing = destination / relative
    if (
        incoming.is_symlink()
        or not incoming.is_file()
        or incoming.stat().st_size != int(row["bytes"])
        or sha256(incoming) != row["sha256"]
    ):
        raise ValueError(f"staged role artifact differs before copy: {incoming}")
    outgoing.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(incoming, outgoing)
    if outgoing.is_symlink() or outgoing.stat().st_size != int(row["bytes"]) or sha256(outgoing) != row["sha256"]:
        raise IOError(f"role-preserving copy verification failed: {outgoing}")
for name in ("role_manifest.tsv", "role_contract.json"):
    incoming, outgoing = source / name, destination / name
    if incoming.is_symlink() or not incoming.is_file():
        raise ValueError(f"staged role metadata is unsafe: {incoming}")
    shutil.copyfile(incoming, outgoing)
shutil.copyfile(source / "SHA256SUMS", destination / "staging_SHA256SUMS")
PY
(
    cd "$working_root"
    sha256sum -c staging_SHA256SUMS >/dev/null
)

resolution=$working_root/input_resolution.tsv
"$python_bin" - "$input_root" "$protocol_root" "$resolution" <<'PY'
import csv
import hashlib
import json
import sys
from pathlib import Path

input_root, protocol_root, output = map(Path, sys.argv[1:])

def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()

with (input_root / "role_contract.json").open(encoding="utf-8") as handle:
    contract = json.load(handle)
if (
    contract.get("schema_version") != "ploidypatch.populus_external_input_stage.v0.4"
    or contract.get("wgd_pairs_enumerated") is not False
    or contract.get("candidate_counts_computed") is not False
    or contract.get("truth_labels_accessed") is not False
    or contract.get("protocol_SHA256SUMS_sha256") != sha256(protocol_root / "SHA256SUMS")
):
    raise ValueError("staged Populus role contract violates frozen pre-truth scope")
with (input_root / "role_manifest.tsv").open(encoding="utf-8", newline="") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))
with (protocol_root / "preflight_input_manifest.tsv").open(
    encoding="utf-8", newline=""
) as handle:
    preflight_rows = list(csv.DictReader(handle, delimiter="\t"))
preflight = {
    (row["role"], row["species_id"], row["artifact"]): row
    for row in preflight_rows
}
if len(preflight) != len(preflight_rows):
    raise ValueError("protocol preflight contains duplicate role/species/artifact rows")
expected = {
    ("target", "Populus_trichocarpa", artifact)
    for artifact in ("genome", "gff3", "protein")
} | {
    ("evaluator_reference", species, artifact)
    for species in ("Manihot_esculenta", "Ricinus_communis")
    for artifact in ("genome", "gff3", "protein")
}
selected = {}
for row in rows:
    key = (row["role"], row["species_id"], row["artifact"])
    if key not in expected:
        continue
    if key in selected:
        raise ValueError(f"duplicate staged evaluator artifact: {key}")
    relative = Path(row["staged_relative_path"])
    if relative.is_absolute() or ".." in relative.parts:
        raise ValueError(f"unsafe staged evaluator path: {relative}")
    allowed = (
        Path("shared_target") if key[0] == "target" and key[2] == "genome"
        else Path("evaluator_only/target_complete") if key[0] == "target"
        else Path("evaluator_only/truth_references")
    )
    if relative.parts[: len(allowed.parts)] != allowed.parts:
        raise ValueError(f"role path does not match evaluator contract: {relative}")
    path = input_root / relative
    if (
        not path.is_file()
        or path.stat().st_size != int(row["bytes"])
        or sha256(path) != row["sha256"]
        or row["staged_sha256"] != row["sha256"]
    ):
        raise ValueError(f"staged evaluator artifact differs from role manifest: {path}")
    frozen = preflight.get(key)
    if (
        frozen is None
        or frozen["bytes"] != row["bytes"]
        or frozen["sha256"] != row["sha256"]
        or frozen["release"] != row["release"]
    ):
        raise ValueError(f"staged evaluator role differs from protocol preflight: {key}")
    selected[key] = path
if set(selected) != expected:
    raise ValueError(f"staged evaluator artifacts incomplete: {expected - set(selected)}")

names = (
    ("target_populus", "target", "Populus_trichocarpa"),
    ("evaluator_manihot", "evaluator_reference", "Manihot_esculenta"),
    ("evaluator_ricinus", "evaluator_reference", "Ricinus_communis"),
)
with output.open("x", encoding="utf-8", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(("bundle", "role", "species_id", "genome", "gff3", "protein"))
    for bundle, role, species in names:
        writer.writerow(
            (
                bundle,
                role,
                species,
                selected[(role, species, "genome")],
                selected[(role, species, "gff3")],
                selected[(role, species, "protein")],
            )
        )
PY

{
    printf 'field\tvalue\n'
    printf 'code_commit\t%s\n' "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'access\tevaluator_only_after_protocol_and_execution_freeze\n'
    printf 'stage\trole_preserving_strict_pre_pair_normalization_only\n'
    printf 'candidate_reference_access\tbyte_verified_copy_only_without_analysis\n'
    printf 'hidden_pair_enumeration\tfalse\n'
    printf 'hidden_event_generation\tfalse\n'
    printf 'truth_label_generation\tfalse\n'
    printf 'provider_gff_compatibility\tstrict_zero_repair\n'
    printf 'protocol_freeze_sha256sums_sha256\t%s\n' \
        "$(sha256sum "$protocol_root/SHA256SUMS" | awk '{print $1}')"
    printf 'execution_freeze_sha256sums_sha256\t%s\n' \
        "$(sha256sum "$execution_root/SHA256SUMS" | awk '{print $1}')"
} > "$working_root/run_contract.tsv"
printf 'bundle\trepair_unescaped_note_semicolons\tdrop_invalid_intron_intervals\tstrip_embedded_fasta\texpected_note_repairs\texpected_dropped_introns\texpected_fasta_directives\n' \
    > "$working_root/provider_compatibility_policy.tsv"

tail -n +2 "$resolution" | while IFS=$'\t' read -r bundle role species genome gff protein; do
    printf '%s\tfalse\tfalse\tfalse\t0\t0\t0\n' "$bundle" \
        >> "$working_root/provider_compatibility_policy.tsv"
    /usr/bin/time -v -o "$working_root/$bundle.provider_compatibility.resource.time.txt" \
        "$python_bin" -m ploidypatch.cli normalize provider-gff3 \
        --gff "$gff" --output-dir "$working_root/provider_compatible/$bundle" \
        > "$working_root/$bundle.provider_compatibility.stdout.json" \
        2> "$working_root/$bundle.provider_compatibility.stderr.log"
    "$python_bin" - "$working_root/provider_compatible/$bundle/manifest.json" <<'PY'
import json
import sys
with open(sys.argv[1], encoding="utf-8") as handle:
    observed = json.load(handle)["observed"]
counts = (
    observed.get("repaired_note_records", 0),
    observed.get("dropped_invalid_intron_intervals", 0),
    observed.get("stripped_embedded_fasta_directives", 0),
)
if counts != (0, 0, 0):
    raise ValueError(f"provider GFF requires an unfrozen compatibility repair: {counts}")
PY
done

commands=$working_root/normalization_commands.tsv
printf 'bundle\tgff3\tgenome\tprimary_seqids\texpected_chromosomes\texpected_genes\texpected_transcripts\texpected_cds\n' > "$commands"
while IFS=$'\t' read -r bundle role species genome gff protein; do
    [[ $bundle == bundle ]] && continue
    case $bundle in
        target_populus)
            primary=$code_root/config/primary_seqids/populus_trichocarpa_v4.0.tsv
            expected_chromosomes=19; expected_genes=34488
            expected_transcripts=52085; expected_cds=308640 ;;
        evaluator_manihot)
            primary=$code_root/config/primary_seqids/manihot_esculenta_v6.tsv
            expected_chromosomes=18; expected_genes=31901
            expected_transcripts=40097; expected_cds=214523 ;;
        evaluator_ricinus)
            primary=$code_root/config/primary_seqids/ricinus_communis_wild_castor.tsv
            expected_chromosomes=10; expected_genes=24585
            expected_transcripts=39624; expected_cds=255805 ;;
        *) echo "unexpected evaluator bundle: $bundle" >&2; exit 1 ;;
    esac
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$bundle" \
        "$working_root/provider_compatible/$bundle/sanitized.gff3" \
        "$genome" "$primary" "$expected_chromosomes" "$expected_genes" \
        "$expected_transcripts" "$expected_cds" >> "$commands"
done < "$resolution"

cd "$code_root"
"$parallel_bin" --colsep '\t' --delay 1 --jobs 3 --header : \
    --joblog "$working_root/parallel.joblog.tsv" \
    '/usr/bin/time -v -o '"$working_root"'/{bundle}.resource.time.txt '"$python_bin"' -m ploidypatch.cli normalize primary-annotation --gff {gff3} --genome {genome} --primary-seqid-table {primary_seqids} --output-dir '"$working_root"'/normalized/{bundle} > '"$working_root"'/{bundle}.stdout.json 2> '"$working_root"'/{bundle}.stderr.log' \
    :::: "$commands"

printf 'bundle\tchromosomes\tgenes\ttranscripts\tcds_segments\n' \
    > "$working_root/annotation_counts.tsv"
while IFS=$'\t' read -r bundle gff genome primary expected_chromosomes expected_genes expected_transcripts expected_cds; do
    [[ $bundle == bundle ]] && continue
    output=$working_root/normalized/$bundle
    for required in "$output/primary_chromosomes.genome.fa" \
                    "$output/primary_chromosomes.genome.fa.fai" \
                    "$output/primary_chromosomes.gff3" "$output/manifest.json"; do
        [[ -s $required ]] || { echo "missing normalized Populus evaluator artifact: $required" >&2; exit 1; }
    done
    chromosomes=$(wc -l < "$output/primary_chromosomes.genome.fa.fai")
    [[ $chromosomes -eq $expected_chromosomes ]] || {
        echo "$bundle chromosome count $chromosomes is not $expected_chromosomes" >&2; exit 1;
    }
    genes=$(grep -Pc '\tgene\t' "$output/primary_chromosomes.gff3")
    transcripts=$(grep -Pc '\t(mRNA|transcript)\t' "$output/primary_chromosomes.gff3")
    cds=$(grep -Pc '\tCDS\t' "$output/primary_chromosomes.gff3")
    [[ $genes -eq $expected_genes && $transcripts -eq $expected_transcripts \
        && $cds -eq $expected_cds ]] || {
        echo "$bundle annotation counts differ from frozen format-only diagnostic" >&2; exit 1;
    }
    printf '%s\t%s\t%s\t%s\t%s\n' "$bundle" "$chromosomes" "$genes" "$transcripts" "$cds" \
        >> "$working_root/annotation_counts.tsv"
done < "$commands"
printf 'field\tvalue\nformal_status\tvalid_pretruth_artifact\nstage\tnormalization\n' \
    > "$working_root/stage_status.tsv"
(
    cd "$working_root"
    find shared_target candidate_only -type f -print0 | sort -z \
        | xargs -0 sha256sum > BLIND_SHA256SUMS
    find shared_target evaluator_only normalized/target_populus \
        normalized/evaluator_manihot normalized/evaluator_ricinus \
        provider_compatible/target_populus provider_compatible/evaluator_manihot \
        provider_compatible/evaluator_ricinus -type f -print0 | sort -z \
        | xargs -0 sha256sum > EVALUATOR_SHA256SUMS
    sha256sum -c BLIND_SHA256SUMS >/dev/null
    sha256sum -c EVALUATOR_SHA256SUMS >/dev/null
)
du -sb "$working_root" > "$working_root/disk_bytes.txt"
(
    cd "$working_root"
    find . -type f ! -path './SHA256SUMS' -print0 \
        | sort -z | xargs -0 sha256sum > SHA256SUMS
    sha256sum -c SHA256SUMS >/dev/null
)
trap - ERR
mv "$working_root" "$result_root"
printf 'normalized Populus evaluator inputs frozen: %s\n' "$result_root"
