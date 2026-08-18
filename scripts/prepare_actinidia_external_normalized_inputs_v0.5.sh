
#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then echo "usage: $0 PROJECT_ROOT" >&2; exit 2; fi
project_root=$(realpath "$1")
parallel_bin=/data/codexli/software/conda/miniforge3/bin/parallel
input_root=$project_root/data/derived/holdout_inputs/actinidia_v0.5
protocol_root=$project_root/results/protocol_freezes/actinidia_external_v0.5
execution_root=$project_root/results/protocol_freezes/actinidia_external_v0.5_execution
code_root=$execution_root/source
export PYTHONPATH="$code_root/src${PYTHONPATH:+:$PYTHONPATH}"
environment_bindings=$execution_root/environment_bindings.tsv
[[ -s $environment_bindings ]] || { echo "missing frozen environment bindings" >&2; exit 1; }
dev_prefix=$(awk -F '\t' '$1 == "ploidypatch-dev" {print $2}' "$environment_bindings")
[[ $dev_prefix == /* ]] || { echo "invalid frozen ploidypatch-dev binding" >&2; exit 1; }
python_bin=$dev_prefix/bin/python
result_root=$project_root/data/derived/external_inputs/actinidia/v0.5
working_root=${result_root}.working
self_relative=scripts/prepare_actinidia_external_normalized_inputs_v0.5.sh

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
                "$protocol_root/SHA256SUMS" "$protocol_root/role_manifest.tsv" \
                "$execution_root/SHA256SUMS" \
                "$execution_root/implementation_manifest.tsv" \
                "${implementation_dependencies[@]/#/$code_root/}" \
                "$code_root/config/primary_seqids/actinidia_chinensis_red5.tsv" \
                "$code_root/config/primary_seqids/rhododendron_simsii_gca_014282245.1.tsv" \
                "$code_root/config/primary_seqids/diospyros_oleifera_v1.0.tsv"; do
    [[ -s $required ]] || { echo "missing Actinidia normalization prerequisite: $required" >&2; exit 1; }
done
(cd "$protocol_root" && sha256sum -c SHA256SUMS >/dev/null)
(cd "$execution_root" && sha256sum -c SHA256SUMS >/dev/null)
(cd "$input_root" && sha256sum -c SHA256SUMS >/dev/null)
for relative in "${implementation_dependencies[@]}"; do verify_implementation "$relative"; done
[[ ! -e $result_root && ! -e $working_root ]] || {
    echo "refusing to overwrite normalized Actinidia evaluator inputs" >&2; exit 1;
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
truth_blind = contract.get("truth_blind", {})
required_false = (
    "wgd_pairs_enumerated_before_protocol_freeze",
    "candidate_counts_computed_before_protocol_freeze",
    "truth_labels_accessed_before_protocol_freeze",
    "candidate_truth_access",
    "candidate_evaluator_reference_access",
)
if (
    contract.get("schema_version") != "ploidypatch.external_holdout_input_stage.v0.5"
    or contract.get("policy_id") != "ploidypatch_actinidia_external_validation_v0.5"
    or any(truth_blind.get(key) is not False for key in required_false)
):
    raise ValueError("staged Actinidia role contract violates frozen pre-truth scope")
with (input_root / "role_manifest.tsv").open(encoding="utf-8", newline="") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))
with (protocol_root / "role_manifest.tsv").open(
    encoding="utf-8", newline=""
) as handle:
    frozen_rows = list(csv.DictReader(handle, delimiter="\t"))
frozen = {
    (row["role"], row["species_id"], row["artifact"]): row
    for row in frozen_rows
}
if len(frozen) != len(frozen_rows) or len(frozen) != 15:
    raise ValueError("protocol role manifest must contain 15 unique frozen artifacts")
staged = {
    (row["role"], row["species_id"], row["artifact"]): row
    for row in rows
}
if len(staged) != len(rows) or set(staged) != set(frozen):
    raise ValueError("staged and protocol role-manifest artifact universes differ")
for key, row in staged.items():
    frozen_row = frozen[key]
    for field in ("release", "bundle_id", "wgdi_prefix", "bytes", "sha256", "staged_relative_path"):
        if row[field] != frozen_row[field]:
            raise ValueError(f"staged artifact differs from protocol role manifest: {key}/{field}")
expected = {
    ("target", "Actinidia_chinensis_Red5", artifact)
    for artifact in ("genome", "gff3", "protein")
} | {
    ("evaluator_reference", species, artifact)
    for species in ("Rhododendron_simsii", "Diospyros_oleifera")
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
    frozen_row = frozen.get(key)
    if (
        frozen_row is None
        or frozen_row["bytes"] != row["bytes"]
        or frozen_row["sha256"] != row["sha256"]
        or frozen_row["release"] != row["release"]
    ):
        raise ValueError(f"staged evaluator role differs from protocol preflight: {key}")
    selected[key] = path
if set(selected) != expected:
    raise ValueError(f"staged evaluator artifacts incomplete: {expected - set(selected)}")

names = (
    ("target_red5", "target", "Actinidia_chinensis_Red5"),
    ("evaluator_rhododendron", "evaluator_reference", "Rhododendron_simsii"),
    ("evaluator_diospyros", "evaluator_reference", "Diospyros_oleifera"),
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
        target_red5)
            primary=$code_root/config/primary_seqids/actinidia_chinensis_red5.tsv
            expected_chromosomes=29; expected_genes=32950
            expected_transcripts=33021; expected_cds=180837 ;;
        evaluator_rhododendron)
            primary=$code_root/config/primary_seqids/rhododendron_simsii_gca_014282245.1.tsv
            expected_chromosomes=13; expected_genes=30804
            expected_transcripts=30419; expected_cds=166134 ;;
        evaluator_diospyros)
            primary=$code_root/config/primary_seqids/diospyros_oleifera_v1.0.tsv
            expected_chromosomes=15; expected_genes=29268
            expected_transcripts=29268; expected_cds=136033 ;;
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
        [[ -s $required ]] || { echo "missing normalized Actinidia evaluator artifact: $required" >&2; exit 1; }
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

# Red5's provider peptide file contains 94 non-chromosome proteins while the
# frozen target annotation is chromosome-only.  Build the only permitted
# target peptide input from the exact CDS protein_id whitelist; never truncate
# or fuzzy-match identifiers.
target_provider_protein=$(awk -F '\t' '$1 == "target_red5" {print $6}' "$resolution")
[[ -s $target_provider_protein ]] || { echo "missing unique Red5 provider protein" >&2; exit 1; }
"$python_bin" - \
    "$working_root/normalized/target_red5/primary_chromosomes.gff3" \
    "$target_provider_protein" \
    "$working_root/normalized/target_red5/chromosome_protein_id_whitelist.tsv" \
    "$working_root/normalized/target_red5/chromosome_whitelist.pep.fa" \
    "$working_root/normalized/target_red5/chromosome_whitelist.manifest.json" <<'PY'
import gzip
import hashlib
import json
import sys
from pathlib import Path

gff, protein, whitelist_path, output, manifest_path = map(Path, sys.argv[1:])
protein_ids = set()
with gff.open(encoding="utf-8") as handle:
    for line in handle:
        if line.startswith("#") or not line.strip():
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) != 9 or fields[2] != "CDS":
            continue
        attrs = {}
        for token in fields[8].split(";"):
            if "=" in token:
                key, value = token.split("=", 1)
                attrs[key] = value
        protein_id = attrs.get("protein_id")
        if not protein_id:
            raise ValueError("Red5 chromosome CDS lacks protein_id")
        protein_ids.add(protein_id)
if len(protein_ids) != 33021:
    raise ValueError(f"Red5 chromosome protein_id whitelist is not 33021: {len(protein_ids)}")
with whitelist_path.open("x", encoding="utf-8") as handle:
    handle.write("protein_id\n")
    for protein_id in sorted(protein_ids):
        handle.write(f"{protein_id}\n")

opener = gzip.open if protein.open("rb").read(2) == b"\x1f\x8b" else open
seen = set()
written = 0
keep = False
with opener(protein, "rt", encoding="utf-8") as incoming, output.open("x", encoding="utf-8") as outgoing:
    for line in incoming:
        if line.startswith(">"):
            identifier = line[1:].split(None, 1)[0]
            keep = identifier in protein_ids
            if keep:
                if identifier in seen:
                    raise ValueError(f"duplicate Red5 peptide identifier: {identifier}")
                seen.add(identifier)
                written += 1
        if keep:
            outgoing.write(line)
if seen != protein_ids or written != 33021:
    raise ValueError("Red5 peptide whitelist does not map exactly and uniquely")
digest = hashlib.sha256(output.read_bytes()).hexdigest()
with manifest_path.open("x", encoding="utf-8") as handle:
    json.dump(
        {
            "schema_version": "ploidypatch.red5_chromosome_peptide_whitelist.v0.5",
            "mapping": "exact_GFF_CDS_protein_id_to_FASTA_first_token",
            "all_provider_proteins": 33115,
            "accepted_chromosome_proteins": written,
            "excluded_nonchromosome_proteins": 94,
            "output_sha256": digest,
        },
        handle,
        indent=2,
        sort_keys=True,
    )
    handle.write("\n")
PY
printf 'field\tvalue\nformal_status\tvalid_pretruth_artifact\nstage\tnormalization\n' \
    > "$working_root/stage_status.tsv"
(
    cd "$working_root"
    find shared_target candidate_only -type f -print0 | sort -z \
        | xargs -0 sha256sum > BLIND_SHA256SUMS
    find shared_target evaluator_only normalized/target_red5 \
        normalized/evaluator_rhododendron normalized/evaluator_diospyros \
        provider_compatible/target_red5 provider_compatible/evaluator_rhododendron \
        provider_compatible/evaluator_diospyros -type f -print0 | sort -z \
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
printf 'normalized Actinidia evaluator inputs frozen: %s\n' "$result_root"
