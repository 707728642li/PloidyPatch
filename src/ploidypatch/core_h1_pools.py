"""Species-neutral, no-ranker primitives for core-H1 candidate pools."""
from __future__ import annotations

from collections import defaultdict
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
from typing import Any, Iterable, Mapping

from .artifact_manifest import sha256_file
from .gff import parse_attributes


def run_command(command: Iterable[str], *, stdout: Path, stderr: Path) -> None:
    stdout.parent.mkdir(parents=True, exist_ok=True)
    with stdout.open("x", encoding="utf-8") as out, stderr.open(
        "x", encoding="utf-8"
    ) as err:
        completed = subprocess.run(list(command), stdout=out, stderr=err, check=False)
    if completed.returncode != 0:
        raise RuntimeError(f"Core-H1 command exited {completed.returncode}")


def load_json_object(path: str | Path) -> dict[str, Any]:
    source = Path(path)
    if not source.is_file() or source.is_symlink() or source.stat().st_size == 0:
        raise ValueError(f"Missing, empty or symlinked JSON: {source}")
    value = json.loads(source.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"Expected JSON object: {source}")
    return value


def extract_candidate_suffix(*, base_gff: Path, adapted_gff: Path, output: Path) -> None:
    """Remove one byte-identical base prefix before reference merging."""

    if output.exists() or output.is_symlink():
        raise FileExistsError(f"Refusing to overwrite candidate suffix: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    with base_gff.open("rb") as base, adapted_gff.open("rb") as adapted, output.open(
        "xb"
    ) as destination:
        last_byte = b""
        for chunk in iter(lambda: base.read(1024 * 1024), b""):
            if adapted.read(len(chunk)) != chunk:
                raise ValueError(
                    f"Adapted candidate does not preserve the base prefix: {adapted_gff}"
                )
            last_byte = chunk[-1:]
        if last_byte not in {b"\n", b"\r"} and adapted.read(1) != b"\n":
            raise ValueError(f"Adapted candidate lacks base separator: {adapted_gff}")
        if adapted.readline() not in {b"###\n", b"###\r\n"}:
            raise ValueError(f"Adapted candidate lacks boundary marker: {adapted_gff}")
        destination.write(b"##gff-version 3\n")
        shutil.copyfileobj(adapted, destination, length=1024 * 1024)


def assemble_method_candidate_gff(
    *, base_gff: Path, merged_candidate_gff: Path, output: Path
) -> None:
    """Write the exact base once followed by merged candidate-only features."""

    if output.exists() or output.is_symlink():
        raise FileExistsError(f"Refusing to overwrite method candidate GFF: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("xb") as destination, base_gff.open("rb") as base:
        last_byte = b""
        for chunk in iter(lambda: base.read(1024 * 1024), b""):
            destination.write(chunk)
            last_byte = chunk[-1:]
        if last_byte not in {b"\n", b"\r"}:
            destination.write(b"\n")
        destination.write(b"###\n")
        with merged_candidate_gff.open("rb") as merged:
            for raw in merged:
                if raw.strip() and not raw.startswith(b"#"):
                    if len(raw.rstrip(b"\r\n").split(b"\t")) != 9:
                        raise ValueError("Malformed merged candidate feature")
                    destination.write(raw)


def collapse_duplicate_exact_chains(
    *, candidate_gff: Path, output: Path
) -> dict[str, int]:
    """Ensure exact phased-CDS chains receive one within-method vote."""

    manifest_path = Path(str(output) + ".manifest.json")
    if any(path.exists() or path.is_symlink() for path in (output, manifest_path)):
        raise FileExistsError(f"Refusing to overwrite exact-chain collapse: {output}")
    transcripts: dict[str, tuple[str, str, str]] = {}
    cds: dict[str, list[tuple[int, int, str]]] = defaultdict(list)
    gene_transcripts: dict[str, set[str]] = defaultdict(set)
    with candidate_gff.open(encoding="utf-8", newline="") as handle:
        for line_number, raw in enumerate(handle, start=1):
            stripped = raw.rstrip("\r\n")
            if not stripped or stripped.startswith("#"):
                continue
            fields = stripped.split("\t")
            if len(fields) != 9 or fields[2] not in {
                "gene", "mRNA", "transcript", "exon", "CDS"
            }:
                raise ValueError(f"Unsupported merged feature at line {line_number}")
            attributes, malformed = parse_attributes(fields[8])
            if malformed:
                raise ValueError(f"Malformed merged attributes at line {line_number}")
            if fields[2] in {"mRNA", "transcript"}:
                transcript_id = attributes.get("ID", "")
                parents = tuple(
                    value for value in attributes.get("Parent", "").split(",") if value
                )
                if not transcript_id or transcript_id in transcripts or len(parents) != 1:
                    raise ValueError(f"Invalid merged transcript at line {line_number}")
                transcripts[transcript_id] = (fields[0], fields[6], parents[0])
                gene_transcripts[parents[0]].add(transcript_id)
            elif fields[2] == "CDS":
                parents = tuple(
                    value for value in attributes.get("Parent", "").split(",") if value
                )
                if len(parents) != 1:
                    raise ValueError(f"Invalid merged CDS parent at line {line_number}")
                try:
                    start, end = int(fields[3]), int(fields[4])
                except ValueError as error:
                    raise ValueError(f"Invalid merged CDS coordinate at line {line_number}") from error
                cds[parents[0]].append((start, end, fields[7]))
    if set(cds) - set(transcripts):
        raise ValueError("Merged CDS references an unknown transcript")
    grouped: dict[tuple[str, str, tuple[tuple[int, int, str], ...]], list[str]] = (
        defaultdict(list)
    )
    for transcript_id, (seqid, strand, _gene_id) in transcripts.items():
        chain = tuple(sorted(cds.get(transcript_id, [])))
        if not chain or any(phase not in {"0", "1", "2"} for _a, _b, phase in chain):
            raise ValueError(f"Merged transcript lacks phased CDS: {transcript_id}")
        grouped[(seqid, strand, chain)].append(transcript_id)
    dropped_transcripts: set[str] = set()
    duplicate_groups = 0
    for values in grouped.values():
        if len(values) > 1:
            duplicate_groups += 1
            dropped_transcripts.update(sorted(values)[1:])
    dropped_genes = {
        gene_id
        for gene_id, values in gene_transcripts.items()
        if values and values <= dropped_transcripts
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    with candidate_gff.open(encoding="utf-8", newline="") as source, output.open(
        "x", encoding="utf-8", newline=""
    ) as destination:
        for raw in source:
            stripped = raw.rstrip("\r\n")
            if not stripped or stripped.startswith("#"):
                destination.write(raw)
                continue
            fields = stripped.split("\t")
            attributes, _malformed = parse_attributes(fields[8])
            identifier = attributes.get("ID", "")
            parents = {
                value for value in attributes.get("Parent", "").split(",") if value
            }
            if identifier in dropped_genes or identifier in dropped_transcripts:
                continue
            if parents & dropped_transcripts:
                continue
            destination.write(raw)
    audit: dict[str, Any] = {
        "schema_version": "ploidypatch.within_method_exact_chain_collapse.v1",
        "input_sha256": sha256_file(candidate_gff),
        "output_sha256": sha256_file(output),
        "input_transcripts": len(transcripts),
        "retained_transcripts": len(transcripts) - len(dropped_transcripts),
        "collapsed_duplicate_transcripts": len(dropped_transcripts),
        "duplicate_chain_groups": duplicate_groups,
        "selection_rule": "lexicographically_smallest_namespaced_transcript_id",
        "truth_access": False,
        "ranker_access": False,
    }
    manifest_path.write_text(
        json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return {
        key: int(audit[key])
        for key in (
            "input_transcripts",
            "retained_transcripts",
            "collapsed_duplicate_transcripts",
            "duplicate_chain_groups",
        )
    }


def canonical_tree_record(root: str | Path) -> dict[str, Any]:
    tree = Path(root)
    if not tree.is_dir() or tree.is_symlink():
        raise ValueError(f"Raw prediction tree is missing or symlinked: {tree}")
    files = sorted(
        path for path in tree.rglob("*") if path.is_file() and not path.is_symlink()
    )
    if not files:
        raise ValueError(f"Raw prediction tree is empty: {tree}")
    digest = hashlib.sha256()
    total = 0
    for path in files:
        relative = path.relative_to(tree).as_posix()
        file_digest = sha256_file(path)
        digest.update(f"{file_digest}  {relative}\n".encode("utf-8"))
        total += path.stat().st_size
    return {"file_count": len(files), "bytes": total, "tree_sha256": digest.hexdigest()}


def build_raw_predictions_manifest(
    *,
    holdout_id: str,
    policy_id: str,
    project_root: str | Path,
    raw_prediction_trees: Mapping[str, str | Path],
    candidate_references: tuple[str, ...],
    input_hashes: Mapping[str, str],
    output: str | Path,
) -> dict[str, Any]:
    required = {
        f"{method}__{bundle}"
        for method in ("miniprot", "gemoma", "lifton")
        for bundle in candidate_references
    }
    if set(raw_prediction_trees) != required:
        raise ValueError("Raw predictions require three methods by all references")
    destination = Path(output)
    if destination.exists() or destination.is_symlink():
        raise FileExistsError(f"Refusing to overwrite raw prediction manifest: {destination}")
    root = Path(project_root).resolve()
    records: dict[str, Any] = {}
    for label, raw_path in sorted(raw_prediction_trees.items()):
        tree = Path(raw_path).resolve()
        try:
            relative = tree.relative_to(root).as_posix()
        except ValueError as error:
            raise ValueError(f"Raw prediction tree escapes project root: {tree}") from error
        records[label] = {"relative_path": relative, **canonical_tree_record(tree)}
    if not input_hashes or any(
        len(value) != 64 or any(character not in "0123456789abcdef" for character in value)
        for value in input_hashes.values()
    ):
        raise ValueError("Raw prediction input hashes are incomplete or malformed")
    payload = {
        "schema_version": "ploidypatch.core_h1_raw_predictions.v1",
        "holdout_id": holdout_id,
        "policy_id": policy_id,
        "truth_access": False,
        "ranker_access": False,
        "method_families": ["miniprot", "gemoma", "lifton"],
        "candidate_references": list(candidate_references),
        "within_method_reference_vote_count": 1,
        "tree_hash_algorithm": "sha256_of_sorted_sha256_two_space_relative_path_newline",
        "input_hashes": dict(sorted(input_hashes.items())),
        "raw_prediction_trees": records,
    }
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return payload


def seal_blind_pool_manifest(
    *,
    holdout_id: str,
    policy_id: str,
    manifest_path: str | Path,
    candidate_gff_path: str | Path,
    decisions_path: str | Path,
    raw_predictions_manifest_path: str | Path,
    policy_arm: str,
    reference_scope: str,
    candidate_reference_count: int,
) -> dict[str, Any]:
    schemas = {
        "retain_distinct_phased_CDS_chains": "ploidypatch.method_candidate_pool.v2",
        "suppress_strongly_overlapping_alternative_chains": "ploidypatch.method_consensus.v1",
    }
    if policy_arm not in schemas:
        raise ValueError(f"Unsupported core-H1 arm: {policy_arm}")
    path = Path(manifest_path)
    manifest = load_json_object(path)
    if manifest.get("schema_version") != schemas[policy_arm]:
        raise ValueError("Core-H1 pool manifest schema differs")
    raw = load_json_object(raw_predictions_manifest_path)
    if raw.get("schema_version") != "ploidypatch.core_h1_raw_predictions.v1":
        raise ValueError("Core-H1 raw-prediction manifest schema differs")
    if raw.get("truth_access") is not False or raw.get("ranker_access") is not False:
        raise ValueError("Raw projections violate the blind firewall")
    inputs = manifest.get("inputs")
    outputs = manifest.get("outputs")
    if not isinstance(inputs, dict) or not isinstance(outputs, dict):
        raise ValueError("Core-H1 pool manifest lacks inputs/outputs objects")
    candidate_output = outputs.get("candidate_gff")
    decision_output = outputs.get("decisions")
    if (
        not isinstance(candidate_output, dict)
        or candidate_output.get("sha256") != sha256_file(candidate_gff_path)
        or not isinstance(decision_output, dict)
        or decision_output.get("sha256") != sha256_file(decisions_path)
    ):
        raise ValueError("Core-H1 pool output hashes differ before sealing")
    inputs["raw_predictions_manifest"] = {
        "file_name": "raw_predictions.manifest.json",
        "sha256": sha256_file(raw_predictions_manifest_path),
    }
    manifest.update(
        {
            "holdout_id": holdout_id,
            "policy_id": policy_id,
            "policy_arm": policy_arm,
            "reference_scope": reference_scope,
            "truth_access": False,
            "ranker_access": False,
            "automatic_approval": False,
            "within_method_reference_vote_count": 1,
            "candidate_references_per_method": candidate_reference_count,
        }
    )
    temporary = path.with_name(path.name + ".working")
    if temporary.exists() or temporary.is_symlink():
        raise FileExistsError(f"Refusing stale pool-manifest working file: {temporary}")
    temporary.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    os.replace(temporary, path)
    return manifest

