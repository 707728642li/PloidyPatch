from __future__ import annotations

from dataclasses import asdict
from io import BytesIO
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import tarfile

import pytest

from ploidypatch.artifact_manifest import sha256_file, verify_sha256sums
from ploidypatch.holdout_contract import (
    CORE_H1_MODEL_VERSION,
    FIXED_CORE_H1_SCIENTIFIC_PARAMETERS,
    SCHEMA_VERSION,
    TRUTH_BLIND_DECLARATIONS,
    load_holdout_contract,
)


ROOT = Path(__file__).parents[1]
SCRIPT = ROOT / "scripts/stage_external_holdout_inputs_v0.5.py"
SPEC = importlib.util.spec_from_file_location("walnut_stage_v08", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
STAGE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = STAGE
SPEC.loader.exec_module(STAGE)


REFERENCE_SPECS = (
    ("target", "Juglans_regia", "target_walnut2", "jre"),
    ("candidate_reference", "Juglans_mandshurica", "candidate_mand", "jma"),
    ("candidate_reference", "Carya_illinoinensis", "candidate_carya", "cil"),
    ("evaluator_reference", "Corylus_avellana", "evaluator_cory", "cav"),
    ("evaluator_reference", "Castanea_mollissima", "evaluator_cast", "cmo"),
)


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def write_single_member_tar(path: Path, member_name: str, content: bytes) -> None:
    with tarfile.open(path, "w:gz", format=tarfile.GNU_FORMAT) as archive:
        member = tarfile.TarInfo(member_name)
        member.size = len(content)
        member.mode = 0o644
        archive.addfile(member, BytesIO(content))


def make_contract_and_sources(
    tmp_path: Path,
) -> tuple[Path, Path, dict[str, object], bytes, Path]:
    source_root = tmp_path / "source"
    references: list[dict[str, object]] = []
    tar_content = b">chr1\nACGTN\n>chr2\nRYSWKMBDHV\n"
    tar_path = source_root / "raw/Juglans_mandshurica/genome.fa.tar.gz"
    for role, species, bundle, prefix in REFERENCE_SPECS:
        artifacts: dict[str, object] = {}
        for artifact, suffix in (
            ("genome", "fa.gz"), ("gff3", "gff3.gz"), ("protein", "faa.gz")
        ):
            relative = Path("raw") / species / f"{artifact}.{suffix}"
            path = source_root / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            if species == "Juglans_mandshurica" and artifact == "genome":
                relative = Path("raw/Juglans_mandshurica/genome.fa.tar.gz")
                path = source_root / relative
                write_single_member_tar(path, "Juglans_mandshurica.genome.fa", tar_content)
                artifacts[artifact] = {
                    "source_relative_path": relative.as_posix(),
                    "bytes": path.stat().st_size,
                    "sha256": sha256_file(path),
                    "container": {
                        "format": "tar.gz",
                        "member_name": "Juglans_mandshurica.genome.fa",
                        "member_bytes": len(tar_content),
                        "member_sha256": sha256_bytes(tar_content),
                    },
                }
            else:
                path.write_bytes(f"{species}\t{artifact}\n".encode())
                artifacts[artifact] = {
                    "source_relative_path": relative.as_posix(),
                    "bytes": path.stat().st_size,
                    "sha256": sha256_file(path),
                }
        references.append(
            {
                "role": role,
                "species_id": species,
                "release": "frozen_release_v1",
                "bundle_id": bundle,
                "wgdi_prefix": prefix,
                "primary_seqid_table": f"config/primary_seqids/{bundle}.tsv",
                "artifacts": artifacts,
            }
        )
    payload: dict[str, object] = {
        "schema_version": SCHEMA_VERSION,
        "holdout_id": "walnut_walnut2_v0.8",
        "policy_id": "ploidypatch_walnut_external_core_h1_v0.8",
        "test_role": "untouched_confirmatory_external_species",
        "model_version": CORE_H1_MODEL_VERSION,
        "references": references,
        "seeds": {"truth_sampler": 20260821, "h1_bootstrap": 20260822},
        "target_resolved_parameters": {
            "primary_chromosome_count": 16,
            "minimum_target_chromosomes_fraction": 0.75,
            "minimum_target_chromosomes": 12,
        },
        "scientific_parameters": asdict(FIXED_CORE_H1_SCIENTIFIC_PARAMETERS),
        "truth_blind": dict(TRUTH_BLIND_DECLARATIONS),
    }
    contract = tmp_path / "contract.json"
    contract.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return contract, source_root, payload, tar_content, tar_path


def manifest_rows(path: Path) -> list[dict[str, str]]:
    import csv
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_generic_stage_extracts_hash_bound_tar_member_into_candidate_role(
    tmp_path: Path,
) -> None:
    contract_path, source_root, _, content, tar_path = make_contract_and_sources(tmp_path)
    contract = load_holdout_contract(contract_path)
    output = tmp_path / "stage"
    STAGE.stage_inputs(
        contract_path=contract_path,
        source_root=source_root,
        output_dir=output,
        code_commit="a" * 40,
        copy_workers=4,
    )
    destination = (
        output
        / "candidate_only/Juglans_mandshurica/Juglans_mandshurica.genome.fa"
    )
    assert destination.read_bytes() == content
    assert destination.stat().st_size != tar_path.stat().st_size
    rows = manifest_rows(output / "role_manifest.tsv")
    row = next(
        row for row in rows
        if row["species_id"] == "Juglans_mandshurica" and row["artifact"] == "genome"
    )
    assert row["source_relative_path"].endswith("genome.fa.tar.gz")
    assert row["staged_relative_path"].endswith("Juglans_mandshurica.genome.fa")
    assert row["bytes"] == str(len(content))
    assert row["sha256"] == row["staged_sha256"] == sha256_bytes(content)
    source_artifact = contract.references[1].genome
    assert source_artifact.bytes == tar_path.stat().st_size
    assert source_artifact.sha256 == sha256_file(tar_path)
    assert source_artifact.staged_bytes == len(content)
    verify_sha256sums(output, ignore_checksum_file=True)


def test_tar_member_hash_failure_cleans_atomic_stage(tmp_path: Path) -> None:
    contract_path, source_root, payload, _, _ = make_contract_and_sources(tmp_path)
    payload["references"][1]["artifacts"]["genome"]["container"][  # type: ignore[index]
        "member_sha256"
    ] = "0" * 64
    contract_path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    output = tmp_path / "stage"
    with pytest.raises(ValueError, match="Member SHA-256"):
        STAGE.stage_inputs(
            contract_path=contract_path,
            source_root=source_root,
            output_dir=output,
            code_commit="b" * 40,
            copy_workers=1,
        )
    assert not output.exists()
    assert not Path(str(output) + ".partial").exists()


def test_tar_outer_hash_failure_happens_before_partial_publication(
    tmp_path: Path,
) -> None:
    contract_path, source_root, payload, _, _ = make_contract_and_sources(tmp_path)
    payload["references"][1]["artifacts"]["genome"]["sha256"] = "0" * 64  # type: ignore[index]
    contract_path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    output = tmp_path / "stage"
    with pytest.raises(ValueError, match="Frozen source SHA-256"):
        STAGE.stage_inputs(
            contract_path=contract_path,
            source_root=source_root,
            output_dir=output,
            code_commit="c" * 40,
        )
    assert not output.exists()
    assert not Path(str(output) + ".partial").exists()
