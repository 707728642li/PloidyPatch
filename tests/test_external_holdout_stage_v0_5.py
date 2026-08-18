from __future__ import annotations

from dataclasses import asdict
import csv
import importlib.util
import json
from pathlib import Path

import pytest

from ploidypatch.artifact_manifest import sha256_file, verify_sha256sums
from ploidypatch.holdout_contract import (
    FIXED_SCIENTIFIC_PARAMETERS,
    SCHEMA_VERSION as CONTRACT_SCHEMA,
    TRUTH_BLIND_DECLARATIONS,
)


ROOT = Path(__file__).parents[1]
MODULE_PATH = ROOT / "scripts" / "stage_external_holdout_inputs_v0.5.py"
SPEC = importlib.util.spec_from_file_location("external_holdout_stage_v05", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


REFERENCE_SPECS = (
    ("target", "Actinidia_chinensis", "target_red5", "ach"),
    ("candidate_reference", "Actinidia_eriantha", "candidate_aeriantha", "aer"),
    ("candidate_reference", "Actinidia_rufa", "candidate_arufa", "aru"),
    ("evaluator_reference", "Rhododendron_simsii", "evaluator_rsimsii", "rsi"),
    ("evaluator_reference", "Diospyros_oleifera", "evaluator_doleifera", "dol"),
)


def make_contract_and_sources(tmp_path: Path) -> tuple[Path, Path, dict[str, object]]:
    source_root = tmp_path / "sources"
    references: list[dict[str, object]] = []
    for role, species, bundle, prefix in REFERENCE_SPECS:
        artifacts: dict[str, object] = {}
        for artifact, suffix in (
            ("genome", "fa.gz"),
            ("gff3", "gff3.gz"),
            ("protein", "faa.gz"),
        ):
            relative = Path("raw") / species / f"{artifact}.{suffix}"
            path = source_root / relative
            path.parent.mkdir(parents=True, exist_ok=True)
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
        "schema_version": CONTRACT_SCHEMA,
        "holdout_id": "actinidia_red5_v0.5",
        "policy_id": "ploidypatch_actinidia_external_validation_v0.5",
        "test_role": "target_level_predeclared_untouched_secondary_replication",
        "model_version": "PloidyPatch_ranker_v0.4",
        "references": references,
        "seeds": {
            "truth_sampler": 20261101,
            "h1_bootstrap": 20261102,
            "h2_bootstrap": 20261103,
            "guard_v03_bootstrap": 20261104,
        },
        "target_resolved_parameters": {
            "primary_chromosome_count": 29,
            "minimum_target_chromosomes_fraction": 0.75,
            "minimum_target_chromosomes": 22,
        },
        "scientific_parameters": asdict(FIXED_SCIENTIFIC_PARAMETERS),
        "truth_blind": dict(TRUTH_BLIND_DECLARATIONS),
    }
    contract = tmp_path / "contract.json"
    contract.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return contract, source_root, payload


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_stage_publishes_disjoint_manifest_driven_role_tree(tmp_path: Path) -> None:
    contract, sources, _ = make_contract_and_sources(tmp_path)
    output = tmp_path / "stage"
    result = MODULE.stage_inputs(
        contract_path=contract,
        source_root=sources,
        output_dir=output,
        code_commit="a" * 40,
        copy_workers=3,
    )
    assert result == output
    assert not Path(str(output) + ".partial").exists()
    rows = read_tsv(output / "role_manifest.tsv")
    assert len(rows) == 15
    assert len({row["staged_relative_path"] for row in rows}) == 15
    assert {
        row["species_id"]
        for row in rows
        if row["staged_relative_path"].startswith("candidate_only/")
    } == {"Actinidia_eriantha", "Actinidia_rufa"}
    assert {
        row["species_id"]
        for row in rows
        if row["staged_relative_path"].startswith(
            "evaluator_only/truth_references/"
        )
    } == {"Rhododendron_simsii", "Diospyros_oleifera"}
    shared = [
        row for row in rows if row["staged_relative_path"].startswith("shared_target/")
    ]
    assert [(row["species_id"], row["artifact"]) for row in shared] == [
        ("Actinidia_chinensis", "genome")
    ]
    target_private = [
        row
        for row in rows
        if row["staged_relative_path"].startswith(
            "evaluator_only/target_complete/Actinidia_chinensis/"
        )
    ]
    assert {row["artifact"] for row in target_private} == {"gff3", "protein"}
    for row in rows:
        staged = output.joinpath(*Path(row["staged_relative_path"]).parts)
        assert staged.is_file() and not staged.is_symlink()
        assert staged.stat().st_size == int(row["bytes"])
        assert sha256_file(staged) == row["sha256"] == row["staged_sha256"]
    assert not any(path.is_symlink() for path in output.rglob("*"))
    verify_sha256sums(output, ignore_checksum_file=True)


def test_stage_contract_freezes_truth_blind_and_target_resolved_claims(
    tmp_path: Path,
) -> None:
    contract, sources, _ = make_contract_and_sources(tmp_path)
    output = tmp_path / "stage"
    MODULE.stage_inputs(
        contract_path=contract,
        source_root=sources,
        output_dir=output,
        code_commit="b" * 40,
    )
    stage = json.loads((output / "role_contract.json").read_text(encoding="utf-8"))
    assert stage["schema_version"] == MODULE.SCHEMA_VERSION
    assert stage["counts"] == {
        "artifacts": 15,
        "bytes": sum(path.stat().st_size for path in sources.rglob("*") if path.is_file()),
        "candidate_references": 2,
        "evaluator_references": 2,
        "references": 5,
        "target_references": 1,
    }
    assert stage["truth_blind"] == TRUTH_BLIND_DECLARATIONS
    assert stage["target_resolved_parameters"] == {
        "primary_chromosome_count": 29,
        "minimum_target_chromosomes_fraction": 0.75,
        "minimum_target_chromosomes": 22,
    }
    assert stage["role_boundaries"]["candidate_evaluator_species_overlap"] is False
    assert stage["contract"]["sha256"] == sha256_file(contract)


def test_stage_rejects_source_mutation_without_partial_publication(
    tmp_path: Path,
) -> None:
    contract, sources, payload = make_contract_and_sources(tmp_path)
    relative = payload["references"][0]["artifacts"]["genome"][  # type: ignore[index]
        "source_relative_path"
    ]
    (sources / relative).write_bytes(b"mutated\n")
    output = tmp_path / "stage"
    with pytest.raises(ValueError, match="byte count|SHA-256"):
        MODULE.stage_inputs(
            contract_path=contract,
            source_root=sources,
            output_dir=output,
            code_commit="c" * 40,
        )
    assert not output.exists()
    assert not Path(str(output) + ".partial").exists()


def test_stage_copy_failure_is_atomic_and_cleans_partial(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    contract, sources, _ = make_contract_and_sources(tmp_path)
    output = tmp_path / "stage"
    real_copy = MODULE.shutil.copyfile
    calls = 0

    def fail_second(*args: object, **kwargs: object) -> object:
        nonlocal calls
        calls += 1
        if calls == 2:
            raise OSError("synthetic copy failure")
        return real_copy(*args, **kwargs)

    monkeypatch.setattr(MODULE.shutil, "copyfile", fail_second)
    with pytest.raises(OSError, match="synthetic"):
        MODULE.stage_inputs(
            contract_path=contract,
            source_root=sources,
            output_dir=output,
            code_commit="d" * 40,
            copy_workers=1,
        )
    assert not output.exists()
    assert not Path(str(output) + ".partial").exists()


@pytest.mark.parametrize("existing", ("output", "partial"))
def test_stage_refuses_overwrite(
    tmp_path: Path, existing: str
) -> None:
    contract, sources, _ = make_contract_and_sources(tmp_path)
    output = tmp_path / "stage"
    target = output if existing == "output" else Path(str(output) + ".partial")
    target.mkdir()
    with pytest.raises(FileExistsError, match="overwrite"):
        MODULE.stage_inputs(
            contract_path=contract,
            source_root=sources,
            output_dir=output,
            code_commit="e" * 40,
        )


def test_stage_rejects_symlinked_source(tmp_path: Path) -> None:
    contract, sources, payload = make_contract_and_sources(tmp_path)
    relative = Path(
        payload["references"][0]["artifacts"]["genome"][  # type: ignore[index]
            "source_relative_path"
        ]
    )
    source = sources / relative
    real = source.with_name("real.fa.gz")
    source.rename(real)
    try:
        source.symlink_to(real.name)
    except OSError:
        pytest.skip("Host does not permit symlink creation")
    with pytest.raises(ValueError, match="symlink"):
        MODULE.stage_inputs(
            contract_path=contract,
            source_root=sources,
            output_dir=tmp_path / "stage",
            code_commit="f" * 40,
        )


def test_stage_requires_full_lowercase_git_sha(tmp_path: Path) -> None:
    contract, sources, _ = make_contract_and_sources(tmp_path)
    with pytest.raises(ValueError, match="full lowercase Git SHA"):
        MODULE.stage_inputs(
            contract_path=contract,
            source_root=sources,
            output_dir=tmp_path / "stage",
            code_commit="not-a-commit",
        )
