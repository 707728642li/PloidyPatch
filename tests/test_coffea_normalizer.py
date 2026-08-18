from __future__ import annotations

import csv
import gzip
import importlib.util
import json
from pathlib import Path
import shutil
import sys

import pytest

from ploidypatch.artifact_manifest import sha256_file, verify_sha256sums, write_sha256sums
from ploidypatch.holdout_contract import load_holdout_contract, staged_relative_path


ROOT = Path(__file__).parents[1]
CONTRACT = ROOT / "config/holdouts/coffea_et39_v1.0/contract.json"


def load_script(relative: str, name: str):
    path = ROOT / relative
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


NORMALIZER = load_script(
    "scripts/prepare_coffea_external_normalized_inputs_v1.0.py",
    "coffea_normalizer_v10",
)


def primary_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_synthetic_stage(stage: Path, execution: Path) -> None:
    contract = load_holdout_contract(CONTRACT)
    source_root = execution / "source"
    aliases = {
        row["genome_seqid"]: row["gff_seqid"]
        for row in primary_rows(ROOT / "config/coffea_gardenia_seqid_aliases_v1.0.tsv")
    }
    rows: list[dict[str, object]] = []
    for reference in contract.references:
        source_primary = ROOT / reference.primary_seqid_table
        frozen_primary = source_root / reference.primary_seqid_table
        frozen_primary.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source_primary, frozen_primary)
        seqids = [row["seqid"] for row in primary_rows(source_primary)]
        for artifact_name, artifact in reference.artifact_items():
            relative = staged_relative_path(reference, artifact_name)
            path = stage.joinpath(*relative.parts)
            path.parent.mkdir(parents=True, exist_ok=True)
            if artifact_name == "genome":
                content = "".join(
                    f">{seqid} provider header\nACGTACGTACGT\n" for seqid in seqids
                )
            elif artifact_name == "gff3":
                gff_seqids = [aliases.get(seqid, seqid) for seqid in seqids]
                content = (
                    "##gff-version 3\n"
                    + "".join(
                        f"{seqid}\ttest\tgene\t1\t12\t.\t+\t.\tID={reference.wgdi_prefix}_{index}\n"
                        for index, seqid in enumerate(gff_seqids, start=1)
                    )
                )
            else:
                content = ">P1\nMAAA\n"
            if path.suffix == ".gz":
                with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
                    handle.write(content)
            else:
                path.write_text(content, encoding="utf-8")
            rows.append(
                {
                    "role": reference.role,
                    "species_id": reference.species_id,
                    "release": reference.release,
                    "bundle_id": reference.bundle_id,
                    "wgdi_prefix": reference.wgdi_prefix,
                    "artifact": artifact_name,
                    "bytes": path.stat().st_size,
                    "sha256": sha256_file(path),
                    "source_relative_path": artifact.source_relative_path.as_posix(),
                    "staged_relative_path": relative.as_posix(),
                    "staged_sha256": sha256_file(path),
                }
            )

    frozen_alias = source_root / "config/coffea_gardenia_seqid_aliases_v1.0.tsv"
    frozen_alias.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(
        ROOT / "config/coffea_gardenia_seqid_aliases_v1.0.tsv", frozen_alias
    )
    with (stage / "role_manifest.tsv").open("x", encoding="utf-8", newline="") as handle:
        fields = (
            "role",
            "species_id",
            "release",
            "bundle_id",
            "wgdi_prefix",
            "artifact",
            "bytes",
            "sha256",
            "source_relative_path",
            "staged_relative_path",
            "staged_sha256",
        )
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    (stage / "role_contract.json").write_text("{}\n", encoding="utf-8")
    write_sha256sums(stage)
    write_sha256sums(execution)


@pytest.mark.skipif(
    sys.platform == "win32",
    reason="Synthetic five-role tree exceeds legacy Windows path-length limits",
)
def test_coffea_normalizer_enforces_role_separation_and_exact_gardenia_aliases(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    stage = tmp_path / "stage"
    execution = tmp_path / "execution"
    stage.mkdir()
    execution.mkdir()
    write_synthetic_stage(stage, execution)
    project = tmp_path / "project"
    project.mkdir()
    monkeypatch.setenv("PLOIDYPATCH_STAGED_INPUT_ROOT", str(stage))
    monkeypatch.setenv("PLOIDYPATCH_EXECUTION_FREEZE", str(execution))
    monkeypatch.setenv("PLOIDYPATCH_HOLDOUT_CONTRACT", str(CONTRACT))

    assert NORMALIZER.main([str(project)]) == 0
    output = project / "data/derived/external_inputs/coffea/v1.0"
    verify_sha256sums(output, ignore_checksum_file=True)
    gardenia = (
        output
        / "evaluator_only/normalized/truth_references/evaluator_gardenia/primary_chromosomes.gff3"
    ).read_text(encoding="utf-8")
    assert "CM023095.1" in gardenia
    assert "Gardenia1" not in gardenia
    blind_names = {path.name for path in (output / "blind").rglob("*")}
    assert "target_complete" not in blind_names
    assert "evaluator_gardenia" not in blind_names
    assert "evaluator_ophiorrhiza" not in blind_names
    assert (
        output / "blind/shared_target/target_et39/primary_chromosomes.genome.fa"
    ).is_file()
    manifest = json.loads((output / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["firewall"] == {
        "blind_complete_target_annotation": False,
        "blind_evaluator_references": False,
        "blind_truth": False,
        "blind_labels": False,
        "blind_homoeolog_group_table": False,
    }

    with pytest.raises(FileExistsError, match="overwrite"):
        NORMALIZER.main([str(project)])
