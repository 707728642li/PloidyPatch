from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from ploidypatch.catalog import write_missing_gene_candidate_catalog
from ploidypatch.perturb import generate_missing_gene_benchmark
from ploidypatch.sampling import sample_candidate_catalog


GFF = (
    "##gff-version 3\n"
    "chr1\ttest\tgene\t1\t90\t.\t+\t.\tID=gene:G1\n"
    "chr1\ttest\tmRNA\t1\t90\t.\t+\t.\tID=transcript:T1;Parent=gene:G1\n"
    "chr1\ttest\texon\t1\t30\t.\t+\t.\tParent=transcript:T1\n"
    "chr1\ttest\tCDS\t1\t30\t.\t+\t0\tParent=transcript:T1\n"
    "chr1\ttest\tgene\t101\t190\t.\t-\t.\tID=gene:G2\n"
    "chr1\ttest\tmRNA\t101\t190\t.\t-\t.\tID=transcript:T2;Parent=gene:G2\n"
    "chr1\ttest\texon\t161\t190\t.\t-\t.\tParent=transcript:T2\n"
    "chr1\ttest\tCDS\t161\t190\t.\t-\t0\tParent=transcript:T2\n"
)


def _write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8", newline="")
    return path


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_stratified_sample_is_deterministic_and_drives_perturbation(
    tmp_path: Path,
) -> None:
    source = _write(tmp_path / "source.gff3", GFF)
    catalog = tmp_path / "catalog.tsv"
    write_missing_gene_candidate_catalog(source, catalog)
    plan = _write(
        tmp_path / "plan.tsv", "transcript_count_bin\tsample_count\n1\t1\n"
    )
    first = tmp_path / "first.tsv"
    second = tmp_path / "second.tsv"

    manifest1 = sample_candidate_catalog(
        catalog_path=catalog, plan_path=plan, output_tsv_path=first, seed=17
    )
    manifest2 = sample_candidate_catalog(
        catalog_path=catalog, plan_path=plan, output_tsv_path=second, seed=17
    )

    assert first.read_bytes() == second.read_bytes()
    assert manifest1["selection"] == manifest2["selection"]
    assert manifest1["output"]["rows"] == 1
    selected_gene = _rows(first)[0]["gene_id"]

    run = tmp_path / "run"
    truth_dir = tmp_path / "evaluator_truth"
    perturb_manifest = generate_missing_gene_benchmark(
        source,
        run,
        count=None,
        seed=17,
        selection_tsv_path=first,
        truth_output_dir=truth_dir,
    )
    truth = json.loads(
        (truth_dir / "hidden_truth.json").read_text(encoding="utf-8")
    )
    assert truth["events"][0]["target"]["gene_id"] == selected_gene
    assert not (run / "hidden_truth.json").exists()
    assert (
        perturb_manifest["outputs"]["hidden_truth"]["storage"]
        == "separate_evaluator_directory"
    )
    assert (
        perturb_manifest["perturbation"]["selection_algorithm"]
        == "evaluator_selection_tsv_v1"
    )


def test_stratified_sample_fails_closed_on_shortage_and_modified_selection(
    tmp_path: Path,
) -> None:
    source = _write(tmp_path / "source.gff3", GFF)
    catalog = tmp_path / "catalog.tsv"
    write_missing_gene_candidate_catalog(source, catalog)
    shortage_plan = _write(
        tmp_path / "shortage.tsv",
        "transcript_count_bin\tsample_count\n1\t3\n",
    )
    with pytest.raises(ValueError, match="has 2 candidates but requests 3"):
        sample_candidate_catalog(
            catalog_path=catalog,
            plan_path=shortage_plan,
            output_tsv_path=tmp_path / "shortage_sample.tsv",
            seed=1,
        )

    plan = _write(
        tmp_path / "plan.tsv", "transcript_count_bin\tsample_count\n1\t1\n"
    )
    selection = tmp_path / "selection.tsv"
    sample_candidate_catalog(
        catalog_path=catalog, plan_path=plan, output_tsv_path=selection, seed=1
    )
    with pytest.raises(ValueError, match="different seed"):
        generate_missing_gene_benchmark(
            source,
            tmp_path / "wrong_seed_run",
            count=None,
            seed=2,
            selection_tsv_path=selection,
        )
    with selection.open("a", encoding="utf-8") as handle:
        handle.write("tampered\n")
    with pytest.raises(ValueError, match="checksum"):
        generate_missing_gene_benchmark(
            source,
            tmp_path / "run",
            count=None,
            seed=1,
            selection_tsv_path=selection,
        )


def test_stratified_sample_can_exclude_a_frozen_prior_selection(
    tmp_path: Path,
) -> None:
    source = _write(tmp_path / "source.gff3", GFF)
    catalog = tmp_path / "catalog.tsv"
    write_missing_gene_candidate_catalog(source, catalog)
    plan = _write(
        tmp_path / "plan.tsv", "transcript_count_bin\tsample_count\n1\t1\n"
    )
    prior = tmp_path / "prior.tsv"
    sample_candidate_catalog(
        catalog_path=catalog,
        plan_path=plan,
        output_tsv_path=prior,
        seed=1,
    )
    held_out = tmp_path / "held_out.tsv"
    manifest = sample_candidate_catalog(
        catalog_path=catalog,
        plan_path=plan,
        output_tsv_path=held_out,
        seed=2,
        exclude_tsv_paths=(prior,),
    )

    assert _rows(prior)[0]["candidate_id"] != _rows(held_out)[0]["candidate_id"]
    assert manifest["exclusions"]["unique_candidates"] == 1
    assert manifest["selection"]["strata"][0]["excluded_candidates"] == 1
