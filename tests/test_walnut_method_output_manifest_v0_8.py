from __future__ import annotations

import csv
import importlib.util
from pathlib import Path

import pytest

from ploidypatch.artifact_manifest import sha256_file


ROOT = Path(__file__).parents[1]
SPEC = importlib.util.spec_from_file_location(
    "verify_published_method_output_v0_8",
    ROOT / "scripts/verify_published_method_output_v0.8.py",
)
assert SPEC is not None and SPEC.loader is not None
VERIFY = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(VERIFY)
POOL_SPEC = importlib.util.spec_from_file_location(
    "build_walnut_h1_candidate_pools_v0_8",
    ROOT / "scripts/build_walnut_h1_candidate_pools_v0.8.py",
)
assert POOL_SPEC is not None and POOL_SPEC.loader is not None
POOL = importlib.util.module_from_spec(POOL_SPEC)
POOL_SPEC.loader.exec_module(POOL)


def make_output(tmp_path: Path) -> Path:
    working = tmp_path / "method.working"
    upstream = working / "upstream"
    upstream.mkdir(parents=True)
    output = upstream / "result.gff3"
    output.write_text("##gff-version 3\n", encoding="utf-8")
    with (working / "output_manifest.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("role", "bytes", "sha256", "path"))
        writer.writerow(
            ("result.gff3", output.stat().st_size, sha256_file(output), str(output))
        )
    final = tmp_path / "method"
    working.rename(final)
    return final


def test_published_manifest_survives_atomic_working_rename(tmp_path: Path) -> None:
    output = make_output(tmp_path)
    VERIFY.verify(output, {"result.gff3"})


def test_published_manifest_survives_exact_blind_namespace_remount(
    tmp_path: Path,
) -> None:
    project = tmp_path / "blind-run" / "project"
    project.mkdir(parents=True)
    output = make_output(project)
    manifest = output / "output_manifest.tsv"
    rows = list(csv.DictReader(manifest.open(newline=""), delimiter="\t"))
    rows[0]["path"] = "/run/blind-run/project/method.working/upstream/result.gff3"
    with manifest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("role", "bytes", "sha256", "path"),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)
    VERIFY.verify(output, {"result.gff3"})


@pytest.mark.parametrize(
    "recorded",
    (
        "/foreign/project/method.working/upstream/result.gff3",
        "/run/blind-run/project/other/method.working/upstream/result.gff3",
        "/run/blind-run/project/method/upstream/result.gff3",
        "run/blind-run/project/method.working/upstream/result.gff3",
    ),
)
def test_published_manifest_rejects_noncanonical_namespace_lineage(
    tmp_path: Path, recorded: str
) -> None:
    project = tmp_path / "blind-run" / "project"
    project.mkdir(parents=True)
    output = make_output(project)
    manifest = output / "output_manifest.tsv"
    rows = list(csv.DictReader(manifest.open(newline=""), delimiter="\t"))
    rows[0]["path"] = recorded
    with manifest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("role", "bytes", "sha256", "path"),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)
    with pytest.raises(ValueError, match="working-path lineage differs"):
        VERIFY.verify(output, {"result.gff3"})


def test_published_manifest_rejects_role_or_byte_drift(tmp_path: Path) -> None:
    output = make_output(tmp_path)
    with pytest.raises(ValueError, match="roles differ"):
        VERIFY.verify(output, {"different.gff3"})
    (output / "upstream/result.gff3").write_text("changed\n", encoding="utf-8")
    with pytest.raises(ValueError, match="output differs"):
        VERIFY.verify(output, {"result.gff3"})


def test_pool_builder_retains_a_failed_atomic_tree(tmp_path: Path) -> None:
    output = tmp_path / "pool"
    with pytest.raises((FileNotFoundError, ValueError)):
        POOL.build_pools(
            project_root=tmp_path / "project",
            base_gff=tmp_path / "perturbed.gff3",
            output=output,
            include_raw_manifest=False,
        )
    assert not output.exists()
    assert output.with_name("pool.invalid_run").is_dir()


def test_candidate_suffix_merge_restores_the_exact_base_once(tmp_path: Path) -> None:
    base = tmp_path / "base.gff3"
    base.write_text(
        "##gff-version 3\nchr1\tbase\tgene\t1\t10\t.\t+\t.\tID=base1\n",
        encoding="utf-8",
        newline="",
    )
    candidate_line = "chr1\ttool\tgene\t20\t30\t.\t+\t.\tID=candidate1\n"
    adapted = tmp_path / "adapted.gff3"
    adapted.write_bytes(base.read_bytes() + b"###\n" + candidate_line.encode())
    candidate_only = tmp_path / "candidate_only.gff3"
    POOL.extract_candidate_suffix(
        base_gff=base, adapted_gff=adapted, output=candidate_only
    )
    assert candidate_only.read_text(encoding="utf-8") == (
        "##gff-version 3\n" + candidate_line
    )
    output = tmp_path / "method.gff3"
    POOL.assemble_method_candidate_gff(
        base_gff=base, merged_candidate_gff=candidate_only, output=output
    )
    assert output.read_bytes() == base.read_bytes() + b"###\n" + candidate_line.encode()


def test_candidate_suffix_rejects_a_nonidentical_base_prefix(tmp_path: Path) -> None:
    base = tmp_path / "base.gff3"
    adapted = tmp_path / "adapted.gff3"
    base.write_text("##gff-version 3\n", encoding="utf-8")
    adapted.write_text("##gff-version 2\n", encoding="utf-8")
    with pytest.raises(ValueError, match="does not preserve base prefix"):
        POOL.extract_candidate_suffix(
            base_gff=base,
            adapted_gff=adapted,
            output=tmp_path / "candidate_only.gff3",
        )


def test_candidate_suffix_requires_the_adapter_boundary(tmp_path: Path) -> None:
    base = tmp_path / "base.gff3"
    adapted = tmp_path / "adapted.gff3"
    base.write_text("##gff-version 3\n", encoding="utf-8")
    adapted.write_bytes(base.read_bytes() + b"not-a-boundary\n")
    with pytest.raises(ValueError, match="lacks boundary marker"):
        POOL.extract_candidate_suffix(
            base_gff=base,
            adapted_gff=adapted,
            output=tmp_path / "candidate_only.gff3",
        )


def test_within_method_exact_chain_collapse_is_deterministic(tmp_path: Path) -> None:
    source = tmp_path / "merged.gff3"
    source.write_text(
        "##gff-version 3\n"
        "chr1\tx\tgene\t1\t9\t.\t+\t.\tID=b_g\n"
        "chr1\tx\tmRNA\t1\t9\t.\t+\t.\tID=b_t;Parent=b_g\n"
        "chr1\tx\tCDS\t1\t9\t.\t+\t0\tParent=b_t\n"
        "chr1\tx\tgene\t1\t9\t.\t+\t.\tID=a_g\n"
        "chr1\tx\tmRNA\t1\t9\t.\t+\t.\tID=a_t;Parent=a_g\n"
        "chr1\tx\tCDS\t1\t9\t.\t+\t0\tParent=a_t\n"
        "chr1\tx\tgene\t20\t28\t.\t+\t.\tID=c_g\n"
        "chr1\tx\tmRNA\t20\t28\t.\t+\t.\tID=c_t;Parent=c_g\n"
        "chr1\tx\tCDS\t20\t28\t.\t+\t0\tParent=c_t\n",
        encoding="utf-8",
        newline="",
    )
    output = tmp_path / "deduplicated.gff3"
    counts = POOL.collapse_duplicate_exact_chains(candidate_gff=source, output=output)
    observed = output.read_text(encoding="utf-8")
    assert "ID=a_t" in observed and "ID=b_t" not in observed
    assert "ID=a_g" in observed and "ID=b_g" not in observed
    assert "ID=c_t" in observed
    assert counts == {
        "input_transcripts": 3,
        "retained_transcripts": 2,
        "collapsed_duplicate_transcripts": 1,
        "duplicate_chain_groups": 1,
    }
