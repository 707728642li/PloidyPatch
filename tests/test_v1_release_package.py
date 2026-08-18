from __future__ import annotations

import hashlib
import importlib.util
import json
import re
import tomllib
from pathlib import Path

from ploidypatch import __version__


ROOT = Path(__file__).resolve().parents[1]
PUBLIC_EMAIL = "litaishan@caf.ac.cn"
PRIVATE_ACCOUNT_EMAIL = "litaishan910706@gmail.com"
REPOSITORY = "https://github.com/707728642li/PloidyPatch"


def test_v1_metadata_is_complete_and_consistent() -> None:
    metadata = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    project = metadata["project"]
    assert project["version"] == __version__ == "1.0.0"
    assert project["authors"] == [{"name": "Taishan Li", "email": PUBLIC_EMAIL}]
    assert project["maintainers"] == [{"name": "Taishan Li", "email": PUBLIC_EMAIL}]
    assert project["license"] == "BSD-3-Clause"
    assert project["license-files"] == ["LICENSE"]
    assert project["urls"]["Repository"] == REPOSITORY
    assert project["scripts"] == {"ploidypatch": "ploidypatch.cli:main"}
    assert not project.get("dependencies")


def test_license_citation_and_public_contact_are_release_ready() -> None:
    license_text = (ROOT / "LICENSE").read_text(encoding="utf-8")
    citation = (ROOT / "CITATION.cff").read_text(encoding="utf-8")
    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    assert "BSD 3-Clause License" in license_text
    assert "Copyright (c) 2026, Taishan Li" in license_text
    assert 'version: "1.0.0"' in citation
    assert 'doi: "10.5281/zenodo.21875561"' in citation
    assert 'license: "BSD-3-Clause"' in citation
    assert f'repository-code: "{REPOSITORY}"' in citation
    assert PUBLIC_EMAIL in citation and PUBLIC_EMAIL in readme
    assert "https://zenodo.org/records/21875561" in readme
    assert PRIVATE_ACCOUNT_EMAIL not in citation
    assert PRIVATE_ACCOUNT_EMAIL not in readme


def test_professional_github_surface_and_user_documentation_exist() -> None:
    required = (
        ".github/assets/ploidypatch-workflow.png",
        ".github/workflows/test.yml",
        ".github/workflows/release.yml",
        ".github/workflows/codeql.yml",
        ".github/dependabot.yml",
        ".github/pull_request_template.md",
        "CHANGELOG.md",
        "CODE_OF_CONDUCT.md",
        "CONTRIBUTING.md",
        "SECURITY.md",
        "SUPPORT.md",
        "docs/INSTALLATION.md",
        "docs/VERSIONING_AND_PAPER_RELEASE_v1.0.md",
        "docs/RELEASE_READINESS_v1.0.0.md",
        "docs/TUTORIAL.md",
        "docs/CLI_PARAMETERS.md",
        "docs/TESTING.md",
        "scripts/run_public_test_suite_v1.py",
        "tests/data/README.md",
    )
    assert all((ROOT / relative).is_file() for relative in required)
    assert not (ROOT / ".github/assets/ploidypatch-banner.svg").exists()
    assert not (ROOT / ".github/assets/workflow.svg").exists()
    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    assert readme.startswith('<h1 align="center">PloidyPatch</h1>')
    workflow_markup = '<img src=".github/assets/ploidypatch-workflow.png"'
    assert workflow_markup in readme
    assert readme.index(workflow_markup) < readme.index("## Why PloidyPatch?")
    for section in (
        "## Installation",
        "## Five-minute verified tutorial",
        "## Command-line interface",
        "## Typical reviewed-patch workflow",
        "## Documentation",
        "## Citation",
        "## License and contact",
    ):
        assert section in readme
    assert "python -m pip install ploidypatch==1.0.0" in readme
    assert "conda install -c conda-forge -c bioconda ploidypatch=1.0.0" in readme
    assert "ploidypatch-1.0.0-py_0.conda" in readme


def test_release_workflow_uses_pypi_trusted_publishing() -> None:
    workflow = (ROOT / ".github/workflows/release.yml").read_text(encoding="utf-8")
    assert "publish-pypi:" in workflow
    assert "id-token: write" in workflow
    assert "environment:" in workflow and "name: pypi" in workflow
    assert "https://pypi.org/project/ploidypatch/" in workflow
    assert "pypa/gh-action-pypi-publish@" in workflow
    assert "packages-dir: dist/" in workflow


def test_checksum_bound_example_has_checkout_stable_line_endings() -> None:
    attributes = (ROOT / ".gitattributes").read_text(encoding="utf-8")
    assert "examples/minimal_reviewed_patch/*.gff3 text eol=lf" in attributes
    assert "examples/minimal_reviewed_patch/*.json text eol=lf" in attributes
    for pattern in ("*.cff", "*.css", "*.html", "*.js", "LICENSE"):
        assert f"{pattern} text eol=lf" in attributes
    assert "config/**/*.json text eol=lf" in attributes
    assert "manuscript/** -text -diff -eol" in attributes
    assert "manuscript/**/*.md text eol=lf" in attributes


def test_generated_parameter_reference_covers_exact_cli_inventory() -> None:
    inventory = json.loads(
        (ROOT / "docs/CLI_COMMAND_INVENTORY_v0.1.json").read_text(encoding="utf-8")
    )
    reference = (ROOT / "docs/CLI_PARAMETERS.md").read_text(encoding="utf-8")
    headings = re.findall(r"^## `ploidypatch (.+)`$", reference, flags=re.MULTILINE)
    assert inventory["package_version"] == "1.0.0"
    assert inventory["leaf_command_count"] == 87
    assert headings == inventory["leaf_commands"]
    assert "| Argument | Required | Type | Default | Choices | Description |" in reference


def test_conda_recipes_cover_local_example_and_container_safe_smoke() -> None:
    recipe = (ROOT / "packaging/conda/meta.yaml").read_text(encoding="utf-8")
    for expected in (
        '{% set version = "1.0.0" %}',
        "noarch: python",
        "ploidypatch = ploidypatch.cli:main",
        "python >=3.11",
        "license: BSD-3-Clause",
        "license_file: LICENSE",
        "examples/minimal_reviewed_patch/run_example.py",
        "- 707728642li",
    ):
        assert expected in recipe
    assert "numpy" not in recipe
    assert "scikit-learn" not in recipe

    bioconda = (ROOT / "packaging/bioconda/meta.yaml").read_text(encoding="utf-8")
    assert '{% set version = "1.0.0" %}' in bioconda
    assert "releases/download/v{{ version }}" in bioconda
    recipe_hashes = re.findall(r"^\s*sha256:\s*([0-9a-f]{64})\s*$", bioconda, re.MULTILINE)
    assert len(recipe_hashes) == 1
    assert recipe_hashes[0] != "9a65b2f752b099f87d1eddc3b75115d79ded819e54ceaf390b22c452590fbaef"
    assert "license_family: BSD" in bioconda
    assert "run_exports:" in bioconda
    assert '{{ pin_subpackage(name, max_pin="x") }}' in bioconda
    assert "ploidypatch patch --help" in bioconda
    assert "pip check" in bioconda
    assert "source_files:" not in bioconda
    assert "examples/minimal_reviewed_patch" not in bioconda

    manifest = (ROOT / "MANIFEST.in").read_text(encoding="utf-8")
    assert "prune packaging/bioconda" in manifest


def test_checksum_bound_tutorial_fixture_is_small_and_exact() -> None:
    root = ROOT / "examples/minimal_reviewed_patch"
    rows = {}
    for line in (root / "input_manifest.tsv").read_text(encoding="utf-8").splitlines()[1:]:
        name, byte_count, digest = line.split("\t")
        rows[name] = (int(byte_count), digest)
    assert set(rows) == {
        "source.gff3",
        "candidate.gff3",
        "decisions.tsv",
        "candidate.gff3.manifest.json",
        "review_decisions.tsv",
        "report_copy_features.tsv",
        "report_scores.tsv",
        "report_topology.tsv",
    }
    assert sum(size for size, _ in rows.values()) < 10_000
    for name, (expected_size, expected_digest) in rows.items():
        payload = (root / name).read_bytes()
        assert len(payload) == expected_size
        assert hashlib.sha256(payload).hexdigest() == expected_digest


def test_public_test_runner_executes_every_shipped_test_without_exclusions() -> None:
    script = ROOT / "scripts/run_public_test_suite_v1.py"
    spec = importlib.util.spec_from_file_location("public_test_suite_v1", script)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    command = module.command(["-x"])
    assert command[1:] == ["-m", "pytest", "-q", "-x"]
    assert "--ignore" not in command
