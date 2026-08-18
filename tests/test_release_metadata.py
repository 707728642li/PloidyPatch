from pathlib import Path


ROOT = Path(__file__).parents[1]


def test_test_extra_covers_numeric_evaluators() -> None:
    text = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
    assert "[project.optional-dependencies]" in text
    assert '"numpy>=1.26"' in text
    assert '"scipy>=1.12"' in text
    assert '"scikit-learn>=1.4"' in text
    assert '"pytest>=8"' in text


def test_ci_runs_all_supported_python_minors() -> None:
    text = (ROOT / ".github" / "workflows" / "test.yml").read_text(
        encoding="utf-8"
    )
    assert 'python-version: ["3.11", "3.12", "3.13"]' in text
    assert 'python -m pip install -e ".[test]"' in text
    assert "python scripts/run_public_test_suite_v1.py" in text
    assert "fail-fast: false" in text


def test_source_distribution_includes_research_reproduction_assets() -> None:
    text = (ROOT / "MANIFEST.in").read_text(encoding="utf-8")
    assert "recursive-include config *" in text
    assert "recursive-include docs *.md" in text
    assert "recursive-include examples *" in text
    assert "recursive-include scripts *.py *.sh" in text
    assert "include LICENSE" in text
    assert "include CITATION.cff" in text
    assert "recursive-include packaging *.yaml *.md" in text
    workflow = (ROOT / ".github" / "workflows" / "test.yml").read_text(
        encoding="utf-8"
    )
    assert "Verify source distribution research assets" in workflow
    assert "examples/minimal_reviewed_patch/run_example.py" in workflow
    assert "config/holdouts/coffea_et39_v1.0/contract.json" in workflow
    assert "/docs/CLI_PARAMETERS.md$" in workflow
    assert "/docs/examples/report_preview/index.html$" in workflow
    assert "/examples/minimal_reviewed_patch/report_scores.tsv$" in workflow
    assert "/packaging/conda/meta.yaml$" in workflow


def test_ci_requires_reproducible_wheel_and_sdist_smoke() -> None:
    workflow = (ROOT / ".github" / "workflows" / "test.yml").read_text(
        encoding="utf-8"
    )
    assert "scripts/build_release_candidate_v0.3.py" in workflow
    assert '--work-root "$RUNNER_TEMP/ploidypatch-release"' in workflow
    assert 'report["wheel_byte_identical"] is True' in workflow
    assert 'report["sdist_byte_identical"] is True' in workflow
    assert 'manifest["installed_smoke"]["wheel"]["pip_check"] == "passed"' in workflow
    assert 'manifest["installed_smoke"]["sdist"]["pip_check"] == "passed"' in workflow
    assert 'example_report_state"] == "validated_reversible_run"' in workflow
    assert "ploidypatch-release-evidence/*" in workflow


def test_readme_installation_preserves_environment_and_release_boundaries() -> None:
    text = (ROOT / "README.md").read_text(encoding="utf-8")
    assert "## Installation" in text
    assert "conda create -p envs/ploidypatch" in text
    assert 'python -m pip install -e ".[test]"' in text
    assert "does not silently install external gene predictors" in text
    assert "Bioconda-ready recipe" in text
    assert "After the community recipe is merged" in text
    assert "packaging/conda/meta.yaml" in (
        ROOT / "packaging" / "bioconda" / "README.md"
    ).read_text(encoding="utf-8")
