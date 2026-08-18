# Changelog

All notable changes to PloidyPatch are documented here. The project follows
[Semantic Versioning](https://semver.org/).

## [1.0.0] - 2026-08-10

First stable public release.

### Added

- Chain-preserving integration of structurally distinct plant gene-model
  candidates.
- Duplication-aware evidence modules with explicit applicability and
  abstention behavior.
- Human-review ledgers and immutable, reversible GFF3 patches.
- Checksum-bound external-validation, truth-custody, and release-audit
  primitives.
- A seven-family command-line interface with a machine-checked inventory.
- A distributable five-minute example with byte-identical patch reversion.
- Python 3.11–3.13 CI, deterministic wheel/sdist builds, a non-root container,
  and Conda/Bioconda packaging assets.
- A checksum-bound, self-contained HTML review workbench with ranked triage,
  exact method intersections, searchable candidate tables, shared-axis
  phased-CDS comparisons, topology explanations and an editable review ledger.
- Full-universe JSON/TSV exports, deterministic report generation, accessible
  light/dark presentation and explicit workflow-attention gates.
- PyPI trusted-publishing automation and a tested Bioconda recipe.

### Safety boundary

- PloidyPatch never converts a ranking into an automatic annotation approval.
- Missing evidence causes abstention or baseline behavior, not a gene-absence
  claim.
- Source annotations are never modified in place.

[1.0.0]: https://github.com/707728642li/PloidyPatch/releases/tag/v1.0
