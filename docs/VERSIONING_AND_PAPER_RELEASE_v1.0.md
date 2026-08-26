# PloidyPatch public versioning

The public software, package registries and manuscript use the single product
identity **PloidyPatch v1.0**. Package registries use the machine-readable
versions `1.0.0` and `1.0.1`. The `1.0.1` patch corrects public authorship and
installation metadata without changing the scientific methods or command
behavior. GitHub uses exact release tags such as `v1.0` and `v1.0.1`.

Identifiers such as `v0.3`, `v0.4`, `v0.5` and `v0.9` in archived protocol,
schema or evidence paths are immutable scientific provenance identifiers. They
are not public product releases. Renaming them would break checksum lineage and
obscure which preregistered estimator or protocol produced an artifact.

Accordingly:

- current PyPI and Bioconda distribution packages use version `1.0.1`;
- the manuscript evaluations and the immutable Zenodo archive remain fixed at
  version `1.0.0`;
- reader-facing text uses **PloidyPatch v1.0**;
- machine-readable protocol, schema and evidence identifiers retain their
  original immutable values;
- future public software changes follow Semantic Versioning and receive a new
  release only after tests, documentation and archival checks pass.

This convention separates the public product identity from reproducibility
metadata without rewriting analysis history.

PyPI and Bioconda currently provide PloidyPatch 1.0.1. This release contains
metadata and distribution updates only and does not change the scientific
methods or command behavior relative to the archived 1.0.0 release.
