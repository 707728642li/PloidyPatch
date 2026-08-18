# PloidyPatch public versioning

The public software, package registries and manuscript use the single product
identity **PloidyPatch v1.0**. The machine-readable Semantic Versioning tag is
`v1.0.0`.

Identifiers such as `v0.3`, `v0.4`, `v0.5` and `v0.9` in archived protocol,
schema or evidence paths are immutable scientific provenance identifiers. They
are not public product releases. Renaming them would break checksum lineage and
obscure which preregistered estimator or protocol produced an artifact.

Accordingly:

- the command-line package, container, documentation and registry records use
  version `1.0.0`;
- reader-facing text uses **PloidyPatch v1.0**;
- machine-readable protocol, schema and evidence identifiers retain their
  original immutable values;
- future public software changes follow Semantic Versioning and receive a new
  release only after tests, documentation and archival checks pass.

This convention separates the public product identity from reproducibility
metadata without rewriting analysis history.
