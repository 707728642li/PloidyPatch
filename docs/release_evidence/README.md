# Release-evidence timeline

The subdirectories in this folder are immutable evidence snapshots, not claims
about the current source tree.

| Evidence | Software identity | Scope |
|---|---|---|
| `wheel_core_v0.2/` | development `0.1.0a0` | Historical first wheel/sdist smoke; predates the public report command |
| `wheel_core_v0.3/` | development `0.1.0a0` | Historical reproducible-build smoke; predates the public report command |
| `container_core_v0.1/` | development `0.1.0a0` | Historical non-root container smoke |

The `report` command and its public HTML example were added after the frozen
historical wheel evidence. Those old manifests are deliberately not rewritten
to imply otherwise. The next tagged release must be built by
`scripts/build_release_candidate_v0.3.py`; its evidence must show eight command
families, the full reviewed example, a `validated_reversible_run` report, and
byte-identical wheel and canonical source-distribution builds before upload.
Final registry state is verified against the live GitHub, PyPI, Bioconda and
Zenodo services rather than represented by a stale repository snapshot.

The engineering and literature decisions for the stable release are recorded
in `docs/RELEASE_READINESS_v1.0.md`. The project concept DOI may be used in
release-independent citation metadata; a version DOI records one immutable
release.

This separation preserves provenance: historical evidence remains true, and
new functionality receives new evidence rather than a retroactive claim.
