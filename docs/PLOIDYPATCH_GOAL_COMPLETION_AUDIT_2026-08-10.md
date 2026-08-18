# PloidyPatch goal-completion audit

> **v1.0 release update:** the public repository has been created, the license
> is BSD-3-Clause, and the stable package version is `1.0.0`. The completion
> analysis below is retained as an immutable pre-release audit. The archive DOI,
> complete manuscript authorship/declarations, author approval, and journal
> submission remain owner- or journal-controlled steps.

Audit date: 2026-08-10
Audited repository commit: `3cf82c0348978a6f6e60398773679da713afcdd3`

## Decision

The scientifically authorized internal work is complete. PloidyPatch has a
tested plant-specific implementation, checksum-bound validation across formal
external systems and natural data, a reproducible release candidate, and a
visually inspected manuscript package. No additional species, topology rescue,
ranker tuning, or endpoint change is scientifically justified before
submission; those actions are explicitly prohibited by the frozen decisions.

The overall publication goal is **externally blocked, not scientifically or
computationally blocked**. The remaining work consists only of author-owned
metadata and legal decisions, creation and archival of the public repository,
final author approval, and journal-portal submission and review. These actions
cannot be inferred or completed truthfully by the analysis pipeline.

## Requirement-by-requirement evidence

| Requirement | Immutable evidence | Audit status | Remaining action |
|---|---|---|---|
| Named, usable tool | PloidyPatch source, CLI and isolated wheel/sdist smoke tests in `docs/release_evidence/wheel_core_v0.3/`; seven command families exercised | **Proved internally** | None for scientific use |
| Plant-specific method | Chain-preserving candidate construction, phased-CDS identity, WGD-aware optional evidence, conflict-winner guard, abstention and byte-identical restoration are implemented and documented | **Proved internally** | None |
| Systematic software validation | Full regression at the audited state: 782 passed, 13 platform-dependent skips, 0 failed; deterministic and fail-closed tests cover manifests, isolation, custody, sentinels and evaluation | **Proved internally** | None |
| Formal external validation | Populus formal positive; Actinidia formal negative because one preregistered topology-coverage safety gate failed despite positive H1/H2 components; Coffea formal-positive H1; Walnut retained as an immutable invalid execution with no biological estimate | **Proved and reported without outcome selection** | None; do not retry or tune on revealed systems |
| Natural biological evidence | Apple long-read case evidence and its source-data lineage in `results/apple_flnc_v0.1_figure/` | **Proved internally** | None |
| Scale and reversibility | Core scaling source data plus zero automatic approvals, conflict-winner preservation and byte-identical restoration sentinels | **Proved internally** | None |
| Publication-grade figures | Six main figures plus nine supplementary figures (36 supplementary panels), real-data source tables, manifests and checksums; rendered manuscript and supplement inspected page by page | **Proved internally** | None unless the journal requests changes |
| Complete scientific manuscript | 10,302-word Software manuscript, cover letter, supplementary figures, source data and reviewer reproducibility archive in `manuscript/submission_v0.8/` | **Proved internally** | Populate owner metadata, render once more and approve |
| Reproducible public software release | Two isolated builds produced byte-identical wheel and sdist; wheel SHA256 `b14592beea9735091e26eb7184c8a38a846616398ac329eb2b20da92eb66a246`, sdist SHA256 `d4be89e48d7768a9449fe14aa0cf0ac0be32856dd5184dff7b26e35d32c0cf01` | **Release candidate proved; public release externally blocked** | Authors choose SPDX license/version, create HTTPS repository and archive release to obtain DOI |
| Portal-ready submission | Submission package is scientifically complete and claim-bounded; owner metadata validator reports no malformed fields | **Externally blocked** | Complete the 16 missing owner fields, obtain all-author approval, materialize final files and upload |
| Publication in *Genome Biology* | A journal-quality submission can be prepared from the frozen package | **External outcome; cannot be proved internally** | Editorial assessment, peer review, revision and acceptance |

## Frozen scientific outcomes

- Populus evaluation: `results/populus_external_v0.4_reveal/evaluation/evaluation.json`,
  SHA256 `0e51f5bde4eb78a11c74179215431254eab938f21684ff2dae2838f619800f37`.
- Actinidia summary: `manuscript/source_data/actinidia_external_v0.5_summary.json`,
  SHA256 `19bc0abc274cf9cc11371d76118e49bf5b21a2b3c4fb28c23e76c067217d2aa1`.
- Coffea evaluation: `manuscript/source_data/coffea_external_v1.0/evaluation.json`,
  SHA256 `4fb69a282f2bdc96cb5eaa27d5c5a29f7543588dcb23f1d5a2fccf44fdaf1764`.
- Walnut invalid-run status: `manuscript/source_data/walnut_external_v0.8_invalid_run7/status.json`,
  SHA256 `0de4945d14865c73422266dc5f66711860acbebe791a62dea21d762562ad76de`.
- Apple natural evidence: `results/apple_flnc_v0.1_figure/figure5_source_data.tsv`,
  SHA256 `bc776988482d5e18a89d4c630d8dadaf7330f4227a19b7f9b87634ae2ee76c08`.
- Scaling source data: `results/core_scaling_v0.1_summary/figure_source_data.tsv`,
  SHA256 `4307905b729afee12cc24fee5bfc8706a603649b3d2adba6f293f56a17908cc7`.

The decisions in `docs/DECISIONS.md` freeze these interpretations. In
particular, v0.9 is retired, Walnut v0.8 must not be retried, and Coffea closes
algorithm development before submission. More post-reveal optimization would
weaken rather than strengthen the evidentiary design.

## Exact owner-controlled blockers

`config/submission_owner_metadata_v0.2.json` remains intentionally fail-closed.
Its 16 unresolved fields are:

1. `authors`
2. `affiliations`
3. `declarations.competing_interests`
4. `declarations.funding`
5. `declarations.author_contributions`
6. `declarations.acknowledgements`
7. `software_release.version`
8. `software_release.license_spdx`
9. `software_release.repository_url`
10. `software_release.archive_doi`
11. `software_release.nonacademic_restrictions`
12. `references.final_export_verified`
13. `submission_date`
14. `cover_letter.all_authors_approved`
15. `cover_letter.original_and_not_under_consideration_elsewhere`
16. `cover_letter.policy_issues_statement`

Until these are supplied, the truthful machine states remain
`formal_public_release=false`, `owner_metadata_complete=false`, and
`portal_upload_ready=false`.

## Minimal completion path

1. Complete `config/submission_owner_metadata_v0.2.json` and set its status to
   `ready`.
2. Create the public repository, select an OSI-approved license, publish the
   frozen release, and archive it to obtain the DOI entered above.
3. Run the fail-closed materializer exactly as documented in
   `manuscript/OWNER_SUBMISSION_HANDOFF_v0.1.md`.
4. Generate the owner-populated DOCX, perform one final native-Word page review,
   obtain all-author approval, and upload the package to *Genome Biology*.

No additional algorithmic experiment is a prerequisite for these four steps.
