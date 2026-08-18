# PloidyPatch v1.0.0 release readiness

Date: 2026-08-11
Scope: the single public software identity used by the paper and package registries.

## Release decision

PloidyPatch v1.0.0 is the stable software and paper release. It includes the
review workbench, packaging, documentation and engineering improvements without
retraining a ranker, changing a scientific threshold,
reinterpret a frozen external result, or make candidate approval automatic.
The release keeps the evidence hierarchy explicit: candidate discovery, ranking,
human review, immutable patch application, and byte-exact reversion remain
separate stages.

## Frontier review translated into design decisions

| Recent work or standard | Relevant frontier | PloidyPatch decision |
|---|---|---|
| [TOGA2 (bioRxiv, 2026)](https://www.biorxiv.org/content/10.64898/2026.06.30.735536v1.full) | Exon-level orthology, gene-tree reconciliation, splice prediction, multiple references and explicit sensitivity to reference quality | Keep orthology/topology as provenance-rich optional evidence; never treat it as universally applicable or sufficient for automatic repair |
| [LiftOn (Genome Research, 2025)](https://genome.cshlp.org/content/early/2025/01/31/gr.279620.124) | Joint DNA/protein evidence improves annotation lift-over | Preserve independent method provenance and report method intersections rather than collapsing support to an opaque consensus |
| [Cross-species annotation benchmark (Genome Research, 2025)](https://genome.cshlp.org/content/early/2025/04/10/gr280377124) | Best method depends on species and available evidence | Report per-method arms, applicability and failure states; avoid a universal-winner claim |
| [Miniprot](https://academic.oup.com/bioinformatics/article/39/1/btad014/6989621) | Fast protein-to-genome alignment is strong upstream evidence but is not a complete annotation or conflict-resolution system | Use Miniprot as an upstream source; perform explicit overlap reconciliation and retain an audit ledger |
| [AEGIS (Bioinformatics, 2026)](https://academic.oup.com/bioinformatics/article/42/6/btag363/8704544) | Emerging annotation systems integrate multiple lift-over engines and overlap analyses | Make overlap/intersection structure first-class in the report and expose the originating models for every consensus candidate |
| [MultiQC](https://github.com/MultiQC/MultiQC) | A single portable interactive report plus machine-readable data is a successful bioinformatics usability pattern | Ship a self-contained HTML review workbench together with JSON/TSV ledgers and checksums |
| [JBrowse 2](https://jbrowse.org/jb2/docs/) | Linked linear and synteny views make genome evidence interpretable | Provide local gene-chain, conflict-set and coordinate context; keep an explicit future integration boundary for full genome browsing |
| [WCAG 2.2](https://www.w3.org/TR/WCAG22/) | Keyboard focus, non-text contrast and reduced-motion behavior are baseline interface quality | Test focus visibility, light/dark contrast, semantic controls and reduced-motion CSS |
| [GitHub secure use guidance](https://docs.github.com/en/actions/reference/security/secure-use) | Full commit SHAs are the immutable way to pin Actions | Pin every third-party workflow action to a 40-character commit SHA and scan Python plus report JavaScript |

## Engineering gates

- source version, CLI version, citation metadata, container label, Conda recipes,
  examples and release documentation agree on `1.0.0`;
- all third-party GitHub Actions are pinned to immutable full commit SHAs;
- Linux test matrix covers Python 3.11--3.13; a macOS smoke job exercises the
  offline reversible example; CodeQL scans Python and JavaScript/TypeScript;
- report JavaScript passes `node --check`; focused Python code passes Ruff and
  Mypy; report tests enforce at least 80% module coverage;
- report assets are packaged in both wheel and source distribution;
- the release builder must produce byte-identical artifacts in two independent
  builds and must smoke-test installation from both wheel and sdist;
- `twine check` must pass for both distributions;
- the reviewed example must finish with `validated_reversible_run` and restore
  the input GFF3 byte-for-byte;
- the Bioconda recipe SHA must equal the final canonical v1.0.0 source archive;
- the PyPI and Bioconda packages must be built from the same immutable v1.0.0
  source archive.

## Reporting architecture

The report is intentionally more than a dashboard. It combines:

1. run status and integrity sentinels;
2. discovery-to-review funnel and exact multi-method intersections;
3. ranked review queues with explicit decisions;
4. linked conflict-set, transcript-chain and phase-aware gene-structure views;
5. reproducibility metadata and downloadable machine-readable ledgers;
6. light/dark, print, keyboard, mobile and reduced-motion presentation modes.

Only a bounded candidate subset is embedded for interactivity; full evidence
remains in TSV/JSON outputs. This keeps large runs usable without hiding the
complete candidate universe.

## External-state boundary

The exact v1.0.0 archive DOI is `10.5281/zenodo.21875561`; the project-level
concept DOI is `10.5281/zenodo.21875560`. PyPI or Bioconda availability is
claimed only after the corresponding registry independently exposes the
package. Internal protocol and artifact-schema identifiers remain immutable
and are not software release numbers.

## Final verification record

The 2026-08-11 clean-snapshot verification produced the following results:

- public release suite: **800 passed, 14 skipped, 0 failed** in 90.74 s;
- complete repository suite: 810 passed and 14 skipped; seven failures were
  confined to the predeclared internal archival boundary (historical manuscript
  checksum trees and upstream result paths not shipped in the public source
  package), with no v1.0.0 tool, CLI, report or patch failure;
- report suite: 17 passed and one environment-specific skip, with measured
  `ploidypatch.review_report` coverage of **87.25%** against an 80% gate;
- report JavaScript: Node.js 22 syntax check passed; light/dark contrast,
  keyboard focus, reduced motion and packaged-asset tests passed;
- Python quality: whole-repository fatal Ruff rules passed; the stricter report
  Ruff subset passed; Mypy reported no issue in the report module;
- workflow audit: every workflow parsed as YAML and every third-party action
  reference was a full 40-character commit SHA;
- two independent PEP 517 builds produced byte-identical wheel and canonical
  sdist archives, stable example outputs and normalized command logs;
- installations from both wheel and sdist passed `pip check`, exposed all eight
  command families, produced `validated_reversible_run`, made no automatic
  approval, and reverted the example GFF3 byte-for-byte;
- `twine check` passed for both distribution formats;
- the canonical sdist excluded the external Bioconda recipe while retaining the
  packaged HTML/CSS/JavaScript report assets, removing checksum self-reference.

Artifact hashes are kept in the external release-candidate evidence generated
by `scripts/build_release_candidate_v0.3.py`, because embedding an archive's own
hash inside that archive would make the claim self-referential. The Bioconda
recipe is outside the sdist and therefore can bind the final canonical sdist
hash exactly.
