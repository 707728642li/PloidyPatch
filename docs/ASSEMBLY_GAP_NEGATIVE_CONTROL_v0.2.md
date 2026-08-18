# Assembly-gap negative control v0.2 — 2026-08-07

## Claim boundary

This is a paired synthetic assembly-gap control for the frozen 800-event
Da-Ae annotation-missing benchmark. It asks whether a repair method abstains
when the genomic sequence at a hidden gene locus has been replaced by `N`
without changing coordinates.

It is not evidence about real presence/absence variation (PAV), structural
deletion, or population-level loss. Those require independent assemblies and
remain a separate validation gate.

The positive and masked arms must be executed in isolated blind environments.
Giving one method both genomes would reveal the hidden target intervals by a
simple sequence diff.

## Why v0.1 was not accepted as the final control

The first coordinate-stable mask covered all 800 selected gene spans. A
post-generation audit found that 134 spans (16.75%) also intersected an exon or
CDS belonging to a gene retained in the perturbed annotation. This is expected
for nested and overlapping plant genes, but it confounds a locus-specific
negative control.

The v0.2 generator therefore uses the blind perturbed GFF3 as a background
exclusion mask. Any target span intersecting a retained exon or CDS is recorded
in evaluator truth and excluded before sequence masking. No candidate output or
hidden WGDI label is used in this decision.

## Frozen construction

| Item | Value |
|---|---:|
| requested paired events | 800 |
| eligible masked events | 666 |
| excluded background-overlap events | 134 |
| chromosomes | 19 |
| assembly length | 816,166,126 bp |
| union masked length | 1,749,428 bp |
| source bases already `N` inside masks | 0 bp |

Key SHA-256 values:

- source genome: `0c24d792c9ee15028ed65e2fe2fc7fbe9627c5fde5d25fb97f5aac5a72d9a8aa`;
- blind perturbed GFF3: `2a4085f95a3d10e9a54b075e478e29bce86a583e81ad95c74dcb2ee96e0c9bbe`;
- evaluator hidden truth: `28e3278b1bc2587625f9b1747d1d235c927034d5bd0dfa476f9bcd12f735e371`;
- v0.2 masked genome: `7143cb6b14debd2573a8286f998fbf28e35694982d694bff965294c02c345364`;
- v0.2 mask truth: `946af5be84002ad3e69717ee3a7d6f51e7c08ade7b95534c105912d557b01a36`.

The independent byte-level audit passed all seven invariants: sequence IDs,
headers, and lengths match; every declared sequence is present; every masked
base is `N`; all bases outside masks are identical; and observed union length
matches evaluator truth. Exactly 1,749,428 bases changed.

All original 24 combinations of subgenome × WGDI stratum × transcript-count
bin remain represented after exclusion. Per-cell retention ranges from 68% to
100%; counts are frozen in the evaluator-only selection summary rather than
silently treating the 666 events as a balanced resample.

## Evaluator sentinels

The scorer defines a false repair as a novel candidate exon or CDS that claims
one or more bases inside a declared all-`N` interval. A transcript whose span
crosses the interval only as an intron is reported separately and is not called
a feature-level false repair.

| Sentinel | Feature-level false repair | Abstention |
|---|---:|---:|
| no-op perturbed annotation | 0/666 | 100% |
| source-annotation oracle | 666/666 | 0% |

The source oracle restores 1,240 hidden transcript structures at these 666
loci and is intentionally treated as maximally wrong on the masked genome.

## Frozen miniprot result

The positive-arm miniprot command, external Brassica proteins, adapter rules,
and thresholds were reused without change. The v0.2 masked arm produced
259,324 raw projections; the frozen adapter accepted 22,579 and rejected
236,745.

| Metric | Result |
|---|---:|
| masked events | 666 |
| events with exon/CDS claim inside mask | 0 |
| feature-level false-repair rate | 0% |
| feature-level abstention rate | 100% |
| unique candidate feature bases inside mask | 0 bp |
| events with transcript-span crossing only | 8 (1.20%) |

This shows that the projection baseline respects absent sequence at the feature
level. It does not solve the positive task: on the unmasked benchmark its
strict precision was 0.0261% and only 5/800 genes were completely recovered.
The combined result strengthens the need for a calibrated decision layer that
can use plant synteny to prioritize real annotation gaps while abstaining on
sequence-unsupported loci.

## Reproducibility

- clean-mask generator commit: `e7b1b0ad8bec429b5710e87c93992362e355e973`;
- byte-level audit commit: `9c503efc29b95635cbeaf5dcc29ebaf720a4b647`;
- selection-summary commit: `22fd859583952de590f35d60505fdfd800c8b406`;
- server benchmark root:
  `/data/codexli/projects/PloidyPatch/benchmark/formal/v0.1/bna_daae_annotation_missing_gene_seed20260807/`;
- server baseline root:
  `/data/codexli/projects/PloidyPatch/results/baselines/miniprot_brassica_v0.1/negative_controls/assembly_gap_v0.2/`.

The blind, evaluator, and baseline directories each contain a `SHA256SUMS`
file. Large genome, index, projection, and candidate files remain on the
server. v0.1 is retained only as an auditable record of the collision finding;
v0.2 is the accepted synthetic assembly-gap control.

