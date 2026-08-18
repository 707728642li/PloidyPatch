# Controlled plant annotation-structure benchmark protocol v0.2

Date frozen: 2026-08-07  
Status: implementation and two-lineage empirical sampling contract frozen

## Scope

Version 0.2 extends the v0.1 missing-gene primitive to five annotation-error
families while preserving the original leakage boundary, immutable source
checksum, evaluator-only truth, and exact byte restoration:

| Event | Blind mutation | Required correction |
|---|---|---|
| `annotation_missing_internal_exon` | remove one internal exon and its exactly coextensive CDS segment | restore the exon/CDS chain |
| `annotation_boundary_shift` | shrink one exact terminal exon/CDS and matching transcript/gene boundary by 15 bp | restore the original terminal boundary |
| `annotation_split_gene` | replace one single-transcript gene by two coordinate-consistent gene hierarchies | reunify the complete source structure |
| `annotation_fused_gene` | join adjacent same-strand single-transcript genes separated by 1–20 kb | recover both source genes and grouping |
| `annotation_copy_collapse` | remove one member of an externally declared copy pair | restore the collapsed copy without confusing it with the retained partner |

These are annotation perturbations. In particular,
`annotation_copy_collapse` does not claim to simulate an assembly sequence
collapse. The distinct `assembly_tandem_copy_collapse` ontology entry remains
unimplemented until genome sequence, annotation coordinates, and copy-depth
evidence can be perturbed together.

## Eligibility and plant-specific safeguards

All events start from coding hierarchies with unique IDs, exactly one parent
per descendant, at least one exon and CDS, and a single transcript for the
v0.2 structural transformations. This deliberately avoids ambiguous
alternative-isoform and shared-feature cases in the first implementation.

- Missing-exon events require at least three exons and CDS segments. The
  removed feature must be internal, with exactly coextensive exon/CDS
  coordinates and no overlapping retained child feature.
- Boundary events require a terminal exon and CDS with equal coordinates, a
  matching gene/transcript boundary, at least 18 bp of length, and no retained
  child extending into the 15-bp trimmed interval. The shift is divisible by
  three so the test does not introduce an unrelated phase rotation.
- Split events require at least four exons and CDS segments. Every direct child
  must lie wholly on one side of a real inter-exon gap, and both products must
  retain CDS.
- Fusion events use genes that are adjacent among all annotated genes, not
  merely adjacent among eligible genes. They must be non-overlapping, on the
  same sequence and strand, with a 1–20 kb intergenic gap.
- Copy-collapse events require an evaluator-supplied TSV with
  `gene_id_a` and `gene_id_b`. PloidyPatch does not infer duplicate truth from
  names or proximity. For WGD tests this file must come from an independently
  frozen synteny/homeology analysis and is never exposed to candidate methods.

Candidate options are ranked by SHA-256 over seed, event type, and stable
source identity. Selection greedily excludes overlapping genes, making output
deterministic for fixed source bytes, event type, count, and seed.

## Reversibility contract

Each hidden event records source line number, source line bytes and SHA-256,
and the exact zero, one, or multiple replacement lines plus their hashes.
Deletion, coordinate rewrite, parent rewrite, one-to-many split, and
many-to-one fusion therefore share one line-edit representation.

`ploidypatch benchmark restore` validates every declared replacement, restores
the original line order, and then requires the complete decompressed source
text SHA-256. Any retained-line modification is detected by the final source
checksum even if it is outside an edited locus. Existing v0.1 missing-gene
truth remains readable through a backward-compatible deletion mapping.
New truth records also freeze the complete perturbed-text SHA-256, which the
v5 scorer checks before building any structure index.

## Scoring contract v5

The evaluator identifies, without using feature IDs for correctness:

1. source transcript structures hidden by each event;
2. erroneous transcript structures introduced into the blind annotation;
3. unaffected blind structures that must be preserved;
4. candidate structures novel relative to the blind input; and
5. reproducible background structures from an optional complete-annotation
   control.

Exact transcript and phased-CDS-chain precision/recall remain the primary
structure metrics. Exon segments, splice junctions, phased CDS segments, and
CDS nucleotide coverage remain secondary metrics. A structural event counts
as completely recovered only when all source structures are present and all
introduced erroneous structures are absent. Exact gene grouping additionally
requires every source gene’s identifier-independent transcript set to occur as
a candidate gene group. Fusion therefore requires two recovered source groups;
split requires one reunited source group.

Collateral loss excludes the event’s intentionally erroneous blind
structures before testing whether unrelated blind structures disappeared.
The report separately exposes the number of events with introduced errors,
complete error removals, and error-removal recall. Missing-gene and copy-
collapse events have no introduced transcript structure, so that denominator
is intentionally undefined rather than reported as perfect.

## CLI examples

```bash
ploidypatch benchmark perturb \
  --gff source.gff3 \
  --output-dir blind/missing_exon_seed17 \
  --truth-dir evaluator/missing_exon_seed17 \
  --event-type annotation_missing_internal_exon \
  --count 100 --seed 17

ploidypatch benchmark perturb \
  --gff source.gff3 \
  --output-dir blind/copy_collapse_seed17 \
  --truth-dir evaluator/copy_collapse_seed17 \
  --event-type annotation_copy_collapse \
  --pair-tsv evaluator/frozen_wgd_pairs.tsv \
  --count 100 --seed 17
```

The v0.2 implementation is a benchmark generator and evaluator, not yet the
PloidyPatch repair engine.

## Frozen empirical matrix v0.1

`config/structure_benchmark_v0.1.tsv` freezes four independently generated
200-event sets in development `B. napus` Da-Ae (seed 20260810) and the same
four sets in held-out `G. max` v2.1 (seed 20260811): missing internal exon,
boundary shift, split gene, and fused gene. Each run must retain source and
code hashes, hidden truth outside the blind directory, a byte-identical
restoration, and no-op/oracle score sentinels before publication.

Copy-collapse is intentionally excluded from this first matrix. It will be
added only after species-specific gene pairs are frozen from independent
WGD/synteny evidence. Inferring pairs from gene names, sequence similarity
alone, or physical proximity would leak the intended plant-duplication answer
into the benchmark.
