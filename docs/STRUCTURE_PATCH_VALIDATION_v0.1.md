# Exact-topology structure patch validation v0.1

Date frozen: 2026-08-07

## Claim boundary

This checkpoint validates an immutable, exactly reversible coding-model repair
step. It does not claim that protein projection reconstructs UTRs or complete
isoforms. The automatic policy is limited to Tier A fused-gene and split-gene
hypotheses because those two classes had no typed hypotheses in either
complete-annotation control. Terminal-boundary and missing-internal-CDS
hypotheses remain review-only pending independent translation, RNA, splice, or
synteny evidence.

The development target is public TAIR10. The untouched holdout is public
*Oryza sativa* IRGSP-1.0, using the previously frozen phylogenetic support
groups *O. rufipogon*, *O. glaberrima*, and *O. brachyantha*. No rice threshold
was retuned.

## Truth-blind compiler contract

`ploidypatch patch compile-structure` reads only the perturbed annotation,
candidate GFF3, and typed-hypothesis table. Hidden truth is supplied later to
the evaluator and is inaccessible to compilation. Every compiled hypothesis
must satisfy all of the following:

1. its event type was explicitly enabled by the caller;
2. at least two declared support groups reproduce its phased CDS topology;
3. the exact fused/split set relation still holds in the current inputs;
4. referenced gene and transcript IDs still map to the declared hierarchy;
5. each target is a complete single-transcript gene hierarchy without a shared
   external parent; and
6. no selected hypothesis overlaps the source lines of another selected
   hypothesis.

The compiler emits stable `PloidyPatchRepair` gene, mRNA, exon, and CDS lines
plus hypothesis ID, event type, candidate group, support count, and support
group provenance. The generic patch engine then freezes original lines and
their checksums, applies the replacement, and can restore the original GFF3
byte for byte.

## Complete-repair result

Each event class contains 200 hidden perturbations. “Compiled” is the number
selected by the blind Tier A policy. Complete CDS repair requires all hidden
CDS chains, removal of every introduced error structure, and exact CDS-to-gene
grouping.

| Target | Event | Compiled | Complete CDS repair | Event recall | Strict CDS precision | Baseline transcript losses | Exact full-transcript recovery |
|---|---|---:|---:|---:|---:|---:|---:|
| TAIR10 development | fused gene | 42 | 42 / 42 | 21.0% | 100% | 0 | 0 / 42 |
| TAIR10 development | split gene | 98 | 98 / 98 | 49.0% | 100% | 0 | 1 / 98 |
| IRGSP-1.0 holdout | fused gene | 21 | 21 / 21 | 10.5% | 100% | 0 | 0 / 21 |
| IRGSP-1.0 holdout | split gene | 62 | 62 / 62 | 31.0% | 100% | 0 | 6 / 62 |

Thus, every automatically applied hypothesis is a correct complete CDS repair
in both species, with no strict-CDS false positives and no loss of baseline
transcript structures outside the perturbed targets. Recall remains modest by
design: the policy abstains unless independent groups agree on an exact
topology.

Full-transcript recovery is intentionally poor. Replacement exon extents are
the supported CDS spans, so a protein-only method cannot normally reproduce
the source UTR boundaries. Coding-model and isoform-model confidence must
remain separate outputs.

## Reversibility and audit checks

All four score quality gates pass. Every file listed in each frozen
`SHA256SUMS` verifies, all four reverted GFF3 files are byte-identical to their
blind perturbed inputs, and the scorer reports zero
`baseline_transcript_structures_missing_from_candidate`.

The first real-data compilation attempt also exposed a standards-relevant
edge case: TAIR10 and rice legitimately reuse one CDS `ID` across multiple CDS
segments. The hierarchy walker was corrected to require uniqueness only when
an ID is itself a parent, and the regression test now uses a segmented,
repeated CDS ID. The complete local and server suites pass 95/95 tests after
the fix. The frozen runs use code commit `d1c0fc7`.

## Frozen artifacts

- runner:
  `scripts/run_structure_patch_one_v0.1.sh`;
- server result root:
  `/data/codexli/projects/PloidyPatch/results/structure_repairs/exact_topology_v0.1`;
- per-arm inputs, compiler report, edit specification, immutable patch,
  repaired and reverted GFF3, event-level score, resource log, and checksums
  are retained below that root.

## Next validation gates

1. integrate WGD/homeolog and subgenome evidence so a correct structural patch
   cannot collapse retained plant duplicates;
2. add an explicit copy-collapse repair arm on polyploid and paleopolyploid
   targets;
3. add independent transcript evidence before promoting boundary, missing-
   internal, UTR, or complete-isoform repairs; and
4. compare the final event-aware policy with matched mature baselines rather
   than treating projection as the contribution.
