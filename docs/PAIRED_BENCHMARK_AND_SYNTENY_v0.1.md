# Paired benchmark and blind synteny checkpoint v0.1

Date: 2026-08-07

## Scope

This checkpoint tests whether progenitor protein projection and blind WGDI
gaps can recover 800 removed protein-coding genes in the Da-Ae *Brassica
napus* annotation. The candidate method sees the perturbed annotation,
assembly, and public *B. rapa*/*B. oleracea* evidence. Hidden truth remains in
the evaluator tree.

The formal candidate was frozen before truth-based scoring. Its checksum is
`744f529a6f09aed04245d3f59884d52a3f6edc5562f98ff0448ff4df19987754`.
The candidate-side freeze checksum is
`631a0b71f7c8ee2bcac6a0f61c13f64a81b1aadee523b9ca0168593224fcdce4`.

## Why paired evaluation is required

The first single-candidate score treated every model absent from the NCBI
source annotation as a synthetic false positive. This is not identifiable on
a real reference: many of those calls are reproducible disagreements between
the source annotation and progenitor proteins, and could include genuine
source-annotation omissions.

Scorer v4 therefore runs the same method on the complete annotation as an
evaluator-only control. It subtracts only identical novel structures produced
in both runs. Source structures themselves are not subtracted as background,
so a projection triggered by a synthetic deletion remains a differential
candidate. The complete control is an evaluation device, not an input that a
deployed repair method may use.

## Frozen formal results

| Method | Differential CDS chains | Exact TP | Precision | Recall | F1 |
|---|---:|---:|---:|---:|---:|
| intact miniprot | 564 | 326 | 0.5780 | 0.3010 | 0.3959 |
| miniprot plus strict local synteny-gap gate | 213 | 104 | 0.4883 | 0.0960 | 0.1605 |

The independent 800-event development set gave the same qualitative result:
miniprot yielded 322/589 exact CDS chains (precision 0.5467, recall 0.3035,
F1 0.3903), whereas the strict synteny gate yielded 114/206 (precision 0.5534,
recall 0.1074, F1 0.1800).

At nucleotide resolution, formal miniprot recovered 63.8% of hidden CDS bases
at 96.7% differential precision. The strict synteny gate retained 21.4% at
89.7% precision. Exact full-transcript recovery remained near zero because
protein alignment does not infer UTR boundaries; exact phased CDS chain is the
appropriate strict endpoint for this component.

## What blind local synteny did and did not achieve

Before protein selection, 10,809 distinct blind gap loci contained 29.5% of
the hidden gene spans at 2.16% locus precision. Protein selection retained
4,383 gap loci, with 21.4% containment recall and 3.86% precision. The local
gate therefore enriched positions but discarded too many true responses and
did not improve paired exact-CDS precision on the formal set.

Decision: strict local synteny is rejected as a mandatory repair gate. It is
retained as an event-context feature and auditable explanation.

## Plant-specific stratification that did replicate

Formal miniprot exact-CDS-chain recall by predeclared WGDI stratum was:

| WGDI stratum | Exact CDS chains | Recall |
|---|---:|---:|
| expected and cross-progenitor support | 141 / 265 | 0.5321 |
| cross-only support | 97 / 273 | 0.3553 |
| expected-only support | 73 / 275 | 0.2655 |
| no WGDI block | 15 / 270 | 0.0556 |

A and C subgenome recalls were 0.3168 and 0.2852, respectively. Single-
transcript genes were easier at the exact-chain level than multi-isoform genes
(0.4042 versus 0.2477 for two transcripts and 0.1848 for three or more).

This supports a plant-specific design in which WGD/subgenome evidence defines
confidence and evolutionary state, while the protein projection remains the
structure generator. It does not support using one local WGDI-gap predicate as
a universal filter.

## Multi-progenitor support exploration

Adapter v2/v3 preserves overlap provenance and links models rejected only as
redundant to the accepted model that absorbed them. In the disjoint development
set, 1,880 accepted projections had support from both progenitor sources.
Among paired differential responses on held-out chromosomes, requiring two
sources increased exact-CDS precision from 0.486 to 0.596 while retaining
0.486 recall. The effect is useful but insufficient for a stand-alone
high-confidence classifier; it remains an evidence column rather than a hard
decision rule.

## Audit artifacts

- Formal paired-v4 freeze:
  `/data/codexli/projects/PloidyPatch/benchmark/formal/v0.1/bna_daae_annotation_missing_gene_seed20260807/evaluator/validation/freeze/paired_v4/SHA256SUMS`
- Freeze checksum:
  `220405d00d6466bec493762f72f55b9c1ec2fae9741f47004150a973af21b36c`
- Formal miniprot paired report:
  `score_miniprot_gap.paired_v4.json`
- Formal stratified report:
  `score_miniprot_gap.paired_stratified_v4.json`
- Formal strict-synteny report:
  `score_synteny_gap_v0.1.paired_v4.json`
- Disjoint development selection checksum:
  `8557df4fdb5354af99bbdeea2671f352a9a5a20d7a679b98f0e6ed06c838a24f`

All 13 files in the paired-v4 freeze passed `sha256sum -c`.

## Remaining validity gates

These results establish a reproducible component baseline, not biological
truth or publication readiness. Required next steps are an independent
species/release, transcript-supported structure refinement, explicit
split/fusion and copy-collapse states, real PAV/loss validation, and matched
comparisons to existing plant annotation-polishing tools.
