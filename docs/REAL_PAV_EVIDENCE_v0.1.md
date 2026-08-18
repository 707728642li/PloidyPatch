# Real assembly PAV evidence checkpoint v0.1

## Outcome

The first real-data sequence-absence arm is complete for *Brassica napus*
Da-Ae (target/reference) versus Darmor-bzh v4.1 (query). The strict funnel
retains **24 gene-level candidates in 13 structural events**. These are
high-confidence, WGD-anchored **sequence-absence candidates**, not validated
biological PAV truth.

The main result is methodological: a generic whole-assembly protein screen
would have rejected 233 of 255 candidates because an intact homolog existed
somewhere in the allotetraploid genome. The locus-aware plant screen instead
separates local contradictions, unlocalized-scaffold uncertainty, and expected
off-locus homeolog/paralog evidence.

## Why this is not a duplicate of existing plant PAV software

[GET_PANGENES](https://pmc.ncbi.nlm.nih.gov/articles/PMC10552430/) already
uses whole-genome alignments and gene projection to construct plant pangene
sets, so whole-genome-alignment-based PAV calling is not a novel claim.
[SynGAP](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-024-03359-8)
already uses synteny gaps and protein-to-genome evidence to repair missing
plant annotations. PloidyPatch therefore treats these capabilities as
baselines/components.

This checkpoint tests the narrower project claim: copy-specific absence in a
polyploid must be evaluated jointly with named WGD/subgenome context and must
not be negated by an expected homeolog elsewhere. It also makes unlocalized
scaffold evidence an abstention state instead of silently assigning presence
or absence.

## Frozen inputs

### Da-Ae target bundle

- 19 primary chromosomes, 816,166,126 bp;
- 99,791 genes; 85,091 conservative coding genes with WGDI evidence;
- GFF3 SHA-256:
  `4cd8fb0244ef51d7055d5866f11086c1fee927f5cfbb41509f118f6fb0f29dad`;
- complete genome hash remains in the normalization manifest.

### Darmor-bzh v4.1 query bundle

- PlantGarden chromosome assembly plus all random/unlocalized sequences,
  41 FASTA records in the protein screen;
- genome SHA-256:
  `593ea249c856099f5477a8ba84452e8ff6f8fb5b161d84973287cbe398f8342e`;
- annotation SHA-256:
  `82a5a05c873fd0afac7dc7ac3bc7921be78084a7a9f258fd8aadd7c6d8696b4f`;
- FAI SHA-256:
  `528ec5becf863457d62abae92a34aa1278e02b667502700a6417bcedcb2406b9`.

The public release is documented by
[PlantGarden](https://plantgarden.jp/en/list/t3708/genome/t3708.G001). The
[Ensembl Plants Darmor page](https://plants.ensembl.org/Brassica_napus/Info/Annotation/)
describes the older assembly as incomplete relative to the estimated genome
size. This is why assembly uncertainty is a first-class state here.

## Frozen evidence policy

1. Align Darmor to Da-Ae with minimap2 2.31-r1302 using both `asm10` and
   `asm20`, `--secondary=no`, `-c`, and `--cs=short`. The
   [minimap2 manual](https://lh3.github.io/minimap2/minimap2.html) defines
   these assembly presets and PAF/CIGAR behavior.
2. Retain only primary alignments on the declared homologous chromosome pair,
   MAPQ at least 20.
3. Extract target-only CIGAR deletions at least 100 bp with immediate aligned
   anchors of at least 1 kb on both sides.
4. Require 5 kb of query sequence on each side of the breakpoint with no
   non-ACGT base; truncated or ambiguous flanks are rejected.
5. Require one deletion to cover at least 95% of a gene span and 100% of its
   CDS union.
6. Require the gene candidate in both asm10 and asm20 runs.
7. Search each representative Da-Ae protein against the complete 41-sequence
   Darmor assembly with miniprot 0.18. A substantial hit requires identity at
   least 0.5 and query coverage at least 0.2; intact presence requires coverage
   at least 0.5 and no frameshift or in-frame stop.
8. A substantial hit within 100 kb of the WGA breakpoint is a local
   contradiction or abstention. A hit on an unlocalized sequence forces
   abstention. A hit on another primary chromosome or farther away does not
   negate copy-specific absence.
9. Promote a retained sequence-absence candidate only when WGDI supports the
   expected diploid progenitor for that subgenome.
10. Collapse genes sharing the same source-specific deletion boundaries and
    query breakpoint into one structural event.

The WGA component is deliberately modular. Future matched comparisons should
include [SyRI](https://schneebergerlab.github.io/syri/) and/or
[AnchorWave](https://pmc.ncbi.nlm.nih.gov/articles/PMC8740769/) rather than
assuming minimap2 alone is authoritative.

## Evidence funnel

| Gate | asm10 / count | asm20 / count | Combined interpretation |
|---|---:|---:|---|
| Homologous-pair alignments passing MAPQ/primary filters | 28,347 | 25,655 | Alignment rows, not events |
| Anchored target-only deletions before query-flank rejection | 2,845 | 2,841 | Derived from pass + rejected counts |
| Rejected for ambiguous query flanks | 699 | 690 | Assembly-uncertain |
| Strict clean-flank deletions | 2,146 | 2,151 | Sequence-level interval candidates |
| Coding genes meeting 95% span / 100% CDS coverage | 258 | 257 | Gene-level candidates |
| Supported by both presets |  |  | 255 genes; 5 preset-sensitive genes rejected |

The locus-aware full-assembly protein screen then classified the 255 genes:

| Decision | Genes |
|---|---:|
| Exclude: intact local sequence present | 23 |
| Abstain: disrupted local sequence | 4 |
| Abstain: substantial hit on unlocalized/random sequence | 170 |
| Retain before WGD gate: off-locus homolog only or no substantial homolog | 58 |

The WGD gate classified those 58 retained candidates:

| WGD result | Genes |
|---|---:|
| Expected-progenitor supported; retain high-confidence candidate | 24 |
| Cross-progenitor only; abstain | 1 |
| No WGDI block; abstain | 33 |

## Event-level result

The 24 retained genes collapse to **13 structural events**:

- 4 A-subgenome events and 9 C-subgenome events;
- 9 single-gene, 2 two-gene, 1 three-gene, and 1 eight-gene event;
- deletion sizes 7,539--61,291 bp, median 18,488 bp;
- total nonoverlapping candidate event span 289,579 bp;
- asm10 and asm20 agree exactly on all 13 event boundaries and query
  breakpoints;
- all 24 genes have an intact off-locus homolog, consistent with the expected
  homeolog/paralog-rich polyploid background rather than evidence against the
  copy-specific event.

The largest event is a 61,291 bp C5 candidate containing eight genes. This and
the C-subgenome enrichment are diagnostic flags, not biological conclusions;
both could reflect release-specific assembly behavior.

## Artifacts and freeze

Server result root:
`/data/codexli/projects/PloidyPatch/results/evidence/pav_darmor_daae_v0.1/`.

Primary outputs:

- gene audit table:
  `final/pav_candidates.wgd_annotated.v2.tsv` (255 rows, 92 columns), SHA-256
  `4d4a08e54dd3fed9271f22417c7763f843cb1a9a50d56687e9af1e8da6eeb72b`;
- event catalog:
  `final/pav_events.high_confidence.tsv` (13 rows), SHA-256
  `2934932c92c0505c80d62f4fec6fdffe9733482e6a05de920897b4835a6cac3c`;
- freeze list: `config/pav_evidence_freeze_v0.1.txt`;
- 26-file checksum set: `freeze/v0.1/SHA256SUMS`, SHA-256
  `fc1d1c1547ab41f76b5353a67190344ae9ed48c275c4f8036d523feca085c16f`.

All 26 frozen files pass `sha256sum -c`. The analysis implementation and
freeze machinery are present through Git commit `d111978`.

## What remains before a biological PAV claim

- reproduce the 13 events with an independent WGA/SV method;
- test a newer Darmor release or a third *B. napus* assembly;
- seek read-depth, long-read, or assembly-graph support at both breakpoints;
- compare against GET_PANGENES with matched assemblies and thresholds;
- manually review the 13 event loci, especially the C5 eight-gene event;
- validate that off-locus hits are the expected homeolog/paralog using an
  explicit event graph or gene tree.

Until those gates pass, downstream code and documentation must use
`high_confidence_sequence_absence_candidate`, never `true_PAV` or
`validated_loss`.
