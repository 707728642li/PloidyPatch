# Public model-plant structure benchmark matrix v0.1

Date frozen: 2026-08-07

This checkpoint validates the multi-event perturbation and evaluator harness
across public plant lineages. It is not a repair-method performance result.

## Frozen design

Four annotations are evaluated independently: *Arabidopsis thaliana* TAIR10,
*Oryza sativa* IRGSP-1.0, *Brassica napus* Da-Ae, and *Glycine max* v2.1. Each
contributes 200 deterministic events for each of four structural errors:
missing internal exon, boundary shift, split gene, and fused gene. The complete
matrix therefore contains 3,200 events.

TAIR10 is a small curated dicot development source. IRGSP-1.0 is the public
monocot holdout. Da-Ae stresses a recent allopolyploid, while soybean provides
an independently held-out paleopolyploid. Copy-collapse remains deferred until
independent WGD/synteny copy pairs are frozen; it is not simulated by treating
arbitrary nearby paralogs as truth.

The TAIR10 and IRGSP inputs are decompressed-only derivatives of the staged,
checksum-verified Ensembl Plants release 62 files. Their exact uncompressed
SHA-256 values are:

| Dataset | SHA-256 |
|---|---|
| TAIR10 GFF3 | `dd1e304a0b3fccb95f5f95ac990b4746ca751ab95e90960f4b40d1d029c96440` |
| IRGSP-1.0 GFF3 | `af797c7eb659ee262e345dc728c4537c1299f43e064a68110952fff18d60b0c5` |

## Integrity result

All 16 dataset/event runs pass their evaluator quality gate. In every run:

- the no-op candidate recovers 0/200 complete CDS chains;
- the oracle source annotation recovers 200/200 complete CDS chains;
- restoration is byte-identical to the input GFF3;
- the expected 200 hidden events are present;
- no non-empty stderr log is produced.

All eight checksum sets in the public TAIR10/rice arm and all eight checksum
sets in the Da-Ae/soybean arm verify successfully. This establishes a usable
cross-lineage test harness but does not establish that PloidyPatch or any
baseline can repair these events.

## TAIR10 missing-annotation development baseline

A separate 800-event TAIR10 missing-annotation run projects only *A. lyrata*
and *A. halleri* proteins; no TAIR10-derived reference protein is supplied.
The paired complete-annotation control and all evaluator gates pass. The
conservative miniprot gap adapter obtains 0.7034 exact-CDS precision, 0.4240
recall, and 0.5291 F1, with 341/800 complete CDS-chain events and 24/800
complete transcript events. CDS nucleotide coverage F1 is 0.9423, again
showing that local coding recovery is substantially easier than full isoform
reconstruction.

Indexing takes 3.1 s, projection 16.4 s, the two adapter arms about 8.3 s each,
and paired scoring 57.8 s on the server. Projection peaks at 1.79 GB RSS and
scoring at 3.45 GB. This is a public model-plant development result for the
candidate layer, not a comparison against the stronger matched methods.

## Frozen artifacts

- public source derivation:
  `/data/codexli/projects/PloidyPatch/data/derived/structure_sources_v0.1`;
- TAIR10/rice matrix:
  `/data/codexli/projects/PloidyPatch/benchmark/structure/public_models_v0.1`;
- Da-Ae/soybean matrix:
  `/data/codexli/projects/PloidyPatch/benchmark/structure/v0.1`;
- public model configuration:
  `config/structure_benchmark_public_models_v0.1.tsv`.
- TAIR10 missing-annotation benchmark:
  `/data/codexli/projects/PloidyPatch/benchmark/public_models_missing_gene/v0.1/ath_tair10_seed20260814`;
- TAIR10 miniprot result:
  `/data/codexli/projects/PloidyPatch/results/baselines/miniprot_v0.18/ath_tair10_missing_gene_v0.1`.

The next gate is to run matched candidate methods on these frozen blind inputs,
then report event-wise sensitivity, false repairs, calibration, and runtime.
