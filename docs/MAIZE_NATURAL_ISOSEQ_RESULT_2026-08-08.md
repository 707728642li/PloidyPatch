# Maize natural Iso-Seq validation result (2026-08-08)

## Result

The RNA-blind candidate and ranking freeze was completed before accessing
validation sequence. From 75,118 published IsoPhase transcript models, the
predeclared pure-B73 count rule (`EM1 + R1 + END1 >= 2`) selected 35,733
transcripts comprising 81,658,050 bases. After de novo alignment to NAM-5.0,
35,257 passed query coverage, identity, mapping-quality, near-best secondary,
and splice-junction filters.

Among 14,853 frozen natural candidates, independent evidence states were:

| Evidence state | Candidates |
| --- | ---: |
| full multi-exon CDS chain supported by one transcript | 198 |
| all junctions supported without single-transcript connectivity | 20 |
| partial junction support | 65 |
| single-exon span support | 281 |
| no qualifying observation | 14,249 |
| not assessable because of assembly N in locus/5-kb flank | 40 |

Only the first state is the primary positive. Non-observation is unlabeled,
not a false positive, and single-exon transcription does not establish gene
identity.

## Frozen top-K yield

| Ranker | K=25 | K=50 | K=100 (primary) | K=200 |
| --- | ---: | ---: | ---: | ---: |
| paired baseline | **5** | **8** | **20** | **38** |
| v0.2 topology | 2 | 7 | 19 | 37 |

For the baseline top 100, the chromosome-stratified random-ranking mean was
1.315 full-chain-supported candidates; the observed 20 represents 15.21-fold
enrichment. None of 20,000 random replicates reached the observed yield
(add-one empirical p = 1/20,001). Enrichment was similarly strong at every
declared baseline budget (12.23- to 16.15-fold).

The natural topology-direction gate failed: topology recovered 19 versus the
baseline's 20 in the top 100. The observed difference was -1 supported locus
per 100 reviews, with a 20,000-replicate target-chromosome bootstrap interval
of -4 to +2. This agrees with the synthetic maize result that v0.2 topology
does not improve the portable baseline.

## Interpretation and remaining audits

This is strong evidence that PloidyPatch's baseline ranking concentrates real,
independently observed multi-exon coding structures in a current plant genome.
It is not yet evidence that all 20 loci are annotation errors or safe patches.
Before a correction claim, each supported locus still requires target-ORF,
same-strand annotation collision, repeat/TE, candidate-conflict, multi-locus,
and leave-one-candidate-out homeolog audits. Manuscript case studies must also
be rechecked against raw E-MTAB-7837 reads or deposited per-read alignments.

The result therefore supports a useful review-prioritization claim and raises
the publication ceiling, while leaving automatic annotation correction and
v0.2 topology superiority explicitly unsupported.
