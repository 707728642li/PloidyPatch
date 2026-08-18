# Apple external validation protocol v0.3

Date frozen: 2026-08-08  
Policy ID: `ploidypatch_apple_external_validation_v0.3`

## Confirmatory question

This is the first untouched external test of the frozen support-conditioned
PloidyPatch v0.3 ranker. The target is apple GDDH13 v1.1. Soybean and Brassica
were used for development; cotton and maize labels were already viewed and are
retrospective diagnostics. No apple duplicate pair, perturbation event,
candidate label, model score or endpoint was inspected before this protocol.
Only file availability, primary sequence names, annotation syntax and feature
counts were checked.

The confirmatory claim has two ordered parts:

1. preserving distinct, strongly overlapping phased-CDS chains increases the
   event-level recovery ceiling over the legacy suppressing union; and
2. conditional on that frozen candidate pool, the frozen v0.3
   support-conditioned topology correction improves candidate average
   precision over its frozen non-topology baseline.

The tool remains a review-prioritization system. No score is a calibrated
probability and no copy addition is automatically approved.

## Frozen data roles

| Role | Assembly/annotation | Permitted use |
|---|---|---|
| target | *Malus x domestica* GDDH13 v1.1 | complete evaluator source, blind target genome and annotation |
| candidate reference | *Pyrus communis* BartlettDH v2.0 | miniprot, GeMoMa and LiftOn candidate generation |
| candidate reference | *Prunus persica* NCBIv2, Ensembl Plants 62 | miniprot, GeMoMa and LiftOn candidate generation |
| evaluator reference | *Rosa chinensis* RchiOBHm-V2, Ensembl Plants 62 | duplicate-pair truth support only |
| evaluator reference | *Fragaria vesca* v4.0.a1 | duplicate-pair truth support only |

Candidate references may not enter the truth definition. Evaluator references
may not generate or rank candidates. These roles cannot be changed after the
freeze, even if pair yield or model performance is disappointing.

The copied public bundle and normalized primary-chromosome bundle are
identified by the SHA-256 values in
`config/apple_external_validation_policy_v0.3.tsv`. Provider compatibility was
performed before truth construction and was limited to three escaped apple
`Note` continuations, five dropped invalid pear `intron` records and one
stripped strawberry `##FASTA` section. No gene, transcript, exon or CDS
coordinate was repaired or removed.

## Evaluator-only duplicate truth

### Protein and WGDI preparation

Proteins for apple, rose and strawberry are reconstructed from their frozen
complete primary-chromosome GFF3 and genome FASTA with gffread 0.12.9. One
representative protein per operable gene is prepared by the existing
PloidyPatch WGDI input command. Preparation does not enumerate pairs and does
not read candidate outputs.

DIAMOND 2.2.2 and WGDI 0.75 use the already established parameters:
`evalue=1e-5`, at most 20 DIAMOND targets, `--more-sensitive`, WGDI
`score=100`, `grading=50,40,25`, `mg=40,40`, `pvalue=0.2`,
`repeat_number=20`, gene-order positions, and 64 processes. All commands,
inputs, outputs and versions are hashed.

### Independent evidence requirements

An apple pair must pass all three evidence streams:

1. apple self-WGDI: cross-chromosome, at least 20 gene pairs in the supporting
   block, and reciprocal-unique after filtering;
2. apple versus rose: one rose counterpart is linked in qualifying collinear
   blocks to exactly two apple genes; and
3. apple versus strawberry: the analogous exactly-one-to-two relation.

The outgroup inference requires both evaluator species as independent support
groups, cross-chromosome apple genes, and reciprocal uniqueness. The final
eligible truth set is the exact unordered apple feature-ID pair intersection
between accepted self-WGDI pairs and accepted two-outgroup pairs. It is a
conservative syntenic duplicate truth set associated with the Maleae WGD, not
a gene-tree proof of homeology.

No pear or peach relationship can support a truth pair. The pair table and all
pair decisions remain evaluator-only and are inaccessible to candidate
generation, blind self-synteny, feature construction and ranking.

## Hidden-event sampling and sentinels

The evaluator requests 800 events with seed `20260831`. If fewer than 800
eligible pairs remain, all eligible pairs are used and the shortfall is
reported. Sampling uses global target-chromosome round robin, then round robin
over one, two-to-three, four-to-six and seven-plus CDS-segment bins, then a
SHA-256 rank. No gene can occur in two selected pairs. The same deterministic
hash rule chooses exactly one pair member to hide while the assembly sequence
remains unchanged.

The test cannot make a confirmatory pass claim unless there are at least 500
events, at least 15 target chromosomes, and at least 20 events in every CDS
complexity bin. These are evaluability gates, not reasons to change references
or relax pair rules.

Before any baseline is evaluated:

- the blind annotation must recover zero complete hidden phased-CDS chains;
- the complete annotation oracle must recover every one;
- applying evaluator restoration must reproduce the complete source GFF3
  byte-for-byte; and
- blind and complete arms must use byte-identical target genome sequence.

Failure of a sentinel invalidates the run rather than counting as poor model
performance.

## Truth-blind candidate generation

The three method families are miniprot 0.18-r281, GeMoMa 1.9 and LiftOn
1.0.11. Pear and peach are run separately for GeMoMa and LiftOn and merged
within method family before voting; two references never count as two
algorithmic votes. Miniprot projects a namespaced combined protein set.

Upstream parameters and adapters are frozen to those used in prior species.
Miniprot uses `-I -N4 --outn 4 --outs 0.8 --outc 0.5`. GeMoMa uses the
MMseqs-backed pipeline, no RNA evidence, `GeMoMa.Score=ReAlign`, and no
AnnotationFinalizer redundancy filter. LiftOn uses gene-only, locus-pipeline
and output-validation modes. The PloidyPatch adapters require identity at
least 0.5, query coverage at least 0.5, intact projections, existing-CDS
overlap at most 0.2, and candidate redundancy overlap at most 0.5.

Every raw prediction is adapted independently against the blind annotation
and paired complete-annotation control. Differential evaluation subtracts
complete-control background; the complete target annotation cannot influence
the blind upstream prediction, adapter threshold or candidate rank.

Exact chains are keyed by target seqid, strand and ordered phased CDS
coordinates. The primary pool deduplicates exact chains but retains distinct
strongly overlapping alternatives in explicit conflict sets. The legacy
comparator uses the same raw predictions and thresholds but suppresses such
alternatives. Both pools are frozen before truth is read.

## Blind plant-specific features and frozen scoring

Blind WGD context is recomputed from the blind candidate union and target
genome. It cannot read the complete annotation or pre-perturbation truth pair
table. Candidate-to-candidate pairs are rejected as circular evidence; only a
unique existing-gene partner can propagate within a conflict set.

Copy features, topology features, candidate GFF, conflict sets and their
checksums are frozen before the portable model is invoked. The model SHA-256
is `9dcafcc3294bfd85ccfbe824c33debd0e967578b89cf7f08b56076d32efd4390`.
It emits the frozen baseline logit and the support-conditioned primary score.
Topology-unavailable candidates must receive exactly the baseline score.
Scores and ranks are frozen before evaluator labels are joined.

## Ordered primary hypotheses

### H1: chain-preserving recovery ceiling

For each hidden event, success means that at least one differential candidate
has the exact hidden phased CDS chain. The paired difference is
`retain_distinct_chains - suppress_overlapping`. A 20,000-replicate paired
event bootstrap reports a two-sided percentile 95% interval. H1 passes only
when the point difference and interval lower bound are both greater than zero.

### H2: support-conditioned ranking

Candidate rows in the frozen chain-preserving pool are positive only when
their exact phased CDS chain matches hidden truth after paired-control
subtraction. The primary metric is average precision of the v0.3 primary
score minus average precision of its frozen baseline score. Uncertainty uses
20,000 target-chromosome grouped bootstrap replicates with seed `20260901`;
all candidates from a sampled chromosome remain together. H2 is tested only
after H1 passes and passes only when the point difference and 95% interval
lower bound are both greater than zero.

The confirmatory algorithmic claim requires both H1 and H2. This fixed-order
procedure prevents choosing whichever endpoint is favorable after reveal.

## Safety gates and secondary endpoints

The result must also satisfy:

- at least 70% of positive candidate rows have blind topology available;
- among conflict sets with exactly one positive chain, primary top-1 accuracy
  is not below baseline top-1 accuracy;
- at the top 1% review budget, the primary score returns at least as many true
  candidates as baseline; and
- all candidate arms retain every baseline transcript structure.

Secondary, non-substitutable results include top 0.5%, 1%, 2%, 100, 250 and
500 review yields; per-chromosome and CDS-complexity results; individual
miniprot, GeMoMa and LiftOn recovery; exact support-two and support-three
consensus; method-support and maximum-quality heuristics; the frozen v0.2
topology ranker; precision-recall curves; candidate ceiling; and failure-mode
strata. Ties use candidate digest as the deterministic secondary key.

All candidate counts, positives, topology coverage, conflicts, coefficients,
software versions, checksums and bootstrap seeds must be reported regardless
of direction. No threshold, feature, pair rule, reference role, sampler,
endpoint or gate may change after labels are visible. Failure is retained as
a formal negative external result; another attempt requires a new model
version and a different untouched external system.
