# Homeolog topology coherence v0.2

Date: 2026-08-07  
Status: implemented and evaluated as post-holdout development; review ranking
only, with automatic approval disabled.

## Biological hypothesis

A proposed missing-copy repair should resemble its independently inferred,
already annotated WGD partner in coding topology. This is more specific to the
plant copy-collapse mechanism than generic projection quality or a binary
`in_synteny` feature.

For every candidate with a blind self-WGD link to an existing gene, the
truth-blind feature builder compares every coding isoform pair and retains the
best deterministic pair. It reports:

- CDS-length ratio;
- CDS-segment-count ratio;
- longest-common-subsequence similarity of CDS phases;
- similarity of intron positions expressed as cumulative coding fractions;
- coding-span ratio; and
- their unweighted mean as an auditable topology-coherence summary.

Candidate generation, self-WGD inference, partner selection, and topology
features read only the perturbed annotation and candidate-side evidence.
Hidden truth is joined only by the evaluator after the feature artifacts are
frozen.

## Coverage and raw separation

| Species | All candidates | Topology available | Hidden positives | Positives with topology | Positive/negative topology-score median |
|---|---:|---:|---:|---:|---:|
| soybean | 16,638 | 1,076 | 619 | 590 | 0.9846 / 0.7434 |
| Brassica napus | 57,072 | 1,910 | 607 | 515 | 0.9957 / 0.8079 |
| cotton | 36,849 | 2,095 | 553 | 523 | 0.9983 / 0.7695 |

All five component features separate positives and negatives in the same
direction in all three species. Topology coverage is intentionally restricted
to candidates with an existing blind WGD partner; the system abstains rather
than inventing a relationship for the remaining candidates.

## Leakage correction

An initial exploratory evaluator stacked topology features on a soybean model
score that had itself been fitted on soybean labels. That makes soybean OOF and
the reverse Brassica-to-soybean transfer invalid. The process was stopped, and
its working/frozen artifacts are retained under names containing
`invalid_stacking_leakage`; none of those estimates support a claim.

The corrected evaluator never trains on a precomputed in-sample model score.
For every chromosome fold or training species it refits the original full copy
feature contract, missing-value transforms, categorical expansion, and
logistic model using training rows only. It then compares four fixed models:

1. original v0.1 copy features;
2. original features plus topology availability and within-species normalized
   WGD block context;
3. original features plus topology availability and the five topology
   components; and
4. the topology and normalized-WGD additions together.

The same L2 logistic model with fixed `C=1` is used for every arm. No cotton
label is used for fitting. Uncertainty is a 2,000-replicate chromosome-group
bootstrap implemented by exact group multiplicity weights; its vectorized
result was checked against scikit-learn weighted average precision.

The corrected frozen result is:

`results/copy_collapse/model_development/homeolog_topology_nested_evaluation_v0.2`

## Corrected results

### Chromosome-grouped development OOF

| Species | Original AP | + normalized WGD AP | + topology AP | Topology AP delta (95% CI) |
|---|---:|---:|---:|---:|
| soybean | 0.8608 | 0.8596 | 0.8897 | +0.0289 (0.0109, 0.0500) |
| Brassica napus | 0.6438 | 0.6439 | 0.7071 | +0.0633 (0.0471, 0.0788) |

Normalized WGD context alone has no supported OOF benefit. The topology arm
improves AP in both development species, so the gain is not explained by
merely restating whether a WGD partner exists.

### Cross-species transfer

| Train | Test | Original AP | + topology AP | AP delta (95% CI) |
|---|---|---:|---:|---:|
| soybean | Brassica napus | 0.6007 | 0.6871 | +0.0864 (0.0656, 0.1046) |
| Brassica napus | soybean | 0.8150 | 0.8171 | +0.0021 (-0.0195, 0.0291) |
| soybean + Brassica napus | cotton diagnostic | 0.7579 | 0.8007 | +0.0428 (0.0317, 0.0533) |

Adding normalized WGD context to topology raises the cotton diagnostic AP to
0.8046, a +0.0467 delta (95% CI 0.0359-0.0568), but it does not improve the
Brassica-to-soybean direction. The defensible conclusion is that explicit
homeolog topology supplies a transferable ranking increment in two directions
and no established increment in one direction. It is not a universal gain.

Cotton is a post-holdout diagnostic for v0.2 because its truth had already been
opened after the v0.1 policy freeze. It cannot serve as the untouched external
validation of the new feature family.

## Quality and claim boundary

This result upgrades the algorithmic contribution from a method-vote wrapper
to a plant-mechanistic relationship model. It does not authorize automatic
annotation repair:

- the result evaluates ranking, not a portable precision-guaranteed threshold;
- controlled copy-collapse perturbations dominate the current truth sets;
- candidates without a reliable existing WGD partner remain unresolved; and
- matched SynGAP and natural-error validation are still pending.

The next formal gate is to freeze the topology schema, dependency policy,
model training, and abstention rule, then evaluate them unchanged on an
untouched plant species. A natural-error set supported by long-read
transcripts, assembly comparison, or blinded expert review remains mandatory
for a broad methods claim.
