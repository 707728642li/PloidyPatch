# Matched Glycine baseline result v0.1

Date frozen: 2026-08-07

GeMoMa 1.9, LiftOn 1.0.11, and SynGAP 1.2.5 receive the same public *Glycine
max* v2.1 target, *G. soja* reference, 800 hidden missing-annotation events,
and paired complete-annotation control. All results pass the evaluator quality
gate and use event-identical paired scoring.

## SynGAP compatibility contract

SynGAP completed both upstream `genblastg` arms successfully. Its added CDS
records write phase as `.`. The strict adapter correctly rejects these by
default; that raw compatibility result is retained separately and is not used
as the comparative score.

The matched v0.2 result enables an explicit normalization that applies only
when every CDS segment of a model has missing phase. It assumes a complete CDS,
sets the first 5-prime segment to phase zero, and derives later phases from
cumulative CDS length. Models with mixed valid/invalid phases remain rejected.
The blind arm normalizes 2,637 models and accepts 1,784 gap candidates; the
complete control normalizes 2,100 and accepts 1,237. This transformation is
recorded per model and in the manifest.

## Strict paired result

| Method | Exact CDS precision | Exact CDS recall | Exact CDS F1 | Complete CDS events | Complete transcript events |
|---|---:|---:|---:|---:|---:|
| GeMoMa 1.9 | 0.7952 | 0.4654 | 0.5872 | 347/800 | 10/800 |
| LiftOn 1.0.11 | 0.8539 | 0.4721 | 0.6080 | 353/800 | 1/800 |
| SynGAP 1.2.5 | 0.5949 | 0.2739 | 0.3751 | 210/800 | 10/800 |

Protein-coding segment recovery is much higher than full-transcript recovery
for every method. This reinforces the scope boundary: comparative projection
is useful candidate evidence, not complete isoform reconstruction.

## Paired event bootstrap

Ten thousand deterministic, event-paired bootstrap replicates give the
following 95% percentile intervals for complete CDS-chain recovery:

| Method | Rate | 95% interval |
|---|---:|---:|
| GeMoMa | 0.4338 | 0.3987–0.4675 |
| LiftOn | 0.4413 | 0.4062–0.4763 |
| SynGAP | 0.2625 | 0.2325–0.2938 |

The LiftOn minus GeMoMa difference is 0.0075, with interval -0.0050 to 0.0200;
these two are not separated by this experiment. GeMoMa minus SynGAP is 0.1713
(0.1437–0.1988), and LiftOn minus SynGAP is 0.1788 (0.1512–0.2075). These are
effect-size intervals, not claims of universal method ranking.

This initial interval resamples paired events within event type. Because the
current benchmark contains only missing-annotation events, a chromosome- or
orthogroup-block bootstrap remains required for the final publication analysis.

## Runtime note

The two SynGAP upstream arms ran concurrently and took 2:37:30 and 2:32:43,
with about 2.17 GB peak RSS each. Phase-normalized adaptation took about 52 s
per arm; strict scoring took 1:41 and 5.75 GB peak RSS. Runtime accounting for
the other methods remains in their frozen resource records.

## Frozen artifacts

- SynGAP upstream and normalized evaluation:
  `/data/codexli/projects/PloidyPatch/results/baselines/syngap_v1.2.5/glycine_v0.1/genblastg`;
- three-method bootstrap:
  `/data/codexli/projects/PloidyPatch/results/baselines/matched_bootstrap/glycine_v0.1`;
- GeMoMa result:
  `/data/codexli/projects/PloidyPatch/results/baselines/gemoma_v1.9/glycine_v0.1/evaluation`;
- LiftOn result:
  `/data/codexli/projects/PloidyPatch/results/baselines/lifton_v1.0.11/glycine_v0.1/evaluation`.
