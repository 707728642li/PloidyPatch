# Populus untouched external validation result v0.4

Date: 2026-08-08

## Decision

The preregistered Populus trichocarpa v4.1 external validation is complete and
the formal outcome is **confirmatory pass**. Both hypotheses passed in their
fixed order, all data and sentinel gates passed, and every safety gate passed.
No rule, reference set, threshold, endpoint, or review budget was changed after
the species and protocol were frozen.

This is the first untouched external confirmation of the composite PloidyPatch
v0.4 policy: the frozen v0.3 review ranker plus the label-independent conflict
winner guard.

## Isolation and audit trail

The blind runner was executed in a bubblewrap namespace with network access
disabled. It mounted only the perturbed target, the two candidate-only Salix
references, frozen source, frozen environments, frozen model and protocol.
The complete Populus annotation, evaluator-only Manihot and Ricinus references,
truth labels, `/nas_data`, and evaluator directories were not mounted.

The blind output tree was checksum-frozen before evaluator access. Independent
verification of all 294 listed files passed. The evaluator then generated the
complete-annotation controls and labels, ran the fixed-sequence tests, and
atomically froze the reveal result.

- protocol `SHA256SUMS`: `1334dc35b709cb047f02eca63cbaea8fc5c77bce5229d9e09e401274b07ffdb8`
- composite model `SHA256SUMS`: `90a5d90b6056f883df32af5556d906849054b0cc1880537cd4c35bd07a3a5797`
- execution patch2 `SHA256SUMS`: `9415b727d1d1cf93cfdbc98d9982742ef12833dd37de39de5274dbe6d27807b9`
- blind tree `SHA256SUMS`: `c1fef1541a9361ec42ae5220487e11400fbecd9345477f14f7c35edb0b7d3f4b`
- custody manifest: `edb1df94710525fdcc80f872de98ad6386e98f3d9102e35c90f6e65d3ade9955`
- reveal `SHA256SUMS`: `796b660950ee68f2e097f94a06cdd164b5fae1bbdc0ab257e63e6a376c899785`
- final `evaluation.json`: `0e51f5bde4eb78a11c74179215431254eab938f21684ff2dae2838f619800f37`

The two failed execution attempts are retained read-only. They failed before
candidate scores or labels existed: the first exposed missing `/bin/sh` and a
provider GFF without exon rows; the second exposed LiftOn exon-ID collisions.
Patch2 fixed the adapter deterministically without changing the scientific
protocol. The successful run used commit
`d4241a01ba35b1bac23c4aab568d7e525ef855f0`.

## Evaluability

The evaluator requested and obtained 800 WGD copy-collapse events across all
19 target chromosomes. The four CDS-complexity bins contained 202, 213, 188,
and 197 events, respectively. No pair or reference rule was relaxed.

All sentinels passed:

- blind no-op exact recovery: 0;
- complete-annotation oracle exact recovery: 800/800;
- restored GFF byte-identical to the source annotation;
- blind and complete target genome SHA256 identical;
- pair-truth audit passed.

## Fixed-sequence confirmatory results

| Test | Comparator | PloidyPatch v0.4 | Observed difference | 95% bootstrap interval | Decision |
| --- | ---: | ---: | ---: | ---: | --- |
| H1 exact phased-CDS event recall | legacy suppress-overlap pool: 663/800 (0.82875) | chain-preserving pool: 731/800 (0.91375) | +0.08500 | +0.06625 to +0.10500 | pass |
| H2 candidate AP | frozen non-topology baseline: 0.845752 | v0.4 guard: 0.853194 | +0.007442 | +0.004210 to +0.010761 | pass |

H1 used 20,000 paired event-bootstrap replicates with seed 20261001. H2
was tested only after H1 passed and used 20,000/20,000 valid target-chromosome
bootstrap replicates with seed 20261002.

The chain-preserving pool recovered 68 additional exact structures. Relative
to the old pool, that is a 10.26% increase in recovered events. Its unranked
strict precision was 0.6521, which is expected for a recall-oriented candidate
ceiling and confirms that ranking and manual review remain essential.

## Ranking and review behavior

There were 8,732 candidates: 731 exact-CDS positives and 8,001 negatives.
The frozen v0.3 primary ranker reached AP 0.853219 and the v0.4 guard reached
0.853194. The guard therefore retained 99.67% of the v0.3 AP gain over baseline;
its descriptive difference from v0.3 was -0.000025 with a 95% interval spanning
zero. The safety guard removed essentially none of the transferable efficacy.

| Review budget | Reviewed | Exact positives | Precision | Positive-candidate recall |
| --- | ---: | ---: | ---: | ---: |
| top 0.5% | 44 | 42 | 0.9545 | 0.0575 |
| top 1% | 88 | 86 | 0.9773 | 0.1176 |
| top 2% | 175 | 172 | 0.9829 | 0.2353 |
| top 100 | 100 | 98 | 0.9800 | 0.1341 |
| top 250 | 250 | 246 | 0.9840 | 0.3365 |
| top 500 | 500 | 472 | 0.9440 | 0.6457 |

These numbers describe review prioritization on the copy-collapse benchmark;
the score is not a calibrated biological probability and does not authorize
automatic annotation changes.

## Safety gates

All preregistered gates passed:

- every one of 1,783 conflict-set winners was byte-for-byte identical to the
  deterministic baseline winner (`mismatch_count = 0`);
- the top-1% exact-positive yield was noninferior to baseline;
- topology-positive coverage was 0.7510, above the frozen 0.70 floor;
- 99.67% of the positive v0.3 AP gain was retained, above the 90% floor;
- primary, legacy, miniprot, GeMoMa, LiftOn, support2 and support3 arms each
  had zero collateral loss of baseline transcript structures;
- 20,000/20,000 H2 bootstrap replicates were valid;
- automatic approvals were zero.

Only 9 candidates required topology abstention. On the 226 conflict sets with
exactly one positive, baseline, v0.3 and v0.4 all selected 138 correct winners
and had mean reciprocal rank 0.7994. Thus this holdout did not need the guard
to repair a realized conflict winner, but the frozen guard guarantee remained
active and passed its exact mapping audit.

## Scientific interpretation

The result confirms two distinct claims on an untouched plant lineage:

1. preserving distinct overlapping phased CDS chains materially raises the
   ceiling for recovering collapsed WGD/homeolog copies; and
2. the frozen plant-specific evidence ranker transfers enough signal to improve
   review prioritization, while the v0.4 guard enforces the conflict-safety
   invariant with negligible loss of AP.

Together with the earlier apple v0.3 external efficacy pass, the retrospective
cotton comparator, development transfers, and maize natural-read enrichment,
this is now strong methods-paper evidence. The remaining gap to a complete
Genome Biology submission is no longer external algorithmic confirmation. It
is biological breadth and presentation: fully audited raw-read natural cases,
one persuasive end-to-end user workflow, release-quality packaging, figures,
and a manuscript that keeps the claim boundary explicit.

## Frozen server artifacts

- blind run: `/data/codexli/projects/PloidyPatch/results/blind_runs/populus_external_v0.4`
- evaluator reveal inputs: `/data/codexli/projects/PloidyPatch/results/evaluator/populus/v0.4/reveal_inputs`
- final reveal: `/data/codexli/projects/PloidyPatch/results/copy_collapse/external/populus_v0.4_reveal`

