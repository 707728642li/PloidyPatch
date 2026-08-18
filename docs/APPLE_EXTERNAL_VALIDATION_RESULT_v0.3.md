# Apple external validation result v0.3

Date: 2026-08-08

## Decision

The frozen apple GDDH13 v1.1 external evaluation is formally evaluable and
both ordered efficacy hypotheses passed. The omnibus confirmatory decision is
nevertheless **fail**, because one predeclared conflict-set safety gate was
not met. This negative gate result is retained without apple retuning.

The result supports chain-preserving candidate generation and a modest but
portable average-precision improvement. It does not support automatic copy
addition or a claim that the support-conditioned correction is uniformly
safer than its baseline inside mutually exclusive candidate conflicts.

## Frozen external test

- Strict evaluator truth contained 5,291 unordered apple duplicate pairs: the
  exact intersection of 8,631 apple self-WGDI pairs and 6,273 pairs supported
  by both rose and strawberry.
- The frozen sampler selected 800 events on all 17 chromosomes. CDS-complexity
  bins contained 204 one-segment, 204 two-to-three-segment, 204 four-to-six-
  segment and 188 seven-plus-segment events.
- All formal evaluability gates passed. The no-op recovered 0/800 chains, the
  complete oracle recovered 800/800, restoration was byte-identical, and the
  pair-truth audit passed.
- The chain-preserving pool contained 30,511 candidates and 4,903 conflict
  sets. The legacy suppressing pool contained 24,629 candidates.
- Exact hidden-chain positives numbered 545; 29,966 candidates were negative.

## Ordered hypotheses

### H1: chain-preserving candidate ceiling

The primary pool recovered 545/800 events (recall 0.68125), versus 438/800
(0.54750) for the suppressing pool. The paired gain was +0.13375, with a
20,000-replicate event-bootstrap 95% interval of 0.1100 to 0.1575. All 438
legacy recoveries were retained and 107 additional events were recovered.
H1 passed.

The unranked primary pool is intentionally not an approved patch set: its
strict-chain precision was 0.3972 and F1 was 0.5018, versus 0.5516 and 0.5496
for the legacy pool. The gain is an information ceiling that requires ranking.

### H2: frozen v0.3 ranking

The frozen non-topology baseline reached AP 0.48016 and the primary
support-conditioned ranker reached AP 0.49142. The difference was +0.01127,
with a 20,000-replicate target-chromosome bootstrap 95% interval of 0.00275 to
0.01920. All 20,000 replicates were valid. H2 was tested after H1 passed and
also passed.

The primary ranker returned 73 exact chains in its top 100, 191 in the top 1%
(306 reviews), and 280 in its top 500. The baseline returned 72, 188 and 272,
respectively. At top 0.5%, the primary returned 103 versus 105 for baseline;
this was secondary and not the frozen review-budget safety gate.

The secondary frozen v0.2 topology ranker reached AP 0.48807. Method-support
count reached 0.22267 and maximum raw method quality reached 0.08438. These
secondary comparisons do not replace H2.

## Safety gates

| Gate | Result | Evidence |
| --- | --- | --- |
| formal evaluability | pass | 800 events, 17 chromosomes, every complexity bin >=188 |
| H1 candidate ceiling | pass | +0.13375; CI 0.1100 to 0.1575 |
| H2 AP in fixed sequence | pass | +0.01127; CI 0.00275 to 0.01920 |
| topology-positive coverage | pass | 393/545 = 0.7211, above 0.70 |
| top-1% review noninferiority | pass | 191 versus 188 exact candidates |
| zero collateral loss | pass | zero for all seven candidate arms |
| no automatic approval | pass | zero approved candidates |
| conflict-set top-1 noninferiority | **fail** | 163/269 versus baseline 165/269 |

The conflict failure was narrow but directional. Among 269 conflict sets with
exactly one positive, both estimators were correct in 163 and both were wrong
in 104. Baseline alone was correct in two sets; primary alone was correct in
none. A non-preregistered paired exact McNemar test gives P=0.5, so the data do
not establish a meaningful degradation, but the frozen zero-margin gate was
literal and therefore failed. No inferential substitution is made.

## Candidate-arm diagnostics

| Arm | Exact events | Recall | Precision | F1 |
| --- | ---: | ---: | ---: | ---: |
| Miniprot | 338 | 0.4225 | 0.4721 | 0.4459 |
| GeMoMa | 478 | 0.5975 | 0.6105 | 0.6039 |
| LiftOn | 374 | 0.4675 | 0.5027 | 0.4845 |
| support >=2 | 401 | 0.5013 | 0.6914 | 0.5812 |
| support =3 | 244 | 0.3050 | 0.8385 | 0.4473 |
| legacy union | 438 | 0.5475 | 0.5516 | 0.5496 |
| chain-preserving union | 545 | 0.6813 | 0.3972 | 0.5018 |

Topology was available for 1,418/30,511 candidates, including 393/545
positives. Within topology-available candidates, AP increased from 0.59518 to
0.60819. Topology-unavailable candidates received identical primary and
baseline scores by construction. AP increased on 11 of 17 chromosomes.

## Execution amendments and audit trail

No apple pair or label existed when execution v0.3 first encountered the
provider `gene:` prefix versus exact `Name` mismatch. Execution v0.3.1 added a
strict evaluator-only alias derivative. After truth, events and the candidate
pool were checksum-frozen, the blind mixed provider/PloidyPatch GFF exposed the
same identifier issue for candidate `ID` values. Execution v0.3.2 generalized
the adapter to a unique exact match among `ID`, `gene_id`, `locus_tag` and
`Name`. It mapped 73,199/73,199 blind WGDI representatives without coordinate
change. Both failed working directories and earlier execution freezes were
retained. No candidate label or model score existed at either algorithmic
amendment.

The first reveal attempt stopped before evaluation because the model Python
environment lacked the project source path. The failed directory was retained
and the same frozen scripts were rerun with an explicit `PYTHONPATH`; no code,
score, label or endpoint changed.

Final hashes:

- outer reveal `SHA256SUMS`: `605ed025b9bd70098218820b4f6f2057158566237b8e83c652c2c9a9e0fc42c9`
- evaluator `SHA256SUMS`: `7f15f58b0e4ba7431fe1aaa063ad61f5644eca4300491ef64e8da09f3281d2d3`
- `evaluation.json`: `f1df1c523ce7e4a5303ca8b09fe4fd8789cffd288e82ad8a1cb3e24292425458`

## Publication interpretation

The fresh species confirms the two central efficacy mechanisms: preserving
alternative exact chains substantially raises recoverable recall, and the
frozen plant-WGD correction improves global candidate ranking. The strict
omnibus outcome remains negative because the correction changed two conflict
winners in the wrong direction. The defensible current claim is therefore a
strong plant-specific review-prioritization framework with independently
replicated efficacy, not a fully safety-confirmed automatic repair system.

A next estimator must be versioned as v0.4 and tested on another untouched
external species. A natural design is to retain v0.3 global ordering while
constraining within-conflict ordering to the proven baseline or abstaining
when topology changes a conflict winner. Apple can be development evidence for
that future version but cannot be reused as its confirmatory test.
