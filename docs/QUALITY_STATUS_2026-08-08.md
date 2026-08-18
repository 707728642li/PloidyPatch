# PloidyPatch quality status (2026-08-08)

## Claim boundary

PloidyPatch is currently supported as a plant-specific, review-oriented system
for recovering collapsed WGD/homeolog coding copies. It is not yet supported as
a general gene annotator, a calibrated probability model, or an automatic
annotation-approval system. The primary structural endpoint is an exact phased
CDS chain after subtraction of projection background measured on a paired
complete-annotation control.

## Engineering and leakage controls

- The local test suite contains 163 passing tests. The apple evaluator was
  additionally exercised in the server model environment with NumPy 2.4.6
  and scikit-learn 1.9.0.
- Candidate generation cannot read hidden truth or the pre-perturbation WGD
  pair table. Candidate artifacts are checksum-frozen before evaluator access.
- Every repair score uses a method-matched complete-annotation control and
  requires zero loss of baseline transcript structures.
- The v0.2 model is a ranker only: no portable probability, threshold, or
  automatic approval is claimed.
- The maize external protocol froze model, source hashes, pair rules, sampler,
  metrics, review budgets, 20,000-replicate chromosome bootstrap, and failure
  policy before hidden labels were produced.

## Established benchmark results

### Frozen v0.1 model on the untouched cotton holdout (800 events)

The soybean-trained rank model reached AP 0.7173 and ROC AUC 0.9888. At its
frozen review decision it recovered 414/800 exact CDS chains with 74 false
positives: precision 0.8484, recall 0.5175, and F1 0.6429. It was substantially
more precise but lower-recall than the best exact two-of-three method consensus
(precision 0.7734, recall 0.6613, F1 0.7129). This is why the absolute soybean
threshold was retired and v0.2 is explicitly rank-only.

### External comparator on the same cotton holdout

SynGAP 1.2.5 (two cotton diploid references, zero retuning) recovered 300/800
exact CDS chains with 237 false positives: precision 0.5587, recall 0.3750, and
F1 0.4488. The frozen PloidyPatch rank decision exceeded SynGAP by 0.1941 F1
(95% bootstrap interval 0.1489-0.2387), 0.2897 precision
(0.2368-0.3424), and 0.1425 recall (0.0935-0.1908).

### Corrected topology evaluation

All transforms and estimators were refitted inside each development fold or
training-species partition. Adding explicit candidate-to-existing-homeolog
coding topology changed AP as follows:

| Evaluation | Baseline AP | Topology AP | Delta AP | 95% group-bootstrap interval |
| --- | ---: | ---: | ---: | ---: |
| Soybean chromosome-grouped OOF | 0.8608 | 0.8897 | +0.0289 | 0.0109 to 0.0500 |
| Brassica chromosome-grouped OOF | 0.6438 | 0.7071 | +0.0633 | 0.0471 to 0.0788 |
| Soybean to Brassica transfer | - | - | +0.0864 | 0.0656 to 0.1046 |
| Brassica to soybean transfer | - | - | +0.0021 | interval includes zero |
| Pooled development to cotton diagnostic | 0.7579 | 0.8007 | +0.0428 | 0.0317 to 0.0533 |

The gain is real but asymmetric. Normalized WGD-block context by itself was
null, so it was excluded from the primary frozen v0.2 estimator.

## Formal maize v0.2 result

The untouched maize B73 NAM-5.0 test is complete. Among 15,636 candidates,
the paired baseline ranker reached AP 0.6529 and the topology ranker reached
AP 0.6511. Topology-minus-baseline AP was -0.0017, with a predeclared 20,000-
replicate chromosome-bootstrap 95% interval of -0.0089 to 0.0045. Topology
coverage among exact positives was 0.7473 and passed its 0.70 floor, but the
two performance requirements failed. This is the retained formal v0.2 result;
maize will not be used to retune or replace it.

GeMoMa was the strongest single maize method (precision 0.7808, recall 0.7612,
F1 0.7709). Exact 2-of-3 consensus was more precise but less sensitive
(precision 0.8317, recall 0.6300, F1 0.7169). The old support-at-least-one
selector had F1 0.7137 and, because it suppressed overlapping alternative CDS
chains, discarded 87 exact structures recovered by at least one method. This
pre-ranking information loss is now the primary v0.3 engineering target.

The natural NAM-5.0 candidate and review ranking arm was checksum-frozen before
Iso-Seq access. Of 14,853 candidates, 198 have independent single-transcript
full multi-exon chain support. The baseline top 100 contains 20 such candidates
versus a chromosome-stratified random mean of 1.315 (15.21-fold enrichment;
20,000-replicate add-one empirical p = 1/20,001). Topology contains 19; its
natural difference from baseline is -1 with chromosome-bootstrap interval
-4 to +2, independently confirming that v0.2 topology does not improve maize.

The v0.3 chain-preserving pool restores exact-chain recall relative to the old
suppressing union in soybean (+0.0138), Brassica (+0.0575), cotton (+0.1075),
and maize (+0.1088); every paired 95% event-bootstrap interval excludes zero.
This raises the candidate ceiling but lowers unranked precision, so conflict-
aware ranking and abstention remain mandatory.

### Audited v0.3 conflict-aware ranker

The preregistered support-conditioned offset ranker passed all development
retention criteria after true repeated chromosome-grouped evaluation. Soybean
AP increased from 0.8534 to 0.8766 (delta +0.0231; 95% chromosome-bootstrap
interval 0.0106 to 0.0389), and Brassica AP increased from 0.6136 to 0.6605
(delta +0.0469; interval 0.0324 to 0.0597). Conflict-set top-1 accuracy also
increased in both species. All five distinct chromosome partitions per species
had positive AP differences. Pooled development AP increased by 0.0331
(interval 0.0205 to 0.0448).

Cross-species development transfers were positive in both directions. The
pooled development model improved retrospective cotton AP by 0.0376, but only
improved maize AP by 0.0027 and slightly reduced exact counts at fixed top
review budgets. The v0.3 estimator is therefore retained for a newly
preregistered external freeze, not declared portable or confirmatory.

## Formal apple v0.3 external result

The untouched apple GDDH13 v1.1 test is complete and formally evaluable. A
strict intersection of apple self-WGDI and independent rose plus strawberry
support yielded 5,291 duplicate pairs; the frozen sampler selected 800 events
across all 17 chromosomes and all sentinel tests passed.

The chain-preserving pool recovered 545/800 exact chains versus 438/800 for
the legacy suppressing pool, a recall gain of 0.13375 with 20,000-replicate
paired-bootstrap interval 0.1100 to 0.1575. The frozen v0.3 primary ranker
reached AP 0.49142 versus 0.48016 for its baseline, a gain of 0.01127 with
20,000-replicate chromosome-bootstrap interval 0.00275 to 0.01920. Thus both
ordered efficacy hypotheses passed on the fresh species.

The omnibus confirmatory decision is nevertheless **fail**. Conflict-set
top-1 accuracy was 163/269 (0.60595) for primary versus 165/269 (0.61338) for
baseline, failing the literal zero-margin noninferiority gate. The other safety
gates passed: topology-positive coverage was 0.7211, primary returned 191
versus 188 positives at the top 1% review budget, all seven candidate arms had
zero collateral transcript loss, and no candidate was automatically approved.
This failure is retained without apple retuning.

## v0.4 conflict winner guard development result

The label-independent v0.4 guard falls back to baseline scores and abstains
from topology only when v0.3 would change the top-ranked member of a conflict
set. It therefore makes every conflict-set winner exactly identical to the
baseline while retaining v0.3 scores everywhere else. All 163 tests pass and
automatic approval remains disabled.

The fixed development gates passed in soybean, Brassica and label-seen apple.
Guard AP was 0.87638 versus baseline 0.85344 in soybean (delta +0.02294; 95%
chromosome-bootstrap interval 0.01024 to 0.03876), 0.65959 versus 0.61365 in
Brassica (+0.04594; 0.03208 to 0.05888), and 0.49138 versus 0.48016 in apple
(+0.01122; 0.00265 to 0.01912). It retained 99.18%, 98.00% and 99.62% of the
respective v0.3 AP gains. Top-1% exact-positive yields were not reduced.

On apple, only 4/4,903 conflict sets and 8/30,511 candidates were guarded;
conflict top-1 returned from 163/269 under v0.3 to the baseline 165/269. The
same exact baseline-winner guarantee held on retrospective cotton and maize.
Cotton retained 96.09% and maize 93.82% of the v0.3 AP gain, although the
maize guard-minus-baseline interval still included zero. These diagnostics
cannot select the policy.

## Formal Populus v0.4 confirmatory result

The untouched Populus trichocarpa v4.1 evaluation is complete. Blind candidate
generation ran without network, `/nas_data`, complete target annotation,
evaluator references, or truth labels; all 294 blind artifacts were checksum-
verified before evaluator access. The formal evaluator obtained 800 events
across all 19 target chromosomes and all four CDS-complexity bins. No frozen
pair, reference, endpoint, threshold, or gate was relaxed.

H1 passed: the chain-preserving pool recovered 731/800 exact phased CDS chains
versus 663/800 for the legacy suppress-overlap pool, a gain of 0.0850 with a
20,000-replicate paired-bootstrap interval of 0.06625 to 0.1050. H2 was then
tested in the fixed sequence and passed: v0.4 AP was 0.85319 versus 0.84575 for
the frozen baseline, a gain of 0.00744 with a 20,000/20,000-valid chromosome-
bootstrap interval of 0.00421 to 0.01076.

All safety gates passed. Every one of 1,783 conflict-set winners matched the
baseline mapping, top-1% review yield was noninferior, topology-positive
coverage was 0.7510, 99.67% of the v0.3 AP gain was retained, all seven arms
had zero collateral transcript loss, and automatic approvals were zero. At
the fixed top-100 review budget, 98/100 candidates were exact positives; at
top 250 the yield was 246/250. These are benchmark review-priority results,
not calibrated probabilities or automatic annotation decisions.

## Publication interpretation

The evidence now supports a strong plant-genomics methods paper: the baseline
ranker has a compelling evidence-blind natural-genome enrichment, v0.3 fixes a
reproducible information-loss mechanism across four retrospective species and
apple, and both frozen v0.3 efficacy endpoints transferred to a fresh species.

The required new-species v0.4 confirmation has now also succeeded. On untouched
Populus, the chain-preserving pool recovered 731/800 exact structures versus
663/800 for the legacy pool (delta +0.0850; 95% interval +0.06625 to +0.1050).
The v0.4 ranker reached AP 0.85319 versus 0.84575 for baseline (delta +0.00744;
95% chromosome-bootstrap interval +0.00421 to +0.01076). All preregistered
evaluability, conflict-safety, review-yield, topology-coverage, gain-retention,
collateral-loss, bootstrap-validity and no-automatic-approval gates passed.

This resolved the largest algorithmic and leakage-control gap for a Genome
Biology-level submission. Subsequent work added untouched Actinidia and Coffea
replications, raw-read apple full-chain evidence, matched SynGAP comparison,
real-workload scaling, a concise reversible end-to-end workflow, independently
reproducible wheel/sdist builds, a non-root container, six audited figures and
the complete editable manuscript/submission package. The remaining release and
portal tasks require owner-supplied authorship, declarations, license,
repository URL and archive DOI; they are not unresolved scientific experiments.
