# Conflict-aware ranker v0.3 development protocol

Date frozen: 2026-08-08  
Candidate contract: `ploidypatch.method_candidate_pool.v2`

## Hypothesis

The v0.2 failure analysis indicates that refitting every baseline coefficient
can cancel a useful topology add-on and that topology effects differ by method-
support regime. The v0.3 primary estimator therefore keeps a fold-specific
baseline logistic score as a fixed offset and learns only a regularized
topology correction. The correction contains:

- the five declared coding-topology components;
- a topology-availability intercept within each of seven exact method-support
  patterns; and
- mean topology coherence within each support pattern.

When topology is unavailable, every correction feature is zero and the score
is exactly the baseline score. Strongly overlapping alternative CDS chains
remain in explicit conflict sets; the ranker may prioritize but never approve
one automatically.

## Fixed estimators

1. `baseline`: L2 logistic regression on the existing full copy-feature
   contract, `C=1`, `lbfgs`, random state 0.
2. `global_refit`: the v0.2-style joint refit of baseline plus topology terms,
   retained as a mechanistic comparator.
3. `offset_global`: baseline logit plus one L2-regularized global topology
   correction.
4. `offset_support_conditioned` (primary): baseline logit plus the fixed
   support-conditioned correction above.

No hyperparameter grid, feature search, probability calibration, portable
threshold, or automatic approval is allowed. Offset corrections minimize
summed logistic loss plus `0.5 * ||beta||^2` (`C=1`) with analytic gradients.

## Data boundary

Soybean and Brassica are development species. Candidate features are frozen
before evaluator labels. Model assessment uses five-seed, five-fold
target-chromosome-grouped out-of-fold predictions plus both cross-species
transfer directions.

Cotton and maize labels have already been viewed. They are retrospective
diagnostics only and may not choose the estimator, features, regularization,
or success rule. A different preregistered species is required for the formal
v0.3 external test.

## Primary development criteria

The primary estimator is retained for external freezing only if:

1. AP is higher than baseline in chromosome-grouped OOF evaluation in both
   soybean and Brassica;
2. its pooled development AP gain is positive; and
3. it does not reduce top-1 selection accuracy among evaluable conflict sets
   containing exactly one true CDS chain in either development species.

Twenty-thousand-replicate chromosome-grouped bootstrap intervals are reported
for AP deltas. Intervals and the two cross-species directions determine the
strength of evidence but are not used to silently substitute another estimator
after results are visible.

## Reporting boundary

All four estimators, every seed, candidate counts, positive counts, topology
coverage, conflict-set metrics, review budgets, coefficients, checksums, and
software versions must be published. Failure leaves the chain-preserving pool
valid but blocks freezing a v0.3 ranking claim.
