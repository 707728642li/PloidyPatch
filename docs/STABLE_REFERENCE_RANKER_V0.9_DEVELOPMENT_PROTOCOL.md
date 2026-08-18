# Candidate-independent reference ranker v0.9 development protocol

Date frozen: 2026-08-09

## Purpose and boundary

The v0.5 ranker has useful review signal, but both its copy-level WGD fields
and its explicit homeolog-topology fields were derived from a self-synteny
graph containing the entire candidate pool. An unrelated candidate could
therefore change an incumbent candidate's features. Target-only fixed-backbone
and bounded local projections removed that dependence but failed the frozen
Actinidia coverage and ranking gates. Directly applying the old model to
candidate-reference-anchored features also caused a semantic distribution
shift.

Version 0.9 tests one correction only: refit the existing simple ranker on a
fully candidate-independent feature meaning. It is a new development model,
not a reinterpretation of any v0.5-v0.8 formal or development result.

## Frozen feature meaning

All method-quality, method-support, span, and CDS features retain their v1
definitions. The three copy-level WGD fields and all homeolog-topology fields
must derive from the same candidate-reference-anchored selection:

1. build each declared candidate reference's primary-chromosome self-WGD
   graph before target candidate access;
2. retain only the predeclared lineage event and reciprocal-unique homeolog
   pairs;
3. trace every target candidate to exact upstream reference genes;
4. map each source homeolog through an already frozen full-CDS target
   projection to exactly one existing target gene;
5. require at least one evidence arm, accept missing arms, and abstain when any
   observed arm names a different partner.

Candidate-candidate edges, target candidate-union self-synteny, peer
inheritance, label-derived thresholds, fuzzy identifier repair, and winner
selection among multiple partners are forbidden. The feature builder must
record `truth_access=false`, exact input hashes, candidate-universe hash, and
per-reason abstention counts.

## Development species and leakage boundary

Actinidia and Populus are label-seen development species. Their stable feature
trees must be checksum-frozen before labels are joined for this evaluation.
Neither species can provide a new external claim. Walnut is the metadata-
selected untouched target and remains unavailable to model selection,
training, thresholds, feature changes, and reference changes.

## Fixed estimators

The comparison uses the existing deterministic feature transforms and only
two rank-only estimators:

1. `stable_copy_baseline`: L2 logistic regression on the full copy-feature
   contract with candidate-reference WGD semantics, `C=1`, `lbfgs`, maximum
   5,000 iterations, random state 0, no calibration or threshold.
2. `stable_reference_topology`: the stable copy baseline logit plus the v0.3
   support-conditioned topology correction, with the same seven exact method-
   support patterns, L2 penalty and `C=1`.

No feature search, hyperparameter grid, probability calibration, automatic
approval, or alternative estimator may be selected after results are visible.
For a pooled final artifact, every species contributes equal total sample
weight; no class balancing is applied. Feature imputation and scaling must use
the same species-balanced weights. Implementation tests must reproduce the
unweighted contract when all weights are equal.

The exact pooled row weight is `N_total / (2 * N_species)`: each of the two
species therefore contributes `N_total / 2` total mass and mean row weight is
one. This preserves the declared `C=1` regularization scale instead of
silently strengthening it by normalizing the complete training set to mass
two.

## Evaluation

The frozen evaluation contains:

- five-seed, five-fold target-chromosome-grouped out-of-fold predictions in
  each species, with a hard assertion that repeated partitions differ;
- Actinidia-to-Populus and Populus-to-Actinidia zero-retuning transfer;
- 20,000-replicate target-chromosome bootstrap intervals for the AP delta;
- the six fixed review budgets: 0.5%, 1%, 2%, 100, 250, and 500;
- conflict winner mapping hashes, automatic-approval count, topology coverage,
  availability by label, and abstention reasons.

The five OOF seeds are `20260901, 20260911, 20260921, 20260931, 20260941`.
The 20,000-replicate Actinidia and Populus chromosome bootstraps use seeds
`20260951` and `20260952`, respectively. The production conflict-winner guard
is applied after averaging the five OOF logits and independently to each
cross-species transfer. Review ordering is descending score followed by the
candidate digest.

## Retention gates

The stable topology estimator may be pooled and frozen for Walnut only if all
of the following hold:

1. OOF AP exceeds the stable copy baseline in both Actinidia and Populus, with
   both 20,000-replicate bootstrap lower bounds above zero;
2. both cross-species transfer AP deltas are positive;
3. top-1% true positives do not decrease in either OOF species or transfer
   direction;
4. conflict-winner mismatch and automatic approvals are zero after the frozen
   guard;
5. every input, feature, prediction, partition, model, and evaluation artifact
   passes its exact manifest.

Failure retires v0.9 as a Walnut ranker and leaves the chain-preserving
candidate workflow and its H1 test intact. It must not trigger threshold
relaxation, a different estimator, a different Walnut reference panel, or a
same-version second attempt after Walnut labels are opened.
