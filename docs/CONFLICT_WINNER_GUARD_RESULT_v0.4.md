# Conflict winner guard development result v0.4

Date: 2026-08-08

## Decision

The label-independent v0.4 conflict winner guard passed every predeclared
development gate. It is suitable for a versioned freeze followed by testing
on a new untouched plant species. It is not externally confirmed: apple
labels were already visible when v0.4 was designed, and cotton and maize are
retrospective diagnostics.

The guard retains the v0.3 score unless v0.3 changes the deterministic top
candidate within a conflict set relative to the frozen non-topology baseline.
In such a set, all members fall back to their baseline logits and are marked
as topology-abstained/manual review. Ties are resolved by descending score and
then candidate digest. No candidate is automatically approved.

## Fixed development gates

For soybean, Brassica and label-seen apple, all of the following were required:

- guard-minus-baseline AP greater than zero with a positive lower bound from
  20,000 chromosome-bootstrap replicates;
- retention of at least 90% of the v0.3 AP gain;
- top-1% exact-positive review yield no lower than baseline;
- every conflict-set top-1 candidate exactly identical to baseline;
- zero automatic approvals.

All gates passed in all three development datasets. Cotton and maize were
reported but could not select the policy.

## Ranking and safety results

| Dataset | Role | Baseline AP | v0.3 AP | v0.4 AP | v0.4 minus baseline (95% CI) | v0.3 gain retained | Conflict top-1 baseline / v0.3 / v0.4 |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Soybean | development | 0.85344 | 0.87657 | 0.87638 | +0.02294 (0.01024 to 0.03876) | 99.18% | 33 / 35 / 33 of 43 |
| Brassica | development | 0.61365 | 0.66053 | 0.65959 | +0.04594 (0.03208 to 0.05888) | 98.00% | 117 / 124 / 117 of 153 |
| Apple | seen external development | 0.48016 | 0.49142 | 0.49138 | +0.01122 (0.00265 to 0.01912) | 99.62% | 165 / 163 / 165 of 269 |
| Cotton | retrospective diagnostic | 0.70991 | 0.74753 | 0.74606 | +0.03615 (0.02837 to 0.04388) | 96.09% | 138 / 147 / 138 of 211 |
| Maize | retrospective diagnostic | 0.59499 | 0.59765 | 0.59749 | +0.00250 (-0.00333 to 0.00760) | 93.82% | 174 / 175 / 174 of 273 |

The top-1% exact-positive counts for baseline, v0.3 and v0.4 were respectively
187/193/193 in soybean, 461/483/484 in Brassica, and 188/191/191 in apple.
Thus the safety correction did not reduce the fixed development review yield.

The guard is deliberately sparse. It changed 8 conflict sets and 19 candidates
in soybean, 81 sets and 193 candidates in Brassica, 4 sets and 8 candidates in
apple, 91 sets and 209 candidates in cotton, and 20 sets and 47 candidates in
maize. It restored the baseline winner in every conflict set by construction
while leaving all other v0.3 scores unchanged.

The small AP cost relative to v0.3 was statistically distinguishable in
Brassica (-0.00094; interval -0.00193 to -0.00017) and cotton (-0.00147;
interval -0.00219 to -0.00080), but both retained at least 96% of the larger
v0.3 gain and remained clearly above baseline. The apple, soybean and maize
v0.4-minus-v0.3 intervals included zero.

## Audit trail

The first checksum list accidentally included its own initially empty file.
Both substantive artifacts passed their listed hashes; the invalid
self-referential list was retained as `SHA256SUMS.invalid_self_reference`, the
generator was corrected to exclude the checksum file, and a valid checksum
list was regenerated. No prediction or evaluation content changed.

Server result:

`/data/codexli/projects/PloidyPatch/results/copy_collapse/model_development/conflict_winner_guard_v0.4_evaluation`

Final artifact hashes:

- `evaluation.json`: `cf726cb6cae6c4bdd8836e959fbbf3ba0bd7509e73ab1c99462ede7ecedb6fa2`
- `predictions.tsv`: `4206624d8923e9891b1ce0951dc61e8077d9b71321efb1607caa90d6d7f4454a`
- corrected `SHA256SUMS`: `161057e6c695db5e144d0629bc0fbc91d728bcb224afa885c08e8ae05eb6ae25`

## Interpretation

This result resolves the observed apple v0.3 safety failure at the algorithmic
level without erasing the cross-species ranking benefit. The defensible next
step is to freeze code, model and policy before labels on a new species are
generated, then require efficacy, review-yield and conflict-safety gates on
that untouched species. Until that test and the natural-locus biological
audits are complete, v0.4 remains promising development evidence rather than
confirmatory evidence for automatic annotation repair.
