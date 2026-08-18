# Chain-preserving candidate pool v0.3 (2026-08-08)

## Design

The v0.3 pool deduplicates identical `seqid + strand + phased CDS chain`
records across method families but does not discard a structurally distinct
chain merely because more than half of its CDS overlaps another candidate.
Strongly overlapping alternatives are retained in a stable connected conflict
set and marked mutually exclusive for ranking and review. The v1 suppressing
policy remains the default for backward compatibility; no frozen v0.1/v0.2
artifact was changed.

This is an information-preserving candidate-universe change, not an automatic
repair policy. Its direct endpoint is exact-chain recall; precision and F1 of
the unranked pool are expected to decline.

## Paired benchmark result

All comparisons use the same frozen individual-method candidates, paired
complete-annotation background subtraction, and 800 exact copy-collapse
events per species. Soybean and Brassica are development data; cotton and
maize are retrospective diagnostics whose labels were already visible and
cannot serve as new v0.3 holdouts.

| Species | Old recovered | v0.3 recovered | Recall gain | 95% paired event-bootstrap CI | v0.3 candidates | Conflict sets |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Soybean | 619 | 630 | +0.0138 | 0.0062 to 0.0225 | 20,504 | 3,182 |
| Brassica napus | 607 | 653 | +0.0575 | 0.0425 to 0.0737 | 76,434 | 16,011 |
| Cotton | 553 | 639 | +0.1075 | 0.0863 to 0.1300 | 44,493 | 6,509 |
| Maize | 566 | 653 | +0.1088 | 0.0875 to 0.1312 | 18,606 | 2,448 |

The recall gains are supported in all four species. In maize, v0.3 reaches the
653/800 oracle ceiling already present among the three individual method
outputs, showing that all 87 structures lost by the old cross-method overlap
suppression are restored.

## Precision boundary

Unranked v0.3 F1 is lower than the old suppressing union: 0.7697 versus 0.7972
in soybean, 0.7010 versus 0.7597 in Brassica, 0.6530 versus 0.6900 in cotton,
and 0.6399 versus 0.7137 in maize. This is not a regression in the declared
candidate-pool objective; it demonstrates why conflict-aware ranking is
mandatory before review. The v0.3 pool must never be interpreted as a set of
approved edits.

## Next confirmatory step

Ranking development is restricted to soybean and Brassica. The next model
should rank alternatives within conflict sets, use method-support regime as an
explicit interaction, and abstain at unresolved conflicts. Cotton and maize
may be used for transparent retrospective diagnostics only. A new species and
frozen protocol are required for a confirmatory v0.3 portability claim.
