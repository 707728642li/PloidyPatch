# Coffea execution patch 2: single-reference descriptive scopes

The second isolated blind attempt completed all six raw projection jobs and
successfully constructed the preregistered combined-reference primary pools.
It then failed before blind custody and before truth-label access while
starting the descriptive `bua_only` scope. The generic candidate-GFF
namespacing primitive rejected an exact one-reference input with:

```text
ValueError: At least two reference GFF inputs are required for merging
```

This was an execution-interface restriction, not a biological or statistical
failure. The frozen protocol explicitly requires the combined two-reference
scope as primary and each one-reference scope as descriptive. Namespacing one
reference is well-defined, retains one method-family vote, and uses the same
identifier transformation and provenance schema as namespacing two or more
references.

Patch 2 changes only that minimum input cardinality from two to one and adds
tests for exact single-reference namespacing and empty-input rejection. It does
not change adapters, candidate coordinates, method support, overlap thresholds,
pool policies, target events, evaluator truth, H1, bootstrap rules, or safety
gates. No label or complete target annotation was visible to the blind runner.

Because the failed attempt had already generated primary combined pools, the
patch freeze records their candidate-GFF and decision-table SHA-256 values.
Final custody must require the rerun's four primary combined outputs to match
those values exactly. Thus the patch can enable the two descriptive
single-reference scopes without changing the already computed primary blind
candidate universe.
