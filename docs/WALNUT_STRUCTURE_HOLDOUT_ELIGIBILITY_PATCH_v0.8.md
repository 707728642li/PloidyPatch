# Walnut structure-holdout eligibility patch v0.8

The first Walnut structure-holdout construction attempt was executed after
the evaluator-only WGD pair tree had been frozen but before candidate
generation, blind execution, scoring, or label access. It stopped at the
pre-candidate no-op sentinel because at least one deterministically collapsed
phased-CDS chain was still represented by another source transcript. The
benchmark was therefore structurally ambiguous and no blind benchmark was
published.

This execution patch preserves that failure and keeps the sentinel unchanged.
Before balanced sampling, it requires the phased-CDS chain of the
deterministically selected collapsed partner to occur exactly once in the
complete source annotation. Non-unique chains receive the explicit rejection
reason `collapsed_phased_cds_chain_not_unique_in_source`. The rule uses only
the source annotation and the already frozen pair-selection seed; it does not
read candidates, scores, performance labels, or review outcomes.

No reference role, WGD pair rule, Ks interval, event target, balance rule,
minimum event gate, bootstrap, hypothesis, candidate method, overlap policy,
or success threshold is changed. If fewer than the preregistered 500 events or
the chromosome/complexity minima remain eligible, Walnut is formally
`not_evaluable`; the eligibility rule and sentinels must not be relaxed.
