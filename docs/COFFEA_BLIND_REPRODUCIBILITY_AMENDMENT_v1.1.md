# Coffea blind reproducibility amendment v1.1

The first complete replay of the Coffea blind candidate workflow stopped at
the frozen custody gate because the primary chain-preserving pool differed
from the preceding blind execution. No target labels, evaluator truth, or
complete target annotation had been opened.

The discrepancy was localized before reveal to nondeterministic LiftOn choices
among duplicated placements. Miniprot and GeMoMa adapter decisions were
byte-identical. LiftOn had seven source-model presence/identity changes and,
under the more conservative complete-row comparison used by this amendment,
19 non-identical decision rows across both references. The biological
references, perturbation,
candidate methods, overlap thresholds, consensus policy, H1 endpoint,
bootstrap seed, evaluability gates, and zero-collateral requirement are not
changed.

For the final Coffea execution, every method/reference adapter model must have
the same model ID and the exact complete decision row in both independent
blind runs. Models absent from either run or differing in any decision field
abstain. Stable accepted models are then passed through the frozen
within-method exact-chain collapse and the frozen retain-distinct and
suppress-overlap pool builders. The rule is applied to all methods and both
candidate references without inspecting labels. It is algorithmic and may not
contain a species gene allowlist or denylist.

This amendment is an execution-reproducibility safeguard, not performance
tuning. Both failed attempts, their logs, input hashes, and model-level
stability audit must remain immutable and be checksum-bound before reveal. The
Coffea result will be reported in its resulting direction; failure or
insufficient evaluability will not trigger another species, threshold change,
or ranker development cycle.
