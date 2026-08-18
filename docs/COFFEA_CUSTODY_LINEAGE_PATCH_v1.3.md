# Coffea custody lineage patch v1.3

The formal stable blind run completed under execution patch 3 and its custody
manifest correctly binds that frozen execution checksum. Reveal patch 4 was
created only after blind custody, so requiring the pre-existing custody to bind
the newer reveal-only execution freeze caused a fail-closed stop before a new
reveal authorization, evaluator input access, or metric computation.

This patch records the existing custody execution checksum explicitly in the
new execution manifest and verifies both that checksum and the immutable
custody-manifest checksum during reveal. It does not regenerate blind outputs
or change candidates, references, truth, perturbations, H1 endpoints,
bootstrap seeds, thresholds, or statistical rules. The next evaluator result
is retained in its resulting direction.
