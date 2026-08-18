# Apple external validation execution amendment v0.3.2

Date: 2026-08-08

Execution v0.3.1 successfully froze 5,291 truth pairs, 800 hidden events
and a 30,511-candidate chain-preserving pool. The blind candidate self-WGDI
then completed DIAMOND and WGDI but stopped before producing any inferred
blind pair, topology feature, candidate label or model score. Its source GFF3
contains both provider genes whose WGDI identifiers are stored in `Name` and
PloidyPatch candidate genes whose WGDI identifiers equal their exact `ID`.
The v0.3.1 adapter was intentionally limited to `Name`, so it was unsuitable
for this mixed source and the unchanged core mapper rejected the provider IDs.

Execution v0.3.2 generalizes only the evaluator-free identifier adapter. Each
WGDI representative ID must equal exactly one value among `ID`, `gene_id`,
`locus_tag` or `Name` on exactly one source `gene` record. The adapter adds
that same value as `gene_id` while retaining the original feature `ID`,
seqid, coordinates and strand. Missing, duplicated, conflicting or unsafe
aliases remain fatal. There is no prefix stripping, substring matching,
coordinate matching or biological selection.

The v0.3.1 truth pairs, hidden benchmark and candidate pool remain byte
unchanged and are checksum inputs to the amended execution freeze. No label
or model score existed when this amendment was made. No reference role,
sequence, coordinate, collinearity result, pair rule, event sampler, candidate
policy, model, numerical feature, endpoint, bootstrap setting or success gate
changed. Both failed working directories and prior execution freezes are
retained for audit.
