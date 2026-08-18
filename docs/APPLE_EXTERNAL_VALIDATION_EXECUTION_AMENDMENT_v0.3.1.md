# Apple external validation execution amendment v0.3.1

Date: 2026-08-08

The first v0.3 execution stopped before producing any self-WGDI pair,
outgroup pair, hidden event, candidate label or apple model score. The frozen
WGDI preparation uses unprefixed apple identifiers such as `MD01G1000100`,
whereas the provider GFF3 stores the corresponding feature ID as
`gene:MD01G1000100` and its exact `Name` as `MD01G1000100`. The existing
pair mapper deliberately accepts only exact `ID`, `gene_id` or `locus_tag`
attributes, so it rejected all retained apple WGDI identifiers rather than
guessing that the prefix should be removed.

Execution v0.3.1 adds an evaluator-only, gene-only GFF3 alias map. For every
gene in the already frozen WGDI representative table, the adapter requires
exactly one source `gene` feature whose `Name` equals the WGDI gene ID. It
retains the original feature `ID`, seqid, coordinates and strand and adds the
same exact value as a `gene_id` alias. Missing, duplicated, unsafe or
conflicting aliases are fatal. The complete frozen source GFF3 remains the
benchmark source; the derivative is used only to translate WGDI identifiers
back to those original feature IDs.

No reference role, sequence, coordinate, gene selection, collinearity result,
pair rule, event sampler, model, feature, score, endpoint, bootstrap setting or
success gate changed. The failed working directory and original execution
freeze are retained as invalid audit artifacts. A new execution freeze must be
created before pair inference is retried.
