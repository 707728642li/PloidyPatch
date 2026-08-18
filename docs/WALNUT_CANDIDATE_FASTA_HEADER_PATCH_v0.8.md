# Walnut candidate FASTA header compatibility patch v0.8

The third blind attempt showed that GeMoMa treated the complete descriptive
FASTA header of the *Juglans mandshurica* reference as its sequence identifier,
although the first token and every retained GFF seqid were already identical.
It consequently reported all 16 reference chromosomes as missing and extracted
zero gene models.

This execution-only patch writes candidate-reference primary FASTA headers as
the exact retained seqid. It does not alter sequence bytes, sequence order,
GFF records, coordinates, proteins, candidate policies, pools, truth, labels,
endpoints, or statistical gates. The default normalization behavior remains
unchanged; exact-seqid headers are enabled only by the Walnut blind candidate
adapter and are recorded in its bundle manifest.

The same failed attempt exposed an independent shell-initialization error in
the Miniprot arm. Its output path and working path were declared together, so
`set -u` expanded the not-yet-assigned output variable and stopped both
Miniprot jobs before candidate generation. The declarations are now ordered on
separate lines. No Miniprot option or biological acceptance rule changed.

Because this failure occurred after candidate methods had started but before
any candidate pool or other formal blind output existed, the execution freezer
now records the exact failed stage. It no longer labels every chained patch as
a structure-holdout failure. The freeze continues to require incomplete
candidate generation, no formal scores, and no truth-label access.

The next full blind attempt completed all six projection arms but stopped
before pool construction because the Walnut wrapper expected a Miniprot-style
`SHA256SUMS` file from every publisher. The shared GeMoMa and LiftOn publishers
instead emit `output_manifest.tsv`. A dedicated verifier now requires the
exact expected role set and rechecks every published file's size, SHA-256, and
pre-rename `.working` path lineage. No raw prediction or pool rule changed.
The pool builder calls the same shared verifier rather than assuming all six
publishers have a `SHA256SUMS` file; this closes the same interface mismatch at
both consumers without changing any adapted prediction or consensus policy.
Pool-construction failures are now retained as an immutable `.invalid_run`
working tree instead of deleting the exact command logs needed for diagnosis.
The retained smoke failure showed that reference-specific adapted GFFs were
being merged with two copies of the exact base prefix. The builder now verifies
and removes each byte-identical base prefix, namespaces and merges only the
candidate suffixes, and restores that same base exactly once. Its result root
also now matches the frozen `external/walnut_v0.8_h1` contract literally.

## Candidate-output publication amendment

The sixth isolated blind attempt completed all six raw method arms and both
candidate-pool policies without access to evaluator data, truth labels, the
complete target annotation, `/nas_data`, a model, or the network.  It then
failed with exit status 72 before custody because the pool builder published
under `external/walnut/v0.8_h1`, while the protocol and custody contract require
`external/walnut_v0.8_h1`.  The failed tree and its isolation evidence are
checksum-bound in the superseding execution freeze.  No artifact from that
attempt may be reused.  The only runtime correction is the literal output-root
constant above; candidate methods, adapters, pool policies, thresholds,
references, truth construction, evaluator, endpoints, and gates are unchanged.
An additional smoke gate found exact phased-CDS chains projected by both
references within one method family. These are collapsed by exact chain only,
with a deterministic namespaced-transcript tie-break and an explicit audit;
distinct chains remain distinct and each method family still contributes at
most one vote.
The final seal preserves each arm's native frozen schema: chain-preserving v2
for retain-distinct and consensus v1 for suppress-overlap. The pool sealer and
blind custody validator consume one shared arm-to-schema mapping, so the legacy
comparator is never mislabeled as a chain-preserving pool.
Execution-patch provenance now distinguishes a completed six-method stage from
a subsequent candidate-pool failure, rather than misreporting both as a method
failure.
