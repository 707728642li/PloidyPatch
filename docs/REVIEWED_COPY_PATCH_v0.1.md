# Reviewed copy-addition patch workflow v0.1

PloidyPatch ranks copy-addition candidates for human review; it does not
automatically approve annotation changes. The public reviewed compiler binds
every accepted addition to the frozen chain-preserving candidate pool, its
checksum manifest, and an explicit reviewer decision.

## Review decision contract

Create a tab-separated file with exactly these five columns:

```text
candidate_digest	review_decision	reviewer	reviewed_at_utc	rationale
```

- `candidate_digest` is the 64-character digest from the frozen pool.
- `review_decision` is `accept_addition`, `reject`, or `defer`.
- `reviewer` is a non-empty person or review-team identifier.
- `reviewed_at_utc` is an ISO-8601 UTC timestamp, for example
  `2026-08-08T12:00:00Z`.
- `rationale` is a non-empty audit note. Avoid tabs and newlines inside it.

Unreviewed, rejected, and deferred candidates are excluded. At most one member
of a conflict set may be accepted. The compiler refuses duplicate or unknown
digests, non-UTC timestamps, truth/label columns in pool decisions, changed
input checksums, non-chain-preserving pool policies, existing outputs, and
stale `.working` files.

## Compile, apply, and prove reversibility

Use the original annotation, the complete candidate GFF, and the adjacent
frozen pool artifacts. Do not prefilter or rewrite these three inputs.

```bash
ploidypatch patch compile-reviewed-copy-additions \
  --annotation-gff primary_chromosomes.gff3 \
  --candidate-gff candidate.gff3 \
  --pool-decisions decisions.tsv \
  --pool-manifest candidate.gff3.manifest.json \
  --review-decisions review_decisions.tsv \
  --output-edits-json reviewed_copy_additions.edits.json

ploidypatch patch create \
  --source-gff primary_chromosomes.gff3 \
  --edits-json reviewed_copy_additions.edits.json \
  --output-patch reviewed_copy_additions.patch.json

ploidypatch patch apply \
  --source-gff primary_chromosomes.gff3 \
  --patch reviewed_copy_additions.patch.json \
  --output-gff annotation.with_reviewed_additions.gff3

ploidypatch patch revert \
  --patched-gff annotation.with_reviewed_additions.gff3 \
  --patch reviewed_copy_additions.patch.json \
  --output-gff annotation.reverted.gff3

sha256sum primary_chromosomes.gff3 annotation.reverted.gff3
```

The two final hashes must be identical. The edits artifact embeds the input
hashes, review provenance for every appended hierarchy, decision counts,
`automatic_approval=false`, zero existing-feature deletions, and zero edits to
existing annotation lines.

`compile-copy-additions` remains available as a low-level compatibility
primitive for an already-prefiltered candidate GFF. New reviewed workflows
should use `compile-reviewed-copy-additions` so that selection provenance and
the frozen candidate universe are verified by the CLI.
