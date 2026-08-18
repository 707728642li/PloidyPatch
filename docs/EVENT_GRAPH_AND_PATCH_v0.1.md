# Typed plant event graph and reversible patch contract v0.1

Date frozen: 2026-08-07  
Status: transparent inference MVP; uncalibrated

## Evidence graph

`ploidypatch graph infer` consumes an event-candidate TSV and a typed evidence
TSV. One biological event may have multiple mutually exclusive candidates.
Every candidate records a locus state, source gene IDs, optional relationship
state, named WGD event, subgenome, location, and proposal URI. Every evidence
edge records:

- a globally unique evidence ID and target candidate;
- evidence type and whether it supports, contradicts, or only contextualizes
  that candidate;
- coding, isoform, both, or relationship scope;
- strength and reliability in `[0,1]`;
- an independent-group label and exact source/provenance string; and
- optional details retained verbatim.

Direction is always relative to the candidate hypothesis. For example, an
intact protein projection is `support` for `missing_annotation` but
`contradict` for `true_absence`. Contradictory edges are retained in the JSON
graph and never discarded because another edge is stronger.

The v1 evidence vocabulary covers protein/CDS/domain projection, synteny,
named-WGD quota, subgenome, gene tree, RNA junction/coverage, full-length
transcripts, assembly gaps, depth/copy number, mappability, repeats, plant
splice scores, pseudogene lesions, whole-genome-alignment absence, PAV
deletions, and haplotype alignment. Unknown evidence types fail closed.

## Correlation guard and dual confidence

The MVP uses declared fixed weights and produces separate coding-model,
isoform-model, and relationship logits. Within one independent group, repeated
supporting edges contribute only the strongest support; repeated contradictory
edges contribute only the strongest contradiction. The two are then added, so
conflict remains visible but correlated tools cannot manufacture confidence by
simple vote counting.

Logistic transforms are exposed as `confidence`, but the graph labels them
`transparent_uncalibrated_v1`. They are ranking scores, not yet empirical
probabilities. The 0.8 accept and 0.2 reject cutoffs are engineering defaults
that must later be fitted or calibrated only on development clades, followed
by held-out reliability, Brier, and ECE evaluation. Missing evidence yields
`null`, not 0.5, and forces review when that scope is required.

Coding and isoform confidence are deliberately distinct. Structural repair
states require both; missing annotation requires coding support. This encodes
the soybean result that an exact CDS chain can be credible while complete
transcript/UTR structure remains weak.

## Plant constraints

Hard constraints cause `abstain` regardless of the aggregate score:

- `true_absence` conflicts with strong intact protein/full-length-transcript
  evidence and with assembly-gap uncertainty, and requires strong sequence-
  level WGA/PAV absence evidence;
- `wgd_homeolog` requires a named WGD event and strong synteny support;
- `collapsed_copies` requires independent copy-number or read-depth support;
- `haplotig_duplicate` requires strong haplotype-alignment evidence; and
- split, fusion, exon, boundary, and missing-annotation repairs require strong
  protein, CDS-chain, RNA-junction, or full-length-transcript structure support;
- a multi-locus projection cannot select a repair without independent synteny
  evidence that disambiguates the intended plant duplicate/homeolog locus.

Candidates sharing an `event_id` are mutually exclusive. If two independently
acceptable alternatives differ by less than 0.1 mean required confidence,
both revert to review rather than being broken by an arbitrary tie. With a
larger margin, only the best is selected and the other remains reviewable.
This is the first constraint layer, not yet a full genome-wide optimizer.

## Input columns

Candidate TSV required columns:

```text
event_id  event_type  candidate_id  gene_ids
```

Optional columns are `seqid`, `start`, `end`, `strand`,
`relationship_state`, `wgd_event`, `subgenome`, and `proposal_uri`.

Evidence TSV required columns:

```text
candidate_id  evidence_id  evidence_type  direction  scope  strength
reliability  independent_group  source
```

`details` is optional. Both input hashes, all weights, group-level
contributions, constraints, decisions, and edges are frozen in
`ploidypatch.event_graph.v1` JSON. A compact decision TSV is emitted alongside
it.

```bash
ploidypatch graph infer \
  --candidates event_candidates.tsv \
  --evidence evidence_edges.tsv \
  --output-json event_graph.json \
  --decisions-tsv decisions.tsv
```

## Immutable annotation patches

`ploidypatch.annotation_patch.v1` represents each edit as one source line plus
zero, one, or multiple replacement lines. It therefore supports deletion,
coordinate/attribute replacement, fusion, and split without an in-place write.
The patch stores original and replacement bytes and hashes, complete source and
patched text hashes, line counts, edit-spec provenance, and event IDs.

An edit specification may contain direct `operations`, or `events` with the
same `line_edits` emitted by the controlled perturbation framework:

```json
{
  "event_ids": ["event-17"],
  "operations": [
    {
      "source_line_number": 42,
      "replacement_lines": ["chr1\tPloidyPatch\tgene\t...\n"]
    }
  ]
}
```

```bash
ploidypatch patch create \
  --source-gff original.gff3 --edits-json accepted_edits.json \
  --output-patch accepted.patch.json

ploidypatch patch apply \
  --source-gff original.gff3 --patch accepted.patch.json \
  --output-gff patched.gff3

ploidypatch patch revert \
  --patched-gff patched.gff3 --patch accepted.patch.json \
  --output-gff restored.gff3
```

Apply requires exact source bytes; revert requires exact patched bytes and
reconstructs the exact original decompressed GFF3 checksum. All commands refuse
to overwrite outputs. Patches are a transport and audit primitive: only
reviewed or sufficiently validated graph candidates should be converted into
an accepted patch.
