# PloidyPatch user guide

## What PloidyPatch is for

PloidyPatch is a review-oriented repair layer for plant annotations affected
by whole-genome duplication, homeology, and overlapping projection models. Its
core output is not a replacement annotation. It is a frozen candidate pool,
an evidence-linked review record, and an immutable patch containing only the
changes a reviewer explicitly accepted.

Use PloidyPatch when all of the following are true:

- a source GFF3 must remain recoverable byte for byte;
- candidate structures come from one or more declared upstream tools;
- distinct CDS chains should be reviewed instead of suppressed merely because
  their spans overlap;
- every candidate must retain its method and source-gene provenance; and
- annotation changes require an explicit human decision.

Do not use it as evidence that an unannotated gene is biologically absent, as
a general de novo gene predictor, or as a calibrated probability of gene
validity.

## Install the core

PloidyPatch supports Python 3.11 or newer. The core wheel has no mandatory
third-party Python dependency.

For the stable tagged release and optional dependency groups, see the
[installation guide](INSTALLATION.md).

```bash
conda create -p envs/ploidypatch -c conda-forge python=3.11 -y
conda run -p envs/ploidypatch --no-capture-output \
  python -m pip install .
conda run -p envs/ploidypatch --no-capture-output \
  ploidypatch --version
```

Upstream projection and synteny tools are intentionally separate. Install and
record them in workflow-specific environments rather than adding them to the
core prefix.

## Run the distributable example

```bash
python examples/minimal_reviewed_patch/run_example.py \
  --output-dir work/minimal_reviewed_patch
```

The input directory contains:

| File | Role |
| --- | --- |
| `source.gff3` | Original annotation; never modified in place. |
| `candidate.gff3` | Frozen candidate structures. |
| `decisions.tsv` | Candidate-pool and conflict decisions. |
| `candidate.gff3.manifest.json` | Candidate-universe and policy binding. |
| `review_decisions.tsv` | Explicit accept/reject ledger. |
| `report_scores.tsv` | Illustrative review-priority scores and automatic-approval audit field. |
| `report_copy_features.tsv` | Illustrative exact method-support intersections. |
| `report_topology.tsv` | Illustrative optional topology context. |
| `input_manifest.tsv` | Expected bytes and SHA-256 for every input. |

The runner invokes five public operations: compile reviewed additions, create
the patch, apply the patch, revert it, and generate the offline report. It then checks that the accepted
candidates are present, the rejected conflict is absent, automatic approval is
false, the report reaches `validated_reversible_run`, and the reverted GFF3 is
byte-identical to the source.

Outputs are first written to an adjacent `.working` directory and atomically
renamed only after all checks pass. Existing output paths are never
overwritten.

## Apply the reviewed-patch path to your data

### 1. Freeze the source and candidate pool

Record file sizes and SHA-256 values before review. Candidate identifiers must
be unique and stable. Mutually exclusive candidates must share a declared
conflict group; do not resolve a conflict by silently deleting alternatives.

### 2. Record review decisions

Use the exact five-column review schema demonstrated in
`examples/minimal_reviewed_patch/review_decisions.tsv`. Every candidate in the
frozen pool must have one explicit decision. A review ledger is data: checksum
it, retain it, and do not infer missing decisions as acceptance.

### 3. Compile accepted additions

```bash
ploidypatch patch compile-reviewed-copy-additions \
  --annotation-gff source.gff3 \
  --candidate-gff candidate.gff3 \
  --pool-decisions decisions.tsv \
  --pool-manifest candidate.gff3.manifest.json \
  --review-decisions review_decisions.tsv \
  --output-edits-json reviewed_additions.edits.json
```

The compiler verifies the pool lineage and conflict decisions and emits edits
only for explicit `accept_addition` records. It refuses ambiguous, incomplete,
or mismatched inputs.

### 4. Create and apply an immutable patch

```bash
ploidypatch patch create \
  --source-gff source.gff3 \
  --edits-json reviewed_additions.edits.json \
  --output-patch reviewed_additions.patch.json

ploidypatch patch apply \
  --source-gff source.gff3 \
  --patch reviewed_additions.patch.json \
  --output-gff annotation.with_reviewed_additions.gff3
```

The patch is content-bound to the source. Applying it to different bytes must
fail rather than guess.

### 5. Verify reversibility

```bash
ploidypatch patch revert \
  --patched-gff annotation.with_reviewed_additions.gff3 \
  --patch reviewed_additions.patch.json \
  --output-gff annotation.reverted.gff3

sha256sum source.gff3 annotation.reverted.gff3
```

The two hashes must be identical. Keep the source, patch, review ledger,
command log, and output checksum manifest together.

### 6. Generate the interactive result report

After the run is frozen, generate a report beside the run directory:

```bash
ploidypatch report \
  --candidate-gff candidate.gff3 \
  --pool-decisions decisions.tsv \
  --pool-manifest candidate.gff3.manifest.json \
  --review-decisions review_decisions.tsv \
  --patch-edits reviewed_additions.edits.json \
  --run-summary work/minimal_reviewed_patch/run_summary.json \
  --output-dir work/minimal_reviewed_patch_report
```

Open `index.html` locally. The report combines the current decision state,
ranked triage queues, exact method intersections, explicit safety
interpretation, searchable candidates, shared-axis conflict comparison,
phased CDS-chain diagrams, and checksum provenance. For a new review, edit the
generated `review_ledger_template.tsv` and supply it back through
`--review-decisions`. Use `--fail-on-attention` in automated pipelines and
`--max-embedded-candidates` to bound HTML size; JSON and TSV remain complete.
It remains a report, not a decision engine: scores only prioritize inspection,
topology is optional supporting evidence, and only explicit human decisions
can enter the patch. See the [complete result-report guide](RESULT_REPORT.md).

## Auditing an annotation bundle

```bash
ploidypatch audit \
  --gff annotation.gff3.gz \
  --protein proteins.fa.gz \
  --cds cds.fa.gz \
  --fai genome.fa.fai \
  --checksums \
  --output audit.json
```

An audit checks identifiers, parent-child structure, coordinates, phases,
sequence relationships, and optional coordinate bounds. Passing the audit
means the bundle is internally usable; it does not establish biological
completeness.

## Understanding ranking and topology

Scores prioritize review; they do not approve annotation changes. The frozen
conflict-winner guard can abstain from applying an additive correction, but it
cannot authorize a candidate. Duplication topology is useful in some species
and unavailable or unstable in others, so it is reported with an applicability
and coverage audit. The cross-species v0.9 replacement ranker was retired after
failing its frozen transfer gates.

For new species, keep the chain-preserving candidate workflow as the stable
core. Add topology only under a preregistered, label-free applicability rule,
and retain a topology-free score for every candidate where the evidence is not
available.

## Failure behavior

PloidyPatch follows fail-closed conventions:

- invalid checksum, schema, custody, or role isolation: `invalid_run`;
- insufficient preregistered events or coverage: `not_evaluable`;
- an evaluable run that misses a scientific gate: retained negative result;
- an existing output path: refusal rather than overwrite; and
- missing review decision: no automatic acceptance.

These states must not be merged. An invalid run is not a biological negative,
and a negative result must not be erased by changing thresholds after labels
are visible.

## Next references

- [Step-by-step tutorial](TUTORIAL.md)
- [Complete generated parameter reference](CLI_PARAMETERS.md)
- [CLI reference](CLI_REFERENCE.md)
- [Interactive result report](RESULT_REPORT.md)
- [Reproducibility guide](REPRODUCIBILITY_GUIDE.md)
- [Reviewed-copy patch contract](REVIEWED_COPY_PATCH_v0.1.md)
- [Species applicability framework](SPECIES_APPLICABILITY_FRAMEWORK_v0.6.md)
- [Submission readiness and claim boundary](SUBMISSION_READINESS_2026-08-10.md)
