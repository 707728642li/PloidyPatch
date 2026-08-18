# PloidyPatch result report

`ploidypatch report` converts a checksum-bound candidate-pool run into an
offline review package for project leads, annotators, and methods reviewers.
It summarizes the run; it does not make annotation decisions. A high score
means **inspect earlier**, never **accept automatically**.

## Report package

Every successful run creates a new, non-overwriting directory:

| File | Purpose |
|---|---|
| `index.html` | Self-contained interactive report; no server or network is required |
| `report.json` | Full machine-readable candidate universe, triage queues, conflicts, policy, and provenance |
| `candidates.tsv` | Full flat export for R, Python, or a spreadsheet |
| `review_ledger_template.tsv` | Five-column ledger ready for human decisions or resuming an existing review |
| `SHA256SUMS` | Exact-universe hashes for the other four files |

The HTML provides:

- a run-state banner and next safe action;
- exact candidate, conflict, review, patch, and automatic-approval audit state;
- ranked triage queues for unreviewed, deferred, and unresolved conflict sets;
- exact Miniprot/GeMoMa/LiftOn intersection counts;
- all-sequence genomic distribution and keyboard-readable score histograms;
- a searchable, keyboard-operable candidate table with readable rank position,
  percentile band, and raw score;
- phased CDS-chain diagrams with strand direction;
- a shared-axis, side-by-side conflict workbench with boundary, phase, and
  overlap differences;
- topology availability, abstention reason, partner, and coherence; and
- input file names, byte counts, and SHA-256 provenance.

## Verified example

Run the complete reviewed-patch example from the repository root:

```bash
python examples/minimal_reviewed_patch/run_example.py \
  --output-dir work/minimal-reviewed-patch
```

The runner verifies all eight inputs, compiles and applies the reviewed patch,
proves byte-identical reversion, and writes the interactive report to
`work/minimal-reviewed-patch/report/index.html`. The expected report state is
`validated_reversible_run`.

The equivalent explicit report command is:

```bash
ploidypatch report \
  --candidate-gff examples/minimal_reviewed_patch/candidate.gff3 \
  --pool-decisions examples/minimal_reviewed_patch/decisions.tsv \
  --pool-manifest examples/minimal_reviewed_patch/candidate.gff3.manifest.json \
  --review-decisions examples/minimal_reviewed_patch/review_decisions.tsv \
  --scores examples/minimal_reviewed_patch/report_scores.tsv \
  --copy-features examples/minimal_reviewed_patch/report_copy_features.tsv \
  --topology-features examples/minimal_reviewed_patch/report_topology.tsv \
  --patch-edits work/minimal-reviewed-patch/reviewed_copy_additions.edits.json \
  --run-summary work/minimal-reviewed-patch/run_summary.json \
  --title "Minimal reviewed patch" \
  --output-dir work/minimal-reviewed-patch-report \
  --fail-on-attention
```

## Required inputs

`--candidate-gff` is the chain-preserving candidate GFF3. Candidate genes must
carry unique `consensus_digest` values. Transcript and CDS features drive the
structure diagrams.

`--pool-decisions` must contain `consensus_digest`, `status`,
`conflict_set_digest`, and `conflict_member_count`.

`--pool-manifest` must use schema `ploidypatch.method_candidate_pool.v2`. The
report verifies the candidate GFF3 and decisions hashes, accepted count, pool
parameters, and exact candidate universe before creating output.

`--output-dir` must not exist. A same-name `.working` directory is also refused.

## Optional inputs and joins

| Argument | Adds | Join requirement |
|---|---|---|
| `--review-decisions` | Completion state, reviewers, decisions, rationales | Strict subset or full candidate universe |
| `--scores` | Rank, percentile, score histogram, automatic-approval audit | Exact candidate universe |
| `--copy-features` | Exact method intersections | Exact candidate universe |
| `--topology-features` | Applicability, reason, partner, coherence | Exact candidate universe |
| `--patch-edits` | Patch-event count and review/patch consistency | Event digests must be a pool subset |
| `--run-summary` | Byte-identical reversion evidence | Must bind accepted count/digests, hashes, and disabled automatic approval |
| `--title` | Display title | One non-empty line, at most 160 characters |
| `--max-embedded-candidates` | HTML payload limit; default 5,000 | Full JSON/TSV/ledger outputs are never truncated |
| `--fail-on-attention` | Exit status 2 for `attention_required` | Useful in CI and pipelines |

Candidate tables containing evaluator-style fields such as `truth`, `label`,
`ground_truth`, `oracle_label`, or `y_true` are rejected. Generic scientific
field names are not blocked solely by substring heuristics.

## Review-ledger workflow

1. Generate the report without `--review-decisions`.
2. Copy or edit `review_ledger_template.tsv`.
3. Set each `review_decision` to `accept_addition`, `reject`, or `defer`; record
   reviewer, UTC timestamp, and rationale.
4. Regenerate the report with `--review-decisions` pointing to that ledger.
5. Resolve every conflict set so no more than one alternative is accepted.
6. Compile the immutable patch, apply it to a new file, revert it, and require
   exact source-byte recovery.

The report preserves existing ledger rows and emits blank rows for candidates
not yet reviewed. The triage queue changes inspection order only; it never
creates an acceptance recommendation.

## Status interpretation

| State | Meaning | Next action |
|---|---|---|
| `review_in_progress` | One or more candidates lack a decision | Continue human review |
| `review_complete_patch_pending` | All candidates are decided; no matching patch supplied | Compile the reviewed patch |
| `ready_to_apply` | Accepted reviews exactly match patch events | Apply and prove reversion |
| `validated_reversible_run` | Patch and byte-identical reversion are verified | Archive the checksum-bound run |
| `attention_required` | Conflict collision, automatic flag, or patch/run mismatch | Resolve all findings before use |

If a score table is absent, the report says automatic-approval status was **not
verified**; it does not claim that zero flags were checked. Topology remains
optional, applicability-gated evidence. Low topology coverage is reported, not
silently treated as missing biological truth.

## Scale and reproducibility

Only the highest-priority `--max-embedded-candidates` are included in the HTML;
complete conflict sets are added when a selected member would otherwise be
orphaned. `report.json`, `candidates.tsv`, and `review_ledger_template.tsv`
always contain the full candidate universe.

The workbench follows system light/dark preference, exposes visible keyboard
focus, includes a skip link and semantic table controls, and disables animation
when reduced motion is requested. Release tests measure key WCAG 2.2 contrast
pairs and parse the JavaScript independently with Node.js. Printing suppresses
interactive controls and keeps the evidence panels readable on paper or PDF.

Set `SOURCE_DATE_EPOCH` to obtain a deterministic `generated_at_utc` value for
reproducible packaging. Then verify the package:

```bash
cd work/minimal-reviewed-patch-report
sha256sum -c SHA256SUMS
```

The HTML uses an explicit Content Security Policy and contains no external
JavaScript, stylesheet, font, image, telemetry, or network dependency.
