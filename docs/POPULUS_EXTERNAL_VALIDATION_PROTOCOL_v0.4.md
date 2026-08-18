# Populus external validation protocol v0.4

Date frozen: 2026-08-08  
Policy ID: `ploidypatch_populus_external_validation_v0.4`

## Scope and selection freeze

This is the first untouched confirmatory test of the composite PloidyPatch v0.4
ranker. The target is *Populus trichocarpa* JGI assembly v4.0 with annotation
v4.1. No Populus WGD pair, truth event, candidate, label, model score or
performance endpoint was enumerated before this protocol. The permitted
preflight was limited to source-file SHA-256, compressed-file parsing, sequence
headers, feature counts, primary-chromosome coverage and exact identifier-alias
coverage. Its machine-readable report explicitly records
`wgd_pairs_enumerated=false`, `candidate_counts_computed=false` and
`truth_labels_accessed=false`.

The contamination registry excludes every species whose labels, target-specific
structure, or protocol suitability were previously viewed. Three metadata-only
eligible systems remained: Populus, Actinidia and Helianthus. The deterministic
selection rule chooses the smallest chromosome-scale target assembly that has a
lineage-specific WGD and two completely separated candidate and evaluator
reference sets. Populus is therefore the unique primary target. Actinidia is a
predeclared secondary replication and Helianthus a future repeat-rich scale
stress; neither may replace a failed or unevaluable Populus result under v0.4.
Pair yield, candidate count and model performance are forbidden selection
variables.

## Frozen reference roles

| Role | Frozen release | Permitted use |
|---|---|---|
| target | *P. trichocarpa* JGI v4.0, annotation v4.1 | complete evaluator source; blind target genome and perturbed annotation |
| candidate reference | *Salix purpurea* JGI v5.0/v5.1 | miniprot, GeMoMa and LiftOn candidate generation only |
| candidate reference | *Salix suchowensis* GCA_017552425.1 | miniprot, GeMoMa and LiftOn candidate generation only |
| evaluator reference | *Manihot esculenta* v6 | evaluator-only duplicate support |
| evaluator reference | *Ricinus communis* Wild castor | evaluator-only duplicate support |

Candidate references can never enter truth construction. Evaluator references
can never generate, feature or rank candidates. Role changes are forbidden even
if the event yield or accuracy is disappointing. All source paths, releases,
primary-sequence regexes, byte counts and SHA-256 values are frozen by the input
source table and metadata preflight.

Identifier normalization permits only an exact unique alias found in the
declared GFF3 fields `ID`, `Name`, `gene_id`, `locus_tag`, `protein_id`,
`orig_protein_id` or `transcript_id`, or the corresponding exact protein-header
field. Prefix stripping, fuzzy matching, nearest-coordinate repair and manual
alias exceptions are forbidden after the freeze.

## Evaluator-only salicoid duplicate truth

The target event is the salicoid WGD retained in Salicaceae after its split
from the evaluator Euphorbiaceae lineages. One target pair must pass three
independent streams:

1. Populus self-WGDI, with target genes on different primary chromosomes, at
   least 20 gene pairs in the supporting block and reciprocal uniqueness after
   filtering;
2. Populus-Manihot collinearity, in which one exact Manihot counterpart links
   to exactly two Populus genes in qualifying blocks; and
3. Populus-Ricinus collinearity with the same exact 1:2 multiplicity.

The final truth is the exact unordered Populus feature-ID pair intersection of
self-WGDI and both evaluator outgroups. Salix relationships cannot support a
truth pair. Candidate and evaluator references must use separately prepared
directories and manifests. The evaluator pair tables, decisions and complete
target annotation are never mounted inside the blind runner.

DIAMOND, WGDI, representative-protein selection and collinearity parameters
remain those frozen for the apple study: DIAMOND `evalue=1e-5`, no more than 20
targets and `--more-sensitive`; WGDI score 100, grading `50,40,25`, `mg=40,40`,
`pvalue=0.2`, `repeat_number=20`, gene-order positions and at least 20 block
pairs. Versions, commands, inputs and outputs are hashed.

## Event sampling and invalidation sentinels

The evaluator requests 800 nonoverlapping events with sampler seed `20260930`.
It uses global target-chromosome round robin, then round robin over one,
two-to-three, four-to-six and seven-plus CDS-segment bins, followed by SHA-256
rank. Exactly one partner per pair is deterministically hidden; the genome is
unchanged.

The run is formally evaluable only with at least 500 events, at least 15 target
chromosomes represented and at least 20 events in every complexity bin. A
shortfall is `not_evaluable_without_rule_relaxation`; references, pair rules and
endpoints remain unchanged.

Before candidate performance is evaluated, all sentinels must pass:

- the blind no-op arm recovers zero complete hidden phased-CDS chains;
- the complete oracle recovers every event;
- evaluator restoration recreates the complete source GFF3 byte-for-byte; and
- blind and complete target genome FASTA files have identical SHA-256.

Any sentinel failure is an invalid run, not a negative performance result.

## Blind generation, controls and custody

The blind runner is an unprivileged bubblewrap mount namespace created with all
namespaces separated and network access disabled. It receives only read-only
mounts for frozen code/model/environment, blind target genome/GFF3 and the two
candidate-reference bundles, plus a new writable output directory. It receives
no mount for `/nas_data`, evaluator references, truth/pair directories, labels,
or the complete target annotation. Mount information, runner identity and input
and output hashes are captured in a custody manifest before evaluator reveal.

Miniprot 0.18-r281, GeMoMa 1.9 and LiftOn 1.0.11 are run independently. Two
Salix references still contribute only one vote per method family. Frozen
adapter thresholds are identity 0.5, query coverage 0.5, intact projection,
existing-CDS overlap at most 0.2 and redundancy overlap at most 0.5. The primary
pool retains distinct phased-CDS chains and records conflict sets; the legacy
arm suppresses strongly overlapping alternatives from the identical raw inputs.

Complete-annotation control adaptation is evaluator-owned and is generated only
after blind raw predictions are checksummed. Thus a blind upstream tool cannot
condition on the complete Populus GFF3. Differential candidates subtract
complete-control background. Every primary, legacy, individual-tool, support-2
and support-3 arm must lose zero baseline transcript structures.

Blind self-WGD, features, v0.3 scores and the production v0.4 conflict guard are
computed before labels. The composite artifact is frozen by SHA-256 and contains
the v0.3 model, immutable guard policy, source archive, explicit conda
environment, development evaluation and production replay. The v0.4 output is
an uncalibrated review rank, never a probability or automatic approval.

## Ordered hypotheses

### H1: chain-preserving candidate ceiling

An event is recovered when a differential candidate exactly matches the hidden
seqid, strand and ordered phased CDS chain. The endpoint is recall of
`retain_distinct_chains - suppress_overlapping`. A paired 20,000-replicate event
bootstrap with seed `20261001` reports a two-sided percentile 95% interval. H1
passes only when the observed difference and lower interval bound are positive.

### H2: v0.4 review ranking

H2 is tested only if Populus is formally evaluable and H1 passes. A candidate is
positive only when its exact phased-CDS chain matches hidden truth after paired
complete-control subtraction. The endpoint is average precision of the frozen
v0.4 guard minus the frozen v0.3 baseline. A 20,000-replicate target-chromosome
grouped bootstrap with seed `20261002` gives a two-sided percentile 95%
interval. H2 passes only when its point difference and lower bound are positive.
At least 19,000 chromosome bootstrap replicates must be valid.

The v0.4-guard versus v0.3-primary bootstrap uses seed `20261003` and is
descriptive only. It cannot replace H1 or H2.

## Safety gates and secondary outputs

Every confirmatory pass additionally requires:

- the complete per-conflict v0.4 winner mapping hash equals the frozen baseline
  mapping hash and mismatch count is zero;
- v0.4 top-1% true positives are not below baseline;
- at least 70% of positive candidates retain effective topology evidence after
  guard abstention;
- the v0.3 AP gain is positive and v0.4 retains at least 90% of it;
- all candidate/control arms have zero collateral transcript loss;
- all automatic-approval fields are zero; and
- the bootstrap-validity requirement passes.

Review budgets are exactly ceil(0.5%, 1%, 2%) and min(N, 100, 250, 500), with
ties resolved by descending score then candidate digest. Each budget reports
reviewed count, true positives, precision, positive-candidate recall and the
ordered selected-digest SHA-256. Secondary outputs include ROC AUC, PR curves,
v0.3, v0.2, method support, maximum quality, individual method arms,
support-2/support-3, per-chromosome, complexity and conflict MRR. None can
substitute for the ordered primary hypotheses.

## Outcome and failure policy

Sentinel, checksum, manifest or custody violation means `invalid_run`.
Insufficient fixed-rule data means `not_evaluable_without_rule_relaxation`.
Failure of H1, H2 or any safety gate is retained as a
`formal_negative_external_result`. No tuning, reference replacement, looser
truth, alternative endpoint or repeated v0.4 Populus attempt is permitted after
reveal. A scientific or code change requires v0.5 and a different untouched
system. Before reveal, a pure execution crash may restart only from identical
frozen bytes; any code change creates a new patch freeze and preserves the old
attempt.

All counts, hashes, versions, event strata, AP/ROC values, bootstrap deltas,
invalid replicates, review yields, collateral results and failure reasons are
reported regardless of direction.
