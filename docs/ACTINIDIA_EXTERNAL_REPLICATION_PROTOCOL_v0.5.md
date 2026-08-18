# Actinidia external replication protocol v0.5

Date preregistered: 2026-08-08
Policy ID: `ploidypatch_actinidia_external_validation_v0.5`
Framework: `PloidyPatch_external_validation_framework_v0.5`
Frozen estimator: `PloidyPatch_ranker_v0.4`

## Scope and claim boundary

This is the predeclared untouched secondary external replication of
PloidyPatch in *Actinidia chinensis* cv. Red5. The target and the four reference
roles were recorded before the Populus result was revealed. The exact
Actinidia-specific event rules, independent seeds and v0.5 execution framework
are being preregistered after the Populus result but before any Red5 WGD pair,
truth event, perturbed annotation, candidate, label, score, performance endpoint
or event/candidate count is enumerated. The correct claim is therefore
`target_level_predeclared_untouched_secondary_replication`; it is not a claim
that every Actinidia analysis detail predates the Populus result.

Framework v0.5 is a species-aware protocol, evaluator and custody freeze. It
does not retrain, recalibrate or otherwise change the estimator. Candidate
ranking uses the byte-frozen PloidyPatch v0.4 composite artifact. Its output is
an uncalibrated review priority, never a probability or an automatic gene-copy
approval.

Permitted preflight before the formal protocol freeze is limited to public
release metadata, exact source paths, file byte counts and SHA-256 values,
compressed-file parsing, sequence headers, primary-chromosome coverage, feature
counts and exact identifier-alias coverage. It must record all of the following
as false:

- `wgd_pairs_enumerated`;
- `truth_events_enumerated`;
- `candidate_counts_computed`;
- `truth_labels_accessed`;
- `model_scores_computed`; and
- `performance_endpoints_computed`.

Pair yield, candidate yield, class balance and model performance are forbidden
inputs to source selection, protocol adjustment or the decision to execute.

## Frozen reference roles

| Role | Frozen release | Permitted use |
|---|---|---|
| target | *A. chinensis* Red5, `Red5_PS1_1.69.0`, GCA_003024255.1, Ensembl Plants 45 annotation | complete evaluator source; blind genome and deterministically perturbed annotation |
| candidate reference | *A. eriantha* cv. White v1.0, QOVS01000000 | Miniprot, GeMoMa and LiftOn candidate generation only |
| candidate reference | *A. rufa* `ARU_r1.0`, GCA_030159155.1 | Miniprot, GeMoMa and LiftOn candidate generation only |
| evaluator reference group 1 | *Rhododendron simsii* v1.0, GCA_014282245.1 | evaluator-only Ad-alpha duplicate support |
| evaluator reference group 2 | *Diospyros oleifera* v1.0 | evaluator-only Ad-alpha duplicate support with Dd-alpha-aware counterpart handling |

Candidate references can never enter self-WGD inference, truth construction,
event sampling or label generation. Evaluator references can never generate,
filter, feature or rank candidates and are never mounted in the blind runner.
Reference roles cannot change after the formal freeze even if event or candidate
yield is disappointing.

Two candidate references still contribute at most one vote per method family.
Thus annotation density in either candidate reference cannot create extra
cross-method support. Individual-reference yields are reported only as
descriptive diagnostics and cannot change the primary pool.

The formal input source table must bind exact source paths, releases, primary
sequence regular expressions, file sizes and SHA-256 values. Identifier
normalization permits only an exact unique alias in the declared GFF3 or protein
header fields. Prefix stripping, fuzzy matching, coordinate-nearest repair and
manual post-freeze exceptions are forbidden.

## Evaluator-only Actinidia Ad-alpha truth

The target event is the younger Actinidia-specific whole-genome duplication
`Ad-alpha`. The older `Ad-beta`, core-eudicot gamma, tandem and small-scale
duplications, assembly haplotigs and allelic pairs are outside the target truth.
Both evaluator lineages are outside Ad-alpha. *Diospyros* has an independent
`Dd-alpha` duplication, so evaluator support is defined per counterpart rather
than by requiring one globally unique Diospyros gene.

One unordered Red5 feature-ID pair enters strict truth only when the identical
pair passes all three independently prepared evaluator streams:

1. Red5 self-WGDI supports the two genes in a qualifying block containing at
   least 20 gene pairs. The genes are on different Red5 primary chromosomes and
   the target pair is reciprocal-unique after the frozen block filters.
2. At least one exact *R. simsii* counterpart in a qualifying target-evaluator
   block links to exactly two Red5 genes, and those genes are exactly the
   unordered self-WGDI pair.
3. At least one exact *D. oleifera* counterpart in a qualifying
   target-evaluator block links to exactly two Red5 genes, and those genes are
   exactly the same unordered self-WGDI pair.

Every evaluator counterpart is judged separately. A counterpart linked to
fewer than two qualifying target genes does not support an event. A counterpart
linked to more than two qualifying target genes is rejected as ambiguous and
cannot be truncated, locally optimized or assigned to a preferred pair.
Multiple Diospyros counterparts may support one Red5 pair, as expected under
independent Dd-alpha duplication, but every accepted counterpart must resolve
to that same unordered target pair. Pair-consistent reciprocal uniqueness
requires that no accepted counterpart or target member supports a different
target pairing after frozen filters. If an evaluator group has no remaining
exact 1:2, pair-consistent counterpart, the target pair fails strict truth.

The final truth is the exact unordered Red5 feature-ID pair intersection of
self-WGDI, evaluator group 1 and evaluator group 2. No partial-support tier can
enter the primary event universe. Counts of rejected 1:1, greater-than-1:2,
discordant-pair and identifier-ambiguous cases are reported after execution,
regardless of direction, but cannot relax the rule.

DIAMOND, representative-protein selection and WGDI settings remain those used
by the Populus confirmatory test: DIAMOND `evalue=1e-5`, at most 20 targets and
`--more-sensitive`; WGDI score 100, grading `50,40,25`, `mg=40,40`,
`pvalue=0.2`, `repeat_number=20`, gene-order positions and at least 20 block
pairs. Exact software builds, commands, inputs and outputs are hashed in the
formal implementation freeze.

## Event sampling and evaluability

The evaluator requests 800 nonoverlapping events using sampler seed
`20261010`. It applies global target-chromosome round robin, then round robin
over one, two-to-three, four-to-six and seven-plus CDS-segment bins, followed by
SHA-256 rank. Exactly one deterministically selected partner is hidden per pair;
the genome sequence is unchanged.

Red5 has 29 primary chromosomes. The Populus coverage rule is generalized as
`ceil(0.75 * primary_chromosome_count)`, which resolves before execution to 22
Red5 chromosomes. The run is formally evaluable only with:

- at least 500 strict events;
- at least 22 target primary chromosomes represented; and
- at least 20 events in each of the four CDS-complexity bins.

A fixed-rule shortfall is `not_evaluable_without_rule_relaxation`. It is not a
negative result and cannot trigger a looser event definition, replacement
reference or same-target retry.

Before candidate performance is evaluated, all sentinels must pass:

- the blind no-op arm recovers zero complete hidden phased-CDS chains;
- the complete oracle recovers every event;
- evaluator restoration recreates the complete Red5 source GFF3 byte-for-byte;
  and
- blind and complete Red5 genome FASTA files have identical SHA-256 values.

Any sentinel failure is `invalid_run`, not a performance result.

## Blind generation, controls and custody

The blind runner is an unprivileged bubblewrap mount namespace with all
namespaces separated and network disabled. It receives read-only mounts for the
frozen code, v0.4 model, environment, blind Red5 genome/GFF3 and the two
candidate-reference bundles, plus a new writable output directory. It receives
no mount for `/nas_data`, evaluator references, truth or pair directories,
labels, complete Red5 annotation or any prior result.

Miniprot 0.18-r281, GeMoMa 1.9 and LiftOn 1.0.11 run independently. Adapter
thresholds remain identity 0.5, query coverage 0.5, intact projection, existing
CDS overlap at most 0.2 and redundancy overlap at most 0.5. The primary policy
retains distinct phased-CDS chains and records conflict sets. The legacy arm
suppresses strongly overlapping alternatives from the identical raw inputs.

Complete-annotation controls are evaluator-owned and can be generated only
after blind raw predictions, pools, scores, manifests and custody records are
checksummed. Differential candidates subtract complete-control background.
Primary, legacy, individual-tool, support-2 and support-3 arms must each lose
zero baseline transcript structures.

Blind self-topology, features, frozen v0.3 scores and the production v0.4
conflict guard are computed before labels. Custody records runner identity,
mounts, network state and hashes. Any truth, complete annotation, evaluator
reference, NAS or network access inside the blind runner invalidates the run.

## Ordered hypotheses

### H1: chain-preserving candidate ceiling

An event is recovered only by an exact match of seqid, strand and ordered
phased CDS chain. The endpoint is event recall of
`retain_distinct_chains - suppress_overlapping`. A paired 20,000-replicate
event bootstrap with seed `20261011` produces a two-sided percentile 95%
interval. H1 passes only when the observed difference and interval lower bound
are both positive.

### H2: frozen v0.4 review ranking

H2 is confirmatorily tested only if the run is formally evaluable and H1
passes. A candidate is positive only when its exact phased-CDS chain matches
hidden truth after complete-control subtraction. The endpoint is average
precision of frozen v0.4 guard ranking minus the frozen v0.3 baseline. A
20,000-replicate target-chromosome grouped bootstrap with seed `20261012`
produces a two-sided percentile 95% interval. H2 passes only when its point
difference and interval lower bound are both positive. At least 19,000
chromosome-bootstrap replicates must be valid.

The v0.4-guard versus v0.3-primary bootstrap uses seed `20261013` and is
descriptive only. Populus-Actinidia effect comparison or meta-analysis is also
secondary and cannot replace either ordered hypothesis.

## Safety gates and secondary outputs

Every confirmatory pass additionally requires all Populus safety gates without
relaxation:

- the full per-conflict v0.4 winner mapping hash equals the frozen baseline
  mapping hash and mismatch count is zero;
- v0.4 top-1% true positives are not below baseline;
- at least 70% of positive candidates retain effective topology evidence after
  guard abstention;
- the v0.3 AP gain is positive and v0.4 retains at least 90% of that gain;
- every candidate/control arm has zero collateral transcript loss;
- every automatic-approval field is zero; and
- the bootstrap-validity requirement passes.

Review budgets are exactly ceil(0.5%, 1%, 2%) and min(N, 100, 250, 500), with
ties resolved by descending score then candidate digest. Each budget reports
reviewed count, true positives, precision, positive-candidate recall and the
ordered selected-digest SHA-256. ROC AUC, PR curves, v0.3, v0.2, method support,
maximum method quality, individual methods/references, support-2/support-3,
chromosome, complexity, conflict MRR and rejected-counterpart audits are
secondary outputs only.

## Outcome and failure policy

Checksum, manifest, sentinel, role-firewall or custody violation means
`invalid_run`. Insufficient fixed-rule events mean
`not_evaluable_without_rule_relaxation`. Failure of H1, H2 or any safety gate is
retained as `formal_negative_external_result`.

No post-reveal tuning, reference replacement, event relaxation, endpoint
substitution or repeated Actinidia attempt is permitted. Before reveal, a pure
execution crash can restart only from identical frozen bytes. Any code change
requires a new patch freeze that preserves the old attempt. Any scientific
change to model, feature, threshold, reference role, pair rule, sampler,
endpoint or gate requires a new protocol version and a different untouched
target.

All counts, hashes, versions, event strata, accepted and rejected counterpart
classes, AP/ROC values, bootstrap deltas, invalid replicates, review yields,
collateral results and failure reasons are reported regardless of direction.

## Required formal freeze before execution

This text and its machine-readable policy, event definition and contamination
audit are preregistration inputs. Execution is forbidden until an immutable
formal freeze additionally binds:

- exact genome, GFF3 and protein source paths, sizes and SHA-256 values;
- metadata-only preflight report and its `SHA256SUMS`;
- the byte-frozen PloidyPatch v0.4 composite artifact;
- source archive, full Git commit and explicit conda environment;
- species-specific pair, holdout, blind, custody and evaluator implementations;
- a run contract stating `post_metadata_pre_pair_pre_candidate_pre_label`; and
- a checksum manifest covering every protocol and implementation artifact.

The freeze builder must refuse to run if any Red5 pair, truth, holdout,
candidate, score, label or performance artifact already exists.
