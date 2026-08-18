# Coffea external core H1 protocol v1.0

Status: untouched post-metadata, pre-pair, pre-candidate and pre-label
protocol draft, 2026-08-09. It is not an execution or result record and cannot
be frozen until all 15 input artifacts pass exact byte, syntax and ID checks.

## Confirmatory claim

The sole confirmatory hypothesis is that retaining distinct exact phased-CDS
chains recovers more deliberately hidden Arabica C/E homoeolog events than
suppressing strongly overlapping alternative chains. The estimand is the
paired event-level exact-chain recall difference. A 20,000-replicate paired
event bootstrap with literal seed `20260912` supplies the two-sided percentile
95% interval. H1 succeeds only when both the point estimate and interval lower
bound are greater than zero.

There is no H2, ranker, topology correction, probability, threshold tuning or
automatic approval. Stable-reference ranker v0.9 remains retired. Candidate
scores, event counts and labels cannot influence this protocol.

## Event definition

ET-39 has 11 predeclared homoeolog groups, each with one C-derived (`c`) and
one E-derived (`e`) primary chromosome. A formal event must satisfy all of the
following:

1. the two target genes occur in a target self-WGDI block containing at least
   20 gene pairs;
2. they are on different primary chromosomes, belong to the same frozen group
   1..11, contain exactly one `c` and one `e` member, and are reciprocal-unique
   after all frozen filters;
3. at least one exact Gardenia counterpart maps to exactly these two target
   genes and no other qualifying target gene;
4. at least one exact Ophiorrhiza counterpart independently maps to exactly
   the same unordered target pair; and
5. every accepted counterpart in each evaluator group is pair-consistent.

Zero/one-target counterparts provide no support. More-than-two multiplicity,
different target pairs, alias ambiguity, same-subgenome pairs and discordant
homoeolog groups are rejected without truncation, best-hit rescue or manual
override. YN00 Ks is reported descriptively but never selects events: the
known subgenome labels and two independent event-outside evaluator groups are
the preregistered event discriminator. This avoids fitting a target-specific
Ks peak after seeing the data.

Candidate-only Bu-A and Mauritiana are forbidden from truth construction.

## Annotation and protein compatibility

The normalized full primary-chromosome GFF remains unchanged and is the only
annotation used for projection and complete-annotation controls. Protein-
dependent WGDI additionally receives a separate immutable gene universe. A
gene enters that universe only when a provider protein maps by an exact GFF
CDS `protein_id`, transcript relation or explicit header gene relation. Missing
proteins are reported and the corresponding genes abstain from WGDI; no
protein is inferred, translated, prefix-matched or fuzzily rescued.

Provider-wide duplicate GFF rows may be collapsed only when normalized ID,
feature type, seqid, coordinates, strand, phase and normalized Parent set are
identical. Gene or transcript IDs occurring at different structures invalidate
the input. This accommodates exact database duplication without turning
annotation multiplicity into artificial syntenic genes. The subset, exclusions,
representative proteins, mapping modes and all source/output hashes are frozen
before any WGDI pair enumeration.

Parentless CDS fragments are not gene models: they are excluded with exact
counts and never assigned a parent from name or proximity. Conversely, a CDS
that declares a Parent which cannot be resolved uniquely invalidates the
protein universe. This distinction prevents isolated provider debris from
discarding an otherwise valid chromosome annotation without weakening the
hierarchy checks on real gene models.

## Role firewall

The exact roles are one target (ET-39), two candidate references (Bu-A and
Mauritiana), and two evaluator-only references (Gardenia and Ophiorrhiza).
Candidate execution receives only:

- the target primary-chromosome genome;
- the perturbed target annotation;
- candidate-reference genome/GFF/protein bundles; and
- a sealed truth-free benchmark contract.

It cannot mount evaluator references, complete target GFF/protein, truth,
labels, NAS, or a network. Evaluator execution cannot see candidate scores
until raw predictions, both pool manifests and custody checksums are frozen.
Complete-annotation controls are evaluator-generated only after that barrier.

## Candidate experiment

Miniprot 0.18-r281, GeMoMa 1.9 and LiftOn 1.0.11 each run against both
candidate references. Multiple references within one method family still
count as one method-family vote. The primary pool retains distinct
`seqid+strand+exact phased-CDS chain` candidates; the fixed comparator
suppresses strongly overlapping alternatives. Adapter thresholds remain
identity >=0.5, query coverage >=0.5, intact model required, existing-CDS
overlap <=0.2 and redundancy overlap 0.5.

The combined-reference pool is primary. Bu-A-only and Mauritiana-only arms
are frozen before labels and reported descriptively so that a positive result
cannot hide dependence on a single reference.

## Evaluability and controls

The evaluator requests 800 nonoverlapping events using seed `20260911` and a
global 11-group round robin followed by the four fixed CDS-complexity bins and
SHA-256 rank. Formal evaluation requires:

- at least 500 events;
- at least 17 of 22 target chromosomes and 9 of 11 homoeolog groups;
- at least 20 events in each of `one`, `two_to_three`, `four_to_six`, and
  `seven_plus` CDS-segment bins;
- no-op recovery exactly zero and oracle recovery exactly one;
- byte-identical restoration of the complete GFF;
- identical complete/blind target-genome SHA-256; and
- zero baseline transcript-structure loss in every method, reference and pool
  arm.

A sentinel, checksum, namespace, custody or role violation is an invalid run.
Insufficient events are `not_evaluable`; no rule may be relaxed. H1 or any
zero-collateral gate failure is retained as a formal negative. After truth is
authorized, the same target/version cannot be retried; subsequent method
changes require a new version and another untouched system.

## Freeze boundary

The formal freeze must bind the contract, policy, event definition, 15 source
artifacts, five primary-seqid tables, ET-39 homoeolog table, Gardenia alias
table, protocol/code/environment manifests, three pipeline entries, blind
outputs, source archive and `SHA256SUMS`. Freeze is permitted only after
metadata/syntax preflight and before any WGDI pair enumeration, candidate
count, truth label or performance computation.
