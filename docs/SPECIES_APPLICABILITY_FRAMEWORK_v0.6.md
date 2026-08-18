# PloidyPatch species applicability framework v0.6

## Purpose

PloidyPatch is not assumed to be equally informative for every plant genome.
The chain-preserving candidate pool can still be useful when evolutionary
context is sparse, but the topology-aware ranker requires a chromosome-scale,
internally coherent target annotation and a resolvable retained-duplication
backbone. This framework separates biological inapplicability from algorithmic
failure before any new holdout candidate, label or performance statistic is
visible.

The framework was written after the frozen Actinidia v0.5 result was revealed.
Actinidia is therefore development evidence for v0.6 and can never be reused as
an untouched v0.6 confirmation.

## Three applicability states

1. `full_topology_evaluable`: assembly, incumbent annotation and the base-only
   fixed self-synteny backbone all pass. The frozen topology ranker may be
   evaluated.
2. `chain_preservation_only`: assembly and annotation pass, but the fixed
   backbone fails one or more biological resolvability gates. Candidate
   generation and the non-topology baseline remain usable; topology is disabled
   for every candidate and is not advertised for that target.
3. `not_evaluable_input_quality`: assembly or annotation integrity fails. The
   target is excluded from formal evaluation; parameters must not be relaxed to
   rescue it.

These states are determined before candidate generation for a prospective
holdout. Candidate yield, labels, truth pairs and ranking performance are
forbidden inputs.

## Assembly gates

- The declared primary chromosome/linkage-group set must match the published
  haploid karyotype or the assembly provider's chromosome table exactly.
- At least 95% of assembled non-gap bases must occur on the declared primary
  chromosomes.
- At least 95% of primary-chromosome bases must be non-`N` sequence.
- Primary seqids must be unique, non-empty and byte-identical across FASTA,
  GFF3 and every downstream normalized artifact.
- Assembly BUSCO completeness must be at least 95% using the preregistered
  lineage dataset. In a known polyploid, duplicated BUSCOs are reported rather
  than treated as assembly contamination.
- If read-backed QV is publicly available, QV must be at least 30. Absence of
  raw reads does not silently pass a QV claim; the metric is recorded as
  unavailable and no QV-specific statement is made.

## Incumbent-annotation gates

- At least 95% of protein-coding genes must lie on primary chromosomes.
- At least 98% of coding transcripts must have an exact, unique GFF-to-protein
  identifier mapping after the frozen alias rules.
- At least 98% of coding transcripts must pass hierarchy, coordinate, strand,
  phase and parent-bound checks.
- At least 95% of representative proteins must be non-empty, start-to-stop
  translatable and free of internal stops under the declared genetic code.
- Protein-mode BUSCO completeness must be at least 90%. A target failing this
  gate may not be called a negative test of missing-copy recovery because its
  incumbent annotation is itself too incomplete to define a stable baseline.
- No fuzzy, prefix or edit-distance identifier repair is permitted.

## Fixed-backbone gates

The backbone is built from the sealed incumbent target annotation and target
genome only. For controlled holdouts this is the perturbed annotation, never
the complete annotation. Candidate projections, evaluator references, truth
pairs and labels are not mounted.

- The target WGD or retained-duplication context must have independent
  biological support from at least two primary publications or one primary
  publication plus an authoritative genome resource.
- WGDI blocks must use at least 20 base-only anchor pairs and cross two distinct
  primary seqids.
- Backbone cells must cover at least 60% of incumbent primary-chromosome coding
  gene midpoints.
- At least 75% of declared primary chromosomes must each contain at least 20
  directed backbone cells.
- Among covered incumbent gene midpoints, at least 70% must resolve to exactly
  one partner chromosome after reverse-block canonicalization. This is a
  property of the target's retained-duplication architecture, not a candidate
  score.
- Backbone construction must be byte-deterministic under input-row permutation
  and must produce the same backbone ID when reverse duplicate WGDI blocks are
  added or removed.

The numerical thresholds above are development policy for v0.6. They may be
audited on already label-seen systems, but must be frozen before metadata for a
new prospective target is staged. They cannot be changed after a new target's
base-only backbone is inspected.

## Candidate-level stability rules

- A candidate never participates in WGDI, the fixed protein database, backbone
  reciprocal uniqueness or another candidate's topology decision.
- Direct overlap with more than one incumbent coding locus is ambiguous and
  causes abstention.
- Gap projection requires two bracketing fixed anchors; block extrapolation is
  prohibited.
- Protein homology is searched exhaustively against the fixed incumbent
  proteome. The set of block-compatible partner gene IDs must contain exactly
  one gene. Multiple partners cause abstention; a bitscore winner is not chosen
  post hoc.
- A missing or ambiguous topology leaves the frozen non-topology score exactly
  unchanged. Topology never grants automatic annotation approval.

## Reporting and failure policy

Every species report must publish all assembly, annotation and backbone gate
values, including unavailable metrics. `chain_preservation_only` is a valid
scope boundary, not a topology success. `not_evaluable_input_quality` is not a
negative biological result. Once a prospective target is selected, failing an
applicability gate cannot be repaired by changing its assembly version,
annotation release, primary seqid set, WGD event, references or thresholds.
Any redesigned test requires a new version and another untouched species.

