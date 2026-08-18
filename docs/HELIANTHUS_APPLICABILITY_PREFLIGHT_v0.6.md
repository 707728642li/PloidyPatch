# Helianthus topology-applicability preflight v0.6

Status: development protocol written after the formal Actinidia v0.5 reveal
and before any Helianthus candidate generation, candidate count, truth-pair
enumeration, label, score or performance inspection.

## Purpose and selection status

`Helianthus annuus` HanXRQr2.0-SUNRISE is the already registered rank-3
repeat-rich scale-stress system in
`config/v0.4_external_species_eligibility_registry.tsv`. It is not a replacement
for Populus or Actinidia and this preflight does not promote it to a formal
confirmatory holdout. The only permitted question at this stage is whether its
assembly, annotation and target-only self-synteny are sufficiently stable for
v0.6 topology features.

The target release and three input archives are byte-locked in
`config/helianthus_applicability_input_sources_v0.6.tsv`. Primary chromosomes
are exactly `1` through `17`, as locked in
`config/primary_seqids/helianthus_annuus_hanxrq_r2.0.tsv`.

## Frozen quality assessment

The decision policy is `config/species_applicability_policy_v0.6.tsv`. BUSCO is
fixed to version 6.1.0 and `eudicotyledons_odb12.2` (1,990 profiles) in both
genome and protein modes. Genome mode uses BUSCO's default eukaryotic Miniprot
pipeline. Dataset bytes, software explicit lock, commands and output hashes
must be retained.

Only the following target-only evidence is permitted:

- archive bytes and SHA-256, FASTA/GFF syntax and exact identifier mapping;
- primary chromosome and non-N assembly fractions;
- genome- and protein-mode BUSCO completeness;
- coding hierarchy, representative-protein and translation consistency;
- a fixed target self-synteny backbone built without candidates, evaluator
  references, truth pairs or labels;
- the label-free backbone coverage and ambiguity metrics preregistered in the
  species-applicability policy.

No candidate reference, candidate count, evaluator reference, truth pair,
label, ranking score or endpoint is allowed in this preflight.

## Decision and failure policy

- All assembly, annotation and backbone gates pass:
  `full_topology_evaluable`; a later frozen v0.6 execution may be designed.
- Assembly/annotation pass but the fixed backbone fails:
  `chain_preservation_only`; Helianthus remains a structure and scale stress
  test and topology features are disabled.
- Assembly or annotation fails: `not_evaluable_input_quality`; no biological
  negative claim is made and no topology result is reported.

Thresholds, input release and outcome meanings may not be changed after the
preflight is inspected. A failed topology-applicability gate may not be rescued
by adding candidates, using candidate-candidate edges, lowering hit or block
thresholds, extrapolating beyond anchors, selecting a best member from an
ambiguous partner set, or substituting another species into this v0.6 test.

## Claim boundary

Passing this preflight establishes only that target-only topology projection is
technically evaluable. It does not establish annotation accuracy, a specific
WGD event, gene-tree orthology, candidate validity or publication-level
performance. Failing it means that the topology component is not suitable for
this release under the frozen rules, not that the plant lineage is unsuitable
for PloidyPatch's chain-preservation component.
