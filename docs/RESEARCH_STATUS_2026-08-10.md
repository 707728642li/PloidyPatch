# PloidyPatch research status (2026-08-10)

## Supported claim

PloidyPatch is supported as a plant-specific, review-oriented workflow for
recovering coding-copy structures suppressed by overlap handling in duplicated
genomes. The stable contribution is the chain-preserving candidate workflow,
explicit evidence provenance, human review boundary, and reversible patch.
Duplication topology is optional applicability-gated evidence. No automatic
annotation approval is claimed.

## Prospective external evidence

| System | Scope | Chain-preserving result | Formal outcome |
| --- | --- | --- | --- |
| *Populus trichocarpa* | Full H1/H2 confirmation | 731/800 versus 663/800; delta 0.0850, 95% interval 0.06625 to 0.1050 | Positive; all scientific and safety gates passed. |
| *Actinidia chinensis* | Full composite replication | 406/800 versus 333/800; delta 0.09125, 95% interval 0.07125 to 0.11125 | Composite negative; topology coverage 62.07% was below the frozen 70% floor. |
| *Coffea arabica* | Core H1-only replication | 672/800 versus 597/800; delta 0.09375, 95% interval 0.07375 to 0.11375 | Positive H1-only; zero collateral transcript loss in all six arms. |

The three prospective systems span a diploid tree, a serial-WGD asterid, and
an allotetraploid. All three support the core claim that suppressing overlapping
models discards exact phased CDS chains that a chain-preserving pool retains.
Only Populus supports the complete frozen ranking composite.

## Natural and engineering evidence

- Golden Delicious apple: 433 strict multi-exon candidate structures have
  full-length long-read support after the candidate pool and review ranks were
  frozen; 23 cases were selected by checksum-bound case-study rules.
- Scaling: the largest of six real plant candidate universes contained 76,434
  candidates and 16,011 conflict sets and completed in 70.48 seconds with
  1.16 GiB peak resident memory.
- The reviewed-patch example creates two explicitly accepted additions,
  records `automatic_approval=false`, and reverts byte-identically.
- Candidate construction is compared with matched miniprot, GeMoMa, LiftOn,
  single-method, method-consensus, legacy overlap-suppression, and complete
  annotation control arms where applicable.

## Retained negative and invalid results

- Actinidia remains formally negative under its composite protocol despite
  positive H1 and H2 component estimates.
- The v0.9 stable-reference ranker is retired because only five of seven frozen
  transfer/retention gates passed. It must not be presented as a general
  replacement ranker.
- Walnut is retained as an invalid execution with no biological estimate; it
  is not evidence for or against the method.
- Earlier maize and development results remain in their versioned protocol and
  result documents and are not substituted for the prospective primary tests.

## Evidence map

- Populus protocol and result:
  [protocol](POPULUS_EXTERNAL_VALIDATION_PROTOCOL_v0.4.md),
  [result](POPULUS_EXTERNAL_VALIDATION_RESULT_v0.4.md)
- Actinidia protocol and result:
  [protocol](ACTINIDIA_EXTERNAL_REPLICATION_PROTOCOL_v0.5.md),
  [result](ACTINIDIA_EXTERNAL_VALIDATION_RESULT_v0.5.md)
- Coffea protocol and execution amendments:
  [core protocol](COFFEA_EXTERNAL_CORE_H1_PROTOCOL_v1.0.md),
  [reproducibility amendment](COFFEA_BLIND_REPRODUCIBILITY_AMENDMENT_v1.1.md)
- Natural long-read evidence:
  [apple protocol](APPLE_NATURAL_ONT_VALIDATION_PROTOCOL_v0.4.md)
- Stable applicability interpretation:
  [species applicability framework](SPECIES_APPLICABILITY_FRAMEWORK_v0.6.md)
- Submission-level synthesis:
  [submission readiness](SUBMISSION_READINESS_2026-08-10.md)

Quantitative manuscript statements are bound to tracked compact evidence in
`manuscript/source_data/` and verified by the claim registry tests. Large raw
and intermediate datasets remain outside the repository.
