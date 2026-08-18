# Actinidia reference-design rationale v0.5

Date: 2026-08-08
Applies to: `ploidypatch_actinidia_external_validation_v0.5`

## Decision

The formal reference roles remain exactly those named in the external-species
eligibility registry before the Populus reveal:

- target: *Actinidia chinensis* cv. Red5;
- candidate-only references: *A. eriantha* and *A. rufa*;
- evaluator-only group 1: *Rhododendron simsii*; and
- evaluator-only group 2: *Diospyros oleifera*.

No reference was selected or retained using Red5 WGD-pair yield, candidate
yield, labels or model performance. Neither evaluator can enter candidate
generation or ranking, and neither candidate reference can enter truth.

## Prospective status

Commit `1d44fe442ac0e7a759b18a4bd0ea3c0a7252a555` recorded Actinidia as the
secondary target and recorded all four reference roles at 2026-08-08
07:09:49+09:00. The Populus confirmatory result was recorded later in commit
`fc7b8c4de0f296734adf3ca622a5a68f8f721b58` at 09:45:32+09:00. Subsequent
auditing found no project Actinidia target artifacts or prior Red5 labels,
scores, event/candidate counts or performance results.

The exact Ad-alpha logic below was specified after the Populus result but before
any Red5 pair or candidate was enumerated. We therefore describe this as a
`target-level predeclared untouched secondary replication`, not as an analysis
whose every rule predates Populus.

## Why the target event is Ad-alpha

Kiwifruit comparative genomics distinguishes a younger Actinidia-specific
duplication, Ad-alpha, from older Ad-beta and core-eudicot gamma events. The
primary truth is explicitly limited to Ad-alpha. Both formal evaluator lineages
are outside Actinidia-specific Ad-alpha, so a retained evaluator locus can
support the expected one-counterpart-to-two-target relationship.

Red5 is a diploid, chromosome-scale 29-pseudomolecule target with a public gene
annotation. Its older assembly and annotation are not asserted to be natural
ground truth. The experiment instead creates a masked-annotation pseudo-truth
with exact restoration and oracle sentinels. Natural biological plausibility is
a separate claim from this replication.

## Why retain the two candidate references

*A. eriantha* White v1.0 and *A. rufa* ARU_r1.0 are chromosome-scale congeners
that permit all three upstream method families to project independent gene
models into Red5. They provide a realistic close-relative annotation-transfer
setting without entering evaluator truth.

The public ARU_r1.0 annotation contains substantially more gene models than the
other Actinidia references. That metadata is a reason to preserve, not relax,
the original safeguards: references receive no independent votes within a
method family; candidates must satisfy the same fixed adapters; complete-control
subtraction is mandatory; and individual-reference behavior is descriptive.
Replacing *A. rufa* after seeing candidate yield would compromise the
prospective design.

## Why retain Rhododendron and Diospyros as formal evaluators

*R. simsii* and *D. oleifera* are chromosome-scale Ericales genomes from two
families outside the Actinidia-specific Ad-alpha event. Requiring the identical
unordered Red5 self-WGD pair to be supported by both evaluator groups reduces
dependence on a single assembly, annotation or lineage.

This does not assume identical genome histories. In particular, Diospyros has
an independent lineage-specific Dd-alpha duplication. Published comparisons
show that a Diospyros region can relate to multiple Actinidia regions after the
two Actinidia duplications, and that an Actinidia region can relate to two
Diospyros regions after Dd-alpha. A globally unique Diospyros counterpart would
therefore be biologically inappropriate.

The conservative v0.5 rule handles this per counterpart:

1. every accepted evaluator counterpart must link to exactly two qualifying
   Red5 genes;
2. those two genes must be exactly the same unordered cross-chromosome,
   reciprocal-unique Red5 self-WGDI pair;
3. the identical target pair must have support from both evaluator groups;
4. multiple Diospyros counterparts may support the same target pair, reflecting
   Dd-alpha, but all must be pair-consistent; and
5. any counterpart linked to more than two qualifying Red5 genes is rejected as
   ambiguous, never split, truncated or assigned post hoc.

This is deliberately a high-specificity truth definition. It may reduce the
number or breadth of evaluable events. A shortfall is reported as
`not_evaluable_without_rule_relaxation`; it cannot justify changing references
or multiplicity after enumeration.

## Alternatives considered and rejected

Replacing Diospyros with Camellia could simplify the evolutionary narrative,
but it would change a reference role that was explicitly recorded before the
Populus reveal. Because the exact per-counterpart rule addresses the known
Dd-alpha complication without inspecting Red5 pair yield, v0.5 preserves the
predeclared evaluator design.

Allowing a Diospyros counterpart to link to three or four target genes was also
rejected. Such a rule would not identify a unique hidden target pair without an
additional data-dependent choice. Accepting the best-scoring two of four,
choosing the nearest target genes or selecting a pair by downstream recovery is
forbidden.

No Ks window is tuned on the Red5 pair distribution. The primary event rule is
the exact intersection of self-WGDI and both evaluator groups under the frozen
1:2 counterpart filter. Age distributions may be reported as evaluator-only
descriptive diagnostics but cannot add, remove or relabel events.

## Public evidence consulted before pair enumeration

- Red5 assembly and annotation: Pilkington et al., 2018,
  <https://pmc.ncbi.nlm.nih.gov/articles/PMC5902842/>.
- Ad-alpha and Ad-beta placement in Actinidia and related Ericales: Shi et al.,
  2010, <https://pmc.ncbi.nlm.nih.gov/articles/PMC2924827/>.
- Rhododendron genome duplication and Ericales comparison: Soza et al., 2019,
  <https://pmc.ncbi.nlm.nih.gov/articles/PMC6907397/>.
- *R. simsii* chromosome-scale assembly: Yang et al., 2020,
  <https://pmc.ncbi.nlm.nih.gov/articles/PMC7572368/>.
- Diospyros Dd-alpha and Actinidia/Diospyros syntenic depth: Akagi et al., 2020,
  <https://pmc.ncbi.nlm.nih.gov/articles/PMC7048303/>.
- *D. oleifera* chromosome-scale assembly: Suo et al., 2020,
  <https://pmc.ncbi.nlm.nih.gov/articles/PMC6964648/>.
- *A. eriantha* White chromosome-scale assembly: Tang et al., 2019,
  <https://pmc.ncbi.nlm.nih.gov/articles/PMC6446220/>.
- ARU_r1.0 public assembly metadata: PlantGARDEN,
  <https://plantgarden.jp/en/list/t165716/genome/t165716.G001>.

These sources were used only to define releases, phylogenetic roles and the
event model. No Red5 project pair, candidate or endpoint was inspected.
