# Algorithm quality strategy v0.1

Date: 2026-08-07  
Status: design freeze for the post-cotton development iteration; it must not
alter the active cotton zero-retuning holdout.

## Executive decision

PloidyPatch should not become another plant annotation-transfer ensemble. Its
defensible algorithmic core is the change in plant evolutionary consistency
caused by a proposed edit: whether adding, splitting, fusing, or rejecting a
gene model makes independently inferred WGD, subgenome, orthology, dosage, and
sequence evidence more coherent without erasing genuine PAV or duplicate
divergence.

The current exact 2-of-3 topology consensus is a strong operational baseline,
not a novelty claim. The current event graph is an auditable schema and
constraint MVP, not yet a genome-wide optimizer. The next algorithm iteration
must produce measurable value beyond both.

## Evidence-backed diagnosis

The frozen soybean copy ranker has strong chromosome-grouped nested ranking
performance, but its largest final coefficients are WGD presence, WGD block
count, and WGD-by-method interactions. Removing WGD features lowers nested OOF
average precision from 0.8641 to 0.4721. Removing method-quality features has
almost no effect (0.8641 to 0.8551) and slightly improves the thresholded OOF
F1. This shows that the model is predominantly recognizing the synthetic
homeolog-recovery setting, not learning a general gene-model quality function.

That signal is not yet portable enough for automatic action. On Brassica, the
same frozen model retains ranking information (AP 0.5879; ROC AUC 0.9804), but
the end-to-end review policy falls to F1 0.6708, below GeMoMa and exact 2-of-3
consensus. No development threshold satisfies the predeclared precision lower
bound for automatic approval. The untouched cotton result confirms the same
boundary: AP is 0.7173, but the unchanged end-to-end policy reaches only
precision 0.8484, recall 0.5175, and F1 0.6429. A zero-refit WGD mask lowers
cotton AP to 0.4492, while a method-quality mask lowers it only to 0.7033.
Thus WGD carries portable ranking information, but the absolute score scale and
decision threshold do not transfer safely.

The immediate implications are:

1. raw WGD membership is context, not a repair gate;
2. absolute block counts are vulnerable to species, annotation-density, and
   candidate-universe shift;
3. method agreement is useful for candidate generation but cannot itself be
   the scientific contribution; and
4. an edit must be evaluated by the relationship coherence it creates, not
   only by features attached to the proposed transcript; and
5. cross-species deployment must separate rank portability from calibration
   portability and fail closed when an external precision guarantee is absent.

## Post-cotton algorithm: event-coherence gain

For each proposed annotation state `h`, compute a before/after evidence delta:

```text
gain(h) = structural_fit_delta
        + evolutionary_relationship_delta
        + expected_dosage_delta
        + independent_reference_delta
        - duplicate_collision_risk
        - assembly_repeat_uncertainty
        - evidence_dependence_penalty
```

The terms are typed evidence summaries rather than an assumed linear truth
model. A calibrated learner may rank the resulting deltas, but hard biological
constraints and abstention remain outside the learner.

### 1. Leave-one-candidate-out WGD evidence

- Infer or reuse blocks from existing genes and external anchors.
- Remove the focal candidate and its candidate-derived alignments before
  computing block membership and quota restoration.
- Attach every homeolog relation to a named WGD event rather than a generic
  `in_synteny` flag.
- Normalize block support within event, chromosome pair, and annotation density
  instead of transferring raw block counts between species.
- Record reciprocal uniqueness, anchor distance, block-edge position, and the
  number of alternative partners.

This prevents the proposed annotation from manufacturing the evidence used to
approve itself and reduces cross-species scale shift.

### 2. Homeolog-partner topology coherence

For the independently inferred existing partner, compare:

- protein alignment coverage and identity;
- CDS length and exon-count ratios;
- intron phase pattern and conserved splice boundaries;
- conserved domains and premature-stop/frameshift status; and
- whether the edit restores a missing relationship without creating a second
  equally plausible partner.

These features are closer to the biological mechanism of a collapsed or
missing homeolog than generic candidate identity scores.

### 3. Subgenome-to-reference lineage consistency

In allopolyploids, preserve reference provenance within each method family.
Method votes remain family-level to avoid pseudo-replication, but the graph
records whether an A-subgenome candidate is independently supported by an
A-lineage reference and likewise for D/C/B subgenomes. Cross-lineage support is
retained as weaker or contradictory context rather than silently merged.

The active cotton holdout cannot use this as a selection feature because it is
absent from the frozen soybean schema. Cotton may only report it as an
evaluator-side stratification. It becomes a predeclared feature in v0.2.

### 4. True absence and assembly uncertainty

A missing annotation and a biological loss must compete as explicit
alternatives. Automatic repair remains prohibited when any of the following is
unresolved:

- an assembly gap or low-mappability interval at the expected locus;
- graph/PAV evidence that the sequence is genuinely absent in the target
  haplotype;
- repeat-driven multi-locus projection without a unique syntenic anchor;
- read depth inconsistent with the expected copy number; or
- mutually incompatible protein, RNA, and whole-genome-alignment evidence.

BUSCO gain alone is never truth because it rewards restoring conserved genes
even when a real loss occurred.

### 5. Dependency-aware evidence aggregation

Every evidence edge must declare both an algorithm family and a biological
source group. Multiple references processed by the same program contribute one
algorithm-family vote; multiple projections derived from the same sequence
alignment cannot be treated as independent. The graph caps correlated support
within each group while retaining the strongest contradiction.

### 6. Genome-wide constrained selection

The next event-graph schema should support explicit conflict and dependency
edges. Selection must enforce:

- at most one mutually exclusive state per biological event;
- no incompatible edits to the same source gene or feature hierarchy;
- all required relationship endpoints exist in the proposed annotation;
- WGD/homeolog assignments satisfy declared event and subgenome constraints;
- no accepted edit introduces an undeclared same-strand CDS collision; and
- ambiguous near-optimal solutions are returned for review rather than broken
  by arbitrary ordering.

Connected conflict components can be solved exactly as small binary programs;
oversized or numerically unstable components fail closed to review.

### 7. Risk-controlled output

The system produces three outputs, not a single thresholded list:

- `auto_patch`: only when a predeclared one-sided precision lower bound is met
  on development data and independently confirmed on a held-out species;
- `review_ranked`: useful candidates below the automatic safety guarantee; and
- `abstain`: hard constraint violation, unresolved absence/assembly state, or
  incompatible near-tied solutions.

Calibration is evaluated by clade, event type, subgenome, repeat context, and
method-support pattern. A high ROC AUC does not authorize automatic repair.

## Validation sequence

1. Finish the active cotton experiment exactly under its frozen v0.1 policy.
2. Use cotton only to diagnose transfer failure modes; do not tune v0.1 to its
   hidden truth.
3. Develop event-coherence features on soybean and Brassica with chromosome-
   grouped nested validation and explicit counterfactuals.
4. Freeze the complete v0.2 feature schema, solver, thresholds, and code hashes.
5. Use an untouched model plant or crop lineage (maize first for duplication
   modes; wheat later for scale) as the next formal external holdout.
6. Add a natural-error truth subset supported by long-read transcript evidence,
   assembly-to-assembly sequence, or blinded expert curation. Controlled
   perturbations remain necessary but are insufficient by themselves.

## Executable v0.2 primitive

The first event-coherence primitive is implemented in the event graph without
changing the frozen cotton v0.1 model. Evidence tables may add three columns:

- `evidence_phase`: `before_edit`, `after_edit`, or `shared`;
- `algorithm_family`: the implementation family that generated the edge; and
- `biological_source_group`: the underlying reference, alignment, library, or
  other biological source from which the observation derives.

For each coding, isoform, and relationship scope, the graph connects edges
that share an independent group, algorithm family, or biological source group.
Within each connected dependency component it retains only the strongest
support and strongest contradiction. Coherence gain is the dependency-adjusted
after-edit logit minus the before-edit logit; raw and adjusted gains are both
reported so evidence inflation is auditable.

This output is intentionally not a fitted classifier. A positive gain enters
`review_ranked`; missing paired state evidence, a non-positive gain, or a hard
constraint enters `abstain`. `automatic_approval` is always false. Legacy graph
decisions remain for compatibility, but the risk output is authoritative for
v0.2 development until external calibration is completed.

The first concrete relationship feature now compares a candidate CDS topology
with its independently inferred existing WGD partner. Correct fold-local
evaluation shows supported AP gains in soybean and Brassica and in
soybean-to-Brassica transfer; Brassica-to-soybean has no supported gain.
Normalized WGD context alone does not explain the improvement. Full methods,
counts, intervals, and the invalidated stacking analysis are documented in
`docs/HOMEOLOG_TOPOLOGY_v0.2.md`.

## Go/no-go gates

The algorithmic paper claim proceeds only if all of the following hold:

- the joint event-coherence method improves the precision-recall frontier over
  the best single mature tool and exact method-family consensus on an untouched
  species;
- the gain is present in more than one plant evolutionary regime and is not
  driven only by WGD membership or synthetic event construction;
- zero collateral loss and exact patch reversibility remain hard gates;
- calibration and abstention are reported with uncertainty, including a valid
  external precision lower bound for any automatic tier; and
- a natural-error set supports the same direction of effect.

If cotton shows only that exact consensus transfers, the defensible outcome is
a rigorous benchmark/workflow or software paper. If event-coherence gain later
survives an untouched species and natural validation, the project has a much
stronger methods-paper claim.
