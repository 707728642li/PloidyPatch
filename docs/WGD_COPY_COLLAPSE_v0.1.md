# WGD-aware annotation copy-collapse checkpoint v0.1

## Scope and claim boundary

This checkpoint tests whether an intact genomic copy that is absent only from
the annotation can be localized and restored in duplication-rich plant
genomes. It covers two deliberately different systems:

- recent allopolyploid *Brassica napus* Da-Ae, with A/C subgenomes and two
  extant progenitor-lineage references; and
- paleopolyploid *Glycine max* v2.1, with conservative within-genome syntenic
  duplicate pairs and *G. soja* as an external reference.

The event is `annotation_copy_collapse`: one annotated member of a declared
duplicate pair is hidden while its assembly sequence remains intact. This is
not an assembly copy collapse. A true assembly collapse changes sequence copy
number and requires depth, graph, k-mer, or phased-assembly evidence; it cannot
be established or repaired by a GFF-only patch.

This checkpoint does **not** establish an automatic PloidyPatch copy-collapse
repair. It establishes a reversible benchmark, leakage-controlled WGD pair
sets, matched mature baselines, and the performance gap that a plant-specific
method must close.

## Pair definitions and leakage control

### Recent allopolyploid pairs

Independent Da-Ae versus *B. rapa* and Da-Ae versus *B. oleracea* WGDI runs
contain 291,618 gene-pair rows. A counterpart must map exactly one A-subgenome
and one C-subgenome Da-Ae gene within each reference comparison. A pair is
accepted only when both independent reference groups support it and the
passing relation is reciprocal-unique.

| Stage | Count |
|---|---:|
| candidate A-C pairs | 10,006 |
| accepted two-group reciprocal-unique pairs | 6,449 |
| below support threshold | 3,553 |
| multi-partner among support-passing pairs | 4 |
| promiscuous anchors rejected upstream | 47,457 |

The frozen pair TSV SHA-256 is
`5aa59261c56f0f8d0f34f300f544a37cac090744a375109f67e0eff1f59725a1`.

### Paleopolyploid pairs

The *G. max* self-WGDI run contains 4,167 blocks and 85,779 collinear pair
rows. Pairs must be cross-chromosome, occur in a block with at least 20 gene
pairs, and be reciprocal-unique after all filters. This yields 10,060
conservative syntenic duplicate pairs.

The first inference exposed a real namespace hazard: WGDI used
`GLYMA_...`, whereas operable GFF feature IDs were `gene:GLYMA_...`. The failed
zero-eligible benchmark is retained, not overwritten. Version 0.2 maps every
one of 55,589 WGDI gene IDs uniquely through the source GFF and retains both
the WGDI and feature-ID namespaces. Pair counts are unchanged, and the mapped
pair TSV SHA-256 is
`dd9d343f8beeaf3e8e14af2004a4f3475cbf1ffeef02c0a103d5ad7d21d15223`.

These are syntenic duplicate candidates associated with the soybean WGD; they
are not claimed to be gene-tree-proven homeologs. In both species, pair TSVs
are evaluator-only benchmark inputs. Candidate run contracts explicitly set
`pair_tsv_candidate_access=false` and `candidate_truth_access=false`.

## Reversible benchmark sentinels

Each species has 800 deterministic held-out events. The no-op and oracle arms
are scored before a baseline is allowed to run.

| Target | No-op complete CDS recovery | Oracle complete CDS recovery | Restoration |
|---|---:|---:|---|
| *B. napus* Da-Ae | 0 / 800 | 800 / 800 | byte-identical |
| *G. max* v2.1 | 0 / 800 | 800 / 800 | byte-identical |

The failed pre-mapping soybean attempt is retained at
`annotation_copy_collapse_seed20260815.failed_gene_id_namespace_v0.1`; the
formal mapped holdout is
`annotation_copy_collapse_seed20260815`.

## Frozen single-source miniprot transfer

Both blind runs use the existing adapter thresholds without copy-collapse or
soybean tuning: identity at least 0.5, query coverage at least 0.5, intact
alignment required, existing-CDS overlap at most 0.2, and candidate redundancy
overlap at most 0.5. Candidate generation uses only the perturbed annotation,
external protein projection, and its protein map. A same-method complete-
annotation control is generated separately and subtracted only by the
evaluator.

| Target | Reference proteins | Accepted raw projections | Differential CDS candidates | TP | FP | Precision | Recall | F1 | Exact CDS grouping | Full transcript | Collateral loss |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Da-Ae | *B. rapa* + *B. oleracea* | 23,241 | 779 | 580 | 199 | 0.7445 | 0.7250 | 0.7346 | 580 / 800 | 6 / 800 | 0 |
| *G. max* | *G. soja* | 6,331 | 668 | 577 | 91 | 0.8638 | 0.7213 | 0.7861 | 577 / 800 | 38 / 800 | 0 |

The coding result transfers: event recovery is almost identical in the two
systems. Specificity does not transfer at the same level, and is substantially
higher in soybean. Therefore one global precision claim would be misleading.

Independent event bootstrap with 20,000 replicates gives Da-Ae recall 0.7250
(95% percentile CI 0.6938-0.7550) and soybean recall 0.7213
(0.6900-0.7525). The Da-Ae-minus-soybean recall difference is 0.00375 with CI
-0.0400 to 0.0475. An independent confusion-outcome bootstrap estimates the
precision difference as -0.1192 with CI -0.1595 to -0.0787, and the F1
difference as -0.0515 with CI -0.0854 to -0.0179. The report identifies the
resampling unit and does not pretend the cross-species events are paired.

## Plant-synteny gates that were rejected

Two plausible plant-specific filters were tested before being accepted as
design features.

| Da-Ae arm | Differential candidates | TP | FP | Precision | Recall | F1 | Decision |
|---|---:|---:|---:|---:|---:|---:|---|
| single-source miniprot baseline | 779 | 580 | 199 | 0.7445 | 0.7250 | 0.7346 | candidate baseline |
| two independent source groups | 679 | 526 | 153 | 0.7747 | 0.6575 | 0.7113 | evidence feature only |
| blind WGDI gap hard gate | 423 | 318 | 105 | 0.7518 | 0.3975 | 0.5200 | reject as hard gate |

Independent source agreement modestly raises precision but sacrifices enough
recall to lower F1. A blind local-WGDI-gap requirement does not materially
raise precision and removes nearly half the recoverable events. Synteny and
subgenome evidence therefore remain contextual features; neither is a
mandatory repair gate.

## Matched mature transfer baselines in soybean

GeMoMa 1.9 and LiftOn 1.0.11 do not consume the target annotation in their
upstream reference-to-target prediction, so their already validated *G. soja*
to *G. max* predictions are reused. They are newly adapted against the
copy-collapse blind annotation and the complete control with the same frozen
overlap thresholds. No pair table or hidden truth is available before
candidate hashes are frozen.

| Method | Differential CDS candidates | TP | FP | Precision | Recall | F1 | Exact CDS grouping | Full transcript | Collateral loss |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| miniprot adapter | 668 | 577 | 91 | 0.8638 | 0.7213 | 0.7861 | 577 / 800 | 38 / 800 | 0 |
| GeMoMa 1.9 | 747 | 612 | 135 | 0.8193 | 0.7650 | 0.7912 | 612 / 800 | 43 / 800 | 0 |
| LiftOn 1.0.11 | 692 | 607 | 85 | 0.8772 | 0.7588 | 0.8137 | 607 / 800 | 6 / 800 | 0 |
| SynGAP 1.2.5 | 514 | 355 | 159 | 0.6907 | 0.4438 | 0.5403 | 355 / 800 | 21 / 800 | 0 |

Paired event bootstrap intervals for complete CDS-chain recovery are:

| Method | Recall | 95% CI |
|---|---:|---:|
| miniprot | 0.7213 | 0.6900-0.7525 |
| GeMoMa | 0.7650 | 0.7350-0.7938 |
| LiftOn | 0.7588 | 0.7288-0.7875 |
| SynGAP | 0.4438 | 0.4100-0.4788 |

Miniprot minus GeMoMa is -0.0438 (CI -0.0613 to -0.0263), and miniprot minus
LiftOn is -0.0375 (-0.0538 to -0.0225). GeMoMa minus LiftOn is 0.00625
(-0.00875 to 0.02125). Thus the simple PloidyPatch miniprot candidate adapter
does not beat mature transfer tools. LiftOn is the current single-tool F1
leader, while GeMoMa has the highest recall. The target-annotation-aware SynGAP
run is also separated from 2-of-3 consensus: SynGAP minus consensus recall is
-0.30875 (paired 95% CI -0.34250 to -0.27500). SynGAP nevertheless recovers
CDS nucleotides at precision 0.9745 and recall 0.7330. This distinction makes
local coding coverage a secondary endpoint rather than evidence of complete
phased-CDS-chain repair.

## Exact method-family consensus baseline

A non-novel ensemble baseline groups accepted miniprot, GeMoMa, and LiftOn
models by exact `seqid + strand + phased CDS chain`. One method family receives
one vote even if it emits duplicate upstream models. The selector verifies
that every input candidate preserves the same base GFF byte prefix, rejects
overlapping alternative consensus chains, and emits CDS-bounded exons with
`coding_only=true`. Both predeclared support tiers were frozen before either
was scored.

| Tier | Blind support-passing chains | Differential chains | TP | FP | Precision | Recall | F1 | Exact CDS grouping | Full transcript |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| at least 2 of 3 methods | 4,466 | 665 | 602 | 63 | 0.9053 | 0.7525 | 0.8218 | 602 / 800 | 38 / 800 |
| all 3 methods | 3,049 | 605 | 564 | 41 | 0.9322 | 0.7050 | 0.8028 | 564 / 800 | 36 / 800 |

All three methods recover 564 truth events, exactly two recover another 38,
and only one recovers 28; all three miss 170. The 2-of-3 tier therefore
recovers exactly the 602 majority-supported truth events. Its recall is 0.0125
below GeMoMa (paired 95% CI 0.00125-0.0250) and is not separated from LiftOn
(LiftOn-minus-consensus difference 0.00625, CI -0.0025 to 0.0150). Its observed
F1 is the highest current value, but the small difference from LiftOn is not
treated as a novel plant-specific result. It is a stronger baseline that the
event-aware method must beat.

## Blind candidate self-WGD reanchoring

To test WGD context without leaking the pre-perturbation pair table, the
2-of-3 candidate GFF was combined with the target genome separately in the
blind and complete-control arms. In each arm, gffread reconstructed proteins,
PloidyPatch rebuilt representative gene-level inputs, and DIAMOND plus WGDI
recomputed self-synteny. Reciprocal pairs were inferred from the mode-specific
extended annotation. A candidate was retained only if its unique partner was
an existing annotated gene; candidate-to-candidate pairs were rejected as
circular support.

| Measurement | Blind | Complete control |
|---|---:|---:|
| genes with representative protein | 59,255 | 59,390 |
| self-WGDI blocks | 4,152 | 4,158 |
| collinear pair rows | 85,103 | 85,409 |
| accepted reciprocal pairs | 10,008 | 10,143 |
| input consensus candidates | 4,466 | 3,801 |
| candidate-existing-gene WGD pairs retained | 750 | 106 |
| candidate-candidate pairs rejected | 31 | 31 |

After paired background subtraction, 645 CDS chains remain: 589 TP and 56 FP,
for precision 0.9132, recall 0.7363, and F1 0.8152. Exact CDS grouping is
589/800, exact full-transcript recovery is 37/800, CDS-nucleotide precision is
0.9946, and collateral loss is zero.

This is not a useful hard gate. Relative to 2-of-3 consensus, reanchoring
removes seven false positives but also thirteen true positives. Recall falls
by 0.01625 (paired 95% CI 0.0075-0.0250), and F1 falls from 0.8218 to 0.8152.
It also has lower precision than the simpler all-three-method tier. Blind
self-WGD reanchoring is therefore retained as an event-graph feature and an
interpretable locus context, not an approval rule.

## Calibrated copy-candidate ranking

A frozen logistic ranking model was trained on 16,638 soybean candidate
chains, with chromosomes kept intact across five outer folds and four inner
folds. It combines method-family support and upstream-quality summaries with
blindly recomputed WGD context. The nested grouped out-of-fold average
precision is 0.8641, ROC AUC is 0.9923, Brier score is 0.00885, and
equal-frequency expected calibration error is 0.00033. The leakage-controlled
cross-fitted review policy has precision 0.8406, recall 0.8772, and F1 0.8585
at the candidate-row level.

The ablations clarify the plant-specific contribution but also its limit.
Average precision falls from 0.8641 to 0.4721 without WGD features and to
0.1655 with method support alone, whereas removing method-quality features
retains 0.8551. WGD context is therefore a real ranking signal in soybean, not
a successful automatic repair rule. No threshold selected from development
folds meets the predeclared one-sided 95% Wilson precision lower bound of 0.90
with at least 30 selected candidates.

For an engineering end-to-end check, the full-fit review rule yields 553 TP
and 40 FP after paired background subtraction: precision 0.9325, recall
0.6913, and F1 0.7940. It sacrifices 49 true event recoveries relative to the
simpler 2-of-3 consensus. Consensus recall minus model recall is 0.06125 with
a paired 95% interval of 0.04500-0.07875. The model is consequently shipped
only as a portable scorer/ranker; it cannot label any copy addition as
automatically approved.

## Brassica multi-reference development and model transfer

The two extant progenitor-lineage references were merged within each method
family before voting, so multiple reference files never count as independent
algorithmic votes. All arms use the same 800 Da-Ae events and a paired
complete-annotation control.

| Method | Differential chains | TP | FP | Precision | Recall | F1 |
|---|---:|---:|---:|---:|---:|---:|
| GeMoMa, merged references | 793 | 618 | 175 | 0.7793 | 0.7725 | 0.7759 |
| LiftOn, merged references | 776 | 540 | 236 | 0.6959 | 0.6750 | 0.6853 |
| method union | 798 | 607 | 191 | 0.7607 | 0.7588 | 0.7597 |
| exact support at least 2/3 | 726 | 594 | 132 | 0.8182 | 0.7425 | 0.7785 |
| exact support 3/3 | 559 | 491 | 68 | 0.8784 | 0.6138 | 0.7226 |

The support-2 tier has the highest observed F1, but it is not a decisive
improvement over GeMoMa. GeMoMa recall exceeds support-2 by 0.0300 with a
paired 95% interval of 0.0125-0.04875. Union recall exceeds support-2 by
0.01625 (0.00875-0.0250). The result supports exact method agreement as a
precision-oriented operational layer, not a claim that the ensemble replaces
the strongest mature predictor.

The frozen soybean logistic ranker was then applied without parameter,
feature, threshold, or calibration changes. Candidate-level average precision
fell from 0.8641 in grouped soybean development to 0.5879 in Brassica; ROC AUC
was 0.9804, Brier score 0.00688, and equal-frequency ECE 0.00748. Thus useful
ordering survives, but the frozen decision threshold does not. After paired
background subtraction, the review arm retains 506 chains: 438 TP and 68 FP,
for precision 0.8656, recall 0.5475, and F1 0.6708. CDS-nucleotide precision is
0.9948 and recall is 0.6447, with zero collateral transcript loss. Support-2
recall exceeds the ranker by 0.1950 (paired 95% interval 0.1675-0.2225).

This Brassica transfer is a secondary portability diagnostic, not a pristine
external holdout. It confirms that ranking metrics alone are insufficient for
a repair decision and freezes the model as review-only.

## Scaling correction

The first full Brassica union exposed a quadratic redundancy scan over 76,434
accepted chains. The selector now indexes accepted models in fixed genomic
bins keyed by sequence and strand while preserving original acceptance order.
On the real support-2 preflight, runtime fell to 12.85 seconds at approximately
241 MB peak RSS. The candidate GFF and decision TSV SHA-256 values were
byte-identical to the original implementation, so the change is a semantics-
preserving scalability correction rather than a method change.

## Product decision

`annotation_copy_collapse` remains review-only. A future automatic policy must
demonstrate held-out value beyond LiftOn and GeMoMa and must add information
that they do not model jointly:

1. locus-specific duplicate quota and reciprocal homeolog/ohnolog identity;
2. subgenome-aware alternative-copy competition rather than a homolog-anywhere
   rule;
3. explicit separation of annotation absence, sequence absence, assembly
   collapse, haplotig redundancy, and biological PAV/loss;
4. calibrated abstention when projections, synteny, depth, or transcript
   evidence contradict; and
5. a reversible copy-addition patch that cannot delete or merge the retained
   partner and that reports coding and isoform confidence separately.

For annotation-only collapse, cross-method agreement can be a confidence
feature but is not itself the novelty claim. For assembly collapse, a distinct
benchmark and sequence/read-depth evidence are mandatory.

Reviewed additions are compiled with
`ploidypatch patch compile-reviewed-copy-additions`. The compiler binds the
complete candidate GFF to the frozen chain-preserving pool manifest and
decisions, then selects only candidates carrying an explicit, provenance-rich
human `accept_addition` decision. It validates complete candidate hierarchies,
identifier non-collision, and at most one accepted alternative per conflict
set. It creates one EOF-only edit that never changes or deletes an existing
feature. The generic patch engine applies the immutable patch and restores the
original annotation byte for byte. The edit artifact explicitly records
`automatic_approval=false`. The exact decision schema and end-to-end command
example are in `docs/REVIEWED_COPY_PATCH_v0.1.md`.

## Reproducibility

Key reusable CLI entry points are:

```bash
# Recent named-WGD pairs from independent reference groups.
ploidypatch evidence infer-homeolog-pairs \
  --gene-evidence gene_evidence.tsv \
  --collinearity bra=bra.collinearity \
  --collinearity bol=bol.collinearity \
  --source-group-map source_groups.tsv \
  --wgd-event brassica_napus_AC_allopolyploidy \
  --subgenome A --subgenome C \
  --min-support-group-count 2 \
  --output-pairs homeolog_pairs.tsv \
  --decisions-tsv homeolog_pair_decisions.tsv

# Conservative within-genome paleopolyploid pairs with GFF ID mapping.
ploidypatch evidence infer-self-wgd-pairs \
  --query-wgdi-gff query.gff \
  --collinearity self.collinearity \
  --source-gff annotation.gff3 \
  --wgd-event glycine_max_paleopolyploid_WGD \
  --min-block-pairs 20 \
  --output-pairs syntenic_duplicate_pairs.tsv \
  --decisions-tsv syntenic_duplicate_decisions.tsv

# Different species/holdouts: independent event resampling.
ploidypatch benchmark bootstrap-independent-events \
  --score brassica=brassica.score.json \
  --score glycine=glycine.score.json \
  --output-json portability.json

# The same hidden events evaluated by different methods: paired resampling.
ploidypatch benchmark bootstrap-events \
  --score miniprot=miniprot.score.json \
  --score gemoma=gemoma.score.json \
  --score lifton=lifton.score.json \
  --output-json matched_methods.json

# Compile only explicit human acceptances from an immutable candidate pool.
ploidypatch patch compile-reviewed-copy-additions \
  --annotation-gff annotation.gff3 \
  --candidate-gff candidate.gff3 \
  --pool-decisions decisions.tsv \
  --pool-manifest candidate.gff3.manifest.json \
  --review-decisions review_decisions.tsv \
  --output-edits-json copy_additions.edits.json
```

- Da-Ae homeolog pairs:
  `/data/codexli/projects/PloidyPatch/results/evidence/homeolog_pairs/brassica_v0.1`;
- soybean mapped self-WGD pairs:
  `/data/codexli/projects/PloidyPatch/results/evidence/wgdi/glycine_self_v0.2`;
- copy-collapse benchmarks:
  `/data/codexli/projects/PloidyPatch/benchmark/structure/copy_collapse_v0.1`;
- miniprot results:
  `/data/codexli/projects/PloidyPatch/results/copy_collapse/miniprot_brassica_v0.1`
  and `miniprot_glycine_v0.1`;
- mature soybean results:
  `/data/codexli/projects/PloidyPatch/results/copy_collapse/gemoma_glycine_v0.1`
  and `lifton_glycine_v0.1`, plus the matched SynGAP result in
  `syngap_glycine_v0.1`;
- method consensus and blind self-WGD reanchoring:
  `/data/codexli/projects/PloidyPatch/results/copy_collapse/consensus_glycine_v0.1`
  and `wgd_reanchor_glycine_v0.1`;
- frozen 20,000-replicate statistics:
  `/data/codexli/projects/PloidyPatch/results/copy_collapse/statistics/portability_and_transfer_trio_v0.1`
  and `glycine_full_method_comparison_v0.1`;
- grouped model training and paired engineering evaluation:
  `/data/codexli/projects/PloidyPatch/results/copy_collapse/model_development/glycine_copy_model_v0.1`
  and `glycine_copy_model_evaluation_v0.1`.

Every listed result has a verified `SHA256SUMS` set. Candidate hashes are
written before evaluator truth is read. The synchronized local and server
suites pass 129/129 tests.
