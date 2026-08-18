# Decision log

## 2026-08-07 — working name

Decision: use **PloidyPatch** as the working project/tool name.

Rationale: it is concise, plant/polyploid-relevant, and emphasizes reversible
repair. Exact-name searches of literature, GitHub, PyPI, and Bioconda returned
no match. This is a preliminary collision check, not formal trademark advice.

Rejected working names:

- PhytoMend: registered agricultural trademark.
- PhytoGRAFT: used commercially and less specific to ploidy-aware inference.
- PlantOGA: too derivative of TOGA and risks misrepresenting the method as a
  plant port.

## 2026-08-07 — product boundary

Decision: build an annotation-repair and orthology-reconciliation layer, not a
new ab initio annotator or graph aligner.

Rationale: GeMoMa, LiftOn/Liftoff, BRAKER3/Helixer, GET_PANGENES, PAP,
Pandagma/GENESPACE, and OrthoFinder already cover major component functions.
The defensible gap is closed-loop plant event inference with calibrated
loss/PAV uncertainty.

Update after full landscape review: SynGAP already provides automated
synteny-based plant gene-structure polishing and OMGene already optimizes gene
models using orthogroup conservation. Therefore neither synteny-guided repair
nor orthology-guided optimization is independently novel. The project proceeds
only around explicit WGD/subgenome event constraints, multi-class duplicate and
assembly-artifact states, calibrated abstention, true loss/PAV discrimination,
and feedback from repaired structures into relationship inference.

## 2026-08-07 — source data handling

Decision: keep the 630 GB source collection in place and access it read-only.

Rationale: copying raw/clean RNA-seq and BAM files would waste storage and
weaken provenance. Project directories will contain manifests, indexes where
appropriate, code, logs, and downstream results only.

## Pending decisions

- Obtain and inspect the complete PAP paper, supplementary material, and code.
- Select the open-source license after auditing all linked dependencies.
- Choose the first independent public benchmark releases and freeze their
  train/tune/test roles before download.

## 2026-08-07 — corrected public-data source

Decision: treat `/nas_data/NFS/Public_genome_data/` as the authoritative local
public-genome pool. Retain the previously inspected hawthorn directory as a
separate discovery cohort.

Rationale: mixing the private/project hawthorn cohort with the public benchmark
pool would obscure data provenance and could create circular validation.

## 2026-08-07 — server project root

Decision: use `/data/codexli/projects/PloidyPatch` as the server project root,
with code, prefix environments, staged data, logs, benchmark artifacts, and
results kept in separate subdirectories.

Rationale: the user explicitly authorized creation of a project directory on
the server. This location follows the shared server convention and isolates all
project writes from both the public source pool and the hawthorn source cohort.

## 2026-08-07 — environment provisioning warning

Observation: while creating `envs/ploidypatch-dev`, mamba reported one cached
repodata shard fetch failure caused by an unexpected TLS EOF. The transaction
completed successfully from available metadata/packages, and the installed
environment passed the project test suite.

Decision: record the warning and retain the configured channels. Do not switch
package sources unless the network/TLS failure recurs and the user directs a
source change.

## 2026-08-07 — copy-specific PAV evidence

Decision: treat whole-assembly protein matches as location-dependent evidence.
A local match contradicts sequence absence, an unlocalized-scaffold match
forces abstention, and an off-locus primary-chromosome match is retained as
homeolog/paralog evidence rather than being allowed to erase a copy-specific
absence hypothesis.

Rationale: in the Da-Ae/Darmor allotetraploid test, a location-naive screen
found intact homologs for 233 of 255 candidates. Most were not at the inferred
missing locus. A generic "homolog anywhere means present" rule therefore
conflates WGD-derived copies with locus presence and is unsuitable for plants.

Decision: report PAV at both gene and structural-event levels. The first 24
WGD-anchored gene candidates collapse to 13 candidate events, including one
eight-gene interval. Event counts, not gene counts, are the primary unit for
structural PAV evaluation.

## 2026-08-07 — paired perturbation evaluation

Decision: use a same-method, complete-annotation control when measuring
controlled perturbation precision. Subtract only novel structures reproduced
in both the complete and perturbed runs; never expose the complete control to
candidate generation.

Rationale: the real NCBI annotation is not guaranteed complete. Counting every
reproducible method/source disagreement as a synthetic false positive confounds
the perturbation response with natural background discoveries. Paired scoring
reduced 23,025 formal miniprot calls to 564 deletion-responsive structures and
produced reproducible exact-CDS metrics in a disjoint development set.

Decision: exact phased CDS chain is the primary strict endpoint for protein-
projection components. Full transcript structure remains reported, but it is
not a fair primary endpoint for a method that does not infer UTRs.

## 2026-08-07 — role of synteny evidence

Decision: reject strict local WGDI-gap membership as a mandatory repair gate.
Retain WGDI/subgenome state and multi-progenitor support as confidence and
event-interpretation features.

Rationale: the frozen formal local gate reduced exact-CDS recall from 30.1% to
9.60% and precision from 57.8% to 48.8%. In contrast, independent gene-level
WGDI strata were highly informative: exact-CDS recall was 53.2% for genes with
expected and cross-progenitor support versus 5.56% for genes with no block.

## 2026-08-07 — independent soybean holdout

Decision: freeze `Glycine max` v2.1 as the target and `Glycine soja`
ASM419377v2 proteins as the external reference before scoring; reuse all
Brassica miniprot and adapter thresholds without soybean tuning.

Rationale: a whole-lineage holdout tests portability better than another
random split within `B. napus`. The paired result reached 79.4% precision and
43.0% recall for exact phased CDS chains, so external projection remains a
valid candidate source. Exact full-transcript recall was only 6.3%, however.

Decision: expose coding-model and isoform-model confidence as separate output
dimensions. Do not promote a CDS-complete projection to a complete gene-model
repair unless independent transcript-structure evidence supports it.

Rationale: the soybean holdout recovered 86.0% of hidden CDS nucleotides and
81.2% of splice junctions but exactly restored only 8 of 800 complete gene
structures. Combining those endpoints into one label would be scientifically
misleading.

## 2026-08-07 - typed structure audit tiers

Decision: expose two exact-topology audit tiers. Tier A requires at least two
independent support groups; Tier B requires at least one. Both require a unique
exact CDS-chain topology and remain hypotheses until a reversible patch is
compiled and validated.

Rationale: on 800 TAIR10 development events, Tier B recovered 385 events and
Tier A recovered 227. Their paired differential precisions were 98.2% and
99.6%, respectively. The event-level recall difference was 19.75 percentage
points (95% stratified bootstrap CI 17.0-22.4). A single threshold would discard
either useful review sensitivity or a defensible high-confidence output.

Decision: count explicit phylogenetic/support groups rather than raw reference
files. Multiple cultivars or accessions from one species may be collapsed to a
single group by a required, checksummed mapping.

Rationale: IR64, MH63, and ZS97 are correlated *O. sativa* cultivar evidence.
The formal no-retuning rice test instead uses grouped *O. rufipogon*,
*O. glaberrima*, and *O. brachyantha* references. It obtains 156/157 correct
Tier A differential hypotheses and 260/266 Tier B hypotheses, while preserving
the source and support-group identities in every candidate artifact.

## 2026-08-07 - automatic exact-topology patch policy

Decision: permit automatic structure patches only for Tier A fused-gene and
split-gene hypotheses. Require at least two explicit support groups, a unique
exact CDS topology, successful topology revalidation at compile time, complete
single-transcript target-gene hierarchies, no shared external parents, and no
overlapping source-line edits. Missing-internal and terminal-boundary events
remain review-only.

Rationale: complete-control annotations produced no Tier A fused/split typed
hypotheses in either development or holdout, while the blind patch benchmark
converted every selected hypothesis into a correct complete CDS repair:
42/42 fused and 98/98 split in TAIR10; 21/21 fused and 62/62 split in rice.
All four arms had strict-CDS precision 1.0, complete introduced-error removal
for every compiled hypothesis, and zero baseline transcript-structure losses.

Decision: label these outputs as coding-model repairs, not complete isoform
repairs. Replacement exon spans are derived from supported CDS segments;
protein projection does not establish UTR boundaries.

Rationale: exact full-transcript recovery was 0/42 and 1/98 in TAIR10 and
0/21 and 6/62 in rice, despite perfect CDS repair among the selected events.
Combining the two confidence dimensions would materially overstate the result.

## 2026-08-07 - WGD pair identity and annotation copy collapse

Decision: keep recent allopolyploid homeolog candidates and paleopolyploid
syntenic duplicate candidates as separate evidence types. Require
reciprocal-unique relations and retain both upstream evidence IDs and operable
GFF feature IDs in every mapped pair artifact.

Rationale: the first soybean self-WGDI artifact used `GLYMA_...`, while the
source annotation used `gene:GLYMA_...`; an otherwise valid benchmark had zero
eligible genes. Version 0.2 maps all 55,589 WGDI IDs uniquely without changing
the 10,060 inferred pairs. The failed attempt is retained as an identifier-
namespace sentinel. Syntenic duplicate pairs are not mislabeled as gene-tree-
proven homeologs.

Decision: distinguish annotation copy collapse from assembly copy collapse.
The annotation benchmark hides one member of a declared duplicate pair while
leaving assembly sequence intact. Assembly collapse requires copy-number,
depth, graph, k-mer, or phased-assembly evidence and is outside a GFF-only
repair claim.

Rationale: external protein projection can recover an intact unannotated
locus, but cannot show that two genomic copies were collapsed into one
sequence. Combining these failure modes would make both evaluation and repair
claims uninterpretable.

Decision: reject local WGDI-gap membership and two-source agreement as hard
automatic copy-repair gates. Retain them as event-graph features, and keep
annotation-copy-collapse output review-only until a plant-specific policy
beats matched mature transfer tools on held-out data.

Rationale: on the Da-Ae 800-event holdout, the single-source baseline reached
precision 0.7445, recall 0.7250, and F1 0.7346. A two-source tier raised
precision to 0.7747 but lowered F1 to 0.7113; a blind WGDI-gap gate reached
only F1 0.5200. On the soybean holdout, GeMoMa reached F1 0.7912 and LiftOn
0.8137, versus 0.7861 for the miniprot adapter. A projection wrapper or strict
synteny filter is therefore neither the best method nor a defensible novelty
claim.

Decision: use independent resampling for cross-species benchmark comparisons
and paired event resampling only when methods share the exact hidden events.
Report the resampling unit in every statistical artifact.

Rationale: Da-Ae and soybean have distinct event IDs and biological contexts,
whereas miniprot, GeMoMa, and LiftOn are evaluated on the same 800 soybean
events. Treating the former as paired or the latter as independent would
discard the actual experimental design.

Decision: treat exact method-family consensus as a mandatory strong baseline,
not as the plant-specific contribution. Count one vote per independently
implemented method family and require exact phased-CDS topology agreement.

Rationale: on the soybean copy-collapse holdout, 2-of-3 consensus reached
precision 0.9053, recall 0.7525, and F1 0.8218; all-three consensus reached
precision 0.9322, recall 0.7050, and F1 0.8028. The former exceeds the observed
F1 of any single tool, and the latter supplies a simpler high-precision review
tier. Any WGD-aware claim must be compared with both.

Decision: reject recomputed blind self-WGD membership as a hard approval gate.
Permit it as provenance-rich event context and a feature for later calibrated
models.

Rationale: the mode-specific pipeline reconstructed proteins and self-synteny
without accessing the pre-perturbation pair table. It retained 589 TP and 56
FP after paired subtraction (precision 0.9132, recall 0.7363, F1 0.8152).
Relative to 2-of-3 consensus it removed seven FP but thirteen TP; the paired
recall loss was 0.01625 with 95% CI 0.0075-0.0250. This is interpretable WGD
evidence, but not a superior classifier.

Decision: permit reviewed annotation-copy additions to compile into immutable
patches, but never equate patchability with automatic biological approval.
The copy-addition compiler may append complete candidate hierarchies only; it
must not change or delete an existing feature.

Rationale: review workflows still need an auditable delivery format. The
EOF-only edit applies to exactly the reviewed candidate GFF and restores the
source annotation byte-for-byte, while `automatic_approval=false` preserves
the current scientific status.

## 2026-08-07 - copy-ranking and transfer freeze

Decision: keep the calibrated copy model as a review ranker and prohibit
automatic approval. Preserve the simpler exact 2-of-3 method consensus as the
primary operational candidate tier.

Rationale: chromosome-grouped nested validation shows strong ranking signal
(average precision 0.8641), and removing blind WGD features lowers average
precision to 0.4721. However, no development threshold meets the predeclared
one-sided 95% precision lower bound of 0.90. The full-fit engineering score
reaches precision 0.9325 but only recall 0.6913 and F1 0.7940, below consensus
F1 0.8218. Ranking utility and automatic-repair safety are separate claims.

Decision: treat WGD evidence as candidate-universe dependent. Do not infer
that a WGD feature is generally useful because it improves one prefiltered
candidate set.

Rationale: self-WGD reanchoring of 2-of-3 consensus gives F1 0.8152, while
applying the analogous gate to the broader method union gives precision
0.8287, recall 0.7375, and F1 0.7804. The evidence must be evaluated jointly
with the candidate generator and cannot be advertised as a portable hard
filter.

Decision: when multiple reference species feed the same implementation,
namespace their feature IDs and merge them before voting. Count the merged
output as one method-family vote, not one vote per reference.

Rationale: *B. rapa* and *B. oleracea* projections are biologically useful
independent sources but their GeMoMa outputs remain one algorithmic family,
as do their LiftOn outputs. Counting reference files as independent methods
would inflate consensus support and make Brassica incomparable with soybean.

Decision: freeze allotetraploid cotton as the first external zero-retuning
test, but keep its evaluator truth sealed until the Brassica/Glycine feature
schema, model artifact, calibration, thresholds, and method-vote policy are
final.

Rationale: *Gossypium hirsutum* supplies explicit A/D subgenomes plus two
diploid lineage references and therefore tests the intended plant-specific
mechanism. Delaying truth generation prevents cotton results from influencing
the rule that cotton is meant to test. Hawthorn remains excluded.

## 2026-08-07 - cross-species ranker boundary and cotton rule freeze

Decision: freeze exact 2-of-3 method-family CDS topology agreement as the
primary cotton candidate policy. Keep 3-of-3 as a secondary high-precision
tier; retain WGD features and the soybean model for context and review ranking
only. Disable automatic copy-addition approval.

Rationale: in the matched Brassica development comparison, support-2 reaches
precision 0.8182, recall 0.7425, and F1 0.7785, narrowly above GeMoMa F1
0.7759. The frozen soybean model transfers useful ranking (AP 0.5879; ROC AUC
0.9804) but its end-to-end review policy falls to precision 0.8656, recall
0.5475, and F1 0.6708. Support-2 recall exceeds it by 0.1950 with paired 95%
CI 0.1675-0.2225. A cross-species ranker is therefore not an automatic repair
policy.

Decision: generalize evaluator WGDI input handling to declared subgenomes and
root-transcript annotations. Infer subgenomes from the longest declared prefix
of a checksummed chromosome label, and treat root transcripts as gene-level
WGDI units only when an annotation contains no gene records at all.

Rationale: cotton uses A/D rather than Brassica A/C, and the public *G.
arboreum* GFF has top-level mRNA records without gene parents. Hard-coding A/C
or silently producing zero genes would invalidate the external test. The
fallback is explicit in the WGDI manifest and leaves normal gene-parent
annotations unchanged.

## 2026-08-07 - first cotton zero-retuning result

Decision: retain exact 2-of-3 method-family topology agreement as the primary
review tier, but make no automatic-repair or algorithm-superiority claim from
the cotton mature-method comparison.

Rationale: the frozen cotton arm reaches precision 0.7734, recall 0.6613, and
the highest observed F1 of 0.7129, narrowly above GeMoMa F1 0.7096. GeMoMa has
significantly higher paired event recovery by 0.0475 (95% CI
0.0263-0.0688). The ensemble's value is a consistent precision/recall trade-off
across soybean, Brassica, and cotton, not a large performance leap.

Decision: keep the 3-of-3 cotton tier review-only.

Rationale: its precision is 0.8729 with a bootstrap interval of
0.8416-0.9021, below the predeclared precision lower-bound requirement for an
automatic patch tier. Observing fewer false candidates does not override the
frozen safety rule.

Decision: retain WGD context in the post-cotton ranker, but retire the frozen
soybean probability threshold as a candidate for cross-species automatic
repair. Treat v0.1 model output as review ranking only.

Rationale: on the untouched cotton holdout, the full frozen model has ranking
AP 0.7173 and ROC AUC 0.9888, but its end-to-end precision 0.8484, recall
0.5175, and F1 0.6429 are below exact 2-of-3 consensus F1 0.7129. Masking WGD
context without refitting lowers AP to 0.4492, whereas masking method-quality
fields lowers it only to 0.7033. WGD is therefore informative, but an absolute
threshold learned in soybean is not a portable safety guarantee. The paired
event-recovery deficit versus 2-of-3 consensus is 0.1438 (95% CI
0.1200-0.1688); the independent confusion-bootstrap F1 deficit is 0.0701 (95%
CI 0.0298-0.1102).

## 2026-08-07 - homeolog topology coherence

Decision: promote explicit candidate-to-existing-homeolog coding topology into
the v0.2 feature schema. Keep normalized WGD block context as an ablation and
secondary context, not the claimed improvement.

Rationale: with all feature transforms and models refitted inside training
folds, topology raises chromosome-grouped OOF AP by 0.0289 in soybean (95% CI
0.0109-0.0500) and 0.0633 in Brassica (0.0471-0.0788). Normalized WGD context
alone changes AP by -0.0013 and +0.0001, respectively. Soybean-to-Brassica
transfer gains 0.0864 (0.0656-0.1046), while the reverse direction gains only
0.0021 with an interval spanning zero. Pooled development-to-cotton diagnostic
gains 0.0428 (0.0317-0.0533). This is a real but asymmetric ranking gain.

Decision: invalidate any topology evaluation that stacks on an in-sample
soybean model score. Preserve the artifacts for audit but exclude their values
from every scientific summary.

Rationale: a score trained on soybean labels is not an admissible input to a
soybean OOF test or a nominal Brassica-to-soybean transfer test. The corrected
pipeline refits the original feature contract and estimator inside every
training fold or training-species partition and never fits on cotton labels.

## 2026-08-07 - pooled v0.2 ranker and maize external test

Decision: freeze the pooled soybean-plus-Brassica topology estimator as the
primary v0.2 ranker. Include the paired original-feature estimator in the same
artifact as the external comparator. Do not include normalized WGD block
context in the primary estimator, calibrate its sigmoid output as a
probability, set a portable threshold, or permit automatic approval.

Rationale: topology has a supported OOF gain in both development species,
whereas normalized WGD context alone is null and its transfer contribution is
direction-dependent. A rank-only paired artifact isolates the new mechanistic
feature on every external candidate universe without repeating the brittle
v0.1 probability-threshold claim. The frozen model SHA-256 is
`e8d0227628bc0688504f1132cca3cd1e418b87a6dbdff39705d35c9c8f933958`.

Decision: use maize B73 NAM-5.0 as the untouched v0.2 external species, with
Sorghum bicolor NCBIv3 and Setaria italica v2.0 as the two declared outgroups.
The primary gate is a positive topology-minus-baseline AP delta whose
20,000-replicate chromosome-group bootstrap 95% interval excludes zero.

Rationale: maize is a distant monocot with lineage-specific ancient
tetraploidy, so it tests a different evolutionary regime from development and
the post-development cotton diagnostic. The model, exact source hashes,
candidate rules, evaluator-only truth construction, review budgets, topology
coverage floor, and failure policy are frozen before a maize hidden pair or
label is generated. Maize-specific retuning invalidates the formal holdout.

Decision: before maize perturbation, replace the first outgroup-pair artifact
with schema v2 that maps each normalized WGDI gene ID to exactly one source-GFF
gene feature ID. Preserve the v1 artifact and the failed empty sampler under
names containing `invalid_unmapped_wgdi_ids`; do not change collinearity,
support groups, pair evidence, filters, seed, or model.

Rationale: all 3,312 v1 pairs used normalized IDs such as
`Zm00001eb000340`, while the operable Ensembl feature ID is
`gene:Zm00001eb000340`. Consequently the reversible sampler rejected every
pair before generating hidden truth. Exact mapping through the source GFF
`ID`/`gene_id` attributes is a representation correction already used by the
self-WGD module, not maize performance tuning. No hidden event, label, or
model score existed when this correction was made.

Decision: classify LiftOn tracebacks only when they are paired one-for-one
with LiftOn's caught per-gene processing error marker; record their count and
continue only if the final GFF passes LiftOn's own validator plus the frozen
target-coordinate and nonempty hierarchy gates. Segmentation faults, killed
processes, missing commands, or any unclassified traceback remain fatal.

Rationale: the Setaria-to-maize LiftOn run completed, emitted a validator-pass
GFF, and reported one caught empty-protein alignment error for a malformed
reference locus, but the outer wrapper rejected the whole run on the literal
word `Traceback`. The completed working artifact was finalized through the
same post-run validation path using `--resume-completed`; it was not recomputed
or inspected against hidden truth. This guard correction occurred while the
maize candidate pipeline was still blind and before any maize label or model
score existed. The independently running Sorghum arm later completed with two
errors of the same classified type and was finalized by the identical path;
its live shell had additionally encountered an end-of-file parse error because
the shared wrapper was replaced while that shell was still running, after the
LiftOn output and validator report had already been written. Both final
artifacts retain their warning counts and original resource/log files.

## 2026-08-08 - maize natural residual-annotation validation

Decision: use the unmodified official B73 NAM-5.0 annotation for the primary
natural arm. Generate and rank candidates without RNA or pan-genome data,
freeze the complete candidate artifact, and only then validate against
B73-specific full-length Iso-Seq observations and NAM-founder recurrence.
Measure supported loci per fixed review budget rather than treating absent
expression as a negative label. Keep every proposed edit review-only.

Rationale: controlled copy-collapse benchmarks establish recovery under known
truth but cannot show that the system finds useful errors in a real current
annotation. E-MTAB-7837/ERP114624 and Zenodo 2611319 provide processed
full-length transcript sequences plus per-transcript counts that distinguish
pure B73 embryo, root, and endosperm observations. The 25 consistently
annotated NAM founders provide a separate recurrence context, although their
shared pipeline prevents counting them as independent truth. A strict evidence
firewall prevents either source from selecting candidates or tuning the frozen
v0.2 model.

Decision: relegate the archived preliminary `Zm00001e.1` to an optional
temporal silver-standard analysis. It may not be the headline natural result.

Rationale: it offers a same-assembly before/after annotation comparison, but
MaizeGDB explicitly marks it as preliminary and unsuitable for current use.
Success on that deliberately retired annotation would be informative only as
a retrospective revision-recovery stress test and could otherwise inflate the
apparent biological value of the method.

## 2026-08-08 - retain the failed maize v0.2 gate

Decision: retain the maize topology result as a failed formal external gate.
Do not refit, recalibrate, select features, change candidate rules, or replace
the result using maize labels.

Rationale: topology-minus-baseline AP is -0.00173 and its predeclared 20,000-
replicate chromosome-bootstrap 95% interval is -0.00887 to 0.00452. Positive
topology coverage is adequate at 0.7473, so missing feature availability alone
does not excuse the failure. Post-hoc analyses are clearly labeled diagnostic
and cannot convert this test into development data.

Decision: prioritize a chain-preserving candidate-universe contract for v0.3.
Deduplicate only identical phased CDS chains; retain distinct overlapping
chains as mutually incompatible review candidates rather than suppressing one
before scoring. Keep the existing selector behavior as the v0.1-compatible
default and preserve all formal v0.2 artifacts byte-for-byte.

Rationale: at least one maize method recovered the exact hidden structure for
653/800 events, but the current support-at-least-one selector retained only
566. Its CDS-overlap suppression lost 87 recoverable structures, 66 of them
from GeMoMa. Candidate recall cannot be restored by a downstream ranker once a
structure has been removed.

Decision: investigate support-conditioned topology correction on soybean and
Brassica development data only. Use maize solely to explain v0.2 behavior and
require a newly preregistered external species for any v0.3 confirmatory claim.

Rationale: maize topology helps sparse GeMoMa-only candidates but slightly
harms several multi-method strata, while joint refitting of all terms nearly
cancels the biological topology contribution. A conditioned correction is a
mechanistically motivated next hypothesis, but selecting it on maize would
erase the independence of the failed holdout.

## 2026-08-08 - chain-preserving pool and natural Iso-Seq result

Decision: adopt `retain_distinct_chains` as the v0.3 development candidate
contract while retaining the v1 suppressing behavior as the software default
for backward compatibility. Treat every strong-overlap component as a
mutually exclusive review conflict and prohibit automatic approval.

Rationale: paired exact-chain recall increases in soybean, Brassica, cotton,
and maize, with every 20,000-replicate event-bootstrap interval above zero.
The unranked pool has lower precision, as expected, so this result supports a
higher candidate ceiling but not wholesale candidate acceptance.

Decision: retain the natural maize baseline top-K result as positive evidence
for review prioritization and retain the topology comparison as a failed
directional test.

Rationale: 20/100 baseline candidates have independent single-transcript full
multi-exon chain support versus a chromosome-stratified random mean of 1.315.
Topology has 19/100 and topology-minus-baseline is -1 with chromosome-bootstrap
interval -4 to +2. The result validates useful baseline prioritization without
overriding either formal topology failure or pending locus-level safety audits.

## 2026-08-08 - retain v0.3 ranker for a new external freeze

Decision: retain the preregistered support-conditioned offset estimator as the
v0.3 primary ranker. Freeze it only for a new external species; do not use
cotton or maize as confirmatory tests and do not introduce an automatic
approval threshold.

Rationale: in audited repeated chromosome-grouped OOF evaluation, primary AP
increased by 0.0231 in soybean and 0.0469 in Brassica, with both 20,000-
replicate chromosome-bootstrap intervals above zero. Conflict-set top-1
accuracy increased in both species, all five genuinely distinct chromosome
partitions per species gave positive AP differences, and cross-species
development transfers were positive in both directions. Retrospective cotton
improved strongly, whereas maize AP was nearly neutral and its top-review
counts were slightly worse; that limitation is retained rather than explained
away.

Decision: require repeated grouped evaluators to verify that partitions differ
across seeds. Preserve the first v0.3 artifact carrying identical partitions
as invalid for repeated-CV claims.

Rationale: `StratifiedGroupKFold(shuffle=True)` produced the same chromosome
partition for all five requested seeds despite accepting different random
states. The corrected evaluator uses randomized `GroupKFold` and asserts
distinct partitions, no chromosome overlap, and both classes in every fold.
The estimator, features, fixed regularization, species roles, and retention
criteria were unchanged.

## 2026-08-08 - normalize provider syntax without altering structural truth

Decision: keep strict GFF3 parsing as the default and expose provider
compatibility as a separate atomic normalization command. Allow only explicitly
requested repairs with a versioned manifest and exact expected-count checks.
Never drop or alter gene, transcript, exon, or CDS coordinates under a
compatibility option.

Rationale: the frozen Rosaceae inputs contain three unescaped free-text
semicolons, five impossible but redundant intron intervals, and one embedded
FASTA terminator. These defects prevent standards-oriented software from
reading otherwise usable community annotations, but none defines the phased
CDS-chain endpoint. A narrow audited compatibility layer improves real-world
portability while preserving the evidence firewall and making every source
transformation reviewable.

## 2026-08-09 - make topology candidate-independent and species-conditional

Decision: replace candidate-union self-synteny as the v0.6 development
topology source with an immutable target-only backbone. Build reciprocal
target edges and directional anchor cells before candidate access, then
project each candidate independently against an exhaustive target-only protein
search. Candidate-candidate edges, conflict-peer inheritance, extrapolation
outside fixed cells, and winner selection among multiple compatible partners
are forbidden. Multiple partners produce an abstention.

Rationale: the formally retained Actinidia v0.5 result exposed a structural
instability in the earlier implementation. Adding 70,097 candidate genes to
the self-WGDI graph changed reciprocal uniqueness among incumbent target
genes. Of 154 positive candidates without usable topology, 74 had a unique
candidate-to-base edge that was rejected only because other candidates created
additional partners. A topology feature should not change when an unrelated
candidate enters or leaves the review pool. The v0.5 formal-negative result is
immutable; the fixed-backbone implementation is evaluated only as subsequent
development.

Decision: gate topology by a three-state, pre-candidate species-applicability
assessment: `full_topology_evaluable`, `chain_preservation_only`, or
`not_evaluable_input_quality`. Do not lower an applicability threshold, change
species, or inspect candidate/truth statistics to rescue a failed gate.

Rationale: plant genomes differ substantially in assembly continuity,
annotation completeness, retained WGD structure, fractionation, and nested
duplication history. Treating topology as universally available would turn
missing or ambiguous biological context into a ranking artifact. In the
predeclared Helianthus scale-stress system, all assembly and annotation gates
passed (genome BUSCO 98.4%, protein BUSCO 97.3%, exact translation 100%), and
all 17 primary chromosomes contained fixed cells. Nevertheless, those cells
covered only 37.31% of primary-gene midpoints versus the frozen 60% minimum.
Helianthus is therefore retained for chain-preservation and repeat-rich scale
stress, not genome-wide topology ranking. This boundary is an applicability
result, not a negative biological result and not evidence of poor input
quality.

Decision: reject the strict fixed-cell projection as a replacement for the
Actinidia v0.5 topology feature. Retain the implementation as a stable
diagnostic primitive and preserve the preregistered v0.5 formal-negative
result unchanged.

Rationale: in a checksum-frozen post-reveal development replay, the fixed
backbone gave topology to only 138/406 positives (33.99%) and failed the
unchanged 70% coverage gate. Its guarded AP was 0.39495 versus 0.42602 for
v0.5; the chromosome-bootstrap v0.6-minus-v0.5 difference was -0.03107 with a
20,000-replicate 95% interval of [-0.04281, -0.01913]. Conflict-winner
mismatch and automatic approvals remained zero, demonstrating engineering
stability but insufficient biological coverage. Of the 268 unavailable
positives, 258 lacked a strict two-sided anchor bracket; the bottleneck was
backbone spatial coverage rather than genome/annotation quality or homology
search failure. Any subsequent extension must remain candidate-independent,
derive its spatial limits from the target-only backbone before label access,
abstain on ambiguity, and pass the same coverage, ranking, review-budget, and
safety gates before consideration.

Decision: reject bounded single-flank and canonical-block-corridor extensions
as v0.7 topology replacements. Stop extending target-only spatial rules after
this registered comparison rather than weakening ambiguity or extrapolation
limits.

Rationale: the conflict-safe v0.7 union raised Actinidia positive coverage only
from 138/406 to 145/406 (35.71%). Its 264 newly available candidates contained
seven positives and 257 negatives; guarded AP was 0.39651 versus 0.42602 for
formal v0.5, and every fixed review budget was no better. The block-corridor
arm added 124 negatives and no positives. Winner mismatch and automatic
approvals remained zero. Stable target-only spatial evidence is therefore too
sparse in this serial-WGD system, and further geometric relaxation would trade
auditable abstention for label-informed instability.

Decision: retain an exact topology-unavailable fallback, but do not remove
topology from the frozen ranker. Treat homeolog topology as optional additive
evidence whose applicability, coverage, conditional contribution, and safety
are reported separately from the core chain-preserving candidate workflow.

Rationale: forcing every topology feature to the frozen unavailable state made
v0.3 and v0.4 scores exactly equal to baseline in both post-label ablations.
Actinidia AP fell from 0.42602 to 0.38320, and Populus AP fell from 0.85319 to
0.84575, with zero guards, automatic approvals, or winner mismatches. Thus a
complete removal is safe but loses real ranking signal, whereas low topology
coverage should not retroactively negate the separately positive
chain-preservation result. Candidate-reference-anchored topology is the next
candidate-independent mechanism to audit; the rejected target-only rules are
not to be revived by threshold tuning.

## 2026-08-09 - retire the v0.9 stable-reference ranker

Decision: retire the v0.9 stable-reference ranker and retain the
chain-preserving candidate workflow. Reference-anchored topology may be used
only as applicability-gated optional evidence; it is not a universal core
dependency and does not authorize automatic approval. Do not tune the v0.9
ranker on these revealed labels or rescue it by changing thresholds, features,
regularization, seeds, event windows, or species weights.

Rationale: five of seven frozen retention gates passed. Actinidia OOF AP
improved from 0.37318 to 0.40185, with a 20,000-replicate interval for the
difference of 0.01875 to 0.04162. Populus OOF AP changed from 0.89713 to
0.89750, but its interval included zero (-0.00345 to 0.00331), and the frozen
Actinidia-to-Populus transfer reduced AP from 0.82582 to 0.82200. The reverse
Populus-to-Actinidia transfer was positive, every top-1% yield was
non-decreasing, and all winner-mismatch and automatic-approval counts were
zero. This is useful but species-conditional topology signal, not evidence for
a universally retained ranker. The checksum-bound aggregate evidence is in
`manuscript/source_data/stable_reference_ranker_v0.9/`.

## 2026-08-09 - retain Walnut v0.8 as an invalid external run

Decision: classify Walnut run 7 as an invalid external run, not a biological
negative or a not-evaluable system. Do not patch and retry Walnut v0.8, do not
calculate an ad hoc H1 result from its blind pools, and do not use those pools
to tune another rule. Correct the namespace/host published-path verifier only
for a newly frozen untouched system.

Rationale: all six blind method/reference arms, both candidate pools, custody,
and 242/242 blind-root hashes completed successfully. The reveal builder then
rejected the immutable GeMoMa manifest because its correct
`/run/blind-run/project/...working` path was compared directly with the host
`run7/...working` path. The evaluator was not invoked and no H1 estimate was
computed. Although `hidden_truth.json` was not passed to the scorer, truth
reveal had been authorized and the evaluator-only checksum tree had already
been read; the strict no-retry policy therefore applies. The checksum-bound
invalid evidence is in
`manuscript/source_data/walnut_external_v0.8_invalid_run7/`.

## 2026-08-10 - retain the Coffea formal-positive H1 result and close development

Decision: retain the untouched *Coffea arabica* ET-39 core-H1 outcome as a
formal positive external result, without post-reveal tuning or a ranking claim.
Freeze chain-preserving candidate construction as the cross-species core and
stop algorithm development before submission. Do not require another species,
topology replacement or ranker revision as a prerequisite for the paper.

Rationale: in the prespecified combined-reference arm, retaining distinct
phased CDS chains recovered 672/800 events versus 597/800 after overlap
suppression. The paired 20,000-replicate difference was 0.09375 with a 95%
interval of 0.07375 to 0.11375. Bu-A-only and Mauritiana-only arms agreed in
direction, all six arms had zero collateral transcript loss, all evaluability
and sentinel checks passed, and no rule was relaxed. The protocol deliberately
executed no ranker, topology correction or H2 endpoint, so the result cleanly
replicates the candidate-universe mechanism already positive in Populus and
Actinidia. The checksum-bound evidence is in
`manuscript/source_data/coffea_external_v1.0/`.
