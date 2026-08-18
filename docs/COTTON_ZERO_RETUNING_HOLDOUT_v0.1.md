# Cotton zero-retuning holdout v0.1

Date frozen: 2026-08-07

Status: candidate policy, normalized inputs, evaluator-only homeolog pairs,
the 800-event benchmark, mature-method candidates, self-WGD context, and the
zero-retuning soybean-model transfer are frozen. No cotton result altered the
rule described below. The independent SynGAP comparator is still running and
will be added as a separately frozen extension.

The machine-readable freeze is
`config/copy_collapse_zero_retuning_policy_v0.1.tsv`. It fixes exact 2-of-3
method-family CDS topology agreement as the primary candidate layer, 3-of-3 as
a secondary high-precision layer, WGD as context rather than a hard gate, and
the soybean model as a review ranker. Automatic copy-addition approval is
disabled. Any post-truth change invalidates the formal holdout claim.

## Selection

The next external annotation-copy-collapse benchmark will use allotetraploid
Gossypium hirsutum as the target. The staged evidence panel contains:

| Role | Species | Biological use |
|---|---|---|
| target | G. hirsutum | A01-A13 and D01-D13 target assembly and annotation |
| A-genome reference | G. arboreum | independent A-lineage projection and collinearity |
| D-genome reference | G. raimondii | independent D-lineage projection and collinearity |
| allotetraploid comparator | G. barbadense | mature-tool comparison and sensitivity analysis only |

All four filesets are the local PGCP/Phytozome exports under
`/nas_data/NFS/Public_genome_data/PGCP/data`. Their original compressed files
are copied into the PloidyPatch project and frozen by byte count, SHA-256, and
`gzip -t`; the public tree remains read-only.

Cotton was chosen because its named A and D subgenomes permit a recent
allopolyploid test that is mechanistically closer to annotation copy collapse
than maize paleopolyploidy. It also supplies two biologically independent
diploid reference lineages. Arachis hypogaea is deferred because the available
PGCP diploid progenitor annotations are tied to highly fragmented sequence
names, whereas the target Tifrunner assembly is newer and chromosome-level.
Zea mays remains a later model-plant portability set, and hexaploid wheat is a
later scale and ambiguity stress test rather than the first external holdout.

## Leakage boundary

Before the candidate rule is frozen, permitted work is limited to:

1. staging and checksum verification;
2. lossless decompression and primary-chromosome normalization;
3. format audits that do not enumerate target homeolog pairs or score hidden
   events; and
4. implementation tests on synthetic fixtures and the existing Brassica and
   Glycine development sets.

After rule freeze, the complete G. hirsutum annotation may be used by the
evaluator process to infer conservative A-D homeolog candidates and generate
the perturbation. The accepted pair table, removed feature IDs, hidden truth,
and complete annotation must remain evaluator-only. Candidate generation may
read only the target genome, perturbed annotation, declared external reference
genomes/annotations, and evidence recomputed from those blind inputs.

No hyperparameter, threshold, feature transformation, reference weighting, or
method-vote policy may be changed after the cotton perturbation is generated.
Any later change creates a new development iteration and requires a different,
untouched external species for a formal claim.

## Frozen truth construction

The evaluator defined A-D pairs using the existing independent-reference
contract:

- target genes are assigned to subgenome A or D from the explicit 26-row
  chromosome manifest;
- each external reference counterpart must connect to exactly one A and one D
  target gene in a WGDI collinear context;
- evidence is collapsed by reference species, with both G. arboreum and
  G. raimondii required;
- target genes must have reciprocal-unique partners across surviving evidence;
- selected pairs are coding-model-scorable and sampled deterministically; and
- the pair artifact is evaluator-only and is not evidence available to the
  repair algorithm.

The two target-reference WGDI runs produced 98,119 (*G. arboreum*) and 127,311
(*G. raimondii*) collinear pair rows. Cross-reference reconciliation produced
15,767 candidate homeolog relationships, of which 8,580 met the evidence
threshold and 8,576 were reciprocal-unique. Only four passing candidates had
multiple partners and were excluded. The frozen accepted-pair TSV SHA-256 is
`79e7f66865ae510a0d83cccafdf88881a53b74b9c938388a15995d3ed5a9650c`.

The deterministic sampler selected the predeclared 800 events. Perturbation is
byte-reversible: no-op recovers 0/800 and the oracle recovers 800/800 exact CDS
chains. The benchmark is
`benchmark/structure/copy_collapse_v0.1/ghi_ad/annotation_copy_collapse_seed20260817`.

## Automated post-candidate evaluation

Once the independent miniprot, GeMoMa, and LiftOn upstreams finish, the frozen
workflow will:

1. freeze all six single-method/consensus candidate and decision hashes before
   truth scoring;
2. recompute blind and complete-control self-WGD context independently, using
   64 threads per arm in parallel;
3. apply the frozen soybean model and frozen review threshold with no cotton
   parameter or threshold tuning, plus zero-refit counterfactual masks that
   separately remove WGD context or method-quality fields; and
4. run paired, event-stratified 20,000-replicate bootstrap intervals for event
   recovery and strict CDS precision/recall/F1.

The model-transfer workflow writes a second candidate freeze before it first
hashes or reads hidden truth. Automatic approval remains disabled regardless of
the observed cotton score. These counterfactual masks are not newly trained
ablation models: they deliberately keep the same fitted coefficients and
threshold, so they quantify deployment-time dependence on an evidence group.
Subgenome identity was not present in the frozen soybean feature schema and
therefore cannot be claimed or retrofitted in this formal holdout.

## Frozen mature-method result

The miniprot, GeMoMa, and LiftOn candidate arms and their exact method-family
consensus outputs were cryptographically frozen before evaluator scoring at:

`results/copy_collapse/holdout/cotton_method_trio_v0.1`

All arms retain every baseline transcript structure. On the 800 hidden events:

| Arm | TP | FP | Precision | Recall | F1 |
|---|---:|---:|---:|---:|---:|
| miniprot | 522 | 215 | 0.7083 | 0.6525 | 0.6792 |
| GeMoMa | 567 | 231 | 0.7105 | 0.7088 | 0.7096 |
| LiftOn | 498 | 288 | 0.6336 | 0.6225 | 0.6280 |
| union (1/3) | 553 | 250 | 0.6887 | 0.6913 | 0.6900 |
| exact support 2/3 | 529 | 155 | 0.7734 | 0.6613 | **0.7129** |
| exact support 3/3 | 419 | 61 | **0.8729** | 0.5238 | 0.6547 |

Twenty-thousand-replicate event-paired bootstrap gives a 2/3 recovery interval
of 0.6288-0.6938. GeMoMa recovers 0.0475 more events than 2/3 consensus (95%
CI 0.0263-0.0688), while the consensus removes enough false candidates to have
the highest observed F1. Its independently resampled strict-CDS intervals are
precision 0.7417-0.8039 and F1 0.6864-0.7384. These intervals support the
precision/recall trade-off but do not establish a statistically paired F1 gain
over GeMoMa.

The predeclared 2/3 rule therefore generalizes as the highest observed F1 on a
third WGD species, but only narrowly. The 3/3 precision interval is
0.8416-0.9021, so even this tier fails the frozen one-sided precision safety
criterion for automatic repair. It remains review-only.

## Frozen model-transfer result

The soybean model, feature transforms, review threshold, and candidate policy
were applied without cotton refitting or threshold tuning. Candidate features
and decisions were frozen before evaluator truth was read at:

`results/copy_collapse/holdout/cotton_glycine_model_transfer_v0.1`

Across 36,849 candidates, 553 match an exact hidden CDS truth. The full model
retains useful ranking performance (average precision 0.7173, ROC AUC 0.9888,
Brier score 0.00696), but the complete end-to-end policy selects 488 novel CDS
chains: 414 TP and 74 FP, for precision 0.8484, recall 0.5175, and F1 0.6429.
It causes zero collateral loss but is substantially below the predeclared
2-of-3 review tier (F1 0.7129). The 20,000-replicate paired event bootstrap
places the 2-of-3 minus model recovery difference at 0.1438 (95% CI
0.1200-0.1688). In the separately declared independent confusion-count
bootstrap, the F1 difference is 0.0701 (95% CI 0.0298-0.1102). The latter is
not a paired-event F1 test, but both analyses reject a chance explanation for
the model's observed deployment deficit under their stated designs.

Zero-refit feature counterfactuals expose the transfer mechanism:

| Frozen-model input | Ranking AP | ROC AUC | Brier score | End-to-end selections |
|---|---:|---:|---:|---:|
| full | 0.7173 | 0.9888 | 0.00696 | 488 |
| WGD context masked | 0.4492 | 0.9622 | 0.01474 | 0 |
| method quality masked | 0.7033 | 0.9838 | 0.00842 | 0 |

The masks keep the fitted coefficients and threshold unchanged; they are not
refitted ablations and do not estimate an independent causal effect. They show
that WGD context carries most of the transferable ranking increment, while the
absolute soybean probability/decision scale is brittle under species shift.
The maximum masked probabilities also fall below the frozen threshold, so an
empty counterfactual selection is expected rather than a scoring failure.

The engineering conclusion is deliberately asymmetric: retain WGD as a
plant-specific ranking signal, reject the current cross-species threshold as
an automatic decision rule, and develop species-normalized event-coherence
features plus abstention in v0.2. Cotton truth may diagnose this failure but
may not be used to retune the frozen v0.1 claim.

## Frozen comparison set

The cotton benchmark must report, on identical hidden events:

- no-op and oracle sentinels;
- miniprot, GeMoMa, LiftOn, and matched SynGAP where technically applicable;
- exact 2-of-3 and 3-of-3 method-family consensus;
- the frozen PloidyPatch calibrated plant-copy model;
- ablations without WGD context, subgenome identity, and method-family
  agreement; and
- strict CDS precision/recall/F1, exact CDS grouping, full-transcript recovery,
  CDS-nucleotide coverage, collateral loss, calibration, abstention, runtime,
  peak memory, and reversible-patch validation.

The primary claim is zero-retuning transfer from the Brassica/Glycine
development regime. G. barbadense projections are reported separately because
their close allotetraploid relationship can inflate transfer performance and
must not be conflated with independent diploid evidence.
