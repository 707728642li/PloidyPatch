# PloidyPatch project charter

Date frozen: 2026-08-07
Status: research design v0.1

## 1. Scientific objective

Develop and rigorously validate a plant-specific system that jointly:

1. detects inconsistent or damaged protein-coding gene models across multiple
   assemblies;
2. proposes reversible structural repairs;
3. reconciles allele, ortholog, WGD-homeolog, tandem, proximal, transposed, and
   dispersed-duplicate relationships;
4. separates true gene loss/PAV from missing annotation, assembly gaps,
   haplotig redundancy, and collapsed sequence;
5. returns calibrated confidence and a complete evidence trail.

The target publication level requires a new inference problem and convincing
biological evidence, not only a convenient workflow wrapper.

## 2. Central hypothesis

Plant gene structure and evolutionary relationship are mutually informative.
Current workflows usually predict or transfer gene models first and infer
orthology afterwards. Annotation errors therefore fragment orthogroups and
inflate lineage-specific genes, while incorrect orthology can also select the
wrong repair template.

We hypothesize that an iterative, plant-event-aware model will improve both
gene structures and evolutionary relationships relative to sequential
pipelines, especially after WGD, subgenome fractionation, tandem expansion,
and pangenome PAV.

## 3. Scope and non-goals

### In scope

- Nuclear protein-coding genes in chromosome-scale or draft plant assemblies.
- Population/close-relative and cross-species/clade operating modes.
- Existing GFF3 annotations plus optional RNA-seq, Iso-Seq, protein, read-depth,
  assembly-gap, and repeat evidence.
- Auditing and patching split, fusion, missing-exon, boundary, missing-gene,
  collapsed-copy, redundant-haplotig, and likely-pseudogene cases.
- Explicit WGD/subgenome and local-duplication relationship states.

### Out of scope for v1

- Replacing full ab initio annotation pipelines such as BRAKER or Helixer.
- Building a new whole-genome or pangenome aligner.
- Comprehensive non-coding RNA annotation.
- Functional GO/pathway annotation.
- Fully automatic replacement of original annotations without review.

## 4. Plant event ontology

### Locus structural state

- `accepted`: source model is consistent with available evidence.
- `split_error`: one biological locus is represented by multiple models.
- `fusion_error`: multiple biological loci are represented by one model.
- `missing_exon`: conserved coding segment is absent from the model.
- `boundary_error`: splice, start, or stop boundary is inconsistent.
- `missing_annotation`: intact coding locus exists but is not annotated.
- `collapsed_copies`: multiple copies are represented by one assembly/model.
- `haplotig_duplicate`: alternative haplotypes are represented as extra genes.
- `pseudogene`: disabled coding copy with supported disabling lesions.
- `true_absence`: sequence-level evidence supports loss/PAV.
- `assembly_uncertain`: gap, low coverage, or unresolved repeat prevents a call.
- `uncertain`: evidence is conflicting or insufficient.

### Relationship state

- `allele`
- `speciation_ortholog`
- `wgd_homeolog`
- `tandem_duplicate`
- `proximal_duplicate`
- `transposed_duplicate`
- `dispersed_duplicate`
- `lineage_specific`
- `uncertain_relationship`

WGD relationships must be attached to an inferred or user-supplied event/node,
not used as a generic synonym for paralogy.

## 5. System architecture

### 5.1 Input normalization

- Validate FASTA/GFF3 consistency, coordinate bounds, phases, parent-child
  relationships, translated CDS, duplicate IDs, and sequence-name mappings.
- Preserve every source feature and identifier in an immutable normalized
  schema.
- Select canonical transcripts only for analyses that require them; retain all
  isoforms for structural evidence.

### 5.2 Candidate generation

Use adapters around mature methods rather than copying them:

- protein-to-genome projection: miniprot and/or GeMoMa;
- close-assembly annotation transfer: Liftoff/LiftOn where applicable;
- ab initio orthogonal evidence: BRAKER3/Helixer/GeneCAD outputs;
- synteny anchors: GENESPACE/MCScanX/WGDI-compatible tables;
- gene trees and duplication nodes: OrthoFinder-compatible outputs;
- population graph or WGA evidence: minigraph-cactus/Cactus/HAL or PAP output
  when available.

Candidate generation must maximize recall. It never makes a final biological
decision.

### 5.3 Comparative locus graph

Represent each locus and transcript hypothesis as nodes. Typed edges carry
independent evidence:

- genomic projection and alignment identity/coverage;
- exon/intron-chain concordance;
- protein similarity and domain architecture;
- local and chromosome-scale synteny;
- gene-tree/species-tree reconciliation;
- known WGD quota and subgenome assignment;
- RNA junction, coverage, and full-length transcript support;
- assembly gaps, read depth, and mappability;
- TE overlap and repeat class;
- plant splice-site probability.

The graph must retain contradictory edges so that uncertainty is visible.

### 5.4 Joint inference

Development will proceed in three levels:

1. transparent rule-based candidate filters to validate the event ontology;
2. calibrated per-event scoring models with held-out clade evaluation;
3. constrained structured inference that jointly selects compatible locus
   structures and relationship states.

The final objective will penalize contradictory states such as simultaneous
`true_absence` and strong intact protein projection, or WGD-homeolog calls
without a compatible synteny/event context. It will not force one-to-one
orthology in polyploids.

### 5.5 Iteration and convergence

1. infer provisional relationships from source annotations;
2. detect structural contradictions and produce repair hypotheses;
3. score hypotheses using evidence not used to generate them where possible;
4. update accepted candidate models;
5. rerun affected relationship components only;
6. stop when no high-confidence state changes occur or a fixed iteration limit
   is reached.

Every iteration is versioned and reversible.

### 5.6 Outputs

- patch GFF3, never a silent in-place rewrite;
- machine-readable event and relationship tables;
- per-call confidence, calibration group, and evidence provenance;
- before/after proteins and transcript structures;
- gene-tree, WGD, synteny, and PAV summaries;
- JBrowse/Apollo-compatible review tracks;
- standardized HTML report and run manifest.

## 6. Benchmark design

### Tier A: controlled perturbation

Start from curated annotations and inject realistic errors:

- adjacent split and fusion;
- random and conserved-exon deletion;
- donor/acceptor shifts;
- start/stop truncation;
- tandem copy collapse;
- haplotig duplication;
- pseudogene-like disabling mutations;
- assembly gaps covering full or partial loci.

Perturbations must preserve a hidden truth table and be stratified by gene
length, exon count, expression, repeat content, family size, and synteny.

### Tier B: natural population/close-relative inconsistency

- public *Arabidopsis* reference and close-relative assemblies;
- public *Oryza sativa* IRGSP-1.0, IR64, MH63, and ZS97 assemblies;
- maize B73/W22/PH207 or equivalent multi-reference set;
- another high-quality public crop pangenome if licensing and evidence permit.

### Tier C: plant duplication stress tests

- paleopolyploid: soybean or Brassica lineage;
- recent/allopolyploid: cotton or Brassica napus;
- high-complexity hexaploid: bread wheat as a late-stage stress test;
- monocot/dicot and rosid/asterid coverage to test generalization.

### Leakage prevention

- Freeze species/clade partitions before model fitting.
- Hold out entire species and orthogroups, not random transcripts only.
- Never score a repair with the same projected model used as its sole source.
- Record database release dates and remove target-derived proteins from
  reference pools in strict tests.
- Keep any user/private cohort outside model selection and primary claims
  unless it is intentionally designated later under a separate protocol.

## 7. Baselines and fair comparison

Required baseline families:

- SynGAP and OMGene;
- GeMoMa;
- Liftoff and LiftOn;
- BRAKER3, Helixer, and GeneCAD where input regimes match;
- GET_PANGENES and PAP if runnable code becomes available;
- PlantGeneAnn as a plant-specific candidate source where GPU inference is
  feasible;
- OrthoFiller and the published maize split-gene method where reproducible;
- OrthoFinder plus GENESPACE and/or Pandagma;
- simple majority-vote and single-best-reference ablations.

All methods receive matched assemblies, annotations, protein evidence, and
RNA evidence. Report default and carefully tuned settings separately. Include
wall time, CPU hours, peak memory, disk use, and failure rate.

## 8. Primary evaluation metrics

### Gene structure

- exact gene/transcript match precision, recall, and F1;
- exact exon and donor/acceptor boundary precision/recall;
- translated CDS identity, completeness, and frame validity;
- event-level precision-recall and PR-AUC by error type.

### Evolutionary relationship

- ortholog/homeolog/duplicate-state accuracy against curated or consensus
  references;
- duplication-node placement accuracy;
- orthogroup fragmentation and erroneous merger rates;
- WGD quota and syntenic-retention consistency.

### Loss and PAV

- false PAV rate;
- true-loss recall under independent sequence/read evidence;
- discrimination among loss, gap, missing annotation, and copy collapse.

### Reliability and utility

- Brier score, expected calibration error, and reliability diagrams;
- number of manual reviews needed per accepted correction;
- effect on downstream family-size, PAV, and expression conclusions;
- scaling with genome size, genome count, ploidy, and repeat content.

BUSCO is a secondary completeness descriptor only.

## 9. Statistical analysis

- Predeclare primary metrics and primary datasets before final runs.
- Use paired comparisons on the same perturbed or natural loci.
- Bootstrap by orthogroup or genomic block to avoid treating correlated genes
  as independent replicates.
- Report effect sizes and 95% confidence intervals, not only P values.
- Correct multiple event/dataset comparisons using an explicit FDR procedure.
- Perform sensitivity analyses for reference choice, WGD priors, transcript
  availability, and assembly quality.
- Publish all exclusions and failed runs.

## 10. Engineering plan

- Python CLI and library with a stable normalized schema.
- Workflow orchestration through a versioned Snakemake pipeline.
- Conda prefix environment under the project; containers only after the
  dependency graph stabilizes.
- Unit tests for parsers and event logic, golden integration fixtures, and
  end-to-end benchmark smoke tests.
- Deterministic seeds, checksummed inputs, machine-readable run manifests, and
  structured logs.
- Parquet/SQLite for large tabular state; graph storage chosen only after
  profiling the MVP.
- Optional compiled acceleration only for demonstrated bottlenecks.

## 11. Milestones and decision gates

### M0 — foundation

- final novelty audit, including full PAP paper/code;
- repository, environments, manifests, and source-data QC;
- frozen event ontology and benchmark protocol.

Gate: no unaddressed direct competitor that already performs the same joint
inference. In particular, the proposed contribution must be demonstrably beyond
SynGAP-style synteny-gap polishing and OMGene-style conservation optimization.

### M1 — benchmark harness

- strict FASTA/GFF normalization;
- perturbation generator with hidden truth;
- baseline adapters and reproducible reports.

Gate: all primary baselines reproduce expected behavior on small fixtures.

### M2 — structure-repair MVP

- split, fusion, missing exon, and boundary candidates;
- evidence scoring and reversible patches;
- held-out controlled and public rice/maize evaluation.

Gate: statistically significant improvement over the best matched baseline on
at least two independent datasets, without unacceptable false corrections.

### M3 — plant evolutionary model

- WGD/homeolog, subgenome, tandem/TE, copy-collapse, and loss/PAV states;
- joint iterative reconciliation.

Gate: gains concentrate in biologically relevant difficult strata and survive
ablation of individual evidence types.

### M4 — publication-grade validation

- broad cross-clade tests, independent long-read evidence, runtime scaling,
  calibration, biological case studies, and manual audit.

Gate: all claims are supported by held-out evidence and reproducible public
artifacts.

### M5 — release and manuscript

- documented CLI/API, Bioconda/container candidate, tutorials, archived data
  manifests, versioned benchmark, manuscript, and editorial checklist.

## 12. Genome Biology-level evidence package

A competitive manuscript should contain all of the following:

1. a clearly new joint inference formulation;
2. broad and independent benchmark data, not only one crop;
3. strong current baselines with matched evidence;
4. controlled truth plus independent natural truth;
5. rigorous uncertainty calibration and error analysis;
6. plant-specific gains attributable to WGD/subgenome/tandem/TE modeling;
7. at least one biological conclusion materially changed by corrected models;
8. scalable, documented, tested, openly reproducible software.

The journal target is an ambition and quality gate, not a guaranteed outcome.

## 13. Principal risks

- PAP may overlap more deeply than its public abstract indicates.
- Gold-standard WGD/homeolog and true-loss labels are scarce.
- RNA-seq absence is tissue-dependent and cannot prove gene loss.
- graph alignment can become prohibitive for repeat-rich polyploids.
- evidence sources are correlated and can create false confidence.
- a wrapper-only implementation will not meet the publication objective.

Mitigations are novelty verification, explicit uncertainty states, independent
evidence partitions, clade holdouts, modular alignment backends, and early
go/no-go gates.
