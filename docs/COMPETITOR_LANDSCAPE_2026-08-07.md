# Competitive landscape and novelty boundary — 2026-08-07

## Bottom line

A generic “plant TOGA2,” “synteny-based plant annotation polisher,” or
“orthology-guided gene-model optimizer” would duplicate published work. The
project remains defensible only if its primary contribution is an explicit,
calibrated plant event model that jointly reconciles gene structure,
allele/ortholog/homeolog/duplicate relationships, and loss/PAV uncertainty in
polyploid and multi-assembly settings.

This is a go/no-go constraint, not marketing language. If controlled and
held-out tests do not show value beyond SynGAP, OMGene, current transfer tools,
and sequential orthology/synteny pipelines, the project should stop or pivot.

## Direct and near-direct prior work

| Method | Existing contribution | Material overlap | Boundary for PloidyPatch |
|---|---|---|---|
| [TOGA2](https://www.biorxiv.org/content/10.64898/2026.06.30.735536v1.full) | Exon-wise annotation from whole-genome alignment, splice-model correction, gene-tree reconciliation, loss/duplication/retrogene calls, multi-reference integration; demonstrated on 2,162 vertebrate assemblies | Joint annotation and orthology inference is already established | Do not port or copy it. Plants require explicit WGD nodes, subgenomes, alleles, local duplicate modes, graph/PAV evidence, and plant-trained splice evidence. TOGA2 also acknowledges reference dependence and inability to discover lineage-specific genes/exons from homology alone. |
| [SynGAP](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-024-03359-8) | Plant/animal synteny-gap detection followed by genBlastG/miniprot model recovery; dual, master, and triple modes; published in Genome Biology | Automated synteny-based gene-structure polishing in plants is not novel | Must exceed pair/triple gap filling: model WGD quota/subgenome context, split/fusion/collapse/haplotig states, contradictory evidence, true absence, calibrated abstention, and feedback into relationship inference. |
| [OMGene](https://pmc.ncbi.nlm.nih.gov/articles/PMC5923031/) | Optimizes existing gene models by maximizing amino-acid alignment agreement within orthogroups | Orthology-guided model optimization is not novel | Avoid “maximize conservation” as the core objective; it can erase real plant divergence. Preserve event alternatives and require independent sequence/RNA/assembly evidence. |
| [OrthoFiller](https://pmc.ncbi.nlm.nih.gov/articles/PMC5437544/) | Uses multi-species orthogroups to find missing gene annotations | Missing-gene recovery from orthology is not novel | Treat as a baseline/candidate generator, not a method contribution. |
| [maize split-gene method](https://pubmed.ncbi.nlm.nih.gov/32264824/) | Detects split/merged annotations across B73, PH207, and W22 and uses ten-tissue RNA-seq plus simulations to select the supported state | Plant split/fusion correction and expression validation are published | Benchmark directly; generalize only through a broader event model and orthogonal plant evolutionary constraints. |
| [PAP](https://www.sciencedirect.com/science/article/pii/S2468014126001792) | Graph-pangenome homology annotation with variant-aware gene models across six representative species | Graph-assisted, population-scale structural annotation is direct overlap | Keep graph construction/projection modular; focus on evolutionary event reconciliation, not graph annotation itself. Full method/code audit is still required because the publisher page presented a CAPTCHA and no definitive repository was located on 2026-08-07. |
| [GET_PANGENES](https://pmc.ncbi.nlm.nih.gov/articles/PMC10552430/) | Calls plant pangenes from whole-genome alignments and projects clustered cDNA/CDS with GMAP; tested on Arabidopsis, rice, wheat, and barley | Consistent plant pangene sets and PAV representation already exist | Use as a pangene/PAV baseline; the new claim must distinguish biological PAV from annotation, assembly, phasing, and copy-collapse artifacts. |
| [Pandagma](https://pmc.ncbi.nlm.nih.gov/articles/PMC11377846/) | Uses homology, synteny, and synonymous divergence to build pan-gene/family sets while accommodating WGD | Plant WGD-aware pan-gene grouping already exists | Do not claim WGD-aware grouping alone; infer explicit repair/event states and calibrated uncertainty. |
| [GENESPACE](https://elifesciences.org/articles/78526) | Synteny-constrained orthology and copy-number tracking across complex genomes, designed with plant WGD in mind | Multi-genome plant synteny/orthology backbone already exists | Use as an evidence adapter/baseline. Novelty must come from altering gene hypotheses and re-evaluating evolutionary relationships. |

## Component tools that should be reused or benchmarked

| Function | Mature methods | Planned role |
|---|---|---|
| Annotation transfer | [Liftoff](https://pmc.ncbi.nlm.nih.gov/articles/PMC8289374/), [LiftOn](https://pmc.ncbi.nlm.nih.gov/articles/PMC11118573/), GeMoMa | Candidate models and matched baselines |
| Evidence-based/de novo annotation | [BRAKER3](https://pmc.ncbi.nlm.nih.gov/articles/11216308/), Helixer, GeneCAD, [PlantGeneAnn](https://www.biorxiv.org/content/10.64898/2026.06.25.733695v1) | Independent candidates; PlantGeneAnn is especially relevant as plant-specific splice/structure evidence, not as the reconciliation engine |
| Orthology and duplication nodes | [OrthoFinder](https://www.nature.com/articles/s41592-026-03126-6) | Initial gene trees/species tree and sequential baseline |
| WGD/subgenome context | [WGDI](https://www.sciencedirect.com/science/article/pii/S1674205222003707), [SubPhaser](https://pubmed.ncbi.nlm.nih.gov/35460274/), DupGen_finder | User-supplied or inferred WGD event, chromosome/subgenome, and duplicate-mode evidence |
| Annotation quality control | [OMArk](https://pmc.ncbi.nlm.nih.gov/articles/PMC11738984/), BUSCO/compleasm, gFACs, AGAT | Input/output diagnostics; never sole truth |
| Manual review | [GSAman](https://doi.org/10.1016/j.xinn.2026.101471), Apollo/JBrowse | Review candidates and build blinded truth subsets; PloidyPatch should export compatible tracks, not reproduce an editor |

## Evaluation of the reference TOGA2 study

### What is strong

- The method solves annotation and orthology together instead of assuming a
  finished annotation.
- Exon-level orthology reduces average memory 513-fold and runtime 6.1-fold
  versus TOGA1 while retaining special handling for intron loss and divergent
  exons.
- It integrates splice prediction, gene-tree reconciliation, reference-absent
  introns, terminal-exon repair, UTRs, retrogenes, loss status, and
  multi-reference annotations.
- The scale demonstration is unusually strong: 9,376 runs over 2,162
  vertebrate assemblies with downloadable annotations and alignments.

### What should not be copied or over-interpreted

- The implementation and biological assumptions are vertebrate-centric;
  SpliceAI was trained on human and the paper validates it across vertebrates,
  not plants.
- The reference-query chain architecture and conventional 1:1/1:many labels do
  not explicitly represent plant WGD events, subgenomes, homeologs, tandem
  expansions, or polyploid dosage constraints.
- Homology transfer cannot discover lineage-specific exons or genes; the
  authors state this limitation directly.
- Reference annotation quality remains a major determinant, while
  multi-reference merging mitigates but does not eliminate reference bias.
- Much of the broad completeness comparison relies on compleasm/BUSCO. Exact
  exon concordance to RefSeq is useful but is not independent truth for all
  model changes. Controlled perturbations, precision of proposed repairs,
  calibrated abstention, and independent long-read validation are still
  needed.

## Frozen novelty claim to test

PloidyPatch will test whether an explicit plant evolutionary event graph can
improve both gene models and relationships by jointly selecting among:

- accepted, split, fusion, missing-exon/boundary, missing-annotation,
  collapsed-copy, haplotig-duplicate, pseudogene, true-absence,
  assembly-uncertain, and abstain states; and
- allele, speciation ortholog, WGD-homeolog tied to a named event,
  tandem/proximal/transposed/dispersed duplicate, lineage-specific, and
  uncertain relationship states.

The claim survives only if the joint model beats the best sequential
combination of SynGAP/LiftOn/GeMoMa plus GENESPACE/Pandagma/OrthoFinder on
held-out clades, particularly WGD, subgenome, tandem-repeat, and PAV strata.

## Immediate design consequences

1. Add SynGAP, OMGene, PlantGeneAnn, GSAman, and PAP to the mandatory audit and
   baseline list.
2. Make a sequential “best components” pipeline the main baseline, not a weak
   single-reference transfer.
3. Freeze Glycine and wheat as held-out WGD stress tests before fitting.
4. Treat uncertainty calibration and true-loss discrimination as primary
   outcomes, not optional reports.
5. Do not begin a large model implementation until a small perturbation study
   demonstrates a measurable gap left by SynGAP and sequential baselines.
