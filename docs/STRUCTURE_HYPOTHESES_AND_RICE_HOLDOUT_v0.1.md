# Typed structure hypotheses and rice holdout v0.1

Date frozen: 2026-08-07

This checkpoint evaluates a plant-aware audit layer that converts independent
protein projections into typed gene-structure hypotheses. The subsequent
Tier A fused/split patch validation is reported separately in
`STRUCTURE_PATCH_VALIDATION_v0.1.md`; boundary and missing-internal hypotheses
remain review-only at this checkpoint.

## Novelty boundary and method

PloidyPatch does not reimplement miniprot or claim protein projection as the
novel component. The new layer groups exact phased-CDS chains across declared
evidence sources, rejects chains already present in the annotation, records
the annotation-overlap contradiction, and infers an error type only when one
of four exact, unique topologies is present:

| Hypothesis | Truth-blind topology rule | Proposed future action |
|---|---|---|
| missing internal CDS segment | annotation chain plus exactly one internal CDS segment | replace transcript CDS chain |
| terminal boundary error | one terminal CDS segment differs at exactly one boundary | replace transcript CDS chain |
| split gene | candidate chain equals the disjoint union of two annotated-gene chains | merge split gene models |
| fused gene | two disjoint candidate chains exactly partition one annotated chain | split fused gene model |

Every candidate, rejected topology, emitted hypothesis, input checksum, and
parameter is retained. Hidden truth is accessible only to the independent
hypothesis scorer.

## Frozen evidence tiers

- **Tier A** requires at least two independent support groups and a unique
  exact topology. It is the high-confidence tier.
- **Tier B** requires at least one support group and the same unique exact
  topology. It is an expanded manual-review tier.

`source` and `support group` are deliberately distinct. Multiple accessions or
cultivars can be mapped to one group, preventing correlated same-species
evidence from masquerading as independent replication. With no explicit map,
each source is its own group and the manifest declares `support_unit=source`.

All thresholds were frozen on TAIR10: miniprot identity and query coverage
0.5, intact alignments only, minimum gene-overlap fraction 0.1, exact phased
CDS chains, and unique exact topology. No threshold was changed after seeing
rice results.

## TAIR10 development result

The development pool uses *Arabidopsis lyrata* and *A. halleri* proteins and
excludes TAIR10-derived reference proteins. Each of four structure-error types
contains 200 frozen events.

| Tier | Complete-control projection disagreements | Complete-control typed hypotheses | True / differential hypotheses | Differential precision | Complete events recovered | Event recall (95% bootstrap CI) | CDS chains recovered / 1,000 |
|---|---:|---:|---:|---:|---:|---:|---:|
| B, >=1 source | 27,292 | 815 | 385 / 392 | 98.2% | 385 / 800 | 48.1% (45.0-51.3%) | 473 |
| A, >=2 sources | 3,438 | 195 | 227 / 228 | 99.6% | 227 / 800 | 28.4% (25.6-31.3%) | 269 |

The paired Tier B-minus-Tier A event-recall difference is 19.75 percentage
points (95% stratified bootstrap CI 17.0-22.4; 10,000 replicates). This
supports keeping two output tiers rather than selecting a single global
threshold.

## Rice zero-retuning holdout

The target is public *Oryza sativa* IRGSP-1.0. The formal reference pool uses
three public genomes: Asian wild rice *O. rufipogon* OR_W1943, independently
domesticated African rice *O. glaberrima* V1, and the more distant FF-genome
species *O. brachyantha* v1.4b. The explicit source-group map is
`config/source_groups/osa_phylo_v0.1.tsv` with SHA-256
`a9bf95303978a9d24ec89dab63e7f05bd95cbe730e728842fbe2e5d1b696606c`.

The choice spans distinct Oryza lineages rather than three cultivars from one
species. Published Oryza phylogenomics places AA lineages within roughly 4.8
million years and the AA-FF split around 35.3 million years
([AA-genome analysis](https://pmc.ncbi.nlm.nih.gov/articles/PMC4246335/)); the
*O. glaberrima* genome supports an independent African domestication
([Wang et al. 2014](https://www.nature.com/articles/ng.3044)); and Ensembl
describes *O. brachyantha* as a distant, basal FF-genome Oryza lineage
([Ensembl Plants](https://plants.ensembl.org/Oryza_brachyantha/Info/Annotation/)).

| Tier | Complete-control projection disagreements | Complete-control typed hypotheses | True / differential hypotheses | Differential precision | Complete events recovered | Event recall (95% bootstrap CI) | CDS chains recovered / 1,000 |
|---|---:|---:|---:|---:|---:|---:|---:|
| B, >=1 group | 70,576 | 2,333 | 260 / 266 | 97.7% | 260 / 800 | 32.5% (29.5-35.6%) | 300 |
| A, >=2 groups | 6,011 | 441 | 156 / 157 | 99.4% | 156 / 800 | 19.5% (16.9-22.3%) | 177 |

The paired Tier B-minus-Tier A recall difference is 13.0 percentage points
(95% stratified bootstrap CI 10.8-15.4). Per-event formal rice results are:

| Event type | Tier B true / differential | Tier B recall | Tier A true / differential | Tier A recall |
|---|---:|---:|---:|---:|
| terminal boundary | 31 / 35 | 15.5% | 19 / 20 | 9.5% |
| fused gene | 40 / 40 | 20.0% | 21 / 21 | 10.5% |
| missing internal CDS | 87 / 89 | 43.5% | 54 / 54 | 27.0% |
| split gene | 102 / 102 | 51.0% | 62 / 62 | 31.0% |

All eight formal rice score quality gates pass. Both result checksum sets
verify, all stderr logs are empty, and the grouped rerun reproduces every
scientific metric from the source-count run while changing the manifest to
`support_unit=explicit_group` and writing full group labels into the GFF.

## Same-species cultivar control

IR64, MH63, and ZS97 were also projected as an exploratory pan-genome control.
Counting them naively as three independent sources recovered 206 Tier B and
174 Tier A events, but left 1,404 and 877 complete-control typed hypotheses.
Because all three are *O. sativa* cultivars, these figures are not used as the
formal independent-evidence result. Under
`config/source_groups/osa_cultivars_v0.1.tsv`, all three collapse to one
support group and cannot satisfy Tier A by themselves.

## Interpretation and limitations

1. Exact topology is highly specific to an induced structural contradiction:
   formal rice Tier A yields only one differential false hypothesis in 157.
2. Raw proposal precision is intentionally much lower because complete public
   annotations contain hundreds of reproducible projection disagreements.
   Those are background review hypotheses, not proven false annotations.
3. Boundary hypotheses are the weakest class in both species and should
   require translation, RNA, synteny, or comparative splice evidence before
   automatic application.
4. Protein projection does not establish UTR or full isoform correctness.
5. The separate patch checkpoint now validates reversible fused/split coding-
   model repair. Boundary and missing-internal events still require additional
   evidence, and protein projection still cannot validate UTR reconstruction.

## Audit incident

The first phylogenetic-pool run was accurately recorded with
`support_unit=source`: a CLI argument was connected to the legacy miniprot-gap
branch instead of the new structure branch. Commit `8b87ab2` corrected the
routing and added an end-to-end CLI regression test. The original directory
`osa_irgsp10_phylo_v0.1` was not overwritten or relabeled. Formal conclusions
use the clean rerun `osa_irgsp10_phylo_grouped_v0.1`, code commit `0c68812`.

## Frozen artifacts

- TAIR10 candidates:
  `/data/codexli/projects/PloidyPatch/results/structure_candidates/miniprot_multisource/ath_tair10_v0.1`;
- TAIR10 topology and evaluation:
  `/data/codexli/projects/PloidyPatch/results/structure_hypotheses`;
- formal rice projection:
  `/data/codexli/projects/PloidyPatch/results/projections/miniprot_v0.18/osa_irgsp10_phylo_v0.1`;
- formal grouped rice result:
  `/data/codexli/projects/PloidyPatch/results/heldout_structure/osa_irgsp10_phylo_grouped_v0.1`;
- event bootstraps:
  `/data/codexli/projects/PloidyPatch/results/structure_hypotheses/bootstrap_v0.1`.

The immutable topology-to-patch gate is complete for Tier A fused/split CDS
models. WGD/homeolog, subgenome, dosage, copy-collapse, and independent
transcript-structure evidence remain separate required validation arms.
