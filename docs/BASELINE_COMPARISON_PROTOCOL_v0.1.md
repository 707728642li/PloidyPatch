# Matched baseline comparison protocol v0.1

Date frozen: 2026-08-07  
Status: pre-execution contract

## Purpose

This protocol defines the comparison that PloidyPatch must pass before a new
plant event model can be claimed as useful. SynGAP is the direct competitor;
GeMoMa and LiftOn are strong homology-transfer controls. The existing miniprot
gap adapter remains the simple component baseline.

The protocol is frozen before installing or scoring the three additional
tools. Parameter changes after aggregate development results are visible must
be declared as tuned analyses and cannot replace the frozen default result.

## Frozen software

| Baseline | Version | Immutable source identity | Role |
|---|---|---|---|
| miniprot-gap | `0.18-r281` | installed project prefix | simple external-protein projection |
| SynGAP | `1.2.5` | Git tag commit `7196d1cb554f089a53b6bf094cb8115893afb7fb` | direct synteny-gap polishing competitor |
| GeMoMa | `1.9` | official `GeMoMa-1.9` distribution | intron-position-aware homology prediction |
| LiftOn | `1.0.11` | Git tag commit `c623f0bddc5b5051a3670a3b0064c52cc61bb719`; PyPI sdist SHA-256 `c2125db9bede3640e13c6aa6a0672887aaf2611710442f9bf69017197d188f88` | DNA plus protein annotation transfer |

The machine-readable registry is `config/baseline_registry_v0.1.tsv`. Each
tool is installed in its own project prefix under `envs/`; no dependency is
added to base or to another analysis environment.

License identity follows the frozen release, not GitHub's current default
branch badge. In particular, the `v1.2.5` SynGAP tag and Bioconda record both
declare `CC-BY-NC-SA-4.0`, although the current repository landing page is
classified as GPL-3.0. SynGAP is therefore an external non-commercial
comparator and no source code is incorporated into PloidyPatch.

## Primary datasets

### Independent primary comparison

- Target: `Glycine max` v2.1, 20 primary chromosomes.
- Reference: `Glycine soja` ASM419377v2.
- Blind input: the frozen 800-gene annotation-missing perturbation with seed
  `20260809`.
- Complete control: the unchanged `G. max` Ensembl Plants 62 annotation.
- Target-derived proteins are forbidden.

This is the primary comparison because every method can receive the same one
reference genome, reference annotation, and target genome without inventing a
multi-reference merge policy.

### Polyploid stress comparison

- Target: chromosome-scale Da-Ae `Brassica napus`.
- References: `B. rapa` and `B. oleracea` Ensembl Plants 62 bundles.
- Blind input: the frozen formal 800-event perturbation with seed `20260807`.

Single-reference runs are reported separately for the A and C progenitors.
Any SynGAP master or other multi-reference merge is a distinct secondary
analysis. It must not be silently compared with a single-reference baseline.

## Information-equivalence rules

Every candidate run receives only:

1. the target genome;
2. the blind target annotation when the method supports annotation polishing;
3. the declared external reference genome, annotation, CDS, and protein files;
4. tool-default sequence databases bundled with the frozen distribution.

It never receives the selection table, hidden coordinates, hidden truth,
complete target annotation, target proteins derived from removed models, or
evaluator-only WGD strata. The complete target annotation is used only for
same-method background control. Methods that consume target annotation
(SynGAP) require a separate complete upstream run. Pure reference-to-target
transfer methods (GeMoMa and LiftOn) receive no target annotation, so their
single frozen raw prediction is adapted independently against blind and
complete target annotations. Running the identical upstream transfer twice
would add no information; the two adapter products are the matched candidate
and background-control arms.

Target and reference inputs are structurally normalized with
`ploidypatch normalize primary-annotation`. This retains every annotation
record on GFF-declared primary chromosomes and excludes unplaced scaffolds;
it does not inspect protein/CDS FASTA files or remove translation warnings.
The resulting manifest freezes source and output hashes, chromosome labels,
record counts, and genome size. Blind and complete target bundles must use the
same genome bytes after normalization.

## Candidate normalization

Tool outputs are not scored directly if that would discard unchanged blind
models. Conventional GFF3 outputs use `ploidypatch baseline adapt-gff`;
miniprot uses its evidence-aware adapter. Each adapter must:

- preserve every blind source record byte-for-byte;
- append only structurally valid protein-coding hypotheses;
- reject models above the same 0.2 existing-CDS overlap threshold;
- collapse equivalent/redundant hypotheses deterministically;
- preserve upstream scores, source IDs, and rejection reasons;
- emit an auditable decisions table and manifest.

The adapter may normalize identifiers and parent links but may not alter model
coordinates or choose a hypothesis using hidden truth. Raw and adapted outputs
are both frozen before evaluation.

## Frozen run modes

### SynGAP 1.2.5

- Primary Glycine mode: `dual`, author biological defaults, default
  `genblastg` prediction backend.
- Secondary Glycine mode: the same `dual` comparison with the documented
  `miniprot` backend and its unchanged SynGAP defaults.
- Primary Brassica modes: two independent `dual` runs, one per progenitor.
- Secondary Brassica mode: `triple` with target plus both progenitors, author
  biological defaults.
- Use `--threads 32`; thread count is treated as a compute setting and must not
  alter predictions.
- Do not use `master`: its separately distributed plant database contains
  `G. max`, which would leak target information into the Glycine holdout, and
  it is not the two-progenitor workflow.
- The polished annotation is adapted only if SynGAP does not already retain
  all blind target records.
- Execute through `scripts/run_syngap_dual.sh`. SynGAP 1.2.5 constructs many
  unquoted `os.system` child commands and does not consistently check their
  return codes, so the wrapper rejects shell-unsafe paths and validates nested
  result directories, clean/full GFF outputs, anchors, completion messages,
  fatal log signatures, and exact preservation of the target GFF as the prefix
  of the full polished output. A zero top-level exit code alone is insufficient.

### GeMoMa 1.9

- Use `GeMoMaPipeline` without RNA-seq for the primary information-matched
  comparison.
- Retain the upstream recommendations `GeMoMa.Score=ReAlign`,
  `AnnotationFinalizer.r=NO`, and `o=true`.
- Use one external reference in Glycine; run A and C Brassica references
  separately before any combined GAF analysis.
- An RNA-supported GeMoMa run is a later evidence-regime experiment, not a
  replacement for this primary baseline.
- Produce one raw target prediction per reference. Apply the frozen GFF
  adapter separately to blind and complete target annotations for paired
  differential scoring.

### LiftOn 1.0.11

- Use current defaults with `--gene-only` for the protein-coding benchmark and
  the documented multi-threaded locus pipeline.
- Current-default and publication-compatible behavior are reported separately
  because versions 1.0.10/1.0.11 changed rescue and merge defaults.
- One reference is used per primary run; multi-reference union is secondary.
- Produce one raw target prediction per reference. Apply the frozen GFF
  adapter separately to blind and complete target annotations for paired
  differential scoring.

### Miniprot-gap 0.18-r281

The already frozen command and adapter thresholds remain unchanged. It is not
rerun or tuned after seeing other baseline results.

## Evaluation

The PloidyPatch evaluator, not an upstream self-evaluator, computes all primary
metrics. The same-method complete-annotation control subtracts reproducible
background proposals before strict scoring.

Primary endpoint:

- precision, recall, and F1 of identifier-independent exact phased CDS chains.

Required secondary endpoints:

- complete-CDS event recovery;
- exact transcript structure;
- CDS nucleotide, phased CDS segment, exon, and splice-junction precision and
  recall;
- collateral loss of blind source models;
- false claims in assembly-gap negative controls;
- results by transcript count, exon count, WGD context, and chromosome/block;
- wall time, CPU time, peak RSS, disk use, warning count, and failed loci.

Confidence intervals are block/chromosome bootstrap intervals. Paired
comparisons use the same hidden events. A method is not credited for a model
that is already produced in its complete-annotation control.

## Go/no-go rule

The project does not claim a new structure-repair method if its event graph
cannot improve the precision-recall or calibrated-review tradeoff over the
best matched existing baseline on at least two independent plant systems.
Improvements restricted to one random split, one reference, or BUSCO counts do
not pass the gate.

## Acquisition incident

The first official GitHub `v1.2.5` SynGAP checkout and tag-archive download
were incomplete because the server connection ended with OpenSSL
`unexpected EOF while reading`. The incomplete checkout and 7.2 MB archive
are retained under the server project's `third_party/_incomplete` and
`third_party/downloads` directories. No mirror or alternate package source was
substituted. A retry through the same configured Bioconda channel succeeded;
the resulting environment and bundled resources are independently audited in
`docs/SYNGAP_ENVIRONMENT_v1.2.5.md`.

The first LiftOn 1.0.11 installation from the verified PyPI source archive
failed while compiling its `mappy` dependency because `zlib.h` was absent.
The failed stdout, stderr, and resource log are retained. Adding `zlib` to the
same project prefix and retrying the same source archive succeeded; no source
or LiftOn version changed. GeMoMa and LiftOn environment identities are frozen
in `docs/GEMOMA_LIFTON_ENVIRONMENTS_v0.1.md`.
