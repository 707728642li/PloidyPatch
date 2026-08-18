# Public plant-genome catalog — 2026-08-07

Authoritative source examined read-only:

`/nas_data/NFS/Public_genome_data/`

## Available collections

The directory contains Ensembl Plants, a normalized `genome_database`, GDR,
Phytozome, PLAZA, PlantGIR, species metadata, and other comparative-genomics
resources. The initial benchmark should use one internally consistent source
where possible; mixing providers is reserved for explicit annotation-method
robustness experiments.

The Ensembl-derived collection inspected here contains release 62 GFF3 files
and matched genome, CDS, cDNA, and peptide FASTA files with checksum manifests.

## Staging status

The ten priority-1 bundles in the pilot, population, and Brassica WGD tiers
were staged under `/data/codexli/projects/PloidyPatch/data/public/` on
2026-08-07. Sixty files were copied and source/destination SHA-256 values were
verified. All ten bundles subsequently passed PloidyPatch input gate v0.3; see
`INPUT_QC_2026-08-07.md`.

Glycine, maize, and wheat remain deliberately unstaged because they are held-out
or late stress-test tiers and must not enter early parameter fitting.

## Recommended benchmark ladder

| Stage | Genome set | Purpose | Role |
|---|---|---|---|
| Pilot | *Arabidopsis thaliana* TAIR10, *A. lyrata* v1.0, *A. halleri* Ahal2.2 | small, curated/near-curated dicot clade; parser and perturbation development | development |
| Population | *Oryza sativa* IRGSP-1.0, IR64, MH63, ZS97 | same-species multi-assembly annotation consistency and PAV | tuning/validation split by accession |
| WGD/allopolyploid | *Brassica rapa* A, *B. oleracea* C, *B. napus* AC | subgenome-aware homeology and copy-loss stress test | tuning then frozen test strata |
| Paleopolyploid | *Glycine max* v2.1 and *G. soja* ASM419377v2 | independent WGD/PAV validation | held-out test |
| TE-rich crop | maize B73 NAM5 plus W22/PH207 to be located or downloaded | split/fusion, tandem family, and TE-associated errors | held-out test |
| Hexaploid stress | wheat IWGSC RefSeq v2.1 with A/D progenitors and tetraploid relatives | scalability and A/B/D homeolog assignment | final stress test |
| Optional external cohort | user-supplied assemblies, only when intentionally designated | engineering regression or later biological case study | excluded from the primary benchmark by default |

## Selection rationale

- The pilot is deliberately small enough for rapid full-pipeline iteration.
- Rice provides many local assemblies under one provider and separates
  population behavior from dicot-specific assumptions.
- Brassica directly tests the primary claimed innovation: explicit
  subgenome/WGD-aware repair.
- Glycine remains outside early parameter fitting so that the WGD claim has an
  independent test.
- Wheat is postponed until correctness and memory scaling are demonstrated;
  beginning with it would slow methodological iteration and complicate error
  diagnosis.

## Copy policy

Only matched unmasked top-level genome FASTA, whole-genome GFF3, all CDS, all
peptides, README, and checksum files should be copied into the project. Soft-
masked and per-chromosome duplicates are not copied unless required by a
specific experiment. Each copied file must retain its provider filename and be
recorded in a manifest with source path, release, size, checksum, role, and
copy date.

## Metadata still required

- assembly accessions and provider release provenance;
- known ploidy and named WGD/subgenome events;
- canonical chromosome mappings;
- structural-annotation method and evidence sources;
- independent curated or long-read truth suitable for each test;
- license/citation obligations.
