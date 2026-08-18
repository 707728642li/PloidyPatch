# Initial hawthorn discovery-cohort audit — 2026-08-07

Source examined read-only:

`/nas_data/NFS/Shanzha/shanzha_analysis/cp_genome/lijiajia_genome`

This path was initially presented as the public-genome collection and was
corrected on 2026-08-07. It is retained as the hawthorn discovery cohort. The
actual public-genome pool is inventoried separately in
`PUBLIC_DATA_CATALOG_2026-08-07.md`.

## Inventory summary

- Total directory footprint: approximately 630 GB.
- RNA-seq subtree: approximately 622 GB.
- Raw RNA-seq: 244 GB; cleaned FASTQ: 213 GB; BAM: 164 GB.
- Nine top-level assembly/annotation sets are listed in `genome_list.txt`.
- All nine have assembly FASTA, GFF3, CDS, and protein files.
- Most chromosome-scale assemblies contain 17 sequences; Black assemblies and
  DJX primary are more fragmented.

| Dataset | Sequences | Assembly bp | Longest sequence bp | Genes | Transcripts/proteins | Initial note |
|---|---:|---:|---:|---:|---:|---|
| Black_Hap1 | 1,735 | 762,509,813 | 46,597,164 | 36,469 | 41,755 | fragmented, multi-isoform |
| Black_Hap2 | 386 | 692,329,401 | 45,871,392 | 35,812 | 41,186 | smaller assembly, multi-isoform |
| Black_Primary | 1,699 | 762,936,038 | 46,561,652 | 36,243 | 41,578 | fragmented, multi-isoform |
| Cp_v1.0 | 17 | 779,308,685 | 55,841,475 | 38,131 | 38,131 | chromosome-scale, EVM IDs |
| DJX_Hap1 | 17 | 777,131,881 | 53,747,536 | 39,194 | 39,194 | chromosome-scale |
| DJX_Hap2 | 17 | 797,623,330 | 55,914,249 | 39,893 | 39,893 | chromosome-scale |
| DJX_Primary | 71 | 822,230,708 | 55,920,689 | 41,206 | 41,206 | modest fragmentation |
| JRY | 17 | 809,118,831 | 55,878,510 | 56,407 | 56,407 | strong gene-count inflation candidate |
| RZ | 17 | 760,537,831 | 51,667,533 | 55,457 | 55,457 | strong gene-count inflation candidate |

## Immediate scientific observations

1. JRY and RZ contain roughly 14,000–18,000 more annotated genes than the other
   chromosome-scale sets. This is too large to treat as biological CNV/PAV
   without sequence-level validation.
2. Source annotations are structurally heterogeneous. Cp_v1.0 uses EVM-style
   models; JRY contains both PASA- and EVM-sourced models. Black annotations
   contain multiple isoforms, while most other sets currently expose one
   protein per gene.
3. Assembly fragmentation differs greatly. Apparent gene absence or truncation
   in Black_Hap1/Black_Primary requires explicit contig-edge and gap states.
4. Existing annotation logs inspected so far are eggNOG functional-annotation
   logs, not sufficient provenance for the structural gene models.
5. The available RNA-seq is potentially valuable independent evidence, but
   sample-to-genotype/tissue metadata and mapping-reference consistency must be
   audited before use.

## Suitability

The collection is well suited to the PloidyPatch close-relative/population
mode and to a real annotation-harmonization case study. It is not sufficient
as the sole publication benchmark because it does not by itself provide:

- broad monocot/dicot generalization;
- recent allopolyploid and hexaploid stress tests;
- curated truth for every error class;
- guaranteed independent long-read transcript evidence;
- complete public accession/version/provenance metadata.

## Required next checks

- obtain species/cultivar, ploidy, assembly technology, accession, release, and
  structural-annotation pipeline metadata for every set;
- validate GFF3 parent-child relationships and FASTA sequence-name matching;
- translate CDS and quantify internal stops, missing start/stop, and phase
  errors;
- quantify single-exon, short-protein, TE-overlap, and contig-edge models;
- calculate assembly N50/gap content and gene-space completeness;
- audit RNA sample metadata, BAM headers, mapping references, and junction
  coverage;
- identify truly independent Iso-Seq/full-length cDNA data;
- preserve all source files read-only and reference them through checksummed
  manifests rather than copying 630 GB.
