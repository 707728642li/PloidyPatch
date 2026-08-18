# Actinidia v0.5 syntax-only preflight

Date: 2026-08-08  
Scope: public reference metadata and file syntax only, before protocol/execution freeze and before any pair, event, holdout, candidate, label, score, or endpoint enumeration.

No WGDI, DIAMOND, homology search, copy-pair inference, candidate generation, or truth/event enumeration was run for this audit. Counts below are annotation syntax counts, not pair or candidate counts.

## Frozen-source candidates inspected

All 20 gzip members passed `gzip -t`. Bytes are compressed-file bytes and hashes are SHA-256 of the compressed bytes.

| Species | Role | File | Bytes | SHA-256 |
|---|---|---|---:|---|
| Red5 | genome | `/nas_data/NFS/Public_genome_data/plantgarden_database/raw_data/t3625/Actinidia_chinensis.Red5_PS1_1.69.0.dna.toplevel.fa.gz` | 165032199 | `b8a50d53bf9b18fa879da3ade5e8162e4d896740afddf0fab1967e23b9c11533` |
| Red5 | GFF3 | `/nas_data/NFS/Public_genome_data/plantgarden_database/raw_data/t3625/Actinidia_chinensis.Red5_PS1_1.69.0.45.chr.gff3.gz` | 7400242 | `41a81f0a85256d58b4fd2fde4032b6326c8295592d4ea8f0341d37575b2d8ca9` |
| Red5 | protein | `/nas_data/NFS/Public_genome_data/plantgarden_database/raw_data/t3625/Actinidia_chinensis.Red5_PS1_1.69.0.pep.all.fa.gz` | 10391929 | `3256699bebf59ef25c6811008d1bec44c95dfb85e3d7de4abac70264cd98663d` |
| Red5 | CDS | `/nas_data/NFS/Public_genome_data/plantgarden_database/raw_data/t3625/Actinidia_chinensis.Red5_PS1_1.69.0.cds.all.fa.gz` | 15469291 | `1f3d577adef843ed72f45ae3da4e66a4751d8ccb7367d3eae9556bb8b4ad7822` |
| A. eriantha | genome | `/nas_data/NFS/Public_genome_data/plantgarden_database/raw_data/t165200/Actinidia_eriantha.chr.fasta.gz` | 199769654 | `259a69329f7a2e07d45a1fea2918b4eba5f3f33b639cc434c93586615fbc5024` |
| A. eriantha | GFF3 | `/nas_data/NFS/Public_genome_data/plantgarden_database/raw_data/t165200/Actinidia_eriantha.gff3.gz` | 3994059 | `b26fca387536811110270d2029f7ceb851723fa15a7a380860216d22d8a066c5` |
| A. eriantha | protein | `/nas_data/NFS/Public_genome_data/plantgarden_database/raw_data/t165200/Actinidia_eriantha.pep.fasta.gz` | 8853995 | `45d89c5d022d7b2ee3f9992010d3bd50b4bcbbee323b35d8c28088ac10c28f73` |
| A. eriantha | CDS | `/nas_data/NFS/Public_genome_data/plantgarden_database/raw_data/t165200/Actinidia_eriantha.cds.fasta.gz` | 13220138 | `162e2340ff3044d558374d4943ee88679dcc137c9c839ea311e408c3c4d1cb87` |
| A. rufa | genome | `/nas_data/NFS/Public_genome_data/plantgarden_database/raw_data/t165716/ARU_r1.0.pmol.fasta.gz` | 187899618 | `a275b3f5f2e3d92ba822906fd8003d205fb6489a58e2ee48f307fd0fd0cd97f0` |
| A. rufa | GFF | `/nas_data/NFS/Public_genome_data/plantgarden_database/raw_data/t165716/ARU1.0.genes.gff.gz` | 7072736 | `a1dc19871b24fb64b44c4ff2a8ecf95daf96c24f5505cc1cc86442d6752b26fe` |
| A. rufa | protein | `/nas_data/NFS/Public_genome_data/plantgarden_database/raw_data/t165716/ARU1.0.proteins.fasta.gz` | 13719497 | `050d96e551c8b0b2b800132364851ab1fa29d5046236237ff13d1884f9291cef` |
| A. rufa | CDS | `/nas_data/NFS/Public_genome_data/plantgarden_database/raw_data/t165716/ARU1.0.cds.fasta.gz` | 21826341 | `4d24c623eab6cef480330126cc01f4cff53afeeb5937ed73c1fb88f9d7d7d21f` |
| R. simsii | genome | `/nas_data/NFS/Public_genome_data/plantgarden_database/raw_data/t118357/GCA_014282245.1_ASM1428224v1_genomic.fna.gz` | 175086237 | `0101fc0d9d7e27a9cc59687cb61cf7014e66372f40558a5b5ab4935aba640aab` |
| R. simsii | GFF3 | `/nas_data/NFS/Public_genome_data/plantgarden_database/raw_data/t118357/GCA_014282245.1_ASM1428224v1_genomic.gff.gz` | 7027134 | `681c5bbb24f58260318315bcf3094b2a20c2d020d2ca3b7adbd730647ebbffff` |
| R. simsii | protein | `/nas_data/NFS/Public_genome_data/plantgarden_database/raw_data/t118357/GCA_014282245.1_ASM1428224v1_protein.faa.gz` | 8870740 | `ff714ad8d26e9e2294d74c193c48be136a1f23469881e06fe515856ef9cb5aef` |
| R. simsii | CDS | `/nas_data/NFS/Public_genome_data/plantgarden_database/raw_data/t118357/GCA_014282245.1_ASM1428224v1_cds_from_genomic.fna.gz` | 15381559 | `5549e057a026c3043880a93e4eb9b54e57fad91a398ff95892a7bd8bfc50f8d0` |
| D. oleifera | genome | `/nas_data/NFS/Public_genome_data/plantgarden_database/raw_data/t227308/Dol.genome.fa.gz` | 220674718 | `f6c6c409198c3dd5854ea7862f37e7ea8644657972dc3eaf642cede4cb7ca200` |
| D. oleifera | GFF3 | `/nas_data/NFS/Public_genome_data/plantgarden_database/raw_data/t227308/Dol.gff3.gz` | 4391416 | `6820b46cd41e098962416d35b98ebbfb7922a36475287e10d9289a0d96f5de68` |
| D. oleifera | protein | `/nas_data/NFS/Public_genome_data/plantgarden_database/raw_data/t227308/Dol.pep.fa.gz` | 6597993 | `6111890b204e9575d83874c4c3bc4a9947c91dd622970ef8af9297fbf7d5adb2` |
| D. oleifera | CDS | `/nas_data/NFS/Public_genome_data/plantgarden_database/raw_data/t227308/Dol.cds.fa.gz` | 9955838 | `96cf7501b2f5f0bef05e172d40042bab92927504c7ea9770282599e03d555c71` |

## Primary-seqid contract

The exact primary sets are:

- Red5: `LG1` through `LG29` (29 exact IDs).
- A. eriantha: `Chr01` through `Chr29` (29 exact IDs). `Chr00` is an unplaced bin and is not primary.
- A. rufa: `ARU1.0ch01` through `ARU1.0ch29` (29 exact IDs).
- R. simsii: `CM024953.1`, `CM024954.1`, `CM024955.1`, `CM024956.1`, `CM024957.1`, `CM024958.1`, `CM024959.1`, `CM024960.1`, `CM024961.1`, `CM024962.1`, `CM024963.1`, `CM024964.1`, `CM024965.1` (13 exact IDs).
- D. oleifera: `Chr1` through `Chr15` (15 exact IDs).

Literal GFF-to-genome seqid coverage was complete in every inspected source: Red5 29/29 annotated seqids, A. eriantha 30/30 including `Chr00`, A. rufa 29/29, R. simsii 552/552, and D. oleifera 368/368. Red5's chromosome-only GFF intentionally leaves 1205 top-level FASTA sequences unused. D. oleifera's GFF uses 368 of 2533 genome records. R. simsii contains the 13 chromosomes plus 539 non-primary records.

## Annotation and exact-ID mapping audit

| Species | Genes | Transcripts | Protein FASTA records | Exact gene-to-protein coverage | Recommended frozen mapping |
|---|---:|---:|---:|---|---|
| Red5 | 32950 | 33021 | 33115 all-sequence proteins | 32950/32950 genes linked; 33021/33021 chromosome transcripts/proteins; 94 all-pep proteins absent from chromosome GFF | Whitelist peptide first tokens by chromosome-GFF CDS `protein_id`; map peptide header `transcript:` to GFF mRNA `transcript_id`/ID and header `gene:` to GFF gene `gene_id`/ID. Strip only the literal GFF type namespace when required. |
| A. eriantha | 42988 | 42988 | 42988 | 42988/42988 exact | Protein/CDS first token equals both gene and mRNA ID. Treat the shared gene/mRNA identifier as a provider-specific alias, not as globally unique GFF3 ID. |
| A. rufa | 63947 | 0 explicit mRNA | 63947 | 63947/63947 exact | Protein/CDS first token equals gene ID; children point directly to gene. A normalized synthetic transcript layer is required for transcript-oriented tools. |
| R. simsii | 33466 | 32999 | 32999 | 32999 protein-coding genes linked one-to-one; 467 genes have no protein mRNA | Protein first token equals CDS `protein_id`/`Name`; follow CDS `Parent` to mRNA and mRNA `Parent` to gene. The CDS FASTA first token is not the protein ID; use its `[protein_id=...]` tag. |
| D. oleifera | 30530 | 30530 | 30530 | 30530/30530 one-to-one | Protein/CDS first token equals mRNA ID; follow mRNA `Parent` to the gene ID. |

Red5 has 71 genes with multiple proteins/transcripts, but every chromosome protein resolves to exactly one gene. R. simsii and D. oleifera have zero protein-to-multiple-gene aliases and zero multi-protein genes in the supplied annotation. All protein and CDS FASTA first-token IDs are unique and nonempty.

## GFF feature counts and parser risks

| Species | Feature counts | Required handling |
|---|---|---|
| Red5 | chromosome 29; gene 32950; mRNA 33021; CDS 180837; exon 180837; five_prime_UTR 32805; three_prime_UTR 32682; biological_region 10846 | CDS IDs are repeated across segments for 26546 multipart CDS values. Filter the all-peptide FASTA by the 33021 chromosome-GFF `protein_id` whitelist before WGDI preparation. |
| A. eriantha | gene 42988; mRNA 42988; CDS 217784; five_prime_UTR 34391; three_prime_UTR 35796 | Gene and mRNA share all 42988 ID values. Multipart CDS and UTR IDs also repeat. Use provider normalization before any strict GFF3 graph parser. Exclude `Chr00` from primary analyses. |
| A. rufa | gene 63947; CDS 297002; exon 297078 | No mRNA features. Child `Parent` values name the gene directly. |
| R. simsii | region 552; gene 33466; mRNA 32999; CDS 178107; exon 185665; rRNA 54; tRNA 413 | CDS IDs repeat across segments for 28031 multipart CDS values. Use exact NCBI `protein_id` aliases. |
| D. oleifera | gene 30530; mRNA 30530; CDS 140908; exon 146752; five_prime_UTR 16567; three_prime_UTR 16365 | CDS IDs repeat across segments for 23040 multipart CDS values. |

Across all five GFF files: malformed-column rows = 0, invalid coordinates = 0, invalid strand = 0, invalid phase = 0, and unresolved `Parent` values = 0. Repeated multipart feature IDs and A. eriantha's cross-feature shared IDs mean a global "ID occurs once" assertion would be incorrect for these provider files; normalization must preserve segment grouping while creating unambiguous internal feature identities.

## Formal-use boundary

This report establishes only that the predeclared source files and deterministic identifier mappings are suitable inputs to a later frozen execution. It does not establish Ad-alpha membership, evaluator 1:2 relationships, self-WGDI pairs, event counts, candidate counts, truth labels, or tool performance. Those remain unopened until both the Actinidia protocol freeze and execution-implementation freeze bind the exact source bytes, mapping rules, code, environments, and model artifact.
