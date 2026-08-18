# Glycine matched-baseline inputs v0.1

Date materialized: 2026-08-07  
Server root: `/data/codexli/projects/PloidyPatch`  
Status: inputs validated; candidate methods not yet run

## Contract

All three bundles were produced with `ploidypatch normalize
primary-annotation`. The operation retains every GFF record on chromosomes
declared by a GFF chromosome feature, an NCBI `genome=chromosome` region, or
an explicit `Chr*` region alias. It performs no protein/CDS translation
filtering. No sequence was selected by length or accession heuristics.

The bundles are under
`data/derived/syngap_glycine_v0.1/{gma_complete,gma_blind,gso_reference}`.
Each directory contains an uncompressed GFF3, genome FASTA, FAI, and manifest.

## Frozen identities

| Bundle | Chromosomes | Genome bp | GFF records | Genes | GFF SHA-256 | Genome SHA-256 |
|---|---:|---:|---:|---:|---|---|
| `gma_complete` | 20 | 949,183,385 | 1,318,728 | 55,589 | `dfb99e2564f314be2e78d6e85477f1779478b44326994dbc8c17270998ee0ceb` | `628e61de1b126631a224956ea233e51ddfd1d8fe4189cbc56e8494b305677c02` |
| `gma_blind` | 20 | 949,183,385 | 1,291,815 | 54,789 | `3858858a686f61eab69f6783781931b22863987ad16a333b1d6f571a00a0fb05` | `628e61de1b126631a224956ea233e51ddfd1d8fe4189cbc56e8494b305677c02` |
| `gso_reference` | 20 | 969,795,186 | 1,272,228 | 46,722 | `8ec0641a924f85caf9249de5e400b35d295e7c9f17fad937243387d1058ce260` | `233306df1fd930121c4326defb851ff3821b776a4ffb10758baa9bfbffebe6e7` |

The complete and blind target FASTA and FAI hashes are identical. Their gene
count differs by exactly the 800 frozen annotation-missing-gene events. The
reference chromosome mapping was recovered from paired aliases such as
`Alias=...,1,Chr01,...`; `Super-Scaffold_*`, `MT`, and `Pltd` regions do not
meet that explicit rule and were excluded.

Source hashes are frozen inside each manifest. The target source genome is
`62a2c3f917f099a504752796d8c709ad3c8466da518874bd282616eb554e5113`;
the blind source GFF is
`fc8c83440805a2a7d0babf46915648ee5ed4dc1606f943bcf5ccccbadc5f28f8`.

## Resource observations

The three materializations completed in 32.60, 39.35, and 36.75 seconds.
Peak RSS was 1.94, 1.97, and 2.04 GB, respectively. Complete timing and stderr
logs are frozen under `logs/baseline/input_normalization_glycine_v0.1`.

## SynGAP input-contract preflight

`scripts/preflight_syngap_annotation.sh` reproduces the upstream preprocessing
commands in isolation: JCVI GFF-to-BED with `--type=mRNA --key=ID
--parent_key=Parent --primary_only`, gffread CDS/protein extraction, and seqkit
selection by the JCVI primary IDs. The final frozen results are under
`results/baselines/syngap_v1.2.5/preflight`.

| Bundle | JCVI primary IDs | Primary CDS | Primary protein | Missing CDS/protein |
|---|---:|---:|---:|---:|
| `gma_blind` | 54,789 | 54,789 | 54,789 | 0 / 0 |
| `gma_complete` | 55,589 | 55,589 | 55,589 | 0 / 0 |
| `gso_reference` | 46,722 | 46,722 | 46,722 | 0 / 0 |

JCVI reported 31,552, 32,388, and 22,007 non-primary isoforms skipped,
respectively; this is the declared upstream `--primary_only` behavior, not a
PloidyPatch filter. Both gffread operations had empty stderr for all three
bundles. JCVI, gffread, and seqkit all exited zero, and every selected ID had
one CDS and protein FASTA record. The three parallel preflights completed in
65--68 seconds; individual gffread steps peaked below 274 MB RSS.

## Materialization incident

After successful bundle creation, a read-only validation command was left
running when its SSH wrapper timed out. PowerShell-to-SSH quote removal changed
the intended FASTA-header pattern `'^>'` into a shell redirection. The orphaned
`grep` process therefore truncated only
`gma_complete/primary_chromosomes.genome.fa` to zero bytes at 08:09.

The exact process group was identified through `fuser`, `lsof`, and `ps` and
terminated before it advanced to another bundle. The zero-byte file is retained
at `third_party/_incomplete/gma_complete_genome_zero_20260807_0809.fa`. Because
the complete and blind bundles deliberately use the same source genome and had
already frozen the same output hash, the verified blind FASTA was copied to a
temporary path, rehashed, and atomically installed as the complete FASTA. Both
files now reproduce the manifest hash
`628e61de1b126631a224956ea233e51ddfd1d8fe4189cbc56e8494b305677c02`.
No baseline candidate run had begun, so this incident cannot affect predictions
or hidden-truth scoring.
