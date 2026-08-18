# Brassica allopolyploid evidence checkpoint — 2026-08-07

## Scope

This checkpoint establishes a QC-passed, chromosome-scale AC allopolyploid
benchmark source and a leakage-controlled external synteny stratum. It does
not evaluate a repair method and is not evidence that PloidyPatch outperforms
an existing tool.

## Dataset decision and immutable staging

The old Ensembl Plants `B. napus` assembly `AST_PRJEB5043_v1` is retained only
as a large, fragmented format sentinel. Its thousands of sequence IDs make it
unsuitable as the core chromosome-scale WGD validation source.

The formal validation source is NCBI `GCF_020379485.1`, cultivar Da-Ae,
Annotation Release 102. Five source files were copied from
`/nas_data/NFS/Public_genome_data/` into the isolated server project and
verified against their source SHA-256 values:

| Role | Bytes | SHA-256 |
|---|---:|---|
| genome | 1,014,840,906 | `67c4dc1c246e021a55e8477eaebb27e07de835e1ccb08ffbcb83f72346f34727` |
| FAI | 117,553 | `fb3ee2c1b712482afd24299ccefb494d6b68f503d5a9c9d47bc3f1681393b316` |
| GFF3 | 548,371,256 | `1816b53f16cd8b07a9108c8ddac2b90e9ca2197613e0c28d1782d3f6b50596ca` |
| protein | 54,635,253 | `9450f2ccb26668374704fcf3d5341bc65ee9af337b467dd2e8f5f95824daf54d` |
| CDS | 159,339,821 | `106e9b0b59e98d1a6eff0e6123b195ca66fa8b9f73f6af56f41e3f67b14ebdb1` |

Staging manifest:
`/data/codexli/projects/PloidyPatch/logs/staging/refseq_bna_daae_gcf0203794851_20260806T174412Z.tsv`.

## Raw-bundle diagnosis and normalization rule

The unmodified full bundle fails the strict gate because of five duplicate
gene/transcript IDs, two invalid and three out-of-bounds organellar features,
and 328 proteins containing the illegal character `.`. The five coordinate/ID
problems occur on the chloroplast record. On the 19 nuclear chromosomes, 292
translation records contain `.` and their corresponding CDS sequences contain
an internal stop codon.

PloidyPatch does not replace `.` with an amino acid or silently accept those
models. `ploidypatch normalize ncbi-primary` applies a declared transformation:

1. select NCBI `region` features with `genome=chromosome`;
2. validate protein/CDS characters, internal stops, and translated length;
3. remove each invalid transcript and its descendant features as one hierarchy;
4. remove a gene only when it has no retained transcript;
5. subset protein and CDS FASTA by the same retained transcript relation;
6. rewrite only the 19 selected genome records and generate a matching FAI;
7. emit `excluded_translations.tsv`, checksums, counts, and reason codes.

The run excluded 292 translation transcripts spanning 288 genes. Every one had
both `invalid_protein_character` and `cds_internal_stop`; no residue was
invented. Runtime was 89.69 s with 4,334,236 KiB peak RSS. The normalized
manifest is at:
`/data/codexli/projects/PloidyPatch/data/derived/normalized_bundles/v0.1/bna_daae_primary/manifest.json`.

## Independent normalized-bundle QC

An independent `ploidypatch audit --checksums` run produced grade `pass`:

| Measurement | Result |
|---|---:|
| Primary chromosomes | 19 |
| Assembly span | 816,166,126 bp |
| Genes | 99,791 |
| Transcripts | 118,340 |
| Matched protein/CDS translations | 113,154 |
| Translation-length mismatches | 0 |
| Invalid protein/CDS records | 0 / 0 |
| Internal-stop protein/CDS records | 0 / 0 |
| Missing or extra protein/CDS relations | 0 |
| GFF errors | 0 |
| Quality-gate warnings | 0 |

Ninety CDS lengths have a frame remainder but all 90 are protein-length
consistent; they are recorded as observations, not silently changed. The audit
report SHA-256 is
`0cb2c10d3695c36f510dc6875f5ff1578e7d6ecfd036b841ab5e0820eac1de1b`.

## External WGDI evidence

The query uses one valid representative protein per Da-Ae gene. The two
external references are `B. rapa` (A lineage) and `B. oleracea` (C lineage).
DIAMOND used `--more-sensitive`, E-value `1e-5`, at most 20 targets, and 48
threads per comparison. WGDI used `multiple=2`, `process=32`, `score=100`,
`evalue=1e-5`, `grading=50,40,25`, `mg=40,40`, `pvalue=0.2`, and `repeat=20`.

| Comparison | DIAMOND hits | WGDI blocks | WGDI gene pairs | Collinearity SHA-256 |
|---|---:|---:|---:|---|
| Da-Ae × B. rapa | 954,707 | 4,001 | 144,935 | `fb30ee1a070b90671cc6ee2b3ffbe35aa77c4bae3d9dde509a9f91a340d3e753` |
| Da-Ae × B. oleracea | 1,107,975 | 4,027 | 146,683 | `2f61fdb32684b5862860e2a78988f5adda75b48537560a593fd744b1221e0ed9` |

The per-gene evidence table has 85,091 rows and SHA-256
`76a8621f6359e22398e2944ed7af12f58ff89413149347a99e39fab2f7cf4072`.
It is evaluator-only because it participates in benchmark stratification.

## Formal candidate catalog

The normalized GFF3 yields 85,091 conservative
`annotation_missing_gene` candidates. All are structurally unique under the
strict identifier-independent signature and all join exactly once to external
WGDI evidence; there are no missing or extra gene IDs.

| Stratum | A subgenome | C subgenome | Total |
|---|---:|---:|---:|
| expected and cross-progenitor blocks | 31,755 | 31,854 | 63,609 |
| expected-progenitor block only | 1,752 | 3,677 | 5,429 |
| cross-progenitor block only | 1,310 | 1,212 | 2,522 |
| no retained block | 4,329 | 9,202 | 13,531 |
| Total | 39,146 | 45,945 | 85,091 |

The catalog SHA-256 is
`56414fb0421ca27db1890f2a8854980e67e755c2d9eba8f02dc72c371e8235fb`.

Subgenome identity is assigned from the assembly chromosome labels (`A1–A10`,
`C1–C9`), not inferred from which reference receives the stronger hit. The
large `expected_and_cross` class is biologically plausible because both
Brassica progenitor lineages share older mesohexaploidy. Consequently, WGDI is
used as evidence availability/difficulty, not as a naive subgenome classifier.

## Reproducibility and next gate

- Code snapshot: `582e99b1e7ef68b35bdb0864056c8d7c8d64e861`
- Snapshot archive SHA-256:
  `d95e721c6a3b00117b797eb54bd087eefdd4b725ec4c38a44d44c885b2674907`
- Server project: `/data/codexli/projects/PloidyPatch`
- Local ignored mirrors of the small manifests:
  `results/brassica_wgd_v0.1/`

The next scientific gate is not a large PloidyPatch model. It is a frozen,
stratified perturbation sample plus a simple sequence/projection baseline,
negative controls, and partial-credit structure metrics. This tests whether
the external plant-WGD evidence adds held-out value beyond a sequential
baseline before investing in the joint event graph.
