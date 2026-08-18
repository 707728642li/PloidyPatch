# Controlled benchmark smoke validation — 2026-08-07

## Scope

This report validates benchmark mechanics only. It does not measure a repair
method and must not be cited as biological performance. Two sentinels were
used for the scorer:

- **no-op**: return the perturbed GFF3 unchanged; expected recall is zero;
- **oracle**: return the immutable source annotation; expected strict scores
  are one with no collateral changes.

Generation used the v0.1 `annotation_missing_gene` primitive and seed 17.
Scoring used code snapshot `8e87821`; 15/15 tests passed on the server.

## Datasets and generation

| Field | Arabidopsis thaliana TAIR10 | Brassica napus AST_PRJEB5043_v1 |
|---|---:|---:|
| Ploidy regime | diploid development sentinel | allopolyploid AC sentinel |
| Source compressed SHA-256 | `9f2701d7151b4467941b5dae4171b8955ea4260c4d78f41be4fb3ff44782d2fd` | `13944fe24c2dd4375c39f9f14cf6034a9d14ea344423d606524099d84a7751a6` |
| Decompressed GFF3 lines | 824,417 | 1,565,491 |
| Conservative eligible genes | 27,628 | 101,040 |
| Generated events | 100 | 100 |
| Removed records | 3,017 | 1,224 |
| Distinct hidden transcript structures | 180 | 100 |
| Perturbed GFF3 SHA-256 | `dcd34fa23fd088ccac60413f10edbd116c91cb13e1d489202ceed13e04f028a3` | `b8e65cff666cd73f691bfb3f576d00a821c4115d267d21b84f2cbb7b5d4f1026` |
| Hidden truth SHA-256 | `fcfdab6faba0ebceaed2d95d2dbc5645f356193364bab0b25b19c91bc96520b8` | `0fd0583d56d675e9d88110e264146c9e74a0a0fe301c7ed761895841eebcfecc` |
| Generation wall time | 10.91 s | 20.91 s |
| Generation peak RSS | 1,299,936 KiB | 2,395,844 KiB |

The TAIR10 run was independently regenerated after a performance-only hashing
change. `perturbed.gff3`, `hidden_truth.json`, and `manifest.json` matched the
first run byte for byte. The redundant reproduction files were then removed.

## Reversibility

Both perturbed annotations were restored using evaluator-only truth. The
restored decompressed text hashes exactly matched the original inputs:

- TAIR10: `dd1e304a0b3fccb95f5f95ac990b4746ca751ab95e90960f4b40d1d029c96440`;
- B. napus: `ea6a2d95b785b5ab1f830b2f3c4795da5d23bce59219ec4a7e206395cc8221e3`.

Restoration took 1.50 s / 413,696 KiB for TAIR10 and 2.11 s / 704,512
KiB for B. napus. The large restored GFF3 files were removed after verification
because they are exactly reproducible; restoration reports were retained.

## Scorer sentinels

| Dataset | Candidate | Strict TP | FP | FN | Precision | Recall | Exact events | Collateral losses |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| TAIR10 | no-op | 0 | 0 | 180 | undefined | 0.0 | 0/100 | 0 |
| TAIR10 | oracle | 180 | 0 | 0 | 1.0 | 1.0 | 100/100 | 0 |
| B. napus | no-op | 0 | 0 | 100 | undefined | 0.0 | 0/100 | 0 |
| B. napus | oracle | 100 | 0 | 0 | 1.0 | 1.0 | 100/100 | 0 |

Both runs passed the source/truth checksum gate and the structural-ambiguity
gate: no hidden transcript signature remained in the perturbed annotation.

The Ensembl `AST_PRJEB5043_v1` assembly is highly fragmented (thousands of
sequence IDs). It remains a useful large-file and provider-format sentinel but
was demoted from biological WGD validation before model fitting. The
chromosome-scale Da-Ae `GCF_020379485.1` bundle is the formal AC allopolyploid
validation source; see `BRASSICA_WGD_EVIDENCE_2026-08-07.md`.

## Artifact locations

Server smoke runs:

- `/data/codexli/projects/PloidyPatch/benchmark/smoke/ath_tair10_annotation_missing_gene_seed17_v0.1`
- `/data/codexli/projects/PloidyPatch/benchmark/smoke/bna_ac_annotation_missing_gene_seed17_v0.1`

The retained server directories contain the perturbed GFF3, hidden truth,
manifest, restoration report, and no-op/oracle score reports. Small JSON
reports are mirrored under local ignored `results/benchmark_smoke/`.

## Interpretation and next gate

The file contracts, deterministic selection, reversibility, strict scoring,
and basic diploid/allopolyploid GFF3 compatibility are working. This does not
yet establish realistic error distributions or plant-specific superiority.
The next gate is a non-leaking stratified candidate catalog followed by a
simple projection baseline and false-positive challenge cases.
