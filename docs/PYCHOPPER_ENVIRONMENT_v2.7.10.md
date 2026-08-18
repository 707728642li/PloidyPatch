# Pychopper 2.7.10 environment audit

Audit date: 2026-08-08.  Server prefix:
`/data/codexli/projects/PloidyPatch/envs/ploidypatch-pychopper`.

## Purpose

This environment is restricted to canonical full-length cDNA preprocessing for
the Golden Delicious natural-data validation.  It is separate from the
PloidyPatch model, development, projection, and synteny environments.

## Resolved core software

| Package | Version |
|---|---:|
| Python | 3.13.14 |
| Pychopper | 2.7.10 |
| HMMER | 3.4 |
| NumPy | 2.5.1 |
| pandas | 2.3.3 |
| pytz | 2026.3.post1 |
| pysam | 0.24.0 |
| parasail-python | 1.3.4 |
| python-edlib | 1.3.9.post1 |
| Biopython | 1.88 |
| tqdm | 4.70.0 |

The formal result captures both `conda list --explicit` and the complete conda
JSON package inventory.  It also captures `pip freeze --all`; those run-local
locks, rather than this human-readable table, are authoritative.

## Compatibility finding

Pychopper 2.7.10 completes read classification under pandas 3, but its PDF
reporting code calls `float()` on a one-row pandas Series and fails under
pandas 3.0.5.  The validation runtime therefore requires pandas 2.x.  The
installed conda-forge pandas 2.3.3 artifact is retained at
`work/software_packages/pandas-2.3.3-py313h08cd8bf_2.conda` with SHA-256
`b998c30e7ff13fc966220891dc0a8318b0a6730933280d76ffa5be46ff928af5`.
Both `pytz` and the Python `tzdata` distribution are installed so pandas has a
complete runtime and metadata dependency set.

Pychopper's Python package metadata pins `tqdm==4.26.0`, whereas the maintained
Bioconda recipe resolves a modern tqdm.  The environment intentionally follows
the curated Bioconda runtime (`tqdm 4.70.0`) rather than forcing a 2018 tqdm
release into Python 3.13.  `pip check` is therefore expected to report exactly
this upstream metadata mismatch and no missing runtime dependency.

## Determinism audit

Pychopper's cutoff autotuning samples reads with NumPy without exposing a CLI
seed.  `scripts/run_pychopper_seeded.py` sets seed `20261005` before invoking
the unmodified Pychopper main function.  The launcher additionally fixes
`PYTHONHASHSEED=0`, `SOURCE_DATE_EPOCH=0`, BLAS thread counts, the HMMER `PATH`,
and a run-private Matplotlib configuration directory.

Two clean executions on the same 1,000-read real Golden Delicious subset gave:

| Artifact | SHA-256 in both runs |
|---|---|
| oriented primary FLNC FASTQ.gz | `900a4ecd660d8e89e3e55cab1b6cbf2cdc28aacb5c38387405a97f5a6ecf1833` |
| Pychopper statistics TSV | `0268f3ed2f7dfd707e83ff0cd2b683c1d958178617a5311a6a24b005557047b6` |
| Pychopper report PDF | `3d7605b4f48f21d1694272ecd66ff2c8767793e10f8b902ac40e48201e6e0afc` |

The subset contained 843 reads classified with one primer-pair configuration,
but 842 materialized FLNC records because one trimmed fragment failed the
50-nt output-length requirement.  The formal summarizer reports both counts
and uses the materialized FASTQ count for downstream denominators.
