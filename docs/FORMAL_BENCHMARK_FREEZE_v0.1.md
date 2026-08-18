# Formal Da-Ae benchmark freeze v0.1 — 2026-08-07

## Claim boundary

This freeze validates target sampling, data separation, reversibility, and
metric sentinels. It contains no method result and makes no superiority claim.

## Frozen source and sampling design

- Source normalized GFF3 SHA-256:
  `4cd8fb0244ef51d7055d5866f11086c1fee927f5cfbb41509f118f6fb0f29dad`
- Candidate catalog: 85,091 rows, SHA-256
  `56414fb0421ca27db1890f2a8854980e67e755c2d9eba8f02dc72c371e8235fb`
- Sampling plan: `config/bna_daae_sampling_plan_v0.1.tsv`, SHA-256
  `54c109e50f191cff2cbf7076bfd53554bad44491ee71541ed832e2978b72114a`
- Sampling seed: `20260807`
- Deterministic algorithm: `sha256_candidate_rank_v1`

The design crosses two assembly-defined subgenomes, four WGDI evidence strata,
and three transcript-count bins. There are 24 exact strata. Each
subgenome × WGDI cell contributes 60 one-transcript, 25 two-transcript, and 15
three-or-more-transcript genes, for 100 events per cell and 800 total. This is
an intentionally balanced stress-test sample, not a prevalence estimate.

Every requested cell passed its availability gate. The smallest cell contained
78 eligible candidates for a requested 15. The evaluator-only sample has 800
rows and SHA-256
`cdbd2b2a889ea8e2f84a3cdca4bd99ce4ea584893058550d4541fc19281222b8`.

## Separation and immutable artifacts

Server root:
`/data/codexli/projects/PloidyPatch/benchmark/formal/v0.1/bna_daae_annotation_missing_gene_seed20260807/`.

- `blind/` contains only `perturbed.gff3` and the non-identifying public
  manifest.
- `evaluator/selection/` contains target identities and its custody manifest.
- `evaluator/truth/` contains exact removed records.
- `evaluator/validation/` contains only aggregate/restoration validation
  reports.

The evaluator directory is mode `0700`. Candidate runs must additionally use a
container or account that mounts only `blind/` and declared public reference
evidence; same-owner filesystem permissions alone are not considered a final
leakage barrier.

| Artifact | SHA-256 |
|---|---|
| blind `perturbed.gff3` | `2a4085f95a3d10e9a54b075e478e29bce86a583e81ad95c74dcb2ee96e0c9bbe` |
| evaluator `hidden_truth.json` | `28e3278b1bc2587625f9b1747d1d235c927034d5bd0dfa476f9bcd12f735e371` |
| public `manifest.json` | `74433a018b970338858b7d03ff6e7ee4681fb77b9f5c6dce3d9bf7f9d33a0a58` |

Generation used code `8bfeedfd9a463de8412e46dfbc8f4c863dae0921`.
Its archive SHA-256 is
`603ccc47643d2e4bef66f464143bff13d929a133c7a575571dfbfc6761c37755`.
Generation took 121.11 s and 4,675,716 KiB peak RSS.

## Reversibility and truth size

The 800 events remove 17,943 GFF3 records and contain 1,488 distinct hidden
transcript structures. Evaluator-only restoration reproduced the normalized
source byte for byte with SHA-256
`4cd8fb0244ef51d7055d5866f11086c1fee927f5cfbb41509f118f6fb0f29dad`.
The 507 MB restored copy was then deleted because it is exactly reproducible;
the checksum-bearing restoration report remains.

## Scorer v2 sentinels

Scorer v2 used code `ce1361a5b20d52dca9118ed50b0444b6efac380e`;
archive SHA-256
`54a717132744d88ca09ec9a496388b35bf8c9c02ff07e72613d89155050b5fe2`.

| Metric | No-op | Oracle | Truth denominator |
|---|---:|---:|---:|
| Strict transcript recall | 0.0 | 1.0 | 1,488 structures |
| Exact complete-event recall | 0.0 | 1.0 | 800 events |
| Exact gene-group recall | 0.0 | 1.0 | 800 genes |
| Exon-segment recall | 0.0 | 1.0 | 4,816 segments |
| Splice-junction recall | 0.0 | 1.0 | 3,638 junctions |
| Phase-aware CDS-segment recall | 0.0 | 1.0 | 3,796 segments |
| CDS nucleotide-coverage recall | 0.0 | 1.0 | 926,236 bp |
| Collateral baseline losses | 0 | 0 | — |

Both sentinel quality gates passed: the source checksum matches truth and no
hidden transcript structure remains in the blind annotation. Scorer v1 reports
were retained under versioned filenames; the unversioned reports now use schema
`ploidypatch.annotation_repair_score.v2`.

## Next executable gate

The deliberately simple, non-PloidyPatch projection baseline is now frozen;
see `MINIPROT_BASELINE_v0.1.md`. The remaining gate is a negative control that
proposes no repair at assembly gaps and biologically genuine absence/PAV loci,
followed by a blind synteny-gap candidate generator. Only then should the
project measure whether WGD-aware evidence improves calibrated precision or
recall over the sequential baseline.
