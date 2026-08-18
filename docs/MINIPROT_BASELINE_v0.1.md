# Miniprot-gap baseline v0.1 — 2026-08-07

## Claim boundary

This is a deliberately simple external-protein projection baseline. It is not
PloidyPatch and was frozen before evaluator access. Its purpose is to establish
what a mature protein-to-genome aligner can recover without the proposed plant
event graph, calibrated uncertainty, or WGD-gap reasoning.

## Leakage controls and fixed inputs

The baseline used only:

- the 19-chromosome Da-Ae genome;
- the blind perturbed GFF3;
- one valid representative protein per `B. rapa` and `B. oleracea` gene.

It did not use Da-Ae proteins, the 800-row target selection, hidden truth, or
the evaluator-only per-gene WGDI strata. The two external protein sets were
prefixed and combined with a complete mapping table:

| Source | Protein records | Input SHA-256 |
|---|---:|---|
| `B. rapa` (`bra_a`) | 39,607 | `b12e39028c7a07b2c47bfb2a4d2410b14013b0f6ba087736257f8d5b00bc6568` |
| `B. oleracea` (`bol_c`) | 54,761 | `03707bb8548a9c2c581042aceaa9778ddfba7b27a80d96f9a496490d79b3dcc9` |

Combined FASTA: 94,368 records, SHA-256
`06c66ad19810947821a1c8f2405f6c928552f43589b06a5653f9c82e6a8aadf4`.

## Frozen aligner and adapter

Miniprot `0.18-r281` was linked from the existing server conda cache into the
isolated `ploidypatch-baseline` prefix; no package was downloaded. Explicit
environment SHA-256:
`fdb9e6dff9138826e99273bffb97373089e2e151253ad835612bdb7261b2d5bb`.

Genome indexing used 64 threads. The 1.7 GB index has SHA-256
`edea2c64d3e05227105ec8a155e3482f9c1a536bd189c7dab5c63f63cbd35646`;
indexing took 34.89 s and 8,747,228 KiB peak RSS.

The projection command was frozen as:

```bash
miniprot -I -t64 --gff-only -N4 --outn 4 --outs 0.8 --outc 0.5 \
  -P BMP target.mpi brassica_progenitors.pep.fa
```

It produced 259,911 models in 1,729,854 GFF3 lines. The raw GFF3 SHA-256 is
`b7cf7b9c598e495f9fc030bee7695a3e11e12becbc7e32b1e5d9fcadeed54438`.
Projection took 64.43 s and 6,419,540 KiB peak RSS.

The blind adapter applied these predeclared rules:

- amino-acid identity at least 0.5;
- query coverage at least 0.5;
- reject any frameshift or in-frame stop;
- projected CDS overlap with same-strand blind CDS at most 0.2;
- reject a lower-ranked projection when CDS overlap with an accepted projection
  exceeds 0.5 of the smaller projection.

Miniprot emits mRNA/CDS/stop-codon features, not full UTR exon models. The
adapter therefore defines coding exons as merged CDS intervals. This is an
explicit baseline limitation and makes strict full-transcript recovery harder
than CDS coverage recovery.

| Decision | Models |
|---|---:|
| accepted | 23,025 |
| overlaps blind CDS | 184,550 |
| non-intact alignment | 43,718 |
| redundant projection | 8,370 |
| identity below threshold | 248 |

The candidate GFF3 SHA-256 is
`dfba0a9ed53206c35ab1ff3f1d975f4bb46067821b3a2b4487d73e1ba5f3f753`.
After replacing a linear chromosome scan with an equivalent binary-search
interval index, adaptation took 33.27 s and 527,852 KiB peak RSS. The first
slow attempt was terminated before producing any artifact; thresholds and
decisions were not changed.

## Frozen formal result

No threshold was changed after evaluator access. Aggregate scorer v2 results:

| Metric | Precision | Recall | F1 |
|---|---:|---:|---:|
| exact transcript structure | 0.000261 | 0.004032 | 0.000490 |
| exact exon segments | 0.024967 | 0.308970 | 0.046200 |
| splice junctions | 0.051129 | 0.514294 | 0.093011 |
| phase-aware CDS segments | 0.036524 | 0.573498 | 0.068675 |
| CDS nucleotide coverage | 0.045467 | 0.638437 | 0.084889 |

Six of 1,488 hidden transcript structures were exact. Five of 800 events were
fully recovered with exact gene grouping. There were 23,019 strict novel-model
false positives and no collateral loss of baseline transcript structures.

The final anonymous stratified report SHA-256 is
`a7265160b530e9ec6f7e29ab4ba2fb130d50bcc38b828417d0cdcc4531253410`.

## Plant/WGD difficulty signal

Each WGDI stratum contains exactly 200 hidden events by design:

| Independent WGDI stratum | Strict transcript recall | Complete-event recall | CDS bp recall | Junction recall |
|---|---:|---:|---:|---:|
| expected and cross-progenitor blocks | 0.0000 | 0.0000 | 0.9713 | 0.7650 |
| expected-progenitor block only | 0.0076 | 0.0150 | 0.6859 | 0.4905 |
| cross-progenitor block only | 0.0092 | 0.0100 | 0.7056 | 0.5602 |
| no retained block | 0.0000 | 0.0000 | 0.1716 | 0.1284 |

The WGDI labels were used only after the global baseline was frozen and scored.
They therefore diagnose difficulty rather than steer predictions. The sharp
drop in `no_block` recall supports syntenic context as useful plant-specific
evidence, while the uniformly poor strict precision shows that projection
alone is not a viable repair policy.

Transcript-complexity marginals also show the expected limitation: complete
event recall is 1.04% for single-transcript genes and zero for two- or
three-plus-transcript genes. CDS bp recall remains 59.52%, 71.23%, and 66.89%,
respectively, because protein projection can locate coding sequence without
reconstructing the complete isoform set.

## Reproducibility and next gate

- Protein preparation code: `b9591e2b5dd7925906dadf08a20385aae8f47ae3`
- Adapter code: `1c0062fc557a82467d30bdfd3a039d5a50d25dd6`
- Stratified scorer code: `346d292689c00a71cbd9a81ad90388946cca7be1`
- Server results:
  `/data/codexli/projects/PloidyPatch/results/baselines/miniprot_brassica_v0.1/`

The next gate is a blind, synteny-gap candidate generator that infers missing
positions from retained neighboring anchors rather than receiving the hidden
gene-level WGDI table. It must rank or abstain among the 23,025 projection-like
possibilities and will be compared with this frozen v0.1 baseline without
retuning on these 800 events.
