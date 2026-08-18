# Independent soybean heldout benchmark v0.1

Date: 2026-08-07  
Status: frozen cross-lineage result

## Question and leakage boundary

This experiment asks whether the already frozen external-protein projection
layer transfers from the Brassica development system to an independent plant
lineage. It is a test of the candidate-generation layer, not evidence that the
full PloidyPatch event graph already outperforms existing tools.

- Target: `Glycine max`, assembly `Glycine_max_v2.1`, Ensembl Plants 62.
- Projection reference: `Glycine soja`, assembly `ASM419377v2`, Ensembl Plants
  62.
- Only the 20 primary `G. max` chromosomes were indexed.
- No target-derived proteins were placed in the reference set.
- The aligner and adapter thresholds were copied unchanged from the frozen
  Brassica experiment. No soybean result was used to tune them.
- Candidate GFF3 files and their complete-annotation control were checksummed
  before the evaluator read hidden truth.

The `G. max` target bundle passed the input gate: 55,897 genes, 88,412
transcripts, and no translation discrepancies. The `G. soja` annotation was
accepted as a projection-reference warning case: all 69,277 protein records
were structurally readable, while 38 transcript CDS/protein lengths disagreed
and 55 CDS records contained internal stops. This warning does not alter the
provided protein sequences, but it is retained in provenance and prevents the
reference annotation from being treated as truth.

## Frozen controlled perturbation

The primary-chromosome catalog contains 55,589 structurally unique eligible
genes. A deterministic 800-event sample was frozen with seed `20260809` over
12 `transcript_count_bin` by `max_exons_bin` cells:

- 400 one-transcript genes, 100 from each exon-complexity bin;
- 200 two-transcript genes, 50 from each exon-complexity bin;
- 200 genes with three or more transcripts, 50 from each exon-complexity bin.

The perturbation removes complete gene hierarchies from the annotation while
leaving the genome unchanged. Evaluator-only restoration recovered 26,913
records and reproduced the decompressed source text hash exactly.

Key custody hashes:

| Artifact | SHA-256 |
|---|---|
| candidate catalog | `9cd7819557eebcc26d181849e4d460ba73c0532208870afe2d78d5a2d582a794` |
| 800-row selection | `cd6d406299fe4df545cf94e540b152cef81e4bb1f3b7cd479131de2cc9f6e52a` |
| blind GFF3 | `fc8c83440805a2a7d0babf46915648ee5ed4dc1606f943bcf5ccccbadc5f28f8` |
| hidden truth | `56d4fa0aa7d13bcb018419ffec42f1c1cbc431588221f8a1b217cb5ea56bdbaa` |
| pre-score freeze list | `4dc39acf824122c7b1573d27cd5bb5b09b1b152160bd59beab40c947fb4e4730` |
| post-score freeze list | `6fc11a8ca8735b2726e03b46f4c39a6549904d1e6d7899be4c7a8b621dbc25b1` |

Server benchmark root:
`/data/codexli/projects/PloidyPatch/benchmark/heldout/v0.1/gma_v21_annotation_missing_gene_seed20260809`.

## Frozen projection baseline

Miniprot `0.18-r281` used the same Brassica settings:

```text
-I -t64 --gff-only -N4 --outn 4 --outs 0.8 --outc 0.5
```

The adapter also retained the frozen thresholds: minimum identity 0.5,
minimum query coverage 0.5, maximum existing-CDS overlap 0.2, maximum
projection redundancy overlap 0.5, and intact alignments only.

The index covers 949,183,385 bp in 20 chromosomes. Index construction took
9.5 s inside miniprot with 9.01 GB peak RSS. Projection of 69,277 proteins took
31.3 s inside miniprot with 5.30 GB peak RSS and emitted 147,274 models.

The blind adapter accepted 6,314 models and the complete-annotation control
accepted 5,663. Paired evaluation subtracts the 5,663 reproducible background
disagreements, leaving 651 deletion-responsive candidates. This subtraction is
performed only by the evaluator and is not available during candidate
generation.

## Heldout results

| Endpoint | Precision | Recall | F1 |
|---|---:|---:|---:|
| exact phased CDS chain | 79.42% | 43.05% | 55.83% |
| exact full transcript structure | 14.59% | 6.28% | 8.78% |
| CDS nucleotide coverage | 99.47% | 86.04% | 92.27% |
| phased CDS segments | 96.52% | 76.83% | 85.56% |
| splice junctions | 98.45% | 81.24% | 89.02% |
| exon segments | 79.83% | 58.07% | 67.23% |

There were 517 exact phased CDS true positives among 651 differential
candidates, with 134 false positives. At the event level, at least one complete
CDS chain was recovered for 323 of 800 hidden genes (40.38%). Only 8 of 800
complete multi-isoform gene structures were restored exactly.

Strict CDS-chain recall by transcript-count bin was:

| Hidden gene complexity | Recovered / truth CDS chains | Recall |
|---|---:|---:|
| one transcript | 208 / 400 | 52.00% |
| two transcripts | 151 / 324 | 46.60% |
| three or more transcripts | 158 / 477 | 33.12% |

The result is consistent with the known boundary of protein-to-genome
projection: it recovers coding sequence and splice junctions well, but does not
reconstruct complete plant isoform sets or untranslated exon boundaries. That
is why exact phased CDS chain remains the primary endpoint for this component,
while full transcript structure is reported as a separate, harder task.

## Interpretation

The cross-lineage result is a positive engineering gate: the frozen projection
and paired-evaluation layer generalizes to a paleopolyploid legume and improves
from the Brassica formal result (57.8% precision and 30.1% recall for exact CDS
chains) without soybean tuning. It also makes the remaining gap more precise:

1. projection is a credible high-recall candidate source, not a final repair
   policy;
2. coding-region confidence and complete-isoform confidence must be distinct;
3. multi-transcript genes require RNA/full-length transcript evidence and
   explicit alternative-structure inference;
4. natural source disagreements must remain separate from perturbation-induced
   responses;
5. the plant-specific contribution still has to beat SynGAP/GeMoMa/LiftOn and
   sequential synteny-plus-orthology baselines on matched inputs.

The next scientific gate is therefore not more tuning on these 800 genes. It
is a matched baseline comparison plus natural-locus validation using
independent transcript, assembly, and copy-specific WGD evidence.
