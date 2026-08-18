# PloidyPatch v1.0 command and parameter reference

This file is generated from `ploidypatch.cli.build_parser` for package
version `1.0.0`. It documents 87 executable leaf commands.
Do not edit it by hand; run `scripts/export_cli_parameters_v1.py`.

Every output path is non-overwriting unless a command explicitly documents
otherwise. Built-in `--help` remains the authoritative runtime view.

## Contents

- [`ploidypatch audit`](#ploidypatch-audit)
- [`ploidypatch baseline adapt-gff`](#ploidypatch-baseline-adapt-gff)
- [`ploidypatch baseline adapt-miniprot`](#ploidypatch-baseline-adapt-miniprot)
- [`ploidypatch baseline adapt-miniprot-structure`](#ploidypatch-baseline-adapt-miniprot-structure)
- [`ploidypatch baseline infer-structure-hypotheses`](#ploidypatch-baseline-infer-structure-hypotheses)
- [`ploidypatch baseline merge-candidate-gffs`](#ploidypatch-baseline-merge-candidate-gffs)
- [`ploidypatch baseline prepare-proteins`](#ploidypatch-baseline-prepare-proteins)
- [`ploidypatch baseline select-method-consensus`](#ploidypatch-baseline-select-method-consensus)
- [`ploidypatch benchmark audit-mask`](#ploidypatch-benchmark-audit-mask)
- [`ploidypatch benchmark bootstrap-confusion`](#ploidypatch-benchmark-bootstrap-confusion)
- [`ploidypatch benchmark bootstrap-events`](#ploidypatch-benchmark-bootstrap-events)
- [`ploidypatch benchmark bootstrap-independent-events`](#ploidypatch-benchmark-bootstrap-independent-events)
- [`ploidypatch benchmark bootstrap-structure-hypotheses`](#ploidypatch-benchmark-bootstrap-structure-hypotheses)
- [`ploidypatch benchmark calibrate-synteny-tiers`](#ploidypatch-benchmark-calibrate-synteny-tiers)
- [`ploidypatch benchmark catalog`](#ploidypatch-benchmark-catalog)
- [`ploidypatch benchmark label-copy-features`](#ploidypatch-benchmark-label-copy-features)
- [`ploidypatch benchmark label-synteny-models`](#ploidypatch-benchmark-label-synteny-models)
- [`ploidypatch benchmark mask-genome`](#ploidypatch-benchmark-mask-genome)
- [`ploidypatch benchmark perturb`](#ploidypatch-benchmark-perturb)
- [`ploidypatch benchmark restore`](#ploidypatch-benchmark-restore)
- [`ploidypatch benchmark sample-catalog`](#ploidypatch-benchmark-sample-catalog)
- [`ploidypatch benchmark sample-copy-pairs`](#ploidypatch-benchmark-sample-copy-pairs)
- [`ploidypatch benchmark score`](#ploidypatch-benchmark-score)
- [`ploidypatch benchmark score-abstention`](#ploidypatch-benchmark-score-abstention)
- [`ploidypatch benchmark score-copy-ranking`](#ploidypatch-benchmark-score-copy-ranking)
- [`ploidypatch benchmark score-localization`](#ploidypatch-benchmark-score-localization)
- [`ploidypatch benchmark score-structure-hypotheses`](#ploidypatch-benchmark-score-structure-hypotheses)
- [`ploidypatch benchmark summarize-catalog`](#ploidypatch-benchmark-summarize-catalog)
- [`ploidypatch benchmark summarize-mask-selection`](#ploidypatch-benchmark-summarize-mask-selection)
- [`ploidypatch evidence aggregate-junctions`](#ploidypatch-evidence-aggregate-junctions)
- [`ploidypatch evidence aggregate-reference-anchored`](#ploidypatch-evidence-aggregate-reference-anchored)
- [`ploidypatch evidence annotate-natural-assembly`](#ploidypatch-evidence-annotate-natural-assembly)
- [`ploidypatch evidence annotate-pav-wgdi`](#ploidypatch-evidence-annotate-pav-wgdi)
- [`ploidypatch evidence apply-conflict-winner-guard`](#ploidypatch-evidence-apply-conflict-winner-guard)
- [`ploidypatch evidence audit-natural-candidates`](#ploidypatch-evidence-audit-natural-candidates)
- [`ploidypatch evidence bootstrap-isoseq-review-yield`](#ploidypatch-evidence-bootstrap-isoseq-review-yield)
- [`ploidypatch evidence build-copy-features`](#ploidypatch-evidence-build-copy-features)
- [`ploidypatch evidence build-fixed-target-backbone`](#ploidypatch-evidence-build-fixed-target-backbone)
- [`ploidypatch evidence build-homeolog-topology-features`](#ploidypatch-evidence-build-homeolog-topology-features)
- [`ploidypatch evidence catalog-pav-genes`](#ploidypatch-evidence-catalog-pav-genes)
- [`ploidypatch evidence collapse-pav-events`](#ploidypatch-evidence-collapse-pav-events)
- [`ploidypatch evidence combine-chromosome-pafs`](#ploidypatch-evidence-combine-chromosome-pafs)
- [`ploidypatch evidence discover-natural`](#ploidypatch-evidence-discover-natural)
- [`ploidypatch evidence evaluate-species-applicability`](#ploidypatch-evidence-evaluate-species-applicability)
- [`ploidypatch evidence export-natural-candidate-cds`](#ploidypatch-evidence-export-natural-candidate-cds)
- [`ploidypatch evidence extract-bam-junctions`](#ploidypatch-evidence-extract-bam-junctions)
- [`ploidypatch evidence extract-paf-deletions`](#ploidypatch-evidence-extract-paf-deletions)
- [`ploidypatch evidence filter-candidate-query-paf`](#ploidypatch-evidence-filter-candidate-query-paf)
- [`ploidypatch evidence freeze-homeolog-review-rankings`](#ploidypatch-evidence-freeze-homeolog-review-rankings)
- [`ploidypatch evidence infer-homeolog-pairs`](#ploidypatch-evidence-infer-homeolog-pairs)
- [`ploidypatch evidence infer-outgroup-duplicated-pairs`](#ploidypatch-evidence-infer-outgroup-duplicated-pairs)
- [`ploidypatch evidence infer-self-wgd-pairs`](#ploidypatch-evidence-infer-self-wgd-pairs)
- [`ploidypatch evidence infer-synteny-gaps`](#ploidypatch-evidence-infer-synteny-gaps)
- [`ploidypatch evidence intersect-copy-pair-evidence`](#ploidypatch-evidence-intersect-copy-pair-evidence)
- [`ploidypatch evidence join-isoseq-review-rankings`](#ploidypatch-evidence-join-isoseq-review-rankings)
- [`ploidypatch evidence prepare-b73-isoseq`](#ploidypatch-evidence-prepare-b73-isoseq)
- [`ploidypatch evidence prepare-natural-graph`](#ploidypatch-evidence-prepare-natural-graph)
- [`ploidypatch evidence prepare-wgdi`](#ploidypatch-evidence-prepare-wgdi)
- [`ploidypatch evidence project-fixed-backbone`](#ploidypatch-evidence-project-fixed-backbone)
- [`ploidypatch evidence propagate-wgd-conflict-partners`](#ploidypatch-evidence-propagate-wgd-conflict-partners)
- [`ploidypatch evidence reconcile-pav-catalogs`](#ploidypatch-evidence-reconcile-pav-catalogs)
- [`ploidypatch evidence score-copy-candidates`](#ploidypatch-evidence-score-copy-candidates)
- [`ploidypatch evidence score-homeolog-copy-candidates`](#ploidypatch-evidence-score-homeolog-copy-candidates)
- [`ploidypatch evidence score-support-conditioned-candidates`](#ploidypatch-evidence-score-support-conditioned-candidates)
- [`ploidypatch evidence screen-pav-proteins`](#ploidypatch-evidence-screen-pav-proteins)
- [`ploidypatch evidence select-projection-support`](#ploidypatch-evidence-select-projection-support)
- [`ploidypatch evidence select-scored-copy-candidates`](#ploidypatch-evidence-select-scored-copy-candidates)
- [`ploidypatch evidence select-synteny-gap-models`](#ploidypatch-evidence-select-synteny-gap-models)
- [`ploidypatch evidence select-wgd-supported-candidates`](#ploidypatch-evidence-select-wgd-supported-candidates)
- [`ploidypatch evidence subset-pav-proteins`](#ploidypatch-evidence-subset-pav-proteins)
- [`ploidypatch evidence summarize-natural`](#ploidypatch-evidence-summarize-natural)
- [`ploidypatch evidence summarize-projection-support`](#ploidypatch-evidence-summarize-projection-support)
- [`ploidypatch evidence summarize-wgdi`](#ploidypatch-evidence-summarize-wgdi)
- [`ploidypatch evidence validate-isoseq-candidates`](#ploidypatch-evidence-validate-isoseq-candidates)
- [`ploidypatch evidence validate-natural-rna`](#ploidypatch-evidence-validate-natural-rna)
- [`ploidypatch evidence validate-natural-secondary-rna`](#ploidypatch-evidence-validate-natural-secondary-rna)
- [`ploidypatch graph infer`](#ploidypatch-graph-infer)
- [`ploidypatch normalize ncbi-primary`](#ploidypatch-normalize-ncbi-primary)
- [`ploidypatch normalize primary-annotation`](#ploidypatch-normalize-primary-annotation)
- [`ploidypatch normalize provider-gff3`](#ploidypatch-normalize-provider-gff3)
- [`ploidypatch patch apply`](#ploidypatch-patch-apply)
- [`ploidypatch patch compile-copy-additions`](#ploidypatch-patch-compile-copy-additions)
- [`ploidypatch patch compile-reviewed-copy-additions`](#ploidypatch-patch-compile-reviewed-copy-additions)
- [`ploidypatch patch compile-structure`](#ploidypatch-patch-compile-structure)
- [`ploidypatch patch create`](#ploidypatch-patch-create)
- [`ploidypatch patch revert`](#ploidypatch-patch-revert)
- [`ploidypatch report`](#ploidypatch-report)

## `ploidypatch audit`

ploidypatch audit

```text
usage: ploidypatch audit [-h] --gff GFF --protein PROTEIN --cds CDS
                         [--fai FAI] [--output OUTPUT] [--checksums]
```

| Argument | Required | Type | Default | Choices | Description |
|---|---:|---|---|---|---|
| `-h`<br>`--help` | no | `string` | `—` | — | show this help message and exit |
| `--gff` | yes | `string` | `None` | — | GFF3 or GFF3.gz path |
| `--protein` | yes | `string` | `None` | — | protein FASTA path |
| `--cds` | yes | `string` | `None` | — | CDS FASTA path |
| `--fai` | no | `string` | `None` | — | optional FASTA index for coordinate checks |
| `--output` | no | `string` | `None` | — | write JSON report to this path |
| `--checksums` | no | `flag` | `false` | — | calculate SHA-256 for audited inputs |

## `ploidypatch baseline adapt-gff`

ploidypatch baseline adapt-gff

```text
usage: ploidypatch baseline adapt-gff [-h] --perturbed-gff PERTURBED_GFF
                                      --candidate-gff CANDIDATE_GFF --source
                                      SOURCE --output-gff OUTPUT_GFF
                                      --decisions-tsv DECISIONS_TSV
                                      [--score-attribute SCORE_ATTRIBUTE]
                                      [--max-existing-cds-overlap MAX_EXISTING_CDS_OVERLAP]
                                      [--max-redundancy-overlap MAX_REDUNDANCY_OVERLAP]
                                      [--infer-missing-cds-phase]
```

| Argument | Required | Type | Default | Choices | Description |
|---|---:|---|---|---|---|
| `-h`<br>`--help` | no | `string` | `—` | — | show this help message and exit |
| `--perturbed-gff` | yes | `string` | `None` | — | — |
| `--candidate-gff` | yes | `string` | `None` | — | — |
| `--source` | yes | `string` | `None` | — | — |
| `--output-gff` | yes | `string` | `None` | — | — |
| `--decisions-tsv` | yes | `string` | `None` | — | — |
| `--score-attribute` | no | `string` | `None` | — | — |
| `--max-existing-cds-overlap` | no | `float` | `0.2` | — | — |
| `--max-redundancy-overlap` | no | `float` | `0.5` | — | — |
| `--infer-missing-cds-phase` | no | `flag` | `false` | — | infer phases only when every CDS phase in a model is missing; assumes a complete CDS whose first phase is zero |

## `ploidypatch baseline adapt-miniprot`

ploidypatch baseline adapt-miniprot

```text
usage: ploidypatch baseline adapt-miniprot [-h] --perturbed-gff PERTURBED_GFF
                                           --miniprot-gff MINIPROT_GFF
                                           --protein-map PROTEIN_MAP
                                           --output-gff OUTPUT_GFF
                                           --decisions-tsv DECISIONS_TSV
                                           [--min-identity MIN_IDENTITY]
                                           [--min-query-coverage MIN_QUERY_COVERAGE]
                                           [--max-existing-cds-overlap MAX_EXISTING_CDS_OVERLAP]
                                           [--max-redundancy-overlap MAX_REDUNDANCY_OVERLAP]
                                           [--allow-disrupted]
```

| Argument | Required | Type | Default | Choices | Description |
|---|---:|---|---|---|---|
| `-h`<br>`--help` | no | `string` | `—` | — | show this help message and exit |
| `--perturbed-gff` | yes | `string` | `None` | — | — |
| `--miniprot-gff` | yes | `string` | `None` | — | — |
| `--protein-map` | yes | `string` | `None` | — | — |
| `--output-gff` | yes | `string` | `None` | — | — |
| `--decisions-tsv` | yes | `string` | `None` | — | — |
| `--min-identity` | no | `float` | `0.5` | — | — |
| `--min-query-coverage` | no | `float` | `0.5` | — | — |
| `--max-existing-cds-overlap` | no | `float` | `0.2` | — | — |
| `--max-redundancy-overlap` | no | `float` | `0.5` | — | — |
| `--allow-disrupted` | no | `flag` | `false` | — | retain alignments with frameshifts or in-frame stops |

## `ploidypatch baseline adapt-miniprot-structure`

ploidypatch baseline adapt-miniprot-structure

```text
usage: ploidypatch baseline adapt-miniprot-structure [-h] --annotation-gff
                                                     ANNOTATION_GFF
                                                     --miniprot-gff
                                                     MINIPROT_GFF
                                                     --protein-map PROTEIN_MAP
                                                     [--source-group-map SOURCE_GROUP_MAP]
                                                     --output-gff OUTPUT_GFF
                                                     --decisions-tsv
                                                     DECISIONS_TSV
                                                     [--min-identity MIN_IDENTITY]
                                                     [--min-query-coverage MIN_QUERY_COVERAGE]
                                                     [--min-source-support MIN_SOURCE_SUPPORT]
                                                     [--min-gene-overlap-fraction MIN_GENE_OVERLAP_FRACTION]
                                                     [--allow-disrupted]
```

| Argument | Required | Type | Default | Choices | Description |
|---|---:|---|---|---|---|
| `-h`<br>`--help` | no | `string` | `—` | — | show this help message and exit |
| `--annotation-gff` | yes | `string` | `None` | — | — |
| `--miniprot-gff` | yes | `string` | `None` | — | — |
| `--protein-map` | yes | `string` | `None` | — | — |
| `--source-group-map` | no | `string` | `None` | — | optional TSV mapping each protein source to an independent support group |
| `--output-gff` | yes | `string` | `None` | — | — |
| `--decisions-tsv` | yes | `string` | `None` | — | — |
| `--min-identity` | no | `float` | `0.5` | — | — |
| `--min-query-coverage` | no | `float` | `0.5` | — | — |
| `--min-source-support` | no | `int` | `1` | — | — |
| `--min-gene-overlap-fraction` | no | `float` | `0.1` | — | — |
| `--allow-disrupted` | no | `flag` | `false` | — | retain alignments with frameshifts or in-frame stops |

## `ploidypatch baseline infer-structure-hypotheses`

ploidypatch baseline infer-structure-hypotheses

```text
usage: ploidypatch baseline infer-structure-hypotheses [-h] --annotation-gff
                                                       ANNOTATION_GFF
                                                       --candidate-gff
                                                       CANDIDATE_GFF
                                                       --output-tsv OUTPUT_TSV
                                                       --candidate-topology-tsv
                                                       CANDIDATE_TOPOLOGY_TSV
```

| Argument | Required | Type | Default | Choices | Description |
|---|---:|---|---|---|---|
| `-h`<br>`--help` | no | `string` | `—` | — | show this help message and exit |
| `--annotation-gff` | yes | `string` | `None` | — | — |
| `--candidate-gff` | yes | `string` | `None` | — | — |
| `--output-tsv` | yes | `string` | `None` | — | — |
| `--candidate-topology-tsv` | yes | `string` | `None` | — | — |

## `ploidypatch baseline merge-candidate-gffs`

ploidypatch baseline merge-candidate-gffs

```text
usage: ploidypatch baseline merge-candidate-gffs [-h] --candidate CANDIDATE
                                                 --output-gff OUTPUT_GFF
                                                 --provenance-tsv
                                                 PROVENANCE_TSV
```

| Argument | Required | Type | Default | Choices | Description |
|---|---:|---|---|---|---|
| `-h`<br>`--help` | no | `string` | `—` | — | show this help message and exit |
| `--candidate` | yes | `string` | `None` | — | repeatable REFERENCE_SOURCE=GFF3 input |
| `--output-gff` | yes | `string` | `None` | — | — |
| `--provenance-tsv` | yes | `string` | `None` | — | — |

## `ploidypatch baseline prepare-proteins`

ploidypatch baseline prepare-proteins

```text
usage: ploidypatch baseline prepare-proteins [-h] --protein PROTEIN
                                             --output-fasta OUTPUT_FASTA
                                             --output-map OUTPUT_MAP
```

| Argument | Required | Type | Default | Choices | Description |
|---|---:|---|---|---|---|
| `-h`<br>`--help` | no | `string` | `—` | — | show this help message and exit |
| `--protein` | yes | `string` | `None` | — | repeatable SOURCE=PATH input |
| `--output-fasta` | yes | `string` | `None` | — | — |
| `--output-map` | yes | `string` | `None` | — | — |

## `ploidypatch baseline select-method-consensus`

ploidypatch baseline select-method-consensus

```text
usage: ploidypatch baseline select-method-consensus [-h] --base-gff BASE_GFF
                                                    --candidate CANDIDATE
                                                    --output-gff OUTPUT_GFF
                                                    --decisions-tsv
                                                    DECISIONS_TSV
                                                    [--min-method-support MIN_METHOD_SUPPORT]
                                                    [--max-redundancy-overlap MAX_REDUNDANCY_OVERLAP]
                                                    [--redundancy-policy {suppress_overlapping,retain_distinct_chains}]
```

| Argument | Required | Type | Default | Choices | Description |
|---|---:|---|---|---|---|
| `-h`<br>`--help` | no | `string` | `—` | — | show this help message and exit |
| `--base-gff` | yes | `string` | `None` | — | — |
| `--candidate` | yes | `string` | `None` | — | repeatable METHOD=adapted_candidate.gff3 input |
| `--output-gff` | yes | `string` | `None` | — | — |
| `--decisions-tsv` | yes | `string` | `None` | — | — |
| `--min-method-support` | no | `int` | `2` | — | — |
| `--max-redundancy-overlap` | no | `float` | `0.5` | — | — |
| `--redundancy-policy` | no | `string` | `suppress_overlapping` | `suppress_overlapping`, `retain_distinct_chains` | suppress strongly overlapping chains (v1 default), or retain every distinct phased CDS chain and label conflict sets |

## `ploidypatch benchmark audit-mask`

ploidypatch benchmark audit-mask

```text
usage: ploidypatch benchmark audit-mask [-h] --source-genome SOURCE_GENOME
                                        --masked-genome MASKED_GENOME
                                        --mask-truth MASK_TRUTH
```

| Argument | Required | Type | Default | Choices | Description |
|---|---:|---|---|---|---|
| `-h`<br>`--help` | no | `string` | `—` | — | show this help message and exit |
| `--source-genome` | yes | `string` | `None` | — | — |
| `--masked-genome` | yes | `string` | `None` | — | — |
| `--mask-truth` | yes | `string` | `None` | — | — |

## `ploidypatch benchmark bootstrap-confusion`

ploidypatch benchmark bootstrap-confusion

```text
usage: ploidypatch benchmark bootstrap-confusion [-h] --score SCORE
                                                 --output-json OUTPUT_JSON
                                                 [--section {strict_cds_chain,strict_transcript_structure}]
                                                 [--replicates REPLICATES]
                                                 [--seed SEED] [--alpha ALPHA]
```

| Argument | Required | Type | Default | Choices | Description |
|---|---:|---|---|---|---|
| `-h`<br>`--help` | no | `string` | `—` | — | show this help message and exit |
| `--score` | yes | `string` | `None` | — | repeatable LABEL=score.json input |
| `--output-json` | yes | `string` | `None` | — | — |
| `--section` | no | `string` | `strict_cds_chain` | `strict_cds_chain`, `strict_transcript_structure` | — |
| `--replicates` | no | `int` | `10000` | — | — |
| `--seed` | no | `int` | `20260807` | — | — |
| `--alpha` | no | `float` | `0.05` | — | — |

## `ploidypatch benchmark bootstrap-events`

ploidypatch benchmark bootstrap-events

```text
usage: ploidypatch benchmark bootstrap-events [-h] --score SCORE --output-json
                                              OUTPUT_JSON
                                              [--metric {complete_cds_chain_recovery,complete_error_removal,complete_transcript_recovery,exact_cds_gene_grouping,exact_gene_grouping}]
                                              [--replicates REPLICATES]
                                              [--seed SEED] [--alpha ALPHA]
```

| Argument | Required | Type | Default | Choices | Description |
|---|---:|---|---|---|---|
| `-h`<br>`--help` | no | `string` | `—` | — | show this help message and exit |
| `--score` | yes | `string` | `None` | — | repeatable LABEL=score.json input |
| `--output-json` | yes | `string` | `None` | — | — |
| `--metric` | no | `string` | `complete_cds_chain_recovery` | `complete_cds_chain_recovery`, `complete_error_removal`, `complete_transcript_recovery`, `exact_cds_gene_grouping`, `exact_gene_grouping` | — |
| `--replicates` | no | `int` | `10000` | — | — |
| `--seed` | no | `int` | `20260807` | — | — |
| `--alpha` | no | `float` | `0.05` | — | — |

## `ploidypatch benchmark bootstrap-independent-events`

ploidypatch benchmark bootstrap-independent-events

```text
usage: ploidypatch benchmark bootstrap-independent-events [-h] --score SCORE
                                                          --output-json
                                                          OUTPUT_JSON
                                                          [--metric {complete_cds_chain_recovery,complete_error_removal,complete_transcript_recovery,exact_cds_gene_grouping,exact_gene_grouping}]
                                                          [--replicates REPLICATES]
                                                          [--seed SEED]
                                                          [--alpha ALPHA]
```

| Argument | Required | Type | Default | Choices | Description |
|---|---:|---|---|---|---|
| `-h`<br>`--help` | no | `string` | `—` | — | show this help message and exit |
| `--score` | yes | `string` | `None` | — | repeatable LABEL=score.json input |
| `--output-json` | yes | `string` | `None` | — | — |
| `--metric` | no | `string` | `complete_cds_chain_recovery` | `complete_cds_chain_recovery`, `complete_error_removal`, `complete_transcript_recovery`, `exact_cds_gene_grouping`, `exact_gene_grouping` | — |
| `--replicates` | no | `int` | `10000` | — | — |
| `--seed` | no | `int` | `20260807` | — | — |
| `--alpha` | no | `float` | `0.05` | — | — |

## `ploidypatch benchmark bootstrap-structure-hypotheses`

ploidypatch benchmark bootstrap-structure-hypotheses

```text
usage: ploidypatch benchmark bootstrap-structure-hypotheses
       [-h] --score SCORE --output-json OUTPUT_JSON [--replicates REPLICATES]
       [--seed SEED] [--alpha ALPHA]
```

| Argument | Required | Type | Default | Choices | Description |
|---|---:|---|---|---|---|
| `-h`<br>`--help` | no | `string` | `—` | — | show this help message and exit |
| `--score` | yes | `string` | `None` | — | repeatable LABEL=score.json; a label may span event files |
| `--output-json` | yes | `string` | `None` | — | — |
| `--replicates` | no | `int` | `10000` | — | — |
| `--seed` | no | `int` | `20260807` | — | — |
| `--alpha` | no | `float` | `0.05` | — | — |

## `ploidypatch benchmark calibrate-synteny-tiers`

ploidypatch benchmark calibrate-synteny-tiers

```text
usage: ploidypatch benchmark calibrate-synteny-tiers [-h] --labels LABELS
                                                     --output OUTPUT
                                                     [--label-column LABEL_COLUMN]
                                                     [--eligible-column ELIGIBLE_COLUMN]
                                                     [--feature-set {synteny,projection_support}]
                                                     [--high-precision-floor HIGH_PRECISION_FLOOR]
                                                     [--high-precision-min-selected HIGH_PRECISION_MIN_SELECTED]
```

| Argument | Required | Type | Default | Choices | Description |
|---|---:|---|---|---|---|
| `-h`<br>`--help` | no | `string` | `—` | — | show this help message and exit |
| `--labels` | yes | `string` | `None` | — | — |
| `--output` | yes | `string` | `None` | — | — |
| `--label-column` | no | `string` | `label_exact_cds_chain` | — | — |
| `--eligible-column` | no | `string` | `None` | — | optional binary column restricting rows eligible for calibration |
| `--feature-set` | no | `string` | `synteny` | `synteny`, `projection_support` | — |
| `--high-precision-floor` | no | `float` | `0.1` | — | — |
| `--high-precision-min-selected` | no | `int` | `20` | — | — |

## `ploidypatch benchmark catalog`

ploidypatch benchmark catalog

```text
usage: ploidypatch benchmark catalog [-h] --gff GFF --output-tsv OUTPUT_TSV
                                     [--external-strata EXTERNAL_STRATA]
                                     [--external-strata-prefix EXTERNAL_STRATA_PREFIX]
                                     [--primary-chromosomes-only]
```

| Argument | Required | Type | Default | Choices | Description |
|---|---:|---|---|---|---|
| `-h`<br>`--help` | no | `string` | `—` | — | show this help message and exit |
| `--gff` | yes | `string` | `None` | — | — |
| `--output-tsv` | yes | `string` | `None` | — | — |
| `--external-strata` | no | `string` | `None` | — | optional TSV keyed by exact gene_id, for WGD/subgenome labels |
| `--external-strata-prefix` | no | `string` | `` | — | prefix added to every imported strata column |
| `--primary-chromosomes-only` | no | `flag` | `false` | — | retain NCBI region features marked genome=chromosome |

## `ploidypatch benchmark label-copy-features`

ploidypatch benchmark label-copy-features

```text
usage: ploidypatch benchmark label-copy-features [-h] --features FEATURES
                                                 --truth TRUTH --output-tsv
                                                 OUTPUT_TSV
```

| Argument | Required | Type | Default | Choices | Description |
|---|---:|---|---|---|---|
| `-h`<br>`--help` | no | `string` | `—` | — | show this help message and exit |
| `--features` | yes | `string` | `None` | — | — |
| `--truth` | yes | `string` | `None` | — | — |
| `--output-tsv` | yes | `string` | `None` | — | — |

## `ploidypatch benchmark label-synteny-models`

ploidypatch benchmark label-synteny-models

```text
usage: ploidypatch benchmark label-synteny-models [-h] --source-gff SOURCE_GFF
                                                  --candidate-gff
                                                  CANDIDATE_GFF --selection
                                                  SELECTION
                                                  --baseline-decisions
                                                  BASELINE_DECISIONS --truth
                                                  TRUTH --output-tsv
                                                  OUTPUT_TSV
                                                  [--control-candidate-gff CONTROL_CANDIDATE_GFF]
```

| Argument | Required | Type | Default | Choices | Description |
|---|---:|---|---|---|---|
| `-h`<br>`--help` | no | `string` | `—` | — | show this help message and exit |
| `--source-gff` | yes | `string` | `None` | — | — |
| `--candidate-gff` | yes | `string` | `None` | — | — |
| `--selection` | yes | `string` | `None` | — | — |
| `--baseline-decisions` | yes | `string` | `None` | — | — |
| `--truth` | yes | `string` | `None` | — | — |
| `--output-tsv` | yes | `string` | `None` | — | — |
| `--control-candidate-gff` | no | `string` | `None` | — | — |

## `ploidypatch benchmark mask-genome`

ploidypatch benchmark mask-genome

```text
usage: ploidypatch benchmark mask-genome [-h] --genome GENOME --truth TRUTH
                                         --output-genome OUTPUT_GENOME
                                         --output-mask-truth OUTPUT_MASK_TRUTH
                                         [--background-gff BACKGROUND_GFF]
```

| Argument | Required | Type | Default | Choices | Description |
|---|---:|---|---|---|---|
| `-h`<br>`--help` | no | `string` | `—` | — | show this help message and exit |
| `--genome` | yes | `string` | `None` | — | — |
| `--truth` | yes | `string` | `None` | — | — |
| `--output-genome` | yes | `string` | `None` | — | — |
| `--output-mask-truth` | yes | `string` | `None` | — | — |
| `--background-gff` | no | `string` | `None` | — | exclude target spans overlapping retained exon/CDS features in this GFF3 |

## `ploidypatch benchmark perturb`

ploidypatch benchmark perturb

```text
usage: ploidypatch benchmark perturb [-h] --gff GFF --output-dir OUTPUT_DIR
                                     [--truth-dir TRUTH_DIR]
                                     [--event-type {annotation_missing_gene,annotation_missing_internal_exon,annotation_boundary_shift,annotation_split_gene,annotation_fused_gene,annotation_copy_collapse}]
                                     (--count COUNT | --selection-tsv SELECTION_TSV)
                                     [--seed SEED] [--pair-tsv PAIR_TSV]
```

| Argument | Required | Type | Default | Choices | Description |
|---|---:|---|---|---|---|
| `-h`<br>`--help` | no | `string` | `—` | — | show this help message and exit |
| `--gff` | yes | `string` | `None` | — | source GFF3 or GFF3.gz |
| `--output-dir` | yes | `string` | `None` | — | new or empty benchmark output directory |
| `--truth-dir` | no | `string` | `None` | — | separate evaluator-only directory for hidden_truth.json |
| `--event-type` | no | `string` | `annotation_missing_gene` | `annotation_missing_gene`, `annotation_missing_internal_exon`, `annotation_boundary_shift`, `annotation_split_gene`, `annotation_fused_gene`, `annotation_copy_collapse` | — |
| `--count` | no | `int` | `None` | — | — |
| `--selection-tsv` | no | `string` | `None` | — | evaluator-only catalog sample with an adjacent verified manifest |
| `--seed` | no | `int` | `1` | — | — |
| `--pair-tsv` | no | `string` | `None` | — | evaluator-only gene_id_a/gene_id_b table required for annotation_copy_collapse |

## `ploidypatch benchmark restore`

ploidypatch benchmark restore

```text
usage: ploidypatch benchmark restore [-h] --perturbed-gff PERTURBED_GFF
                                     --truth TRUTH --output-gff OUTPUT_GFF
```

| Argument | Required | Type | Default | Choices | Description |
|---|---:|---|---|---|---|
| `-h`<br>`--help` | no | `string` | `—` | — | show this help message and exit |
| `--perturbed-gff` | yes | `string` | `None` | — | — |
| `--truth` | yes | `string` | `None` | — | — |
| `--output-gff` | yes | `string` | `None` | — | — |

## `ploidypatch benchmark sample-catalog`

ploidypatch benchmark sample-catalog

```text
usage: ploidypatch benchmark sample-catalog [-h] --catalog CATALOG --plan PLAN
                                            --output-tsv OUTPUT_TSV --seed
                                            SEED
                                            [--exclude-selection EXCLUDE_SELECTION]
```

| Argument | Required | Type | Default | Choices | Description |
|---|---:|---|---|---|---|
| `-h`<br>`--help` | no | `string` | `—` | — | show this help message and exit |
| `--catalog` | yes | `string` | `None` | — | — |
| `--plan` | yes | `string` | `None` | — | — |
| `--output-tsv` | yes | `string` | `None` | — | — |
| `--seed` | yes | `int` | `None` | — | — |
| `--exclude-selection` | no | `string` | `[]` | — | repeatable prior selection TSV whose candidate IDs are excluded |

## `ploidypatch benchmark sample-copy-pairs`

ploidypatch benchmark sample-copy-pairs

```text
usage: ploidypatch benchmark sample-copy-pairs [-h] --source-gff SOURCE_GFF
                                               --pairs PAIRS --count COUNT
                                               --seed SEED --output-pairs
                                               OUTPUT_PAIRS --decisions-tsv
                                               DECISIONS_TSV
```

| Argument | Required | Type | Default | Choices | Description |
|---|---:|---|---|---|---|
| `-h`<br>`--help` | no | `string` | `—` | — | show this help message and exit |
| `--source-gff` | yes | `string` | `None` | — | — |
| `--pairs` | yes | `string` | `None` | — | — |
| `--count` | yes | `int` | `None` | — | — |
| `--seed` | yes | `int` | `None` | — | — |
| `--output-pairs` | yes | `string` | `None` | — | — |
| `--decisions-tsv` | yes | `string` | `None` | — | — |

## `ploidypatch benchmark score`

ploidypatch benchmark score

```text
usage: ploidypatch benchmark score [-h] --source-gff SOURCE_GFF
                                   --perturbed-gff PERTURBED_GFF
                                   --candidate-gff CANDIDATE_GFF --truth TRUTH
                                   [--include-event-details]
                                   [--control-candidate-gff CONTROL_CANDIDATE_GFF]
                                   [--event-strata EVENT_STRATA]
                                   [--stratum-column STRATUM_COLUMN]
```

| Argument | Required | Type | Default | Choices | Description |
|---|---:|---|---|---|---|
| `-h`<br>`--help` | no | `string` | `—` | — | show this help message and exit |
| `--source-gff` | yes | `string` | `None` | — | — |
| `--perturbed-gff` | yes | `string` | `None` | — | — |
| `--candidate-gff` | yes | `string` | `None` | — | — |
| `--truth` | yes | `string` | `None` | — | — |
| `--include-event-details` | no | `flag` | `false` | — | include target IDs; keep disabled for aggregate-only tuning |
| `--control-candidate-gff` | no | `string` | `None` | — | optional same-method candidate generated from the complete annotation; subtracts background proposals in controlled perturbation evaluation |
| `--event-strata` | no | `string` | `None` | — | evaluator-only TSV keyed by hidden gene_id |
| `--stratum-column` | no | `string` | `[]` | — | repeatable marginal/joint recall column from --event-strata |

## `ploidypatch benchmark score-abstention`

ploidypatch benchmark score-abstention

```text
usage: ploidypatch benchmark score-abstention [-h] --perturbed-gff
                                              PERTURBED_GFF --candidate-gff
                                              CANDIDATE_GFF --mask-truth
                                              MASK_TRUTH
                                              [--include-event-details]
```

| Argument | Required | Type | Default | Choices | Description |
|---|---:|---|---|---|---|
| `-h`<br>`--help` | no | `string` | `—` | — | show this help message and exit |
| `--perturbed-gff` | yes | `string` | `None` | — | — |
| `--candidate-gff` | yes | `string` | `None` | — | — |
| `--mask-truth` | yes | `string` | `None` | — | — |
| `--include-event-details` | no | `flag` | `false` | — | include evaluator-only event identifiers and coordinates |

## `ploidypatch benchmark score-copy-ranking`

ploidypatch benchmark score-copy-ranking

```text
usage: ploidypatch benchmark score-copy-ranking [-h] --scores SCORES
                                                --labeled-features
                                                LABELED_FEATURES --output-json
                                                OUTPUT_JSON
```

| Argument | Required | Type | Default | Choices | Description |
|---|---:|---|---|---|---|
| `-h`<br>`--help` | no | `string` | `—` | — | show this help message and exit |
| `--scores` | yes | `string` | `None` | — | — |
| `--labeled-features` | yes | `string` | `None` | — | — |
| `--output-json` | yes | `string` | `None` | — | — |

## `ploidypatch benchmark score-localization`

ploidypatch benchmark score-localization

```text
usage: ploidypatch benchmark score-localization [-h] --gaps GAPS --selection
                                                SELECTION --truth TRUTH
                                                [--include-event-details]
```

| Argument | Required | Type | Default | Choices | Description |
|---|---:|---|---|---|---|
| `-h`<br>`--help` | no | `string` | `—` | — | show this help message and exit |
| `--gaps` | yes | `string` | `None` | — | repeatable upstream gap TSV |
| `--selection` | yes | `string` | `None` | — | — |
| `--truth` | yes | `string` | `None` | — | — |
| `--include-event-details` | no | `flag` | `false` | — | — |

## `ploidypatch benchmark score-structure-hypotheses`

ploidypatch benchmark score-structure-hypotheses

```text
usage: ploidypatch benchmark score-structure-hypotheses [-h] --source-gff
                                                        SOURCE_GFF
                                                        --perturbed-gff
                                                        PERTURBED_GFF
                                                        --candidate-gff
                                                        CANDIDATE_GFF
                                                        --hypotheses-tsv
                                                        HYPOTHESES_TSV --truth
                                                        TRUTH
                                                        [--control-hypotheses-tsv CONTROL_HYPOTHESES_TSV]
                                                        [--include-event-details]
```

| Argument | Required | Type | Default | Choices | Description |
|---|---:|---|---|---|---|
| `-h`<br>`--help` | no | `string` | `—` | — | show this help message and exit |
| `--source-gff` | yes | `string` | `None` | — | — |
| `--perturbed-gff` | yes | `string` | `None` | — | — |
| `--candidate-gff` | yes | `string` | `None` | — | — |
| `--hypotheses-tsv` | yes | `string` | `None` | — | — |
| `--truth` | yes | `string` | `None` | — | — |
| `--control-hypotheses-tsv` | no | `string` | `None` | — | — |
| `--include-event-details` | no | `flag` | `false` | — | — |

## `ploidypatch benchmark summarize-catalog`

ploidypatch benchmark summarize-catalog

```text
usage: ploidypatch benchmark summarize-catalog [-h] --catalog CATALOG --output
                                               OUTPUT --column COLUMN
                                               [--cross CROSS]
```

| Argument | Required | Type | Default | Choices | Description |
|---|---:|---|---|---|---|
| `-h`<br>`--help` | no | `string` | `—` | — | show this help message and exit |
| `--catalog` | yes | `string` | `None` | — | — |
| `--output` | yes | `string` | `None` | — | — |
| `--column` | yes | `string` | `None` | — | repeatable column to summarize; empty values are counted explicitly |
| `--cross` | no | `string` | `[]` | — | repeatable COLUMN_A,COLUMN_B joint count |

## `ploidypatch benchmark summarize-mask-selection`

ploidypatch benchmark summarize-mask-selection

```text
usage: ploidypatch benchmark summarize-mask-selection [-h] --mask-truth
                                                      MASK_TRUTH
                                                      --hidden-truth
                                                      HIDDEN_TRUTH
                                                      --strata-tsv STRATA_TSV
                                                      --column COLUMN
```

| Argument | Required | Type | Default | Choices | Description |
|---|---:|---|---|---|---|
| `-h`<br>`--help` | no | `string` | `—` | — | show this help message and exit |
| `--mask-truth` | yes | `string` | `None` | — | — |
| `--hidden-truth` | yes | `string` | `None` | — | — |
| `--strata-tsv` | yes | `string` | `None` | — | — |
| `--column` | yes | `string` | `None` | — | — |

## `ploidypatch evidence aggregate-junctions`

ploidypatch evidence aggregate-junctions

```text
usage: ploidypatch evidence aggregate-junctions [-h] --input-dir INPUT_DIR
                                                --primary-sample
                                                PRIMARY_SAMPLE --output-tsv
                                                OUTPUT_TSV
                                                [--min-reads-per-sample MIN_READS_PER_SAMPLE]
                                                [--min-supporting-samples MIN_SUPPORTING_SAMPLES]
                                                [--sample-groups SAMPLE_GROUPS]
                                                [--min-samples-per-group MIN_SAMPLES_PER_GROUP]
                                                [--min-secondary-groups MIN_SECONDARY_GROUPS]
```

| Argument | Required | Type | Default | Choices | Description |
|---|---:|---|---|---|---|
| `-h`<br>`--help` | no | `string` | `—` | — | show this help message and exit |
| `--input-dir` | yes | `string` | `None` | — | repeatable directory containing *.junctions.tsv plus manifests |
| `--primary-sample` | yes | `string` | `None` | — | — |
| `--output-tsv` | yes | `string` | `None` | — | — |
| `--min-reads-per-sample` | no | `int` | `2` | — | — |
| `--min-supporting-samples` | no | `int` | `2` | — | — |
| `--sample-groups` | no | `string` | `None` | — | optional explicit sample/filename_stem_group TSV |
| `--min-samples-per-group` | no | `int` | `2` | — | — |
| `--min-secondary-groups` | no | `int` | `2` | — | — |

## `ploidypatch evidence aggregate-reference-anchored`

ploidypatch evidence aggregate-reference-anchored

```text
usage: ploidypatch evidence aggregate-reference-anchored [-h]
                                                         --candidate-source-provenance
                                                         CANDIDATE_SOURCE_PROVENANCE
                                                         --accepted-wgd-pairs
                                                         ACCEPTED_WGD_PAIRS
                                                         --homeolog-base-evidence
                                                         HOMEOLOG_BASE_EVIDENCE
                                                         --output-dir
                                                         OUTPUT_DIR
                                                         [--evidence-type EVIDENCE_TYPE]
```

| Argument | Required | Type | Default | Choices | Description |
|---|---:|---|---|---|---|
| `-h`<br>`--help` | no | `string` | `—` | — | show this help message and exit |
| `--candidate-source-provenance` | yes | `string` | `None` | — | — |
| `--accepted-wgd-pairs` | yes | `string` | `None` | — | repeatable REFERENCE=PATH accepted directed reciprocal pair table |
| `--homeolog-base-evidence` | yes | `string` | `None` | — | repeatable REFERENCE=PATH arm/homeolog/base evidence table |
| `--output-dir` | yes | `string` | `None` | — | — |
| `--evidence-type` | no | `string` | `reference_anchored_projection` | — | — |

## `ploidypatch evidence annotate-natural-assembly`

ploidypatch evidence annotate-natural-assembly

```text
usage: ploidypatch evidence annotate-natural-assembly [-h] --validation-tsv
                                                      VALIDATION_TSV --genome
                                                      GENOME --fai FAI
                                                      --output-tsv OUTPUT_TSV
                                                      [--flank-bp FLANK_BP]
                                                      [--max-ambiguous-fraction MAX_AMBIGUOUS_FRACTION]
```

| Argument | Required | Type | Default | Choices | Description |
|---|---:|---|---|---|---|
| `-h`<br>`--help` | no | `string` | `—` | — | show this help message and exit |
| `--validation-tsv` | yes | `string` | `None` | — | — |
| `--genome` | yes | `string` | `None` | — | — |
| `--fai` | yes | `string` | `None` | — | — |
| `--output-tsv` | yes | `string` | `None` | — | — |
| `--flank-bp` | no | `int` | `5000` | — | — |
| `--max-ambiguous-fraction` | no | `float` | `0.0` | — | — |

## `ploidypatch evidence annotate-pav-wgdi`

ploidypatch evidence annotate-pav-wgdi

```text
usage: ploidypatch evidence annotate-pav-wgdi [-h] --protein-screen
                                              PROTEIN_SCREEN --wgdi-genes
                                              WGDI_GENES --output-tsv
                                              OUTPUT_TSV
```

| Argument | Required | Type | Default | Choices | Description |
|---|---:|---|---|---|---|
| `-h`<br>`--help` | no | `string` | `—` | — | show this help message and exit |
| `--protein-screen` | yes | `string` | `None` | — | — |
| `--wgdi-genes` | yes | `string` | `None` | — | — |
| `--output-tsv` | yes | `string` | `None` | — | — |

## `ploidypatch evidence apply-conflict-winner-guard`

ploidypatch evidence apply-conflict-winner-guard

```text
usage: ploidypatch evidence apply-conflict-winner-guard [-h] --v03-scores
                                                        V03_SCORES
                                                        --pool-decisions
                                                        POOL_DECISIONS
                                                        --pool-manifest
                                                        POOL_MANIFEST
                                                        --output-tsv
                                                        OUTPUT_TSV
```

| Argument | Required | Type | Default | Choices | Description |
|---|---:|---|---|---|---|
| `-h`<br>`--help` | no | `string` | `—` | — | show this help message and exit |
| `--v03-scores` | yes | `string` | `None` | — | — |
| `--pool-decisions` | yes | `string` | `None` | — | — |
| `--pool-manifest` | yes | `string` | `None` | — | — |
| `--output-tsv` | yes | `string` | `None` | — | — |

## `ploidypatch evidence audit-natural-candidates`

ploidypatch evidence audit-natural-candidates

```text
usage: ploidypatch evidence audit-natural-candidates [-h] --candidate-gff
                                                     CANDIDATE_GFF --base-gff
                                                     BASE_GFF --genome-fasta
                                                     GENOME_FASTA
                                                     --review-rankings
                                                     REVIEW_RANKINGS
                                                     [--isoseq-evidence ISOSEQ_EVIDENCE]
                                                     [--self-map-paf SELF_MAP_PAF]
                                                     [--repeat-gff REPEAT_GFF]
                                                     [--repeat-seqid-map REPEAT_SEQID_MAP]
                                                     [--repeat-flank-bp REPEAT_FLANK_BP]
                                                     [--minimum-full-length-read-support MINIMUM_FULL_LENGTH_READ_SUPPORT]
                                                     [--review-budget REVIEW_BUDGET]
                                                     [--minimum-query-coverage MINIMUM_QUERY_COVERAGE]
                                                     [--minimum-identity MINIMUM_IDENTITY]
                                                     [--near-equal-score-fraction NEAR_EQUAL_SCORE_FRACTION]
                                                     --output-tsv OUTPUT_TSV
                                                     --output-summary
                                                     OUTPUT_SUMMARY
```

| Argument | Required | Type | Default | Choices | Description |
|---|---:|---|---|---|---|
| `-h`<br>`--help` | no | `string` | `—` | — | show this help message and exit |
| `--candidate-gff` | yes | `string` | `None` | — | — |
| `--base-gff` | yes | `string` | `None` | — | — |
| `--genome-fasta` | yes | `string` | `None` | — | — |
| `--review-rankings` | yes | `string` | `None` | — | — |
| `--isoseq-evidence` | no | `string` | `None` | — | — |
| `--self-map-paf` | no | `string` | `None` | — | — |
| `--repeat-gff` | no | `string` | `None` | — | — |
| `--repeat-seqid-map` | no | `string` | `None` | — | — |
| `--repeat-flank-bp` | no | `int` | `2000` | — | — |
| `--minimum-full-length-read-support` | no | `int` | `1` | — | — |
| `--review-budget` | no | `int` | `None` | — | — |
| `--minimum-query-coverage` | no | `float` | `0.9` | — | — |
| `--minimum-identity` | no | `float` | `0.98` | — | — |
| `--near-equal-score-fraction` | no | `float` | `0.95` | — | — |
| `--output-tsv` | yes | `string` | `None` | — | — |
| `--output-summary` | yes | `string` | `None` | — | — |

## `ploidypatch evidence bootstrap-isoseq-review-yield`

ploidypatch evidence bootstrap-isoseq-review-yield

```text
usage: ploidypatch evidence bootstrap-isoseq-review-yield [-h] --evidence
                                                          EVIDENCE
                                                          --review-rankings
                                                          REVIEW_RANKINGS
                                                          [--review-budget REVIEW_BUDGET]
                                                          [--replicates REPLICATES]
                                                          [--seed SEED]
                                                          [--alpha ALPHA]
                                                          [--comparator-estimator COMPARATOR_ESTIMATOR]
                                                          [--primary-estimator PRIMARY_ESTIMATOR]
                                                          --output-json
                                                          OUTPUT_JSON
```

| Argument | Required | Type | Default | Choices | Description |
|---|---:|---|---|---|---|
| `-h`<br>`--help` | no | `string` | `—` | — | show this help message and exit |
| `--evidence` | yes | `string` | `None` | — | — |
| `--review-rankings` | yes | `string` | `None` | — | — |
| `--review-budget` | no | `int` | `None` | — | — |
| `--replicates` | no | `int` | `20000` | — | — |
| `--seed` | no | `int` | `20260808` | — | — |
| `--alpha` | no | `float` | `0.05` | — | — |
| `--comparator-estimator` | no | `string` | `baseline` | — | — |
| `--primary-estimator` | no | `string` | `topology` | — | — |
| `--output-json` | yes | `string` | `None` | — | — |

## `ploidypatch evidence build-copy-features`

ploidypatch evidence build-copy-features

```text
usage: ploidypatch evidence build-copy-features [-h] --consensus-decisions
                                                CONSENSUS_DECISIONS
                                                --method-decisions
                                                METHOD_DECISIONS
                                                [--wgd-selection WGD_SELECTION]
                                                --output-tsv OUTPUT_TSV
```

| Argument | Required | Type | Default | Choices | Description |
|---|---:|---|---|---|---|
| `-h`<br>`--help` | no | `string` | `—` | — | show this help message and exit |
| `--consensus-decisions` | yes | `string` | `None` | — | — |
| `--method-decisions` | yes | `string` | `None` | — | repeatable METHOD=PATH decision table |
| `--wgd-selection` | no | `string` | `None` | — | — |
| `--output-tsv` | yes | `string` | `None` | — | — |

## `ploidypatch evidence build-fixed-target-backbone`

ploidypatch evidence build-fixed-target-backbone

```text
usage: ploidypatch evidence build-fixed-target-backbone [-h] --base-gff
                                                        BASE_GFF
                                                        --query-wgdi-gff
                                                        QUERY_WGDI_GFF
                                                        --base-only-collinearity
                                                        BASE_ONLY_COLLINEARITY
                                                        --primary-seqid-table
                                                        PRIMARY_SEQID_TABLE
                                                        [--min-block-pairs MIN_BLOCK_PAIRS]
                                                        --output-dir
                                                        OUTPUT_DIR
```

| Argument | Required | Type | Default | Choices | Description |
|---|---:|---|---|---|---|
| `-h`<br>`--help` | no | `string` | `—` | — | show this help message and exit |
| `--base-gff` | yes | `string` | `None` | — | — |
| `--query-wgdi-gff` | yes | `string` | `None` | — | — |
| `--base-only-collinearity` | yes | `string` | `None` | — | — |
| `--primary-seqid-table` | yes | `string` | `None` | — | — |
| `--min-block-pairs` | no | `int` | `20` | — | — |
| `--output-dir` | yes | `string` | `None` | — | — |

## `ploidypatch evidence build-homeolog-topology-features`

ploidypatch evidence build-homeolog-topology-features

```text
usage: ploidypatch evidence build-homeolog-topology-features
       [-h] --copy-features COPY_FEATURES --wgd-selection WGD_SELECTION
       --candidate-gff CANDIDATE_GFF --base-gff BASE_GFF --output-tsv
       OUTPUT_TSV
```

| Argument | Required | Type | Default | Choices | Description |
|---|---:|---|---|---|---|
| `-h`<br>`--help` | no | `string` | `—` | — | show this help message and exit |
| `--copy-features` | yes | `string` | `None` | — | — |
| `--wgd-selection` | yes | `string` | `None` | — | — |
| `--candidate-gff` | yes | `string` | `None` | — | — |
| `--base-gff` | yes | `string` | `None` | — | — |
| `--output-tsv` | yes | `string` | `None` | — | — |

## `ploidypatch evidence catalog-pav-genes`

ploidypatch evidence catalog-pav-genes

```text
usage: ploidypatch evidence catalog-pav-genes [-h] --gff GFF --deletions
                                              DELETIONS --output-tsv
                                              OUTPUT_TSV
                                              [--min-gene-coverage MIN_GENE_COVERAGE]
                                              [--min-cds-coverage MIN_CDS_COVERAGE]
```

| Argument | Required | Type | Default | Choices | Description |
|---|---:|---|---|---|---|
| `-h`<br>`--help` | no | `string` | `—` | — | show this help message and exit |
| `--gff` | yes | `string` | `None` | — | — |
| `--deletions` | yes | `string` | `None` | — | — |
| `--output-tsv` | yes | `string` | `None` | — | — |
| `--min-gene-coverage` | no | `float` | `0.95` | — | — |
| `--min-cds-coverage` | no | `float` | `1.0` | — | — |

## `ploidypatch evidence collapse-pav-events`

ploidypatch evidence collapse-pav-events

```text
usage: ploidypatch evidence collapse-pav-events [-h] --annotated-genes
                                                ANNOTATED_GENES --output-tsv
                                                OUTPUT_TSV
```

| Argument | Required | Type | Default | Choices | Description |
|---|---:|---|---|---|---|
| `-h`<br>`--help` | no | `string` | `—` | — | show this help message and exit |
| `--annotated-genes` | yes | `string` | `None` | — | — |
| `--output-tsv` | yes | `string` | `None` | — | — |

## `ploidypatch evidence combine-chromosome-pafs`

ploidypatch evidence combine-chromosome-pafs

```text
usage: ploidypatch evidence combine-chromosome-pafs [-h] --input-dir INPUT_DIR
                                                    --chromosome-pairs
                                                    CHROMOSOME_PAIRS
                                                    --output-paf OUTPUT_PAF
```

| Argument | Required | Type | Default | Choices | Description |
|---|---:|---|---|---|---|
| `-h`<br>`--help` | no | `string` | `—` | — | show this help message and exit |
| `--input-dir` | yes | `string` | `None` | — | — |
| `--chromosome-pairs` | yes | `string` | `None` | — | — |
| `--output-paf` | yes | `string` | `None` | — | — |

## `ploidypatch evidence discover-natural`

ploidypatch evidence discover-natural

```text
usage: ploidypatch evidence discover-natural [-h] --target-gff TARGET_GFF
                                             --miniprot-gff MINIPROT_GFF
                                             --protein-map PROTEIN_MAP
                                             --output-tsv OUTPUT_TSV
                                             [--min-identity MIN_IDENTITY]
                                             [--min-query-coverage MIN_QUERY_COVERAGE]
                                             [--max-existing-cds-overlap MAX_EXISTING_CDS_OVERLAP]
                                             [--min-boundary-extension-bp MIN_BOUNDARY_EXTENSION_BP]
                                             [--near-best-score-fraction NEAR_BEST_SCORE_FRACTION]
```

| Argument | Required | Type | Default | Choices | Description |
|---|---:|---|---|---|---|
| `-h`<br>`--help` | no | `string` | `—` | — | show this help message and exit |
| `--target-gff` | yes | `string` | `None` | — | — |
| `--miniprot-gff` | yes | `string` | `None` | — | — |
| `--protein-map` | yes | `string` | `None` | — | — |
| `--output-tsv` | yes | `string` | `None` | — | — |
| `--min-identity` | no | `float` | `0.8` | — | — |
| `--min-query-coverage` | no | `float` | `0.8` | — | — |
| `--max-existing-cds-overlap` | no | `float` | `0.1` | — | — |
| `--min-boundary-extension-bp` | no | `int` | `30` | — | — |
| `--near-best-score-fraction` | no | `float` | `0.95` | — | — |

## `ploidypatch evidence evaluate-species-applicability`

ploidypatch evidence evaluate-species-applicability

```text
usage: ploidypatch evidence evaluate-species-applicability [-h] --metrics
                                                           METRICS --policy
                                                           POLICY
                                                           --output-json
                                                           OUTPUT_JSON
```

| Argument | Required | Type | Default | Choices | Description |
|---|---:|---|---|---|---|
| `-h`<br>`--help` | no | `string` | `—` | — | show this help message and exit |
| `--metrics` | yes | `string` | `None` | — | — |
| `--policy` | yes | `string` | `None` | — | — |
| `--output-json` | yes | `string` | `None` | — | — |

## `ploidypatch evidence export-natural-candidate-cds`

ploidypatch evidence export-natural-candidate-cds

```text
usage: ploidypatch evidence export-natural-candidate-cds [-h] --candidate-gff
                                                         CANDIDATE_GFF
                                                         --genome-fasta
                                                         GENOME_FASTA
                                                         --output-fasta
                                                         OUTPUT_FASTA
```

| Argument | Required | Type | Default | Choices | Description |
|---|---:|---|---|---|---|
| `-h`<br>`--help` | no | `string` | `—` | — | show this help message and exit |
| `--candidate-gff` | yes | `string` | `None` | — | — |
| `--genome-fasta` | yes | `string` | `None` | — | — |
| `--output-fasta` | yes | `string` | `None` | — | — |

## `ploidypatch evidence extract-bam-junctions`

ploidypatch evidence extract-bam-junctions

```text
usage: ploidypatch evidence extract-bam-junctions [-h] --bam BAM --sample
                                                  SAMPLE --samtools SAMTOOLS
                                                  --output-tsv OUTPUT_TSV
                                                  [--threads THREADS]
                                                  [--min-mapq MIN_MAPQ]
```

| Argument | Required | Type | Default | Choices | Description |
|---|---:|---|---|---|---|
| `-h`<br>`--help` | no | `string` | `—` | — | show this help message and exit |
| `--bam` | yes | `string` | `None` | — | — |
| `--sample` | yes | `string` | `None` | — | — |
| `--samtools` | yes | `string` | `None` | — | — |
| `--output-tsv` | yes | `string` | `None` | — | — |
| `--threads` | no | `int` | `4` | — | — |
| `--min-mapq` | no | `int` | `20` | — | — |

## `ploidypatch evidence extract-paf-deletions`

ploidypatch evidence extract-paf-deletions

```text
usage: ploidypatch evidence extract-paf-deletions [-h] --paf PAF --output-tsv
                                                  OUTPUT_TSV
                                                  [--min-deletion-bp MIN_DELETION_BP]
                                                  [--min-flanking-anchor-bp MIN_FLANKING_ANCHOR_BP]
                                                  [--min-mapq MIN_MAPQ]
                                                  [--chromosome-pairs CHROMOSOME_PAIRS]
                                                  [--query-genome QUERY_GENOME]
                                                  [--query-fai QUERY_FAI]
                                                  [--query-flank-bp QUERY_FLANK_BP]
                                                  [--max-query-flank-ambiguous-fraction MAX_QUERY_FLANK_AMBIGUOUS_FRACTION]
                                                  [--retain-uncertain-query-flanks]
                                                  [--allow-secondary]
```

| Argument | Required | Type | Default | Choices | Description |
|---|---:|---|---|---|---|
| `-h`<br>`--help` | no | `string` | `—` | — | show this help message and exit |
| `--paf` | yes | `string` | `None` | — | — |
| `--output-tsv` | yes | `string` | `None` | — | — |
| `--min-deletion-bp` | no | `int` | `100` | — | — |
| `--min-flanking-anchor-bp` | no | `int` | `1000` | — | — |
| `--min-mapq` | no | `int` | `20` | — | — |
| `--chromosome-pairs` | no | `string` | `None` | — | optional query_seqid/target_seqid TSV for homologous-pair filtering |
| `--query-genome` | no | `string` | `None` | — | uncompressed query FASTA for assembly-ambiguity checks at breakpoints |
| `--query-fai` | no | `string` | `None` | — | samtools FAI for --query-genome |
| `--query-flank-bp` | no | `int` | `5000` | — | — |
| `--max-query-flank-ambiguous-fraction` | no | `float` | `0.0` | — | — |
| `--retain-uncertain-query-flanks` | no | `flag` | `false` | — | report rather than exclude truncated or ambiguous query flanks |
| `--allow-secondary` | no | `flag` | `false` | — | retain non-primary PAF rows |

## `ploidypatch evidence filter-candidate-query-paf`

ploidypatch evidence filter-candidate-query-paf

```text
usage: ploidypatch evidence filter-candidate-query-paf [-h] --candidate-gff
                                                       CANDIDATE_GFF
                                                       --paf-input PAF_INPUT
                                                       [--alignment-strand-source {query_orientation,minimap2_ts}]
                                                       --output-paf OUTPUT_PAF
                                                       --output-counts
                                                       OUTPUT_COUNTS
                                                       --output-summary
                                                       OUTPUT_SUMMARY
                                                       --output-manifest
                                                       OUTPUT_MANIFEST
```

| Argument | Required | Type | Default | Choices | Description |
|---|---:|---|---|---|---|
| `-h`<br>`--help` | no | `string` | `—` | — | show this help message and exit |
| `--candidate-gff` | yes | `string` | `None` | — | — |
| `--paf-input` | yes | `string` | `None` | — | repeatable ACCESSION=PATH PAF input |
| `--alignment-strand-source` | no | `string` | `query_orientation` | `query_orientation`, `minimap2_ts` | — |
| `--output-paf` | yes | `string` | `None` | — | — |
| `--output-counts` | yes | `string` | `None` | — | — |
| `--output-summary` | yes | `string` | `None` | — | — |
| `--output-manifest` | yes | `string` | `None` | — | — |

## `ploidypatch evidence freeze-homeolog-review-rankings`

ploidypatch evidence freeze-homeolog-review-rankings

```text
usage: ploidypatch evidence freeze-homeolog-review-rankings
       [-h] --scores SCORES [--review-budget REVIEW_BUDGET] --output-tsv
       OUTPUT_TSV
```

| Argument | Required | Type | Default | Choices | Description |
|---|---:|---|---|---|---|
| `-h`<br>`--help` | no | `string` | `—` | — | show this help message and exit |
| `--scores` | yes | `string` | `None` | — | — |
| `--review-budget` | no | `int` | `None` | — | repeatable top-K review budget; defaults to 25,50,100,200 |
| `--output-tsv` | yes | `string` | `None` | — | — |

## `ploidypatch evidence infer-homeolog-pairs`

ploidypatch evidence infer-homeolog-pairs

```text
usage: ploidypatch evidence infer-homeolog-pairs [-h] --gene-evidence
                                                 GENE_EVIDENCE --collinearity
                                                 COLLINEARITY
                                                 [--source-group-map SOURCE_GROUP_MAP]
                                                 --wgd-event WGD_EVENT
                                                 --subgenome SUBGENOME
                                                 [--min-support-group-count MIN_SUPPORT_GROUP_COUNT]
                                                 [--allow-multiple-partners]
                                                 --output-pairs OUTPUT_PAIRS
                                                 --decisions-tsv DECISIONS_TSV
```

| Argument | Required | Type | Default | Choices | Description |
|---|---:|---|---|---|---|
| `-h`<br>`--help` | no | `string` | `—` | — | show this help message and exit |
| `--gene-evidence` | yes | `string` | `None` | — | — |
| `--collinearity` | yes | `string` | `None` | — | repeatable SOURCE=PATH value |
| `--source-group-map` | no | `string` | `None` | — | — |
| `--wgd-event` | yes | `string` | `None` | — | — |
| `--subgenome` | yes | `string` | `None` | — | repeat exactly twice in the desired A/B output order |
| `--min-support-group-count` | no | `int` | `2` | — | — |
| `--allow-multiple-partners` | no | `flag` | `false` | — | disable the default reciprocal-unique partner gate |
| `--output-pairs` | yes | `string` | `None` | — | — |
| `--decisions-tsv` | yes | `string` | `None` | — | — |

## `ploidypatch evidence infer-outgroup-duplicated-pairs`

ploidypatch evidence infer-outgroup-duplicated-pairs

```text
usage: ploidypatch evidence infer-outgroup-duplicated-pairs
       [-h] --query-wgdi-gff QUERY_WGDI_GFF [--source-gff SOURCE_GFF]
       --collinearity COLLINEARITY [--source-group-map SOURCE_GROUP_MAP]
       --wgd-event WGD_EVENT
       [--min-support-group-count MIN_SUPPORT_GROUP_COUNT]
       [--min-block-pairs MIN_BLOCK_PAIRS] [--allow-same-seqid]
       [--allow-multiple-partners] --output-pairs OUTPUT_PAIRS --decisions-tsv
       DECISIONS_TSV
```

| Argument | Required | Type | Default | Choices | Description |
|---|---:|---|---|---|---|
| `-h`<br>`--help` | no | `string` | `—` | — | show this help message and exit |
| `--query-wgdi-gff` | yes | `string` | `None` | — | — |
| `--source-gff` | no | `string` | `None` | — | map normalized WGDI IDs uniquely to operable source gene feature IDs |
| `--collinearity` | yes | `string` | `None` | — | repeatable SOURCE=PATH value |
| `--source-group-map` | no | `string` | `None` | — | — |
| `--wgd-event` | yes | `string` | `None` | — | — |
| `--min-support-group-count` | no | `int` | `2` | — | — |
| `--min-block-pairs` | no | `int` | `20` | — | — |
| `--allow-same-seqid` | no | `flag` | `false` | — | — |
| `--allow-multiple-partners` | no | `flag` | `false` | — | — |
| `--output-pairs` | yes | `string` | `None` | — | — |
| `--decisions-tsv` | yes | `string` | `None` | — | — |

## `ploidypatch evidence infer-self-wgd-pairs`

ploidypatch evidence infer-self-wgd-pairs

```text
usage: ploidypatch evidence infer-self-wgd-pairs [-h] --query-wgdi-gff
                                                 QUERY_WGDI_GFF --collinearity
                                                 COLLINEARITY
                                                 [--source-gff SOURCE_GFF]
                                                 --wgd-event WGD_EVENT
                                                 [--min-block-pairs MIN_BLOCK_PAIRS]
                                                 [--allow-same-seqid]
                                                 [--allow-multiple-partners]
                                                 --output-pairs OUTPUT_PAIRS
                                                 --decisions-tsv DECISIONS_TSV
```

| Argument | Required | Type | Default | Choices | Description |
|---|---:|---|---|---|---|
| `-h`<br>`--help` | no | `string` | `—` | — | show this help message and exit |
| `--query-wgdi-gff` | yes | `string` | `None` | — | — |
| `--collinearity` | yes | `string` | `None` | — | — |
| `--source-gff` | no | `string` | `None` | — | map WGDI IDs uniquely to operable source gene feature IDs |
| `--wgd-event` | yes | `string` | `None` | — | — |
| `--min-block-pairs` | no | `int` | `20` | — | — |
| `--allow-same-seqid` | no | `flag` | `false` | — | — |
| `--allow-multiple-partners` | no | `flag` | `false` | — | — |
| `--output-pairs` | yes | `string` | `None` | — | — |
| `--decisions-tsv` | yes | `string` | `None` | — | — |

## `ploidypatch evidence infer-synteny-gaps`

ploidypatch evidence infer-synteny-gaps

```text
usage: ploidypatch evidence infer-synteny-gaps [-h] --query-wgdi-gff
                                               QUERY_WGDI_GFF
                                               --target-wgdi-gff
                                               TARGET_WGDI_GFF --collinearity
                                               COLLINEARITY --source-label
                                               SOURCE_LABEL --output-tsv
                                               OUTPUT_TSV
                                               [--expected-chromosome-pairs EXPECTED_CHROMOSOME_PAIRS]
                                               [--max-query-intervening-genes MAX_QUERY_INTERVENING_GENES]
                                               [--min-target-excess-genes MIN_TARGET_EXCESS_GENES]
                                               [--max-target-gap-genes MAX_TARGET_GAP_GENES]
                                               [--max-query-locus-bp MAX_QUERY_LOCUS_BP]
```

| Argument | Required | Type | Default | Choices | Description |
|---|---:|---|---|---|---|
| `-h`<br>`--help` | no | `string` | `—` | — | show this help message and exit |
| `--query-wgdi-gff` | yes | `string` | `None` | — | — |
| `--target-wgdi-gff` | yes | `string` | `None` | — | — |
| `--collinearity` | yes | `string` | `None` | — | — |
| `--source-label` | yes | `string` | `None` | — | — |
| `--output-tsv` | yes | `string` | `None` | — | — |
| `--expected-chromosome-pairs` | no | `string` | `None` | — | optional query/target/source TSV for subgenome-aware filtering |
| `--max-query-intervening-genes` | no | `int` | `0` | — | — |
| `--min-target-excess-genes` | no | `int` | `1` | — | — |
| `--max-target-gap-genes` | no | `int` | `5` | — | — |
| `--max-query-locus-bp` | no | `int` | `500000` | — | — |

## `ploidypatch evidence intersect-copy-pair-evidence`

ploidypatch evidence intersect-copy-pair-evidence

```text
usage: ploidypatch evidence intersect-copy-pair-evidence [-h] --pairs PAIRS
                                                         --pair-set-label
                                                         PAIR_SET_LABEL
                                                         [--allow-multiple-partners]
                                                         --output-pairs
                                                         OUTPUT_PAIRS
                                                         --decisions-tsv
                                                         DECISIONS_TSV
```

| Argument | Required | Type | Default | Choices | Description |
|---|---:|---|---|---|---|
| `-h`<br>`--help` | no | `string` | `—` | — | show this help message and exit |
| `--pairs` | yes | `string` | `None` | — | repeatable LABEL=PATH pair table |
| `--pair-set-label` | yes | `string` | `None` | — | — |
| `--allow-multiple-partners` | no | `flag` | `false` | — | — |
| `--output-pairs` | yes | `string` | `None` | — | — |
| `--decisions-tsv` | yes | `string` | `None` | — | — |

## `ploidypatch evidence join-isoseq-review-rankings`

ploidypatch evidence join-isoseq-review-rankings

```text
usage: ploidypatch evidence join-isoseq-review-rankings [-h] --evidence
                                                        EVIDENCE
                                                        --review-rankings
                                                        REVIEW_RANKINGS
                                                        [--review-budget REVIEW_BUDGET]
                                                        [--comparator-estimator COMPARATOR_ESTIMATOR]
                                                        [--primary-estimator PRIMARY_ESTIMATOR]
                                                        --output-tsv
                                                        OUTPUT_TSV
                                                        --output-summary
                                                        OUTPUT_SUMMARY
```

| Argument | Required | Type | Default | Choices | Description |
|---|---:|---|---|---|---|
| `-h`<br>`--help` | no | `string` | `—` | — | show this help message and exit |
| `--evidence` | yes | `string` | `None` | — | — |
| `--review-rankings` | yes | `string` | `None` | — | — |
| `--review-budget` | no | `int` | `None` | — | — |
| `--comparator-estimator` | no | `string` | `baseline` | — | — |
| `--primary-estimator` | no | `string` | `topology` | — | — |
| `--output-tsv` | yes | `string` | `None` | — | — |
| `--output-summary` | yes | `string` | `None` | — | — |

## `ploidypatch evidence prepare-b73-isoseq`

ploidypatch evidence prepare-b73-isoseq

```text
usage: ploidypatch evidence prepare-b73-isoseq [-h] --fasta FASTA --counts-csv
                                               COUNTS_CSV
                                               [--minimum-b73-reads MINIMUM_B73_READS]
                                               --output-fasta OUTPUT_FASTA
                                               --output-counts OUTPUT_COUNTS
```

| Argument | Required | Type | Default | Choices | Description |
|---|---:|---|---|---|---|
| `-h`<br>`--help` | no | `string` | `—` | — | show this help message and exit |
| `--fasta` | yes | `string` | `None` | — | — |
| `--counts-csv` | yes | `string` | `None` | — | — |
| `--minimum-b73-reads` | no | `int` | `2` | — | — |
| `--output-fasta` | yes | `string` | `None` | — | — |
| `--output-counts` | yes | `string` | `None` | — | — |

## `ploidypatch evidence prepare-natural-graph`

ploidypatch evidence prepare-natural-graph

```text
usage: ploidypatch evidence prepare-natural-graph [-h] --validation-tsv
                                                  VALIDATION_TSV
                                                  --output-candidates
                                                  OUTPUT_CANDIDATES
                                                  --output-evidence
                                                  OUTPUT_EVIDENCE
```

| Argument | Required | Type | Default | Choices | Description |
|---|---:|---|---|---|---|
| `-h`<br>`--help` | no | `string` | `—` | — | show this help message and exit |
| `--validation-tsv` | yes | `string` | `None` | — | — |
| `--output-candidates` | yes | `string` | `None` | — | — |
| `--output-evidence` | yes | `string` | `None` | — | — |

## `ploidypatch evidence prepare-wgdi`

ploidypatch evidence prepare-wgdi

```text
usage: ploidypatch evidence prepare-wgdi [-h] --gff GFF --protein PROTEIN
                                         --fai FAI --output-dir OUTPUT_DIR
                                         --prefix PREFIX
                                         [--min-genes-per-seqid MIN_GENES_PER_SEQID]
                                         [--primary-chromosomes-only]
```

| Argument | Required | Type | Default | Choices | Description |
|---|---:|---|---|---|---|
| `-h`<br>`--help` | no | `string` | `—` | — | show this help message and exit |
| `--gff` | yes | `string` | `None` | — | — |
| `--protein` | yes | `string` | `None` | — | — |
| `--fai` | yes | `string` | `None` | — | — |
| `--output-dir` | yes | `string` | `None` | — | — |
| `--prefix` | yes | `string` | `None` | — | — |
| `--min-genes-per-seqid` | no | `int` | `1` | — | — |
| `--primary-chromosomes-only` | no | `flag` | `false` | — | retain NCBI region features marked genome=chromosome |

## `ploidypatch evidence project-fixed-backbone`

ploidypatch evidence project-fixed-backbone

```text
usage: ploidypatch evidence project-fixed-backbone [-h] --backbone-dir
                                                   BACKBONE_DIR
                                                   --candidate-catalog
                                                   CANDIDATE_CATALOG
                                                   --fixed-base-hits
                                                   FIXED_BASE_HITS
                                                   --fixed-base-hits-manifest
                                                   FIXED_BASE_HITS_MANIFEST
                                                   [--minimum-query-coverage MINIMUM_QUERY_COVERAGE]
                                                   [--maximum-evalue MAXIMUM_EVALUE]
                                                   --output-dir OUTPUT_DIR
```

| Argument | Required | Type | Default | Choices | Description |
|---|---:|---|---|---|---|
| `-h`<br>`--help` | no | `string` | `—` | — | show this help message and exit |
| `--backbone-dir` | yes | `string` | `None` | — | — |
| `--candidate-catalog` | yes | `string` | `None` | — | — |
| `--fixed-base-hits` | yes | `string` | `None` | — | — |
| `--fixed-base-hits-manifest` | yes | `string` | `None` | — | — |
| `--minimum-query-coverage` | no | `float` | `0.5` | — | — |
| `--maximum-evalue` | no | `float` | `1e-05` | — | — |
| `--output-dir` | yes | `string` | `None` | — | — |

## `ploidypatch evidence propagate-wgd-conflict-partners`

ploidypatch evidence propagate-wgd-conflict-partners

```text
usage: ploidypatch evidence propagate-wgd-conflict-partners
       [-h] --base-gff BASE_GFF --candidate-gff CANDIDATE_GFF --pool-decisions
       POOL_DECISIONS --prior-wgd-selection PRIOR_WGD_SELECTION
       --output-selection OUTPUT_SELECTION
```

| Argument | Required | Type | Default | Choices | Description |
|---|---:|---|---|---|---|
| `-h`<br>`--help` | no | `string` | `—` | — | show this help message and exit |
| `--base-gff` | yes | `string` | `None` | — | — |
| `--candidate-gff` | yes | `string` | `None` | — | — |
| `--pool-decisions` | yes | `string` | `None` | — | — |
| `--prior-wgd-selection` | yes | `string` | `None` | — | — |
| `--output-selection` | yes | `string` | `None` | — | — |

## `ploidypatch evidence reconcile-pav-catalogs`

ploidypatch evidence reconcile-pav-catalogs

```text
usage: ploidypatch evidence reconcile-pav-catalogs [-h] --catalog CATALOG
                                                   --output-tsv OUTPUT_TSV
                                                   [--min-support MIN_SUPPORT]
```

| Argument | Required | Type | Default | Choices | Description |
|---|---:|---|---|---|---|
| `-h`<br>`--help` | no | `string` | `—` | — | show this help message and exit |
| `--catalog` | yes | `string` | `None` | — | repeatable LABEL=PATH catalog input |
| `--output-tsv` | yes | `string` | `None` | — | — |
| `--min-support` | no | `int` | `2` | — | — |

## `ploidypatch evidence score-copy-candidates`

ploidypatch evidence score-copy-candidates

```text
usage: ploidypatch evidence score-copy-candidates [-h] --features FEATURES
                                                  --model-json MODEL_JSON
                                                  [--mask-feature-group {wgd_context,method_quality}]
                                                  --output-tsv OUTPUT_TSV
```

| Argument | Required | Type | Default | Choices | Description |
|---|---:|---|---|---|---|
| `-h`<br>`--help` | no | `string` | `—` | — | show this help message and exit |
| `--features` | yes | `string` | `None` | — | — |
| `--model-json` | yes | `string` | `None` | — | — |
| `--mask-feature-group` | no | `string` | `[]` | `wgd_context`, `method_quality` | repeatable counterfactual evidence mask; applies the same frozen model and thresholds without refitting |
| `--output-tsv` | yes | `string` | `None` | — | — |

## `ploidypatch evidence score-homeolog-copy-candidates`

ploidypatch evidence score-homeolog-copy-candidates

```text
usage: ploidypatch evidence score-homeolog-copy-candidates [-h]
                                                           --copy-features
                                                           COPY_FEATURES
                                                           --topology-features
                                                           TOPOLOGY_FEATURES
                                                           --model-json
                                                           MODEL_JSON
                                                           --output-tsv
                                                           OUTPUT_TSV
```

| Argument | Required | Type | Default | Choices | Description |
|---|---:|---|---|---|---|
| `-h`<br>`--help` | no | `string` | `—` | — | show this help message and exit |
| `--copy-features` | yes | `string` | `None` | — | — |
| `--topology-features` | yes | `string` | `None` | — | — |
| `--model-json` | yes | `string` | `None` | — | — |
| `--output-tsv` | yes | `string` | `None` | — | — |

## `ploidypatch evidence score-support-conditioned-candidates`

ploidypatch evidence score-support-conditioned-candidates

```text
usage: ploidypatch evidence score-support-conditioned-candidates
       [-h] --copy-features COPY_FEATURES --topology-features
       TOPOLOGY_FEATURES --model-json MODEL_JSON --output-tsv OUTPUT_TSV
```

| Argument | Required | Type | Default | Choices | Description |
|---|---:|---|---|---|---|
| `-h`<br>`--help` | no | `string` | `—` | — | show this help message and exit |
| `--copy-features` | yes | `string` | `None` | — | — |
| `--topology-features` | yes | `string` | `None` | — | — |
| `--model-json` | yes | `string` | `None` | — | — |
| `--output-tsv` | yes | `string` | `None` | — | — |

## `ploidypatch evidence screen-pav-proteins`

ploidypatch evidence screen-pav-proteins

```text
usage: ploidypatch evidence screen-pav-proteins [-h] --catalog CATALOG
                                                --miniprot-gff MINIPROT_GFF
                                                --protein-map PROTEIN_MAP
                                                --chromosome-pairs
                                                CHROMOSOME_PAIRS --output-tsv
                                                OUTPUT_TSV
                                                [--min-identity MIN_IDENTITY]
                                                [--min-query-coverage MIN_QUERY_COVERAGE]
                                                [--min-partial-query-coverage MIN_PARTIAL_QUERY_COVERAGE]
                                                [--max-local-distance-bp MAX_LOCAL_DISTANCE_BP]
                                                [--max-breakpoint-spread-bp MAX_BREAKPOINT_SPREAD_BP]
```

| Argument | Required | Type | Default | Choices | Description |
|---|---:|---|---|---|---|
| `-h`<br>`--help` | no | `string` | `—` | — | show this help message and exit |
| `--catalog` | yes | `string` | `None` | — | — |
| `--miniprot-gff` | yes | `string` | `None` | — | — |
| `--protein-map` | yes | `string` | `None` | — | — |
| `--chromosome-pairs` | yes | `string` | `None` | — | — |
| `--output-tsv` | yes | `string` | `None` | — | — |
| `--min-identity` | no | `float` | `0.5` | — | — |
| `--min-query-coverage` | no | `float` | `0.5` | — | — |
| `--min-partial-query-coverage` | no | `float` | `0.2` | — | — |
| `--max-local-distance-bp` | no | `int` | `100000` | — | — |
| `--max-breakpoint-spread-bp` | no | `int` | `1000` | — | — |

## `ploidypatch evidence select-projection-support`

ploidypatch evidence select-projection-support

```text
usage: ploidypatch evidence select-projection-support [-h] --candidate-gff
                                                      CANDIDATE_GFF
                                                      --projection-support
                                                      PROJECTION_SUPPORT
                                                      [--source-group-map SOURCE_GROUP_MAP]
                                                      [--min-support-group-count MIN_SUPPORT_GROUP_COUNT]
                                                      --output-gff OUTPUT_GFF
                                                      --selection-tsv
                                                      SELECTION_TSV
```

| Argument | Required | Type | Default | Choices | Description |
|---|---:|---|---|---|---|
| `-h`<br>`--help` | no | `string` | `—` | — | show this help message and exit |
| `--candidate-gff` | yes | `string` | `None` | — | — |
| `--projection-support` | yes | `string` | `None` | — | — |
| `--source-group-map` | no | `string` | `None` | — | — |
| `--min-support-group-count` | no | `int` | `2` | — | — |
| `--output-gff` | yes | `string` | `None` | — | — |
| `--selection-tsv` | yes | `string` | `None` | — | — |

## `ploidypatch evidence select-scored-copy-candidates`

ploidypatch evidence select-scored-copy-candidates

```text
usage: ploidypatch evidence select-scored-copy-candidates [-h] --base-gff
                                                          BASE_GFF
                                                          --candidate-gff
                                                          CANDIDATE_GFF
                                                          --scores SCORES
                                                          --model-json
                                                          MODEL_JSON
                                                          [--policy {review,high_confidence}]
                                                          --output-gff
                                                          OUTPUT_GFF
                                                          --selection-tsv
                                                          SELECTION_TSV
```

| Argument | Required | Type | Default | Choices | Description |
|---|---:|---|---|---|---|
| `-h`<br>`--help` | no | `string` | `—` | — | show this help message and exit |
| `--base-gff` | yes | `string` | `None` | — | — |
| `--candidate-gff` | yes | `string` | `None` | — | — |
| `--scores` | yes | `string` | `None` | — | — |
| `--model-json` | yes | `string` | `None` | — | — |
| `--policy` | no | `string` | `review` | `review`, `high_confidence` | — |
| `--output-gff` | yes | `string` | `None` | — | — |
| `--selection-tsv` | yes | `string` | `None` | — | — |

## `ploidypatch evidence select-synteny-gap-models`

ploidypatch evidence select-synteny-gap-models

```text
usage: ploidypatch evidence select-synteny-gap-models [-h] --gaps GAPS
                                                      --baseline-decisions
                                                      BASELINE_DECISIONS
                                                      --adapted-candidate-gff
                                                      ADAPTED_CANDIDATE_GFF
                                                      --output-selection
                                                      OUTPUT_SELECTION
                                                      --output-candidate-gff
                                                      OUTPUT_CANDIDATE_GFF
```

| Argument | Required | Type | Default | Choices | Description |
|---|---:|---|---|---|---|
| `-h`<br>`--help` | no | `string` | `—` | — | show this help message and exit |
| `--gaps` | yes | `string` | `None` | — | repeatable filtered synteny-gap TSV |
| `--baseline-decisions` | yes | `string` | `None` | — | — |
| `--adapted-candidate-gff` | yes | `string` | `None` | — | — |
| `--output-selection` | yes | `string` | `None` | — | — |
| `--output-candidate-gff` | yes | `string` | `None` | — | — |

## `ploidypatch evidence select-wgd-supported-candidates`

ploidypatch evidence select-wgd-supported-candidates

```text
usage: ploidypatch evidence select-wgd-supported-candidates
       [-h] --base-gff BASE_GFF --candidate-gff CANDIDATE_GFF --pairs PAIRS
       --output-gff OUTPUT_GFF --selection-tsv SELECTION_TSV
```

| Argument | Required | Type | Default | Choices | Description |
|---|---:|---|---|---|---|
| `-h`<br>`--help` | no | `string` | `—` | — | show this help message and exit |
| `--base-gff` | yes | `string` | `None` | — | — |
| `--candidate-gff` | yes | `string` | `None` | — | — |
| `--pairs` | yes | `string` | `None` | — | — |
| `--output-gff` | yes | `string` | `None` | — | — |
| `--selection-tsv` | yes | `string` | `None` | — | — |

## `ploidypatch evidence subset-pav-proteins`

ploidypatch evidence subset-pav-proteins

```text
usage: ploidypatch evidence subset-pav-proteins [-h] --catalog CATALOG
                                                --proteins PROTEINS
                                                --output-fasta OUTPUT_FASTA
```

| Argument | Required | Type | Default | Choices | Description |
|---|---:|---|---|---|---|
| `-h`<br>`--help` | no | `string` | `—` | — | show this help message and exit |
| `--catalog` | yes | `string` | `None` | — | — |
| `--proteins` | yes | `string` | `None` | — | — |
| `--output-fasta` | yes | `string` | `None` | — | — |

## `ploidypatch evidence summarize-natural`

ploidypatch evidence summarize-natural

```text
usage: ploidypatch evidence summarize-natural [-h] --validation-tsv
                                              VALIDATION_TSV --output-json
                                              OUTPUT_JSON
```

| Argument | Required | Type | Default | Choices | Description |
|---|---:|---|---|---|---|
| `-h`<br>`--help` | no | `string` | `—` | — | show this help message and exit |
| `--validation-tsv` | yes | `string` | `None` | — | — |
| `--output-json` | yes | `string` | `None` | — | — |

## `ploidypatch evidence summarize-projection-support`

ploidypatch evidence summarize-projection-support

```text
usage: ploidypatch evidence summarize-projection-support [-h] --decisions
                                                         DECISIONS
                                                         --output-tsv
                                                         OUTPUT_TSV
```

| Argument | Required | Type | Default | Choices | Description |
|---|---:|---|---|---|---|
| `-h`<br>`--help` | no | `string` | `—` | — | show this help message and exit |
| `--decisions` | yes | `string` | `None` | — | — |
| `--output-tsv` | yes | `string` | `None` | — | — |

## `ploidypatch evidence summarize-wgdi`

ploidypatch evidence summarize-wgdi

```text
usage: ploidypatch evidence summarize-wgdi [-h] --query-gff QUERY_GFF
                                           --query-input-manifest
                                           QUERY_INPUT_MANIFEST --collinearity
                                           COLLINEARITY --expected-source
                                           EXPECTED_SOURCE --output-gene-tsv
                                           OUTPUT_GENE_TSV --output-block-tsv
                                           OUTPUT_BLOCK_TSV
```

| Argument | Required | Type | Default | Choices | Description |
|---|---:|---|---|---|---|
| `-h`<br>`--help` | no | `string` | `—` | — | show this help message and exit |
| `--query-gff` | yes | `string` | `None` | — | — |
| `--query-input-manifest` | yes | `string` | `None` | — | — |
| `--collinearity` | yes | `string` | `None` | — | repeatable SOURCE=PATH value |
| `--expected-source` | yes | `string` | `None` | — | repeatable SUBGENOME=SOURCE mapping |
| `--output-gene-tsv` | yes | `string` | `None` | — | — |
| `--output-block-tsv` | yes | `string` | `None` | — | — |

## `ploidypatch evidence validate-isoseq-candidates`

ploidypatch evidence validate-isoseq-candidates

```text
usage: ploidypatch evidence validate-isoseq-candidates [-h] --candidate-gff
                                                       CANDIDATE_GFF --paf PAF
                                                       --selected-counts
                                                       SELECTED_COUNTS
                                                       --genome-fasta
                                                       GENOME_FASTA
                                                       --output-evidence
                                                       OUTPUT_EVIDENCE
                                                       [--minimum-query-coverage MINIMUM_QUERY_COVERAGE]
                                                       [--minimum-identity MINIMUM_IDENTITY]
                                                       [--minimum-mapq MINIMUM_MAPQ]
                                                       [--maximum-secondary-score-fraction MAXIMUM_SECONDARY_SCORE_FRACTION]
                                                       [--minimum-candidate-cds-coverage MINIMUM_CANDIDATE_CDS_COVERAGE]
                                                       [--flank-bp FLANK_BP]
                                                       [--alignment-strand-source {query_orientation,minimap2_ts}]
```

| Argument | Required | Type | Default | Choices | Description |
|---|---:|---|---|---|---|
| `-h`<br>`--help` | no | `string` | `—` | — | show this help message and exit |
| `--candidate-gff` | yes | `string` | `None` | — | — |
| `--paf` | yes | `string` | `None` | — | — |
| `--selected-counts` | yes | `string` | `None` | — | — |
| `--genome-fasta` | yes | `string` | `None` | — | — |
| `--output-evidence` | yes | `string` | `None` | — | — |
| `--minimum-query-coverage` | no | `float` | `0.9` | — | — |
| `--minimum-identity` | no | `float` | `0.98` | — | — |
| `--minimum-mapq` | no | `int` | `20` | — | — |
| `--maximum-secondary-score-fraction` | no | `float` | `0.95` | — | — |
| `--minimum-candidate-cds-coverage` | no | `float` | `0.9` | — | — |
| `--flank-bp` | no | `int` | `5000` | — | — |
| `--alignment-strand-source` | no | `string` | `query_orientation` | `query_orientation`, `minimap2_ts` | — |

## `ploidypatch evidence validate-natural-rna`

ploidypatch evidence validate-natural-rna

```text
usage: ploidypatch evidence validate-natural-rna [-h] --candidates CANDIDATES
                                                 --junctions JUNCTIONS
                                                 --output-tsv OUTPUT_TSV
```

| Argument | Required | Type | Default | Choices | Description |
|---|---:|---|---|---|---|
| `-h`<br>`--help` | no | `string` | `—` | — | show this help message and exit |
| `--candidates` | yes | `string` | `None` | — | — |
| `--junctions` | yes | `string` | `None` | — | — |
| `--output-tsv` | yes | `string` | `None` | — | — |

## `ploidypatch evidence validate-natural-secondary-rna`

ploidypatch evidence validate-natural-secondary-rna

```text
usage: ploidypatch evidence validate-natural-secondary-rna [-h]
                                                           --primary-validation
                                                           PRIMARY_VALIDATION
                                                           --grouped-junctions
                                                           GROUPED_JUNCTIONS
                                                           --output-tsv
                                                           OUTPUT_TSV
```

| Argument | Required | Type | Default | Choices | Description |
|---|---:|---|---|---|---|
| `-h`<br>`--help` | no | `string` | `—` | — | show this help message and exit |
| `--primary-validation` | yes | `string` | `None` | — | — |
| `--grouped-junctions` | yes | `string` | `None` | — | — |
| `--output-tsv` | yes | `string` | `None` | — | — |

## `ploidypatch graph infer`

ploidypatch graph infer

```text
usage: ploidypatch graph infer [-h] --candidates CANDIDATES --evidence
                               EVIDENCE --output-json OUTPUT_JSON
                               --decisions-tsv DECISIONS_TSV
```

| Argument | Required | Type | Default | Choices | Description |
|---|---:|---|---|---|---|
| `-h`<br>`--help` | no | `string` | `—` | — | show this help message and exit |
| `--candidates` | yes | `string` | `None` | — | — |
| `--evidence` | yes | `string` | `None` | — | — |
| `--output-json` | yes | `string` | `None` | — | — |
| `--decisions-tsv` | yes | `string` | `None` | — | — |

## `ploidypatch normalize ncbi-primary`

ploidypatch normalize ncbi-primary

```text
usage: ploidypatch normalize ncbi-primary [-h] --gff GFF --protein PROTEIN
                                          --cds CDS --genome GENOME
                                          --output-dir OUTPUT_DIR
```

| Argument | Required | Type | Default | Choices | Description |
|---|---:|---|---|---|---|
| `-h`<br>`--help` | no | `string` | `—` | — | show this help message and exit |
| `--gff` | yes | `string` | `None` | — | — |
| `--protein` | yes | `string` | `None` | — | — |
| `--cds` | yes | `string` | `None` | — | — |
| `--genome` | yes | `string` | `None` | — | — |
| `--output-dir` | yes | `string` | `None` | — | — |

## `ploidypatch normalize primary-annotation`

ploidypatch normalize primary-annotation

```text
usage: ploidypatch normalize primary-annotation [-h] --gff GFF --genome GENOME
                                                --output-dir OUTPUT_DIR
                                                [--primary-seqid-table PRIMARY_SEQID_TABLE]
```

| Argument | Required | Type | Default | Choices | Description |
|---|---:|---|---|---|---|
| `-h`<br>`--help` | no | `string` | `—` | — | show this help message and exit |
| `--gff` | yes | `string` | `None` | — | — |
| `--genome` | yes | `string` | `None` | — | — |
| `--output-dir` | yes | `string` | `None` | — | — |
| `--primary-seqid-table` | no | `string` | `None` | — | optional strict TSV with seqid and chromosome_label columns; use for provider GFF files without chromosome declaration records |

## `ploidypatch normalize provider-gff3`

ploidypatch normalize provider-gff3

```text
usage: ploidypatch normalize provider-gff3 [-h] --gff GFF --output-dir
                                           OUTPUT_DIR
                                           [--repair-unescaped-note-semicolons]
                                           [--drop-invalid-intron-intervals]
                                           [--strip-embedded-fasta]
```

| Argument | Required | Type | Default | Choices | Description |
|---|---:|---|---|---|---|
| `-h`<br>`--help` | no | `string` | `—` | — | show this help message and exit |
| `--gff` | yes | `string` | `None` | — | — |
| `--output-dir` | yes | `string` | `None` | — | — |
| `--repair-unescaped-note-semicolons` | no | `flag` | `false` | — | percent-escape bare semicolon continuations only inside Note values |
| `--drop-invalid-intron-intervals` | no | `flag` | `false` | — | drop negative-length intron records without changing exon or CDS records |
| `--strip-embedded-fasta` | no | `flag` | `false` | — | stop at ##FASTA when a separate genome FASTA is supplied downstream |

## `ploidypatch patch apply`

ploidypatch patch apply

```text
usage: ploidypatch patch apply [-h] --source-gff SOURCE_GFF --patch PATCH
                               --output-gff OUTPUT_GFF
```

| Argument | Required | Type | Default | Choices | Description |
|---|---:|---|---|---|---|
| `-h`<br>`--help` | no | `string` | `—` | — | show this help message and exit |
| `--source-gff` | yes | `string` | `None` | — | — |
| `--patch` | yes | `string` | `None` | — | — |
| `--output-gff` | yes | `string` | `None` | — | — |

## `ploidypatch patch compile-copy-additions`

ploidypatch patch compile-copy-additions

```text
usage: ploidypatch patch compile-copy-additions [-h] --annotation-gff
                                                ANNOTATION_GFF --candidate-gff
                                                CANDIDATE_GFF
                                                --output-edits-json
                                                OUTPUT_EDITS_JSON
```

| Argument | Required | Type | Default | Choices | Description |
|---|---:|---|---|---|---|
| `-h`<br>`--help` | no | `string` | `—` | — | show this help message and exit |
| `--annotation-gff` | yes | `string` | `None` | — | — |
| `--candidate-gff` | yes | `string` | `None` | — | — |
| `--output-edits-json` | yes | `string` | `None` | — | — |

## `ploidypatch patch compile-reviewed-copy-additions`

ploidypatch patch compile-reviewed-copy-additions

```text
usage: ploidypatch patch compile-reviewed-copy-additions [-h] --annotation-gff
                                                         ANNOTATION_GFF
                                                         --candidate-gff
                                                         CANDIDATE_GFF
                                                         --pool-decisions
                                                         POOL_DECISIONS
                                                         --pool-manifest
                                                         POOL_MANIFEST
                                                         --review-decisions
                                                         REVIEW_DECISIONS
                                                         --output-edits-json
                                                         OUTPUT_EDITS_JSON
```

| Argument | Required | Type | Default | Choices | Description |
|---|---:|---|---|---|---|
| `-h`<br>`--help` | no | `string` | `—` | — | show this help message and exit |
| `--annotation-gff` | yes | `string` | `None` | — | — |
| `--candidate-gff` | yes | `string` | `None` | — | — |
| `--pool-decisions` | yes | `string` | `None` | — | — |
| `--pool-manifest` | yes | `string` | `None` | — | — |
| `--review-decisions` | yes | `string` | `None` | — | — |
| `--output-edits-json` | yes | `string` | `None` | — | — |

## `ploidypatch patch compile-structure`

ploidypatch patch compile-structure

```text
usage: ploidypatch patch compile-structure [-h] --annotation-gff
                                           ANNOTATION_GFF --candidate-gff
                                           CANDIDATE_GFF --hypotheses-tsv
                                           HYPOTHESES_TSV --output-edits-json
                                           OUTPUT_EDITS_JSON --event-type
                                           {annotation_boundary_shift,annotation_fused_gene,annotation_missing_internal_exon,annotation_split_gene}
                                           [--min-support-group-count MIN_SUPPORT_GROUP_COUNT]
```

| Argument | Required | Type | Default | Choices | Description |
|---|---:|---|---|---|---|
| `-h`<br>`--help` | no | `string` | `—` | — | show this help message and exit |
| `--annotation-gff` | yes | `string` | `None` | — | — |
| `--candidate-gff` | yes | `string` | `None` | — | — |
| `--hypotheses-tsv` | yes | `string` | `None` | — | — |
| `--output-edits-json` | yes | `string` | `None` | — | — |
| `--event-type` | yes | `string` | `None` | `annotation_boundary_shift`, `annotation_fused_gene`, `annotation_missing_internal_exon`, `annotation_split_gene` | — |
| `--min-support-group-count` | no | `int` | `2` | — | — |

## `ploidypatch patch create`

ploidypatch patch create

```text
usage: ploidypatch patch create [-h] --source-gff SOURCE_GFF --edits-json
                                EDITS_JSON --output-patch OUTPUT_PATCH
```

| Argument | Required | Type | Default | Choices | Description |
|---|---:|---|---|---|---|
| `-h`<br>`--help` | no | `string` | `—` | — | show this help message and exit |
| `--source-gff` | yes | `string` | `None` | — | — |
| `--edits-json` | yes | `string` | `None` | — | — |
| `--output-patch` | yes | `string` | `None` | — | — |

## `ploidypatch patch revert`

ploidypatch patch revert

```text
usage: ploidypatch patch revert [-h] --patched-gff PATCHED_GFF --patch PATCH
                                --output-gff OUTPUT_GFF
```

| Argument | Required | Type | Default | Choices | Description |
|---|---:|---|---|---|---|
| `-h`<br>`--help` | no | `string` | `—` | — | show this help message and exit |
| `--patched-gff` | yes | `string` | `None` | — | — |
| `--patch` | yes | `string` | `None` | — | — |
| `--output-gff` | yes | `string` | `None` | — | — |

## `ploidypatch report`

ploidypatch report

```text
usage: ploidypatch report [-h] --candidate-gff CANDIDATE_GFF --pool-decisions
                          POOL_DECISIONS --pool-manifest POOL_MANIFEST
                          --output-dir OUTPUT_DIR
                          [--review-decisions REVIEW_DECISIONS]
                          [--scores SCORES] [--copy-features COPY_FEATURES]
                          [--topology-features TOPOLOGY_FEATURES]
                          [--patch-edits PATCH_EDITS]
                          [--run-summary RUN_SUMMARY] [--title TITLE]
                          [--max-embedded-candidates MAX_EMBEDDED_CANDIDATES]
                          [--fail-on-attention]
```

| Argument | Required | Type | Default | Choices | Description |
|---|---:|---|---|---|---|
| `-h`<br>`--help` | no | `string` | `—` | — | show this help message and exit |
| `--candidate-gff` | yes | `string` | `None` | — | chain-preserving consensus GFF3 containing consensus_digest genes |
| `--pool-decisions` | yes | `string` | `None` | — | candidate-pool decisions TSV bound by the pool manifest |
| `--pool-manifest` | yes | `string` | `None` | — | ploidypatch.method_candidate_pool.v2 JSON manifest |
| `--output-dir` | yes | `string` | `None` | — | new non-overwriting directory for HTML, JSON, TSV, and checksums |
| `--review-decisions` | no | `string` | `None` | — | optional explicit human review ledger TSV |
| `--scores` | no | `string` | `None` | — | optional exact-universe rank-score TSV |
| `--copy-features` | no | `string` | `None` | — | optional exact-universe copy-feature TSV for method overlap |
| `--topology-features` | no | `string` | `None` | — | optional exact-universe topology TSV for applicability context |
| `--patch-edits` | no | `string` | `None` | — | optional reviewed patch-edits JSON for accept/event consistency |
| `--run-summary` | no | `string` | `None` | — | optional run summary proving byte-identical reversion |
| `--title` | no | `string` | `Annotation repair review` | — | human-readable project or sample title |
| `--max-embedded-candidates` | no | `int` | `5000` | — | maximum ranked candidates embedded in the interactive HTML; report.json and candidates.tsv always retain the full universe |
| `--fail-on-attention` | no | `flag` | `false` | — | exit with status 2 when the generated report requires attention |
