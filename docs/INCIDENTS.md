# Reproducibility incidents

## 2026-08-08 - scaling display metadata named the wrong Brassica assembly

The first six-species core-scaling reference attempt used the correct Da-Ae
files, checksums, candidates, and feature matrices, but its human-readable
assembly field incorrectly said Darmor-bzh.  The mismatch was found while the
first workload was being reconstructed, before any formal timed scaling run.
The process was terminated and the partial tree was preserved as
`core_v0.1_references_74e8480.invalid_metadata`.

The corrected registry identifies the target as Da-Ae `GCF_020379485.1` with
NCBI Annotation Release 102 and also replaces the remaining shorthand display
labels with exact source releases.  No biological input, candidate policy,
model, endpoint, or performance result was changed.  A new commit, source
archive, reference freeze, and timed result are required; the partial attempt
is not eligible for reporting.

The next reference reconstruction completed successfully at commit `c4f8cbe`,
but a pre-timing completeness audit found that its registry recorded the
merged pool counts without separately recording each method-family decision
count or the number of conflict sets by size.  That byte-valid tree was kept as
`core_v0.1_references_c4f8cbe.superseded_protocol_metadata`; it is also not a
reportable timing reference.  The successor adds only provenance metadata and
must independently reproduce the same canonical bytes.

## 2026-08-08 - post-download Iso-Seq FASTA truncation

After the first Zenodo 2611319 evidence package had passed byte-count, MD5,
and SHA-256 validation, a remote read-only inspection command was misquoted by
the invoking shell and interpreted a FASTA-header search token as output
redirection. This truncated `F1maize.FINAL.fasta` to zero bytes. No candidate,
score, rank, or validation result had been generated from the damaged file.

The entire affected directory was preserved as
`maize_isophase_zenodo_2611319.invalid_accidental_shell_redirection_20260808T0130`.
The evidence package was then downloaded anew from the frozen URLs and again
validated against Zenodo byte counts and MD5 plus local SHA-256 before use.
The valid FASTA SHA-256 is
`e2deca471fcfe0eaf2b4935203b906abeacd9754060a4fe441ae579c5bcb68c8`.

Two otherwise valid statistic artifacts carrying incorrect manually supplied
code-commit labels were likewise preserved under names containing
`invalid_provenance`; correctly labeled reruns are the canonical artifacts.

## 2026-08-08 - repeated chromosome OOF seeds produced one partition

The first candidate-pool ranker v0.3 evaluation requested five shuffled
`StratifiedGroupKFold` seeds. Post-run audit showed that every seed produced the
same chromosome partition in both development species. The numerical result
was a valid single five-fold grouped evaluation but not a valid five-seed
repeated evaluation.

The complete artifact was preserved as
`candidate_pool_ranker_v0.3_evaluation.invalid_repeated_partition_88da631`.
Before any model or retention rule was changed, the evaluator was corrected to
use randomized `GroupKFold` with explicit assertions for five distinct
partitions, zero train/test chromosome overlap, and both labels in each fold.
The canonical `d1b42c8` rerun passed these audits and all predeclared retention
criteria.

## 2026-08-08 - first apple normalization exposed provider GFF defects

The first pre-truth apple external-input normalization stopped before pair
enumeration or event generation. Strict parsing found three unescaped
semicolon continuations in apple `Note` attributes, five negative-length
provider-supplied `intron` records in pear, and the standard `##FASTA`
terminator in the strawberry GFF3. Peach and rose normalized successfully.
The complete working directory was preserved as
`apple_v0.3.invalid_provider_gff_incompatibility_9f8622c`.

No gene, exon, CDS, pair, event, label, candidate score, or endpoint was
inspected or changed in response. A reusable compatibility normalizer was
added with all repairs disabled by default. The apple rerun may enable only
three narrow, manifest-recorded actions: percent-escape bare continuations of
`Note`, discard invalid redundant `intron` records, and strip an embedded FASTA
section when a separately checksummed genome FASTA is used. The rerun asserts
the exact frozen counts 3, 5, and 1 respectively and fails on any other defect.
