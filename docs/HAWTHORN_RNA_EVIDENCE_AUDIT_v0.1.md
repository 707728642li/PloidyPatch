# Hawthorn RNA evidence audit v0.1

Date frozen: 2026-08-07  
Source (read-only):
`/nas_data/NFS/Shanzha/shanzha_analysis/cp_genome/lijiajia_genome`  
Audit result:
`/data/codexli/projects/PloidyPatch/results/hawthorn/rna_audit_v0.1`

## Result

The collection contains 60 paired-end short-read RNA-seq samples, 60 cleaned
pairs, 60 coordinate-sorted BAMs, and 60 BAM indexes. Every BAM passes
`samtools quickcheck` and has the same ordered 1,699-sequence dictionary as
`Black_Primary.fa`.

The reference and annotation copied into the RNA pipeline are byte-identical
to the top-level Black_Primary bundle:

| Artifact | Bytes | SHA-256 |
|---|---:|---|
| Black_Primary genome | 762,971,615 | `db56ab135a49e09f01eba8f896969b183bc0a29f21492889b6562651c5573edc` |
| Black_Primary GFF3 | 92,929,116 | `6a3aaa90d213f7d75d6a0d4a1809a2e23b4a6e327f8eeec6ef35992d47d36bb6` |

The BAM sequence-dictionary SHA-256 is
`ed112e51a33218d0c3d5ddbc77e1ad95b6572783fa98ed81a4a368ff56837c6c`
for all 60 files. The audit used samtools 1.23.1 from the project-only
`ploidypatch-pav` prefix. Per-BAM bytes, index status, mapped/unmapped alignment
counts, sequence dictionary, header/program counts, and paths are frozen in
`bam_audit.tsv`; full headers are retained under `headers/`.

## Provenance limitations

- Every BAM has zero `@RG` records. Sample identity is therefore recoverable
  only from file names and pipeline logs, not embedded library/sample/tissue
  metadata.
- The supplied pipeline and logs show fastp 0.23.4, HISAT2, samtools, and
  featureCounts. `Black-1` reports 86.92% overall alignment; `JRY_S1-1`
  reports 82.99%. These are examples, not substitutes for the per-BAM audit.
- All samples—including JRY, RZ, Red, Yellow, BM, and XZ names—were mapped to
  Black_Primary. Their alignments cannot be treated as coordinate-valid
  evidence for the other eight assemblies without remapping or a separately
  validated coordinate lift.
- File names strongly suggest that `Black-1`, `Black-2`, and `Black-3` are the
  closest matched replicates for Black_Primary, but genotype and tissue cannot
  be proven from the available files alone. This is an explicit inference.
- A recursive filename audit found no Iso-Seq, FLNC, PacBio, CCS/HiFi, or
  Nanopore artifact. The collection therefore provides no independent
  full-length/long-read transcript validation.
- Short-read non-expression cannot prove absence, and cross-genotype mapping
  can suppress divergent junctions. Negative RNA evidence is never used to
  call true loss.

## Frozen evidence tiers

Natural validation will use Black_Primary coordinates only:

1. **Primary RNA support:** a splice junction observed with at least two
   uniquely usable split reads in at least two of `Black-1..3`.
2. **Secondary population support:** the same junction observed in at least
   two non-Black named samples, reported separately because genotype/reference
   mismatch can bias both presence and absence.
3. **RNA-unsupported:** no qualifying junction in the available tissues;
   this means missing evidence, never contradiction or biological absence.

Candidate discovery must not use the same junction thresholds later used for
validation. Black_Hap1/Hap2 homology projections and cross-assembly structure
discordance are the intended discovery sources; the BAM-derived junction
table is held back until proposals and source hashes are frozen.

## Long-read gate

The requested long-read arm is unavailable in this source collection. It
cannot be silently replaced by short-read assembly or by the annotation used
to generate candidates. A publication-grade natural truth set still requires
independent Iso-Seq/full-length cDNA or blinded manual review supported by new
data. Until then, hawthorn results must be labeled RNA-supported discovery,
not full-transcript ground truth.

## Reproducibility

`scripts/audit_hawthorn_rna.sh` performs the reference hashes, BAM quickcheck,
sequence-dictionary comparison, index/coordinate-sort checks, idxstats,
FASTQ-pair audit, and long-read filename search. It does not write to the
source tree and refuses to overwrite an existing result directory.

`ploidypatch evidence extract-bam-junctions` streams one coordinate BAM
through an explicitly supplied samtools executable. The frozen filters retain
MAPQ >= 20 primary mapped alignments and exclude SAM flags 4, 256, and 2048
(combined mask 2308). Every `N` CIGAR operation flanked by aligned reference
bases contributes one unstranded, one-based `(left_exon_end,
right_exon_start)` observation. The output manifest records the exact command,
reference-header hash, thresholds, counts, software version, elapsed time, and
TSV checksum.

`ploidypatch evidence aggregate-junctions` verifies every per-sample manifest
and checksum before separating the frozen Black primary tier from secondary
samples. The aggregate retains raw read and sample counts alongside Boolean
threshold calls and accepts repeatable input directories so primary evidence
can remain immutable while secondary samples are added later.
`scripts/run_hawthorn_primary_junctions.sh` launches the three
Black replicates concurrently, preserves failed work, refuses overwrite, and
publishes the result only by an atomic directory rename after all checks pass.
`scripts/run_hawthorn_secondary_junctions.sh` uses the audited sample/BAM table,
excludes the three Black files exactly, processes the remaining 57 samples in
bounded groups of eight, and emits a separate combined aggregate. Secondary
support never upgrades the primary genotype-matched tier.
The 57 non-primary files are not treated as 57 independent biological units.
`config/hawthorn_rna_sample_groups_v0.1.tsv` freezes 19 filename-stem groups,
each with three files. These stems are not interpreted as cultivar, tissue, or
time metadata. A separate grouped aggregate requires at least two reads in at
least two files per group and at least two non-Black groups before reporting
secondary group-level recurrence. Absence remains missing evidence rather than
contradiction.
