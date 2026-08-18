# Apple Golden Delicious natural ONT validation protocol v0.4

Status: candidate and ranking freeze completed before RNA-sequence access;
canonical Pychopper FLNC analysis running.  This is a descriptive natural-data
validation, not a new confirmatory model test.

## Question and claim boundary

This analysis asks whether evidence-blind PloidyPatch review ranks enrich for
candidate coding structures that are independently traversed by full-length
long reads in the same Golden Delicious genotype.  It does not use RNA data to
generate candidates, change candidate coordinates, tune the ranker, choose a
score threshold, or grant automatic annotation approval.

The 29,144-candidate complete-annotation control pool and all review ranks were
frozen at source commit `8f0b539921a383b87296dea4148c2738ca2ba08f` with
`RNA_access=false`, `validation_sequence_access=false`, and deterministic order
`descending_score_then_candidate_digest`.  The frozen v0.3 model SHA-256 is
`9dcafcc3294bfd85ccfbe824c33debd0e967578b89cf7f08b56076d32efd4390`.
Consequently, later read evidence cannot alter the candidate universe or its
ordering.  Because apple labels have now been inspected, all apple comparisons
are explicitly descriptive and cannot substitute for the untouched Populus
v0.4 confirmatory result.

## Public validation data

The raw input is GSA experiment CRA021523 from the 2026 Golden Delicious
long-read transcriptome study (DOI `10.3389/fpls.2026.1819201`).  The seven
official ONT cDNA-PCR FASTQ files are:

| Accession | Tissue |
|---|---|
| CRR1429911 | root |
| CRR1429912 | stem |
| CRR1429913 | leaf |
| CRR1429914 | flower |
| CRR1429915 | young fruit, 45 DAB |
| CRR1429916 | expanding fruit, 95 DAB |
| CRR1429917 | mature fruit, 145 DAB |

Every FASTQ is retained under its official MD5 contract and the entire input
tree is SHA-256 frozen.  Alignment uses the same GDDH13 v1.1 primary-chromosome
FASTA used for candidate generation.  Repeat screening uses the official
GDDH13 v1.1 TE GFF3.  The published reference-guided `All.gtf` is not used as
validation evidence because reference inheritance could make structural
agreement circular.

## Canonical FLNC preprocessing

1. Stream-validate FASTQ structure and retain reads of at least 500 nt.
2. Run unmodified Pychopper 2.7.10 with kit `PCS109`, method `phmm`, minimum
   mean quality 7, and minimum output segment length 50.  Quality is computed
   once by Pychopper from mean error probability.  Rescued multi-segment reads
   are counted but excluded from the canonical FLNC output.
3. Autotune the cutoff with NumPy seed `20261005`.  Fix `PYTHONHASHSEED=0`,
   `SOURCE_DATE_EPOCH=0`, and the executable `PATH`.  Pychopper 2.7.10 report
   generation is constrained to pandas 2.x; the conda explicit lock records
   the complete resolved environment.
4. Count the materialized FLNC FASTQ records independently of Pychopper's
   `Primers_found` statistic.  The latter can include a classified fragment
   removed by the 50-nt output filter, so both quantities are reported.

Two independent executions on the same 1,000-read real-data subset produced
identical SHA-256 hashes for FLNC FASTQ, statistics TSV, and report PDF after
the environment and timestamp controls were applied.

## Alignment and full-chain evidence

Each oriented FLNC set is aligned separately with minimap2 using
`-x splice -uf -k14 -G 1000000`, 10 retained secondary alignments, PAF CIGAR,
and long `cs` tags.  Since Pychopper has oriented reads 5-prime to 3-prime,
transcript strand is taken from PAF query orientation; minimap2 `ts:A` is not
used for this FLNC arm.

A read is eligible only if query coverage is at least 0.85, identity at least
0.90, MAPQ at least 20, and no secondary score reaches 0.95 of the selected
score.  Candidate CDS coverage must be at least 0.90 within a 5-kb locus flank.
A multi-exon candidate is `full_chain_supported` only when a single selected
read traverses the complete ordered intron/CDS chain on the correct strand.
Single-exon span support is reported separately and cannot establish complete
structure support.

## Fixed analyses

The pipeline reports:

- FLNC classification and exclusion counts for all seven tissues;
- evidence state and supporting-read/tissue counts for every one of 29,144
  frozen candidates;
- review yield at budgets 100, 146, 250, 292, 500, and 583;
- a 20,000-replicate fixed-seed comparison of v0.4 guard ranks with the frozen
  baseline ranks;
- exact transitions from the previously frozen raw-read/`ts:A` sensitivity
  analysis to the canonical Pychopper FLNC analysis; and
- a biological audit combining full-chain support, ORF/CDS plausibility,
  candidate self-mapping ambiguity, collision checks, and CDS/span/2-kb-flank
  overlap with the official TE annotation.

A manuscript case requires at least two canonical FLNC full-chain reads and a
multi-exon structure, then must also pass the sequence, repeat, collision, and
mapping audit.  Tissue breadth, method-family support, and WGD-partner context
are reported rather than converted into an automatic decision.

## Validity and negative-result policy

All source, ranking, model, candidate, raw-read, raw-alignment, TE, and code
manifests must verify before analysis.  The final universe must contain exactly
seven tissues and 29,144 candidates, and the PAF filter and validator must use
the same retained read universe.  Outputs are written to a non-overwriting
working directory, checksum-verified, made read-only, and atomically renamed.

The full-length counts and fractions reported by the source paper are
descriptive comparators, not gates.  A lower natural yield, no v0.4-over-v0.3
gain, or no case passing the strict audit is retained as a valid negative
result.  Thresholds, candidates, and ranks are not changed after observing
FLNC evidence.  No result from this protocol authorizes an automatic GFF
patch.
