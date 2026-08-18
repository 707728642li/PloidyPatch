# Maize natural residual-annotation validation protocol v0.1

Date frozen: 2026-08-08  
Target assembly: `Zm-B73-REFERENCE-NAM-5.0`  
Target annotation: official `Zm00001eb.1` protein-coding annotation  
Status: design freeze before natural candidates are generated and before any
candidate is intersected with validation evidence

## Purpose and claim boundary

This experiment asks whether PloidyPatch can prioritize plausible residual
missing WGD/homeolog coding copies in an unperturbed, current maize annotation.
It is independent of the 800-event synthetic maize holdout and does not alter
that holdout, the frozen v0.2 model, or its success gates.

The primary natural endpoint is **independent-evidence yield per fixed manual
review budget**, not precision. A candidate without observed expression or a
pan-genome match is unlabeled, not a false positive. A supported candidate is
still a reviewable hypothesis unless its complete coding structure, locus, and
copy identity survive all declared checks. No annotation is patched
automatically.

## Evidence firewall and execution order

The following order is mandatory.

1. Finish and freeze the formal maize synthetic evaluation.
2. Run the already frozen homology-transfer methods on the complete official
   annotation. RNA, Iso-Seq, pan-genome, later annotation, and manual knowledge
   may not be read during candidate generation.
3. Infer self-WGD partners from the complete target annotation and score the
   union with the frozen baseline and topology rankers. No maize labels,
   thresholds, or feature changes are allowed.
4. Freeze every candidate, score, rank, method provenance, command, source
   checksum, code commit, and top-K membership in a write-once directory.
5. Only after step 4 may validation sequences be copied or downloaded,
   aligned, or intersected with candidates.
6. Preserve all unmatched candidates as `not_confirmed_by_available_evidence`;
   never relabel them as negatives.

Any operational correction after step 4 must preserve candidate coordinates,
scores, and ranks byte-for-byte. Otherwise the arm is invalidated and repeated
under a new protocol version with a new validation source.

## Discovery system

The discovery target is the unmodified official `Zm00001eb.1` annotation on
the ten primary NAM-5.0 chromosomes. Candidate models are produced from the
same Sorghum bicolor NCBIv3 and Setaria italica v2 references, program versions,
and commands frozen for the formal maize holdout:

- miniprot;
- LiftOn;
- GeMoMa; and
- exact two-of-three method-family consensus.

Existing target CDS chains are subtracted before ranking. The candidate union
is annotated with complete-annotation self-WGD evidence and scored once by:

- the frozen v0.2 baseline ranker; and
- the frozen v0.2 homeolog-topology ranker.

The natural experiment introduces no score threshold. Deterministic ranking is
descending model score, then candidate ID. The predeclared review budgets are
K = 25, 50, 100, and 200; K = 100 is primary. Individual mature methods and the
exact consensus are comparator rankings. A method's native alignment score is
used only inside that method; cross-method scores are never treated as
calibrated probabilities.

## Primary independent full-length transcript evidence

The first validation source is the published B73/Ki11 Iso-Seq resource
E-MTAB-7837 / ERP114624 and its processed Zenodo record 2611319. Candidate
discovery does not use this resource. The processed files are:

- `F1maize.FINAL.fasta` (nonredundant full-length transcript models);
- `F1maize.FINAL.demux_FL_count.txt` (per-model full-length read counts); and
- `F1maize.FINAL.gff` (B73 RefGen_v4 placement, context only).

Only pure B73 columns `EM1`, `R1`, and `END1` contribute primary biological
support. Ki11 and reciprocal-hybrid columns are excluded from the primary
support count and retained only as context. A transcript is B73-observed when
`EM1 + R1 + END1 >= 2`. This threshold is applied before looking at candidate
overlap.

Observed transcript sequences are aligned de novo to NAM-5.0 with minimap2
using the high-quality spliced-cDNA preset. A usable alignment must satisfy all
of the following:

- query coverage at least 0.90;
- nucleotide identity at least 0.98;
- mapping quality at least 20;
- no secondary alignment with at least 95% of the primary alignment score;
- no noncanonical splice junction unless the same junction has at least two
  independent full-length B73 reads; and
- no assembly `N` in the candidate interval or its predeclared 5-kb flanks.

Candidate evidence states are hierarchical:

1. `full_chain_supported`: one usable B73 transcript covers at least 90% of
   candidate CDS bases and matches every candidate intron boundary exactly;
2. `all_junctions_supported`: the union of usable B73 transcripts matches all
   candidate introns, but no single transcript establishes connectivity;
3. `partial_junction_support`;
4. `single_exon_span_support`: a usable transcript covers at least 90% of a
   single-exon candidate, which validates transcription but not gene identity;
5. `no_qualifying_observation`; or
6. `not_assessable`.

The primary positive state is `full_chain_supported`. Single-exon support and
union-only junction support are reported separately and cannot establish the
complete structure. Iso-Seq non-observation is never negative evidence.

The processed transcript set is a computational derivative of public raw
reads. Every primary supported locus selected for manuscript illustration must
also be checked against the raw ENA reads or a deposited per-read alignment so
that processing artifacts cannot create the case study.

## Independent maize pan-genome recurrence

After the candidate freeze, canonical proteins and PANDAGMA family assignments
from the 25 NAM founder annotations may be used as a second evidence family.
Candidate target CDS is translated from NAM-5.0, not copied from the discovery
reference. For each founder, the best and near-best protein hits are retained.

A candidate has `pan_genome_recurrent` context when:

- an intact target ORF exists without internal stop or frameshift;
- at least 10 non-B73 NAM founders contain a hit with amino-acid identity at
  least 0.95 and bidirectional coverage at least 0.90;
- hits from at least 80% of qualifying founders agree on one PANDAGMA family;
  and
- no alternative family or target locus is within 95% of the best score.

This is recurrence evidence, not truth: the NAM annotations share a pipeline
and biological absence is possible. Pan-genome recurrence is therefore kept
separate from Iso-Seq and never counted as 25 independent confirmations.

## Hard abstentions and locus audit

A candidate is ineligible for a manuscript correction claim if any of the
following remains unresolved:

- multi-locus near-best mapping;
- overlap with an existing same-strand CDS chain inconsistent with a missing
  copy;
- transposable-element/repeat-driven protein or transcript alignment;
- assembly gap or truncated 5-kb chromosome flank;
- non-intact translated ORF;
- incompatible gene grouping or mutually exclusive candidate models;
- the proposed model collides with another top-ranked edit; or
- the putative existing homeolog was itself used to manufacture the candidate
  relation rather than inferred leave-one-candidate-out.

All top-100 topology candidates receive the same automated audit. Manual
reviewers see coordinates and evidence tracks but are blinded to ranker name
and score. Review labels and reasons are frozen before ranker identity is
revealed.

## Primary and secondary statistics

The primary statistic is the difference in `full_chain_supported` loci among
the top 100 topology-ranked and top 100 baseline-ranked candidates. Report:

- supported loci per 100 reviews for each ranker;
- the observed difference;
- a 20,000-replicate chromosome-grouped bootstrap interval; and
- the exact overlap of the two top-100 sets.

The topology natural-direction gate passes only when the observed difference
is positive. A lower bootstrap bound above zero is strong confirmatory
evidence but is not required to retain the natural arm, because absent
expression is unlabeled and the number of observable loci is unknown in
advance.

Secondary analyses report the same yield at K = 25, 50, and 200, full support
plus pan-genome recurrence, results by method provenance, and enrichment over
20,000 chromosome-stratified random candidate rankings. No K or evidence tier
may be selected after seeing results. If fewer than ten full-chain-supported
loci occur in the union, the arm is descriptive and no comparative p-value is
interpreted.

## Retrospective annotation-revision arm

The archived preliminary `Zm00001e.1` annotation and later official
`Zm00001eb.1` annotation share the NAM-5.0 assembly. A later optional analysis
may ask whether RNA-blind PloidyPatch candidates against the preliminary
annotation were subsequently added in the official annotation. This is a
temporal silver standard only: MaizeGDB explicitly marks the preliminary
annotation as unsuitable for current use. It cannot be the headline natural
validation, cannot tune the model, and must be labeled separately from current
residual errors.

## Required frozen artifacts

The completed arm must publish:

- source and validation file contracts with URLs, accessions, sizes, and
  SHA-256 checksums;
- candidate and evidence firewalls with timestamps and code commits;
- complete union scores and top-K lists for every ranker/method;
- full-length transcript alignments and per-candidate evidence table;
- translated target ORFs and pan-genome recurrence table;
- blinded manual-review form and adjudication log;
- chromosome bootstrap replicates/seeds and machine-readable summaries; and
- resource-use records plus a one-command reproduction entry point.

