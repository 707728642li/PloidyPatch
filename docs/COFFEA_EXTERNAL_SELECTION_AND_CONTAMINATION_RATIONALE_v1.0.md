# Coffea external selection and contamination rationale v1.0

Status: metadata-only prospective selection record, 2026-08-09. No Coffea
WGDI pair, candidate pool, truth label, score, event count, or performance
result existed when this record was written.

## Purpose

This system is selected as the next untouched core-H1 replication after the
Walnut reveal was retained as an invalid run. The target is the recent
allotetraploid *Coffea arabica* ET-39 HiFi assembly. The biological event is
the known C/E-subgenome homoeolog pairing created by the origin of Arabica,
not an inferred generic ancient-WGD peak.

The deterministic metadata rule is:

1. chromosome-level target with a complete annotation and at least 95% of
   the full release assembly on primary nuclear chromosomes;
2. two candidate-only annotated references, at least one a diploid progenitor
   or close congener, each passing the same 95% rule;
3. two annotated evaluator-only Rubiaceae references outside the Arabica
   allopolyploid event, also passing 95%;
4. exact release binding for genome, GFF3 and protein, or an explicit
   assembly-report alias map with primary-sequence length identity; and
5. no project-history evidence that any selected species, pair, candidate,
   label, score, or performance result has previously been used.

No species may be selected by pair yield, candidate count, event count,
runtime, label distribution, H1 effect, or any downstream result.

## Selected role design

| Role | Release | Primary bp / full-release bp | Fraction | Rationale |
|---|---|---:|---:|---|
| target | *C. arabica* ET-39 HiFi, GCF_036785885.1 | 1,185,726,978 / 1,198,288,228 | 0.989517 | 22 labeled C/E chromosomes and current RefSeq annotation |
| candidate | *C. eugenioides* Bu-A v1 | 633,862,471 / 643,979,419 | 0.984290 | diploid E-subgenome progenitor; public Coffee Genome Hub annotation |
| candidate | *C. mauritiana* v1 | 406,702,695 / 426,289,761 | 0.954052 | independent diploid Coffea congener; avoids an easy same-cultivar arm |
| evaluator | *Gardenia jasminoides* ASM1310374v1 | 531,426,856 / 535,679,409 | 0.992061 | event-outside Rubiaceae group 1 |
| evaluator | *Ophiorrhiza pumila* v1 | 439,801,894 / 440,319,191 | 0.998825 | event-outside Rubiaceae group 2 |

The Bu-A Coffee Genome Hub pseudomolecules are exactly 100 bp longer per
chromosome than the corresponding GenBank submitted molecules. Therefore the
Hub genome, Hub GFF3 and Hub protein are treated as one inseparable release;
they must not be mixed with the shorter GenBank chromosome bytes. The
GenBank full-assembly total is used only for the conservative 95% metadata
audit.

Gardenia primary GFF3 coordinates and chromosome lengths exactly match the 11
assembly-report sequence names. The normalizer must use the frozen
`Gardenia1..Gardenia11` to `CM023095.1..CM023105.1` alias table. Twelve genes
on two unmatched non-primary contigs are excluded by the primary-only rule;
they cannot be fuzzy-mapped or used in formal truth.

## Rejected or deferred alternatives

- quinoa was rejected because only 90.5% of the reported assembly was on its
  18 chromosomes;
- old *C. arabica* Cara_1.0 was rejected at 90.45%;
- old RefSeq *C. eugenioides* Ceug_1.0 was rejected at 93.83%;
- *C. canephora* DH200-94 v2 was rejected at 87.09% when the denominator was
  correctly taken from the full assembly rather than the chromosome-only Hub
  download;
- *C. humblotiana* was rejected at 88.97%; and
- *C. arabica* Geisha passed the 95% assembly gate (96.76%) but was not
  selected because a same-species cultivar projection would make the primary
  test easier and less phylogenetically informative than *C. mauritiana*.

These exclusions are retained regardless of the eventual Coffea result.

## Project contamination audit

Before any Coffea pair or candidate calculation:

- `git grep` and full Git-history searches found no Coffea, Gardenia,
  Ophiorrhiza, ET-39, Bu-A, Geisha, or their selected accessions;
- path searches of existing server `results`, `data/derived`, `config`, and
  `docs` found no matching prior artifacts after excluding the new
  metadata-download directory; and
- no WGDI, DIAMOND, method-trio, pair, benchmark, label, score, or evaluation
  process has been run for the selected system.

The system is therefore eligible to be described as an untouched prospective
external H1 replication, subject to final byte hashes and syntax/ID preflight.

## Sources

- Arabica, Bu-A and DH200-94 comparative assemblies: [Nature Genetics
  2024](https://www.nature.com/articles/s41588-024-01695-w)
- Public Coffea releases and annotations: [Coffee Genome
  Hub](https://coffee-genome-hub.southgreen.fr/)
- Gardenia reference: [Plant Communications
  2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7302004/)
- Ophiorrhiza reference: [Communications Biology
  2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC7810986/)
