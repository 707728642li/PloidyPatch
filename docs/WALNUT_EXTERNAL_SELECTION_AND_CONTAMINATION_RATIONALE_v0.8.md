# Walnut v0.8 metadata selection and contamination rationale

Status: metadata-stage design record, 2026-08-09. This document does not
freeze a scientific protocol, model, pair universe, event truth set, or
performance endpoint.

## Selection boundary

*Juglans regia* Walnut 2.0 is the next untouched target selected for
reference-anchored development. Selection used only public assembly and
annotation metadata, raw FASTA/GFF3/protein syntax, phylogeny, and physical
role separation. No WGD pair yield, candidate count, truth label, event count,
model score, or performance endpoint was computed or inspected.

The historical 0.95 primary-chromosome base-fraction gate was not relaxed and
was applied to every member of the proposed 1/2/2 role system:

| Role | Species and release | Primary bp / assembly bp | Fraction |
| --- | --- | ---: | ---: |
| target | *J. regia* Walnut 2.0 | 545,778,456 / 572,786,133 | 0.9528485844348267 |
| candidate | *J. mandshurica* GWHBEUN00000000.1 | 528,151,479 / 548,694,591 | 0.9625600245802314 |
| candidate | *Carya illinoinensis* Pawnee GCF_018687715.1 | 674,272,338 / 674,433,157 | 0.9997615493865762 |
| evaluator | *Corylus avellana* CavTom2PMs-1.0 | 369,779,143 / 369,779,143 | 1.0 |
| evaluator | *Castanea mollissima* Nanking HAP2 GCA_053895285.1 | 704,228,712 / 706,692,489 | 0.9965136505080359 |

The first proposed replacements, *J. sigillata* and *Quercus lobata*, were
rejected rather than grandfathered: their exact fractions were
0.9299730706665554 and 0.8367602683815891. The three surveyed *Cucurbita*
targets were also rejected as formal holdouts because their public/NFS
assemblies did not meet 0.95. They remain eligible only for gene-centric
development that makes no chromosome-coverage claim.

Walnut 2.0 is chromosome-scale and has high assembly/gene completeness
reported by its primary publication. Comparative Juglandaceae analyses place
one juglandoid paleopolyploid event before the *Juglans*--*Carya* radiation;
the event is commonly interpreted as ancient allopolyploidy. This makes
*J. mandshurica* and *C. illinoinensis* event-sharing candidate references,
while Betulaceae (*Corylus*) and Fagaceae (*Castanea*) are physically separate
evaluator references outside that event. Public evidence:

- [Walnut 2.0 assembly](https://pmc.ncbi.nlm.nih.gov/articles/PMC7238675/)
- [Juglandoid WGD comparative analysis](https://pmc.ncbi.nlm.nih.gov/articles/PMC6431679/)
- [Juglandaceae paleopolyploid history](https://pmc.ncbi.nlm.nih.gov/articles/PMC9899254/)
- [NCBI Pawnee pecan assembly and annotation](https://www.ncbi.nlm.nih.gov/assembly/10213041)
- [NCBI Nanking chestnut BioProject](https://www.ncbi.nlm.nih.gov/bioproject/1190970)

Ancient allotetraploidy is a topology-risk declaration, not a truth label.
Subgenome fractionation and retained multi-copy regions can make an exact
one-to-two interpretation ambiguous. Any later fixed-target-backbone analysis
must quarantine ambiguous whole blocks, report applicability/coverage, and
abstain when the predeclared coverage gate is not met.

## Contamination audit

Immediately before this metadata record was added, case-insensitive tracked
worktree searches and `git log --all -S` searches returned no occurrence for
`Juglans_regia`, `Juglans_mandshurica`, `Carya_illinoinensis`,
`Corylus_avellana`, or `Castanea_mollissima`. None appears in the v0.4
contamination registry, development datasets, previous candidate/evaluator
roles, revealed labels, or result narratives.

The audit establishes target-level untouched status only as of this metadata
selection boundary. It does not claim that every later implementation detail
was specified before prior holdouts, and it must not be used to describe
Walnut v0.8 as an executed or frozen external validation.

## Source custody and current preflight

The exact 15 source artifacts, bytes, SHA-256 digests, public URLs, container
metadata, and primary-seqid tables are recorded in
`config/walnut_external_input_sources_v0.8.tsv`. The two NCBI bundles were
downloaded into the server project from their exact official directories.
All six gzip streams pass `gzip -t`, and their local MD5 values match the
official `md5checksums.txt` files.

Raw syntax checks found no duplicate FASTA seqids, invalid sequence symbols,
GFF seqids absent from the corresponding genome, malformed rows, invalid
coordinates, or primary nuclear CDS-parent/axis violations in the two NCBI
bundles. Pecan has 54,401 distinct primary-nuclear protein IDs and chestnut has
33,513; every one is present in its corresponding protein FASTA. Provider-wide
protein-ID sets also match exactly.

The *J. mandshurica* genome source is a tar-gzip container and must never be
renamed to `.fa.gz` or unpacked with `extractall`. Its only member is the
regular file `Juglans_mandshurica.genome.fa`, with 554,212,158 bytes and
SHA-256 `7b1be8bfb34096526ac94bc6f1806f0257c7b31293b81590225d283f8949679f`.
The generic fail-closed extractor and its tests are metadata-stage plumbing;
they are not a scientific protocol freeze.

## Explicitly deferred

- reference-anchored topology/backbone construction;
- WGD pair or event enumeration;
- candidate generation and candidate counts;
- truth construction, labels, scoring, and performance evaluation;
- final holdout contract, scientific policy, seeds, model, environment, or
  execution freeze.

Those items may be frozen only after the separately scoped reference-anchored
development work is complete.
