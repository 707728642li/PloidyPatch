# PloidyPatch controlled benchmark protocol v0.1

## Purpose and scope

This protocol tests whether a method can recover deliberately hidden plant gene
models without confusing annotation defects with biological gene loss. It is a
development protocol, not evidence of method superiority by itself.

Version 0.1 freezes the file separation, provenance rules, deterministic
selection rule, and the first smoke-test event
`annotation_missing_gene`. Changes to those contracts require a new protocol
version. The complete event vocabulary is machine-readable in
`config/event_ontology_v0.1.tsv`.

## Four physically separate artifacts

1. **Immutable source**: a QC-passed reference GFF3. It is never edited.
2. **Blind input**: `perturbed.gff3`, containing only the altered annotation.
3. **Public run manifest**: `manifest.json`, containing source/output hashes,
   generator version, event counts, eligibility count, seed, and selection
   algorithm. It does not disclose target loci.
4. **Evaluator-only truth**: `hidden_truth.json`, containing target loci and
   exact removed records. Candidate methods must not be able to read this file.

Training, tuning, or running a candidate in the same directory as the hidden
truth is protocol leakage. Production workflows must place truth in a separate
permission domain or evaluator process. A checksum in the public manifest is
for custody verification and conveys no target identity.

## Implemented event: annotation missing gene

The genome is unchanged and an intact coding gene hierarchy is removed from
the annotation. This must never be labelled as a biological deletion.

Eligibility is conservative in v0.1:

- the gene has a unique `ID`;
- its descendant closure includes at least one transcript, exon, and CDS;
- every descendant feature has exactly one parent;
- removing the hierarchy therefore cannot alter a shared feature.

Eligible genes are ranked by SHA-256 over the fixed seed and stable locus
identity. The lowest ranks are selected. This avoids dependence on input hash
iteration order and makes the same source/seed/count deterministic.

Each hidden-truth event stores a stable event ID, locus, transcript IDs,
feature counts, original line numbers, raw GFF3 records, and per-line hashes.
The inverse operation reconstructs the decompressed source text and must match
its SHA-256 exactly. Reversibility is an implementation invariant, not a
biological success metric.

## Smoke-test command

```bash
ploidypatch benchmark perturb \
  --gff source.gff3.gz \
  --output-dir benchmark_runs/ath_missing_gene_seed17 \
  --event-type annotation_missing_gene \
  --count 100 \
  --seed 17
```

Evaluator-only restoration check:

```bash
ploidypatch benchmark restore \
  --perturbed-gff benchmark_runs/ath_missing_gene_seed17/perturbed.gff3 \
  --truth benchmark_runs/ath_missing_gene_seed17/hidden_truth.json \
  --output-gff benchmark_runs/ath_missing_gene_seed17/restored.gff3
```

Output paths are exclusive: the generator and restorer refuse to overwrite an
existing artifact.

## Plant-specific stratification required before comparative results

The current deterministic sampler is only a software smoke test. Comparative
experiments must first construct a candidate catalog and sample within frozen
strata for:

- retained WGD/homeolog, tandem, proximal, transposed, and dispersed duplicate
  state;
- subgenome and homoeolog dosage where defined;
- syntenic context and local copy density;
- family size and orthogroup ambiguity;
- gene length, transcript count, exon count, and CDS length;
- repeats/TE proximity, mappability, and assembly gaps;
- evidence availability, including RNA-seq and long-read support.

The perturbation generator may use a relationship engine to label strata, but
the evaluated method must not receive hidden target labels. Any evidence used
to create a candidate must be declared, and scoring should preferentially use
independent evidence.

## Leakage and dataset partitions

- Freeze species and clade roles in `config/benchmark_dataset_plan.tsv` before
  fitting thresholds.
- Hold out entire species and orthogroups for final tests.
- Do not project a model into the target and then score it solely against that
  same projection.
- Remove target-derived proteins from strict homology reference pools.
- Keep the hawthorn genomes as a discovery/engineering cohort unless an
  independently validated truth set is available.
- Do not inspect hidden truth while tuning a candidate method; only aggregate
  evaluator metrics may return to the tuning loop.

## Planned evaluation contract

Candidate methods will emit repaired GFF3 plus a machine-readable change table.
The evaluator will report, at minimum:

- exact event recovery precision, recall, and F1;
- exact transcript, exon, CDS, and splice-boundary recovery;
- false repairs outside perturbed loci;
- per-event-type precision-recall curves and calibrated abstention;
- paired effects within WGD/duplicate strata;
- CPU hours, peak memory, disk use, failure rate, and excluded loci.

Confidence intervals must be block- or orthogroup-bootstrap intervals. Final
claims require held-out natural inconsistencies and plant duplication stress
tests in addition to controlled perturbations.

### Implemented strict scorer

The v0.1 scorer compares identifier-independent transcript signatures composed
of sequence ID, strand, transcript span, exact exon intervals, and exact CDS
intervals plus phase. Candidate structures absent from the perturbed input are
treated as proposed changes. Exact hidden structures are true positives; all
other novel structures are false positives. It also counts perturbed baseline
structures that disappear from the candidate as collateral losses.

The scorer fails closed if the source checksum differs from hidden truth or if
a hidden transcript structure is still present in the blind annotation. It
reports aggregate metrics by default. `--include-event-details` exposes target
IDs and is evaluator-only.

```bash
ploidypatch benchmark score \
  --source-gff /evaluator/source.gff3.gz \
  --perturbed-gff benchmark_runs/example_seed17/perturbed.gff3 \
  --candidate-gff /candidate/output.gff3 \
  --truth /evaluator/hidden_truth.json
```

## Version 0.1 limitations and next gate

Only `annotation_missing_gene` is implemented. It does not test split/fusion,
missing exons, splice boundaries, assembly defects, copy collapse, or the joint
WGD relationship model. The next valid engineering gate is:

1. pass real-GFF generation/restoration checks on diploid and polyploid input;
2. emit a non-leaking stratified candidate catalog;
3. extend the implemented exact scorer with boundary-level and translated-CDS
   partial-credit metrics;
4. reproduce expected behavior of at least one simple projection baseline on
   a small fixture before scaling up.

### Implementation checkpoint — 2026-08-07

Gates 1–3 now pass for the chromosome-scale Da-Ae `B. napus` reference.
The formal catalog contains 85,091 structurally unique candidates and joins
one-to-one to evaluator-only WGDI evidence. This checkpoint does not change the
v0.1 protocol and is not a repair-performance result. The remaining next gate
is a real projection baseline and negative controls. Scorer v2 now reports
exact exon segments, splice junctions, phase-aware CDS segments, and unioned
CDS nucleotide coverage in addition to strict transcript recovery.
