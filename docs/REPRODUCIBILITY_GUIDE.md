# PloidyPatch reproducibility guide

## Reproducibility layers

PloidyPatch separates four layers that are often conflated in annotation
evaluation:

1. **source identity**: exact bytes, assembly/annotation release, primary
   sequence identifiers, and aliases;
2. **candidate execution**: upstream software, parameters, candidate universe,
   conflict policy, and blind outputs;
3. **evaluation custody**: truth construction, role isolation, reveal
   authorization, and formal result state; and
4. **reviewed patch**: human decisions, immutable edits, application, and
   byte-identical reversion.

A reproducible patch does not by itself prove biological correctness. A
biologically informative benchmark is not reproducible unless its bytes,
software, roles, and decision rules are bound before labels are opened.

## File identity

Every formal input and output is bound by byte count and SHA-256. Manifests use
canonical relative paths, reject symlinks where role isolation matters, and
cover an exact file universe rather than an arbitrary subset. `SHA256SUMS`
does not include itself.

For a distributable result directory:

```bash
sha256sum -c SHA256SUMS
```

On Windows PowerShell, individual hashes can be checked with
`Get-FileHash -Algorithm SHA256`; formal package tests perform exact-universe
verification in Python.

## Data roles and hidden truth

Prospective evaluations use disjoint roles:

- `shared_target`: target genome and truth-free perturbed annotation;
- `candidate_only`: references available to candidate-generation methods;
- `blind_benchmark`: sealed benchmark metadata without labels; and
- `evaluator_only`: source annotation, outgroup truth evidence, or labels.

Blind runs must not mount evaluator directories, a complete target annotation,
the public-genome store, or the network. The blind output tree is checksummed
and placed under custody before reveal. Evaluator code reads the frozen blind
scores only after authorization.

## Freeze before result-dependent computation

Formal protocols bind, before candidate or truth labels are inspected:

- species and release roles;
- source bytes and primary-sequence tables;
- event definition and multiplicity rules;
- requested event count and evaluability floors;
- endpoints, bootstrap units, replicates, seeds, and ordered hypotheses;
- candidate-pool and conflict policy;
- model or explicit no-ranker sentinel;
- software and environment locks; and
- invalid, not-evaluable, negative, and positive result policies.

Syntax-only checks such as gzip integrity, identifier uniqueness, GFF3 column
shape, coordinate validity, and source checksum verification are allowed before
the freeze. Candidate yield, truth-pair yield, labels, and performance are not
metadata-selection criteria.

## Result states

| State | Meaning | Permitted response |
| --- | --- | --- |
| `invalid_run` | Checksum, custody, isolation, schema, or execution contract failed. | Repair execution under a new recorded attempt; do not interpret biology. |
| `not_evaluable` | Frozen data sufficiency or sentinel criteria were not met. | Report insufficiency; do not relax rules for the same revealed system. |
| formal negative | Run was valid and evaluable but a scientific or safety gate failed. | Retain and report the result; do not tune it away. |
| formal positive | All ordered scientific and safety gates passed. | Report all endpoints, counts, intervals, and audits. |

The Actinidia result demonstrates why these states matter: its core structural
endpoint improved, but its composite outcome remains formally negative because
a preregistered topology-coverage gate failed. Walnut demonstrates a different
case: an orchestration/custody failure yielded no biological estimate and is
retained as invalid rather than negative.

## Environments

The core package is pure Python. External tools are installed in separate,
project-local conda prefixes and captured with explicit package locks and
executable/version manifests. Formal execution archives the exact source tree;
server mirrors without Git metadata are not treated as provenance by
themselves.

The development wheel and source distribution are rebuilt in CI. A clean
wheel-smoke environment installs with `--no-deps`, runs `pip check`, queries
the CLI version, and executes the checksum-bound minimal reviewed-patch
example.

## Reproducing manuscript evidence

The compact evidence needed to audit manuscript claims is tracked under
`manuscript/source_data/`. Each package has its own README or manifest and
checksums. Formal protocols and result interpretations are under `docs/`.
Large genomes, alignments, raw reads, and recomputable intermediate files are
not copied into the repository.

The claim registry maps every quantitative manuscript statement to a tracked
source field. Figure scripts consume tracked source data and produce a figure
manifest and checksums. The submission package binds the manuscript, figures,
cover letter, checklist, and source-data index.

For review and archival deposition, the deterministic archive
`manuscript/PloidyPatch_reviewer_reproducibility_v0.1.zip` collects the compact
evidence into a single 76-file exact universe. After extraction, run
`python verify_bundle.py .` from its root; the standard-library verifier rejects
missing, additional, renamed, or byte-modified files. `file_index.tsv` maps
every included file to its project-relative source. Large public genomes and
reads remain in the cited repositories and are represented by frozen source
manifests rather than duplicated in this archive.

## Minimum release audit

Before distributing a development build:

1. require a clean tracked source tree;
2. run the complete Python suite;
3. derive `SOURCE_DATE_EPOCH` from the clean source commit, independently build
   sdist and wheel twice with PEP 517 isolation, canonicalize sdist archive
   metadata, and require each pair to be byte-identical;
4. inspect the sdist for the example, core protocols, configs, and scripts;
5. install both wheel and sdist without runtime dependencies in separate clean
   prefixes;
6. run `pip check`, `ploidypatch --version`, and `ploidypatch --help` in each;
7. rerun the reviewed-patch example from each installed archive;
8. verify stable outputs agree, path-normalized command logs agree, automatic
   approval remains false, and reversion is byte-identical; and
9. record artifact bytes, SHA-256, canonicalization policy, Python version,
   source commit, and test totals.

The formal public archive records the final author identity, BSD-3-Clause
license, repository URL and persistent DOI in `CITATION.cff`. Future releases
must revalidate these fields rather than infer them from code or results.
