# PloidyPatch CLI reference

## Scope

PloidyPatch exposes eight top-level command families and 87 executable leaf
commands in the current source tree. The exact list is exported from the
production `argparse` parser to
[`CLI_COMMAND_INVENTORY_v0.1.json`](CLI_COMMAND_INVENTORY_v0.1.json) and is
checked for byte-for-byte reproducibility in the test suite.

```text
audit  baseline  benchmark  evidence  graph  normalize  patch  report
```

This document distinguishes the small supported review-and-patch path from the
larger research surface. Exact required options and defaults are authoritative
in command help:

```bash
ploidypatch --help
ploidypatch patch --help
ploidypatch patch compile-reviewed-copy-additions --help
```

The generated [complete parameter reference](CLI_PARAMETERS.md) contains one
section for every executable command, including usage, required arguments,
types, defaults, choices, and help text. It is regenerated directly from the
same parser in CI, so the documentation cannot silently drift from the CLI.

## Supported reviewed-patch path

| Command | Purpose | Decision boundary |
| --- | --- | --- |
| `ploidypatch audit` | Validate an annotation/protein/CDS bundle and optionally bind input checksums. | Reports consistency; does not call biological completeness. |
| `ploidypatch report` | Build a self-contained interactive review report from checksum-bound pool artifacts. | Explains and prioritizes; never converts rank into acceptance. |
| `ploidypatch patch compile-reviewed-copy-additions` | Convert explicit reviewer acceptances from a frozen candidate pool into addition edits. | Never approves a candidate itself. |
| `ploidypatch patch create` | Bind edits to immutable source bytes. | Refuses malformed or source-mismatched edits. |
| `ploidypatch patch apply` | Write a patched annotation to a new path. | Never overwrites the source annotation. |
| `ploidypatch patch revert` | Reverse a patch and recover the source bytes. | Reversion should be verified by SHA-256. |

The complete runnable contract is in
`examples/minimal_reviewed_patch/run_example.py`.

The report package, optional inputs, statuses, and interpretation boundary are
documented in [`RESULT_REPORT.md`](RESULT_REPORT.md).

## Additional patch command

- `patch compile-structure` compiles explicitly permitted exact-topology
  structure hypotheses such as internal-exon recovery, split/fused correction,
  or boundary repair.
- `patch compile-copy-additions` compiles an already prefiltered candidate GFF3
  and is retained for controlled workflows. Public use should prefer
  `compile-reviewed-copy-additions`, which verifies the pool and review ledger.

## Research command families

### `baseline`

Adapters and candidate-pool construction for miniprot and provider GFF3
outputs. These commands normalize source provenance, combine candidate GFF3
files, retain or suppress overlap according to an explicit policy, and infer
structure hypotheses. They require method-specific upstream outputs and are
not packaged substitutes for the upstream tools.

### `benchmark`

Controlled annotation perturbation, evaluator-only restoration and scoring,
catalog sampling, bootstrap intervals, copy-feature labeling, masking, and
abstention audits. Hidden truth must remain outside the candidate environment.
The `benchmark score` and restoration commands are evaluator operations, not
blind-method inputs.

### `evidence`

WGDI preparation and pair inference, projection-support aggregation,
duplication-aware feature construction, rank scoring, conflict-winner guards,
long-read validation, natural-candidate audits, and PAV reconciliation. This
family contains most research commands. Each formal run binds the relevant
protocol, policy, model, reference roles, and upstream executable versions.

### `normalize`

Provider-specific GFF3 and primary-annotation normalization. Normalization is
conservative: aliases must resolve uniquely, invalid translations are not
fabricated, and a normalization report must accompany derived files.

### `graph`

Transparent event-graph inference from candidate and evidence tables. Graph
decisions are evidence records, not automatic patch approval.

## Output conventions

Depending on the command, PloidyPatch writes JSON, TSV, GFF3, FASTA, or patch
JSON. Formal workflows additionally require:

- a manifest containing expected file size and SHA-256;
- a command log with exact arguments;
- a schema or policy identifier;
- explicit `truth_access` and role declarations when hidden evaluation data
  exist;
- adjacent `.working` output followed by atomic publication; and
- non-overwrite behavior.

CLI reports use machine-readable JSON where practical. Do not parse human help
text as a scientific result.

## Exporting the command inventory

From a source checkout:

```bash
PYTHONPATH=src python scripts/export_cli_inventory_v0.1.py \
  --output docs/CLI_COMMAND_INVENTORY_v0.1.json --check
```

Without `--check`, the exporter writes a new path and refuses to overwrite an
existing file. With `--check`, it compares generated bytes with the tracked
inventory and returns nonzero on any drift.
