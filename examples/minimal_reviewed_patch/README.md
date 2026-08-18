# Minimal reviewed-copy workflow

This distributable example exercises the public PloidyPatch safety path from
an immutable candidate-pool manifest and an explicit review ledger through a
reversible annotation patch. It intentionally contains one mutually exclusive
pair: the review accepts candidate `a`, rejects candidate `b`, and accepts the
independent candidate `c`.

From the repository root, run:

```bash
python examples/minimal_reviewed_patch/run_example.py
```

The runner first verifies every input byte against `input_manifest.tsv`, then
invokes the public CLI to:

1. compile only `accept_addition` decisions;
2. create an immutable patch;
3. apply the patch to the source GFF3; and
4. revert it and require byte-identical recovery of the source.

Successful output is written atomically to `output/`. The runner also creates
`output/report/index.html`, the complete offline review report, plus its
machine-readable JSON/TSV exports and an editable review-ledger template. The
summary must report
two accepted additions, `automatic_approval=false`, and identical source and
reverted SHA-256 values. `command_log.jsonl` records the five public CLI calls
and their machine-readable reports; `SHA256SUMS` binds every output file.

This small example covers PloidyPatch's review and reversible-patch core. The
publication workflows additionally run miniprot, GeMoMa, LiftOn, and WGDI in
separate frozen environments before producing the candidate-pool inputs shown
here.

The checksum-bound `report_scores.tsv`, `report_copy_features.tsv`, and
`report_topology.tsv` files demonstrate the interactive result-report views.
They are illustrative review evidence, not benchmark truth. See
`docs/RESULT_REPORT.md` for the complete schema and interpretation guide.
