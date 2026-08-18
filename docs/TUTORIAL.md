# PloidyPatch v1.0 tutorial

This tutorial walks through the complete reviewed-copy workflow using the
small dataset in `examples/minimal_reviewed_patch/`. It requires only Python
and the installed PloidyPatch core.

## Learning objectives

You will:

1. inspect a checksum-bound candidate pool;
2. understand conflict and review records;
3. compile only explicit acceptances;
4. create and apply a source-bound GFF3 patch; and
5. prove byte-identical reversion; and
6. inspect the offline report and resume review from its ledger template.

## 1. Install and inspect the CLI

```bash
python -m pip install .
ploidypatch --version
ploidypatch patch --help
```

The expected version is `1.0.0`.

## 2. Understand the example inputs

| File | Important fields | Purpose |
|---|---|---|
| `source.gff3` | normal GFF3 columns | immutable source annotation |
| `candidate.gff3` | candidate digest in attributes | proposed copy structures |
| `decisions.tsv` | digest, status, conflict group | frozen candidate-pool decisions |
| `candidate.gff3.manifest.json` | policy and SHA-256 lineage | binds the pool universe |
| `review_decisions.tsv` | digest and review decision | explicit human ledger |
| `report_scores.tsv` | digest, rank, percentile, automatic flag | illustrative inspection priority |
| `report_copy_features.tsv` | exact method-presence columns | method intersections |
| `report_topology.tsv` | availability, reason, partner | optional topology context |
| `input_manifest.tsv` | bytes and SHA-256 | detects any changed input |

Candidates `a` and `b` are mutually exclusive. The ledger accepts `a`, rejects
`b`, and independently accepts `c`. A missing decision is never interpreted as
acceptance.

## 3. Run the verified end-to-end example

```bash
python examples/minimal_reviewed_patch/run_example.py \
  --output-dir work/tutorial-reviewed-patch
```

Inspect the summary:

```bash
python -m json.tool work/tutorial-reviewed-patch/run_summary.json
```

The important fields are:

```json
{
  "accepted_additions": 2,
  "automatic_approval": false,
  "byte_identical_reversion": true,
  "selected_candidate_digests": ["aaaa...", "cccc..."]
}
```

## 4. Run each CLI stage manually

### 4.1 Compile reviewed additions

```bash
ploidypatch patch compile-reviewed-copy-additions \
  --annotation-gff examples/minimal_reviewed_patch/source.gff3 \
  --candidate-gff examples/minimal_reviewed_patch/candidate.gff3 \
  --pool-decisions examples/minimal_reviewed_patch/decisions.tsv \
  --pool-manifest examples/minimal_reviewed_patch/candidate.gff3.manifest.json \
  --review-decisions examples/minimal_reviewed_patch/review_decisions.tsv \
  --output-edits-json work/tutorial.edits.json
```

This command checks the candidate universe, manifest lineage, review coverage,
conflict consistency, and decision vocabulary before writing addition events.

### 4.2 Create a content-bound patch

```bash
ploidypatch patch create \
  --source-gff examples/minimal_reviewed_patch/source.gff3 \
  --edits-json work/tutorial.edits.json \
  --output-patch work/tutorial.patch.json
```

The patch stores the expected source digest. Applying it to a different file
must fail rather than guess.

### 4.3 Apply without overwriting the source

```bash
ploidypatch patch apply \
  --source-gff examples/minimal_reviewed_patch/source.gff3 \
  --patch work/tutorial.patch.json \
  --output-gff work/tutorial.patched.gff3
```

Confirm that accepted candidates `PPCONS_gene_a` and `PPCONS_gene_c` are
present and rejected `PPCONS_gene_b` is absent.

### 4.4 Revert and compare hashes

```bash
ploidypatch patch revert \
  --patched-gff work/tutorial.patched.gff3 \
  --patch work/tutorial.patch.json \
  --output-gff work/tutorial.reverted.gff3
```

Linux/macOS:

```bash
sha256sum examples/minimal_reviewed_patch/source.gff3 \
  work/tutorial.reverted.gff3
```

PowerShell:

```powershell
Get-FileHash examples/minimal_reviewed_patch/source.gff3 -Algorithm SHA256
Get-FileHash work/tutorial.reverted.gff3 -Algorithm SHA256
```

The digests must match exactly.

## 5. Audit a real annotation bundle

```bash
ploidypatch audit \
  --gff annotation.gff3.gz \
  --protein proteins.fa.gz \
  --cds cds.fa.gz \
  --fai genome.fa.fai \
  --checksums \
  --output audit.json
```

The audit checks internal identifier, hierarchy, phase, translation, and
coordinate consistency. Passing does not establish biological completeness;
it establishes that the bundle is internally suitable for the declared
workflow.

## 6. Move from the tutorial to a real project

1. Keep original annotation bytes read-only.
2. Normalize each upstream tool independently and preserve source-gene IDs.
3. Build the candidate pool under one explicit overlap policy.
4. Freeze candidate and decision manifests before review.
5. Rank for review, but do not convert score thresholds into approval.
6. Record one explicit decision per candidate.
7. Compile, apply, revert, and checksum the patch.
8. Archive the source, candidate pool, review ledger, patch, command log, and
   output manifest together.

For all 87 commands and every option, see the generated
[parameter reference](CLI_PARAMETERS.md). For formal blind validation and
truth custody, continue with the [reproducibility guide](REPRODUCIBILITY_GUIDE.md).

The verified example now writes `output/report/index.html` automatically. Open
it to use the triage queues, exact method intersections, readable rank,
shared-axis conflict comparison, phased CDS structures, topology explanation,
and checksum provenance. To start a fresh review, edit the adjacent
`review_ledger_template.tsv`, regenerate the report with `--review-decisions`,
and never treat its score as an approval threshold. See the
[result-report guide](RESULT_REPORT.md) for the full workflow and scale limits.
